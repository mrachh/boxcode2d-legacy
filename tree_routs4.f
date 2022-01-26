c
c
c    common routines for generating and processing
c     a level restricted quad tree in 2D
c   
c
c
c
      subroutine tree_refine_boxes(irefinebox,nboxes,
     1  ifirstbox,nbloc,centers,bs,nbctr,nlctr,
     2  ilevel,iparent,nchild,ichild)
      implicit none
      integer nboxes,nbloc,nbctr,nlctr
      real *8 centers(2,nboxes),bs
      integer ilevel(nboxes),iparent(nboxes)
      integer ichild(4,nboxes),nchild(nboxes)
      integer irefinebox(nbloc)
      integer ifirstbox
      integer, allocatable :: isum(:)
      integer ii

      integer i,ibox,nel0,j,l,jbox,nel1,nbl

      allocate(isum(nbloc))
      if(nbloc.gt.0) call cumsum(nbloc,irefinebox,isum)
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,nbl,j,jbox)
      do i = 1,nbloc
        ibox = ifirstbox + i-1
        if(irefinebox(i).eq.1) then
          nbl = nbctr + (isum(i)-1)*4
          
          nchild(ibox) = 4
          do j=1,4
            ii = 2
            if(j.le.2) ii = 1
            jbox = nbl+j
            centers(1,jbox) = centers(1,ibox)+(-1)**j*bs/2
            centers(2,jbox) = centers(2,ibox)+(-1)**ii*bs/2 
            iparent(jbox) = ibox
            nchild(jbox) = 0
            do l=1,4
              ichild(l,jbox) = -1
            enddo
            ichild(j,ibox) = jbox
            ilevel(jbox) = nlctr 
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      

      if(nbloc.gt.0) nbctr = nbctr + isum(nbloc)*4


      return
      end
c
c
c
c
c
c
       subroutine tree_copy(nb,centers,ilevel,iparent,nchild,ichild,
     1              centers2,ilevel2,iparent2,nchild2,ichild2)

       implicit none
       integer nd,nb,npb
       real *8 centers(2,nb),centers2(2,nb)
       integer ilevel(nb),ilevel2(nb)
       integer iparent(nb),iparent2(nb)
       integer nchild(nb),nchild2(nb)
       integer ichild(4,nb),ichild2(4,nb)

       integer i,j,nel


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
       do i=1,nb
         centers2(1,i) = centers(1,i)
         centers2(2,i) = centers(2,i)
         ilevel2(i) = ilevel(i)
         iparent2(i) = iparent(i)
         nchild2(i) = nchild(i)
         ichild2(1,i) = ichild(1,i)
         ichild2(2,i) = ichild(2,i)
         ichild2(3,i) = ichild(3,i)
         ichild2(4,i) = ichild(4,i)
       enddo
C$OMP END PARALLEL DO       
       

       return
       end
c
c
c
c
c

      subroutine computecoll(nlevels,nboxes,laddr,boxsize,
     1                       centers,iparent,nchild,ichild,iper,
     2                       nnbors,nbors)

c     This subroutine computes the colleagues for an adaptive
c     pruned tree. box j is a colleague of box i, if they share a
c     vertex or an edge and the two boxes are at the same
c     level in the tree
c
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c
c     boxsize     in: double precision(0:nlevels)
c                 Array of boxsizes
c 
c     centers     in: double precision(2,nboxes)
c                 array of centers of boxes
c   
c     iparent     in: integer(nboxes)
c                 iparent(i) is the box number of the parent of
c                 box i
c
c     nchild      in: integer(nboxes)
c                 nchild(i) is the number of children of box i
c
c     ichild      in: integer(4,nboxes)
c                 ichild(j,i) is the box id of the jth child of
c                 box i
c
c     iper        in: integer
c                 flag for periodic implementations. 
c                 Currently not used. Feature under construction.
c
c----------------------------------------------------------------
c     OUTPUT
c     nnbors      out: integer(nboxes)
c                 nnbors(i) is the number of colleague boxes of
c                 box i
c
c     nbors       out: integer(9,nboxes)
c                 nbors(j,i) is the box id of the jth colleague
c                 box of box i
c---------------------------------------------------------------
      implicit none
      integer nlevels,nboxes
      integer iper
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(2,nboxes)
      integer iparent(nboxes), nchild(nboxes), ichild(4,nboxes)
      integer nnbors(nboxes)
      integer nbors(9,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,dad
      integer i,j,ifirstbox,ilastbox


c     Setting parameters for level = 0
      nnbors(1) = 1
      nbors(1,1) = 1
      do ilev = 1,nlevels
c        Find the first and the last box at level ilev      
         ifirstbox = laddr(1,ilev)
         ilastbox = laddr(2,ilev)
c        Loop over all boxes to evaluate neighbors, list1 and updating
c        hunglists of targets

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox)
         do ibox = ifirstbox,ilastbox
c           Find the parent of the current box         
            dad = iparent(ibox)
c           Loop over the neighbors of the parent box
c           to find out list 1 and list 2
            do i=1,nnbors(dad)
                jbox = nbors(i,dad)
                do j=1,4
c               ichild(j,jbox) is one of the children of the
c               neighbors of the parent of the current
c               box
                   kbox = ichild(j,jbox)
                   if(kbox.gt.0) then
c               Check if kbox is a nearest neighbor or in list 2
                      if((abs(centers(1,kbox)-centers(1,ibox)).le.
     1                   1.05*boxsize(ilev)).and.
     2                   (abs(centers(2,kbox)-centers(2,ibox)).le.
     3                   1.05*boxsize(ilev))) then
                     
                         nnbors(ibox) = nnbors(ibox)+1
                         nbors(nnbors(ibox),ibox) = kbox
                      endif
                   endif
                enddo
            enddo
c           End of computing colleagues of box i
         enddo
C$OMP END PARALLEL DO         
      enddo

      return
      end
c
c
c
c
c--------------------------------------------------------------------      
      subroutine updateflags(curlev,nboxes,nlevels,laddr,nchild,ichild,
     1                    nnbors,nbors,centers,boxsize,iflag)

c      This subroutine is to check the boxes flagged as flag++
c      and determine which of the boxes need refinement. The flag
c      of the box which need refinement is updated to iflag(box)=1
c      and that of the boxes which do not need refinement is
c      updated to iflag(box) = 0
c
c      INPUT arguments
c      curlev         in: integer
c                     the level for which boxes need to be processed
c
c      nboxes         in: integer
c                     total number of boxes
c
c      nlevels        in: integer
c                     total number of levels
c
c      laddr          in: integer(2,0:nlevels)
c                     boxes from laddr(1,ilev) to laddr(2,ilev)
c                     are at level ilev
c
c      nchild         in: integer(nboxes)
c                     nchild(ibox) is the number of children
c                     of box ibox
c
c      ichild         in: integer(4,nboxes)
c                     ichild(j,ibox) is the box id of the jth
c                     child of box ibox
c
c      nnbors         in: integer(nboxes)
c                     nnbors(ibox) is the number of colleagues
c                     of box ibox
c
c      nbors          in: integer(9,nboxes)
c                     nbors(j,ibox) is the jth colleague of box
c                     ibox
c
c      centers        in: double precision(2,nboxes)
c                     x and y coordinates of the box centers
c
c      boxsize        in: double precision(0:nlevels)
c                     boxsize(i) is the size of the box at level i
c
c      iflag          in/out: integer(nboxes)
c                     iflag(ibox)=3 if it is flag++. iflag(ibox) =1
c                     or 0 at the end of routine depending on
c                     whether box needs to be subdivided or not
c
      implicit none
c     Calling sequence variables
      integer curlev, nboxes, nlevels
      integer laddr(2,0:nlevels),nchild(nboxes),ichild(4,nboxes)
      integer nnbors(nboxes), nbors(9,nboxes)
      integer iflag(nboxes)
      double precision centers(2,nboxes),boxsize(0:nlevels)

c     Temporary variables
      integer i,j,k,l,ibox,jbox,kbox,lbox, ict
      double precision distest,xdis,ydis,zdis

      distest = 1.05d0*(boxsize(curlev) + boxsize(curlev+1))/2.0d0
c     Loop over all boxes at the current level     

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,i,jbox,j,kbox,xdis,ydis)
C$OMP$PRIVATE(zdis,ict)
      do ibox = laddr(1,curlev),laddr(2,curlev)
         if(iflag(ibox).eq.3) then
            iflag(ibox) = 0
c           Loop over colleagues of the current box      
            do i=1,nnbors(ibox)
c              Loop over colleagues of flag++ box        
               jbox = nbors(i,ibox)
              
c              Loop over the children of the colleague box
c              Note we do not need to exclude self from
c              the list of colleagues as a self box which
c              is flag++ does not have any children 
c              and will not enter the next loop
               do j=1,4
                  kbox = ichild(j,jbox)
                  if(kbox.gt.0) then
                     if(nchild(kbox).gt.0) then
                        xdis = centers(1,kbox) - centers(1,ibox)
                        ydis = centers(2,kbox) - centers(2,ibox)
                        ict = 0
                        if(abs(xdis).le.distest) ict = ict + 1
                        if(abs(ydis).le.distest) ict = ict + 1
                        if(ict.eq.2) then
                           iflag(ibox) = 1
                           goto 1111
                        endif
                     endif
                  endif
c                 End of looping over the children of the child
c                 of the colleague box
               enddo
c              End of looping over the children of the colleague box       
            enddo
c           End of looping over colleagues            
 1111       continue        
         endif
c        End of testing if the current box needs to checked for         
      enddo
c     End of looping over boxes at the current level      
C$OMP END PARALLEL DO      

      return
      end
c
c
c
c
c
c

      subroutine tree_refine_boxes_flag(iflag,nboxes,
     1  ifirstbox,nbloc,centers,bs,nbctr,nlctr,
     2  ilevel,iparent,nchild,ichild)
      implicit none
      integer nboxes,nbloc,nbctr,nlctr
      real *8 centers(2,nboxes),bs
      integer ilevel(nboxes),iparent(nboxes)
      integer ichild(4,nboxes),nchild(nboxes)
      integer iflag(nboxes)
      integer ifirstbox
      integer, allocatable :: isum(:),itmp(:)

      integer i,ibox,nel0,j,l,jbox,nel1,nbl
      integer ii

      allocate(isum(nbloc),itmp(nbloc))
      do i=1,nbloc
        ibox = ifirstbox+i-1
        itmp(i) = 0
        if(iflag(ibox).gt.0) itmp(i) = 1
      enddo
      if(nbloc.gt.0) call cumsum(nbloc,itmp,isum)
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,nbl,j,jbox)
      do i = 1,nbloc
        ibox = ifirstbox + i-1
        if(iflag(ibox).gt.0) then
          nbl = nbctr + (isum(i)-1)*4
          
          nchild(ibox) = 4
          do j=1,4
            ii = 2
            if(j.le.2) ii = 1
            jbox = nbl+j
            centers(1,jbox) = centers(1,ibox)+(-1)**j*bs/2
            centers(2,jbox) = centers(2,ibox)+(-1)**ii*bs/2 
            iparent(jbox) = ibox
            nchild(jbox) = 0
            do l=1,4
              ichild(l,jbox) = -1
            enddo
            ichild(j,ibox) = jbox
            ilevel(jbox) = nlctr+1 
            if(iflag(ibox).eq.1) iflag(jbox) = 3
            if(iflag(ibox).eq.2) iflag(jbox) = 0
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      

      if(nbloc.gt.0) nbctr = nbctr + isum(nbloc)*4


      return
      end
c
c
c
c
c
c
      subroutine print_tree(itree,ltree,nboxes,centers,boxsize,nlevels,
     1   iptr,fname)
c
c        this subroutine writes the tree info to a file
c
c        input arguments:
c          itree - integer (ltree)
c             packed array containing tree info
c          ltree - integer
c            length of itree
c          nboxes - integer
c             number of boxes
c          centers - real *8 (2,nboxes)
c             xy coordinates of box centers in tree hierarchy
c          boxsize - real *8 (0:nlevels)
c             size of box at various levels
c          nlevels - integer
c             number of levels
c          iptr - integer(8)
c            pointer to various arrays inside itree
c          fname - character *
c            file name to which tree info is to be written
c 
c          output
c            file with name fname, which contains the tree info
c            file can be plotted using the python script
c              tree_plot.py containted in src/common
c


      implicit real *8 (a-h,o-z)
      integer itree(ltree),ltree,nboxes,nlevels,iptr(8)
      real *8 centers(2,nboxes),boxsize(0:nlevels)
      character (len=*) fname

      open(unit=33,file=trim(fname))
      

 1111 format(10(2x,e11.5))      

      do ibox=1,nboxes
         if(itree(iptr(4)+ibox-1).eq.0) then
           ilev = itree(iptr(2)+ibox-1)
           bs = boxsize(ilev)
           x1 = centers(1,ibox) - bs/2
           x2 = centers(1,ibox) + bs/2

           y1 = centers(2,ibox) - bs/2
           y2 = centers(2,ibox) + bs/2
           
           write(33,1111) x1,x2,x2,x1,x1,y1,y1,y2,y2,y1
         endif
      enddo

      close(33)

      return
      end
c
c
c
c
c

      subroutine computemnlists(nlevels,nboxes,itree,ltree,
     1   iptr,centers,
     1   boxsize,iper,mnlist1,mnlist2,mnlist3,mnlist4)
c
c        determine maximum number of elements in list1,list2,list3,list4
c
c        NOTE in 2D: we use max values
c
      implicit real *8 (a-h,o-z)
      integer ltree
      integer nlevels,nboxes,itree(ltree),iptr(8)
      real *8 centers(2,nboxes),boxsize(0:nlevels)
      integer mnlist1,mnlist2,mnlist3,mnlist4

      mnlist1 = 13
      mnlist2 = 27
      mnlist3 = 20
      mnlist4 = 5

      return
      end
c
c
c
c
c

      subroutine computelists(nlevels,nboxes,itree,ltree,
     1            iptr,centers,
     2            boxsize,iper,mnlist1,nlist1,list1,
     3            mnlist2,nlist2,list2,
     4            mnlist3,nlist3,list3,
     5            mnlist4,nlist4,list4)
c
ck
c     This subroutine computes the various fmm lists of a given tree
c     structure
c     
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     itree       in: integer(ltree)
c                   array containing tree info - see start of file
c                   for documentation
c     ltree       in: integer
c                   length of itree array
c 
c     iptr        in: integer(8)
c                   pointer for various arrays in itree
c
c     centers     in: real *8(2,nboxes)
c                 xy coordinates of centers of boxes
c   
c     boxsize     in: real *8(0:nlevels)
c                 Array of boxsizes
c   
c     iper        in: integer
c                 flag for periodic implementations. Currently not used.
c                 Feature under construction
c 
c     mnlist1     in: integer
c                 max number of boxes in list 1 of a box
c
c     mnlist2     in: integer
c                 max number of boxes in list 2 of a box
c
c     mnlist3     in: integer
c                 max number of boxes in list 3 of a box
c
c     mnlist4     in: integer
c                 max number of boxes in list 4 of a box 
c  
c--------------------------------------------------------------
c     OUTPUT arguments:
c     nlist1      out: integer(nboxes)
c                 nlist1(i) is the number of boxes in list 1 
c                 of box i
c
c     list1       out: integer(mnlist1,nboxes)
c                 list1(j,i) is the box id of the jth box in 
c                 list1 of box i
c                      
c     nlist2      out: integer(nboxes)
c                 nlist2(i) is the number of boxes in the list 2
c                 of box i
c 
c     list2       out: integer(mnlist2,nboxes)
c                 list2(j,i) is the box id of the jth box in 
c                 list2 of box i
c
c     nlist3      out: integer(nboxes)
c                 nlist3(i) is the number of boxes in list 3
c                 of box i
c
c     list3       out: integer(mnlist3,nboxes)
c                 list3(j,i) is the box id of the jth box in 
c                 list3 of box i
c                      
c     nlist4      out: integer(nboxes)
c                 nlist4(i) is the number of boxes in the list 2
c                 of box i
c 
c     list4       out: integer(mnlist4,nboxes)
c                 list4(j,i) is the box id of the jth box in 
c                 list4 of box i
c---------------------------------------------------------------
      implicit none
      integer nlevels,nboxes,ltree
      integer iper
      integer itree(ltree),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 centers(2,nboxes)
      integer mnlist1,mnlist2,mnlist3,mnlist4
      integer nlist1(nboxes), list1(mnlist1,nboxes)
      integer nlist2(nboxes), list2(mnlist2,nboxes)
      integer nlist3(nboxes), list3(mnlist3,nboxes)
      integer nlist4(nboxes), list4(mnlist4,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,dad
      integer i,j,ifirstbox,ilastbox,ii
      real *8 distest,xdis,ydis

      do i=1,nboxes
        nlist1(i) = 0
        nlist2(i) = 0
        nlist3(i) = 0
        nlist4(i) = 0
      enddo


c     Setting parameters for level = 0
      if(itree(iptr(4)).eq.0) then
         nlist1(1) = 1 
         list1(1,1) = 1
      else
         nlist1(1) = 0
      endif
      nlist2(1) = 0
      nlist3(1) = 0
      nlist4(1) = 0

      





      do ilev = 1,nlevels
c        Find the first and the last box at level ilev      
         ifirstbox = itree(2*ilev+1)
         ilastbox = itree(2*ilev+2)

         do ibox = ifirstbox,ilastbox
            dad = itree(iptr(3)+ibox-1)

c           Loop over the neighbors of the parent box
c           to find out list 2

            do i=1,itree(iptr(6)+dad-1)

c
c                 jbox is a colleague of ibox's parent
c
               jbox = itree(iptr(7) + 9*(dad-1)+i-1)
               do j=1,itree(iptr(4)+jbox-1)

c
c                  kbox is a child of a colleague of ibox's parent
c
                  kbox = itree(iptr(5)+4*(jbox-1)+j-1)
                   
                  if((abs(centers(1,kbox)-centers(1,ibox)).ge.
     1                1.05*boxsize(ilev)).or.
     2                (abs(centers(2,kbox)-centers(2,ibox)).ge.
     3                1.05*boxsize(ilev))) then
                     
                     nlist2(ibox) = nlist2(ibox)+1
                     list2(nlist2(ibox),ibox) = kbox
                  endif
               enddo
            enddo

cc            goto 1120

c           Compute list 1 and list 3 of ibox if ibox is childless
            if(itree(iptr(4)+ibox-1).eq.0) then
c              Loop over all colleagues of ibox              
               do i=1,itree(iptr(6)+ibox-1)
                  jbox = itree(iptr(7)+9*(ibox-1)+i-1)

c                 If the colleague box is childless, then
c                 colleague box is in list 1
                  if(itree(iptr(4)+jbox-1).eq.0) then
                     nlist1(ibox) = nlist1(ibox)+1
                     list1(nlist1(ibox),ibox) = jbox

c                 If colleague box is not childless, then
c                 test to see if children of colleague
c                 box are in list1 or list 3. 

                  else
                     distest = 1.05d0*(boxsize(ilev) + 
     1                                 boxsize(ilev+1))/2.0d0
c                    Loop over children of colleague box              
                     do j=1,itree(iptr(4)+jbox-1)
                        kbox = itree(iptr(5)+4*(jbox-1)+j-1)
                        xdis = dabs(centers(1,kbox)-centers(1,ibox))
                        ydis = dabs(centers(2,kbox)-centers(2,ibox))
c                       Test to see if child of colleague box
c                       is in list1
                        if(xdis.lt.distest.and.ydis.lt.distest) then
                           nlist1(ibox) = nlist1(ibox)+1
                           list1(nlist1(ibox),ibox)=kbox

                           nlist1(kbox) = nlist1(kbox)+1
                           list1(nlist1(kbox),kbox) = ibox

c                       If it is not in list 1 of ibox then it 
c                       is in list3
                        else
                           nlist3(ibox) = nlist3(ibox)+1
                           list3(nlist3(ibox),ibox)=kbox

                           nlist4(kbox) = nlist4(kbox)+1
                           list4(nlist4(kbox),kbox)=ibox
                        endif
c                       End of figuring out whether child 
c                       of colleague box is in list 1 or list3
                     enddo
c                    End of looping over of children of colleague
c                    box
                  endif
c                 End of checking if colleague box has children
               enddo
c              End of looping over colleague boxes
            endif 
c           End of checking of current box is childless


 1120      continue
         enddo
c        End of looping over boxes at level ilev         
      enddo
c     End of looping over levels      

      return
      end
c
c
c
c
      subroutine compute_modified_list1(nlevels,npwlevel,nboxes,
     1            itree,ltree,iptr,centers,
     2            boxsize,iper,mnlist1,nlist1,list1)
c
c
c     This subroutine computes the modified list1 of a given tree
c     structure
c
c     Note: new definition of list 1 - (1) for ilev < npwlevel,             
c                  list1 is the same as before; (2) for ilev =
c                  npwlevel, boxes at the finer level in list1 
c                  are replaced by their parent, i.e., the colleague 
c                  of ibox; (3) for ilev > npwlevel, list1 has
c                  no use in the FGT.
c                  
c     Assume that the tree is level-restricted. Then the
c     new list1 of source ibox contains all target boxes 
c     that require the evaluation of direct interaction with ibox 
c     in the point FGT.  
c      
c      
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     npwlevel    in: integer
c                 Cutoff level
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     itree       in: integer(ltree)
c                   array containing tree info - see start of file
c                   for documentation
c     ltree       in: integer
c                   length of itree array
c 
c     iptr        in: integer(8)
c                   pointer for various arrays in itree
c
c     centers     in: real *8(2,nboxes)
c                 xy coordinates of centers of boxes
c   
c     boxsize     in: real *8(0:nlevels)
c                 Array of boxsizes
c   
c     iper        in: integer
c                 flag for periodic implementations. Currently not used.
c                 Feature under construction
c 
c     mnlist1     in: integer
c                 max number of boxes in list 1 of a box
c
c--------------------------------------------------------------
c     OUTPUT arguments:
c     nlist1      out: integer(nboxes)
c                 nlist1(i) is the number of boxes in list 1 
c                 of box i
c
c     list1       out: integer(mnlist1,nboxes)
c                 list1(j,i) is the box id of the jth box in 
c                 list1 of box i.
c                  
c      
c---------------------------------------------------------------
      implicit none
      integer nlevels,npwlevel,nboxes,ltree
      integer iper
      integer itree(ltree),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 centers(2,nboxes)
      integer mnlist1,mnlist2,mnlist3,mnlist4
      integer nlist1(nboxes), list1(mnlist1,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,dad
      integer i,j,ifirstbox,ilastbox,ii
      real *8 distest,xdis,ydis,bs,xdp1,xdm1,ydp1,ydm1

      bs=boxsize(0)
      
      do i=1,nboxes
        nlist1(i) = 0
      enddo


c     Setting parameters for level = 0
      if(itree(iptr(4)).eq.0) then
         nlist1(1) = 1 
         list1(1,1) = 1
      else
         nlist1(1) = 0
      endif

      
      do ilev = 1,nlevels
c        Find the first and the last box at level ilev      
         ifirstbox = itree(2*ilev+1)
         ilastbox = itree(2*ilev+2)

         do ibox = ifirstbox,ilastbox
            dad = itree(iptr(3)+ibox-1)

c           Compute list 1 and list 3 of ibox if ibox is childless
            if(itree(iptr(4)+ibox-1).eq.0) then
c              Loop over all colleagues of ibox              
               do i=1,itree(iptr(6)+ibox-1)
                  jbox = itree(iptr(7)+9*(ibox-1)+i-1)

                  if (ilev .ge. npwlevel) then
                     if (itree(iptr(4)+jbox-1).eq.0) then
                        nlist1(ibox) = nlist1(ibox)+1
                        list1(nlist1(ibox),ibox) = jbox
                     endif
                  else
c                    If the colleague box is childless, then
c                    colleague box is in list 1
                     if(itree(iptr(4)+jbox-1).eq.0) then
                        nlist1(ibox) = nlist1(ibox)+1
                        list1(nlist1(ibox),ibox) = jbox

c                    If colleague box is not childless, then
c                    test to see if children of colleague
c                    box are in list1 or list 3. 

                     else
                        distest = 1.05d0*(boxsize(ilev) + 
     1                                 boxsize(ilev+1))/2.0d0
c                       Loop over children of colleague box              
                        do j=1,itree(iptr(4)+jbox-1)
                           kbox = itree(iptr(5)+4*(jbox-1)+j-1)
                           xdis = dabs(centers(1,kbox)-centers(1,ibox))
                           ydis = dabs(centers(2,kbox)-centers(2,ibox))
                           
                           if (iper .eq. 1) then
                              xdp1 = bs-xdis
                              if (xdp1.lt.xdis) xdis=xdp1
                              ydp1 = bs-ydis
                              if (ydp1.lt.ydis) ydis=ydp1
                           endif
c                       Test to see if child of colleague box
c                       is in list1
                           if(xdis.lt.distest.and.ydis.lt.distest) then
                              nlist1(ibox) = nlist1(ibox)+1
                              list1(nlist1(ibox),ibox)=kbox

                              nlist1(kbox) = nlist1(kbox)+1
                              list1(nlist1(kbox),kbox) = ibox
                           endif
c                       End of figuring out whether child 
c                       of colleague box is in list 1 or list3
                        enddo
c                    End of looping over of children of colleague
c                    box
                     endif
c                    End of checking if colleague box has children
                  endif
               enddo
c              End of looping over colleague boxes
            endif 
c           End of checking of current box is childless


 1120      continue
         enddo
c        End of looping over boxes at level ilev         
      enddo
c     End of looping over levels      

      return
      end
c
c
c
c
      subroutine gt2d_computemnlistpw(nlevels,nboxes,itree,ltree,
     1   iptr,centers,
     1   boxsize,iper,mnlistpw)
c
c     determine maximum number of elements in listsoe and listsx
c
c     NOTE in 2D: we use max values
c
      implicit real *8 (a-h,o-z)
      integer ltree
      integer nlevels,nboxes,itree(ltree),iptr(8)
      real *8 centers(2,nboxes),boxsize(0:nlevels)
      integer mnlistpw

      mnlistpw = 8

      return
      end
c
c
c
c
c
      subroutine gt2d_computelistpw(nlevels,npwlevel,nboxes,
     1    itree,ltree,iptr,centers,boxsize,laddr,
     2    mnlistpw,nlistpw,listpw)
c     this subroutine returns lists of pw interaction boxes
c
c     input parameters:
c
c     nlevels = number of levels, level 0 contains boxes of size 6\sqrt{delta}
c     npwlevel   = recursive SOE picture stops at the level npwlevel
c     nboxes = total number of boxes in the tree
c     itree - laddr: tree stuff
c
c     mnlistpw = maximum number of a particular type of PW expansions
c      
c     output parameters:
c
c     nlistpw (nboxes) = contains the number of boxes for the PW interaction
c                            for each box
c     listpw (mnlistpw,nboxes) = contains the box ID for the PW interaction
c                            for each box
c      
      implicit real *8 (a-h,o-z)
      integer nlevels,nboxes,ltree
      integer itree(ltree),iptr(8),laddr(2,0:nlevels)
      real *8 centers(2,nboxes),boxsize(0:nlevels)
      integer mnlistpw
      integer nlistpw(nboxes), listpw(mnlistpw,nboxes)

      integer ibox

      do ibox=1,nboxes
         nlistpw(ibox)=0
      enddo

      nlevstart=max(npwlevel,0)
      do ilev=npwlevel,npwlevel
        do ibox = laddr(1,ilev),laddr(2,ilev)
           ncoll = itree(iptr(6)+ibox-1)
           do i=1,ncoll
              jbox = itree(iptr(7) + (ibox-1)*9+i-1)
              jlev = itree(iptr(2)+jbox-1)
              if (jbox.ne.ibox .and. ilev.eq.jlev .and. 
     1            (itree(iptr(4)+jbox-1).gt.0
     2            .or. itree(iptr(4)+ibox-1).gt.0)) then
                 nlistpw(ibox)=nlistpw(ibox)+1
                 listpw(nlistpw(ibox),ibox) = jbox
              endif
           enddo
        enddo
c     end of ilev do loop
      enddo

      return
      end
c
c
c
c
c
      subroutine oldtree2newtree(nlev,levelbox,iparentbox,
     2    ichildbox,icolbox,irowbox,nboxes,nblevel,
     3    iboxlev,istartlev,cent0,xsize0,iperiod,
     4    ltree,nlevels,itree,iptr,centers,boxsize)
      implicit real *8 (a-h,o-z)
c
c     convert an old tree to the new tree
c
c     input parameters:
c     nlev - total number of  levels
c     levelbox - an array determining the level of each box
c     iparentbox - the parent of each box
c     ichildbox - the four children of each box
c     icolbox - the column of each box
c     irowbox - the row of each box
c     nboxes - integer
c          number of boxes
c     nblevel - the total number of boxes per level
c     iboxlev - the array in which the boxes are arranged
c     istartlev - the pointer to where each level
c               begins in the iboxlev array
c     cent0 - center of the root box
c     xsize0 - size of the root box
c     iperiod = 0 : free space
c               1: periodic
c     ltree - integer
c          length of tree = 2*(nlevels+1)+17*nboxes
c
c     output:
c     nlevels - integer
c          number of levels
c     itree - integer(ltree)
c          tree info
c     iptr - integer(8)
c          iptr(1) - laddr
c          iptr(2) - ilevel
c          iptr(3) - iparent
c          iptr(4) - nchild
c          iptr(5) - ichild
c          iptr(6) - ncoll
c          iptr(7) - coll
c          iptr(8) - ltree
c     centers - double precision (2,nboxes)
c          xy coordinates of box centers in the oct tree
c     boxsize - double precision (0:nlevels)
c          size of box at each of the levels
      integer *4  levelbox(1)
      integer *4  nlev, nboxes
      integer *4  icolbox(1), irowbox(1)
      integer *4  iparentbox(1), ichildbox(4,1)
      integer *4  nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer nlevels,ltree
      integer iptr(8)
      integer itree(ltree)
      real *8 cent0(2),xsize0
      real *8 centers(2,*),boxsize(0:1)
      integer iboxlevinv(nboxes)
      integer, allocatable :: icolleagbox(:,:)

      allocate(icolleagbox(9,nboxes))
      
      nlevels = nlev
      
      iptr(1) = 1
      iptr(2) = 2*(nlevels+1)+1
      iptr(3) = iptr(2) + nboxes
      iptr(4) = iptr(3) + nboxes
      iptr(5) = iptr(4) + nboxes
      iptr(6) = iptr(5) + 4*nboxes
      iptr(7) = iptr(6) + nboxes
      iptr(8) = iptr(7) + 9*nboxes



      boxsize(0) = xsize0

      centers(1,1) = cent0(1)
      centers(2,1) = cent0(2)
      
cccc      call prinf('iboxlev=*',iboxlev,nboxes)
cccc      call prinf('istartlev=*',istartlev,nlev+1)
      
      do i=1,nlevels
         boxsize(i)=boxsize(i-1)/2
      enddo

      do i=1,nboxes
         ibox=iboxlev(i)
         iboxlevinv(ibox)=i
      enddo
      
      do ilev=0,nlevels
         do i=istartlev(ilev),istartlev(ilev)+nblevel(ilev)-1
            ibox=iboxlev(i)
            icol=icolbox(ibox)
            irow=irowbox(ibox)
         
            bs = boxsize(ilev)
         
            centers(1,i) = cent0(1)-xsize0/2 + (icol-0.5d0) * bs
            centers(2,i) = cent0(2)-xsize0/2 + (irow-0.5d0) * bs
         enddo
      enddo

      call mkcolls(icolbox,
     1      irowbox, icolleagbox, nboxes, nlev,
     2      iparentbox, ichildbox, nblevel,
     3      iboxlev, istartlev, iperiod)

      do ilev=0,nlevels
         itree(iptr(1)+2*ilev)=istartlev(ilev)
         itree(iptr(1)+2*ilev+1)=istartlev(ilev)+nblevel(ilev)-1
      enddo

      do ilev=0,nlevels
        do i=istartlev(ilev),istartlev(ilev)+nblevel(ilev)-1
          ibox=iboxlev(i)
          itree(iptr(2)+i-1) = levelbox(ibox)
cccc          print *, ilev, levelbox(ibox)
          
          jbox=iparentbox(ibox)
          if (jbox.gt.0) then
             itree(iptr(3)+i-1) = iboxlevinv(jbox)
          else
             itree(iptr(3)+i-1) = -1
          endif

          do j=1,4
             itree(iptr(5)+4*(i-1)+j-1)=-1
          enddo
          
          ichild=0
          do j=1,4
             jbox=ichildbox(j,ibox)
             if (jbox.gt.0) then
                itree(iptr(5)+4*(i-1)+ichild)=iboxlevinv(jbox)
                ichild=ichild+1
             endif 
          enddo
          itree(iptr(4)+i-1) = ichild

          do j=1,9
             itree(iptr(7)+9*(i-1)+j-1)=-1
          enddo
          
          icoll=0
          do j=1,9
             jbox=icolleagbox(j,ibox)
             if (jbox.gt.0) then
                itree(iptr(7)+9*(i-1)+icoll)=iboxlevinv(jbox)
                icoll=icoll+1
             endif
          enddo
          itree(iptr(6)+i-1)=icoll
        enddo
      enddo

cccc      call prinf('itree=*',itree,ltree)
      
      return
      end subroutine
c
c
c
c
      subroutine newtree2oldtree(nboxes,nlevels,ltree,itree,iptr,
     1    centers,boxsize,nlev,levelbox,iparentbox,
     2    ichildbox,icolbox,irowbox,nblevel,
     3    iboxlev,istartlev)
      implicit real *8 (a-h,o-z)
c
c     convert a new tree to the old tree used by Frank's box code
c
c     input parameters:
c     nboxes - integer
c          number of boxes
c     nlevels - integer
c          number of levels
c     itree - integer(ltree)
c          tree info
c     iptr - integer(8)
c          iptr(1) - laddr
c          iptr(2) - ilevel
c          iptr(3) - iparent
c          iptr(4) - nchild
c          iptr(5) - ichild
c          iptr(6) - ncoll
c          iptr(7) - coll
c          iptr(8) - ltree
c     centers - double precision (2,nboxes)
c          xy coordinates of box centers in the oct tree
c     boxsize - double precision (0:nlevels)
c          size of box at each of the levels
c     ltree - integer
c          length of tree = 2*(nlevels+1)+17*nboxes
c
c     output:
c     nlev - total number of  levels
c     levelbox - an array determining the level of each box
c     iparentbox - the parent of each box
c     ichildbox - the four children of each box
c     icolbox - the column of each box
c     irowbox - the row of each box
c     nblevel - the total number of boxes per level
c     iboxlev - the array in which the boxes are arranged
c     istartlev - the pointer to where each level
c               begins in the iboxlev array
c     cent0 - center of the root box
c     xsize0 - size of the root box
c     iperiod = 0 : free space
c               1: periodic

      integer nlevels,ltree
      integer iptr(8)
      integer itree(ltree)
      real *8 centers(2,*),boxsize(0:1)

      integer *4  levelbox(1)
      integer *4  nlev, nboxes
      integer *4  icolbox(1), irowbox(1)
      integer *4  iparentbox(1), ichildbox(4,1)
      integer *4  nblevel(0:1), iboxlev(1), istartlev(0:1)
      real *8 cent0(2),xsize0

      integer ichild(4),icolj(4),irowj(4)
      
      nlev = nlevels
      
      do ilev=0,nlevels
         istartlev(ilev) = itree(2*ilev+1)
         nblevel(ilev) = itree(2*ilev+2)-itree(2*ilev+1)+1
      enddo

c      call prinf('istartlev=*',istartlev(0),nlevels+1)
c      call prinf('nblevel=*',nblevel(0),nlevels+1)

      xsize0 = boxsize(0)
     
      cent0(1) = centers(1,1)
      cent0(2) = centers(2,1)

c      call prin2('xsize=*',xsize0,1)
c      call prin2('centers=*',centers,2)
      
      ibox=0
      do ilev=0,nlevels
         do i=istartlev(ilev),istartlev(ilev)+nblevel(ilev)-1
            ibox=ibox+1
c            if (ibox.ne. i) print *, ibox, i
            
            iboxlev(ibox)=ibox
            levelbox(ibox)=ilev
            
            bs = boxsize(ilev)

            icol = (centers(1,ibox)-cent0(1)+xsize0/2)/bs+0.55d0
            irow = (centers(2,ibox)-cent0(2)+xsize0/2)/bs+0.55d0
c            print *, ibox,centers(1,ibox),bs,icol
c            print *, ibox,centers(2,ibox),bs,irow
            
            icolbox(ibox)=icol
            irowbox(ibox)=irow
         enddo
      enddo

c      call prinf('iboxlev=*',iboxlev,nboxes)
c      call prinf('icolbox=*',icolbox,nboxes)
c      call prinf('irowbox=*',irowbox,nboxes)

      ibox=0
      do ilev=0,nlevels
        do i=istartlev(ilev),istartlev(ilev)+nblevel(ilev)-1
          ibox=ibox+1

          iparentbox(ibox)=itree(iptr(3)+i-1)
          
          do j=1,4
             ichild(j)=itree(iptr(5)+4*(i-1)+j-1)
          enddo
c     reorder the children in the clockwise order, starting from the topleft one
          if (ichild(1).gt.0) then
             icolj(1)=icolbox(ichild(1))
             irowj(1)=irowbox(ichild(1))
             icolmin=icolj(1)
             irowmin=irowj(1)
             do j=2,4
                icolj(j)=icolbox(ichild(j))
                irowj(j)=irowbox(ichild(j))
                if (icolj(j).lt.icolmin) icolmin=icolj(j)
                if (irowj(j).lt.irowmin) irowmin=irowj(j)
             enddo
             
             do j=1,4
                if (icolj(j).eq.icolmin.and.irowj(j).eq.irowmin+1) then
                   ichildbox(1,ibox)=ichild(j)
                elseif (icolj(j).eq.icolmin+1 .and.
     1                 irowj(j).eq.irowmin+1) then
                   ichildbox(2,ibox)=ichild(j)
                elseif (icolj(j).eq.icolmin+1 .and.
     1                 irowj(j).eq.irowmin) then
                   ichildbox(3,ibox)=ichild(j)
                elseif (icolj(j).eq.icolmin .and.
     1                 irowj(j).eq.irowmin) then
                   ichildbox(4,ibox)=ichild(j)
                endif
             enddo
          else
             do j=1,4
                ichildbox(j,ibox)=-1
             enddo
          endif
        enddo
      enddo

c      call prinf('iparentbox=*',iparentbox,nboxes)
c      call prinf('ichildbox=*',ichildbox,4*nboxes)

      return
      end subroutine
      
      
