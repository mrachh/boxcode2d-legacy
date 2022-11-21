
      subroutine lbfmm2d(nd,eps,iperiod,nboxes,nlevels,ltree,itree,iptr,
     1    norder,npbox,fvals,centers,boxsize,ifnear,pot,timeinfo)
c
c     This code computes the solution to the Poisson equation in a box
c     This is wrapper to Frank's old code with several caveats.
c     1. norder can only be 8
c     2. iperiod is 0 - other boundary conditions have not been translated;
c     3. the internal flag for eps is iprec that takes values 0-3
c        (3,6,9,12 digits of accuracy);
c     
c     
c       input
c         eps - double precision
c            tolerance requested
c         nboxes - integer
c            number of boxes
c         nlevels - integer
c            number of levels
c         ltree - integer
c            length of array containing the tree structure
c         itree - integer(ltree)
c            array containing the tree structure
c         iptr - integer(8)
c            pointer to various parts of the tree structure
c           iptr(1) - laddr
c           iptr(2) - ilevel
c           iptr(3) - iparent
c           iptr(4) - nchild
c           iptr(5) - ichild
c           iptr(6) - ncoll
c           iptr(7) - coll
c           iptr(8) - ltree
c         norder - integer
c           order of expansions for input coefficients array
c         fvals - double precision (npbox,nboxes)
c           function tabulated on the Chebyshev tensor grid
c         centers - double precision (2,nboxes)
c           xy coordintes of boxes in the tree structure
c         boxsize - double precision (0:nlevels)
c           size of boxes at each of the levels
c         ifnear - integer
c            flag for deciding whether to include list1 interactions
c         npbox - integer
c           number of points per box where potential is to be dumped = (norder**2)
c
c     output:
c         pot - double precision (npbox,nboxes)
c            volume potential on the tree structure (note that 
c            the potential is non-zero only in the leaf boxes)
c
      implicit real *8 (a-h,o-z)
      real *8 eps
      integer nboxes,nlevels,ltree
      integer itree(ltree),iptr(8),norder,ncbox,npbox
      real *8 fvals(nd,npbox,nboxes)
      real *8 pot(nd,npbox,nboxes)
      real *8 centers(2,nboxes)
      real *8 boxsize(0:nlevels)
      real *8 timeinfo(6)

      real *8, allocatable :: fvalstmp(:,:),pottmp(:,:)
      real *8, allocatable :: centers_use(:,:),boxsize_use(:)

      integer, allocatable :: levelbox(:),iparentbox(:),ichildbox(:,:)
      integer, allocatable :: icolbox(:),irowbox(:),nblevel(:)
      integer, allocatable :: iboxlev(:),istartlev(:)
      
      allocate(levelbox(nboxes),iparentbox(nboxes),ichildbox(4,nboxes))
      allocate(icolbox(nboxes),irowbox(nboxes),nblevel(0:nlevels))
      allocate(iboxlev(nboxes),istartlev(0:nlevels))

      allocate(centers_use(2,nboxes),boxsize_use(0:nlevels))
c
c
c  rescale and recenter boxes
c
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox)
      do ibox=1,nboxes
        centers_use(1,ibox) = (centers(1,ibox)-centers(1,1))/boxsize(0)
        centers_use(2,ibox) = (centers(2,ibox)-centers(2,1))/boxsize(0)
      enddo
C$OMP END PARALLEL DO

      do ilev=0,nlevels
        boxsize_use(ilev) = boxsize(ilev)/boxsize(0)
      enddo

      call newtree2oldtree(nboxes,nlevels,ltree,itree,iptr,
     1    centers_use,boxsize_use,nlev,levelbox,iparentbox,
     2    ichildbox,icolbox,irowbox,nblevel,
     3    iboxlev,istartlev)

      allocate(fvalstmp(npbox,nboxes),pottmp(npbox,nboxes))

      if (eps.le.1d-12) then
         iprec=3
      elseif (eps.le.1d-9) then
         iprec=2
      elseif (eps.le.1d-6) then
         iprec=1
      else
         iprec=0
      endif


      do idim=1,nd

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,j)      
        do ibox=1,nboxes
          do j=1,npbox
            fvalstmp(j,ibox) = fvals(idim,j,ibox)
          enddo
        enddo
C$OMP END PARALLEL DO        

        call fmmstart8_wrap(nlev,levelbox,iparentbox,ichildbox,
     1      icolbox,irowbox,nboxes,nblevel,iboxlev,istartlev,
     2      fvalstmp,ifnear,pottmp,iprec)
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,j)      
        do ibox=1,nboxes
          do j=1,npbox
            pot(idim,j,ibox) = pottmp(j,ibox)
          enddo
        enddo
C$OMP END PARALLEL DO        
      enddo
      
      return
      end
c
c
c
