      subroutine chebtens_exps_2d(itype,n,type,x,u,ldu,v,ldv,w)
c                 input parameters:
c
c  itype - the type of the calculation to be performed
c          itype=0 means that only the gaussian nodes are 
c                  to be constructed. 
c          itype=1 means that only the nodes and the weights 
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c          itype=3 only construct u
c          itype=4 only construct v
      
      implicit none
      integer itype, n, ldu, ldv
      character type
      real *8 x(2,*),w(*)
      real *8 u(ldu,*), v(ldv,*)
      real *8 x1d(n), w1d(n), u1d(n,n), v1d(n,n)
      integer i,j,ipt,itype1d,io, jo,ipol
      
      itype1d = 0
      if (itype .ge. 1) then
         itype1d = 1
      endif
      if (itype .ge. 2) then
         itype1d = 2
      endif
      
      call chebexps(itype1d,n,x1d,u1d,v1d,w1d)

      ipt = 0
      do i=1,n
         do j=1,n
            ipt = ipt + 1
            x(1,ipt) = x1d(j)
            x(2,ipt) = x1d(i)               
        enddo
      enddo

      if (itype .ge. 1) then
         ipt = 0
         do i=1,n
            do j=1,n
               ipt = ipt + 1
               w(ipt) = w1d(i)*w1d(j)
            enddo
         enddo
      endif


      if (itype .eq. 2 .or. itype .eq. 3) then
c     construct u from 1d u
         if (type .eq. 'f' .or. type .eq. 'F') then         
            ipt = 0
            do io = 1,n
            do jo = 1,n
               ipt = ipt + 1
               ipol = 0
               do i=1,n
               do j=1,n
                  ipol = ipol + 1
                  u(ipol,ipt) =  u1d(i,io)*u1d(j,jo)
               enddo
               enddo
            enddo
            enddo

         else if (type .eq. 't' .or. type .eq. 'T') then

            ipt = 0
            do io = 1,n
            do jo = 1,n
               ipt = ipt + 1
               ipol = 0
               do i=1,n
               do j=1,n+1-i
                  ipol = ipol + 1
                  u(ipol,ipt) = u1d(i,io)*u1d(j,jo)
               enddo
               enddo
            enddo
            enddo
         endif
      endif
      
      if (itype .eq. 2 .or. itype .eq. 4) then
c     construct v from 1d v
         if (type .eq. 'f' .or. type .eq. 'F') then         
            ipol = 0
            do io = 1,n
            do jo = 1,n
               ipol = ipol + 1
               ipt = 0
               do i=1,n
               do j=1,n
                  ipt = ipt + 1
                  v(ipt,ipol) =  v1d(i,io)*v1d(j,jo)
               enddo
               enddo
            enddo
            enddo

         else if (type .eq. 't' .or. type .eq. 'T') then

            ipol = 0
            do io = 1,n
            do jo = 1,n+1-io
               ipol = ipol + 1
               ipt = 0
               do i=1,n
               do j=1,n
                  ipt = ipt + 1
                  v(ipt,ipol) = v1d(i,io)*v1d(j,jo)
               enddo
               enddo
            enddo
            enddo
         endif
      endif         

      return
      end
c
c
