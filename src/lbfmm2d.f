
      subroutine lbfmm2d(nd,eps,iperiod,nboxes,nlevels,ltree,itree,iptr,
     1    norder,npols,ttype,fvals,centers,boxsize,npbox,
     2    pot,timeinfo)
c
c     This code computes the solution to the Poisson equation in a box
c     This is wrapper to Frank's old code with several caveats.
c     1. nd is always 1 for now; 2. norder can only be 8; 3. iperiod 
c     is 0 - other boundary conditions have not been translated;
c     4. the internal flag for eps is iprec that takes values 0-3
c        (3,6,9,12 digits of accuracy);
c     5. boxsize(0) and centers(:,1) are fixed.
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
c         npols - integer
c           number of coefficients of expansions of functions
c           in each of the boxes
c         ttype - character *1
c            type of coefs provided, total order ('t') or full order('f')
c         fvals - double precision (npbox,nboxes)
c           function tabulated on the Chebyshev tensor grid
c         centers - double precision (2,nboxes)
c           xy coordintes of boxes in the tree structure
c         boxsize - double precision (0:nlevels)
c           size of boxes at each of the levels
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
      character *1 ttype
      real *8 fvals(npbox,nboxes)
      real *8 pot(npbox,nboxes)
      real *8 centers(2,nboxes)
      real *8 boxsize(0:nlevels)
      real *8 timeinfo(6)

      integer, allocatable :: levelbox(:),iparentbox(:),ichildbox(:,:)
      integer, allocatable :: icolbox(:),irowbox(:),nblevel(:)
      integer, allocatable :: iboxlev(:),istartlev(:)
      
      allocate(levelbox(nboxes),iparentbox(nboxes),ichildbox(4,nboxes))
      allocate(icolbox(nboxes),irowbox(nboxes),nblevel(0:nlevels))
      allocate(iboxlev(nboxes),istartlev(0:nlevels))

      call newtree2oldtree(nboxes,nlevels,ltree,itree,iptr,
     1    centers,boxsize,nlev,levelbox,iparentbox,
     2    ichildbox,icolbox,irowbox,nblevel,
     3    iboxlev,istartlev)

      if (eps.le.1d-12) then
         iprec=3
      elseif (eps.le.1d-9) then
         iprec=2
      elseif (eps.le.1d-6) then
         iprec=1
      else
         iprec=0
      endif

      ifnear = 1
      
      call fmmstart8_wrap(nlev,levelbox,iparentbox,ichildbox,
     1    icolbox,irowbox,nboxes,nblevel,iboxlev,istartlev,
     2    fvals,ifnear,pot,iprec)
      
      return
      end
c
c
c
