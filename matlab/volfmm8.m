function [pot] = volfmm8(t,f,eps,opts)
%  This subroutine calls frank's boxcode to evaluate
%  the free space convolution against the laplace green's
%  function but using the new tree interface
%
%  input arguments:
%    - t: tree structure
%        t.nlev: int
%           number of levels in the tree, root level is 0
%        t.nboxes: int 
%           number of boxes in the tree
%        t.itree: int(ltree)
%           array containing tree structure
%        t.iptr: int(8)
%           iptr - integer(8)
%              pointer to various parts of the tree structure
%              iptr(1) - laddr
%              iptr(2) - ilevel
%              iptr(3) - iparent
%              iptr(4) - nchild
%              iptr(5) - ichild
%              iptr(6) - ncoll
%              iptr(7) - coll
%              iptr(8) - ltree
%        t.centers: double(2,nboxes)
%           location of center of boxes in tree structure
%        t.boxsize: double(nlevels+1)
%           boxsize of various levels strating from level0
%    - f: (64,nboxes)
%         right hand side: ordering of chebyshev nodes given by
%         do j=1,8
%           do l=1,8
%              xf(8*(l-1)+j) = x(j) + xshift
%              yf(8*(l-1)+j) = x(j) + yshift
%           enddo
%         enddo
%
%         where x(j) are scaled chebyshev nodes on the box
%    - opts: options structure (optional)
%         opts.ifnear: flag for determining whether to include list1
%           interactions or not, default value, ifnear = 1, 
%           which includes the near interaction
%
%           
%
%
    if(nargin == 3)
        opts = [];
    end
    nlev = t.nlev;
    itree = t.itree(:);
    ltree = length(itree);
    iptr = t.iptr;
    centers = t.centers;
    boxsize = t.boxsize;

    nboxes = t.nboxes;
    nlevp1 = nlev + 1;


    ifnear = 1;
    if(isfield(opts,'ifnear'))
       ifnear = opts.ifnear;
    end
    npbox = 64;
    nd = 1;
    iperiod = 0;
    norder = 8;
    timeinfo = zeros(6,1);


    mex_id_ = 'lbfmm2d(i int[x], i double[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i double[x], i int[x], io double[xx], io double[x])';
[pot, timeinfo] = poisson8_matlab(mex_id_, nd, eps, nboxes, nlev, ltree, itree, iptr, norder, npbox, fvals, centers, boxsize, ifnear, pot, timeinfo, 1, 1, 1, 1, 1, ltree, 8, 1, 1, npbox, nboxes, 2, nboxes, nlevp1, 1, npbox, nboxes, 6);

end

%
%
%
%
%
%


