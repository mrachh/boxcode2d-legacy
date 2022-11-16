function [pot] = volfmm(t,f,eps)
%  This subroutine calls frank's boxcode to evaluate
%  the free space convolution against the laplace green's
%  function
%
%  input arguments:
%    - t: tree structure
%        t.nlev: int
%           number of levels in the tree, root level is 0
%        t.nboxes: int 
%           number of boxes in the tree
%        t.levelbox: (nboxes)
%           levelbox(i) indicates the level of box i
%        t.iparentbox: (nboxes)
%           iparentbox(i) is the parent of box i
%        t.ichildbox: (4,nboxes)
%           ichildbox(1:4,i) are the children of box i,
%           child 1 is top left corner, child 2 is top right corner,
%           child 3 is bottom right corner, child 4 is bottom left corner
%        t.icolbox: (nboxes)
%           column number of box i, assuming a uniform grid of boxes 
%           at the level of box i (1 indexed)
%        t.irowbox: (nboxes)
%           row number of box i, assuming a uniform grid of boxes 
%           at the level of box i (1 indexed)
%        t.nblevel: (nlev+1)
%           number of boxes at each level
%        t.iboxlev: (nboxes)
%           possible permutation of boxes at a fixed level
%        t.istartlev (nlev+1)
%           location in the sorted boxes array where boxes at level
%           i begin
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
%
%           
%
%
    iprec = 1;
    if(eps >= 0.5e-3)
        iprec = 1;
    else if(eps >= 0.5e-6)
        iprec = 2;
    else if(eps >= 0.5e-9)
        iprec = 3;
    else
        iprec = 4;
    end
    nlev = t.nlev;
    levelbox = t.levelbox;
    iparentbox = t.iparentbox;
    ichildbox = t.ichildbox;
    icolbox = t.icolbox;
    irowbox = t.irowbox;
    nboxes = t.nboxes;
    nblevel = t.nblevel;
    iboxlev = t.iboxlev;
    istartlev = t.istartlev;
    nlevp1 = nlev + 1;
    ifnear = 1;
    nd = 64;

    mex_id_ = 'fmmstart8_wrap(i int[x], i int[x], i int[x], i int[xx], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i int[x], io double[xx], i int[x])';
[pot] = poisson8_matlab(mex_id_, nlev, levelbox, iparentbox, ichildbox, icolbox, irowbox, nboxes, nblevel, iboxlev, istartlev, f, ifnear, pot, iprec, 1, nboxes, nboxes, 4, nboxes, nboxes, nboxes, 1, nlevp1, nboxes, nlevp1, nd, nboxes, 1, nd, nboxes, 1);

end
