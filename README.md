# boxcode-legacy
This repo contains a basic Fortran and Matlab wrapper for running the 
8th order box code by Greengard and Ethridge.

To verify successful compilation of the library, run `make test`
after copying over the appropriate version of make.inc

To generate the mex file run `make matlab` after copying over the 
appropriate version of make.inc

The interface on input requires a level restricted adaptive tree, 
and function values tabulated on a 8th order tensor product 
chebyshev grid at all the leaf boxes.

For more detailed documentation for the input format of the level
restricted tree checkout either the subroutine lbfmm2d in src/lbfmm2d.f 
or the matlab file `volfmm8.m` or `volfmm8_legacy.m`

There are two matlab interfaces, one for the legacy tree format, and
one for the new tree format

In fortran on a single core, the code runs at approximately 500 k points 
per second per core for 12 digits of accuracy. 

The legacy version of the fortran code, poisson8.f contains more
complicated options for evalauting the volume potential such as
imposing periodic boundary conditions, or imposing dirichlet/neumann
boundary conditions on the edge of the square. However, the more
complicated routines have not been thoroughly tested in this wrapped
environment and should be used with caution.
