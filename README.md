# boxcode-legacy
This repo contains a basic Matlab wrapper for running the 8th order box code by Greengard and Ethridge.

To generate the mex file on a mac, simply run the bash script ./runsh.sh

The interface on input requires a level restricted adaptive tree, and function values tabulated on a 8th order tensor product chebyshev grid at all the leaf boxes.

In fortran on a single core, the code runs at approximately 500 k points per second per core.
