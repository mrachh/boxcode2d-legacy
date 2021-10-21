gfortran -c tables8.f
gfortran -std=legacy -c chebrouts.f
gfortran -std=legacy -O3 -c poisson8.f
/Applications/MATLAB_R2021a.app/bin/mex poisson8_matlab.c chebrouts.o poisson8.o tables8.o \
    -compatibleArrayDims -DMWF77_UNDERSCORE1 -L/usr/local/lib/gcc/11 -output poisson8_matlab \
    -lm -lstdc++ -ldl -lgfortran

