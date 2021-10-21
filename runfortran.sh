gfortran -c tables8.f
gfortran -std=legacy -c chebrouts.f
gfortran -std=legacy -O3 -c poisson8.f
gfortran -o int2 tables8.o chebrouts.o poisson8.o
./int2


