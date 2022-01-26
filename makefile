
EXEC = int2

#HOST = osx
HOST=linux-gfortran
#HOST=linux-ifort
#HOST=linux-gfortran-prof
#HOST=linux-gfortran-openmp

ifeq ($(HOST),osx)
FC = gfortran
FFLAGS = -O3 -march=native --openmp -funroll-loops -c -w
FLINK = gfortran -w --openmp -o $(EXEC)
FEND = -lopenblas ${LDFLAGS}
endif

ifeq ($(HOST),linux-gfortran)
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -funroll-loops -ftree-vectorize -ffast-math -std=legacy -c -w
#FFLAGS = -fPIC -O3 -march=native -funroll-loops -std=legacy -fcx-limited-range -c -w
FLINK = gfortran -w -o $(EXEC)
FEND = -lblas -llapack
#FEND = -lopenblas -L/usr/local/opt/openblas/lib 
endif

ifeq ($(HOST),linux-gfortran-prof)
FC = gfortran
FFLAGS = -O3 -march=native -pg -g -funroll-loops -ftree-vectorize -ffast-math -c -w  
FLINK = gfortran -w -o $(EXEC) -pg
FEND = -lblas -llapack
endif

ifeq ($(HOST),linux-gfortran-openmp)
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -c --openmp
FLINK = gfortran -w --openmp -o $(EXEC) 
FEND = -lopenblas -L/usr/local/opt/openblas/lib 
endif

ifeq ($(HOST),linux-ifort)
FC = ifort
FFLAGS = -O1 -g -c -w -xW -qopenmp 
FFLAGS = -fast -f77rtl -c
FLINK = ifort -w -qopenmp -o $(EXEC)
FLINK = ifort -w -o $(EXEC)
WITH_SECOND = 0
endif


SRC = .
FMM3D = ../src
COMMON = ./common


.PHONY: all clean list

SOURCES = lbfmm2d.f \
	chebrouts.f \
	tables8.f \
	tree_routs4.f \
	tree_vol_coeffs_cheb.f \
	$(COMMON)/prini_new.f \
	$(COMMON)/legeexps.f \
	$(COMMON)/chebexps.f \
	$(COMMON)/legetens.f \
	$(COMMON)/chebtens.f \
	$(COMMON)/voltab2d.f \
	test_lbfmm2d.f	
#	test_lbfmm2d_old.f 

ifeq ($(WITH_SECOND),1)
SOURCES += $(SRC)/second-r8.f
endif

OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCES)))

#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : %.f
	$(FC) $(FFLAGS) $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@

%.mod : %.f90
	$(FC) $(FFLAGS) $< 

all: $(OBJECTS)
	rm -f $(EXEC)
	$(FLINK) $(OBJECTS) $(FEND)
	./$(EXEC) 2 

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)
	rm -f fort*
	rm -f int*
	rm -f *~

list: $(SOURCES)
	$(warning Requires:  $^)



