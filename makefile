
EXEC = int2

# Set compiler
FC=gfortran

# Set compile flags
FFLAGS=-fPIC -O3 -march=native -funroll-loops -std=legacy -w
FFLAGS2=-fPIC -std=legacy -w

# extra flags for multithreading
OMPFLAGS=-fopenmp
OMPLIBS=-lgomp

# set linking libraries
LIBS = -lm

# flags for MATLAB MEX compilation..
MFLAGS=-compatibleArrayDims -DMWF77_UNDERSCORE1 
MWFLAGS=-c99complex 
MOMPFLAGS = -D_OPENMP

# location of MATLAB's mex compiler
MEX=mex

# For experts, location of Mwrap executable
MWRAP=../../mwrap/mwrap
MEXLIBS=-lm -lstdc++ -ldl -lgfortran

BC_INSTALL_DIR=$(PREFIX)
ifeq ($(PREFIX),)
	BC_INSTALL_DIR = ${HOME}/lib
endif

DYLIBS = $(LIBS)

LIBNAME=libboxcode2dlegacy
DYNAMICLIB = $(LIBNAME).so
STATICLIB = $(LIBNAME).a
LIMPLIB = $(DYNAMICLIB)

LLINKLIB = -lboxcode2dlegacy


-include make.inc

ifneq ($(OMP),OFF)
FFLAGS += $(OMPFLAGS)
MEXLIBS += $(OMPLIBS)

endif

SRC = ./src
TEST = ./test


OBJS = $(SRC)/lbfmm2d.o \
	$(SRC)/chebrouts.o \
	$(SRC)/poisson8.o \
	$(SRC)/tables8.o \
	$(SRC)/tree_routs4.o \
	$(SRC)/tree_vol_coeffs_cheb.o \
	$(SRC)/prini_new.o \
	$(SRC)/legeexps.o \
	$(SRC)/chebexps.o \
	$(SRC)/legetens.o \
	$(SRC)/chebtens.o \
	$(SRC)/voltab2d.o \
	$(SRC)/lapack_f77.o \

.PHONY: usage install lib test matlab

default: usage


usage:
	@echo "Makefile for boxcode2d-legacy. Specify what to make:"
	@echo "  make install - compile and install the main library"
	@echo "  make install PREFIX=(INSTALL_DIR) - compile and install the main library at custom location given by PREFIX"
	@echo "  make lib - compile the main library (in lib/ and lib-static/)"
	@echo "  make test - compile and run validation tests (will take a couple of mins)"
	@echo "  make matlab - compile matlab interfaces"
	@echo "  make objclean - removal all object files, preserving lib & MEX"
	@echo "  make clean - also remove lib, MEX, py, and demo executables"
	@echo "For faster (multicore) making, append the flag -j"
	@echo "  'make [task] OMP=OFF' for single-threaded"

#
# use only the file part of the filename, then manually specify
# the build location
#

$(SRC)/lbfmm2d.o : $(SRC)/lbfmm2d.f
	$(FC) -c $(FFLAGS) $(SRC)/lbfmm2d.f -o $(SRC)/lbfmm2d.o

$(SRC)/chebrouts.o : $(SRC)/chebrouts.f
	$(FC) -c $(FFLAGS) $(SRC)/chebrouts.f -o $(SRC)/chebrouts.o

$(SRC)/poisson8.o : $(SRC)/poisson8.f
	$(FC) -c $(FFLAGS) $(SRC)/poisson8.f -o $(SRC)/poisson8.o

$(SRC)/tables8.o : $(SRC)/tables8.f
	$(FC) -c $(FFLAGS2) $(SRC)/tables8.f -o $(SRC)/tables8.o

$(SRC)/tree_routs4.o : $(SRC)/tree_routs4.f
	$(FC) -c $(FFLAGS) $(SRC)/tree_routs4.f -o $(SRC)/tree_routs4.o

$(SRC)/tree_vol_coeffs_cheb.o : $(SRC)/tree_vol_coeffs_cheb.f
	$(FC) -c $(FFLAGS) $(SRC)/tree_vol_coeffs_cheb.f -o $(SRC)/tree_vol_coeffs_cheb.o

$(SRC)/prini_new.o : $(SRC)/prini_new.f
	$(FC) -c $(FFLAGS) $(SRC)/prini_new.f -o $(SRC)/prini_new.o

$(SRC)/legeexps.o : $(SRC)/legeexps.f
	$(FC) -c $(FFLAGS) $(SRC)/legeexps.f -o $(SRC)/legeexps.o

$(SRC)/chebexps.o : $(SRC)/chebexps.f
	$(FC) -c $(FFLAGS) $(SRC)/chebexps.f -o $(SRC)/chebexps.o

$(SRC)/legetens.o : $(SRC)/legetens.f
	$(FC) -c $(FFLAGS) $(SRC)/legetens.f -o $(SRC)/legetens.o

$(SRC)/chebtens.o : $(SRC)/chebtens.f
	$(FC) -c $(FFLAGS) $(SRC)/chebtens.f -o $(SRC)/chebtens.o

$(SRC)/voltab.o : $(SRC)/voltab.f
	$(FC) -c $(FFLAGS) $(SRC)/voltab.f -o $(SRC)/voltab.o

$(SRC)/lapack_f77.o : $(SRC)/lapack_f77.f
	$(FC) -c $(FFLAGS) $(SRC)/lapack_f77.f -o $(SRC)/lapack_f77.o



# build the library...
lib: $(STATICLIB) $(DYNAMICLIB)
ifneq ($(OMP),OFF)
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, multithread versions"
else
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, single-threaded versions"
endif


install: $(STATICLIB) $(DYNAMICLIB)
	echo $(BC_INSTALL_DIR)
	mkdir -p $(BC_INSTALL_DIR)
	cp -f lib/$(DYNAMICLIB) $(BC_INSTALL_DIR)/
	cp -f lib-static/$(STATICLIB) $(BC_INSTALL_DIR)/
	[ ! -f lib/$(LIMPLIB) ] || cp lib/$(LIMPLIB) $(BC_INSTALL_DIR)/
	@echo "Make sure to include " $(BC_INSTALL_DIR) " in the appropriate path variable"
	@echo "    LD_LIBRARY_PATH on Linux"
	@echo "    PATH on windows"
	@echo "    DYLD_LIBRARY_PATH on Mac OSX (not needed if default installation directory is used"
	@echo " "
	@echo "In order to link against the dynamic library, use -L"$(FMM_INSTALL_DIR) " -lboxcode2dlegacy"


$(STATICLIB): $(OBJS) 
	ar rcs $(STATICLIB) $(OBJS)
	mv $(STATICLIB) lib-static/
$(DYNAMICLIB): $(OBJS) 
	$(FC) -shared -fPIC $(OBJS) -o $(DYNAMICLIB) $(DYLIBS) 
	mv $(DYNAMICLIB) lib/
	[ ! -f $(LIMPLIB) ] || mv $(LIMPLIB) lib/

test: $(STATICLIB) test/lbfmm2d test/poisson8
	cd test; ./runsh.sh

test/lbfmm2d:
	$(FC) $(FFLAGS) test/test_lbfmm2d.f $(OBJS) -o test/int2-lb $(LIBS) 

test/poisson8:
	$(FC) $(FFLAGS) test/test_poisson8.f $(OBJS) -o test/int2-p $(LIBS) 


# matlab ..
MWRAPFILE = poisson8_matlab
GATEWAY = $(MWRAPFILE)


matlab:	$(STATICLIB) matlab/$(GATEWAY).c 
	$(MEX) matlab/$(GATEWAY).c lib-static/$(STATICLIB) $(MFLAGS) \
	-output matlab/poisson8_matlab $(MEXLIBS) 


mex:  $(STATICLIB)
	cd matlab; $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY) -mb $(MWRAPFILE).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY) -c $(GATEWAY).c $(MWRAPFILE).mw;\
	$(MEX) $(GATEWAY).c ../lib-static/$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE) \
	$(MEXLIBS); \

clean: objclean
	rm -f lib-static/*.a lib/*.so lib/*.dll
	rm -f matlab/*.mex*
	rm -f test/fort*
	rm -f test/int*

objclean:
	rm -f $(OBJS)
	rm -f test/*.o

list: $(SOURCES)
	$(warning Requires:  $^)



