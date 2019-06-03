
SRCDIR=trunk
MODDIR=trunk/module

#Default path for FFTW
FFTWINC=/usr/local/lib/fftw-3.3.7/include
FFTWLIB=/usr/local/lib/fftw-3.3.7/lib

FC=gfortran
FFLAGS=-ffree-form -cpp -O3  -fPIC -fno-second-underscore -m64 -Wl,--no-as-needed

vpath %.f90 $(SRCDIR)

%.o:%.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) -L$(FFTWLIB) -I$(FFTWINC)  -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lgomp -lpthread -lm -ldl -lfftw3 -c -o $(MODDIR)/$@ $<	

.PHONY:	all	clean

all: phys_cte.o	algo.o geometry_defs.o geometry_func.o mode_solver.o scan_process.o

clean:
	@echo "Cleaning .o and .mod files in "$(MODDIR)"."; \
	rm -f *.o $(MODDIR)/*.o	\
	rm -f *.o $(MODDIR)/*.mod
	






