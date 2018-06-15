# Makefile.inc example. Works with the GNU compilers & OpenMPI.

# Flags to be passed to compilers (C++, C, Fortran, mex)
# CXXFLAGS is mandatory. Flags for the C++ compiler.
# CFLAGS is optional (used only for the C and Fortran interface).
# FFLAGS is optional (used only for the Fortran interface, or with -DHQR). 
# MEXFLAGS is optional (used only for the Octave/Matlab interface).

MKLROOT = /home/opt/mkl
CXXFLAGS  = -O2 -qopenmp -I$(MKLROOT)/include -mavx 
CFLAGS    = -O2 -qopenmp -I$(MKLROOT)/include -mavx
FFLAGS    = -O2 -qopenmp -I$(MKLROOT)/include -mavx # -DOUTPUT_Mat

#CXXFLAGS  = -O1 
#CFLAGS    = -O1 
#FFLAGS    = -O2
# -traceback

MEXFLAGS  = --mex

# Compilers  (C++, C, Fortran, mex)
# CXX is mandatory. C++ compiler.
# CC is optional (used only for the C and Fortran interface).
# FC is optional (used only for the Fortran interface, or with -DHQR). 
# MEXis optional (used only for the Octave/Matlab interface).
CXX       = mpicxx
CC        = mpicc
FC        = mpif90
MEX       = mkoctfile

# Libraries:
# LIB should provide ScaLAPACK, BLACS, LAPACK, BLAS
# LIBCXX is used for the Fortran interface, to provide C++ libraries to the linker
# LIBMEX is used for the mex interface
#LIB = -L/home/lsg/library/scalapack-2.0.2/ -lscalapack -L/home/lsg/library/lapack-3.6.1 -llapack -lblas -lgfortran -lm

#LIB=$(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a \
#	$(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm
#BLAS = -mkl=parallel

LIB=$(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a \
	$(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm

#LIBCXX    = -lstdc++ -lmpi_cxx
LIBCXX    = -lstdc++ -lmpicxx
LIBMEX    = -L/usr/lib/openmpi -lmpi_cxx -lscalapack -lblas
