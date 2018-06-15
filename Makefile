#
include make.inc
#
#F90=mpiifort -fc=ifort
F90=mpif90 -fc=ifort
F90OPT=$(F90) # -xSSE4.2
FLIBS = $(LIB) $(MPIlib)

ELPALIB=/home/lsg/local/ParEig/lib/lib/libelpa_openmp.a
ELPAMOD=-I/home/lsg/local/ParEig/lib/include/elpa_openmp-2016.05.003/modules/ -I/media/lsg/娱乐/GZ_20180527/GZ_ELPA/elpa-2016.05.003/private_modules
LIn= -qopenmp -I/home/opt/mkl/include/intel64/lp64


exe=test_zolodrv1
#exe=test_pqdwh

# ------------------------------------------------------------------------------
# qdwh 
#
#OBJ    = prepare_matrix.o pqdwh_fall.o pqdwhfacs.o hfunc.o pqdwhsubit.o  mpdgeqrf.o mpdorgqr.o \
          pdmatgen2.o pmatgeninc.o elpa_pdsyevd2.o 
# zolosvd
OBJ = check_routines.o prepare_matrix.o Auxil_Zolo.o test_zolodrv1.o pzolopd1.o zoloqr.o zolochol.o \
      FormX2.o ZoloCHK1.o ZoloCHK.o mpdgeqrf.o mpdorgqr.o elpa_pdsyevd2.o

# mqr
OBJ1 = mpdgeqrf.o mpdorgqr.o test_mpdgeqrf.o

# QR
OBJ_qr = test_qr.o


$(exe): $(OBJ) 
	$(F90) $(FFLAGS) -o $@ $(OBJ) $(ELPAMOD) $(ELPALIB) $(FLIBS) $(LIn)

test_qr: $(OBJ)
	$(F90) $(FFLAGS) -o $@ $(OBJ) $(FLIBS)

test_mpdgeqrf: $(OBJ1)
	$(F90) $(FFLAGS) -o $@ $(OBJ1) $(FLIBS)

elpa_pdsyevd2.o: elpa_pdsyevd2.F90
	$(F90) $(FFLAGS) -c $< $(ELPAMOD) $(LIn)

elpa_pdsyevd.o: elpa_pdsyevd.F90
	$(F90) $(FFLAGS) -c $< $(ELPAMOD) $(LIn)

%.o:%.f90
	$(F90) -c $(FFLAGS) $<

%.o:%.f
	$(F90) -c $(FFLAGS) $<

clean:
	rm -f *.o *.mod test pqdwh_caller
