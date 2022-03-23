
ifneq ($(CUSTOM),yes)
FC = ifort
LIBS = -L $(NETCDF_ROOT)/lib -lnetcdf
ifneq ($(NCCLIB),yes)
LIBS += -lnetcdff
endif
HOST = -xHost
ifeq ($(CASCADELAKE),yes)
HOST = -xCASCADELAKE
endif
ifeq ($(ZEN3),yes)
HOST = -axCORE-AVX2
endif
INC = -I $(NETCDF_ROOT)/include
FFLAGS =  -qopenmp $(HOST) -fp-model precise -traceback
PPFLAG90 = -fpp
PPFLAG77 = -fpp
DEBUGFLAG = -check all -debug all -traceback -fpe0
endif

ifeq ($(GFORTRAN),yes)
FC = gfortran
FFLAGS = -O2 -mtune=native -march=native -I $(NETCDF_ROOT)/include
ifeq ($(ZEN3),yes)
FFLAGS = -O2 -fallow-argument-mismatch -mtune=native -march=native -I $(NETCDF_ROOT)/include
endif
PPFLAG90 = -x f95-cpp-input
PPFLAG77 = -x f77-cpp-input
DEBUGFLAG = -g -Wall -Wextra -fbounds-check -fbacktrace
endif

ifeq ($(CRAY),yes)
FC = ftn
FFLAGS = -h noomp
PPFLAG90 = -eZ
PPFLAG77 = -eZ
DEBUGFLAG =
endif

ifeq ($(MAUI),yes)
FC = ftn
HOST = -xSKYLAKE-AVX512
FFLAGS = -qopenmp $(HOST) -fp-model precise -traceback
PPFLAG90 = -fpp
PPFLAG77 = -fpp
DEBUGFLAG = -check all -debug all -traceback -fpe0
endif

# Testing - I/O and fpmodel
ifeq ($(TEST),yes)
FFLAGS += $(DEBUGFLAG)
endif

ifeq ($(NCCLIB),yes)
FFLAGS += -Dncclib
endif


OBJ2= dryadj.o findxn.o filt.o sintp16.o vidar.o invert.o\
      cdfvidar.o vispl.o esmtrv.o amap.o mslp.o \
      maxmin.o fill.o outcdf.o prt_pan.o \
      setxyz_m.o ccinterp.o jimcc_m.o \
      latltoij_m.o xyzinfo_m.o newmpar_m.o indices_m.o \
      parm_m.o precis_m.o ind_m.o jimco_m.o jim_utils.o nfft_m.o \
      latlong_m.o comsig_m.o cll_m.o sigdata_m.o netcdf_m.o \
      stacklimit.o

cdfvidar : $(OBJ2)
	$(FC) $(FFLAGS) $(OBJ2) $(LIBS) -o cdfvidar

clean:
	rm -f *.o core cdfvidar *.mod

.SUFFIXES:.f90

stacklimit.o: stacklimit.c
	cc -c stacklimit.c
version.h: FORCE
	rm -f brokenver tmpver
	echo "      character(len=*), parameter :: version ='CDFVIDAR r'" > brokenver
	echo "      character(len=*), parameter :: version ='CDFVIDAR r`svnversion .`'" > tmpver
	grep exported tmpver || grep Unversioned tmpver || cmp tmpver brokenver || cmp tmpver version.h || mv tmpver version.h
FORCE:

.f.o:
	$(FC) -c $(FFLAGS) $(INC) $(PPFLAG77) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INC) $(PPFLAG90) $<
%.o : %.mod

cdfvidar.o sintp16.o setxyz.o: latlong_m.o
cdfvidar.o vidar.o : comsig_m.o
cdfvidar.o outcdf.o : cll_m.o
cdfvidar.o vidar.o : sigdata_m.o
cdfvidar.o : ccinterp.o version.h
cdfvidar.o vispl.o dryadj.o : lmax.h
cdfvidar.o vidar.o vispl.o : nplevs.h
cdfvidar.o vidar.o : vidar.h
cdfvidar.o outcdf.o : netcdf_m.o
vidar.o : outcdf.o
utilities.o : utilities.f90
ccinterp.o : ccinterp.f90 setxyz_m.o xyzinfo_m.o latltoij_m.o newmpar_m.o indices_m.o precis_m.o
latltoij_m.o : latltoij_m.f90 xyzinfo_m.o newmpar_m.o precis_m.o
setxyz_m.o : setxyz_m.f90 newmpar_m.o indices_m.o parm_m.o precis_m.o ind_m.o xyzinfo_m.o jimco_m.o jimcc_m.o 
xyzinfo_m.o : xyzinfo_m.f90 precis_m.o
newmpar_m.o : newmpar_m.f90 
precis_m.o : precis_m.f90
indices_m.o : indices_m.f90
parm_m.o : parm_m.f90 precis_m.o 
ind_m.o : ind_m.f90 newmpar_m.o 
jimcc_m.o : jimcc_m.f90 parm_m.o precis_m.o xyzinfo_m.o
jimco_m.o : jimco_m.f90 precis_m.o jim_utils.o nfft_m.o 
jim_utils.o : jim_utils.f90 precis_m.o 
nfft_m.o : nfft_m.f90 precis_m.o 
