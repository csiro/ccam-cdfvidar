CMP=ifort
FFLAGS = -O
INC = -I/park/tools/netcdf/3.6.1/include
LIBS = -L/park/tools/netcdf/3.6.1/lib -lnetcdf

OBJ2= dryadj.o findxn.o filt.o sintp16.o vidar.o invert.o\
	cdfvidar.o vispl.o esmtrv.o amap.o mslp.o lconset.o tot.o tomod.o \
	maxmin.o setxyz.o staguv.o jimcc.o fill.o outcdf.o prt_pan.o 

cdfvidar : $(OBJ2)
	$(CMP) $(FFLAGS) $(OBJ2) $(LIBS) -o cdfvidar

clean:
	rm -f *.o core cdfvidar

.f.o:
	$(CMP) -c $(FFLAGS) $(INC) $<
.f90.o:
	$(CMP) -c $(FFLAGS) $(INC) $<
%.o : %.mod

gribvida.o alphagrib.o pvidar.o sintp16.o tomod.o tot.o vidar.o cdfvidar.o: newmpar.h
cdfvidar.o one.o sintp16.o setxyz.o: latlong.h
cdfvidar.o one.o sintp16.o sintp16x.f: gblparm.h
vispl.o dryadj.o : lmax.h
cdfvidar.o vidar.o vispl.o : nplevs.h
cdfvidar.o vidar.o : vidar.h
utilities.o : utilities.f90 
