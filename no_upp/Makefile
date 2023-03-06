
# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .f .F

COMPILERF90    =       ifort
FREESOURCE     =       -free
F90FLAGS       =       #-fdefault-real-8 -fdefault-double-8

OBJS =	TABLE.o CTLBLK.o LOOKUP.o MASKS_mod.o TTBLEX.o VRBLS2D_mod.o \
        VRBLS3D_mod.o params.o UPP_PHYSICS_minimal.o SPLINE.o TABLEQ.o

all:	cape_driver.exe

.f.o:
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(*).f

.F.o:
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(*).F

cape_driver.exe: $(OBJS) cape_driver.f90
	$(COMPILERF90) -o $(@) $(F90FLAGS) $(FREESOURCE) \
	cape_driver.f90 $(OBJS)

clean:
	rm -f *.o *.exe *.mod
#
# Dependencies:
#
TTBLEX.o:		        CTLBLK.o
UPP_PHYSICS_minimal.o:		TABLE.o CTLBLK.o LOOKUP.o MASKS_mod.o \
                                TTBLEX.o VRBLS2D_mod.o VRBLS3D_mod.o params.o 
TABLE.o:                        SPLINE.o
TABLEQ.o:                       SPLINE.o
