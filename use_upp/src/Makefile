
# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .f .f90

include ../user_build_config

OBJS =	UPP_PHYSICS_minimal.o \
        cape_driver.o

.f.o:
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) -I../build_upp $(NETCDFMOD) $(*).f

.f90.o:
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) -I../build_upp $(NETCDFMOD) $(*).f90

all:	$(OBJS)

clean:
	rm -f *.o *.exe *.mod
#
# Dependencies:
#
UPP_PHYSICS_minimal.o:  ../build_upp/TABLE.o ../build_upp/CTLBLK.o ../build_upp/LOOKUP.o \
                        ../build_upp/MASKS_mod.o ../build_upp/TTBLEX.o ../build_upp/TABLE.o \
                        ../build_upp/VRBLS2D_mod.o ../build_upp/VRBLS3D_mod.o ../build_upp/params.o 
cape_driver.o:          UPP_PHYSICS_minimal.o ../build_upp/TABLE.o ../build_upp/CTLBLK.o \
                        ../build_upp/LOOKUP.o ../build_upp/MASKS_mod.o ../build_upp/TABLEQ.o \
                        ../build_upp/VRBLS2D_mod.o ../build_upp/VRBLS3D_mod.o
