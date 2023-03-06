# Makefile 
#

.SUFFIXES:
.SUFFIXES: .o .f90

include ../user_build_config

OBJS =	../build_upp/TABLE.o \
        ../build_upp/CTLBLK.o \
        ../build_upp/LOOKUP.o \
        ../build_upp/MASKS_mod.o \
        ../build_upp/TTBLEX.o \
        ../build_upp/VRBLS2D_mod.o \
        ../build_upp/VRBLS3D_mod.o \
        ../build_upp/params.o \
        ../build_upp/SPLINE.o \
        ../build_upp/TABLEQ.o \
        ../src/UPP_PHYSICS_minimal.o \
        ../src/cape_driver.o 

all:	cape_driver.exe

cape_driver.exe: $(OBJS)
	$(COMPILERF90) -o $(@) $(F90FLAGS) $(FREESOURCE) \
	-I../src -I../build_upp $(NETCDFMOD) $(OBJS) $(NETCDFLIB)

clean:
	rm -f *.exe

