# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .f .F

include ../user_build_config

OBJS =	TABLE.o \
        CTLBLK.o \
        LOOKUP.o \
        MASKS_mod.o \
        TTBLEX.o \
        VRBLS2D_mod.o \
        VRBLS3D_mod.o \
        params.o \
        SPLINE.o \
        TABLEQ.o
	
all:	$(OBJS)

TABLE.o: $(UPPDIR)/TABLE.f
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(UPPDIR)/TABLE.f

CTLBLK.o: $(UPPDIR)/CTLBLK.f
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(UPPDIR)/CTLBLK.f

LOOKUP.o: $(UPPDIR)/LOOKUP.f
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(UPPDIR)/LOOKUP.f

MASKS_mod.o: $(UPPDIR)/MASKS_mod.f
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(UPPDIR)/MASKS_mod.f

TTBLEX.o: $(UPPDIR)/TTBLEX.f
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(UPPDIR)/TTBLEX.f

VRBLS2D_mod.o: $(UPPDIR)/VRBLS2D_mod.f
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(UPPDIR)/VRBLS2D_mod.f

VRBLS3D_mod.o: $(UPPDIR)/VRBLS3D_mod.f
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(UPPDIR)/VRBLS3D_mod.f

params.o: $(UPPDIR)/params.F
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(UPPDIR)/params.F

SPLINE.o: $(UPPDIR)/SPLINE.f
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(UPPDIR)/SPLINE.f

TABLEQ.o: $(UPPDIR)/TABLEQ.f
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(UPPDIR)/TABLEQ.f

clean:
	rm -f *.o *.mod

#
# Dependencies:
#
TTBLEX.o:		        CTLBLK.o
TABLE.o:                        SPLINE.o
TABLEQ.o:                       SPLINE.o
