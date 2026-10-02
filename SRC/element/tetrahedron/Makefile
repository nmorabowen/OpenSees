include ../../../Makefile.def

OBJS       = FourNodeTetrahedron.o TenNodeTetrahedron.o TenNodeTetrahedronThermal.o ThermalBoundaryTri6.o


all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
