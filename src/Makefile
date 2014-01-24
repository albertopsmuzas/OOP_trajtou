# Include general rules
debugmode=on
include make.inc
# Add here all folder names in src
folderlist= debug constants units surface lapackcontrol PES interpol1d

# Object list. There should be a rule for each object in this list
objectlist= debug.o constants.o lapackcontrol.o units.o surface.o PES.o PES_CRP.o interpol1d.o

# Where to find fortran files:
vpath %.f90 $(folderlist)

# DEPENDENCIES SECTION ###################################################

# Main rules
libtrajtou.a: $(objectlist)
	cd ../bin; ar rvs libtrajtou.a $(objectlist); mv libtrajtou.a ../lib/.
# Each object file dependencies
debug.o:
constants.o:
lapackcontrol.o: 
PES.o:
units.o: debug.o constants.o
surface.o: debug.o units.o constants.o lapackcontrol.o
zinterpol.o: debug.o lapackcontrol.o
PES_CRP.o: debug.o interpol1d.o PES.o surface.o 
# Clean rule
.PHONY: clean
clean:
	rm -f ../bin/*.o ../mod/*.mod ../lib/*.a
