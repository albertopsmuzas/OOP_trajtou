# Include general rules
include make.inc
.PHONY: clean
# Add here all folder names in src
folderlist= debug constants units surface lapackcontrol

# Object list. There should be a rule for each object in this list
objectlist= debug.o constants.o lapackcontrol.o units.o surface.o 

# Where to find fortran files:
vpath %.f90 $(folderlist)

# DEPENDENCIES SECTION ###################################################

# Main rules
testing: $(objectlist)
build: $(objectlist)
	cd ../bin; ar rvs libtrajtou.a $(objectlist); mv libtrajtou.a ../lib/.

# Each object file dependencies
debug.o:
constants.o:
lapackcontrol.o: 
units.o : debug.o constants.o
surface.o: debug.o units.o constants.o lapackcontrol.o

# Clean rule
clean:
	rm -f ../bin/*.o ../mod/*.mod ../lib/*.a
