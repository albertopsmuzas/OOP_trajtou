# Include general rules
debugmode=on
include make.inc
# Add here all folder names in src
folderlist= debug constants units surface lapackcontrol PES interpol1d interpol1d/cubicsplines maths PES/CRP interpol2d interpol2d/fourier

# Object list. There should be a rule for each object in this list
objectlist= debug.o constants.o lapackcontrol.o units.o surface.o PES.o CRP.o interpol1d.o maths.o interpol2d.o \
				fourier2d.o cubicsplines.o
				

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
maths.o:
interpol1d.o: 
interpol2d.o: 
units.o: debug.o constants.o
surface.o: debug.o units.o constants.o lapackcontrol.o
cubicsplines.o: debug.o interpol1d.o
fourier2d.o: debug.o constants.o surface.o lapackcontrol.o interpol2d.o
CRP.o: debug.o maths.o cubicsplines.o fourier2d.o PES.o surface.o 
# Clean rule
.PHONY: clean
clean:
	rm -f ../bin/*.o ../mod/*.mod ../lib/*.a
