# Include general rules
debugmode=on
include make.inc
# Add here all folder names in src
folderlist= debug constants units surface PES interpol1d interpol1d/cubicsplines maths PES/CRP3D interpol2d interpol2d/fourier \
				inicond inicond/atom dynamics dynamics/atom interpol2d/bicsplines PES/CRP6D interpol1d/fourier

# Object list. There should be a rule for each object in this list
objectlist= debug.o constants.o units.o surface.o PES.o CRP3D.o interpol1d.o maths.o fourier2d.o fourier.p4mm.o cubicsplines.o \
				inicond.o initatom.o dynamics.o dynatom.o bicsplines.o CRP6D.o fourier1d.o wyckoff.generic.o wyckoff.p4mm.o \
				extrapol.o interpolgrid2d.o

# Where to find fortran files:
vpath %.f90 $(folderlist)

# DEPENDENCIES SECTION ###################################################

# Main rules
libtrajtou.a: $(objectlist)
	cd ../bin; ar rvs libtrajtou.a $(objectlist); mv libtrajtou.a ../lib/.
# Each object file dependencies
debug.o:
constants.o:
PES.o:
maths.o:
interpol1d.o: 
interpolgrid2d.o: 
fourier2d.o:
inicond.o:
dynamics.o: units.o
bicsplines.o: maths.o constants.o cubicsplines.o interpolgrid2d.o
dynatom.o: debug.o dynamics.o CRP3D.o
initatom.o: debug.o constants.o units.o inicond.o surface.o
units.o: debug.o constants.o
surface.o: debug.o units.o constants.o maths.o 
cubicsplines.o: debug.o interpol1d.o
fourier.p4mm.o: debug.o constants.o surface.o fourier2d.o
CRP3D.o: debug.o units.o maths.o cubicsplines.o fourier.p4mm.o PES.o surface.o 
CRP6D.o: debug.o units.o maths.o CRP3D.o wyckoff.p4mm.o extrapol.o
fourier1d.o: constants.o debug.o interpol1d.o
wyckoff.generic.o: debug.o fourier1d.o bicsplines.o
wyckoff.p4mm.o: wyckoff.generic.o
extrapol.o: units.o debug.o maths.o
# Clean rule
.PHONY: clean
clean:
	rm -f ../bin/*.o ../mod/*.mod ../lib/*.a
