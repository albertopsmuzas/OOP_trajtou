# Include general rules
include make.inc
# Add here all folder names in src
folderlist= debug constants units surface PES interpol1d interpol1d/cubicsplines maths PES/CRP3D interpol2d interpol2d/fourier \
				inicond inicond/atom inicond/molecule dynamics dynamics/atom interpol2d/bicsplines PES/CRP6D interpol1d/fourier functions\
				dynamics/molecule utils analysis analysis/diffraction

# Object list. There should be a rule for each object in this list
objectlist= debug.o constants.o units.o surface.o PES.o CRP3D.o interpol1d.o maths.o fourier2d.o fourier.p4mm.o cubicsplines.o \
				inicond.o initatom.o dynamics.o dynatom.o bicsplines.o CRP6D.o fourier1d.o wyckoff.generic.o wyckoff.p4mm.o \
				extrapol.o interpolgrid2d.o fourier.4mm.o fourier.2.o fourier.mm2.o fourier.e.o fourier.m.o fourier.m45.o fourier.m45m1352.o \
				function1d.o logistic.o one.o link_PES.o link_interpol1d.o initdiatomic.o dyndiatomic.o link_function1d.o linear.o\
				drawtraj.o link_inicond.o diffractionCRP3D.o diffractionCRP6D.o fourier.p6mm.o link_fourier2d.o

# Where to find fortran files:
vpath %.f90 $(folderlist)

# DEPENDENCIES SECTION ###################################################

# Static library section ----------------------------------------------------------------
libtrajtou.a: $(objectlist)
	cd ../bin; ar rvs libtrajtou.a $(objectlist); mv libtrajtou.a ../lib/.
# System / Fundamental section ----------------------------------------------------------
debug.o:
maths.o:
constants.o:
units.o: debug.o constants.o
surface.o: debug.o units.o constants.o maths.o 
# Potential energy surfaces section -----------------------------------------------------
link_PES.o: PES.o CRP3D.o CRP6D.o
PES.o: surface.o
CRP3D.o: debug.o units.o maths.o cubicsplines.o link_fourier2d.o PES.o link_function1d.o
CRP6D.o: debug.o units.o maths.o CRP3D.o wyckoff.p4mm.o extrapol.o link_function1d.o
# Example of new PES implementation:
# newPES.o: PES.o +list of needed modules. Don't forget to add this routine to objectlist.
 
# 1D interpolations ---------------------------------------------------------------------
link_interpol1d.o: interpol1d.o fourier1d.o cubicsplines.o
interpol1d.o: 
fourier1d.o: debug.o interpol1d.o
cubicsplines.o: debug.o constants.o interpol1d.o
# Specific interpolations for CRP6D -----------------------------------------------------
fourier.4mm.o: debug.o fourier1d.o
fourier.mm2.o: debug.o fourier1d.o
fourier.m45m1352.o: debug.o fourier1d.o
fourier.2.o: debug.o fourier1d.o
fourier.e.o: debug.o fourier1d.o
fourier.m.o: debug.o fourier1d.o
fourier.m45.o: debug.o fourier1d.o
# 2D interpolations ----------------------------------------------------------------------
interpolgrid2d.o: 
bicsplines.o: maths.o constants.o cubicsplines.o interpolgrid2d.o
link_fourier2d.o: fourier2d.o fourier.p4mm.o fourier.p6mm.o
fourier2d.o: constants.o surface.o constants.o
fourier.p4mm.o: debug.o fourier2d.o
fourier.p6mm.o: debug.o fourier2d.o
# 1D Functions section -------------------------------------------------------------------
link_function1d.o: function1d.o one.o logistic.o linear.o exponential.o
function1d.o:
exponential.o: function1d.o
logistic.o: function1d.o
one.o: function1d.o
linear.o: function1d.o
# Initial conditions section -------------------------------------------------------------
link_inicond.o: inicond.o initatom.o initdiatomic.o
inicond.o: units.o constants.o link_PES.o 
initatom.o: debug.o inicond.o 
initdiatomic.o: debug.o inicond.o
# Dynamics section -----------------------------------------------------------------------
dynamics.o: units.o link_PES.o inicond.o
dynatom.o: debug.o dynamics.o initatom.o
dyndiatomic.o: debug.o dynamics.o initdiatomic.o
# Analysis -------------------------------------------------------------------------------
#diffraction.o: debug.o surface.o initatom.o CRP3D.o
diffractionCRP3D.o: debug.o surface.o initatom.o CRP3D.o
diffractionCRP6D.o: debug.o surface.o initdiatomic.o CRP6D.o
# Utils ----------------------------------------------------------------------------------
drawtraj.o: surface.o
#-----------------------------------------------------------------------------------------
wyckoff.generic.o: debug.o fourier1d.o bicsplines.o
wyckoff.p4mm.o: wyckoff.generic.o fourier.4mm.o fourier.mm2.o fourier.2.o fourier.e.o fourier.m.o fourier.m45.o fourier.m45m1352.o
extrapol.o: units.o debug.o maths.o
# Clean rule
.PHONY: clean
clean:
	rm -f ../bin/*.o ../mod/*.mod ../lib/*.a
