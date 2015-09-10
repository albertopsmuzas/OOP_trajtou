# Include general rules
include make.inc
# Add here all folder names in src
folderlist= debug units surface PES interpol1d interpol1d/cubicsplines maths PES/CRP3D interpol2d interpol2d/fourier \
				inicond inicond/atom inicond/molecule dynamics dynamics/atom interpol2d/bicsplines PES/CRP6D interpol1d/fourier functions\
				dynamics/molecule utils analysis analysis/diffraction integrator system inicond/atomSurf dynamics/atomSurf PES/Specific\
				PES/Specific/HLiF001 PES/Specific/H2LiF001 interpol3d interpol3d/fourierCombined

# Object list. There should be a rule for each object in this list
objectlist= debug.o units.o surface.o PES.o CRP3D.o interpol1d.o maths.o fourier2d.o fourier.p4mm.o cubicsplines.o \
				inicond.o initatom.o dynamics.o dynatom.o bicsplines.o CRP6D.o fourier1d.o wyckoff.generic.o wyckoff.p4mm.o \
				extrapol.o interpolgrid2d.o fourier.4mm.o fourier.2.o fourier.mm2.o fourier.e.o fourier.m.o fourier.m45.o fourier.m45m1352.o \
				function1d.o logistic.o one.o link_PES.o link_interpol1d.o initdiatomic.o dyndiatomic.o link_function1d.o linear.o\
				drawtraj.o link_inicond.o diffractionCRP3D.o diffractionCRP6D.o link_fourier2d.o exponential.o xexponential.o\
				xplusnexponential.o integrator_o1.o link_integrator.o bulirsch_stoer.o system.o initAtomSurf.o dynAtomSurf.o\
				diffractionAtomSurf.o interpol3d.o fourier3dCombined.o fourier3d.p4mm.o initSurfaceGrowDiatomic.o diffractionGrow.o 
				

# Where to find fortran files:
vpath %.f90 $(folderlist)

default: all
all: libaotus.a libtrajtou.a

# Static library: trajtou ----------------------------------------------------------------
libtrajtou.a: $(objectlist)
	cd ../bin; ar rvs libtrajtou.a $(objectlist); mv libtrajtou.a ../lib/.
# Static library: aotus --------------------------------------------------------------
libaotus.a: 
	cd utils/haraldkl-aotus-707f4da5293b;\
	./waf configure build --check-fortran-compiler=$(F);\
	mv build/*.a ../../../lib/.;\
	mv build/*.mod ../../../mod/.;\
	mv build/lua ../../../bin/lua
# System / Fundamental section ----------------------------------------------------------
debug.o:
maths.o:
units.o: debug.o 
system.o: debug.o units.o surface.o maths.o
# Surface section
surface.o: debug.o units.o maths.o
# Potential energy surfaces section -----------------------------------------------------
link_PES.o: PES.o CRP3D.o CRP6D.o
PES.o: surface.o
CRP3D.o: debug.o units.o maths.o PES.o cubicsplines.o link_fourier2d.o link_function1d.o system.o
CRP6D.o: debug.o units.o maths.o PES.o CRP3D.o wyckoff.p4mm.o extrapol.o link_function1d.o fourier3d.p4mm.o
# Example of new PES implementation:
# newPES.o: PES.o +list of needed modules. Don't forget to add this routine to objectlist.
 
# 1D interpolations ---------------------------------------------------------------------
link_interpol1d.o: interpol1d.o fourier1d.o cubicsplines.o
interpol1d.o: 
fourier1d.o: debug.o interpol1d.o
cubicsplines.o: debug.o units.o interpol1d.o
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
bicsplines.o: maths.o units.o cubicsplines.o interpolgrid2d.o
link_fourier2d.o: fourier2d.o fourier.p4mm.o 
fourier2d.o: 
fourier.p4mm.o: debug.o fourier2d.o maths.o units.o system.o
# 3D interpolations
interpol3d.o: debug.o
fourier3dCombined.o: debug.o interpol3d.o maths.o fourier1d.o fourier2d.o
fourier3d.p4mm.o: debug.o fourier3dCombined.o fourier.4mm.o fourier.p4mm.o  
# 1D Functions section -------------------------------------------------------------------
link_function1d.o: function1d.o one.o logistic.o linear.o exponential.o xexponential.o xplusnexponential.o
function1d.o:
exponential.o: function1d.o
xexponential.o: function1d.o
xplusnexponential.o: function1d.o
logistic.o: function1d.o
one.o: function1d.o
linear.o: function1d.o
# Initial conditions section -------------------------------------------------------------
link_inicond.o: inicond.o initatom.o initdiatomic.o
inicond.o: units.o link_PES.o 
initatom.o: debug.o inicond.o system.o
initdiatomic.o: debug.o inicond.o system.o
initAtomSurf.o: debug.o inicond.o system.o
initSurfaceGrowDiatomic.o: debug.o inicond.o system.o
# Dynamics section -----------------------------------------------------------------------
dynamics.o: units.o link_PES.o inicond.o link_integrator.o
dynatom.o: debug.o dynamics.o initatom.o system.o
dynAtomSurf.o: debug.o dynamics.o initAtomSurf.o system.o
dyndiatomic.o: debug.o dynamics.o initdiatomic.o system.o
# Analysis -------------------------------------------------------------------------------
diffractionCRP3D.o: debug.o  initatom.o CRP3D.o system.o
diffractionCRP6D.o: debug.o  initdiatomic.o CRP6D.o system.o
diffractionAtomSurf.o: debug.o initAtomSurf.o CRP3D.o system.o
diffractionGrow.o: debug.o initSurfaceGrowDiatomic.o system.o
# Integrator -------------------------------------------------------------------------------
link_integrator.o: integrator_o1.o bulirsch_stoer.o
integrator_o1.o:
bulirsch_stoer.o: integrator_o1.o
# Utils ----------------------------------------------------------------------------------
drawtraj.o: debug.o system.o
#-----------------------------------------------------------------------------------------
wyckoff.generic.o: debug.o fourier1d.o bicsplines.o
wyckoff.p4mm.o: wyckoff.generic.o fourier.4mm.o fourier.mm2.o fourier.2.o fourier.e.o fourier.m.o fourier.m45.o fourier.m45m1352.o
extrapol.o: units.o debug.o maths.o
# Clean rule
.PHONY: clean
clean:
	rm -f ../bin/*.o ../mod/*.mod ../lib/*.a
