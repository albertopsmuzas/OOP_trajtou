##################################################
# MAKE_FILE FOR  OOP TRAJTOU #####################
##################################################
include $(OOPTRAJTOUPATH)/make.compiler.inc
include $(OOPTRAJTOUPATH)/include/make.general.inc
.PHONY : trajtouJobs clean trajtouSrc
vpath %.f90 jobs
joblist = test_inicond3dInput.x\
			 trajtouDynamics3D.x\
			 trajtouDynamics6D.x\
			 trajtouGetPot_crp6d.x\
			 CRP6D_gridstats.ZR.x\
			 trajtouMolec2Atom.x\
			 au2angst.x\
			 get_CRP6D_shift.x\
			 test_CRP6D_cuts.x\
			 get_CRP6D_smoothvalue.x\
			 get_CRP6D_rawvalue.x\
			 cheat_CRP6D_cartwheelinput.x\
			 cart2surf.x\
			 test_inicond6dInput.x\
			 drawtraj6d.x anlys_diff3d.x\
			 anlys_diff6d.x\
			 anlys_diffAtomSurf.x\
			 test_systemInput.x\
			 test_crp3dInput.x\
			 trajtouGetPot_crp3d.x\
			 test_crp6dInput.x\
			 test_CRP6D.x\
			 test_dynamics3dInput.x\
			 test_dynamics6dInput.x\
			 trajtouAtom2Molec.x\
			 trajtouNewSeed.x\
			 trajtouZcut_crp6d.x\
			 test_inicondAtomSurfInput.x\
			 trajtouDynamicsAtomSurf.x\
			 trajtouEvaluateEnergyRovibrState.x\
			 test_CRP6D_bugs.x trajtouGetPot.x\
			 trajtouZcut_rawInterpolMolec.x\
			 trajtouPrintPolygon.x\
			 anlys_growMolec.x\
			 test_inicondGrowMolec.x\
			 test_inicondGrowMolecFromFile.x\
			 trajtouGetCut_ZR_crp6d.x
# Rules 
trajtouJobs: trajtouSrc $(joblist) 
trajtouSrc:
	cd src && make default

%.x: %.f90 libtrajtou.a libaotus.a
	$(LF) -o bin/$@ $< lib/libtrajtou.a lib/libaotus.a

%.a:
	cd src && make $@ 
clean:
	rm -f mod/*.mod lib/*.a bin/*.x bin/*.o
