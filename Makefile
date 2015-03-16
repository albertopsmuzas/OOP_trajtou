##################################################
# MAKE_FILE FOR  OOP TRAJTOU #####################
##################################################
include make.inc
.PHONY : build
joblist = test_inicond3dInput.x trajtouDynamics3D.x trajtouDynamics6D.x trajtouGetPot_crp6d.x\
			 CRP6D_gridstats.ZR.x trajtouMolec2Atom.x au2angst.x get_CRP6D_shift.x\
			 test_CRP6D_cuts.x get_CRP6D_smoothvalue.x get_CRP6D_rawvalue.x\
			 cheat_CRP6D_cartwheelinput.x cart2surf.x test_inicond6dInput.x\
			 drawtraj6d.x anlys_diff3d.x anlys_diff6d.x anlys_diffAtomSurf.x test_systemInput.x\
			 test_crp3dInput.x trajtouGetPot_crp3d.x test_crp6dInput.x test_CRP6D.x\
			 test_dynamics3dInput.x test_dynamics6dInput.x trajtouAtom2Molec.x\
			 trajtouNewSeed.x trajtouZcut_crp6d.x test_inicondAtomSurfInput.x trajtouDynamicsAtomSurf.x\
			 trajtouEvaluateEnergyRovibrState.x test_CRP6D_bugs.x trajtouGetPot.x
# Rules 
link_jobs2lib: $(joblist)
lib: 
	cd src && $(MAKE) libtrajtou.a

# Clean rule
.PHONY : clean
clean:
	rm -f mod/*.mod lib/*.a bin/*.x bin/*.o
