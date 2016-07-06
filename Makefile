##################################################
# MAKE_FILE FOR  OOP TRAJTOU #####################
##################################################
ifeq ($(OOPTRAJTOUPATH),) 
	export OOPTRAJTOUPATH=$(PWD)
endif
include $(OOPTRAJTOUPATH)/make.compiler.inc
include $(OOPTRAJTOUPATH)/include/make.general.inc
.PHONY : trajtouJobs clean trajtouSrc
vpath %.f90 jobs
joblist = test_inicond3dInput.x\
			 CRP6D_gridstats.ZR.x\
			 au2angst.x\
			 get_CRP6D_shift.x\
			 test_CRP6D_cuts.x\
			 get_CRP6D_smoothvalue.x\
			 get_CRP6D_rawvalue.x\
			 cheat_CRP6D_cartwheelinput.x\
			 cart2surf.x\
			 drawtraj6d.x\
			 test_inicond6dInput.x\
			 test_systemInput.x\
			 test_crp3dInput.x\
			 test_crp6dInput.x\
			 test_dynamics3dInput.x\
			 test_dynamics6dInput.x\
			 test_inicondGrowMolec.x\
			 test_inicondGrowMolecFromFile.x\
			 test_inicondAtomSurfInput.x\
			 test_CRP6D_bugs.x\
			 trajtouGetPot.x\
			 trajtouAtom2Molec.x\
			 trajtouNewSeed.x\
			 trajtouMolec2Atom.x\
			 trajtouDynamicsAtomSurf.x\
			 trajtouEvaluateEnergyRovibrState.x\
			 trajtouPrintPolygon.x\
			 trajtouDynamics3D.x\
			 trajtouDynamics6D.x\
			 trajtouAnalysis_diffraction_crp3d.x\
			 trajtouAnalysis_diffraction_crp6d.x\
			 trajtouAnalysis_diffraction_grow6d.x\
			 trajtouAnalysis_diffraction_crp3d_withSpring.x\
			 trajtouGetPot_crp6d.x\
			 trajtouGetPot_crp3d.x\
			 trajtouGetInfo_vacuumPot.x\
			 trajtouGetInfo_crp6d_vacuumPot.x\
			 trajtouGetInfo_vacuumPotCoeff.x\
			 trajtouGetGraph_vacuumPot.x\
			 trajtouGetGraph_vacuumPotShifted.x\
			 trajtouGetGraph_XYcut_crp3d.x\
			 trajtouGetGraph_Zcut_crp3d.x\
			 trajtouGetGraph_Zcut_crp6d.x\
			 trajtouGetGraph_ZRcut_crp6d.x\
			 trajtouGenerateCRP3DInput.x\
			 trajtouSymmetrizeCRP3DInput.x\
			 trajtouTest_interpolationAtSitios_crp3d.x\
			 trajtouTest_interpolationAt2dCuts_crp6d.x
# Rules 
trajtouJobs: trajtouSrc $(joblist) 
trajtouSrc:
	cd $(OOPTRAJTOUPATH)/src && make default

%.x: %.f90 libtrajtou.a libaotus.a
	$(LF) -o $(OOPTRAJTOUPATH)/bin/$@ $< $(OOPTRAJTOUPATH)/lib/libtrajtou.a $(OOPTRAJTOUPATH)/lib/libaotus.a

%.a:
	cd $(OOPTRAJTOUPATH)/src && make $@ 
clean:
	rm -f $(OOPTRAJTOUPATH)/mod/*.mod $(OOPTRAJTOUPATH)/lib/*.a $(OOPTRAJTOUPATH)/bin/*.x $(OOPTRAJTOUPATH)/bin/*.o
