##################################################
# MAKE_FILE FOR  OOP TRAJTOU #####################
##################################################
include make.inc
.PHONY : build
joblist = test_inicond3dInput.x trajtouDynamics3D.x trajtouGetPot_crp6d.x\
			 CRP6D_gridstats.ZR.x molec2atom.x au2angst.x get_CRP6D_shift.x\
			 test_CRP6D_cuts.x get_CRP6D_smoothvalue.x get_CRP6D_rawvalue.x\
			 cheat_CRP6D_cartwheelinput.x cart2surf.x test_surface.x test_inicond6D.x \
			 test_dynamics6D.x drawtraj6d.x anlys_diff3d.x anlys_diff6d.x test_systemInput.x\
			 test_crp3dInput.x trajtouGetPot_crp3d.x test_crp6dInput.x
# Rules 
link_jobs2lib: $(joblist)
lib: 
	cd src && $(MAKE) libtrajtou.a

# Clean rule
.PHONY : clean
clean:
	rm -f mod/*.mod lib/*.a bin/*.x bin/*.o
