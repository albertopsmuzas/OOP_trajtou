##################################################
# MAKE_FILE FOR  OOP TRAJTOU #####################
##################################################
include make.inc
.PHONY : build
joblist = test_CRP3D.x test_inicond3D.x test_dynamics3D.x get_CRP3D_value.x test_CRP6D.x\
			 CRP6D_gridstats.ZR.x get_CRP6D_value.x molec2atom.x au2angst.x get_CRP6D_shift.x\
			 test_CRP6D_cuts.x get_CRP6D_smoothvalue.x get_CRP6D_rawvalue.x test_CRP6D_bugs.x\
			 cheat_CRP6D_cartwheelinput.x cart2surf.x test_surface.x test_inicond6D.x \
			 test_dynamics6D.x drawtraj6d.x anlys_diff3d.x anlys_diff6d.x
# Rules 
build: $(joblist)
libtrajtou.a: 
	cd src && $(MAKE) libtrajtou.a

# Clean rule
.PHONY : clean
clean:
	rm -f mod/*.mod lib/*.a bin/*.x bin/*.o
