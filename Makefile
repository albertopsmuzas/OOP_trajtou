##################################################
# MAKE_FILE FOR  OOP TRAJTOU #####################
##################################################
include make.inc
.PHONY : build
joblist = test_CRP3D.x test_3Dinicond.x test_3Ddynamics.x get_CRP3D_value.x test_CRP6D.x\
			 CRP6D_gridstats.ZR.x
hola="adios"
# Rules 
build: $(joblist)
libtrajtou.a: 
	cd src && $(MAKE) libtrajtou.a

# Clean rule
.PHONY : clean
clean:
	rm -f mod/*.mod lib/*.a bin/*.x bin/*.o
