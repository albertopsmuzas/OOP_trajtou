##################################################
# MAKE_FILE FOR  OOP TRAJTOU #####################
##################################################
include make.inc
.PHONY : build
joblist = test_CRP.x test_inicond.x test_dynamics.x get_CRP_value.x test_CRP6D.x
hola="adios"
# Rules 
build: $(joblist)
libtrajtou.a: 
	cd src && $(MAKE) libtrajtou.a

# Clean rule
.PHONY : clean
clean:
	rm -f mod/*.mod lib/*.a bin/*.x bin/*.o
