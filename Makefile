##################################################
# MAKE_FILE FOR  OOP TRAJTOU #####################
##################################################
include make.inc
.PHONY : build
joblist = test.x
hola="adios"
# Rules 
build: $(joblist)
libtrajtou.a: 
	cd src && $(MAKE) libtrajtou.a

# Clean rule
.PHONY : clean
clean:
	rm -f mod/*.mod lib/*.a bin/*.x bin/*.o
