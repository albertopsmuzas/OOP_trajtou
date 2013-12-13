##################################################
# MAKE_FILE FOR  OOP TRAJTOU #####################
##################################################
include make.inc
.PHONY : clean
.PHONY : build
.PHONY : testing
.PHONY : jobs
joblist = test.x

# Rules 
jobs: $(joblist)
build: libtrajtou.a jobs
libtrajtou.a: 
	cd src && $(MAKE) build
testing: 
	cd src && $(MAKE) testing
clean:
	rm -f mod/*.mod lib/*.a bin/*.x bin/*.o
