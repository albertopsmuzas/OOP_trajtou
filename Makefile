##################################################
# MAKE_FILE FOR  OOP TRAJTOU #####################
##################################################
.PHONY : clean
.PHONY : all
.PHONY : build
.PHONY : testing

build :
	cd src && $(MAKE) build
testing : 
	cd src && $(MAKE) testing
clean :
	rm -f mod/*.mod lib/*.a bin/*.x bin/*.o
