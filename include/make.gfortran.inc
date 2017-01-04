export F  = gfortran
export CC = gcc
export CF = $(F) -cpp -J $(OOPTRAJTOUPATH)/mod -I $(OOPTRAJTOUPATH)/mod -c -Wall -Wtabs -O2 -ffree-line-length-none -g
export LF = $(F) -cpp -I $(OOPTRAJTOUPATH)/mod -Wall -Wtabs -O2 -ffree-line-length-none -g
ifeq ($(OOPTRAJTOU_DEBUGMODE),on)
	CF = $(CF) -fcheck=all
	LF = $(LF) -fcheck=all
endif
