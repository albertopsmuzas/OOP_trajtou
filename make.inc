# Compile options

#F = ifort
#CF = $(F) -cpp -I mod -c -O2 -traceback -CB -gen-interfaces -no-wrap-margin -warn truncated-source
#LF = $(F) -cpp -I mod -O2
F = gfortran
CF = $(F) -cpp -J mod -I mod -c -Wall -Wtabs -O2
LF = $(F) -cpp -I mod -Wall -Wtabs -O2

# Some paths to find files 
vpath %.f90 jobs
vpath %.a lib

# General rules

%.x: %.o libtrajtou.a
	$(LF) -o bin/$@ bin/$< lib/libtrajtou.a 

%.o: %.f90
	$(CF) $< -o bin/$@
