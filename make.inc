# Compile options

F = gfortran
CF = $(F) -cpp -J mod -I mod -c -Wall -O2
LF = $(F) -cpp -I mod -Wall -O2

# Some paths to find files 
vpath %.f90 jobs
vpath %.a lib

# General rules

%.x: %.o libtrajtou.a
	$(LF) -o bin/$@ bin/$< lib/libtrajtou.a 

%.o: %.f90
	$(CF) $< -o bin/$@
