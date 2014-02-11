# Compile options

F = gfortran
CF = $(F) -cpp -J mod -I mod -c -Wall
LF = $(F) -cpp -I mod -Wall

# Some paths to find files 
vpath %.f90 jobs
vpath %.a lib

# General rules

%.x: %.o libtrajtou.a
	$(LF) -o bin/$@ bin/$< lib/libtrajtou.a 

%.o: %.f90
	$(CF) $< -o bin/$@
