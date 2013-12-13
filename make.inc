# Compile options

F = gfortran
CF = $(F) -J mod -I mod -c
LF = $(F) -I mod 
Lapack = /usr/lib/liblapack.so.3 

# Some paths to find files 
vpath %.o bin
vpath %.f90 jobs
vpath %.a lib

# General rules
%.o : %.f90
	$(CF) $< -o bin/$@
%.x : %.f90
	$(LF) -o bin/@ $(Lapack) lib/libtrajtou.a $< 
