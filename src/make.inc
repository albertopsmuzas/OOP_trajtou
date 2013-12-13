# Compile options

F = gfortran
CF = $(F) -cpp -J ../mod -I ../mod -c
LF = $(F) -I ../mod 

# Where to find compiled but not yet linked binaries
vpath %.o ../bin

# Some conditionals
ifdef DEBUG
	FFLAGS =  -DDEBUG $(FFLAGS)
endif

# General rules
%.o : %.f90
	$(CF) $< -o ../bin/$@
