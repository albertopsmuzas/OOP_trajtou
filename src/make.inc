# Compile options

F = gfortran
FFLAGS = -cpp -J ../mod -I ../mod -c
LFFLAGS = -I ../mod 

# Where to find compiled but not yet linked binaries
vpath %.o ../bin

# Some conditionals
ifdef DEBUG
	FFLAGS =  -DDEBUG $(FFLAGS)
endif

# General rules
%.o : %.f90
	$(F) $(FFLAGS) $< -o ../bin/$@
