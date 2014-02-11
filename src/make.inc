# Compile options

F = gfortran
CF = $(F) -cpp -J ../mod -I ../mod -c -Wall -fbacktrace -ffixed-line-length-none -Wtabs -W -fbounds-check -O2
LF = $(F) -I ../mod -Wall -fbounds-check -Wtabs -ffixed-line-length-none -fbacktrace -W -O2

# Where to find compiled but not yet linked binaries
vpath %.o ../bin

# General rules
%.o : %.f90
ifeq ($(debugmode),on)
	$(CF) -DDEBUG $< -o ../bin/$@
else
	$(CF) $< -o ../bin/$@
endif
