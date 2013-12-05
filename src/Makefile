# Include general rules
include make.inc

# Add here all folder names in src
folder_list= debug constants units PES surfaces

# Object list. There should be a rule for each object in this list
object_list= debug.o constants.o units.o 

# Where to find fortran files:
vpath %.f90 $(folder_list)

# DEPENDENCIES SECTION ###################################################
testing : $(object_list)
build : $(object_list)

debug.o : 
constants.o : 
units.o : debug.o constants.o
