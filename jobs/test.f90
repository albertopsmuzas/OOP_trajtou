PROGRAM TEST
USE PES_MOD
IMPLICIT NONE
TYPE(PES) :: crap
CALL crap%INITIALIZE()
WRITE(*,*) "Alias:", crap%get_alias()
WRITE(*,*) "Dimensions:", crap%get_dimensions()
WRITE(*,*) "Is initialized?: ", crap%is_initialized()
END PROGRAM TEST
