!##################################################
! PROGRAM: TEST_SURFACE
!> @brief
!! Program for testing purpose
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!##################################################
PROGRAM TEST_SURFACE
   USE SURFACE_MOD
   USE DEBUG_MOD
IMPLICIT NONE
TYPE(Surface):: surf
CHARACTER(LEN=30) :: filename
CALL SET_DEBUG_MODE(.TRUE.) 
filename="INsurface.inp"
CALL surf%INITIALIZE(filename)
CALL surf%PRINT_PATTERN(6,2,"XYZ")

END PROGRAM TEST_SURFACE
