!##################################################
! PROGRAM: au2angst
!> @brief
!! Just a tool to improve workflow
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!##################################################
PROGRAM convert
! Initial declarations
USE UNITS_MOD
IMPLICIT NONE
! Variables
TYPE(Length) :: len
REAL(KIND=8) :: input
! variables
READ(*,*) input
CALL len%READ(input,"au")
CALL len%TO_ANGST()
WRITE(*,*) len%getvalue()
END PROGRAM convert
