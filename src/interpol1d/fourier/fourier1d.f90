!########################################################
! MODULE : FOURIER1D
!
!> @brief
!! Provides tools to perform 1D periodical interpolations
!
!> @warning
!! - Includes Interpol1d_mod in its scope
!########################################################
MODULE FOURIER1D_MOD
USE INTERPOL1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////
! TYPE: FOURIER1D
!> @brief
!! Class to store all information needed for a 1D REAL fourier inteprolation
!------------------------------------------------
TYPE,EXTENDS(Interpol1d):: Fourier1d
   PRIVATE
   LOGICAL :: par
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: coeffpar,coeffodd
   CONTAINS
      PROCEDURE,PUBLIC :: ispar => ispar_fourier1d
END TYPE Fourier1d
!//////////////////////////////////////////////////
CONTAINS
!##################################################
! FUNCTION: ispar_fourier1d
!> @brief
!! Basic enquire function. If the function is par,
!! sinus elements vanish from expansion
!--------------------------------------------------
LOGICAL FUNCTION ispar_fourier1d(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d),INTENT(IN) :: this
   ! Run section ----------------
   ispar_fourier1d = this%par
   RETURN
END FUNCTION ispar_fourier1d
!#################################################
! SUBROUTINE: SET_COEFF_FURIER1D
!> @brief
!! Performs interpolation
!-------------------------------------------------
SUBROUTINE SET_COEFF_FOURIER1D(this,
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
END SUBROUTINE SET_COEFF_FOURIER1D
END MODULE FOURIER1D_MOD
