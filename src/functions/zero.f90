!#########################################################
! MODULE: ZERO_FUNCTION
!> @brief
!! Compendium of routines and types to manage a ZERO function
!##########################################################
MODULE ZERO_FUNCTION_MOD
USE FUNCTION1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: ZERO_func
!> @brief
!! Just zero for all "x" values
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2016
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Function1d) :: Zero_func
CONTAINS
   PROCEDURE,PUBLIC :: getvalue => getvalue_ZERO
   PROCEDURE,PUBLIC :: getderiv => getderiv_ZERO
END TYPE Zero_func
!/////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# FUNCTION: getvalue_ZERO 
!###########################################################
!> @brief 
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getvalue_ZERO(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(ZERO_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   getvalue_ZERO=0.D0
   RETURN
END FUNCTION getvalue_ZERO
!###########################################################
!# FUNCTION: getderiv_ZERO 
!###########################################################
!> @brief 
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_ZERO(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(ZERO_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   getderiv_ZERO=0.D0
   RETURN
END FUNCTION getderiv_ZERO

END MODULE ZERO_FUNCTION_MOD! contains, body
