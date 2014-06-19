!#########################################################
! MODULE: ONE_FUNCTION
!> @brief
!! Compendium of routines and types to manage a ONE function
!##########################################################
MODULE ONE_FUNCTION_MOD
USE FUNCTION1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: ONE_func
!> @brief
!! Function @f$ f(x)=\over{1}{1+e^{ax-b}} @f$, where @b a and @b are
!! parameters
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Function1d) :: ONE_func
CONTAINS
   PROCEDURE,PUBLIC :: getvalue => getvalue_ONE
   PROCEDURE,PUBLIC :: getderiv => getderiv_ONE
END TYPE ONE_func
!/////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# FUNCTION: getvalue_ONE 
!###########################################################
!> @brief 
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getvalue_ONE(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(ONE_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   getvalue_ONE=1.D0
   RETURN
END FUNCTION getvalue_ONE
!###########################################################
!# FUNCTION: getderiv_ONE 
!###########################################################
!> @brief 
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_ONE(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(ONE_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   getderiv_ONE=0.D0
   RETURN
END FUNCTION getderiv_ONE

END MODULE ONE_FUNCTION_MOD! contains, body
