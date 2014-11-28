!#########################################################
! MODULE: EXPONENTIAL_FUNCTION
!> @brief
!! Compendium of routines and types to manage a exponential function
!##########################################################
MODULE EXPONENTIAL_FUNCTION_MOD
USE FUNCTION1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: exponential_func
!> @brief
!! Function @f$ f(x)=ae^{bx} @f$, where @b a and @b b are
!! parameters
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Function1d) :: exponential_func
CONTAINS
   PROCEDURE,PUBLIC :: getvalue => getvalue_exponential
   PROCEDURE,PUBLIC :: getderiv => getderiv_exponential
END TYPE exponential_func
!/////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# FUNCTION: getvalue_exponential 
!###########################################################
!> @brief 
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getvalue_exponential(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(exponential_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: a
   ! Run section
   a=this%getparam()
   getvalue_exponential=a(1)*dexp(a(2)*x)
   RETURN
END FUNCTION getvalue_exponential
!###########################################################
!# FUNCTION: getderiv_exponential 
!###########################################################
!> @brief 
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_exponential(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(exponential_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: a
   ! Run section
   a=this%getparam()
   getderiv_exponential=a(1)*a(2)*dexp(a(2)*x)
   RETURN
END FUNCTION getderiv_exponential

END MODULE EXPONENTIAL_FUNCTION_MOD! contains, body
