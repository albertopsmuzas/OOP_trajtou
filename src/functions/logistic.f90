!#########################################################
! MODULE: LOGISTIC_FUNCTION
!> @brief
!! Compendium of routines and types to manage a logistic function
!##########################################################
MODULE LOGISTIC_FUNCTION_MOD
USE FUNCTION1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: Logistic_func
!> @brief
!! Function @f$ f(x)=\over{1}{1+e^{ax-b}} @f$, where @b a and @b are
!! parameters
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Function1d) :: Logistic_func
CONTAINS
   PROCEDURE,PUBLIC :: getvalue => getvalue_logistic
   PROCEDURE,PUBLIC :: getderiv => getderiv_logistic
END TYPE Logistic_func
!/////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# FUNCTION: getvalue_logistic 
!###########################################################
!> @brief 
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getvalue_logistic(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Logistic_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: a
   ! Run section
   a=this%getparam()
   getvalue_logistic=1.D0/(1.D0+dexp(a(1)*x-a(2)))
   RETURN
END FUNCTION getvalue_logistic
!###########################################################
!# FUNCTION: getderiv_logistic 
!###########################################################
!> @brief 
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_logistic(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Logistic_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: a
   ! Run section
   a=this%getparam()
   getderiv_logistic=-a(1)*dexp(a(1)*x-a(2))/((1.D0+dexp(a(1)*x-a(2)))**2.D0)
   RETURN
END FUNCTION getderiv_logistic

END MODULE LOGISTIC_FUNCTION_MOD! contains, body
