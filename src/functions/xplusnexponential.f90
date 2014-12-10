!#########################################################
! MODULE: XPLUSNEXPONENTIAL_FUNCTION
!> @brief
!! Compendium of routines and types to manage a xplusnexponential function
!##########################################################
MODULE XPLUSNEXPONENTIAL_FUNCTION_MOD
USE FUNCTION1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: xplusnexponential_func
!> @brief
!! Function @f$ f(x)=a(x+n)e^{bx} @f$, where @b a, @b b and @b n are
!! parameters
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Function1d) :: xplusnexponential_func
CONTAINS
   PROCEDURE,PUBLIC :: getvalue => getvalue_xplusnexponential
   PROCEDURE,PUBLIC :: getderiv => getderiv_xplusnexponential
END TYPE xplusnexponential_func
!/////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# FUNCTION: getvalue_xplusnexponential 
!###########################################################
!> @brief 
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getvalue_xplusnexponential(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(xplusnexponential_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),DIMENSION(3) :: a
   ! Run section
   a=this%getparam()
   getvalue_xplusnexponential=a(1)*(x+a(2))*dexp(a(3)*x)
   RETURN
END FUNCTION getvalue_xplusnexponential
!###########################################################
!# FUNCTION: getderiv_xplusnexponential 
!###########################################################
!> @brief 
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_xplusnexponential(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(xplusnexponential_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),DIMENSION(3) :: a
   ! Run section
   a=this%getparam()
   getderiv_xplusnexponential=a(1)*(1.D0+a(2)*a(3)+a(3)*x)*dexp(a(3)*x)
   RETURN
END FUNCTION getderiv_xplusnexponential

END MODULE XPLUSNEXPONENTIAL_FUNCTION_MOD! contains, body
