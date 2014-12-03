!#########################################################
! MODULE: XEXPONENTIAL_FUNCTION
!> @brief
!! Compendium of routines and types to manage a xexponential function
!##########################################################
MODULE XEXPONENTIAL_FUNCTION_MOD
USE FUNCTION1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: xexponential_func
!> @brief
!! Function @f$ f(x)=axe^{bx} @f$, where @b a and @b b are
!! parameters
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Function1d) :: xexponential_func
CONTAINS
   PROCEDURE,PUBLIC :: getvalue => getvalue_xexponential
   PROCEDURE,PUBLIC :: getderiv => getderiv_xexponential
END TYPE xexponential_func
!/////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# FUNCTION: getvalue_xexponential 
!###########################################################
!> @brief 
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getvalue_xexponential(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Xexponential_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: a
   ! Run section
   a=this%getparam()
   getvalue_xexponential=a(1)*x*dexp(a(2)*x)
   RETURN
END FUNCTION getvalue_xexponential
!###########################################################
!# FUNCTION: getderiv_xexponential 
!###########################################################
!> @brief 
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_xexponential(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(xexponential_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: a
   ! Run section
   a=this%getparam()
   getderiv_xexponential=a(1)*(1.D0+a(2)*x)*dexp(a(2)*x)
   RETURN
END FUNCTION getderiv_xexponential

END MODULE XEXPONENTIAL_FUNCTION_MOD! contains, body
