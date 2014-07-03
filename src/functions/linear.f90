!#########################################################
! MODULE: LINEAR_FUNCTION
!> @brief
!! Compendium of routines and types to manage a LINEAR function
!##########################################################
MODULE LINEAR_FUNCTION_MOD
USE FUNCTION1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: LINEAR_func
!> @brief
!! Function @f$ f(x)=\over{1}{1+e^{ax-b}} @f$, where @b a and @b are
!! parameters
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Function1d) :: LINEAR_func
CONTAINS
   PROCEDURE,PUBLIC :: getvalue => getvalue_LINEAR
   PROCEDURE,PUBLIC :: getderiv => getderiv_LINEAR
END TYPE LINEAR_func
!/////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# FUNCTION: getvalue_LINEAR 
!###########################################################
!> @brief 
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getvalue_LINEAR(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(LINEAR_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: a
   ! Run section
   a=this%getparam()
   getvalue_LINEAR=a(1)*x+a(2)
   RETURN
END FUNCTION getvalue_LINEAR
!###########################################################
!# FUNCTION: getderiv_LINEAR 
!###########################################################
!> @brief 
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_LINEAR(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(LINEAR_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: a
   ! Run section
   a=this%getparam()
   getderiv_LINEAR=a(1)
   RETURN
END FUNCTION getderiv_LINEAR

END MODULE LINEAR_FUNCTION_MOD! contains, body
