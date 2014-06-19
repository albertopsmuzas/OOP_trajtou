!#########################################################
! MODULE: FUNCTION1D_MOD
!> @brief
!! Types and procedures to manage a function 
!##########################################################
MODULE FUNCTION1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: Function1d
!> @brief
!! General type for functions
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!----------------------------------------------------------------
TYPE,ABSTRACT :: Function1d
   PRIVATE
   CHARACTER(LEN=10) :: id
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: parameters
CONTAINS
   PROCEDURE,PUBLIC,NON_OVERRIDABLE :: SET_ID => SET_ID_FUNCTION1D
   PROCEDURE(getvalue_function1d),PUBLIC,DEFERRED :: getvalue
   PROCEDURE(getvalue_function1d),PUBLIC,DEFERRED :: getderiv
   PROCEDURE(READ_FUNCTION1D),PUBLIC,DEFERRED :: READ
END TYPE Function1d

ABSTRACT INTERFACE
   !###########################################################
   !# SUBROUTINE: getvalue_function1d 
   !###########################################################
   !-----------------------------------------------------------
   REAL(KIND=8) FUNCTION getvalue_function1d(this,x)
      IMPORT Function1d
      CLASS(Function1d),INTENT(IN)::this
      REAL(KIND=8),INTENT(IN) :: x
   END FUNCTION getvalue_function1d
   !###########################################################
   !# SUBROUTINE: READ_FUNCTION1D 
   !###########################################################
   !-----------------------------------------------------------
   SUBROUTINE READ_FUNCTION1D(this,coeff)
      IMPORT Function1d
      CLASS(Function1d),INTENT(INOUT):: this
      REAL(KIND=8),DIMENSION(:) :: coeff
   END SUBROUTINE READ_FUNCTION1D
END INTERFACE

!//////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_ID_FUNCTION1D 
!###########################################################
!> @brief
!! Common set function. Sets atribute "id"
!-----------------------------------------------------------
SUBROUTINE SET_ID_FUNCTION1D(this,id)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Function1d),INTENT(INOUT):: this
   CHARACTER(LEN=*),INTENT(IN) :: id
   ! Run section
   this%id=id
   RETURN
END SUBROUTINE SET_ID_FUNCTION1D
END MODULE FUNCTION1D_MOD 
