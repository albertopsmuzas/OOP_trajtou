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
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: param
CONTAINS
   ! Initialization block
   PROCEDURE,PUBLIC,NON_OVERRIDABLE :: READ => READ_FUNCTION1D
   ! Set block
   PROCEDURE,PUBLIC,NON_OVERRIDABLE :: SET_ID => SET_ID_FUNCTION1D
   ! Get block
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: getparam => getparam_function1d
   PROCEDURE(getvalue_function1d),PUBLIC,DEFERRED :: getvalue
   PROCEDURE(getvalue_function1d),PUBLIC,DEFERRED :: getderiv
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
END INTERFACE

!//////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: READ_FUNCTION1D
!###########################################################
!> @brief
!! Read parameters.
!-----------------------------------------------------------
SUBROUTINE READ_FUNCTION1D(this,param)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Function1d),INTENT(INOUT):: this
   REAL(KIND=8),DIMENSION(:) :: param
   ! Run section
   ALLOCATE(this%param(size(param)))
   this%param=param
   RETURN
END SUBROUTINE READ_FUNCTION1D
!###########################################################
!# FUNCTION: getparam_function1d 
!###########################################################
!> @brief 
!! Common get function. Gets atribute "param"
!-----------------------------------------------------------
FUNCTION getparam_function1d(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Function1d),INTENT(IN):: this
   ! Local variables
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: getparam_function1d
   ! Run section
   ALLOCATE(getparam_function1d(size(this%param)))
   getparam_function1d=this%param
   RETURN
END FUNCTION getparam_function1d
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
