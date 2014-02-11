!########################################################
! MODULE: INICOND_MOD
!
!> @brief
!! This module provides general routines to define initial conditions
!! for a dynobject, which is a general entity prepared for
!! dynamics simulation
!########################################################
MODULE INICOND_MOD
IMPLICIT NONE
!////////////////////////////////////////////////////////
! TYPE: Dynobject
!> @brief
!! General object prepared to perform dynamics simulation
!
!> @param E - Energy of the object
!> @param init_r - Initial position
!> @param init_p - Initial momentum
!> @param turning_point - Classic turning point
!> @param r - Actual position of the object
!> @param p - Actial momentum of the object
!> @param stat - String that defines the state of this object
!> @param id - Integer ID number
!--------------------------------------------------------
TYPE ::  Dynobject
   INTEGER :: id 
   REAL*8 :: E
   REAL*8,DIMENSION(:),ALLOCATABLE :: init_r
   REAL*8,DIMENSION(:),ALLOCATABLE :: init_p
   REAL*8,DIMENSION(:),ALLOCATABLE :: turning_point
   REAL*8,DIMENSION(:),ALLOCATABLE :: r 
   REAL*8,DIMENSION(:),ALLOCATABLE :: p 
   CHARACTER(LEN=10) :: stat 
   CONTAINS
      PROCEDURE,PUBLIC :: INITIALIZE => INITIALIZE_DYNOBJECT
END TYPE Dynobject
!////////////////////////////////////////////////////////
! TYPE Inicond
!> @brief
!! Stores information to create general initial conditions
!
!> @param alias - Human friendly name
!> @param kind
!> @param input_file - Input file
!> @param output_file - Output file
!> @param ntrajs - Number of trajectories
!> @param nstart - Initial trajectory
!> @param seed - Allocatable integer array to feed random functions
!-------------------------------------------------------
TYPE :: Inicond
	CHARACTER(LEN=30) :: alias
	CHARACTER(LEN=5) :: kind
	CHARACTER(LEN=30) :: input_file
	CHARACTER(LEN=30) :: output_file
	INTEGER :: ntraj ! number of trajectories
	INTEGER :: nstart ! initial trajectory
	INTEGER,DIMENSION(:),ALLOCATABLE :: seed ! Seed for random number generation
END TYPE Inicond
!////////////////////////////////////////////////////////
CONTAINS
!########################################################
! SUBROUTINE: INITIALIZE_DYNOBJECT
!> @brief
!! Allocates and initializes some atributes of a
!! Dynobject type variable
!
!> @param[out] this - Dynobject type variable
!> @param[in] id - Integer ID number to set
!> @param[in] dimens - Number of dimensions to set up
!--------------------------------------------------------
SUBROUTINE INITIALIZE_DYNOBJECT(this,id,dimens)
   IMPLICIT NONE
   CLASS(Dynobject),INTENT(OUT) :: this
   INTEGER(KIND=4),INTENT(IN) :: id, dimens
   ALLOCATE(this%init_r(dimens))
   ALLOCATE(this%init_p(dimens))
   ALLOCATE(this%p(dimens))
   ALLOCATE(this%r(dimens))
   ALLOCATE(this%turning_point(dimens))
   this%id = id
   this%stat="Unknown"
   RETURN
END SUBROUTINE INITIALIZE_DYNOBJECT
END MODULE INICOND_MOD
