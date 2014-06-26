!########################################################
! MODULE: INICOND_MOD
!
!> @brief
!! This module provides general routines to define initial conditions
!! for a dynobject, which is a general entity prepared for
!! dynamics simulation
!########################################################
MODULE INICOND_MOD
   USE LINK_PES_MOD
   USE UNITS_MOD
   USE CONSTANTS_MOD
IMPLICIT NONE
!////////////////////////////////////////////////////////
! TYPE: Dynobject
!> @brief
!! General object prepared to perform dynamics simulation. It can be an atom, molecule, etc.
!
!> @param E - Energy of the object
!> @param init_r - Initial position
!> @param init_p - Initial momentum
!> @param turning_point - Classic turning point
!> @param r - Actual position of the object
!> @param p - Actial momentum of the object
!> @param stat - String that defines the state of this object
!> @param ireb - Number of bouncings in Z direction
!> @param ixyboun - Number of bouncings in XY direction
!> @param id - Integer ID number
!--------------------------------------------------------
TYPE,ABSTRACT ::  Dynobject
   REAL(KIND=8) :: E
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: init_r
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: init_p
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: turning_point
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: r 
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: p 
   INTEGER(KIND=4) :: ireb
   INTEGER(KIND=4) :: ixyboun
   CHARACTER(LEN=10) :: stat 
   CONTAINS
      PROCEDURE(INITIALIZE_DYNOBJECT),DEFERRED,PUBLIC :: INITIALIZE
END TYPE Dynobject
ABSTRACT INTERFACE
   !###########################################################
   !# SUBROUTINE: INITIALIZE_DYNOBJECT 
   !###########################################################
   SUBROUTINE INITIALIZE_DYNOBJECT(this)
      IMPORT Dynobject
      CLASS(Dynobject),INTENT(OUT) :: this
   END SUBROUTINE INITIALIZE_DYNOBJECT
END INTERFACE
!/////////////////////////////////////////////////////////////
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
!-------------------------------------------------------------
TYPE,ABSTRACT :: Inicond
   CHARACTER(LEN=30) :: alias
   CHARACTER(LEN=5) :: kind
   CHARACTER(LEN=30) :: input_file
   CHARACTER(LEN=30) :: output_file
   INTEGER :: ntraj ! number of trajectories
   INTEGER :: nstart ! initial trajectory
   INTEGER,DIMENSION(:),ALLOCATABLE :: seed ! Seed for random number generation
   CLASS(Dynobject),DIMENSION(:),ALLOCATABLE :: trajs
   CONTAINS
      PROCEDURE(INITIALIZE_INICOND),DEFERRED,PUBLIC :: INITIALIZE
      PROCEDURE(GENERATE_TRAJS_INICOND),DEFERRED,PUBLIC :: GENERATE_TRAJS
END TYPE Inicond
ABSTRACT INTERFACE
   !###########################################################
   !# SUBROUTINE: INITIALIZE_INICOND 
   !###########################################################
   !-----------------------------------------------------------
   SUBROUTINE INITIALIZE_INICOND(this,filename)
      IMPORT Inicond
      CLASS(Inicond),INTENT(OUT) :: this
      CHARACTER(LEN=*),INTENT(IN) :: filename
   END SUBROUTINE INITIALIZE_INICOND
   !###########################################################
   !# SUBROUTINE: GENERATE_TRAJS_INICOND 
   !###########################################################
   !-----------------------------------------------------------
   SUBROUTINE GENERATE_TRAJS_INICOND(this,thispes)
      IMPORT Inicond
      IMPORT PES
      CLASS(Inicond),INTENT(INOUT) :: this
      CLASS(PES),INTENT(IN) :: thispes
   END SUBROUTINE GENERATE_TRAJS_INICOND
END INTERFACE
!////////////////////////////////////////////////////////
END MODULE INICOND_MOD
