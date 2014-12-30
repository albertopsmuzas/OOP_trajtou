!#############################################################
! MODULE: DYNAMICS_MOD
!
!> @brief
!! General module for general dynamics schemes
!
!#############################################################
MODULE DYNAMICS_MOD
USE UNITS_MOD
USE LINK_INTEGRATOR_MOD
USE LINK_PES_MOD
USE INICOND_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////
! TYPE: Dynamics
!> @brief
!! Stores common information for dynamics jobs
!
!> @param alias - Human friendly name
!> @param filename - Name fo the input file
!> @param kind - Kind of dynamics
!> @param delta_t - Initial time-step
!> @param max_t - Maximum time 
!> @param nfollow - Number of trajectories to be followed
!> @param followtraj - array with ID of the trajectories to follow
!> @param init - initial conditions
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Feb/2014, Jun/2014
!> @version 1.0, 2.0
!---------------------------------------------------------------
TYPE,ABSTRACT :: Dynamics
   CHARACTER(LEN=30) :: alias
   CHARACTER(LEN=30) :: filename
   CHARACTER(LEN=10) :: kind
   TYPE(Time) :: delta_t
   TYPE(Time) :: max_t 
   INTEGER :: nfollow
   INTEGER,DIMENSION(:),ALLOCATABLE :: followtraj
   CLASS(Integrator_o1),ALLOCATABLE :: thisintegrator
   CLASS(Inicond),ALLOCATABLE :: thisinicond
   CLASS(PES),ALLOCATABLE :: thispes
   CONTAINS
      PROCEDURE(INITIALIZE_DYNAMICS),DEFERRED,PUBLIC:: INITIALIZE
      PROCEDURE(RUN_DYNAMICS),DEFERRED,PUBLIC:: RUN
END TYPE Dynamics
ABSTRACT INTERFACE
   !###########################################################
   !# SUBROUTINE: INITIALIZE_DYNAMICS 
   !###########################################################
   !> @brief
   !! Reads all specifications for a dynamics job
   !-----------------------------------------------------------
   SUBROUTINE INITIALIZE_DYNAMICS(this,filename)
      IMPORT Dynamics
      CLASS(Dynamics),INTENT(OUT):: this
      CHARACTER(LEN=*),INTENT(IN):: filename
   END SUBROUTINE INITIALIZE_DYNAMICS
   !###########################################################
   !# SUBROUTINE: RUN_DYNAMICS 
   !###########################################################
   !> @brief
   !! Starts dynamics of trajectories indicated in @b thisinicond
   !-----------------------------------------------------------
   SUBROUTINE RUN_DYNAMICS(this)
      IMPORT Dynamics
      CLASS(Dynamics),INTENT(INOUT):: this
   END SUBROUTINE RUN_DYNAMICS
END INTERFACE
!/////////////////////////////////////////////////////////////
END MODULE DYNAMICS_MOD
