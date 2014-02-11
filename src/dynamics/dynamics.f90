!#############################################################
! MODULE: DYNAMICS_MOD
!
!> @brief
!! General module for general dynamics schemes
!
!#############################################################
MODULE DYNAMICS_MOD
USE UNITS_MOD
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
!> @param follow_dynobj - array with ID of the trajectories to follow
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 10/Feb/2014
!> @version 1.0
!---------------------------------------------------------------
TYPE :: Dynamics
   CHARACTER(LEN=30) :: alias
   CHARACTER(LEN=30) :: filename
   CHARACTER(LEN=10) :: kind
   TYPE(Time) :: delta_t ! initial time-step
   TYPE(Time) :: max_t  ! upper limit in time for all trajectori
   INTEGER :: nfollow ! number of trajectories to follow
   INTEGER, DIMENSION(:), ALLOCATABLE :: follow_atom ! atoms to be followed (store id)
END TYPE Dynamics
!/////////////////////////////////////////////////////////////
END MODULE DYNAMICS_MOD
