!#########################################################
! MODULE: TRAJDRAW_MOD
!> @brief
!! Module to analyze OUTtrajxxxxx.out files in order to obtain
!! graphic representations
!
!> @todo
!! - Generalize for CRP6D and other formats
!##########################################################
MODULE TRAJDRAW_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: Drawtraj
!> @brief
!! All utils and information to create pgraphic representation of a
!! trajectory
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!----------------------------------------------------------------
TYPE :: Drawtraj
PRIVATE
   CHARACTER(LEN=3) :: format_out
   TYPE(Surface):: surf
   INTEGER(KIND=4) :: npattern
   INTEGER(KIND=4) :: nprojectile
CONTAINS
   ! Initialization block
   PROCEDURE,PUBLIC :: INITIALIZE => INITIALIZE_TRAJDRAW
   ! Tools block
   PROCEDURE,PUBLIC :: DRAW => DRAW_TRAJDRAW
END TYPE Drawtraj
CONTAINS
!###########################################################
!# SUBROUTINE: INITIALIZE_TRAJDRAW 
!###########################################################
!> @brief
!! brief description
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INITIALIZE_TRAJDRAW(this,format_out,filename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Drawtraj),INTENT(OUT):: this
   CHARACTER(LEN=3),INTENT(IN) :: format_out
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables
   ! Run section
   this%format_out=format_out
   CALL this%surf%INITIALIZE(filename)
   npattern=sum(this%atomtype(:)%n)
   nprojectile=1
   RETURN
END SUBROUTINE INITIALIZE_TRAJDRAW
END MODULE TRAJDRAW_MOD
