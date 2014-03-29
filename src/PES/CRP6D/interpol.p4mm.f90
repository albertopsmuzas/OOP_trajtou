!#########################################################
! MODULE: INTERPOL_WYCKOFF_P4MM
!> @brief
!! Provides tools to interpolate through wyckoff sites belonging to
!! p4mm wallpaper symmetry
!##########################################################
MODULE INTERPOL_WYCKOFF_P4MM_MOD
! Initial declarations
USE INTERPOL_WYCKOFF_GENERIC_MOD
USE FOURIER1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: Wyckoffp4mm
!> @brief
!! Subclass of Wyckoffsitio for generic p4mm symmetry.
!! in p4mm symmetry does not matter
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 20/03/2014 
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Wyckoffsitio) :: Wyckoffp4mm
END TYPE Wyckoffp4mm
!/////////////////////////////////////////////////////////////////
CONTAINS
! contains, body
END MODULE INTERPOL_WYCKOFF_P4MM_MOD
