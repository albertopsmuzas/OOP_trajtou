!#########################################################
! RAISON D'ÊTRE:
! - Implementation of those routines needed to define a CRP PES 
! UPDATES:
! - Created: 18/12/2013 ---> Alberto Muzas 
! FUNCTIONALITY:
! - Type definition 
! IDEAS FOR THE FUTURE:
! - None 
!##########################################################
MODULE CRP_MOD
   USE PES_MOD
   USE SURFACE_MOD
! Initial declarations
IMPLICIT NONE
!============================================================================
! Pair potential derived data 
! --------------------------------
TYPE, PUBLIC :: Pair_pot
   CHARACTER*30 :: filename ! Usually: intrep?.dat (?) is a 1 digit numbe
   CHARACTER*30 :: alias ! human-legible alias
	INTEGER(KIND=4) :: id ! Integer identification label
	INTEGER(KIND=4) :: n ! Number of points
   REAL(KIND=8)  :: vasint ! Potential far from the surface
   REAL(KIND=8)  :: dz1, dz2 ! 1st V-derivative in Z(1) and Z(n). Boundary for spline interpolation
	REAL(KIND=8) :: rumpling ! rumpling associated with this pair potential
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: z, v ! Collection of v(zi)
	TYPE(Point_3D), DIMENSION(:), ALLOCATABLE :: vecnet
	TYPE(Interpol_Z) :: interz
END TYPE Pair_pot
!===========================================================
! Type Sitio
TYPE, PUBLIC :: Sitio
   CHARACTER(LEN=30) :: filename !Usually: Sitio?.dat (?) is a 1 digit number
   CHARACTER(LEN=30) :: alias
	CHARACTER(LEN=10) :: units_z
	CHARACTER(LEN=10) :: units_v
	INTEGER :: n ! Number of points
	REAL(KIND=8) :: x, y ! Position on the XY surface plane
	REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: z, v
	REAL(KIND=8) :: dz1, dz2 ! 1st V-derivative in Z(1) and Z(n)
	REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dvdx, dvdy, dvdz
	TYPE(Interpol_Z) :: interz
END TYPE Sitio
!=============================================================================
! Type Compendium of sitios
!---------------------------
TYPE, PUBLIC :: Sitios
	INTEGER :: n ! Number of sitios
	TYPE(Sitio), DIMENSION(:), ALLOCATABLE :: site ! An array of sitios
END TYPE Sitios
!=============================================================================
! Type Compendium of pair potentials
!---------------------------
TYPE, PUBLIC :: Pair_pots
	INTEGER :: n ! Number of sitios
	TYPE(Pair_pot), DIMENSION(:), ALLOCATABLE :: pairpot ! An array of pair potentials
END TYPE Pair_pots
!====================================================================================================
! Type for CRP PES
!------------------
TYPE, EXTENDS(PES) :: CRP
   INTEGER(KIND=4) :: max_order
   TYPE(Pair_pots) :: all_pairpots
   TYPE(Sitios) :: all_sitios
   TYPE(Surface) :: surf
CONTAINS
   PROCEDURE, PUBLIC :: INITIALIZE => INIT_CRP_PES
END TYPE CRP
END MODULE CRP_MOD
