!#########################################################
! MODULE: CRP_MOD
!> @brief
!! Provides all those routines and types needed to define a CRP PES
!
!> @warning
!! - This module uses PES_MOD, SURFACE_MOD and INTERPOL1D_MOD
!
!> @see pes_mod, surface_mod, interpol1d_mod
!##########################################################
MODULE CRP_MOD
   USE PES_MOD
   USE SURFACE_MOD
   USE INTERPOL1D_MOD
! Initial declarations
IMPLICIT NONE
!///////////////////////////////////////////////////////////////////////////////
! TYPE: Symmetric point
! ---------------------
!
!> @brief 
!! For CRP PES. Symmetric point in XY plane. Composed by a scan in Z @f$(Z_{i},F(Z_{i}))@f$
!! variable, perpendicular to the surface.
!
!> @details
!! - Should be chosen inside the irreducible Wegner-Seitz cell of the surface
!> @param filename - File that contains initial information for this point
!> @param alias - Human-friendly alias for the Symmetric point
!> @param n - Number of points in the scan
!> @param x,y - Position in the XY plane
!> @param units_z,units_v - Units in which v(:) and z(:) are stored. Should be compatible
!!                          with units_mod module
!> @param dz1,dz2 - Second derivatives at Z(1) and Z(n)
!> @param interz - One dimensional interpolation for Z variable
!
!> @author A.P. Muzas
!> @version 1.0
!> @date 21/Jan/2014
!------------------------------------------------------------------------------
TYPE :: Symmpoint
   PRIVATE
   CHARACTER(LEN=30) :: filename
   CHARACTER(LEN=30) :: alias
   INTEGER(KIND=4) :: n
   REAL(KIND=8) :: x
   REAL(KIND=8) :: y
   CHARACTER(LEN=10) :: units_z
   CHARACTER(LEN=10) :: units_v
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: z
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: v
   REAL(KIND=8) :: dz1
   REAL(KIND=8) :: dz2
   TYPE(Interpol1d) :: interz 
   CONTAINS
      PROCEDURE,PUBLIC :: READ_RAW => READ_SYMMPOINT_RAW
END TYPE Symmpoint
!///////////////////////////////////////////////////////////////////////
! SUBTYPE: Pair potential  
! -----------------------
!> @brief 
!! For CRP PES. Symmetric point used to create a pair potential to extract corrugation.
!
!> @details
!! - It is needed one pair potentials for each different kind of atom in the surface
!! - Usually, the pair potentian is built on a top site of the surface.
!
!> @param id - Integer number to identify the pair potential
!> @param vasint - Potential far from the surface
!> @param rumpling - Rumpling associated with this pair potential. As pair potentials
!!                      are built on top sites and are related to different kinds of atoms
!!                      in the surface, the position of the coulomb repulsive singularity
!!                      may not be situated at Z=0. This rumpling is the shift so that
!!                      the singularity is plazed at zero.
!> @param vecnet - Collection of neightbour atoms
!
!> @author A.P. Muzas
!> @version 1.0
!> @date 21/Jan/2014
!------------------------------------------------------------------------------
TYPE,EXTENDS(Symmpoint) :: Pair_pot
   PRIVATE
   INTEGER(KIND=4) :: id 
   REAL(KIND=8) :: vasint
	REAL(KIND=8) :: rumpling
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: vecnet
END TYPE Pair_pot
!/////////////////////////////////////////////////////////////////////
! SUBTYPE: Sitio
! --------------
!> @brief 
!! For CRP PES. Symmetric point inside the irreducible Wigner-Seitz cell of the surface.
!! Used just to build up the PES.
!
!> @param dvdx,dvdy,dvdz - Collection of derivatives for each @f$V(Z_{i})@f$
!
!> @author A.P. Muzas
!> @version 1.0
!> @date 21/Jan/2014
!------------------------------------------------------------------------------
TYPE,EXTENDS(Symmpoint) :: Sitio
   PRIVATE
	REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dvdx, dvdy, dvdz 
END TYPE Sitio
!/////////////////////////////////////////////////////////////////////
! TYPE: Compendium of sitios
! --------------------------
!> @brief 
!! For CRP PES. Just a collection of Sitio types
!
!> @param n - Number of sitios
!> @param site - An array of sitios
!
!> @author A.P. Muzas
!> @version 1.0
!> @date 21/Jan/2014
!------------------------------------------------------------------------------
TYPE :: Sitios
   PRIVATE
   INTEGER :: n
   TYPE(Sitio), DIMENSION(:), ALLOCATABLE :: site
END TYPE Sitios
!/////////////////////////////////////////////////////////////////////
! TYPE: Compendium of pair potentials
! -----------------------------------
!> @brief
!! For CRP PES. Just a collection of Sitio types
!
!> @param n - Number of pair potentials
!> @param pairpot - An array of pair potentials
!
!> @author A.P. Muzas
!> @version 1.0
!> @date 21/Jan/2014
!------------------------------------------------------------------------------
TYPE :: Pair_pots
   PRIVATE
   INTEGER :: n
   TYPE(Pair_pot), DIMENSION(:), ALLOCATABLE :: pairpot
END TYPE Pair_pots
!////////////////////////////////////////////////////////////////////////////////
! SUBTYPE: CRP
! ------------
!> @brief 
!! Contains all information that can be extracted or needed during a CRP interpolation
!
!> @details
!! - This is a subtype PES. In principle, should be compatible with all available 
!!   dynamics jobs.
!
!> @param max_order - Integer number. Order of the environment used to extract corrugation
!> @param all_pairpots,all_sitios - Compendium of all pair potentials and sitios used
!!                                  in the interpolation procedure
!> @param surf - Surface parameters in which the CRP surface is defined
!
!> @author A.P. Muzas
!> @version 1.0
!> @date 21/Jan/2014
!------------------------------------------------------------------------------
TYPE, EXTENDS(PES) :: CRP
   PRIVATE
   INTEGER(KIND=4) :: max_order
   TYPE(Pair_pots) :: all_pairpots
   TYPE(Sitios) :: all_sitios
   TYPE(Surface) :: surf
CONTAINS
   PROCEDURE, PUBLIC :: INITIALIZE => INIT_CRP_PES
END TYPE CRP
!///////////////////////////////////////////////////////////////////////////
CONTAINS
!######################################################################
! SUBROUTINE: READ_SYMMPOINT_RAW ######################################
!######################################################################
!> @brief 
!! Initializes all data contained in a generic Symmpoint type variable 
!! from a "raw" input file.
!
!> @details
!! - Data is stored in a.u. automatically despite the units in input file.
!
!> @param[in,out] symmraw - A Symmpoint type variable. All its members will be loaded
!!                          from the input file
!> @param[in] filename - File that contains the "raw" input file
!
!> @author A.P. Muzas
!> @version 1.0
!> @date 22/Jan/2014
!
!> @warning
!! - This routine uses UNITS_MOD and maybe DEBUG_MOD
!! - This is the @b format needed for the "raw" input file:
!!    -# line 1: @b real(kind=8),@b real(kind=8),@b character(len=10) <-- X location, Y location, length units
!!    -# line 2: @b integer(kind=4) <-- Number of points (N) contained in the file
!!    -# line 3: @b character(len=10),@b character(len=10) <-- unit codes for Z(length) and potential(energy)
!!    -# lines 4~4+N: @b real(kind=8),@b real(kind=8) <-- cuples of values @f$Z_{i}, F(Z_{i})@f$. They should
!!                    be consistent with the previous choice of units
!
!> @see units_mod, debug_mod
!---------------------------------------------------------------------- 
SUBROUTINE READ_SYMMPOINT_RAW(symmraw,filename)
   USE UNITS_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
	IMPLICIT NONE
	! I/O Variables ----------------
	CLASS(Symmpoint),INTENT(INOUT) :: symmraw
   CHARACTER(LEN=*),INTENT(IN) :: filename
	! Local variables
	INTEGER :: i ! Counter
   REAL(KIND=8) :: aux_r1, aux_r2
	CHARACTER(LEN=20), PARAMETER :: routinename = "READ_SYMMPOINT_RAW: "
	TYPE(Length) :: x, y, z
   TYPE(Energy) :: v
	CHARACTER(LEN=10) :: units1,units2
	! FIRE IN THE HOLE!
   symmraw%filename=filename
	OPEN(10,FILE=symmraw%filename,STATUS="OLD")
		READ(10,*) aux_r1,aux_r2,units1
      CALL x%READ(aux_r1,units1)
      CALL y%READ(aux_r2,units1)
		symmraw%x = x%getvalue()
		symmraw%y = y%getvalue()
		! ------
		READ(10,*) symmraw%n
		READ(10,*) units1, units2
		ALLOCATE(symmraw%v(symmraw%n))
		ALLOCATE(symmraw%z(symmraw%n))
		DO i=1,symmraw%n
         READ(10,*) aux_r1,aux_r2
         CALL z%READ(aux_r1,units1)
         CALL v%READ(aux_r2,units2)
         CALL z%TO_STD() ! go to standard units (a.u.)
         CALL v%TO_STD() ! go to standard units (a.u.)
			symmraw%z(i) = z%getvalue()
			symmraw%v(i) = v%getvalue()
		END DO
		symmraw%units_z = z%getunits()
		symmraw%units_v = v%getunits()
	CLOSE(10)
	! Write status message
#ifdef DEBUG
	CALL VERBOSE_WRITE(routinename, symmraw%filename, "Done")
#endif
END SUBROUTINE READ_SYMMPOINT_RAW
!###########################################################
!# SUBROUTINE: READ_STANDARD_PAIRPOT #######################
!###########################################################
!> @brief
!! Initializes all data for a Pair_pot subtype variable from
!! an standard input file, i.e. one written in a.u.
!
!> @param[in,out] pairpot - Pair_pot subtype variable to be initialized
!> @param[out] filename - String with the name of the standard pairpot
!!                        file to be loaded
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 24/Jan/2014
!> @version 1.0
!
!> @warning
!! - There is a routine to create standard input files from "raw" ones.
!-----------------------------------------------------------
SUBROUTINE READ_PAIRPOT (pairpot,filename)
      ! Initial declarations
      IMPLICIT NONE
      ! I/O variables ------------------------------
      TYPE(Pair_pot), INTENT(INOUT) :: pairpot
      CHARACTER(LEN=*),INTENT(IN) :: filename
      ! Local variables ----------------------------
      INTEGER :: i
      CHARACTER(LEN=14), PARAMETER :: routinename = "READ_PAIRPOT: "
      ! Run section ---------------------------------
      OPEN(10,FILE=pairpot%filename,STATUS='OLD')
      READ(10,*) ! dummy line
      READ(10,*) ! dummy line
      READ(10,*) pairpot%vasint
      READ(10,*) pairpot%dz1
      READ(10,*) pairpot%dz2
      READ(10,*) pairpot%id,pairpot%rumpling
      READ(10,*) pairpot%n
      ALLOCATE(pairpot%z(1:pairpot%n))
      ALLOCATE(pairpot%v(1:pairpot%n))
      DO i=1, pairpot%n
         READ(10,*) pairpot%z(i), pairpot%v(i)
      END DO
      CLOSE(10)
      ! Status message
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename, pairpot%filename, "Done")
#endif
END SUBROUTINE READ_PAIRPOT

END MODULE CRP_MOD
