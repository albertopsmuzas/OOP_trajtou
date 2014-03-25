!#########################################################
! MODULE: CRP3D_MOD
!> @brief
!! Provides all those routines and types needed to define a CRP3D PES
!
!> @warning
!! - This module uses PES_MOD, SURFACE_MOD and INTERPOL1D_MOD
!
!> @see pes_mod, surface_mod, interpol1d_mod
!##########################################################
MODULE CRP3D_MOD
   USE PES_MOD
   USE SURFACE_MOD
   USE CUBICSPLINES_MOD
   USE FOURIER2D_MOD
! Initial declarations
IMPLICIT NONE
!///////////////////////////////////////////////////////////////////////////////
! TYPE: Symmetric point
! ---------------------
!
!> @brief 
!! For CRP3D PES. Symmetric point in XY plane. Composed by a scan in Z @f$(Z_{i},F(Z_{i}))@f$
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
   TYPE(Csplines),PUBLIC :: interz 
   CONTAINS
      PROCEDURE,PUBLIC :: READ_RAW => READ_SYMMPOINT_RAW
      PROCEDURE,PUBLIC :: GEN_SYMM_RAW => GEN_SYMMETRIZED_RAW_INPUT
      PROCEDURE,PUBLIC :: PLOT_DATA => PLOT_DATA_SYMMPOINT
END TYPE Symmpoint
!///////////////////////////////////////////////////////////////////////
! SUBTYPE: Pair potential  
! -----------------------
!> @brief 
!! For CRP3D PES. Symmetric point used to create a pair potential to extract corrugation.
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
   CONTAINS
      PROCEDURE,PUBLIC :: READ => READ_STANDARD_PAIRPOT
END TYPE Pair_pot
!/////////////////////////////////////////////////////////////////////
! SUBTYPE: Sitio
! --------------
!> @brief 
!! For CRP3D PES. Symmetric point inside the irreducible Wigner-Seitz cell of the surface.
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
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: dvdx, dvdy, dvdz 
   CONTAINS
      PROCEDURE,PUBLIC :: READ => READ_STANDARD_SITIO
END TYPE Sitio
!////////////////////////////////////////////////////////////////////////////////
! SUBTYPE: CRP3D
! ------------
!> @brief 
!! Contains all information that can be extracted or needed during a CRP3D interpolation
!
!> @details
!! - This is a subtype PES. In principle, should be compatible with all available 
!!   dynamics jobs.
!
!> @param max_order - Integer number. Order of the environment used to extract corrugation
!> @param all_pairpots,all_sitios - Compendium of all pair potentials and sitios used
!!                                  in the interpolation procedure
!> @param surf - Surface parameters in which the CRP3D surface is defined
!
!> @warning
!! - Use this type variable and related procedures once you've obtained a reliable set of input
!!   files to read. In order to create such a set, use Newinput variables and member procedures
!!
!> @author A.P. Muzas
!> @version 1.1
!> @date 04/Feb/2014
!
!> @see newinput
!------------------------------------------------------------------------------
TYPE, EXTENDS(PES) :: CRP3D
   INTEGER(KIND=4) :: max_order
   TYPE(Pair_pot),DIMENSION(:),ALLOCATABLE :: all_pairpots
   TYPE(Sitio),DIMENSION(:),ALLOCATABLE :: all_sites
   TYPE(Surface) :: surf
   CONTAINS
      PROCEDURE,PUBLIC :: READ => READ_CRP3D
      PROCEDURE,PUBLIC :: EXTRACT_VASINT => EXTRACT_VASINT_CRP3D
      PROCEDURE,PUBLIC :: SMOOTH => SMOOTH_CRP3D
      PROCEDURE,PUBLIC :: INTERPOL_Z => INTERPOL_Z_CRP3D
      PROCEDURE,PUBLIC :: GET_V_AND_DERIVS => GET_V_AND_DERIVS_CRP3D
      PROCEDURE,PUBLIC :: PLOT_XYMAP => PLOT_XYMAP_CRP3D
      PROCEDURE,PUBLIC :: PLOT_DIRECTION1D => PLOT_DIRECTION1D_CRP3D
      PROCEDURE,PUBLIC :: is_allowed => is_allowed_CRP3D
END TYPE CRP3D
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
!! - This routine uses UNITS_MOD, MATHS_MOD and maybe DEBUG_MOD
!! - This is the @b format needed for the "raw" input file:
!!    -# line 1: @b real(kind=8),@b real(kind=8),@b character(len=10); X location, Y location, length units
!!    -# line 2: @b integer(kind=4); Number of points (N) contained in the file
!!    -# line 3: @b character(len=10),@b character(len=10); unit codes for Z(length) and potential(energy)
!!    -# lines 4~4+N: @b real(kind=8),@b real(kind=8); cuples @f$(Z_{i},V(Z_{i}))@f$. They should
!!                    be consistent with the previous choice of units. There is no need of a specific ordering
!!                    of the data
!
!> @see units_mod, debug_mod, maths_mod
!---------------------------------------------------------------------- 
SUBROUTINE READ_SYMMPOINT_RAW(symmraw,filename)
   USE UNITS_MOD
   USE MATHS_MOD, ONLY: ORDER
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
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: aux_v1,aux_v2
   CHARACTER(LEN=20), PARAMETER :: routinename = "READ_SYMMPOINT_RAW: "
   TYPE(Length) :: x,y,z
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
      ALLOCATE(aux_v1(symmraw%n))
      ALLOCATE(aux_v2(symmraw%n))
      DO i=1,symmraw%n
         READ(10,*) aux_r1,aux_r2
         CALL z%READ(aux_r1,units1)
         CALL v%READ(aux_r2,units2)
         CALL z%TO_STD() ! go to standard units (a.u.)
         CALL v%TO_STD() ! go to standard units (a.u.)
         aux_v1(i) = z%getvalue()
         aux_v2(i) = v%getvalue()
      END DO
   CLOSE(10)
   CALL ORDER(aux_v1,aux_v2) ! Should order Z(i) and V(i)
   symmraw%z = aux_v1
   symmraw%v = aux_v2
   symmraw%units_z = z%getunits()
   symmraw%units_v = v%getunits()
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
!> @param[in] filename - String with the name of the standard pairpot
!!                        file to be loaded
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 24/Jan/2014
!> @version 1.0
!
!> @warning
!! - There is a routine to create standard input files from "raw" ones.
!! - May use DEBUG_MOD
!! - This is the @b format needed for a standard pair potential input:
!!    -# line 1: Dummy line, just for some comments
!!    -# line 2: Dummy line, just for some comments
!!    -# line 3: @b character(len=30); human-frienfly alias for this pair potential
!!    -# line 4: @b real(kind=8); Energy at long dinstance from the surface
!!    -# line 5: @b real(kind=8); @f$\partial^{2}{V(Z_{1})}\over{\partial{Z}^{2}}@f$
!!    -# line 6: @b real(kind=8); @f$\partial^{2}{V(Z_{N})}\over{\partial{Z}^{2}}@f$
!!    -# line 7: @b integer, @b real(kind=8); integer ID number for the pair potential, Z value
!!               in which @f$V(Z)=\pm\infty@f$
!!    -# line 8: @b integer; Number of points (N)
!!    -# line 9~9+N: @b real(kind=8), @b real(kind=8); couples @f$(Z_{i},V(Z_{i})), Z_{i}<Z_{i+1}@f$ in a.u.
!
!> @see 
!! debug_mod
!-----------------------------------------------------------
SUBROUTINE READ_STANDARD_PAIRPOT (pairpot,filename)
      ! Initial declarations
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   ! I/O variables ------------------------------
   CLASS(Pair_pot), INTENT(INOUT) :: pairpot
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables ----------------------------
   INTEGER :: i
   CHARACTER(LEN=14), PARAMETER :: routinename = "READ_PAIRPOT: "
   ! Run section ---------------------------------
   pairpot%filename=filename
   OPEN(10,FILE=pairpot%filename,STATUS='OLD')
   READ(10,*) ! dummy line
   READ(10,*) ! dummy line
   READ(10,*) pairpot%alias
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
   CALL pairpot%interz%READ(pairpot%z,pairpot%v)
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,pairpot%filename,pairpot%alias)
#endif
   RETURN
END SUBROUTINE READ_STANDARD_PAIRPOT
!###########################################################
!# SUBROUTINE: READ_STANDARD_SITIO #########################
!###########################################################
!> @brief
!! Initializes all data for a Sitio subtype variable from
!! an standard input file, i.e. one written in a.u.
!
!> @param[in,out] site - Sitio subtype variable to be initialized
!> @param[in] filename - String with the name of the standard site
!!                        file to be loaded
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 24/Jan/2014
!> @version 1.0
!
!> @warning
!! - There is a routine to create standard input files from "raw" ones.
!! - May use DEBUG_MOD
!! - This is the @b format needed for a standard Sitio input:
!!    -# line 1: Dummy line, just for some comments
!!    -# line 2: Dummy line, just for some comments
!!    -# line 3; @b character(len=30); human_friendly alias for this site
!!    -# line 4: @b real(kind=8),@b real(kind=8); X, Y location in a.u.
!!    -# line 5: @b integer; Number of points (N)
!!    -# line 6: @b real(kind=8); @f$\partial{V(Z_{1})}\over{\partial{Z}}@f$
!!    -# line 7: @b real(kind=8); @f$\partial{V(Z_{N})}\over{\partial{Z}}@f$
!!    -# line 8~8+N: @b real(kind=8),@b real(kind=8); couples @f$(Z_{i},V(Z_{i})), Z_{i}<Z_{i+1}@f$ in a.u.
!
!-----------------------------------------------------------
SUBROUTINE READ_STANDARD_SITIO(site,filename)
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   ! I/O variables
   CLASS(Sitio),INTENT(INOUT) :: site
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables
   INTEGER :: i ! counter
   CHARACTER(LEN=12) , PARAMETER :: routinename = "READ_SITIO: "
   !
   site%filename=filename
   OPEN (10,file=site%filename,status="old")
   READ(10,*) ! dummy line
   READ(10,*) ! dummy line
   READ(10,*) site%alias
   READ(10,*) site%x, site%y
   READ(10,*) site%n
   READ(10,*) site%dz1
   READ(10,*) site%dz2
   ALLOCATE(site%z(1:site%n))
   ALLOCATE(site%v(1:site%n))
   DO i=1, site%n
      READ(10,*) site%z(i), site%v(i)
   END DO
   CLOSE(10)
   CALL site%interz%READ(site%z,site%v)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,site%filename,site%alias)
#endif
   RETURN
END SUBROUTINE READ_STANDARD_SITIO
!####################################################################
! SUBROUTINE:  GEN_SYMMETRIZED_RAW_INPUT ############################
!####################################################################
!> @brief
!! Tool that adds new points to a "raw" input (can be other kind of input)
!! so that the function is symmetric respect to @b zero. These changes
!! are made only for those couples that satisfy @f$V(z_{i})> f_{0}@f$. Changes
!! are stored in a file called @b filename, which is a standard "raw" input.
!
!> @param[in,out] symmraw - Symmpoint variable to be symmetrized
!> @param[in,out] zero - Length quantity, contains @f$x_{0}@f$. After execution, in a.u.
!> @param[in,out] vtop - Energy quantity, contains @f$f_{0}@f$. After execution, in a.u.
!> @param[in] filename - Name of file to generate
!
!> @warning
!! - This subroutine should be used carefully. The only aim is to add new points
!!   to a specific symmpoint instead of performing new DFT calculations. 
!!   If @b vtop is chosen wisely (high enough) we do not
!!   commit any error projecting points in the sorroundings of @b zero. Repulsive 
!!   potentials in the vecinity of a singularity are symmetric.
!! - @b Symmraw, should have been read before applying this routine.
!
!> @author A.s. Muzas - alberto.muzas@uam.es
!> @date 03/Feb/2014
!> @version 1.0
!
!> @see symmetrize, maths_mod, units_mod
!--------------------------------------------------------------------
SUBROUTINE GEN_SYMMETRIZED_RAW_INPUT(symmraw,zero,vtop,filename)
   USE MATHS_MOD
   USE UNITS_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
	IMPLICIT NONE
	! I/O variable ---------------------------
	CLASS(Symmpoint),INTENT(INOUT) :: symmraw
   TYPE(Length),INTENT(INOUT) :: zero
	TYPE(Energy),INTENT(INOUT) :: vtop
	CHARACTER(LEN=*),INTENT(IN) :: filename 
	! Local variable -------------------------
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: x,v
	INTEGER :: n
	INTEGER :: i ! Counter
   CHARACTER(LEN=27),PARAMETER :: routinename="GEN_SYMMETRIZED_RAW_INPUT: "
	! GABBA GABBA HEY! ------------------------
	CALL zero%TO_STD()
	CALL vtop%TO_STD()
	CALL SYMMETRIZE(symmraw%n,symmraw%z,symmraw%v,zero%getvalue(),vtop%getvalue(),n,x,v)
	DEALLOCATE(symmraw%z)
	DEALLOCATE(symmraw%v)
	symmraw%n = n
	ALLOCATE(symmraw%z(1:symmraw%n))
	ALLOCATE(symmraw%v(1:symmraw%n))
	DO i=1, symmraw%n
		symmraw%z(i)=x(i)
		symmraw%v(i)=v(i)
	END DO
	! Store data in filename ----------------
	OPEN(11,FILE=filename,STATUS="replace")
	WRITE(11,*) symmraw%x, symmraw%y, ' au    <----(X,Y) location in a.u.'
	WRITE(11,*) symmraw%n , " SYMM job with V(z) > (a.u.) ", vtop%getvalue()
	WRITE(11,*) "au         au       <---- everything in a.u."
	DO i=1, symmraw%n
		WRITE(11,*) symmraw%z(i), symmraw%v(i)
	END DO
	CLOSE(11)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Symmetrized input generated: ",filename)
#endif
	RETURN
END SUBROUTINE GEN_SYMMETRIZED_RAW_INPUT
!#######################################################################
!# SUBROUTINE: INITIALIZE_CRP3D_PES ######################################
!#######################################################################
!> @brief
!! Read from file useful data to initialize a CRP3D PES files
!
!> @param[in,out] this - CRP3D to be initializated
!> @param[in] filename - Name of the input file
!
!> @warning
!! - Input file structure:
!!    -# line 1: dummy line
!!    -# line 2: character(len=30),character(len=30); surface filename, surface alias
!!    -# line 3: integer(kind=4); max order environment to smooth sitios
!!    -# line 4: integer(kind=4); Number of pair potentials (N)
!!    -# line 5~5+N: character(kind=30),character(len=30); pairpot filename, pairpot alias
!!    -# line 5+N+1: integer(kind=4); number of sitios (M)
!!    -# lines 5N+2~5N+2+M: character(len=30),character(len=30); sitio filename, sitio alias
!-----------------------------------------------------------------------
SUBROUTINE READ_CRP3D(this,filename)
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   USE UNITS_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),INTENT(OUT) :: this 
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables
   INTEGER(KIND=4) :: n_pairpots,n_sites,max_order
   CHARACTER(LEN=30) :: file_surf
   CHARACTER(LEN=30),DIMENSION(:),ALLOCATABLE :: files_sites
   CHARACTER(LEN=30),DIMENSION(:),ALLOCATABLE :: files_pairpots
   CHARACTER(LEN=12),PARAMETER :: routinename="READ_CRP3D: "
   CHARACTER(LEN=2),DIMENSION(1) :: symbol
   CHARACTER(LEN=10) :: units
   REAL(KIND=8):: aux
   REAL(KIND=8),DIMENSION(1) :: aux1
   TYPE(Mass) :: masss
   INTEGER(KIND=4) :: i ! Counter
   ! HEY HO!, LET'S GO!! ------------------
   CALL this%SET_DIMENSIONS(3)
   CALL this%SET_ALIAS("CRP3D PES")
   ! Read input file
   OPEN(11,FILE=filename,STATUS="old")
   READ(11,*) !dummy line
   READ(11,*) symbol,aux,units
   CALL masss%READ(aux,units)
   CALL masss%TO_STD()
   aux1=masss%getvalue()
   CALL this%SET_ATOMS(1,symbol,aux1)
   READ(11,*) file_surf
   READ(11,*) max_order
   READ(11,*) n_pairpots
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename, file_surf)
   CALL VERBOSE_WRITE(routinename,"Max order: ",max_order)
   CALL VERBOSE_WRITE(routinename,"Pairpots: ",n_pairpots)
#endif
   ALLOCATE(files_pairpots(n_pairpots))
   DO i = 1, n_pairpots
      READ(11,*) files_pairpots(i)
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,files_pairpots(i))
#endif
   END DO
   READ(11,*) n_sites
#ifdef DEBUG
   CALL DEBUG_WRITE(routinename,"Sitios: ",n_sites)
#endif
   ALLOCATE(files_sites(n_sites))
   DO i = 1, n_sites
      READ(11,*) files_sites(i)
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,files_sites(i))
#endif
   END DO
   CLOSE(11)
   ! Initialize CRP3D
   ALLOCATE(this%all_pairpots(n_pairpots))
   DO i = 1,n_pairpots
      CALL this%all_pairpots(i)%READ(files_pairpots(i)) 
   END DO
   ALLOCATE(this%all_sites(n_sites))
   DO i = 1, n_sites
     CALL this%all_sites(i)%READ(files_sites(i)) 
   END DO
   CALL this%surf%INITIALIZE(file_surf)
   this%max_order = max_order
   RETURN
END SUBROUTINE READ_CRP3D
!#######################################################################
!# SUBROUTINE: EXTRACT_VASINT_CRP3D ######################################
!#######################################################################
!> @brief
!! Extracts the potential at long distances from every value stored in pairpots
!! and sites that belongs to @b thispes. Basically, this procedure sets
!! our zero potential value at the vacuum.
!
!> @param[in,out] thispes - A CRP3D subtype variable
!
!> @warning 
!! - The value of the potential at long distance should be the same in all
!!   pair potentials, read from files previously
!! - At least one pair potential should've defined in the CRP3D variable
!! - Obviously, @ thispes should contain reliable data
!! - Only modifies things inside @b interz
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 04/Feb/2014
!> @version 1.0
!------------------------------------------------------------------------
SUBROUTINE EXTRACT_VASINT_CRP3D(thispes)
   ! Initial declarations
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),INTENT(INOUT) :: thispes
   ! Local variables
   INTEGER(KIND=4) :: npairpots, nsites
   INTEGER(KIND=4) :: i,j ! counters
   REAL(KIND=8) :: control_vasint
   CHARACTER(LEN=22) :: routinename="EXTRACT_VASINT_CRP3D: "
   ! Run section ------------------------
   npairpots=size(thispes%all_pairpots)
   control_vasint=thispes%all_pairpots(1)%vasint
   DO i = 1, npairpots
      IF (thispes%all_pairpots(1)%vasint/=control_vasint) THEN
         WRITE(0,*) "EXTRACT_VASINT_CRP3D ERR: Incoherences in vasint values found"
         WRITE(0,*) "EXTRACT_VASINT_CRP3D ERR: vasint's value at pairpot",1,control_vasint
         WRITE(0,*) "EXTRACT_VASINT_CRP3D ERR: vasint's value at pairpot",i,control_vasint
         CALL EXIT(1)
      END IF
      DO j = 1, thispes%all_pairpots(i)%n
         thispes%all_pairpots(i)%interz%f(j)=thispes%all_pairpots(i)%interz%f(j)-thispes%all_pairpots(i)%vasint
      END DO
#ifdef DEBUG
      CALL DEBUG_WRITE(routinename,"Vasint extracted from pair potential ",i)
#endif
   END DO
   nsites=size(thispes%all_sites)
   DO i = 1, nsites
      DO j = 1, thispes%all_sites(i)%n
         thispes%all_sites(i)%interz%f(j)=thispes%all_sites(i)%interz%f(j)-thispes%all_pairpots(1)%vasint
      END DO
#ifdef DEBUG
      CALL DEBUG_WRITE(routinename,"Vasint extracted from pair site ",i)
#endif
   END DO
   RETURN
END SUBROUTINE EXTRACT_VASINT_CRP3D
!############################################################
! SUBROUTINE: SMOOTH_CRP3D ####################################
!############################################################
!> @brief
!! Substracts to all sitio potentials the pair-potential interactions 
!! with all atoms of the surface inside an environment of order @b
!! max order
! 
!> @param[in,out] thispes - CRP3D PES variable to be smoothed (only it's sites)
!
!> @warning
!! - Vasint should've  been extracted before
!! - Only afects values stored in @b interz variable
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 05/Feb/2014
!> @version 1.0
!
!> @see extract_vasint_CRP3D, initialize_CRP3D_pes
!-----------------------------------------------------------
SUBROUTINE SMOOTH_CRP3D(thispes)
   ! Initial declaraitons
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   CLASS(CRP3D),INTENT(INOUT) :: thispes
   ! Local variables
   REAL(KIND=8),DIMENSION(3) :: A,aux
   REAL(KIND=8) :: dummy
   INTEGER(KIND=4) :: i,j,k,l ! counters
   INTEGER(KIND=4) :: npairpots,nsites
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: interac,dvdzr ! collects interactions and corrections to
   CHARACTER(LEN=14) :: routinename="SMOOTH_CRP3D: "
   ! Run section ----------
   nsites = size(thispes%all_sites)
   npairpots = size(thispes%all_pairpots)
   ALLOCATE(interac(0:thispes%max_order))
   ALLOCATE(dvdzr(0:thispes%max_order))
   DO i = 1, nsites ! loop over sites
      aux(1) = thispes%all_sites(i)%x
      aux(2) = thispes%all_sites(i)%y
      DO j = 1, thispes%all_sites(i)%n ! loop over pairs v,z
         aux(3)=thispes%all_sites(i)%z(j)
         A=aux
         A(1:2)=thispes%surf%cart2surf(aux(1:2))
         DO l = 1, npairpots ! loop over pairpots
            DO k = 0, thispes%max_order ! loop over environment orders
               CALL INTERACTION_AENV(k,A,thispes%surf,thispes%all_pairpots(l),interac(k),dvdzr(k),dummy,dummy)
            END DO
            thispes%all_sites(i)%interz%f(j)=thispes%all_sites(i)%interz%f(j)-sum(interac)
            IF (j.EQ.1) THEN
               thispes%all_sites(i)%dz1=thispes%all_sites(i)%dz1-SUM(dvdzr) ! correct first derivative
            ELSE IF (j.EQ.thispes%all_sites(i)%n) THEN
               thispes%all_sites(i)%dz2=thispes%all_sites(i)%dz2-SUM(dvdzr) ! correct first derivative
            END IF
         END DO
      END DO
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Site smoothed: ",i)
#endif
   END DO
   RETURN
END SUBROUTINE SMOOTH_CRP3D
!############################################################
! SUBROUTINE: INTERPOL_Z_CRP3D ################################
!############################################################
!> @brief
!! Interpolates all sitios and pairpot potentials that belong
!! to an specific CRP3D PES. Sitio potentials are smoothed
!
!> @param[in,out] thispes - CRP3D PES to be interpolated 
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 05/Feb/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INTERPOL_Z_CRP3D(thispes)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D) :: thispes
   ! Local variables
   INTEGER(KIND=4) :: nsites,npairpots
   REAL(KIND=8) :: dz1,dz2
   INTEGER(KIND=4) :: i ! counters
   ! Run secton------------------------
   nsites=size(thispes%all_sites)
   npairpots=size(thispes%all_pairpots)
   CALL thispes%EXTRACT_VASINT()
   DO i = 1, npairpots ! loop pairpots
      dz1=thispes%all_pairpots(i)%dz1
      dz2=thispes%all_pairpots(i)%dz2
      CALL thispes%all_pairpots(i)%interz%INTERPOL(dz1,1,dz2,1)
   END DO
   CALL thispes%SMOOTH()
   DO i = 1, nsites
      dz1=thispes%all_sites(i)%dz1
      dz2=thispes%all_sites(i)%dz2
      CALL thispes%all_sites(i)%interz%INTERPOL(dz1,1,dz2,1) 
   END DO
   RETURN
END SUBROUTINE INTERPOL_Z_CRP3D
!##################################################################
!# SUBROUTINE: INTERACTION_AP #####################################
!##################################################################
!> @brief
!! Calculates the interaction between two given 3D Points with X and Y
!! in surface coordinates. The potential is extracted from pairpot,
!! which should've been interpolated before
!
!> @param[in] A, P - Pair of 3D points that are interacting through @ pairpot potential
!> @param[in] surf - Periodic surface
!> @param[in] pairpot - Pair potential. Source of @f$V(r)@f$, where @b r is the distance
!!                      between @b A and @b P. @f$V(r)@f$ is just a shifted version of
!!                      @f$V(z)@f$.
!> @param[out] interac - Actual interaction between @b A and @b P
!> @param[out] dvdz_corr - Corrections to first derivatives: 
!!                         @f$\frac {\partial V(r)}{\partial z}@f$
!> @param[out] dvdx_corr - Corrections to first derivatives: 
!!                         @f$\frac {\partial V(r)}{\partial z}@f$
!> @param[out] dvdy_corr - Corrections to first derivatives: 
!!                         @f$\frac {\partial V(r)}{\partial y}@f$
!
!> @warning
!! - Sitio and Pairpot should have been shifted so that the interaction
!!   at large distance from the surface is 0 (extracted vasint)
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 04/Feb/2014
!> @version 1.0
!
!> @see extra documentation
!------------------------------------------------------------------
SUBROUTINE INTERACTION_AP(A,P,surf,pairpot,interac,dvdz_corr,dvdx_corr,dvdy_corr)
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   ! I/O VAriables --------------------------------------------
   REAL(KIND=8),DIMENSION(3),INTENT(IN) :: A, P
   TYPE(Surface),INTENT(IN) :: surf
   TYPE(Pair_pot),INTENT(IN) :: pairpot
   REAL(KIND=8),INTENT(OUT) :: interac,dvdz_corr,dvdx_corr,dvdy_corr
   ! Local variables ------------------------------------------
   REAL(KIND=8) :: r ! distance
   REAL(KIND=8) :: modA ! A modulus
   REAL(KIND=8) :: modP_plus ! P modulus plus other term
   REAL(KIND=8), DIMENSION(3) :: vect, vect2
   REAL(KIND=8) :: aux
   CHARACTER(LEN=16), PARAMETER :: routinename = "INTERACTION_AP: "
   INTEGER :: i ! Counter
   ! GABBA, GABBA HEY! ----------------------------------------
   ! Find the distance between A and P
#ifdef DEBUG
   CALL DEBUG_WRITE(routinename, "A vector: (In surface coordinates)")
   CALL DEBUG_WRITE(routinename,"A_1: ",A(1))
   CALL DEBUG_WRITE(routinename,"A_2: ",A(2))
   CALL DEBUG_WRITE(routinename,"A_3: ",A(3))
   CALL DEBUG_WRITE(routinename, "P vector: (In surface coordinates)")
   CALL DEBUG_WRITE(routinename,"P_1: ",P(1))
   CALL DEBUG_WRITE(routinename,"P_2: ",P(2))
   CALL DEBUG_WRITE(routinename,"P_3: ",P(3))
#endif  
   vect=A
   vect(1:2)=matmul(surf%metricsurf_mtrx,A(1:2))
   modA = dot_product(A,vect)
   FORALL(i=1:3) vect(i) = P(i)-2.D0*A(i)
   vect2=vect
   vect2(1:2)=matmul(surf%metricsurf_mtrx,vect(1:2))
   modP_plus = dot_product(P,vect2)
   r=dsqrt(modA+modP_plus)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename, "Distance: ",r)
#endif
   IF (r.GT.pairpot%z(pairpot%n)) THEN
      interac = 0.D0 ! if distance too high, interaction 0 
      dvdz_corr = 0.D0
      dvdx_corr = 0.D0
      dvdy_corr = 0.D0
#ifdef DEBUG
      CALL DEBUG_WRITE(routinename,"A is really far away from P. Negligible interaction.")
#endif
   ELSE
      vect = A
      vect(1:2) = surf%surf2cart(A(1:2))
      vect2 = P
      vect2(1:2) = surf%surf2cart(P(1:2))
#ifdef DEBUG
      CALL DEBUG_WRITE(routinename, "A vector: (In cartesian coordinates)")
      CALL DEBUG_WRITE(routinename,"A_1: ",vect(1))
      CALL DEBUG_WRITE(routinename,"A_2: ",vect(2))
      CALL DEBUG_WRITE(routinename,"A_3: ",vect(3))
      CALL DEBUG_WRITE(routinename, "P vector: (In cartesian coordinates)")
      CALL DEBUG_WRITE(routinename,"P_1: ",vect2(1))
      CALL DEBUG_WRITE(routinename,"P_2: ",vect2(2))
      CALL DEBUG_WRITE(routinename,"P_3: ",vect2(3))
#endif
      ! We need the shifted version of the potential
      interac = pairpot%interz%getvalue(r,pairpot%rumpling)
      aux=pairpot%interz%getderiv(r,pairpot%rumpling) ! better performance
      dvdz_corr = aux*(vect(3)-vect2(3))/r
      dvdx_corr = aux*(vect(1)-vect2(1))/r
      dvdy_corr = aux*(vect(2)-vect2(2))/r
   END IF
#ifdef DEBUG
   CALL DEBUG_WRITE(routinename,"Repul. interaction: ",interac)
   CALL DEBUG_WRITE(routinename,"Correction to derivative dvdx: ",dvdx_corr)
   CALL DEBUG_WRITE(routinename,"Correction to derivative dvdy: ",dvdy_corr)
   CALL DEBUG_WRITE(routinename,"Correction to derivative dvdz: ",dvdz_corr)
#endif
   RETURN
END SUBROUTINE INTERACTION_AP
!##################################################################
!# SUBROUTINE: INTERACTION_AENV ###################################
!##################################################################
!> @brief
!! Calculates the complete interaction between a given point @b A with an
!! environment of order @b n.
!
!> @param[in] n - Environment order
!> @param[in] A - Point in space whose X and Y coordinates are in surface units
!!                and inside the IWS cell defined by @b surf
!> @param[in] surf - Periodic surface 
!> @param[in] pairpot - Pair potential which defines interaction potential
!> @param[out] interac - Total interaction with the environment defined
!> @param[out] dvdx_term - 
!> @param[out] dvdy_term - 
!> @param[out] dvdz_term - 
!> @param[in] gnp - Optional logical variable which tells this routine to print extra
!!                  output for GNP_MAP_INTERACT in file gnp-data.dat. Just for testing jobs
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 04/Feb/2014
!> @version 1.0
!------------------------------------------------------------------
SUBROUTINE INTERACTION_AENV(n,A,surf,pairpot,interac,dvdz_term,dvdx_term,dvdy_term,gnp)
   IMPLICIT NONE
   ! I/O VAriables ------------------------------------------
   INTEGER,INTENT(IN) :: n
   REAL(KIND=8),DIMENSION(3), INTENT(IN) :: A
   TYPE(Pair_pot),INTENT(IN) :: pairpot
   TYPE(Surface),INTENT(IN) :: surf
   LOGICAL,INTENT(IN),OPTIONAL :: gnp
   REAL(KIND=8),INTENT(OUT) :: interac, dvdz_term, dvdx_term, dvdy_term
   ! Local variables ----------------------------------------
   REAL(KIND=8),DIMENSION(3) :: P
   REAL(KIND=8) :: dummy1, dummy2, dummy3, dummy4
   REAL(KIND=8) :: atomx, atomy
   INTEGER :: pairid
   INTEGER :: i, k ! Counters
   CHARACTER(LEN=18), PARAMETER :: routinename = "INTERACTION_AENV: "
   ! SUSY IS A HEADBANGER !!!! -------------------
   ! Defining some aliases to make the program simpler:
   pairid = pairpot%id
   P(3) = pairpot%rumpling ! rumpling associated with pairpot
   ! Case n = 0
   interac=0.D0
   dvdz_term=0.D0
   dvdx_term=0.D0
   dvdy_term=0.D0
   IF (present(gnp).AND.gnp) OPEN(10,FILE="gnp-data.dat",STATUS="replace")
   IF (n.EQ.0) THEN 
      DO i=1, surf%atomtype(pairid)%n
         P(1) = surf%atomtype(pairid)%atom(i,1)
         P(2) = surf%atomtype(pairid)%atom(i,2)
         CALL INTERACTION_AP(A,P,surf,pairpot,dummy1,dummy2,dummy3,dummy4)
         interac = interac +dummy1
         dvdz_term = dvdz_term + dummy2
         dvdx_term = dvdx_term + dummy3
         dvdy_term = dvdy_term + dummy4
         IF (present(gnp).AND.gnp) WRITE(10,*) P(:),dummy1,dummy2,dummy3,dummy4
      END DO
      IF (present(gnp).AND.gnp) CLOSE(10)
      RETURN
   ELSE IF (n.GT.0) THEN
      DO i=1, surf%atomtype(pairid)%n
         atomx = surf%atomtype(pairid)%atom(i,1)
         atomy = surf%atomtype(pairid)%atom(i,2)
         DO k= -n,n
            P(1) = atomx + DFLOAT(n)
            P(2) = atomy + DFLOAT(k)
            CALL INTERACTION_AP(A,P,surf,pairpot,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
            IF (present(gnp).AND.gnp) WRITE(10,*) P(:),dummy1,dummy2,dummy3,dummy4
            P(1) = atomx + DFLOAT(-n)
            P(2) = atomy + DFLOAT(k)
            CALL INTERACTION_AP(A,P,surf,pairpot,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
            IF (present(gnp).AND.gnp) WRITE(10,*) P(:),dummy1,dummy2,dummy3,dummy4
         END DO
         DO k= -n+1, n-1
            P(1) = atomx + DFLOAT(k)
            P(2) = atomy + DFLOAT(n)
            CALL INTERACTION_AP(A,P,surf,pairpot,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
            IF (present(gnp).AND.gnp) WRITE(10,*) P(:),dummy1,dummy2,dummy3,dummy4
            P(1) = atomx + DFLOAT(k)
            P(2) = atomy + DFLOAT(-n)
            CALL INTERACTION_AP(A,P,surf,pairpot,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
            IF (present(gnp).AND.gnp) WRITE(10,*) P(:),dummy1,dummy2,dummy3,dummy4
         END DO
      END DO
      IF(present(gnp).AND.gnp) CLOSE(10)
      RETURN
   ELSE
      WRITE(0,*) "INTERACTION_AENV ERR: Wrong Environment order."
      CALL EXIT(1)
   END IF
END SUBROUTINE INTERACTION_AENV
!############################################################
!# SUBROUTINE: GET_V_AND_DERIVS_CRP3D #########################
!############################################################
!> @brief
!! Subroutine that calculates the 3D potential for a point A and
!! its derivatives in cartesian coordinates.
!
!> @param[in] thispes - CRP3D PES
!> @param[in] X - Point in space to calculate the potential and it's derivatives. Cartesian's
!> @param[out] v - Value of the potential at X
!> @param[out] dvdu - derivatives, cartesian coordinates 
!
!> @warning
!! - All sitios should have been interpolated in Z direction previously. It
!!   is optinal to extract "vasint".
!! - Corrugation, MUST have been extracted from sitios previously (smoothed) 
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 06/Feb/2014
!> @version 1.0
!------------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_CRP3D(thispes,X,v,dvdu)
#ifdef DEBUG
   USE DEBUG_MOD
#endif
	IMPLICIT NONE
	! I/O variables
   CLASS(CRP3D),TARGET,INTENT(IN) :: thispes
	REAL(KIND=8),DIMENSION(3), INTENT(IN) :: X
	REAL(KIND=8),INTENT(OUT) :: v
	REAL(KIND=8),DIMENSION(3),INTENT(OUT) :: dvdu
	! Local variables
	TYPE(Fourier2d) :: interpolxy
   INTEGER(KIND=4) :: nsites,npairpots
	REAL(KIND=8),DIMENSION(3) :: A, deriv
	REAL(KIND=8) :: pot
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: f, dfdz ! arguments to the xy interpolation
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: xy ! arguments to the xy interpolation
	INTEGER :: i, j,k ! counters
	! Pointers
	REAL(KIND=8), POINTER :: zmax
	TYPE(Pair_pot),POINTER :: pairpot
   CHARACTER(LEN=24),PARAMETER :: routinename="GET_V_AND_DERIVS_CRP3D: "
	zmax => thispes%all_sites(1)%interz%x(thispes%all_sites(1)%n)
   npairpots = size(thispes%all_pairpots)
   nsites = size(thispes%all_sites)
	! GABBA, GABBA HEY! ----------------------
	IF (X(3).GT.zmax) THEN
		dvdu = 0.D0
		v = 0.D0
		RETURN
	END IF
	! For the corrugation addition, we need to project upon IWS cell
	A = X
   A(1:2) = thispes%surf%project_unitcell(X(1:2))
   A(1:2) = thispes%surf%cart2surf(A(1:2)) ! go to surface coordinates, but now, inside the unit cell 
	!
	pot = 0.D0 ! Initialization value
	v = 0.D0 ! Initialization value
	FORALL(i=1:3) 
		dvdu(i) = 0.D0 ! Initialization value
		deriv(i) = 0.D0 ! Initialization value
	END FORALL 
#ifdef DEBUG   
	CALL VERBOSE_WRITE(routinename,"Proceeding to calculate repulsive interactions")
	CALL VERBOSE_WRITE(routinename,"Max_order: ", thispes%max_order)
#endif
	DO j=1, npairpots
		pairpot => thispes%all_pairpots(j)
		DO k=0,thispes%max_order
#ifdef DEBUG
			CALL VERBOSE_WRITE(routinename, "---> Invoking order", k)
#endif
			CALL INTERACTION_AENV(k,A,thispes%surf,pairpot,pot,deriv(3),deriv(1),deriv(2))
			v = v + pot
			FORALL (i=1:3) dvdu(i) = dvdu(i) + deriv(i)
		END DO
#ifdef DEBUG
		CALL VERBOSE_WRITE(routinename, "Repulsive interaction: ", pot)
		CALL VERBOSE_WRITE(routinename, "Cumulative repulsive interaction:", v )
		CALL VERBOSE_WRITE(routinename, "dvdx: ",deriv(1))
		CALL VERBOSE_WRITE(routinename, "dvdy: ",deriv(2))
		CALL VERBOSE_WRITE(routinename, "dvdz: ",deriv(3))
		CALL VERBOSE_WRITE(routinename, "Cumulative dvdx: ",dvdu(1))
		CALL VERBOSE_WRITE(routinename, "Cumulative dvdy: ",dvdu(2))
		CALL VERBOSE_WRITE(routinename, "Cumulative dvdz: ",dvdu(3))
#endif
	END DO
	! Now, we have all the repulsive interaction and corrections to the derivarives
	! stored in v(:) and dvdu(:) respectively.
	! Let's get v and derivatives from xy interpolation of the corrugationless function
   ALLOCATE(f(nsites))
   ALLOCATE(xy(nsites,2))
   ALLOCATE(dfdz(nsites))
   DO i=1,nsites
      xy(i,1)=thispes%all_sites(i)%x
      xy(i,2)=thispes%all_sites(i)%y
      CALL thispes%all_sites(i)%interz%GET_V_AND_DERIVS(X(3),f(i),dfdz(i))
   END DO
	CALL interpolxy%READ(xy,f,dfdz)
   CALL interpolxy%SET_COEFF(thispes%surf)
   CALL interpolxy%GET_F_AND_DERIV(thispes%surf,X,pot,deriv)
	! Corrections from the smoothing procedure
	v = v + pot
	FORALL(i=1:3) dvdu(i) = dvdu(i) + deriv(i)
	RETURN
END SUBROUTINE GET_V_AND_DERIVS_CRP3D
!######################################################################
! SUBROUTINE: PLOT_DATA_SYMMPOINT #####################################
!######################################################################
!> @brief
!! Creates a data file called "filename" with the data z(i) and v(i) 
!! for a Symmpoint type variable
!----------------------------------------------------------------------
SUBROUTINE PLOT_DATA_SYMMPOINT(this,filename)
   IMPLICIT NONE
   ! I/O Variables ----------------------
   CLASS(Symmpoint), INTENT(IN) :: this
   CHARACTER(LEN=*), INTENT(IN) :: filename
   ! Local variables --------------------
   INTEGER :: i ! Counter
   ! GABBA GABBA HEY!! ------------------
   OPEN(11, FILE=filename, STATUS="replace")
   DO i=1,this%n
      WRITE(11,*) this%z(i),this%v(i)
   END DO
   CLOSE (11)
   WRITE(*,*) "PLOT_DATA_SYMMPOINT: ",this%alias,filename," file created"
END SUBROUTINE PLOT_DATA_SYMMPOINT
!#######################################################################
! SUBROUTINE: GRAPH_V3D_POTMAP_XY
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (X,Y) of the PES. 
!
!> @param[in] thispes - CRP3D PES used
!> @param[in] filename - Name of the file to print the output
!> @param[in] init_xyz - Initial position to start the scan (a.u.)
!> @param[in] nxpoints - Number of points in X axis (auxiliar cartesian coordinates)
!> @param[in] nypoints - Number of points in Y axis (auxiliar cartesian coordinates)
!> @param[in] Lx - Length of X axis (a.u.)
!> @param[in] Ly - Length of Y axis (a.u.)
!
!> @warning
!! - Z parameter is taken from @b init_xyz and is constant
!
!> @author A.S. Muzas
!> @date 08/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT_XYMAP_CRP3D(thispes,filename,init_xyz,nxpoints,nypoints,Lx,Ly)
   IMPLICIT NONE
   CLASS(CRP3D),INTENT(IN) :: thispes
   REAL*8,DIMENSION(3),INTENT(IN) :: init_xyz ! Initial position to start the scan (in a.u.)
   INTEGER,INTENT(IN) :: nxpoints, nypoints ! number of points in XY plane
   CHARACTER(LEN=*),INTENT(IN) :: filename ! filename
   REAL*8,INTENT(IN) :: Lx ! Length of X axis 
   REAL*8,INTENT(IN) :: Ly ! Length of X axis 
   ! Local variables
   REAL*8 :: xmin, ymin, xmax, ymax, z
   REAL*8, DIMENSION(3) :: r, dvdu
   REAL*8 :: xdelta, ydelta
   INTEGER :: xinpoints, nxdelta
   INTEGER :: yinpoints, nydelta
   INTEGER :: i, j ! counters
   REAL*8 :: v ! potential
	! GABBA, GABBA HEY! ---------
   xmin = init_xyz(1)
   ymin = init_xyz(2)
   z = init_xyz(3)
   xmax = init_xyz(1)+Lx
   ymax = init_xyz(2)+Ly
   ! For X, grid parameters
   xinpoints=nxpoints-2
   nxdelta=nxpoints-1
   xdelta=Lx/DFLOAT(nxdelta)
   ! For Y, grid parameters
   yinpoints=nypoints-2
   nydelta=nypoints-1
   ydelta=(ymax-ymin)/DFLOAT(nydelta)
   ! Let's go! 
   ! 1st XY point
   OPEN(11,file=filename,status="replace")
   r(1) = xmin
   r(2) = ymin
   r(3) = z
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v
   END DO
   r(2) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v
   ! inpoints in XY
   DO i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      r(3) = z
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v
      DO j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
         WRITE(11,*) r(1), r(2), v
      END DO
      r(2) = ymax
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v
   END DO
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   r(3) = z
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v
   END DO
   r(2) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_XYMAP_CRP3D
!#######################################################################
! SUBROUTINE: PLOT_DIRECTION1D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. To define 
!! the direction, the angle alpha is given. 
!
!> @param[in] thispes - CRP3D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] angle - Angle between the surface vector S1 and the direction of the
!!                    cut. It should be given in degrees.
!> @param[in] z - Distance to the surface. All points have the same z.
!> @param[in] L - Length of the graphic
!
!> @warning
!! - The graph starts always at 0,0
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 09/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT_DIRECTION1D_CRP3D(thispes,filename,npoints,angle,z,L)
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP3D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*), INTENT(IN) :: filename
   REAL*8, INTENT(IN) :: z, angle
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL*8 :: delta,L,v,s,alpha
   REAL*8 :: xmax, xmin, ymax, ymin 
   REAL*8, DIMENSION(3) :: r, dvdu
   INTEGER :: i ! Counter
   CHARACTER(LEN=24), PARAMETER :: routinename = "PLOT_DIRECTION1D_CRP3D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT_DIRECTION1D_CRP3D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   ! Change alpha to radians
   alpha = angle * PI / 180.D0
   !
   xmin = 0.D0
   ymin = 0.D0
   xmax = L*DCOS(alpha)
   ymax = L*DSIN(alpha)
   r(3) = z
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=L/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(1) = xmin
   r(2) = ymin
   s = DSQRT(r(1)**2.D0+r(2)**2.D0)
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) s, v , dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   ! cycle for inpoints
   DO i=1, inpoints
      r(1)=xmin+(DFLOAT(i)*delta)*DCOS(alpha)
      r(2)=ymin+(DFLOAT(i)*delta)*DSIN(alpha)
      s = DSQRT(r(1)**2.D0+r(2)**2.D0)
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) s, v , dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   END DO
   ! Final value
   r(1) = xmax
   r(2) = ymax
   s = DSQRT(r(1)**2.D0+r(2)**2.D0)
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) s, v, dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT_DIRECTION1D_CRP3D
!###########################################################
!# FUNCTION: is_allowed_CRP3D
!###########################################################
!> @brief
!! Determines if the potential can be calculated in this point
!-----------------------------------------------------------
LOGICAL FUNCTION is_allowed_CRP3D(thispes,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),TARGET,INTENT(IN) :: thispes
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),POINTER :: xmin,xmax
   ! Run section
   xmin => thispes%all_sites(1)%z(1)
   xmax => thispes%all_sites(1)%z(thispes%all_sites(1)%n)
   IF (size(x)/=3) THEN
      WRITE(0,*) "is_allowed_CRP3D ERR: checked array without the correct number of dimensions: 3"
      is_allowed_CRP3D=.FALSE.
      CALL EXIT(1)
   ENDIF
   IF ((x(3).LT.xmin).OR.(x(3).GT.xmax)) THEN
     is_allowed_CRP3D=.FALSE. 
   END IF  
   is_allowed_CRP3D=.TRUE.
   RETURN
END FUNCTION is_allowed_CRP3D
END MODULE CRP3D_MOD
