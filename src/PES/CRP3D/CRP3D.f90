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
   use PES_MOD
   use SYSTEM_MOD
   use CUBICSPLINES_MOD
   use LINK_FOURIER2D_MOD
   use LINK_FUNCTION1D_MOD
   use AOTUS_MODULE, only: flu_State, OPEN_CONFIG_FILE, CLOSE_CONFIG, AOT_GET_VAL
   use AOT_TABLE_MODULE, only: AOT_TABLE_OPEN, AOT_TABLE_CLOSE, AOT_TABLE_LENGTH, AOT_TABLE_GET_VAL
   use MATHS_MOD, only: ORDER_VECT, ORDER, SYMMETRIZE
   use UNITS_MOD, only: pi
#ifdef DEBUG
   use DEBUG_MOD, only: VERBOSE_WRITE, DEBUG_WRITE
#endif
! Initial declarations
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: CRP3D_details
!> @brief
!! Type used in jobs to generate good inputs for a CRP3D PES
!
!> @param vasint - potential when the atom is far from the surface (in the vacuum)
!> @param dfin - arbitrary boundary: 1st derivative of all potentials at zgrid(nzgrid)
!> @param nrumpling - number of different rumplings defined
!> @param rumpling - storage of rumplings
!> @param zeropos - position of ith rumpling in the grid
!> @param nzgrid - dimension of zgrid
!> @param zgrid - grid proposed to generate input files
!> @param first, last - first and last points in the grid

!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jul/2014
!> @version 1.0
!----------------------------------------------------------------
TYPE :: CRP3D_details
PRIVATE
   REAL(KIND=8):: vasint 
   REAL(KIND=8):: dfin 
   INTEGER(KIND=4):: nrumpling 
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: rumpling 
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE :: zeropos 
   INTEGER(KIND=4) :: nzgrid 
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: zgrid 
   REAL(KIND=8):: first, last 
   CONTAINS
      PROCEDURE,PUBLIC :: READ => READ_CRP3D_details
END TYPE CRP3D_details
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
   CHARACTER(LEN=:),ALLOCATABLE:: filename
   CHARACTER(LEN=:),ALLOCATABLE:: alias
   INTEGER(KIND=4):: n
   REAL(KIND=8):: x
   REAL(KIND=8):: y
   CHARACTER(LEN=10):: units_z
   CHARACTER(LEN=10):: units_v
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: z
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: v
   REAL(KIND=8):: dz1
   REAL(KIND=8):: dz2
   TYPE(Csplines),PUBLIC:: interz 
   CONTAINS
      PROCEDURE,PUBLIC:: READ_RAW => READ_SYMMPOINT_RAW
      PROCEDURE,PUBLIC:: GET_SYMM_RAW => GET_SYMMETRIZED_RAW_INPUT
      PROCEDURE,PUBLIC:: PLOT_DATA => PLOT_DATA_SYMMPOINT
      PROCEDURE,PUBLIC:: PLOT => PLOT_INTERPOL_SYMMPOINT
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
!!                      the singularity is plazed at zero. Should be a negative number o 0.D0.
!
!> @author A.P. Muzas
!> @version 1.0
!> @date 21/Jan/2014
!------------------------------------------------------------------------------
TYPE,EXTENDS(Symmpoint) :: Pair_pot
   PRIVATE
   INTEGER(KIND=4):: id 
   REAL(KIND=8):: vasint
	REAL(KIND=8):: rumpling
   CONTAINS
      ! Initialization block
      PROCEDURE,PUBLIC:: READ => READ_STANDARD_PAIRPOT
      ! Tools block
      PROCEDURE,PUBLIC:: GET_V_AND_DERIVS => GET_V_AND_DERIVS_PAIRPOT
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
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: dvdx,dvdy,dvdz 
   CONTAINS
      PROCEDURE,PUBLIC:: READ => READ_STANDARD_SITIO
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
TYPE,EXTENDS(PES) :: CRP3D
   INTEGER(KIND=4):: max_order
   TYPE(Pair_pot),DIMENSION(:),ALLOCATABLE:: all_pairpots
   TYPE(Sitio),DIMENSION(:),ALLOCATABLE:: all_sites
   INTEGER(KIND=4),DIMENSION(:,:),ALLOCATABLE:: kList
   character(len=1),dimension(:),allocatable:: parityList
   character(len=2),dimension(:),allocatable:: irrepList
   CLASS(Function1d),ALLOCATABLE:: dampfunc
   CONTAINS
      ! Initialization block
      PROCEDURE,PUBLIC:: READ => READ_CRP3D
      PROCEDURE,PUBLIC:: INITIALIZE => INITIALIZE_CRP3D
      ! Get block 
      PROCEDURE,PUBLIC:: GET_V_AND_DERIVS => GET_V_AND_DERIVS_CRP3D
      PROCEDURE,PUBLIC:: GET_V_AND_DERIVS_CORRECTION => GET_V_AND_DERIVS_CORRECTION_CRP3D
      PROCEDURE,PUBLIC:: GET_V_AND_DERIVS_SMOOTH => GET_V_AND_DERIVS_SMOOTH_CRP3D
      PROCEDURE,PUBLIC:: GET_REPUL_CORRECTIONS => GET_REPUL_CORRECTIONS_CRP3D
      PROCEDURE,PUBLIC:: getpot => getpot_crp3d
      PROCEDURE,PUBLIC:: getrumpling => getrumpling_CRP3D
      ! Enquire block
      PROCEDURE,PUBLIC:: is_allowed => is_allowed_CRP3D
      ! Tools block
      PROCEDURE,PUBLIC:: EXTRACT_VASINT => EXTRACT_VASINT_CRP3D
      PROCEDURE,PUBLIC:: SMOOTH => SMOOTH_CRP3D
      PROCEDURE,PUBLIC:: INTERPOL => INTERPOL_Z_CRP3D
      PROCEDURE,PUBLIC:: RAWINTERPOL => RAWINTERPOL_Z_CRP3D
      PROCEDURE,PUBLIC:: SET_FOURIER_SYMMETRY => SET_FOURIER_SYMMETRY_CRP3D
      ! Plot tools
      PROCEDURE,PUBLIC:: PLOT_XYMAP => PLOT_XYMAP_CRP3D
      PROCEDURE,PUBLIC:: PLOT_XYMAP_SMOOTH => PLOT_XYMAP_CRP3D
      PROCEDURE,PUBLIC:: PLOT_XYMAP_CORRECTION => PLOT_XYMAP_CORRECTION_CRP3D
      PROCEDURE,PUBLIC:: PLOT_DIRECTION1D => PLOT_DIRECTION1D_CRP3D
      PROCEDURE,PUBLIC:: PLOT_DIRECTION1D_SMOOTH => PLOT_DIRECTION1D_SMOOTH_CRP3D
      PROCEDURE,PUBLIC:: PLOT_DIRECTION1D_CORRECTION => PLOT_DIRECTION1D_CORRECTION_CRP3D
      PROCEDURE,PUBLIC:: PLOT_SITIOS => PLOT_SITIOS_CRP3D
      PROCEDURE,PUBLIC:: PLOT_PAIRPOTS => PLOT_PAIRPOTS_CRP3D
      PROCEDURE,PUBLIC:: PLOT_Z => PLOT_Z_CRP3D
END TYPE CRP3D
!///////////////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: READ_CRP3D_details 
!###########################################################
!> @brief
!! Reads from file interesting details to generate new inputs for a CRP3D PES.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jul/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE READ_CRP3D_details(this,filename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D_details),INTENT(OUT):: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! WARNING: unit used to read
   INTEGER(KIND=4) :: runit=233
   ! Local variables
   INTEGER(KIND=4) :: i,k ! counters
   TYPE(Energy):: en
   TYPE(Length):: len
   REAL(KIND=8) :: raux,raux2
   CHARACTER(LEN=10) :: units
   CHARACTER(LEN=20),PARAMETER :: routinename="READ_CRP3D_details: "
   CHARACTER(LEN=4) :: control
   LOGICAL:: exists_zero
   ! Run section
   OPEN (runit,FILE=filename,STATUS="old",ACTION="read")
   READ(runit,*) ! dummy line
   READ(runit,*) raux,units
   CALL en%READ(raux,units)
   CALL en%TO_STD()
   this%vasint=en%getvalue()
   READ(runit,*) this%nrumpling
   ALLOCATE(this%rumpling(this%nrumpling))
   this%dfin = 0.D0
   DO i=1,this%nrumpling
      READ(runit,*) raux,units
      CALL len%READ(raux,units)
      CALL len%TO_STD()
      this%rumpling(i)=len%getvalue()
   END DO
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Rumplings: ",this%rumpling(:))
#endif
   READ(runit,*) this%nzgrid,control,raux,raux2,units
   ALLOCATE(this%zgrid(this%nzgrid))
   CALL len%READ(raux,units)
   CALL len%TO_STD()
   this%first=len%getvalue()
   CALL len%READ(raux2,units)
   CALL len%TO_STD()
   this%last=len%getvalue()
   SELECT CASE(control)
      CASE("MANU") ! manual grid input
         DO i = 1, this%nzgrid
            READ(runit,*) raux
            CALL len%READ(raux,units)
            CALL len%TO_STD()
            this%zgrid(i)=len%getvalue()
         END DO
         CALL ORDER_VECT(this%zgrid)
         ALLOCATE(this%zeropos(this%nrumpling))
         FORALL(i=1:this%nrumpling) this%zeropos(i)=0
         DO i=1,this%nzgrid
            DO k=1,this%nrumpling
               SELECT CASE(this%zgrid(i)==this%rumpling(k))
                  CASE(.TRUE.)
                     this%zeropos(k)=i
                  CASE(.FALSE.)
                     ! do nothing
               END SELECT
            END DO
         END DO
         exists_zero=.TRUE.
         DO k=1,this%nrumpling
            SELECT CASE(this%zeropos(k))
               CASE(0)
                  exists_zero=.FALSE.
               CASE DEFAULT
                  ! do nothing
            END SELECT
         END DO
         SELECT CASE(exists_zero)
            CASE(.FALSE.)
               WRITE(0,*) "READ_CRP3D_details ERR: Manual grid does not have a rumpling point in the grid"
               CALL EXIT(1)
            CASE(.TRUE.)
               ! do nothing
         END SELECT
      !
      CASE DEFAULT
         WRITE(0,*) "READ_CRP3D_details ERR: wrong grid control keyword"
         WRITE(0,*) "You geave: ", control
         WRITE(0,*) "Implemented ones: MANU"
         WRITE(0,*) "Warning: case-sensitive"
         CALL EXIT(1)
   END SELECT

   CLOSE(runit)
   RETURN
END SUBROUTINE READ_CRP3D_details
!###########################################################
!# SUBROUTINE: SET_FOURIER_SYMMETRY_CRP3D
!###########################################################
!> @brief
!! Given a Fourier2d allocatable variable, allocates it with the correct surface
!! symmetry
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jul/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_FOURIER_SYMMETRY_CRP3D(this,interpolxy)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),INTENT(IN):: this
   CLASS(Fourier2d),ALLOCATABLE,INTENT(INOUT):: interpolxy
   ! Run section
   SELECT CASE(system_surface%getsymmlabel())
      CASE("p4mm")
         ALLOCATE(Fourierp4mm::interpolxy)
      CASE DEFAULT
         WRITE(0,*) "SET_FOURIER_SYMMETRY_CRP3D ERR: Incorrect surface symmlabel"
         WRITE(0,*) "Used: ",system_surface%getsymmlabel()
         WRITE(0,*) "Implemented ones: p4mm, p6mm"
         CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE SET_FOURIER_SYMMETRY_CRP3D
!###########################################################
!# SUBROUTINE: INITIALIZE_CRP3D 
!###########################################################
!> @brief
!! Specific implementation of initialize PES from input file
!
!> @param[in] filename - char(len=*),optional: input file. If system_inputfile
!!                       is allocated, it is ignored.
!> @param[in] tablename - char(len=*),optional: name of Lua table where the PES is
!!                        stored. By default, it has the value 'pes'.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014 
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INITIALIZE_CRP3D(this,filename,tablename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),INTENT(OUT):: this
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: filename,tablename
   ! Local variables
   CHARACTER(LEN=:),ALLOCATABLE:: auxstring
   ! Run section
   SELECT CASE(allocated(system_inputfile) .or. .not.present(filename))
      CASE(.TRUE.)
         auxstring=trim(system_inputfile)
      CASE(.FALSE.)
         auxstring=trim(filename)
   END SELECT
   SELECT CASE(present(tablename))
      CASE(.TRUE.) ! present tablename
         CALL this%READ(filename=trim(auxstring),tablename=trim(tablename))
      CASE(.FALSE.) ! not present tablename
         CALL this%READ(filename=trim(auxstring),tablename='pes')
   END SELECT
   CALL this%INTERPOL()
   RETURN
END SUBROUTINE INITIALIZE_CRP3D
!###########################################################
!# SUBROUTINE: GET_REPUL_CORRECTIONS_CRP3D 
!###########################################################
!> @brief
!! Collects interactions 
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_REPUL_CORRECTIONS_CRP3D(this,P,v,dvdz,dvdx,dvdy)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(3),INTENT(IN) :: P
   REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT):: dvdx,dvdy,dvdz ! corrections to the derivatives
   ! Local variables
   INTEGER(KIND=4) :: npairpots
   INTEGER(KIND=4) :: l,k ! counters
   REAL(KIND=8) :: aux1,aux2,aux3,aux4
   ! Run section
   npairpots=size(this%all_pairpots)
   FORALL(l=1:npairpots)
      v(l)=0.D0
      dvdz(l)=0.D0
      dvdx(l)=0.D0
      dvdy(l)=0.D0
   END FORALL
   DO l = 1, npairpots
      DO k = 0, this%max_order
         CALL INTERACTION_AENV(k,P,this%all_pairpots(l),this%dampfunc,aux1,aux2,aux3,aux4)
         v(l)=v(l)+aux1
         dvdz(l)=dvdz(l)+aux2
         dvdx(l)=dvdx(l)+aux3
         dvdy(l)=dvdy(l)+aux4
      END DO
   END DO
   RETURN
END SUBROUTINE GET_REPUL_CORRECTIONS_CRP3D
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_PAIRPOT 
!###########################################################
!> @brief
!! Gets the pairpot potential shifted to 0 (maximum interaction is at 0)
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_PAIRPOT(this,x,v,dvdu)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Pair_pot),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   REAL(KIND=8),INTENT(OUT) :: v,dvdu
   ! Run section
   SELECT CASE(X>this%z(this%n)-this%rumpling)
      CASE(.TRUE.)
         v=0.D0
         dvdu=0.D0
      CASE(.FALSE.)
         CALL this%interz%GET_V_AND_DERIVS(x,v,dvdu,this%rumpling)
   END SELECT
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_PAIRPOT
!###########################################################
!# FUNCTION: getrumpling_CRP3D 
!###########################################################
!> @brief 
!! Common get function. Gets rumpling value of parit potential "toptype"
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date MAy/2014
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getrumpling_CRP3D(this,toptype) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),INTENT(IN):: this
   INTEGER(KIND=4),INTENT(IN) :: toptype
   ! Run section
   getrumpling_CRP3D=this%all_pairpots(toptype)%rumpling
   RETURN
END FUNCTION getrumpling_CRP3D
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
SUBROUTINE READ_STANDARD_PAIRPOT(pairpot,filename)
      ! Initial declarations
   IMPLICIT NONE
   ! I/O variables ------------------------------
   CLASS(Pair_pot),INTENT(INOUT) :: pairpot
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables ----------------------------
   INTEGER :: i
   CHARACTER(LEN=14), PARAMETER :: routinename = "READ_PAIRPOT: "
   character(len=1024):: auxString
   ! Run section ---------------------------------
   pairpot%filename=filename
   OPEN(10,FILE=pairpot%filename,STATUS='OLD')
   READ(10,*) ! dummy line
   READ(10,*) ! dummy line
   READ(10,*) auxString
   pairpot%alias=trim(auxString)
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
   IMPLICIT NONE
   ! I/O variables
   CLASS(Sitio),INTENT(INOUT):: site
   CHARACTER(LEN=*),INTENT(IN):: filename
   ! Local variables
   INTEGER:: i ! counter
   CHARACTER(LEN=*),PARAMETER:: routinename = "READ_SITIO: "
   character(len=1024):: auxString
   !
   site%filename=filename
   OPEN (10,file=site%filename,status="old")
   READ(10,*) ! dummy line
   READ(10,*) ! dummy line
   READ(10,*) auxString
   site%alias=trim(auxString)
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
! SUBROUTINE:  GET_SYMMETRIZED_RAW_INPUT ############################
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
SUBROUTINE GET_SYMMETRIZED_RAW_INPUT(symmraw,zero,vtop,filename)
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
END SUBROUTINE GET_SYMMETRIZED_RAW_INPUT
!###########################################################
!# SUBROUTINE: READ_CRP3D
!###########################################################
!> @brief
!! Read CRP3D input file from Lua config file.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE READ_CRP3D(this,filename,tablename)
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),INTENT(OUT) :: this 
   CHARACTER(LEN=*),INTENT(IN) :: filename
   CHARACTER(LEN=*),INTENT(IN):: tablename
   ! Local variables
   INTEGER(KIND=4) :: n_pairpots,n_sites
   CHARACTER(LEN=1024),DIMENSION(:),ALLOCATABLE:: files_pairpots,files_sites
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: param
   INTEGER(KIND=4):: ierr
   INTEGER(KIND=4) :: i ! counter
   ! Lua-related variables
   TYPE(flu_State):: conf ! Lua state
   INTEGER(KIND=4):: pes_table,pairpot_table,sitio_table,dampfunc_table,param_table,fourier_table ! tables
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE:: subtables
   ! Auxiliar (dummy) variables
   INTEGER(KIND=4):: auxint
   CHARACTER(LEN=1024):: auxstring
   ! Parameters
   CHARACTER(LEN=*),PARAMETER :: routinename="READ_CRP3D: "
   ! HEY HO!, LET'S GO!! ------------------
   ! Open Lua file
   CALL OPEN_CONFIG_FILE(L=conf,filename=filename,ErrCode=ierr)
   SELECT CASE(ierr)
      CASE(0)
         ! do nothing
      CASE DEFAULT
         WRITE(0,*) "READ_CRP3D ERR: error reading Lua config file: ",filename
         CALL EXIT(1)
   END SELECT
   ! Open PES table
   CALL AOT_TABLE_OPEN(L=conf,thandle=pes_table,key=tablename)
   ! Set pestype (kind)
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='kind',val=auxstring)
   CALL this%SET_PESTYPE(trim(auxstring))
   SELECT CASE(trim(auxstring))
      CASE('CRP3D')
         ! do nothing
      CASE DEFAULT
         WRITE(0,*) "READ_CRP3D ERR: wrong type of PES. Expected: CRP3D. Encountered: "//trim(auxstring)
         CALL EXIT(1)
   END SELECT
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'Type of PES: '//trim(auxstring))
#endif
   ! Set alias (name)
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='name',val=auxstring)
   CALL this%SET_ALIAS(trim(auxstring))
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'PES Name: '//trim(auxstring))
#endif
   ! Set dimensions
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='dimensions',val=auxint)
   CALL this%SET_DIMENSIONS(auxint)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'PES dimensions: ',auxint)
#endif
   ! Set max environment
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='maxEnvironment',val=auxint)
   this%max_order=auxint
   ! Set pair potentials
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=pairpot_table,key='pairPotentials')
   n_pairpots=aot_table_length(L=conf,thandle=pairpot_table)
   ALLOCATE(files_pairpots(n_pairpots))
   ALLOCATE(this%all_pairpots(n_pairpots))
   DO i = 1, n_pairpots
      CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pairpot_table,pos=i,val=files_pairpots(i))
      CALL this%all_pairpots(i)%READ(trim(files_pairpots(i)))
   END DO
   CALL AOT_TABLE_CLOSE(L=conf,thandle=pairpot_table)
   ! Set damping function
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=dampfunc_table,key='dampFunction')
   CALL AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=dampfunc_table,key='kind',val=auxstring)
   SELECT CASE(trim(auxstring))
      CASE("Logistic")
         ALLOCATE(Logistic_func::this%dampfunc)
         ALLOCATE(param(2))
         ! open param table
         CALL AOT_TABLE_OPEN(L=conf,parent=dampfunc_table,thandle=param_table,key='param')
         auxint=aot_table_length(L=conf,thandle=param_table)
         SELECT CASE(auxint/=2) 
            CASE(.TRUE.)
               WRITE(0,*) "READ_CRP3D ERR: wrong number of parameters in pes.dampFunc.param table"
               CALL EXIT(1)
            CASE(.FALSE.)
               ! do nothing
         END SELECT
         DO i = 1, 2
            CALL AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=param_table,pos=i,val=param(i))
         END DO
         CALL this%dampfunc%READ(param)
         CALL AOT_TABLE_CLOSE(L=conf,thandle=param_table)
      CASE("fullCRP")
         ALLOCATE(One_func::this%dampfunc)
      CASE("fullRaw")
         ALLOCATE(Zero_func::this%dampfunc)
      CASE DEFAULT
         WRITE(0,*) "READ_CRP3D ERR: dampfunction keyword is not implemented"
         WRITE(0,*) "Implemented ones: Logistic, fullCRP, fullRaw"
         WRITE(0,*) "Case sensitive"
         CALL EXIT(1)
   END SELECT
   CALL AOT_TABLE_CLOSE(L=conf,thandle=dampfunc_table)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'Damping function used: '//trim(auxstring))
   SELECT CASE(allocated(param))
      CASE(.TRUE.)
         CALL VERBOSE_WRITE(routinename,'Damping function parameters: ',param(:))
      CASE(.FALSE.)
         ! do nothing
   END SELECT
#endif
   ! Set sitios
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=sitio_table,key='sitios')
   n_sites=aot_table_length(L=conf,thandle=sitio_table)
   ALLOCATE(files_sites(n_sites))
   ALLOCATE(this%all_sites(n_sites))
   DO i = 1, n_sites
      CALl AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=sitio_table,pos=i,val=files_sites(i))
      CALL this%all_sites(i)%READ(trim(files_sites(i))) 
   END DO
   CALL AOT_TABLE_CLOSE(L=conf,thandle=sitio_table)
   ! Read fourier Kpoints
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=fourier_table,key='fourierKpoints')
   auxint=aot_table_length(L=conf,thandle=fourier_table)
   allocate( this%irrepList(auxInt) )
   allocate( this%parityList(auxInt) )
   this%irrepList(:)='A1'
   this%parityList(:)='+'
   SELECT CASE(auxint/=n_sites)
      CASE(.TRUE.)
         WRITE(0,*) "READ_CRP3D ERR: dimension mismatch between fourierKpoints and number of sitios"
         WRITE(0,*) 'Number of Fourier Kpoints: ',auxint
         WRITE(0,*) 'Number of sitios: ',n_sites
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   ALLOCATE(this%klist(n_sites,2))
   ALLOCATE(subtables(n_sites))
   DO i = 1, n_sites
      CALL AOT_TABLE_OPEN(L=conf,parent=fourier_table,thandle=subtables(i),pos=i)
      CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=subtables(i),pos=1,val=auxint)
      this%klist(i,1)=auxint
      CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=subtables(i),pos=2,val=auxint)
      this%klist(i,2)=auxint
      CALL AOT_TABLE_CLOSE(L=conf,thandle=subtables(i))
   END DO
   CALL AOT_TABLE_CLOSE(L=conf,thandle=fourier_table)
   CALL AOT_TABLE_CLOSE(L=conf,thandle=pes_table)
   CALL CLOSE_CONFIG(conf)
   ! VERBOSE PRINT 
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Maximum environmental order: ",this%max_order)
   CALL VERBOSE_WRITE(routinename,"Number of pair potentials: ",n_pairpots)
   CALL VERBOSE_WRITE(routinename,"Pair potentials input files:")
   DO i = 1, n_pairpots
      CALL VERBOSE_WRITE(routinename,trim(files_pairpots(i)))
   END DO
   CALL VERBOSE_WRITE(routinename,"Number of sitios: ",n_sites)
   CALL VERBOSE_WRITE(routinename,"Sitios input files:")
   DO i = 1, n_sites
      CALL VERBOSE_WRITE(routinename,trim(files_sites(i)))
   END DO
   CALL VERBOSE_WRITE(routinename,"List of Kpoints for Fourier interpolation: ")
   DO i = 1, n_sites
      CALL VERBOSE_WRITE(routinename,this%klist(i,:))
   END DO
#endif
   RETURN
END SUBROUTINE READ_CRP3D
!#######################################################################
!# SUBROUTINE: EXTRACT_VASINT_CRP3D ######################################
!#######################################################################
!> @brief
!! Extracts the potential at long distances from every value stored in pairpots
!! and sites that belongs to @b this. Basically, this procedure sets
!! our zero potential value at the vacuum.
!
!> @param[in,out] this - A CRP3D subtype variable
!
!> @warning 
!! - The value of the potential at long distance should be the same in all
!!   pair potentials, read from files previously
!! - At least one pair potential should've defined in the CRP3D variable
!! - Obviously, @ this should contain reliable data
!! - Only modifies things inside @b interz
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 04/Feb/2014
!> @version 1.0
!------------------------------------------------------------------------
SUBROUTINE EXTRACT_VASINT_CRP3D(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),INTENT(INOUT) :: this
   ! Local variables
   INTEGER(KIND=4) :: npairpots, nsites
   INTEGER(KIND=4) :: i,j ! counters
   REAL(KIND=8) :: control_vasint
   CHARACTER(LEN=22) :: routinename="EXTRACT_VASINT_CRP3D: "
   ! Run section ------------------------
   npairpots=size(this%all_pairpots)
   control_vasint=this%all_pairpots(1)%vasint
   DO i = 1, npairpots
      IF (this%all_pairpots(1)%vasint/=control_vasint) THEN
         WRITE(0,*) "EXTRACT_VASINT_CRP3D ERR: Incoherences in vasint values found"
         WRITE(0,*) "EXTRACT_VASINT_CRP3D ERR: vasint's value at pairpot",1,control_vasint
         WRITE(0,*) "EXTRACT_VASINT_CRP3D ERR: vasint's value at pairpot",i,control_vasint
         CALL EXIT(1)
      END IF
      DO j = 1, this%all_pairpots(i)%n
         this%all_pairpots(i)%interz%f(j)=this%all_pairpots(i)%interz%f(j)-this%all_pairpots(i)%vasint
      END DO
#ifdef DEBUG
      CALL DEBUG_WRITE(routinename,"Vasint extracted from pair potential ",i)
#endif
   END DO
   nsites=size(this%all_sites)
   DO i = 1, nsites
      DO j = 1, this%all_sites(i)%n
         this%all_sites(i)%interz%f(j)=this%all_sites(i)%interz%f(j)-this%all_pairpots(1)%vasint
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
!> @param[in,out] this - CRP3D PES variable to be smoothed (only it's sites)
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
SUBROUTINE SMOOTH_CRP3D(this)
   ! Initial declaraitons
   IMPLICIT NONE
   CLASS(CRP3D),INTENT(INOUT):: this
   ! Local variables
   REAL(KIND=8),DIMENSION(3):: A
   INTEGER(KIND=4):: i,j ! counters
   INTEGER(KIND=4):: npairpots,nsites
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: v,dvdzr,dummy
   CHARACTER(LEN=*),PARAMETER:: routinename="SMOOTH_CRP3D: "
   ! Run section ----------
   nsites = size(this%all_sites)
   npairpots = size(this%all_pairpots)
   ALLOCATE(v(npairpots))
   ALLOCATE(dvdzr(npairpots))
   ALLOCATE(dummy(npairpots))
   DO i = 1, nsites ! loop over sites
      A(1) = this%all_sites(i)%x
      A(2) = this%all_sites(i)%y
      DO j = 1, this%all_sites(i)%n ! loop over pairs v,z
         A(3)=this%all_sites(i)%z(j)
         CALL this%GET_REPUL_CORRECTIONS(A,v,dvdzr,dummy,dummy)
         this%all_sites(i)%interz%f(j)=this%all_sites(i)%interz%f(j)-sum(v)
         IF (j.EQ.1) THEN
            this%all_sites(i)%dz1=this%all_sites(i)%dz1-sum(dvdzr) ! correct first derivative
         ELSE IF (j.EQ.this%all_sites(i)%n) THEN
            this%all_sites(i)%dz2=this%all_sites(i)%dz2-sum(dvdzr) ! correct first derivative
         END IF
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
!> @param[in,out] this - CRP3D PES to be interpolated 
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 05/Feb/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INTERPOL_Z_CRP3D(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),INTENT(INOUT) :: this
   ! Local variables
   INTEGER(KIND=4) :: nsites,npairpots
   REAL(KIND=8) :: dz1,dz2
   INTEGER(KIND=4) :: i ! counters
   ! Run secton------------------------
   nsites=size(this%all_sites)
   npairpots=size(this%all_pairpots)
   CALL this%EXTRACT_VASINT()
   DO i = 1, npairpots ! loop pairpots
      dz1=this%all_pairpots(i)%dz1
      dz2=this%all_pairpots(i)%dz2
      CALL this%all_pairpots(i)%interz%INTERPOL(dz1,1,dz2,1)
   END DO
   CALL this%SMOOTH()
   DO i = 1, nsites
      dz1=this%all_sites(i)%dz1
      dz2=this%all_sites(i)%dz2
      CALL this%all_sites(i)%interz%INTERPOL(dz1,1,dz2,1) 
   END DO
   RETURN
END SUBROUTINE INTERPOL_Z_CRP3D
!############################################################
! SUBROUTINE: RAWINTERPOL_Z_CRP3D ################################
!############################################################
!> @brief
!! Interpolates all sitios and pairpot potentials that belong
!! to an specific CRP3D PES. Sitio potentials aren't smoothed
!
!> @param[in,out] this - CRP3D PES to be interpolated 
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE RAWINTERPOL_Z_CRP3D(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D) :: this
   ! Local variables
   INTEGER(KIND=4) :: nsites,npairpots
   REAL(KIND=8) :: dz1,dz2
   INTEGER(KIND=4) :: i ! counters
   ! Run secton------------------------
   nsites=size(this%all_sites)
   npairpots=size(this%all_pairpots)
   CALL this%EXTRACT_VASINT()
   DO i = 1, npairpots ! loop pairpots
      dz1=this%all_pairpots(i)%dz1
      dz2=this%all_pairpots(i)%dz2
      CALL this%all_pairpots(i)%interz%INTERPOL(dz1,1,dz2,1)
   END DO
   DO i = 1, nsites
      dz1=this%all_sites(i)%dz1
      dz2=this%all_sites(i)%dz2
      CALL this%all_sites(i)%interz%INTERPOL(dz1,1,dz2,1) 
   END DO
   RETURN
END SUBROUTINE RAWINTERPOL_Z_CRP3D
!##################################################################
!# SUBROUTINE: INTERACTION_AP #####################################
!##################################################################
!> @brief
!! Calculates the interaction between two given 3D Points.
!! The potential is extracted from @b pairpot, which should've been interpolated before
!
!> @param[in] A, P - Pair of 3D points that are interacting through a @b pairpot @b potential
!> @param[in] pairpot - Pair potential. Source of @f$V(r)@f$, where @b r is the distance
!!                      between @b A and @b P. @f$V(r)@f$ is just a shifted version of
!!                      @f$V(z)@f$.
!> @param[out] interac - Actual interaction between @b A and @b P
!> @param[out] dvdz_corr - Corrections to first derivatives: 
!!                         @f$\frac {\partial V(r)}{\partial z}@f$
!> @param[out] dvdx_corr - Corrections to first derivatives: 
!!                         @f$\frac {\partial V(r)}{\partial x}@f$
!> @param[out] dvdy_corr - Corrections to first derivatives: 
!!                         @f$\frac{\partial V(r)}{\partial y}@f$
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
SUBROUTINE INTERACTION_AP(A,P,pairpot,dampfunc,interac,dvdz_corr,dvdx_corr,dvdy_corr)
   IMPLICIT NONE
   ! I/O VAriables --------------------------------------------
   REAL(KIND=8),DIMENSION(3),INTENT(IN):: A, P
   TYPE(Pair_pot),INTENT(IN):: pairpot
   CLASS(Function1d),INTENT(IN):: dampfunc
   REAL(KIND=8),INTENT(OUT):: interac,dvdz_corr,dvdx_corr,dvdy_corr
   ! Local variables ------------------------------------------
   REAL(KIND=8):: r ! distance
   REAL(KIND=8):: v,pre
   REAL(KIND=8):: aux ! dv/dr
   CHARACTER(LEN=*),PARAMETER:: routinename = "INTERACTION_AP: "
   ! GABBA, GABBA HEY! ----------------------------------------
   ! Find the distance between A and P, in a.u.
   r=dsqrt((A(1)-P(1))**2.D0+(A(2)-P(2))**2.D0+(A(3)-P(3))**2.D0)
   CALL pairpot%GET_V_AND_DERIVS(r,v,aux)
   interac=v*dampfunc%getvalue(r)
   pre=aux*dampfunc%getvalue(r)+v*dampfunc%getderiv(r)
   dvdz_corr=pre*(A(3)-P(3))/r
   dvdx_corr=pre*(A(1)-P(1))/r
   dvdy_corr=pre*(A(2)-P(2))/r
   RETURN
END SUBROUTINE INTERACTION_AP
!##################################################################
!# SUBROUTINE: INTERACTION_AENV ###################################
!##################################################################
!> @brief
!! Calculates the complete interaction between a given point @b A with an
!! environment of order @b n. Depending on the surface, it can use a hexagonal
!! environment or an octagonal one. If a new surface symmetry is implemented, this routine should be
!! edited as well.
!
!> @param[in] n - Environment order
!> @param[in] A - Point in space, cartesian coordinates
!> @param[in] pairpot - Pair potential which defines interaction potential
!> @param[out] interac - Total interaction with the environment defined
!> @param[out] dvdx_term - 
!> @param[out] dvdy_term - 
!> @param[out] dvdz_term - 
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Feb/2014;Jun/2014
!> @version 2.0
!------------------------------------------------------------------
SUBROUTINE INTERACTION_AENV(n,A,pairpot,dampfunc,interac,dvdz_term,dvdx_term,dvdy_term)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   INTEGER,INTENT(IN):: n
   REAL(KIND=8),DIMENSION(3),INTENT(IN):: A
   TYPE(Pair_pot),INTENT(IN):: pairpot
   CLASS(Function1d),INTENT(IN):: dampfunc
   REAL(KIND=8),INTENT(OUT):: interac, dvdz_term, dvdx_term, dvdy_term
   ! Run section
   SELECT CASE(system_surface%getsymmlabel())
      CASE("p4mm")
         CALL INTERACTION_AENV_OCTA(n,A,pairpot,dampfunc,interac,dvdz_term,dvdx_term,dvdy_term)
      CASE("p6mm")
         CALL INTERACTION_AENV_HEXA(n,A,pairpot,dampfunc,interac,dvdz_term,dvdx_term,dvdy_term)
      CASE DEFAULT
         WRITE(0,*) "INTERACTION_AENV ERR: wrong surface symmlabel or it's not been implemented yet"
         WRITE(0,*) "Surface reads: ",system_surface%getsymmlabel()
         WRITE(0,*) "Implemented ones: p4mm, p6mm"
         CALL EXIT(1)      
   END SELECT
   RETURN
END SUBROUTINE INTERACTION_AENV
SUBROUTINE INTERACTION_AENV_OCTA(n,A,pairpot,dampfunc,interac,dvdz_term,dvdx_term,dvdy_term)
   IMPLICIT NONE
   ! I/O VAriables ------------------------------------------
   INTEGER,INTENT(IN) :: n
   REAL(KIND=8),DIMENSION(3),INTENT(IN) :: A
   TYPE(Pair_pot),INTENT(IN) :: pairpot
   CLASS(Function1d),INTENT(IN) :: dampfunc
   REAL(KIND=8),INTENT(OUT) :: interac, dvdz_term, dvdx_term, dvdy_term
   ! Local variables ----------------------------------------
   REAL(KIND=8),DIMENSION(3) :: P
   REAL(KIND=8),DIMENSION(3) :: ghost_A ! A in cartesians, but inside unitcell
   REAL(KIND=8),DIMENSION(2) :: aux
   REAL(KIND=8) :: dummy1, dummy2, dummy3, dummy4
   REAL(KIND=8) :: atomx, atomy
   INTEGER :: pairid
   INTEGER :: i, k ! Counters
   CHARACTER(LEN=18), PARAMETER :: routinename = "INTERACTION_AENV: "
   ! SUSY IS A HEADBANGER !!!! -------------------
   ! Defining some aliases to make the program simpler:
   pairid = pairpot%id
   interac=0.D0
   dvdz_term=0.D0
   dvdx_term=0.D0
   dvdy_term=0.D0
   ! ghost A definition
   ghost_A(1:2)=system_surface%project_unitcell(A(1:2))
   ghost_A(3)=A(3)

   SELECT CASE(n)
      CASE(0)
         DO i=1, system_surface%atomtype(pairid)%n
            P(:)=system_surface%atomtype(pairid)%atom(i,:)
            CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac=interac+dummy1
            dvdz_term=dvdz_term+dummy2
            dvdx_term=dvdx_term+dummy3
            dvdy_term=dvdy_term+dummy4
         END DO
         RETURN

      CASE(1 :)
         DO i=1, system_surface%atomtype(pairid)%n
            atomx=system_surface%atomtype(pairid)%atom(i,1)
            atomy=system_surface%atomtype(pairid)%atom(i,2)
            P(3)=system_surface%atomtype(pairid)%atom(i,3)
            DO k= -n,n
               aux(1)=dfloat(n)
               aux(2)=dfloat(k)
               aux=system_surface%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac=interac+dummy1
               dvdz_term=dvdz_term+dummy2
               dvdx_term=dvdx_term+dummy3
               dvdy_term=dvdy_term+dummy4
               !
               aux(1)=dfloat(-n)
               aux(2)=dfloat(k)
               aux=system_surface%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
            END DO
            DO k= -n+1, n-1
               aux(1)=dfloat(k)
               aux(2)=dfloat(n)
               aux=system_surface%surf2cart(aux)
               P(1) =atomx+aux(1)
               P(2) =atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
               !
               aux(1)=dfloat(k)
               aux(2)=dfloat(-n)
               aux=system_surface%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
            END DO
         END DO
         RETURN

      CASE DEFAULT
         WRITE(0,*) "INTERACTION_AENV ERR: Wrong environment order."
         CALL EXIT(1)
   END SELECT
END SUBROUTINE INTERACTION_AENV_OCTA
SUBROUTINE INTERACTION_AENV_HEXA(n,A,pairpot,dampfunc,interac,dvdz_term,dvdx_term,dvdy_term)
   IMPLICIT NONE
   ! I/O VAriables ------------------------------------------
   INTEGER,INTENT(IN) :: n
   REAL(KIND=8),DIMENSION(3),INTENT(IN) :: A
   TYPE(Pair_pot),INTENT(IN) :: pairpot
   CLASS(Function1d),INTENT(IN) :: dampfunc
   REAL(KIND=8),INTENT(OUT) :: interac, dvdz_term, dvdx_term, dvdy_term
   ! Local variables ----------------------------------------
   REAL(KIND=8),DIMENSION(3) :: P
   REAL(KIND=8),DIMENSION(3) :: ghost_A ! A in cartesians, but inside unitcell
   REAL(KIND=8),DIMENSION(2) :: aux
   REAL(KIND=8) :: dummy1, dummy2, dummy3, dummy4
   REAL(KIND=8) :: atomx, atomy
   INTEGER :: pairid
   INTEGER :: i, k ! Counters
   CHARACTER(LEN=18), PARAMETER :: routinename = "INTERACTION_AENV: "
   ! SUSY IS A HEADBANGER !!!! -------------------
   ! Defining some aliases to make the program simpler:
   pairid = pairpot%id
   interac=0.D0
   dvdz_term=0.D0
   dvdx_term=0.D0
   dvdy_term=0.D0
   ! ghost A definition
   ghost_A(1:2)=system_surface%project_unitcell(A(1:2))
   ghost_A(3)=A(3)
   SELECT CASE(n)
      CASE(0)
         DO i=1, system_surface%atomtype(pairid)%n
            P(:)=system_surface%atomtype(pairid)%atom(i,:)
            CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac=interac+dummy1
            dvdz_term=dvdz_term+dummy2
            dvdx_term=dvdx_term+dummy3
            dvdy_term=dvdy_term+dummy4
         END DO
         RETURN
      CASE(1)
         DO i=1, system_surface%atomtype(pairid)%n
            atomx=system_surface%atomtype(pairid)%atom(i,1)
            atomy=system_surface%atomtype(pairid)%atom(i,2)
            P(3)=system_surface%atomtype(pairid)%atom(i,3)
            k=0
            !
            aux(1)=dfloat(n)
            aux(2)=dfloat(-k)
            aux=system_surface%surf2cart(aux)
            P(1) =atomx+aux(1)
            P(2) =atomy+aux(2)
            CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
            !
            aux(1)=dfloat(-n)
            aux(2)=dfloat(k)
            aux=system_surface%surf2cart(aux)
            P(1)=atomx+aux(1)
            P(2)=atomy+aux(2)
            CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
            !
            aux(1)=dfloat(-k)
            aux(2)=dfloat(n)
            aux=system_surface%surf2cart(aux)
            P(1)=atomx+aux(1)
            P(2)=atomy+aux(2)
            CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
            !
            aux(1)=dfloat(k)
            aux(2)=dfloat(-n)
            aux=system_surface%surf2cart(aux)
            P(1)=atomx+aux(1)
            P(2)=atomy+aux(2)
            CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
            !
            aux(1)=dfloat(-n)
            aux(2)=dfloat(n)
            aux=system_surface%surf2cart(aux)
            P(1)=atomx+aux(1)
            P(2)=atomy+aux(2)
            CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
            !
            aux(1)=dfloat(n)
            aux(2)=dfloat(-n)
            aux=system_surface%surf2cart(aux)
            P(1)=atomx+aux(1)
            P(2)=atomy+aux(2)
            CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
         END DO
         RETURN
      CASE(2 :)
         DO i=1, system_surface%atomtype(pairid)%n
            atomx=system_surface%atomtype(pairid)%atom(i,1)
            atomy=system_surface%atomtype(pairid)%atom(i,2)
            P(3)=system_surface%atomtype(pairid)%atom(i,3)
            DO k= 1,n-1
               aux(1)=dfloat(k)
               aux(2)=dfloat(n-k)
               aux=system_surface%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac=interac+dummy1
               dvdz_term=dvdz_term+dummy2
               dvdx_term=dvdx_term+dummy3
               dvdy_term=dvdy_term+dummy4
               !
               aux(1)=dfloat(-k)
               aux(2)=dfloat(k-n)
               aux=system_surface%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
            END DO
            DO k= 0, n-1
               aux(1)=dfloat(n)
               aux(2)=dfloat(-k)
               aux=system_surface%surf2cart(aux)
               P(1) =atomx+aux(1)
               P(2) =atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
               !
               aux(1)=dfloat(-n)
               aux(2)=dfloat(k)
               aux=system_surface%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
               !
               aux(1)=dfloat(-k)
               aux(2)=dfloat(n)
               aux=system_surface%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
               !
               aux(1)=dfloat(k)
               aux(2)=dfloat(-n)
               aux=system_surface%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
            END DO
            !
            aux(1)=dfloat(-n)
            aux(2)=dfloat(n)
            aux=system_surface%surf2cart(aux)
            P(1)=atomx+aux(1)
            P(2)=atomy+aux(2)
            CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
            !
            aux(1)=dfloat(n)
            aux(2)=dfloat(-n)
            aux=system_surface%surf2cart(aux)
            P(1)=atomx+aux(1)
            P(2)=atomy+aux(2)
            CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
         END DO
         RETURN

      CASE DEFAULT
         WRITE(0,*) "INTERACTION_AENV ERR: Wrong environment order."
         CALL EXIT(1)
   END SELECT
END SUBROUTINE INTERACTION_AENV_HEXA
!############################################################
!# SUBROUTINE: GET_V_AND_DERIVS_CRP3D #########################
!############################################################
!> @brief
!! Subroutine that calculates the 3D potential for a point A and
!! its derivatives in cartesian coordinates.
!
!> @param[in] this - CRP3D PES
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
SUBROUTINE GET_V_AND_DERIVS_CRP3D(this,X,v,dvdu,errCode)
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),TARGET,INTENT(IN):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN):: X
   REAL(KIND=8),INTENT(OUT):: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT):: dvdu
   integer(kind=1),optional,intent(out):: errCode
   ! Local variables
   INTEGER(KIND=4):: nsites,npairpots
   CLASS(Fourier2d),ALLOCATABLE:: interpolxy
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: pot,dvdz,dvdx,dvdy
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: potarr
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: f,derivarr ! arguments to the xy interpolation
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: xy ! arguments to the xy interpolation
   INTEGER :: i ! counters
   ! Pointers
	REAL(KIND=8),POINTER:: zmax
   CHARACTER(LEN=*),PARAMETER:: routinename="GET_V_AND_DERIVS_CRP3D: "
   zmax => this%all_sites(1)%z(this%all_sites(1)%n)
   npairpots = size(this%all_pairpots)
   nsites = size(this%all_sites)
   ! GABBA, GABBA HEY! ----------------------
   SELECT CASE(X(3)>zmax)
      CASE(.TRUE.)
         dvdu = (/0.D0,0.D0,0.D0/)
         v = 0.D0
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   ! Initialization section
   ALLOCATE(pot(npairpots))
   ALLOCATE(dvdz(npairpots))
   ALLOCATE(dvdx(npairpots))
   ALLOCATE(dvdy(npairpots))
   CALL this%GET_REPUL_CORRECTIONS(X,pot,dvdz,dvdx,dvdy)
   v=sum(pot)
   dvdu(1)=sum(dvdx)
   dvdu(2)=sum(dvdy)
   dvdu(3)=sum(dvdz)
   ! Now, we have all the repulsive interaction and corrections to the derivarives
   ! stored in v(:) and dvdu(:) respectively.
   ! Let's get v and derivatives from xy interpolation of the corrugationless function
   ! f(1,i) ==> v values interpolated for Z at site "i"
   ! f(2,i) ==> dvdz values interpolated for Z at site "i"
   ALLOCATE(f(2,nsites))
   ALLOCATE(xy(nsites,2))
   ! potarr(1) ==> v interpolated at X,Y for a given Z
   ! potarr(2) ==> dvdz interpolated at X,Y for a given Z
   ALLOCATE(potarr(2))
   ! derivarr(1,1:2) ==> dvdx, dvdy interpolated at X,Y
   ! derivarr(2,1:2) ==> d(dvdz)/dx, d(dvdz)/dy interpolated at X,Y. This is not needed
   ALLOCATE(derivarr(2,2))
   DO i=1,nsites
      xy(i,1)=this%all_sites(i)%x
      xy(i,2)=this%all_sites(i)%y
      CALL this%all_sites(i)%interz%GET_V_AND_DERIVS(X(3),f(1,i),f(2,i))
   END DO
   CALL this%SET_FOURIER_SYMMETRY(interpolxy)
   CALL interpolxy%READ( xy=xy,f=f,kList=this%klist,irrepList=this%irrepList(:),parityList=this%parityList )
   call interpolXY%initializeTerms()
   CALL interpolxy%INTERPOL()
   CALL interpolxy%GET_F_AND_DERIVS(X,potarr,derivarr)
   ! Corrections from the smoothing procedure
   v = v + potarr(1)
   dvdu(1)=dvdu(1)+derivarr(1,1)
   dvdu(2)=dvdu(2)+derivarr(1,2)
   dvdu(3)=dvdu(3)+potarr(2)
   call interpolXY%cleanAll()
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_CRP3D
!############################################################
!# SUBROUTINE: GET_V_AND_DERIVS_CORRECTION_CRP3D ############
!############################################################
!> @brief
!! Subroutine that calculates the correction to the 3D PES for a point A and
!! its derivatives in cartesian coordinates.
!
!> @param[in] this - CRP3D PES
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
SUBROUTINE GET_V_AND_DERIVS_CORRECTION_CRP3D(this,X,v,dvdu)
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),TARGET,INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(3), INTENT(IN) :: X
   REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(3),INTENT(OUT) :: dvdu
   ! Local variables
   INTEGER(KIND=4) :: npairpots
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: pot,dvdz,dvdx,dvdy
   INTEGER :: i ! counters
   ! Pointers
	REAL(KIND=8), POINTER :: zmax
   CHARACTER(LEN=24),PARAMETER :: routinename="GET_V_AND_DERIVS_CRP3D: "
   zmax => this%all_sites(1)%interz%x(this%all_sites(1)%n)
   npairpots = size(this%all_pairpots)
   SELECT CASE(size(v)==npairpots+1)
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         WRITE(*,*) "GET_V_AND_DERIVS_CORRECTION_CRP3D ERR: wrong number of dimensions array v"
         CALL EXIT(1)
   END SELECT
   ! GABBA, GABBA HEY! ----------------------
   SELECT CASE(X(3)>zmax)
      CASE(.TRUE.)
         dvdu = (/0.D0,0.D0,0.D0/)
         FORALL(i=1:npairpots+1) v(i) = 0.D0
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   ! Initializing variables
   ALLOCATE(pot(npairpots))
   ALLOCATE(dvdz(npairpots))
   ALLOCATE(dvdx(npairpots))
   ALLOCATE(dvdy(npairpots))
   FORALL(i=1:npairpots+1) v(i) = 0.D0
   ! Compute
   CALL this%GET_REPUL_CORRECTIONS(X,pot,dvdz,dvdx,dvdy)
   v(1)=sum(pot)
   FORALL(i=2:npairpots+1) v(i)=pot(i-1)
   dvdu(1)=sum(dvdx)
   dvdu(2)=sum(dvdy)
   dvdu(3)=sum(dvdz)
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_CORRECTION_CRP3D
!############################################################
!# SUBROUTINE: GET_V_AND_DERIVS_SMOOTH_CRP3D #########################
!############################################################
!> @brief
!! Subroutine that calculates the 3D potential for a point A and
!! its derivatives in cartesian coordinates without any correction to
!! the smooth potential
!
!> @param[in] this - CRP3D PES
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
!> @date Jun/2014
!> @version 1.0
!------------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_SMOOTH_CRP3D(this,X,v,dvdu)
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),TARGET,INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(3), INTENT(IN) :: X
   REAL(KIND=8),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(3),INTENT(OUT) :: dvdu
   ! Local variables
   CLASS(Fourier2d),ALLOCATABLE:: interpolxy
   INTEGER(KIND=4) :: nsites,npairpots
   REAL(KIND=8),DIMENSION(3) :: deriv
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: potarr
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: f,derivarr ! arguments to the xy interpolation
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: xy ! arguments to the xy interpolation
   INTEGER :: i ! counters
   ! Pointers
	REAL(KIND=8), POINTER :: zmax
   CHARACTER(LEN=24),PARAMETER :: routinename="GET_V_AND_DERIVS_CRP3D: "
   zmax => this%all_sites(1)%interz%x(this%all_sites(1)%n)
   npairpots = size(this%all_pairpots)
   nsites = size(this%all_sites)
   ! GABBA, GABBA HEY! ----------------------
   IF (X(3).GT.zmax) THEN
      dvdu = 0.D0
      v = 0.D0
      RETURN
   END IF
   !
   v = 0.D0 ! Initialization value
   FORALL(i=1:3) 
      dvdu(i) = 0.D0 ! Initialization value
      deriv(i) = 0.D0 ! Initialization value
   END FORALL 
   ! Let's get v and derivatives from xy interpolation of the corrugationless function
   ! f(1,i) ==> v values interpolated for Z at site "i"
   ! f(2,i) ==> dvdz values interpolated for Z at site "i"
   ALLOCATE(f(2,nsites))
   ALLOCATE(xy(nsites,2))
   ! potarr(1) ==> v interpolated at X,Y for a given Z
   ! potarr(2) ==> dvdz interpolated at X,Y for a given Z
   ALLOCATE(potarr(2))
   ! derivarr(1,1:2) ==> dvdx, dvdy interpolated at X,Y
   ! derivarr(2,1:2) ==> d(dvdz)/dx, d(dvdz)/dy interpolated at X,Y. It is not needed
   ALLOCATE(derivarr(2,2))
   DO i=1,nsites
      xy(i,1)=this%all_sites(i)%x
      xy(i,2)=this%all_sites(i)%y
      CALL this%all_sites(i)%interz%GET_V_AND_DERIVS(X(3),f(1,i),f(2,i))
   END DO
   CALL this%SET_FOURIER_SYMMETRY(interpolxy)
   CALL interpolxy%READ( xy=xy,f=f,kList=this%klist,irrepList=this%irrepList(:),parityList=this%parityList )
   call interpolXY%initializeTerms()
   CALL interpolxy%INTERPOL()
   CALL interpolxy%GET_F_AND_DERIVS(X,potarr,derivarr)
   ! Corrections from the smoothing procedure
   v = v + potarr(1)
   dvdu(1)=dvdu(1)+derivarr(1,1)
   dvdu(2)=dvdu(2)+derivarr(1,2)
   dvdu(3)=dvdu(3)+potarr(2)
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_SMOOTH_CRP3D
!############################################################
!# FUNCTION: getpot_crp3d ###################################
!############################################################
!> @brief
!! Subroutine that calculates the 3D potential for a point A
!
!> @param[in] this - CRP3D PES
!> @param[in] X - Point in space to calculate the potential. Cartesian's
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
REAL(KIND=8) FUNCTION getpot_crp3d(this,X)
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),TARGET,INTENT(IN) :: this
	REAL(KIND=8),DIMENSION(3), INTENT(IN) :: X
   ! Local variables
   CLASS(Fourier2d),ALLOCATABLE:: interpolxy
   REAL(KIND=8):: v
   INTEGER(KIND=4) :: nsites,npairpots
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: pot,dummy
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: potarr
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: f,derivarr ! arguments to the xy interpolation
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: xy ! arguments to the xy interpolation
   INTEGER :: i ! counters
   ! Pointers
   REAL(KIND=8), POINTER :: zmax
   CHARACTER(LEN=24),PARAMETER :: routinename="GET_V_AND_DERIVS_CRP3D: "
   zmax => this%all_sites(1)%interz%x(this%all_sites(1)%n)
   npairpots = size(this%all_pairpots)
   nsites = size(this%all_sites)
   ! GABBA, GABBA HEY! ----------------------
   ALLOCATE(pot(npairpots))
   ALLOCATE(dummy(npairpots))
   SELECT CASE(X(3)>zmax)
      CASE(.TRUE.)
         getpot_crp3d=0.D0
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   !
   CALL this%GET_REPUL_CORRECTIONS(X,pot,dummy,dummy,dummy)
   v=sum(pot)
   ! Now, we have all the repulsive interaction and corrections to the derivarives
   ! stored in v(:) and dvdu(:) respectively.
   ! Let's get v and derivatives from xy interpolation of the corrugationless function
   ALLOCATE(f(2,nsites))
   ALLOCATE(xy(nsites,2))
   ALLOCATE(potarr(2))
   ALLOCATE(derivarr(2,2))
   DO i=1,nsites
      xy(i,1)=this%all_sites(i)%x
      xy(i,2)=this%all_sites(i)%y
      CALL this%all_sites(i)%interz%GET_V_AND_DERIVS(X(3),f(1,i),f(2,i))
   END DO
   CALL this%SET_FOURIER_SYMMETRY(interpolxy)
   CALL interpolxy%READ( xy=xy,f=f,kList=this%klist,irrepList=this%irrepList(:),parityList=this%parityList )
   call interpolXY%initializeTerms()
   CALL interpolxy%INTERPOL()
   CALL interpolxy%GET_F_AND_DERIVS(X,potarr,derivarr)
   ! Corrections from the smoothing procedure
   getpot_crp3d=v+potarr(1)
   RETURN
END FUNCTION getpot_crp3d
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
! SUBROUTINE: PLOT_XYMAP_CRP3D
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (X,Y) of the PES. 
!
!> @param[in] this - CRP3D PES used
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
SUBROUTINE PLOT_XYMAP_CRP3D(this,filename,init_xyz,nxpoints,nypoints,Lx,Ly)
   IMPLICIT NONE
   CLASS(CRP3D),INTENT(IN) :: this
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
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1),r(2),v,dvdu(:)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v,dvdu(:)
   END DO
   r(2) = ymax
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v,dvdu(:)
   ! inpoints in XY
   DO i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      r(3) = z
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v,dvdu(:)
      DO j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         CALL this%GET_V_AND_DERIVS(r,v,dvdu)
         WRITE(11,*) r(1), r(2), v,dvdu(:)
      END DO
      r(2) = ymax
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v,dvdu(:)
   END DO
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   r(3) = z
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v,dvdu(:)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v,dvdu(:)
   END DO
   r(2) = ymax
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v,dvdu(:)
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_XYMAP_CRP3D
!#######################################################################
! SUBROUTINE: PLOT_XYMAP_CORRECTION_CRP3D
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (X,Y) of the corrections
!! to the PES
!
!> @param[in] this - CRP3D PES used
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
SUBROUTINE PLOT_XYMAP_CORRECTION_CRP3D(this,filename,init_xyz,nxpoints,nypoints,Lx,Ly)
   IMPLICIT NONE
   CLASS(CRP3D),INTENT(IN) :: this
   REAL*8,DIMENSION(3),INTENT(IN) :: init_xyz ! Initial position to start the scan (in a.u.)
   INTEGER,INTENT(IN) :: nxpoints, nypoints ! number of points in XY plane
   CHARACTER(LEN=*),INTENT(IN) :: filename ! filename
   REAL*8,INTENT(IN) :: Lx ! Length of X axis 
   REAL*8,INTENT(IN) :: Ly ! Length of X axis 
   ! Local variables
   REAL*8 :: xmin, ymin, xmax, ymax, z
   REAL*8, DIMENSION(3) :: r, dvdu
   REAL*8 :: xdelta, ydelta
   INTEGER :: xinpoints, nxdelta, npairpots
   INTEGER :: yinpoints, nydelta
   INTEGER :: i, j ! counters
   REAL*8,DIMENSION(:),ALLOCATABLE :: v ! potential
	! GABBA, GABBA HEY! ---------
   npairpots=size(this%all_pairpots)
   ALLOCATE(v(npairpots+1))
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
   CALL this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v(:)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v(:)
   END DO
   r(2) = ymax
   CALL this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v(:)
   ! inpoints in XY
   DO i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      r(3) = z
      CALL this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v(:)
      DO j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         CALL this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
         WRITE(11,*) r(1), r(2), v(:)
      END DO
      r(2) = ymax
      CALL this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v(:)
   END DO
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   r(3) = z
   CALL this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v(:)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v(:)
   END DO
   r(2) = ymax
   CALL this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v(:)
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_XYMAP_CORRECTION_CRP3D
!#######################################################################
! SUBROUTINE: PLOT_XYMAP_SMOOTH_CRP3D
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (X,Y) of the PES without
!! corrections to the potential
!
!> @param[in] this - CRP3D PES used
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
SUBROUTINE PLOT_XYMAP_SMOOTH_CRP3D(this,filename,init_xyz,nxpoints,nypoints,Lx,Ly)
   IMPLICIT NONE
   CLASS(CRP3D),INTENT(IN) :: this
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
   CALL this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v
   END DO
   r(2) = ymax
   CALL this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v
   ! inpoints in XY
   DO i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      r(3) = z
      CALL this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v
      DO j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         CALL this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
         WRITE(11,*) r(1), r(2), v
      END DO
      r(2) = ymax
      CALL this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v
   END DO
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   r(3) = z
   CALL this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v
   END DO
   r(2) = ymax
   CALL this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_XYMAP_SMOOTH_CRP3D
!#######################################################################
! SUBROUTINE: PLOT_DIRECTION1D_CRP3D ###################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. To define 
!! the direction, the angle alpha is given. 
!
!> @param[in] this - CRP3D PES used
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
SUBROUTINE PLOT_DIRECTION1D_CRP3D(this,filename,npoints,angle,z,L)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP3D),INTENT(IN) :: this
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
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) s, v , dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   ! cycle for inpoints
   DO i=1, inpoints
      r(1)=xmin+(DFLOAT(i)*delta)*DCOS(alpha)
      r(2)=ymin+(DFLOAT(i)*delta)*DSIN(alpha)
      s = DSQRT(r(1)**2.D0+r(2)**2.D0)
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) s, v , dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   END DO
   ! Final value
   r(1) = xmax
   r(2) = ymax
   s = DSQRT(r(1)**2.D0+r(2)**2.D0)
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) s, v, dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT_DIRECTION1D_CRP3D
!#######################################################################
! SUBROUTINE: PLOT_DIRECTION1D_CORRECTION_CRP3D ###################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the corrections to the PES. To define 
!! the direction, the angle alpha is given. 
!
!> @param[in] this - CRP3D PES used
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
SUBROUTINE PLOT_DIRECTION1D_CORRECTION_CRP3D(this,filename,npoints,angle,z,L)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP3D),INTENT(IN) :: this
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*), INTENT(IN) :: filename
   REAL*8, INTENT(IN) :: z, angle
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta,npairpots
   REAL*8 :: delta,L,s,alpha
   REAL*8 :: xmax, xmin, ymax, ymin 
   REAL*8, DIMENSION(3) :: r, dvdu
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: v
   INTEGER :: i ! Counter
   CHARACTER(LEN=24), PARAMETER :: routinename = "PLOT_DIRECTION1D_CRP3D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT_DIRECTION1D_CRP3D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   npairpots=size(this%all_pairpots)
   ALLOCATE(v(npairpots+1))
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
   CALL this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
   WRITE(11,*) s, v(:)
   ! cycle for inpoints
   DO i=1, inpoints
      r(1)=xmin+(DFLOAT(i)*delta)*DCOS(alpha)
      r(2)=ymin+(DFLOAT(i)*delta)*DSIN(alpha)
      s = DSQRT(r(1)**2.D0+r(2)**2.D0)
      CALL this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
      WRITE(11,*) s, v(:)
   END DO
   ! Final value
   r(1) = xmax
   r(2) = ymax
   s = DSQRT(r(1)**2.D0+r(2)**2.D0)
   CALL this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
   WRITE(11,*) s, v(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT_DIRECTION1D_CORRECTION_CRP3D
!#######################################################################
! SUBROUTINE: PLOT_Z_CRP3D #############################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut along z direction
!
!> @param[in] this - CRP3D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!!                    cut. It should be given in degrees.
!> @param[in] xyz - Initial point
!> @param[in] L - Length of the graphic
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 09/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT_Z_CRP3D(this,npoints,xyz,L,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP3D),INTENT(IN):: this
   INTEGER,INTENT(IN):: npoints
   CHARACTER(LEN=*),INTENT(IN):: filename
   REAL(KIND=8),DIMENSION(3),INTENT(IN):: xyz
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8):: delta,L
   REAL(KIND=8):: zmax, zmin 
   REAL(KIND=8),DIMENSION(3):: r, dvdu
   REAL(KIND=8):: v
   INTEGER:: i ! Counter
   CHARACTER(LEN=*),PARAMETER:: routinename = "PLOT_DIRECTION1D_CRP3D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT_Z_CRP3D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   r(1)=xyz(1)
   r(2)=xyz(2)
   zmin=xyz(3)
   zmax=xyz(3)+L
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=L/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(3) = zmin
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)
   ! cycle for inpoints
   DO i=1, inpoints
      r(3)=zmin+(DFLOAT(i)*delta)
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(3),v,dvdu(:)
   END DO
   ! Final value
   r(3) = zmax
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT_Z_CRP3D
!#######################################################################
! SUBROUTINE: PLOT_DIRECTION1D_SMOOTH_CRP3D ############################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. To define 
!! the direction, the angle alpha is given. There's not any correction to
!! PES values.
!
!> @param[in] this - CRP3D PES used
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
SUBROUTINE PLOT_DIRECTION1D_SMOOTH_CRP3D(this,filename,npoints,angle,z,L)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP3D),INTENT(IN) :: this
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
   CALL this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) s, v , dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   ! cycle for inpoints
   DO i=1, inpoints
      r(1)=xmin+(DFLOAT(i)*delta)*DCOS(alpha)
      r(2)=ymin+(DFLOAT(i)*delta)*DSIN(alpha)
      s = DSQRT(r(1)**2.D0+r(2)**2.D0)
      CALL this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) s, v , dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   END DO
   ! Final value
   r(1) = xmax
   r(2) = ymax
   s = DSQRT(r(1)**2.D0+r(2)**2.D0)
   CALL this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) s, v, dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT_DIRECTION1D_SMOOTH_CRP3D
!###########################################################
!# SUBROUTINE: PLOT_INTERPOL_SYMMPOINT 
!###########################################################
!> @brief
!! Plots a Z scan along a symmpoint once it's been interpolated
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE PLOT_INTERPOL_SYMMPOINT(this,npoints,filename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Symmpoint),INTENT(IN):: this
   INTEGER(KIND=4),INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Run section
   CALL this%interz%PLOT(npoints,filename)
   RETURN
END SUBROUTINE PLOT_INTERPOL_SYMMPOINT
!###########################################################
!# SUBROUTINE: PLOT_PAIRPOTS_CRP3D 
!###########################################################
!> @brief
!! Plots all pairpot potentials
!
!> @warning
!! - It only works if there are less than 9 pairpotentials defined
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE PLOT_PAIRPOTS_CRP3D(this,npoints)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),INTENT(IN)::this
   INTEGER(KIND=4),INTENT(IN) :: npoints
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   INTEGER(KIND=4) :: npairpots
   CHARACTER(LEN=12) :: stringbase
   CHARACTER(LEN=13) :: filename
   ! Run section
   npairpots=size(this%all_pairpots)
   WRITE(stringbase,'(A12)') "-pairpot.dat"
   DO i = 1, npairpots
      WRITE(filename,'(I1,A12)') i,stringbase
      CALL this%all_pairpots(i)%PLOT(npoints,filename)
   END DO
   RETURN
END SUBROUTINE PLOT_PAIRPOTS_CRP3D
!###########################################################
!# SUBROUTINE: PLOT_SITIOS_CRP3D 
!###########################################################
!> @brief
!! Plots all sitio potentials
!
!> @warning
!! - It only works if there are less than 9 sitio potentials defined
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE PLOT_SITIOS_CRP3D(this,npoints)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),INTENT(IN)::this
   INTEGER(KIND=4),INTENT(IN) :: npoints
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   INTEGER(KIND=4) :: nsitios
   CHARACTER(LEN=10) :: stringbase
   CHARACTER(LEN=11) :: filename
   ! Run section
   nsitios=size(this%all_sites)
   WRITE(stringbase,'(A10)') "-sitio.dat"
   DO i = 1, nsitios
      WRITE(filename,'(I1,A10)') i,stringbase
      CALL this%all_sites(i)%PLOT(npoints,filename)
   END DO
   RETURN
END SUBROUTINE PLOT_SITIOS_CRP3D
!###########################################################
!# FUNCTION: is_allowed_CRP3D
!###########################################################
!> @brief
!! Determines if the potential can be calculated in this point
!-----------------------------------------------------------
LOGICAL FUNCTION is_allowed_CRP3D(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP3D),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8) :: xmin,xmax
   ! Run section
   xmin=this%all_sites(1)%z(1)
   xmax=this%all_sites(1)%z(this%all_sites(1)%n)
   SELECT CASE(size(x)/=3)
      CASE(.TRUE.)
         WRITE(0,*) "is_allowed_CRP3D ERR: array doesn't have 3 dimensions: 3"
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE( x(3)<xmin )
      CASE(.TRUE.)
         is_allowed_CRP3D=.FALSE.
      CASE(.FALSE.)
         is_allowed_CRP3D=.TRUE.
   END SELECT
   RETURN
END FUNCTION is_allowed_CRP3D
END MODULE CRP3D_MOD
