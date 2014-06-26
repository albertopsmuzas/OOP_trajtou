!###############################################
! MODULE: INITATOM_MOD
!> @brief
!! This module provides routines and onjects to create
!! initial conditions for an atom or list of atoms
!###############################################
MODULE INITATOM_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   USE INICOND_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////
! TYPE: Atoms
!> @brief
!! Atom subtype dynamics object
!----------------------------------------------------
TYPE,EXTENDS(Dynobject) ::  Atom
   CONTAINS
      PROCEDURE,PUBLIC :: INITIALIZE => INITIALIZE_ATOM
END TYPE Atom
!/////////////////////////////////////////////////////
! TYPE: INITATOM
!> @brief
!! Sets initial conditions for atoms
!----------------------------------------------------
TYPE,EXTENDS(Inicond) :: Initatom
   LOGICAL :: control_vel, control_posX, control_posY, control_out, control_seed
   REAL(KIND=8) :: impact_x, impact_y
   TYPE(Energy) :: E_norm
   TYPE(Angle) :: vz_angle, vpar_angle
   TYPE(Length) :: init_z ! initial Z value
   CONTAINS
      PROCEDURE,PUBLIC :: INITIALIZE => INITIALIZE_INITATOM
      PROCEDURE,PUBLIC :: GENERATE_TRAJS => GENERATE_TRAJS_INITATOM
END TYPE Initatom
!////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: INITIALIZE_ATOM #############################
!###########################################################
!> @brief
!! Initializes type atom
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INITIALIZE_ATOM(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Atom),INTENT(OUT):: this
   ! Local variables
   INTEGER(KIND=4),PARAMETER :: dimens=3
   ! Run section
   ALLOCATE(this%init_r(dimens))
   ALLOCATE(this%init_p(dimens))
   ALLOCATE(this%p(dimens))
   ALLOCATE(this%r(dimens))
   ALLOCATE(this%turning_point(dimens))
   this%stat="Dummy"
   RETURN
   RETURN
END SUBROUTINE INITIALIZE_ATOM
!##################################################################################
!# SUBROUTINE: INITIALIZE_INITATOM ################################################
!##################################################################################
!> @brief
!! Reads from an input file enough data to generate initial conditions
!! for a batch of atoms. Adapted to atom/surface systems
!
!> @param[out] this - Initial conditions for a batch of atoms
!> @param[in] filename - Name of the input file
!
!> @warning
!! - Input file syntax:
!!    -# line 1: character(len=30); human friendly alias 
!!    -# line 2: character(len=10); kind of input. "Atoms" is the only available label
!!    -# line 3: integer(kind=4); initial trajectory
!!    -# line 4: integer(kind=4); final trajectory
!!    -# line 5: real(kind=8),character(len=10); mass, units
!!    -# line 6: real(kind=8),character(len=10); initial perpendicular energy, units
!!    -# line 7: real(kind=8),character(len=10); angle of velocity respect to surface plane, units
!!    -# line 8: real(kind=8),character(len=10); angle of velocity respect to surface vector S1, units
!!    -# line 9: real(kind=8),character(len=10); initial Z, units
!!    -# line 10: logical,logical; random seed generation in X and Y?
!!    -# line 11: real(kind=8); initial X parameter in surface coordinates 1<= x >= 0
!!    -# line 12: real(kind=8); initial Y parameter in surface coordinates 1<= x >= 0
!!    -# line 13: logical; Create output file?
!!    -# line 14: character(len=30); name of the output file
!!    -# line 15: logical; read seed from INseed.inp?
!
!> @author A.S. Muzas alberto.muzas@uam.es
!> @date 10/Feb/2014
!> @version 1.0
!------------------------------------------------------------------------------
SUBROUTINE INITIALIZE_INITATOM(this,filename)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Initatom),INTENT(OUT) :: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! IMPORTANT: unit used to read info
   INTEGER(KIND=4),PARAMETER :: runit=134
   ! Local variables
   REAL(KIND=8) :: aux
   CHARACTER(LEN=10) :: units
   INTEGER :: i ! counters
   INTEGER :: size_seed, clock
   CHARACTER(LEN=23), PARAMETER :: routinename = "DEFINE_INICOND_SCHEME: "
   ! YIPEE KI YAY !! -------
   this%input_file = filename
   OPEN(runit,FILE=filename,STATUS="old")
   READ(runit,*) ! dumy line
   READ(runit,*) this%alias
   READ(runit,*) this%kind
   READ(runit,*) this%nstart
   READ(runit,*) this%ntraj
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Alias: ",this%alias)
   CALL VERBOSE_WRITE(routinename,"Kind: ",this%kind)
   CALL VERBOSE_WRITE(routinename,"Initial traj: ",this%nstart)
   CALL VERBOSE_WRITE(routinename,"Final traj: ",this%ntraj)
#endif
   IF (this%kind.EQ."Atoms") THEN
      READ(runit,*) aux,units
      CALL this%E_norm%READ(aux,units)
      CALL this%E_norm%TO_STD()

      READ(runit,*) aux,units
      CALL this%vz_angle%READ(aux,units)
      CALL this%vz_angle%TO_STD()
      IF (this%vz_angle%getvalue() > 90.D0*PI/180.D0) THEN
         WRITE(0,*) "DEFINE_INICOND_SCHEME ERR: incidence angle greater than 90.0 deg (1.57079632679D0 rad)"
         WRITE(0,*) "Incidence angle (rad) : ",this%vz_angle%getvalue()
         CALL EXIT(1)
      ELSE IF (this%vz_angle%getvalue() < 0.D0) THEN
         WRITE(0,*) "DEFINE_INICOND_SCHEME ERR: incidence angle lower than 0.0 deg."
         WRITE(0,*) "Incidence angle (rad) : ",this%vz_angle%getvalue()
         CALL EXIT(1)
      END IF

      READ(runit,*) aux,units
      CALL this%vpar_angle%READ(aux,units)
      CALL this%vpar_angle%TO_STD()
      IF (this%vpar_angle%getvalue() > 90.D0*PI/180.D0) THEN
         WRITE(0,*) "DEFINE_INICOND_SCHEME ERR: parallele angle greater than 90.0 deg (1.57079632679D0 rad)"
         WRITE(0,*) "Incidence angle (rad) : ",this%vpar_angle%getvalue()
         CALL EXIT(1)
      ELSE IF (this%vz_angle%getvalue() < 0.D0) THEN
         WRITE(0,*) "DEFINE_INICOND_SCHEME ERR: parallele angle lower than 0.0 deg."
         WRITE(0,*) "Incidence angle (rad) : ",this%vpar_angle%getvalue()
         CALL EXIT(1)
      END IF
      
      READ(runit,*) aux,units
      CALL this%init_z%READ(aux,units)
      CALL this%init_z%TO_STD()
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Initial normal energy in au:",this%E_norm%getvalue())
      CALL VERBOSE_WRITE(routinename,"Incidence angle respect to surface plane in radians: ",this%vz_angle%getvalue())
      CALL VERBOSE_WRITE(routinename,"Angle between trajectory and S1 vector in radians",this%vpar_angle%getvalue())
      CALL VERBOSE_WRITE(routinename,"Initial Z in au: ",this%init_z%getvalue())
#endif
      READ(runit,*) this%control_posX, this%control_posY 
      READ(runit,*) this%impact_x
      READ(runit,*) this%impact_y
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Random X impact parameter?: ",this%control_posX)
      CALL VERBOSE_WRITE(routinename,"Random Y impact parameter?: ",this%control_posY)
      CALL VERBOSE_WRITE(routinename,"Impact param X: ",this%impact_x)
      CALL VERBOSE_WRITE(routinename,"Impact param Y: ",this%impact_x)
#endif
      IF(((this%impact_x > 1.D0).OR.(this%impact_x < 0.D0)).AND.(.NOT.this%control_posX)) THEN
         WRITE(0,*) "DEFINE_INICOND_SCHEME ERR: X impact parameter outside range 0-1"
         CALL EXIT(1)
      ELSE IF(((this%impact_y > 1.D0).OR.(this%impact_y < 0.D0)).AND.(.NOT.this%control_posy))THEN
         WRITE(0,*) "DEFINE_INICOND_SCHEME ERR: Y impact parameter outside range 0-1"
         CALL EXIT(1)
      END IF
      READ(runit,*) this%control_out
      READ(runit,*) this%output_file
      READ(runit,*) this%control_seed
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Output files?: ",this%control_out)
      CALL VERBOSE_WRITE(routinename,"Output file name: ",this%output_file)
      CALL VERBOSE_WRITE(routinename,"Seed read from file?: ",this%control_seed)
#endif
      IF (this%control_seed.EQV..TRUE.) THEN
         CALL RANDOM_SEED(SIZE=size_seed)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Default size for seed array: ",size_seed)
#endif
         ALLOCATE(this%seed(1:size_seed))
         OPEN(12,FILE="INseed.inp",STATUS="old")
         READ(12,*) this%seed
         CLOSE(12)
      ELSE
         CALL RANDOM_SEED(SIZE=size_seed)
         ALLOCATE(this%seed(1:size_seed))
         CALL SYSTEM_CLOCK(COUNT=clock)
         this%seed = clock+ 37*(/ (i - 1, i = 1, size_seed) /)
         CALL RANDOM_SEED(PUT=this%seed)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Seed generated from CPU time: ")
         CALL VERBOSE_WRITE(routinename,"CPU time: ", clock)
#endif
         IF (this%control_out.EQV..TRUE.) THEN
            OPEN(12,FILE="INseed.inp",STATUS="replace")
            WRITE(12,*) this%seed
            CLOSE(12)
         END IF
      END IF
#ifdef DEBUG
      DO i=1,size_seed
         CALL VERBOSE_WRITE(routinename,this%seed(i))
      END DO
#endif
   ELSE
      WRITE(0,*) "DEFINE_INICOND_ATOM ERR: Wrong kind of initial conditions"
      WRITE(0,*) "Available kinds: Atoms"
      WRITE(0,*) "You wrote: ", this%kind
      CALL EXIT(1)
   END IF
   CLOSE(runit)
   RETURN
END SUBROUTINE INITIALIZE_INITATOM
!####################################################################
!# SUBROUTINE: GENERATE_TRAJS_INITATOM ##############################
!####################################################################
!> @brief
!! Creates trajectories of atoms. All data stored as atomic units.
!
!> @param[in] this - Initial conditions for atoms to be used
!> @param[in] surf - Surface used
!
!> @warning
!! - Auxiliar cartesian coordinates.
!! - For atom/surface systems
!
!> @author A.S. Muzas alberto.muzas@uam.es
!> @date Feb/2014, Jun/2014
!> @version 1.0, 1.1
!--------------------------------------------------------------------
SUBROUTINE GENERATE_TRAJS_INITATOM(this,thispes)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Initatom),INTENT(INOUT):: this
   CLASS(PES),INTENT(IN):: thispes
   ! IMPORTANT: unit used to write
   INTEGER(KIND=4),PARAMETER :: wunit=923
   ! Local variables
   INTEGER :: i,j ! counters
   CHARACTER(LEN=22),PARAMETER:: routinename = "GENERATE_TRAJS_ATOMS: "
   REAL(KIND=8),DIMENSION(2):: proj_iws_r
   REAL(KIND=8):: delta,alpha,Enorm,masa
   ! YIPPIEE KI YAY !! -------------------
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"New set of trajectories")
   CALL VERBOSE_WRITE(routinename,"Allocating trajs: ", this%ntraj)
#endif
   ALLOCATE(Atom::this%trajs(this%ntraj))
   IF((this%control_posX).AND.(.NOT.this%control_posY)) THEN
      DO i=1, this%ntraj
         CALL this%trajs(i)%INITIALIZE()
         this%trajs(i)%r(3) = this%init_z%getvalue()
         CALL RANDOM_NUMBER(this%trajs(i)%r(1))
         this%trajs(i)%r(2) = this%impact_y
      END DO
   ELSE IF ((.NOT.this%control_posX).AND.(this%control_posY)) THEN
      DO i=1, this%ntraj
         CALL this%trajs(i)%INITIALIZE()
         this%trajs(i)%r(3) = this%init_z%getvalue()
         this%trajs(i)%r(1) = this%impact_x
         CALL RANDOM_NUMBER(this%trajs(i)%r(2))
      END DO
   ELSE IF ((this%control_posX).AND.(this%control_posY)) THEN
      DO i=1, this%ntraj
         CALL this%trajs(i)%INITIALIZE()
         this%trajs(i)%r(3) = this%init_z%getvalue()
         CALL RANDOM_NUMBER(this%trajs(i)%r(1:2))
      END DO
   ELSE
      DO i=1, this%ntraj
         CALL this%trajs(i)%INITIALIZE()
         this%trajs(i)%r(3) = this%init_z%getvalue()
         this%trajs(i)%r(1) = this%impact_x
         this%trajs(i)%r(2) = this%impact_y
      END DO
   END IF
   DO i=1,this%ntraj
      ! Change to cartesian coordinates (impact parameters are in surface coordinates)
      this%trajs(i)%r(1:2) = thispes%surf%surf2cart(this%trajs(i)%r(1:2))
      ! projectin into IWS cell leads to errors in the dynamics (wrong sampling)
      IF (thispes%surf%units/="au") THEN
         WRITE(0,*) "GENERATE_TRAJS_ATOMS ERR: Surface should be in atomic units!"
         CALL EXIT(1)
      END IF   
      delta = this%vpar_angle%getvalue()
      alpha = this%vz_angle%getvalue()
      masa = thispes%atomdat(1)%getmass()
      Enorm = this%E_norm%getvalue()
      this%trajs(i)%p(3) = -DSQRT(2.D0*masa*Enorm)  ! Z momentum (m*v_z), negative (pointing uppon the surface)
      this%trajs(i)%p(1) = DCOS(delta)*DSQRT(2.D0*masa*Enorm/(DTAN(alpha)**2.D0))
      this%trajs(i)%p(2) = DSIN(delta)*DSQRT(2.D0*masa*Enorm/(DTAN(alpha)**2.D0))
      this%trajs(i)%E = Enorm/(DSIN(alpha)**2.D0)
      ! Setting initial values
      this%trajs(i)%init_r(1) = this%trajs(i)%r(1)
      this%trajs(i)%init_r(2) = this%trajs(i)%r(2)
      this%trajs(i)%init_r(3) = this%trajs(i)%r(3)
      this%trajs(i)%init_p(1) = this%trajs(i)%p(1)
      this%trajs(i)%init_p(2) = this%trajs(i)%p(2)
      this%trajs(i)%init_p(3) = this%trajs(i)%p(3)
   END DO
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Initial trajectory: ",this%nstart)
   CALL VERBOSE_WRITE(routinename,"Final trajectory: ",this%ntraj)
#endif
   ! Print if the option was given
   SELECT CASE(this%control_out)
      CASE(.TRUE.)
         OPEN(wunit,FILE=this%output_file,STATUS="replace")
         WRITE(wunit,*) "# FILE CREATED BY : GENERATE_TRAJS_ATOMS ================================================================="
         WRITE(wunit,*) "# Format:   traj_num      X,Y,Z (a.u.)      Px,Py,Pz(a.u.)   X,Y(Projected IWS, a.u.)"
         WRITE(wunit,*) "# Initial total Energy (a.u.) / (eV) : ", this%trajs(1)%E, " /  ", this%trajs(1)%E*au2ev
         WRITE(wunit,*) "# Perpendicular Energy (a.u.) / (eV) : ", this%E_norm%getvalue(), " / ", this%E_norm%getvalue()*au2ev
         WRITE(wunit,*) "# MASS (a.u.) / proton_mass : ", masa," / ", masa/pmass2au
         WRITE(wunit,*) "# Incidence angle (deg): ", this%vz_angle%getvalue()*180.D0/PI
         WRITE(wunit,*) "# Parallel velocity direction (deg): ", this%vpar_angle%getvalue()*180.D0/PI
         IF ((this%control_posX).AND.(.NOT.this%control_posY)) THEN
            WRITE(wunit,*) "# Random X impact parameter"
            WRITE(wunit,*) "# Seed used: ", this%seed
         ELSE IF ((.NOT.this%control_posX).AND.(this%control_posY)) THEN
            WRITE(wunit,*) "# Random Y impact parameter"
            WRITE(wunit,*) "# Seed used: ", this%seed
         ELSE IF ((this%control_posX).AND.(this%control_posY)) THEN
            WRITE(wunit,*) "# Random X and Y impact parameters"
            WRITE(wunit,*) "# Seed used: ", this%seed
         ELSE
            WRITE(wunit,*) "# X, Y values are not random numbers"
         END IF
         WRITE(wunit,*) "# ======================================================================================================="
         DO i=this%nstart,this%ntraj
            FORALL(j=1:2) proj_iws_r(j) = this%trajs(i)%init_r(j)
            proj_iws_r = thispes%surf%project_iwscell(proj_iws_r)
            WRITE(wunit,'(1X,I10,3(3F20.7))') i,this%trajs(i)%init_r,this%trajs(i)%init_p,proj_iws_r
         END DO
         CLOSE(wunit)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Outputfile generated: ",this%output_file)
#endif
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE GENERATE_TRAJS_INITATOM
END MODULE INITATOM_MOD
