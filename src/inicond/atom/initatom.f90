!###############################################
! MODULE: INITATOM_MOD
!> @brief
!! This module provides routines and objects to create
!! initial conditions for an atom or list of atoms
!###############################################
MODULE INITATOM_MOD
   USE SYSTEM_MOD
   USE INICOND_MOD
   USE AOTUS_MODULE, ONLY: flu_State, OPEN_CONFIG_FILE, CLOSE_CONFIG, AOT_GET_VAL
   USE AOT_TABLE_MODULE, ONLY: AOT_TABLE_OPEN, AOT_TABLE_CLOSE, AOT_TABLE_LENGTH, AOT_TABLE_GET_VAL
#ifdef DEBUG
   USE DEBUG_MOD
#endif
IMPLICIT NONE
!/////////////////////////////////////////////////////
! TYPE: Atoms
!> @brief
!! Atom subtype dynamics object
!----------------------------------------------------
TYPE,EXTENDS(Dynobject)::  Atom
   CONTAINS
      PROCEDURE,PUBLIC :: INITIALIZE => INITIALIZE_ATOM
END TYPE Atom
!/////////////////////////////////////////////////////
! TYPE: INITATOM
!> @brief
!! Sets initial conditions for atoms
!----------------------------------------------------
TYPE,EXTENDS(Inicond) :: Initatom
   LOGICAL :: control_vel, control_posX, control_posY, control_out
   TYPE(Energy) :: E_norm
   TYPE(Angle) :: vz_angle, vpar_angle
   REAL(KIND=8) :: impact_x, impact_y
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
!! - Input file syntax: lua
!
!> @author A.S. Muzas alberto.muzas@uam.es
!> @date 10/Feb/2014
!> @version 1.0
!------------------------------------------------------------------------------
SUBROUTINE INITIALIZE_INITATOM(this,filename)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Initatom),INTENT(OUT):: this
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: filename
   ! Local variables
   INTEGER:: i ! counters
   INTEGER:: size_seed, clock
   CHARACTER(LEN=*), PARAMETER :: routinename = "INITIALIZE_INITATOM: "
   ! Lua variables
   TYPE(flu_State):: conf
   INTEGER(KIND=4):: inicond_table,trajlist_table,magnitude_table,control_table,out_table
   INTEGER(KIND=4):: auxtable
   INTEGER(KIND=4):: ierr
   ! Auxiliary (dummy) variables
   CHARACTER(LEN=1024):: auxstring
   INTEGER(KIND=4):: auxint
   REAL(KIND=8):: auxreal
   ! YIPEE KI YAY !! -------
   SELECT CASE(allocated(system_inputfile) .or. .not.present(filename))
      CASE(.TRUE.)
         auxstring=system_inputfile
      CASE(.FALSE.)
         auxstring=filename
   END SELECT
   ! Open lua conf file
   ALLOCATE(this%input_file,source=trim(auxstring))
   CALL OPEN_CONFIG_FILE(L=conf,ErrCode=ierr,filename=this%input_file)
   ! Open initial conditions table
   CALL AOT_TABLE_OPEN(L=conf,thandle=inicond_table,key='initialConditions')
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=inicond_table,key='kind',val=auxstring)
   this%kind=trim(auxstring)
   SELECT CASE(this%kind)
      CASE('Atoms')
         ! do nothing
      CASE DEFAULT
         WRITE(0,*) "INITIALIZE_IINITATOM ERR: wrong kind of initial conditions"
         CALL EXIT(1)
   END SELECT
   ! get traj list to initialize
   CALL AOT_TABLE_OPEN(L=conf,parent=inicond_table,thandle=trajlist_table,key='trajList')
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=trajlist_table,key='from',val=auxint)
   this%nstart=auxint
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=trajlist_table,key='to',val=auxint)
   this%ntraj=auxint
   CALL AOT_TABLE_CLOSE(L=conf,thandle=trajlist_table)
   ! Some debugging messages
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Kind: ",this%kind)
   CALL VERBOSE_WRITE(routinename,"Initial traj: ",this%nstart)
   CALL VERBOSE_WRITE(routinename,"Final traj: ",this%ntraj)
#endif
   ! get enorm
   CALL AOT_TABLE_OPEN(L=conf,parent=inicond_table,thandle=magnitude_table,key='Enormal')
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=1,val=auxreal)
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=2,val=auxstring)
   CALL this%E_norm%READ(auxreal,trim(auxstring))
   CALL this%E_norm%TO_STD()
   CALL AOT_TABLE_CLOSE(L=conf,thandle=magnitude_table)
   ! get incidence angle
   CALL AOT_TABLE_OPEN(L=conf,parent=inicond_table,thandle=magnitude_table,key='incidenceAngle')
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=1,val=auxreal)
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=2,val=auxstring)
   CALL this%vz_angle%READ(auxreal,trim(auxstring))
   CALL this%vz_angle%TO_STD()
   CALL AOT_TABLE_CLOSE(L=conf,thandle=magnitude_table)
   SELECT CASE(this%vz_angle%getvalue() > 90.D0*PI/180.D0 .or. this%vz_angle%getvalue() < 0.D0)
      CASE(.true.)
         WRITE(0,*) "INITIALIZE_INITATOM ERR: wrong incidence angle: ",this%vz_angle%getvalue()
         CALL EXIT(1)
      CASE(.false.)
         ! do nothing
   END SELECT
   ! get direction angle
   CALL AOT_TABLE_OPEN(L=conf,parent=inicond_table,thandle=magnitude_table,key='directionAngle')
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=1,val=auxreal)
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=2,val=auxstring)
   CALL this%vpar_angle%READ(auxreal,trim(auxstring))
   CALL this%vpar_angle%TO_STD()
   CALL AOT_TABLE_CLOSE(L=conf,thandle=magnitude_table)
   SELECT CASE(this%vz_angle%getvalue() < 0.D0)
      CASE(.true.)
         WRITE(0,*) "INITIALIZE_INITATOM ERR: parallele angle lower than 0.0 deg."
         WRITE(0,*) "Incidence angle (rad) : ",this%vpar_angle%getvalue()
         CALL EXIT(1)
      CASE(.false.)
         ! do nothing
   END SELECT
   ! get initial Z
   CALL AOT_TABLE_OPEN(L=conf,parent=inicond_table,thandle=magnitude_table,key='initialZ')
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=1,val=auxreal)
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=2,val=auxstring)
   CALL this%init_z%READ(auxreal,trim(auxstring))
   CALL this%init_z%TO_STD()
   CALL AOT_TABLE_CLOSE(L=conf,thandle=magnitude_table)
   ! Some debugging info
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Initial normal energy in au:",this%E_norm%getvalue())
   CALL VERBOSE_WRITE(routinename,"Incidence angle respect to surface plane in radians: ",this%vz_angle%getvalue())
   CALL VERBOSE_WRITE(routinename,"Angle between trajectory and S1 vector in radians",this%vpar_angle%getvalue())
   CALL VERBOSE_WRITE(routinename,"Initial Z in au: ",this%init_z%getvalue())
#endif
   ! get control random initial XY position
   CALL AOT_TABLE_OPEN(L=conf,parent=inicond_table,thandle=control_table,key='randomXY')
   CALL AOT_TABLE_OPEN(L=conf,parent=control_table,thandle=auxtable,key='X')
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=auxtable,pos=1,val=this%control_posX)
   SELECT CASE(this%control_posX)
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=auxtable,pos=2,val=this%impact_x)
   END SELECT
   CALL AOT_TABLE_CLOSE(L=conf,thandle=auxtable)
   CALL AOT_TABLE_OPEN(L=conf,parent=control_table,thandle=auxtable,key='Y')
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=auxtable,pos=1,val=this%control_posY)
   SELECT CASE(this%control_posY)
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=auxtable,pos=2,val=this%impact_y)
   END SELECT
   CALL AOT_TABLE_CLOSE(L=conf,thandle=auxtable)
   IF(((this%impact_x > 1.D0).OR.(this%impact_x < 0.D0)).AND.(.NOT.this%control_posX)) THEN
      WRITE(0,*) "INITIALIZE_INITATOM ERR: X impact parameter outside range 0-1"
      CALL EXIT(1)
   ELSE IF(((this%impact_y > 1.D0).OR.(this%impact_y < 0.D0)).AND.(.NOT.this%control_posy))THEN
      WRITE(0,*) "INITIALIZE_INITATOM ERR: Y impact parameter outside range 0-1"
      CALL EXIT(1)
   END IF
   ! Some debugging messages
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Random X impact parameter?: ",this%control_posX)
   CALL VERBOSE_WRITE(routinename,"Random Y impact parameter?: ",this%control_posY)
   IF (this%control_posX.eqv..false. .and. this%control_posY.eqv..false.) THEN
      CALL VERBOSE_WRITE(routinename,"Impact param X: ",this%impact_x)
      CALL VERBOSE_WRITE(routinename,"Impact param Y: ",this%impact_x)
   ELSE IF (.not.this%control_posX) THEN
      CALL VERBOSE_WRITE(routinename,"Impact param X: ",this%impact_x)
   ELSE IF (.not.this%control_posY) THEN
      CALL VERBOSE_WRITE(routinename,"Impact param X: ",this%impact_y)
   END IF
#endif
   ! get output control
   CALL AOT_TABLE_OPEN(L=conf,parent=inicond_table,thandle=out_table,key='outputFile')
   auxint=aot_table_length(L=conf,thandle=out_table)
   SELECT CASE(auxint)
      CASE(0)
         this%control_out=.false.
      CASE DEFAULT
         this%control_out=.true.
         CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=out_table,pos=1,val=auxstring)
         this%output_file=trim(auxstring)
   END SELECT
   CALL AOT_TABLE_CLOSE(L=conf,thandle=out_table)
   ! get seed control
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Output files?: ",this%control_out)
   CALL VERBOSE_WRITE(routinename,"Output file name: ",this%output_file)
#endif
   call generate_seed()
   call random_seed(put=system_iSeed(:))
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'System seed: ',system_iSeed(:))
#endif
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
   INTEGER :: i ! counters
   CHARACTER(LEN=*),PARAMETER:: routinename = "GENERATE_TRAJS_ATOMS: "
   REAL(KIND=8):: delta,alpha,Enorm,masa,v
   REAL(KIND=8),DIMENSION(3) :: dummy
   REAL(KIND=8),DIMENSION(2) :: random_kernel
   ! YIPPIEE KI YAY !! -------------------
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"New set of trajectories")
   CALL VERBOSE_WRITE(routinename,"Allocating trajs: ", this%ntraj)
#endif
   ALLOCATE(Atom::this%trajs(this%ntraj))
   DO i=1,this%ntraj
      CALL RANDOM_NUMBER(random_kernel(:))
      IF((this%control_posX).AND.(.NOT.this%control_posY)) THEN
         CALL this%trajs(i)%INITIALIZE()
         this%trajs(i)%r(1) = random_kernel(1)
         this%trajs(i)%r(2) = this%impact_y
      ELSE IF ((.NOT.this%control_posX).AND.(this%control_posY)) THEN
         CALL this%trajs(i)%INITIALIZE()
         this%trajs(i)%r(1) = this%impact_x
         this%trajs(i)%r(2) = random_kernel(2)
      ELSE IF ((this%control_posX).AND.(this%control_posY)) THEN
         CALL this%trajs(i)%INITIALIZE()
         this%trajs(i)%r(1:2)=random_kernel(1:2)
      ELSE
         CALL this%trajs(i)%INITIALIZE()
         this%trajs(i)%r(1) = this%impact_x
         this%trajs(i)%r(2) = this%impact_y
      END IF
      this%trajs(i)%r(3) = this%init_z%getvalue()
      ! Change to cartesian coordinates (impact parameters are in surface coordinates)
      this%trajs(i)%r(1:2) = system_surface%surf2cart(this%trajs(i)%r(1:2))
      ! projectin into IWS cell leads to errors in the dynamics (wrong sampling)
      IF (system_surface%units/="au") THEN
         WRITE(0,*) "GENERATE_TRAJS_ATOMS ERR: Surface should be in atomic units!"
         CALL EXIT(1)
      END IF   
      delta = this%vpar_angle%getvalue()
      alpha = this%vz_angle%getvalue()
      masa = system_mass(1)
      Enorm = this%E_norm%getvalue()
      this%trajs(i)%p(3) = -DSQRT(2.D0*masa*Enorm)  ! Z momentum (m*v_z), negative (pointing uppon the surface)
      this%trajs(i)%p(1) = DCOS(delta)*DSQRT(2.D0*masa*Enorm/(DTAN(alpha)**2.D0))
      this%trajs(i)%p(2) = DSIN(delta)*DSQRT(2.D0*masa*Enorm/(DTAN(alpha)**2.D0))
      CALL thispes%GET_V_AND_DERIVS(this%trajs(i)%r,v,dummy)
      this%trajs(i)%E = Enorm/(DSIN(alpha)**2.D0)+v
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
         WRITE(wunit,*) "# Format:   traj_num      X,Y,Z (a.u.)      Px,Py,Pz(a.u.)"
         WRITE(wunit,*) "# Initial total Energy (a.u.) / (eV) : ", this%trajs(1)%E, " /  ", this%trajs(1)%E*au2ev
         WRITE(wunit,*) "# Perpendicular Energy (a.u.) / (eV) : ", this%E_norm%getvalue(), " / ", this%E_norm%getvalue()*au2ev
         WRITE(wunit,*) "# MASS (a.u.) / proton_mass : ", masa," / ", masa/pmass2au
         WRITE(wunit,*) "# Incidence angle (deg): ", this%vz_angle%getvalue()*180.D0/PI
         WRITE(wunit,*) "# Parallel velocity direction (deg): ", this%vpar_angle%getvalue()*180.D0/PI
         IF ((this%control_posX).AND.(.NOT.this%control_posY)) THEN
            WRITE(wunit,*) "# Random X impact parameter"
            WRITE(wunit,*) "# Seed used: ",system_iSeed(:)
         ELSE IF ((.NOT.this%control_posX).AND.(this%control_posY)) THEN
            WRITE(wunit,*) "# Random Y impact parameter"
            WRITE(wunit,*) "# Seed used: ",system_iSeed(:)
         ELSE IF ((this%control_posX).AND.(this%control_posY)) THEN
            WRITE(wunit,*) "# Random X and Y impact parameters"
            WRITE(wunit,*) "# Seed used: ",system_iSeed(:)
         ELSE
            WRITE(wunit,*) "# X, Y values are not random numbers"
         END IF
         WRITE(wunit,*) "# ======================================================================================================="
         DO i=this%nstart,this%ntraj
            WRITE(wunit,'(1X,I10,3(3F20.7))') i,this%trajs(i)%init_r,this%trajs(i)%init_p
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
