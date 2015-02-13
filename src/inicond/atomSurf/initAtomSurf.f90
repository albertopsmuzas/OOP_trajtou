!###############################################
! MODULE: INITATOM_MOD
!> @brief
!! This module provides routines and objects to create
!! initial conditions for an atom or list of atoms
!###############################################
MODULE INITATOMSURF_MOD
   use SYSTEM_MOD
   use INICOND_MOD
   use AOTUS_MODULE, only: flu_State, OPEN_CONFIG_FILE, CLOSE_CONFIG, AOT_GET_VAL
   use AOT_TABLE_MODULE, only: AOT_TABLE_OPEN, AOT_TABLE_CLOSE, AOT_TABLE_LENGTH, AOT_TABLE_GET_VAL
#ifdef DEBUG
   use DEBUG_MOD, only: VERBOSE_WRITE, DEBUG_WRITE
#endif
IMPLICIT NONE
!/////////////////////////////////////////////////////
! TYPE: AtomSurf
!> @brief
!! Atom subtype dynamics object
!----------------------------------------------------
TYPE,EXTENDS(DynObject)::  AtomSurf
   CONTAINS
      PROCEDURE,PUBLIC :: INITIALIZE => INITIALIZE_ATOMSURF
END TYPE AtomSurf
!/////////////////////////////////////////////////////
! TYPE: INITATOMSURF
!> @brief
!! Sets initial conditions for atoms
!----------------------------------------------------
TYPE,EXTENDS(IniCond):: InitAtomSurf
   LOGICAL:: control_vel, control_posX, control_posY, control_out, control_seed
   TYPE(Energy):: E_norm
   TYPE(Angle):: vz_angle, vpar_angle
   TYPE(Length):: init_z ! initial Z value
   type(Temperature):: temp
   type(Mass):: surfMass
   real(kind=8),dimension(3):: freqs
   REAL(KIND=8):: impact_x, impact_y
   CONTAINS
      PROCEDURE,PUBLIC:: INITIALIZE => INITIALIZE_INITATOMSURF
      PROCEDURE,PUBLIC:: GENERATE_TRAJS => GENERATE_TRAJS_INITATOMSURF
END TYPE InitAtomSurf
!////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: INITIALIZE_ATOMSURF #############################
!###########################################################
!> @brief
!! Initializes type atomSurf
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Feb/2015
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INITIALIZE_ATOMSURF(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(AtomSurf),INTENT(OUT):: this
   ! Local variables
   INTEGER(KIND=4),PARAMETER :: dimens=6
   ! Run section
   ALLOCATE(this%init_r(dimens))
   ALLOCATE(this%init_p(dimens))
   ALLOCATE(this%p(dimens))
   ALLOCATE(this%r(dimens))
   ALLOCATE(this%turning_point(dimens))
   this%stat="Dummy"
   RETURN
END SUBROUTINE INITIALIZE_ATOMSURF
!##################################################################################
!# SUBROUTINE: INITIALIZE_INITATOMSURF ################################################
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
!> @date 10/Feb/2015
!> @version 1.0
!------------------------------------------------------------------------------
SUBROUTINE INITIALIZE_INITATOMSURF(this,filename)
   IMPLICIT NONE
   ! I/O variables
   CLASS(InitAtomSurf),INTENT(OUT):: this
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: filename
   ! Local variables
   INTEGER:: i ! counters
   INTEGER:: size_seed, clock
   CHARACTER(LEN=*), PARAMETER :: routinename = "INITIALIZE_INITATOMSURF: "
   ! Lua variables
   TYPE(flu_State):: conf
   INTEGER(KIND=4):: inicond_table,trajlist_table,magnitude_table,control_table,out_table,osciSurf_table
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
      CASE('AtomSurfs')
         ! do nothing
      CASE DEFAULT
         WRITE(0,*) "INITIALIZE_IINITATOMSURF ERR: wrong kind of initial conditions"
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
         WRITE(0,*) "INITIALIZE_INITATOMSURF ERR: wrong incidence angle: ",this%vz_angle%getvalue()
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
         WRITE(0,*) "INITIALIZE_INITATOMSURF ERR: parallele angle lower than 0.0 deg."
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
      WRITE(0,*) "INITIALIZE_INITATOMSURF ERR: X impact parameter outside range 0-1"
      CALL EXIT(1)
   ELSE IF(((this%impact_y > 1.D0).OR.(this%impact_y < 0.D0)).AND.(.NOT.this%control_posy))THEN
      WRITE(0,*) "INITIALIZE_INITATOMSURF ERR: Y impact parameter outside range 0-1"
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
   ! get Surface oscillator parameters
   call aot_table_open(L=conf,parent=inicond_table,thandle=osciSurf_table,key='surfaceOscillator')
   call aot_get_val(L=conf,ErrCode=ierr,thandle=osciSurf_table,key='wx',val=this%freqs(1))
   call aot_get_val(L=conf,ErrCode=ierr,thandle=osciSurf_table,key='wy',val=this%freqs(2))
   call aot_get_val(L=conf,ErrCode=ierr,thandle=osciSurf_table,key='wz',val=this%freqs(3))
   select case(product(this%freqs(:))==0.d0)
      case(.true.)
         write(0,*) 'ERR: '//routinename//'One or more of the initial surface frequencies are zero'
         call exit(1)
      case(.false.)
         ! do nothing
   end select
   call aot_table_open(L=conf,parent=osciSurf_table,thandle=magnitude_table,key='temperature')
   call aot_get_val(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=1,val=auxReal)
   call aot_get_val(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=2,val=auxString)
   call this%temp%READ(auxReal,trim(auxString))
   call this%temp%to_std()
   call aot_table_close(L=conf,thandle=magnitude_table)
   call aot_table_open(L=conf,parent=osciSurf_table,thandle=magnitude_table,key='mass')
   call aot_get_val(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=1,val=auxReal)
   call aot_get_val(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=2,val=auxString)
   call this%surfMass%READ(auxReal,trim(auxString))
   call this%surfMass%to_std()
   call aot_table_close(L=conf,thandle=magnitude_table)
   call aot_table_close(L=conf,thandle=osciSurf_table)
#ifdef DEBUG
   call verbose_write(routinename,'Surface mass (au): ',this%surfMass%getValue())
   call verbose_write(routinename,'Surface frequencies (au): ',this%freqs(:))
   call verbose_write(routinename,'Surface temperature (kelvin): ',this%temp%getValue())
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
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=inicond_table,key='seedRead',val=this%control_seed)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Output files?: ",this%control_out)
   CALL VERBOSE_WRITE(routinename,"Output file name: ",this%output_file)
   CALL VERBOSE_WRITE(routinename,"Seed read from file?: ",this%control_seed)
#endif
   SELECT CASE(this%control_seed)
      CASE(.TRUE.)
         CALL RANDOM_SEED(SIZE=size_seed)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Default size for seed array: ",size_seed)
#endif
         ALLOCATE(this%seed(1:size_seed))
         OPEN(12,FILE="INseed.inp",STATUS="old")
         READ(12,*) this%seed
         CLOSE(12)
         CALL RANDOM_SEED(PUT=this%seed)
      CASE(.FALSE.)
         CALL RANDOM_SEED(SIZE=size_seed)
         ALLOCATE(this%seed(1:size_seed))
         CALL SYSTEM_CLOCK(COUNT=clock)
         this%seed = clock+ 37*(/ (i - 1, i = 1, size_seed) /)
         CALL RANDOM_SEED(PUT=this%seed)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Seed generated from CPU time: ",clock)
#endif
         OPEN(12,FILE="INseed.inp",STATUS="replace")
         WRITE(12,*) this%seed
         CLOSE(12)
   END SELECT
#ifdef DEBUG
   DO i=1,size_seed
      CALL VERBOSE_WRITE(routinename,this%seed(i))
   END DO
#endif
   RETURN
END SUBROUTINE INITIALIZE_INITATOMSURF
!####################################################################
!# SUBROUTINE: GENERATE_TRAJS_INITATOMSURF ##############################
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
!> @date Feb/2015
!> @version 1.0
!--------------------------------------------------------------------
SUBROUTINE GENERATE_TRAJS_INITATOMSURF(this,thispes)
   IMPLICIT NONE
   ! I/O variables
   CLASS(InitAtomSurf),INTENT(INOUT):: this
   CLASS(PES),INTENT(IN):: thispes
   ! IMPORTANT: unit used to write
   INTEGER(KIND=4),PARAMETER :: wunit=923
   ! Local variables
   INTEGER :: i ! counters
   CHARACTER(LEN=*),PARAMETER:: routinename = "GENERATE_TRAJS_ATOMS: "
   REAL(KIND=8):: delta,alpha,Enorm,masa,v
   REAL(KIND=8),DIMENSION(3) :: dummy
   REAL(KIND=8),DIMENSION(2) :: random_kernel
   real(kind=8):: beta
   ! YIPPIEE KI YAY !! -------------------
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"New set of trajectories")
   CALL VERBOSE_WRITE(routinename,"Allocating trajs: ", this%ntraj)
#endif
   beta=1.d0/(boltzmann*this%temp%getvalue())
   ALLOCATE(AtomSurf::this%trajs(this%ntraj))
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
      ! Initial surf parameters
      this%trajs(i)%r(4) = normalDistRandom()/(this%freqs(1)*dsqrt(beta*this%surfMass%getvalue()))
      this%trajs(i)%r(5) = normalDistRandom()/(this%freqs(2)*dsqrt(beta*this%surfMass%getvalue()))
      this%trajs(i)%r(6) = normalDistRandom()/(this%freqs(3)*dsqrt(beta*this%surfMass%getvalue()))
      this%trajs(i)%p(4) = normalDistRandom()/dsqrt(beta/this%surfMass%getvalue())
      this%trajs(i)%p(5) = normalDistRandom()/dsqrt(beta/this%surfMass%getvalue())
      this%trajs(i)%p(6) = normalDistRandom()/dsqrt(beta/this%surfMass%getvalue())
      ! Change to cartesian coordinates (impact parameters are in surface coordinates)
      this%trajs(i)%r(1:2) = system_surface%surf2cart(this%trajs(i)%r(1:2))
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
      CALL thispes%GET_V_AND_DERIVS(this%trajs(i)%r(1:3)-this%trajs(i)%r(4:6),v,dummy)
      dummy(:)=dot_product(this%trajs(i)%r(4:6),this%freqs(:))
      this%trajs(i)%E = 0.5d0*this%surfMass%getValue()*norm2(dummy(:))**2.d0+ &                 ! Surface's kinetic term +
                        (0.5d0/this%surfMass%getvalue())*norm2(this%trajs(i)%p(4:6))**2.d0+ &   ! Surface's harmonic potential +
                        Enorm/(DSIN(alpha)**2.D0)+v                                             ! Projectile's kinetic term + pot
      ! Setting initial values
      this%trajs(i)%init_r(:) = this%trajs(i)%r(:)
      this%trajs(i)%init_p(:) = this%trajs(i)%p(:)
   END DO
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Initial trajectory: ",this%nstart)
   CALL VERBOSE_WRITE(routinename,"Final trajectory: ",this%ntraj)
#endif
   ! Print if the option was given
   SELECT CASE(this%control_out)
      CASE(.TRUE.)
         OPEN(wunit,FILE=this%output_file,STATUS="replace")
         WRITE(wunit,*) "# FILE CREATED BY : GENERATE_TRAJS_ATOMSSURF ============================================================="
         WRITE(wunit,*) "# Format:   Id/E(au)/X,Y,Z,Xs,Ys,Zs(au)/Px,Py,Pz,Pxs,Pys,Pzs(au)"
         WRITE(wunit,*) "# Projectile's normal energy (a.u.) / (eV) : ", this%E_norm%getvalue(), " / ", this%E_norm%getvalue()*au2ev
         WRITE(wunit,*) "# Projectile's mass (a.u.) / proton mass : ", masa," / ", masa/pmass2au
         write(wunit,*) "# Surface's mass (au) / proton mass: ",this%surfMass%getValue()," / ",this%surfMass%getValue()/pmass2au
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
            WRITE(wunit,'(1X,I10,F15.7,2(6F20.7))') i,this%trajs(i)%E,this%trajs(i)%init_r(:),this%trajs(i)%init_p(:)
         END DO
         CLOSE(wunit)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Outputfile generated: ",this%output_file)
#endif
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE GENERATE_TRAJS_INITATOMSURF
END MODULE INITATOMSURF_MOD
