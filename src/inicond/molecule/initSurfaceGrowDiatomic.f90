!###################################################################
! MODULE: INITSURFACEGROW_MOD
!> @brief
!! This module provides routines and objects to read
!! initial conditions for a diatomic molecule or a list of them inside
!! surfaceGrow program
!###################################################################
MODULE INITSURFACEGROW_MOD
   use SYSTEM_MOD
   use INICOND_MOD
   use AOTUS_MODULE, only: flu_State, OPEN_CONFIG_FILE, CLOSE_CONFIG, AOT_GET_VAL,aoterr_WrongType
   use AOT_TABLE_MODULE, only: AOT_TABLE_OPEN, AOT_TABLE_CLOSE, AOT_TABLE_LENGTH, AOT_TABLE_GET_VAL
#ifdef DEBUG
   use DEBUG_MOD, only: VERBOSE_WRITE,DEBUG_WRITE
#endif
IMPLICIT NONE
!/////////////////////////////////////////////////////
! TYPE: DiatomicGrow
!> @brief
!! DiatomicGrow subtype dynamics object
!
!> @param   evirot - rovibrational energy
!> @param   init_qn - initial quantum number (v,J)
!> @param   final_qn - final quantum number (v,J)
!----------------------------------------------------
TYPE,EXTENDS(Dynobject) ::  DiatomicGrow
   CONTAINS
      PROCEDURE,PUBLIC :: INITIALIZE => INITIALIZE_DIATOMICGROW
END TYPE DiatomicGrow
!///////////////////////////////////////////////////////////////////
! TYPE: InitSurfaceGrowDiatomic
!> @brief
!! Sets initial conditions for a diatomic molecule. Wrapper for grow
!-------------------------------------------------------------------
TYPE,EXTENDS(Inicond):: InitSurfaceGrowDiatomic
   CHARACTER(LEN=:),ALLOCATABLE::extrapol
   REAL(KIND=8):: period
   REAL(KIND=8):: eps
   LOGICAL:: control_vel,control_posX,control_posY,control_out
   LOGICAL:: is_classic=.FALSE.
   LOGICAL:: fixed_theta
   REAL(KIND=8):: impact_x, impact_y
   INTEGER(KIND=4),DIMENSION(3):: init_qn
   TYPE(Csplines):: pr_t,r_t
   TYPE(Time):: delta_t
   TYPE(Energy):: E_norm, evirot
   TYPE(Angle):: vz_angle, vpar_angle
   TYPE(Length):: init_z ! initial Z value
   TYPE(Vacuumpot):: vibrpot
   TYPE(Vacuumpot):: rovibrpot
   CONTAINS
      ! Initialization block
      PROCEDURE,PUBLIC:: INITIALIZE => INITIALIZE_INITSURFACEGROW
      ! Tools block
      PROCEDURE,PUBLIC:: GENERATE_TRAJS => GENERATE_TRAJS_INITSURFACEGROW
      PROCEDURE,PUBLIC:: GENERATE_TRAJS_FROM_FILE => GENERATE_TRAJS_FROM_FILE_INITSURFACEGROW
      ! Private routines
      PROCEDURE,PRIVATE:: SET_PERIOD => SET_PERIOD_INITSURFACEGROW
      PROCEDURE,PRIVATE:: TIME_DERIVS => TIME_DERIVS_INITSURFACEGROW
      PROCEDURE,PRIVATE:: MMID => MMID_INITSURFACEGROW
      PROCEDURE,PRIVATE:: BSSTEP => BSSTEP_INITSURFACEGROW
      PROCEDURE,PRIVATE:: POLINOM_EXTRAPOL => PZEXTR_INITSURFACEGROW
      PROCEDURE,PRIVATE:: RATIONAL_EXTRAPOL => RZEXTR_INITSURFACEGROW
END TYPE InitSurfaceGrowDiatomic
CONTAINS
!#################################################################################
!# SUBROUTINE: GENERATE_TRAJS_FROM_FILE_INITSURFACEGROW ##########################
!#################################################################################
!> @brief
!! Initializes trajectories Initdiatomic from a Grow output file. (init.raw.dat)
!
!> @param[in] this - Initial conditions for molecules to be used
!--------------------------------------------------------------------
SUBROUTINE GENERATE_TRAJS_FROM_FILE_INITSURFACEGROW(this,fileName)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(InitSurfaceGrowDiatomic),INTENT(INOUT):: this
   CHARACTER(LEN=*),INTENT(IN) :: fileName
   ! IMPORTANT: unit used to write
   INTEGER(KIND=4),PARAMETER:: wunit=923
   integer(kind=4),parameter:: runit=922
   ! Local variables
   INTEGER :: i ! counters
   real(kind=8),dimension(12):: phaseSpaceVect
   REAL(KIND=8):: delta,alpha,Enorm,masa,mu,Eint,Epar
   real(kind=8),dimension(6):: p,r
   ! Parameters
   character(len=*),parameter:: formatFile='("# Format: id/Etot/X,Y,Z,R(a.u.)/THETA,PHI(rad)/Px,Py,Pz,pr,ptheta,pphi(a.u.)")'
   character(len=*),parameter:: routinename = "GENERATE_TRAJS_ATOMS: "
   ! YIPPIEE KI YAY !! -------------------
   delta = this%vpar_angle%getvalue()
   alpha = this%vz_angle%getvalue()
   masa = sum(system_mass(1:2))
   mu = product(system_mass(1:2))/masa
   Enorm = this%E_norm%getvalue()
   Epar = Enorm/(dtan(alpha)**2.d0)
   Eint = this%evirot%getvalue()
   ALLOCATE(DiatomicGrow::this%trajs(this%ntraj))
   open(unit=runit,file=fileName,status='old',action='read')
   do i=1,this%ntraj
      call this%trajs(i)%initialize()
      read(runit,*)
      read(runit,*)
      read(runit,*) 
      read(runit,*)
      read(runit,*) this%trajs(i)%init_r(1:3)
      read(runit,*) this%trajs(i)%init_r(4:6)
      read(runit,*)
      read(runit,*) this%trajs(i)%init_p(1:3) ! velocities!
      read(runit,*) this%trajs(i)%init_p(4:6) ! velocities!
      this%trajs(i)%init_p(1:3)=this%trajs(i)%init_p(1:3)*system_mass(1)/dsqrt(pmass2au)
      this%trajs(i)%init_p(4:6)=this%trajs(i)%init_p(4:6)*system_mass(2)/dsqrt(pmass2au)
      phaseSpaceVect(1:6)=this%trajs(i)%init_r(:)
      phaseSpaceVect(7:12)=this%trajs(i)%init_p(:)
      phaseSpaceVect(:)=from_atomic_to_molecular_phaseSpace(atomcoord=phaseSpaceVect)
      this%trajs(i)%init_r(:)=phaseSpaceVect(1:6)
      this%trajs(i)%init_p(:)=phaseSpaceVect(7:12)
      this%trajs(i)%init_p(1)=dcos(delta)*dsqrt(2.d0*masa*Epar)
      this%trajs(i)%init_p(2)=dsin(delta)*dsqrt(2.d0*masa*Epar)
      this%trajs(i)%init_p(3)=-dsqrt(2.d0*masa*Enorm)
      p(:)=this%trajs(i)%init_p(:)
      r(:)=this%trajs(i)%init_r(:)
      !write(*,*) 0.5d0*(p(4)**2.d0)/mu+0.5d0*(p(5)**2.d0+(p(6)/dsin(r(5)))**2.d0)/(mu*r(4)**2.d0)+this%vibrPot%getPot(r(4))
   enddo
   close(runit)
   ! Print if the option was given
   SELECT CASE(this%control_out)
      CASE(.TRUE.)
         OPEN(wunit,FILE=this%output_file,STATUS="replace")
         WRITE(wunit,'("# FILE CREATED BY : GENERATE_TRAJS_ATOMS =================================================================")')
         WRITE(wunit,formatFile) 
         WRITE(wunit,'("# Perpendicular Energy (a.u.) / (eV): ",F10.5," / ",F10.5)') Enorm,Enorm*au2ev
         WRITE(wunit,'("# Initial internal energy (a.u.) / (eV) : ",F10.5," / ",F10.5)') Eint,Eint*au2ev
         WRITE(wunit,'("# Initial center of mass energy (a.u.) / (eV) : ",F10.5," / ",F10.5)') Enorm/(DSIN(alpha)**2.D0),&
            (Enorm/(DSIN(alpha)**2.D0))*au2ev
         WRITE(wunit,'("# Initial Rovibrational state: ",3(I5,1X))') this%init_qn(:)
         WRITE(wunit,'("# MASS (a.u.) / proton_mass : ",F10.5," / ",F10.5)')  masa,masa/pmass2au
         WRITE(wunit,'("# Reduced MASS (a.u.) / proton_mass : ",F10.5," / ",F10.5)') mu,mu/pmass2au
         WRITE(wunit,'("# Incidence angle (deg): ",F10.5)')  alpha*180.D0/PI
         WRITE(wunit,'("# Parallel velocity direction (deg): ",F10.5)') delta*180.D0/PI
         WRITE(wunit,'("# =======================================================================================================")')
         DO i=this%nstart,this%ntraj
            WRITE(wunit,'(1X,I10,F16.7,2(6F16.7))') i,this%trajs(i)%E,this%trajs(i)%init_r,this%trajs(i)%init_p
         END DO
         CLOSE(wunit)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Outputfile generated: ",this%output_file)
#endif
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE GENERATE_TRAJS_FROM_FILE_INITSURFACEGROW
!###########################################################
!# SUBROUTINE: INITIALIZE_DIATOMICGROW #########################
!###########################################################
!> @brief
!! Initializes type atom
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INITIALIZE_DIATOMICGROW(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(DiatomicGrow),INTENT(OUT):: this
   ! Local variables
   INTEGER(KIND=4),PARAMETER :: dimens=6
   ! Run section
   ALLOCATE(this%init_r(dimens))
   ALLOCATE(this%init_p(dimens))
   ALLOCATE(this%init_qn(3))
   ALLOCATE(this%final_qn(3))
   this%stat="Dummy"
   RETURN
END SUBROUTINE INITIALIZE_DIATOMICGROW
!##################################################################################
!# SUBROUTINE: INITIALIZE_INITDIATOMIC ################################################
!##################################################################################
!> @brief
!! Reads from an input file enough data to generate initial conditions.
!! for a batch of atoms. Adapted to atom/surface systems in surface Grow
!
!> @param[out] this - Initial conditions for a batch of molecules
!> @param[in] filename - Name of the input file
!------------------------------------------------------------------------------
SUBROUTINE INITIALIZE_INITSURFACEGROW(this,filename)
   IMPLICIT NONE
   ! I/O variables
   CLASS(InitSurfaceGrowDiatomic),INTENT(OUT):: this
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: filename
   ! Local variables
   CHARACTER(LEN=*), PARAMETER :: routinename = "INITIALIZE_INITDIATOMIC: "
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
      CASE('Molecules')
         ! do nothing
      CASE DEFAULT
         WRITE(0,*) "INITIALIZE_INITDIATOMIC ERR: wrong kind of initial conditions"
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
   ! get internal state
   CALL AOT_TABLE_OPEN(L=conf,parent=inicond_table,thandle=magnitude_table,key='internalState')
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,key='v',val=auxint)
   this%init_qn(1)=auxint
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,key='J',val=auxint)
   this%init_qn(2)=auxint
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,key='mJ',val=auxint)
   SELECT CASE(btest(ierr,aoterr_WrongType))
      CASE(.false.)
         this%fixed_theta=.true.
         this%init_qn(3)=auxint
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,'Quantum state (v,J,mJ): ',this%init_qn(:))
#endif
      CASE(.true.)
         this%fixed_theta=.false.
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,'Quantum state (v,J,mJ average): ',this%init_qn(1:2))
#endif
   END SELECT
   CALL AOT_TABLE_CLOSE(L=conf,thandle=magnitude_table)
   ! get internal energy
   CALL AOT_TABLE_OPEN(L=conf,parent=inicond_table,thandle=magnitude_table,key='internalEnergy')
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=1,val=auxreal)
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=2,val=auxstring)
   CALL this%evirot%READ(auxreal,trim(auxstring))
   CALL this%evirot%TO_STD()
   CALL AOT_TABLE_CLOSE(L=conf,thandle=magnitude_table)
   SELECT CASE(this%evirot%getvalue()>0.D0+1.d-8)
      CASE(.TRUE.)
         this%is_classic=.FALSE.
      CASE(.FALSE.)
         this%is_classic=.TRUE.
   END SELECT
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Is this classical?: ",this%is_classic)
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
         WRITE(0,*) "INITIALIZE_INITDIATOMIC ERR: wrong incidence angle: ",this%vz_angle%getvalue()
         CALL EXIT(1)
      CASE(.false.)
         ! do nothing
   END SELECT
   ! get vibrational potential
   CALL AOT_TABLE_OPEN(L=conf,parent=inicond_table,thandle=magnitude_table,key='vibrationalFunction')
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,key='kind',val=auxstring)
   SELECT CASE(trim(auxstring))
      CASE('Numerical')
         CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,key='source',val=auxstring)
         CALL this%vibrpot%INITIALIZE(trim(auxstring))
         CALL this%vibrpot%SHIFTPOT()
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,'Numerical Vacuumpot loaded')
         CALL VERBOSE_WRITE(routinename,'Equilibrium distance (au): ',this%vibrpot%getreq())
         CALL VERBOSE_WRITE(routinename,'Potential at minimum, prior shifting (au): ',this%vibrpot%getpotmin())
         CALL VERBOSE_WRITE(routinename,'Force Constant (d2V(r)/dr2 at Req (au): ',this%vibrpot%getForceConstant())
#endif
      CASE DEFAULT
         WRITE(0,*) "INITIALIZE_INITDIATOMIC ERR: vibrational function kind not implemented"
         WRITE(0,*) "Implemented ones: Numerical"
         WRITE(0,*) 'Warning: case sensitive'
         CALL EXIT(1)
   END SELECT
   CALL AOT_TABLE_CLOSE(L=conf,thandle=magnitude_table)
   ! get integration parameters
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=inicond_table,key='extrapolation',val=auxstring)
   this%extrapol=trim(auxstring)
   CALL AOT_TABLE_OPEN(L=conf,parent=inicond_table,thandle=magnitude_table,key='timeStep')
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=1,val=auxreal)
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=2,val=auxstring)
   CALL this%delta_t%READ(auxreal,trim(auxstring))
   CALL this%delta_t%TO_STD()
   CALL AOT_TABLE_CLOSE(L=conf,thandle=magnitude_table)
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=inicond_table,key='precision',val=this%eps)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'Extrapolation used: '//this%extrapol)
   CALL VERBOSE_WRITE(routinename,'Initial time step (a.u.): ',this%delta_t%getvalue())
   CALL VERBOSE_WRITE(routinename,'Precision during integration (dimensionless factor): ',this%eps)
#endif
   ! get direction angle
   CALL AOT_TABLE_OPEN(L=conf,parent=inicond_table,thandle=magnitude_table,key='directionAngle')
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=1,val=auxreal)
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=2,val=auxstring)
   CALL this%vpar_angle%READ(auxreal,trim(auxstring))
   CALL this%vpar_angle%TO_STD()
   CALL AOT_TABLE_CLOSE(L=conf,thandle=magnitude_table)
   SELECT CASE(this%vz_angle%getvalue() < 0.D0)
      CASE(.true.)
         WRITE(0,*) "INITIALIZE_INIDIATOMIC ERR: parallele angle lower than 0.0 deg."
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
      WRITE(0,*) "INITIALIZE_INITDIATOMIC ERR: X impact parameter outside range 0-1"
      CALL EXIT(1)
   ELSE IF(((this%impact_y > 1.D0).OR.(this%impact_y < 0.D0)).AND.(.NOT.this%control_posy))THEN
      WRITE(0,*) "INITIALIZE_INITDIATOMIC ERR: Y impact parameter outside range 0-1"
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
END SUBROUTINE INITIALIZE_INITSURFACEGROW
!####################################################################
!# SUBROUTINE: GENERATE_TRAJS_INITSURFACEGROW #######################
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
SUBROUTINE GENERATE_TRAJS_INITSURFACEGROW(this,thispes)
   IMPLICIT NONE
   ! I/O variables
   CLASS(InitSurfaceGrowDiatomic),INTENT(INOUT):: this
   CLASS(PES),INTENT(IN):: thispes
   ! IMPORTANT: unit used to write
   INTEGER(KIND=4),PARAMETER :: wunit=923
   ! Local variables
   INTEGER :: i ! counters
   REAL(KIND=8):: delta,alpha,Enorm,masa,mu,Eint
   REAL(KIND=8),DIMENSION(9):: random_kernel
   REAL(KIND=8):: rnd_delta,rnd_phi,rnd_theta,rnd_eta,rnd_sign,rnd_mu
   REAL(KIND=8):: eta
   REAL(KIND=8):: L,L_theta
   real(kind=8),dimension(12):: atomCoordPhaseSpace,molecCoordPhaseSpace
   ! Parameters
   character(len=*),parameter:: formatFile='("# Format: id/Etot,Eint/X,Y,Z,R(a.u.)/THETA,PHI(rad)/Px,Py,Pz,pr,ptheta,pphi(a.u.)")'
   character(len=*),parameter:: routinename = "GENERATE_TRAJS_ATOMS: "
   ! YIPPIEE KI YAY !! -------------------
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"New set of trajectories")
   CALL VERBOSE_WRITE(routinename,"Allocating trajs: ", this%ntraj)
#endif
   ALLOCATE(DiatomicGrow::this%trajs(this%ntraj))
   SELECT CASE(this%is_classic)
      CASE(.TRUE.)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Calculating period in R is not meanful")
#endif
         ! body
      CASE(.FALSE.)
         CALL this%SET_PERIOD()
   END SELECT
   delta = this%vpar_angle%getvalue()
   alpha = this%vz_angle%getvalue()
   masa = sum(system_mass(1:2))
   mu = product(system_mass(1:2))/masa
   Enorm = this%E_norm%getvalue()
   Eint = this%evirot%getvalue()
   L=dsqrt(dfloat(this%init_qn(2)*(this%init_qn(2)+1)))
   DO i=1,this%ntraj
      CALL RANDOM_NUMBER(random_kernel(:))
      IF((this%control_posX).AND.(.NOT.this%control_posY)) THEN
         CALL this%trajs(i)%INITIALIZE()
         this%trajs(i)%r(1)=random_kernel(1)
         this%trajs(i)%r(2)=this%impact_y
      ELSE IF ((.NOT.this%control_posX).AND.(this%control_posY)) THEN
         CALL this%trajs(i)%INITIALIZE()
         this%trajs(i)%r(1)=this%impact_x
         this%trajs(i)%r(2)=random_kernel(2)
      ELSE IF ((this%control_posX).AND.(this%control_posY)) THEN
         CALL this%trajs(i)%INITIALIZE()
         this%trajs(i)%r(1:2)=random_kernel(1:2)
      ELSE
         CALL this%trajs(i)%INITIALIZE()
         this%trajs(i)%r(1)=this%impact_x
         this%trajs(i)%r(2)=this%impact_y
      END IF
      ! CENTER OF MASS COORD ---------------------------------------
      ! Set center of mass coordinates
      this%trajs(i)%r(3) = this%init_z%getvalue()
      this%trajs(i)%r(1:2) = system_surface%surf2cart(this%trajs(i)%r(1:2))
      ! center of mass momenta DOFS. For a init.raw.dat output, should be zero.
      this%trajs(i)%p(1)=0.d0
      this%trajs(i)%p(2)=0.d0
      this%trajs(i)%p(3)=0.d0
      ! center of mass energy. For init.raw.dat output, should be zero.
      this%trajs(i)%Ecm = 0.d0
      ! INTERNAL COORD --------------------------------------------
      ! Set random kernel
      rnd_delta=random_kernel(3)
      rnd_phi=random_kernel(4)
      rnd_theta=random_kernel(5)
      rnd_eta=random_kernel(6)
      rnd_mu=random_kernel(7)
      SELECT CASE(random_kernel(8)<0.5d0)
         CASE(.true.)
            rnd_sign=1.d0
         CASE DEFAULT
            rnd_sign=-1.d0
      END SELECT
      SELECT CASE(L==0.d0)
        CASE(.true.)
            L_theta=0.d0
        CASE(.false.)
            L_theta=dacos(dfloat(this%init_qn(3))/L)
      END SELECT
      SELECT CASE(this%is_classic)
         CASE(.TRUE.)
            this%trajs(i)%r(4)=this%vibrpot%getreq()                                 ! r
            this%trajs(i)%p(4)=0.D0                                                  ! pr
            SELECT CASE(this%fixed_theta)
               CASE(.TRUE.)
                  ! Get internal coordinates
                  this%trajs(i)%r(5)=dacos(dsin(L_theta)*dcos(PI*rnd_mu))            ! theta
                  this%trajs(i)%r(6)=2.D0*PI*rnd_phi                                 ! phi
                  ! Get internal momenta
                  eta=dacos(dcos(L_theta)/dsin(this%trajs(i)%r(5)))
                  this%trajs(i)%p(5)=rnd_sign*L*dsin(eta)                                    ! ptheta
                  this%trajs(i)%p(6)=L*dcos(eta)*dsin(this%trajs(i)%r(5))  ! pphi
               CASE(.FALSE.)
                  ! Get internal coordinates
                  this%trajs(i)%r(5)=dacos(2.D0*rnd_theta-1.d0)                      ! theta
                  this%trajs(i)%r(6)=2.D0*PI*rnd_phi                                 ! phi
                  ! Get internal momenta
                  eta=2.D0*PI*rnd_eta
                  this%trajs(i)%p(5)=-L*dsin(eta)                                    ! ptheta
                  this%trajs(i)%p(6)=L*dcos(eta)*dsin(this%trajs(i)%r(5))  ! pphi
            END SELECT
            ! Internal energy
            this%trajs(i)%Eint=((L/this%vibrpot%getreq())**2.D0)/(2.D0*mu)
         CASE(.FALSE.)
            this%trajs(i)%r(4)=this%r_t%getvalue(rnd_delta*this%period)             ! r
            this%trajs(i)%p(4)=this%pr_t%getvalue(rnd_delta*this%period)            ! pr
            SELECT CASE(this%fixed_theta)
               CASE(.TRUE.)
                  ! Get internal coordinates
                  this%trajs(i)%r(5)=dacos(dsin(L_theta)*dcos(PI*rnd_mu))           ! theta
                  this%trajs(i)%r(6)=2.D0*PI*rnd_phi                                ! phi
                  ! Get internal momenta
                  eta=dacos(dcos(L_theta)/dsin(this%trajs(i)%r(5)))
                  this%trajs(i)%p(5)=rnd_sign*L*dsin(eta)                                   ! ptheta
                  this%trajs(i)%p(6)=L*dcos(eta)*dsin(this%trajs(i)%r(5)) ! pphi
               CASE(.FALSE.)
                  ! Get internal coordinates
                  this%trajs(i)%r(5)=dacos(2.D0*rnd_theta-1.d0)                     ! theta
                  this%trajs(i)%r(6)=2.D0*PI*rnd_phi                                ! phi
                  ! Get internal momenta
                  eta=2.D0*PI*rnd_eta
                  this%trajs(i)%p(5)=-L*dsin(eta)                                   ! ptheta
                  this%trajs(i)%p(6)=L*dcos(eta)*dsin(this%trajs(i)%r(5)) ! pphi
            END SELECT
            ! Internal energy
            this%trajs(i)%Eint=(this%trajs(i)%p(4)**2.D0)/(2.D0*mu)+this%rovibrpot%getpot(this%trajs(i)%r(4))
      END SELECT
      ! Total energy
      this%trajs(i)%E=this%trajs(i)%Ecm+this%trajs(i)%Eint
      this%trajs(i)%init_E=this%trajs(i)%Ecm+this%trajs(i)%Eint
      ! STORE INITIAL VALUES .....................................
      this%trajs(i)%init_r = this%trajs(i)%r
      this%trajs(i)%init_p = this%trajs(i)%p
   END DO
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Initial trajectory: ",this%nstart)
   CALL VERBOSE_WRITE(routinename,"Final trajectory: ",this%ntraj)
#endif
   ! Print if the option was given
   SELECT CASE(this%control_out)
      CASE(.TRUE.)
         OPEN(wunit,FILE=this%output_file,STATUS="replace")
         WRITE(wunit,'("# FILE CREATED BY : GENERATE_TRAJS_ATOMS =================================================================")')
         WRITE(wunit,formatFile)
         WRITE(wunit,'("# Perpendicular Energy (a.u.) / (eV): ",F10.5," / ",F10.5)') this%E_norm%getvalue(),&
            this%E_norm%getvalue()*au2ev
         WRITE(wunit,'("# Initial center of mass energy (a.u.) / (eV) : ",F10.5," / ",F10.5)') Enorm/(DSIN(alpha)**2.D0),&
            (Enorm/(DSIN(alpha)**2.D0))*au2ev
         WRITE(wunit,'("# Initial internal energy (a.u.) / (eV) : ",F10.5," / ",F10.5)') Eint,Eint*au2ev
         WRITE(wunit,'("# Initial Rovibrational state: ",3(I5,1X))') this%init_qn(:)
         WRITE(wunit,'("# MASS (a.u.) / proton_mass : ",F10.5," / ",F10.5)')  masa,masa/pmass2au
         WRITE(wunit,'("# Reduced MASS (a.u.) / proton_mass : ",F10.5," / ",F10.5)') mu,mu/pmass2au
         WRITE(wunit,'("# Incidence angle (deg): ",F10.5)')  this%vz_angle%getvalue()*180.D0/PI
         WRITE(wunit,'("# Parallel velocity direction (deg): ",F10.5)') this%vpar_angle%getvalue()*180.D0/PI
         IF ((this%control_posX).AND.(.NOT.this%control_posY)) THEN
            WRITE(wunit,'("# Random X impact parameter")')
            WRITE(wunit,*) "# Seed used: ",system_iSeed(:)
         ELSE IF ((.NOT.this%control_posX).AND.(this%control_posY)) THEN
            WRITE(wunit,'("# Random Y impact parameter")')
            WRITE(wunit,*) "# Seed used: ",system_iSeed(:)
         ELSE IF ((this%control_posX).AND.(this%control_posY)) THEN
            WRITE(wunit,'("# Random X and Y impact parameters")')
            WRITE(wunit,*) "# Seed used: ",system_iSeed(:)
         ELSE
            WRITE(wunit,'("# X, Y values are not random numbers")')
            WRITE(wunit,*) "# Seed used: ",system_iSeed(:)
         END IF
         WRITE(wunit,'("# =======================================================================================================")')
         DO i=this%nstart,this%ntraj
            WRITE(wunit,'(1X,I10,2(F16.7),2(6F16.7))') i,this%trajs(i)%E,this%trajs(i)%Eint,&
               this%trajs(i)%init_r,this%trajs(i)%init_p
         END DO
         CLOSE(wunit)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Outputfile generated: ",this%output_file)
#endif
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   open( unit=wunit,file='init.raw.dat',status='new',action='write' )
   do i=this%nstart,this%ntraj
      write(wunit,*)
      write(wunit,*)
      write(wunit,*)
      write(wunit,*)
      molecCoordPhaseSpace(1:6)=this%trajs(i)%init_r(:)
      molecCoordPhaseSpace(1:6)=this%trajs(i)%init_p(:)
      atomCoordPhaseSpace(:)=from_molecular_to_atomic_phaseSpace( molecCoordPhaseSpace )
      write(wunit,'(F16.7,1X,F16.7,1X,F16.7)') atomCoordPhaseSpace(1:3)
      write(wunit,'(F16.7,1X,F16.7,1X,F16.7)') atomCoordPhaseSpace(4:6)
      write(wunit,*)
      atomCoordPhaseSpace(7:9)=atomCoordPhaseSpace(7:9)*dsqrt(pmass2au)/system_mass(1)
      atomCoordPhaseSpace(10:12)=atomCoordPhaseSpace(10:12)*dsqrt(pmass2au)/system_mass(2)
      write(wunit,'(F16.7,1X,F16.7,1X,F16.7)') atomCoordPhaseSpace(7:9)
      write(wunit,'(F16.7,1X,F16.7,1X,F16.7)') atomCoordPhaseSpace(10:12)
   enddo
   close(wunit)
   return
END SUBROUTINE GENERATE_TRAJS_INITSURFACEGROW
!###############################################################
!# SUBROUTINE : TIME_DERIVS_INITSURFACEGROW ####################
!###############################################################
!> @brief
!! Gives dzdt at z and t values from Hamilton equations of motion
!
!> @param[in] this - Provides some data
!> @param[in] z - array of positions and momenta. z(1) -> position, z(2) -> momenta
!> @param[out] dzdt - time derivatives of position and momenta
!> @param[out] fin - controls errors
!--------------------------------------------------------------
SUBROUTINE TIME_DERIVS_INITSURFACEGROW(this,z,dzdt,fin)
   IMPLICIT NONE
   ! I/O variables
   CLASS(InitSurfaceGrowDiatomic),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: z
   REAL(KIND=8),DIMENSION(2),INTENT(OUT) :: dzdt
   LOGICAL,INTENT(OUT) :: fin
   ! Local variables
   CHARACTER(LEN=*),PARAMETER :: routinename = "TIME_DERIVS_INISURFACEGROW: "
   ! ROCK THE CASBAH ! ---------------------
   SELECT CASE(this%rovibrpot%is_allowed(z(1)))
      CASE(.FALSE.)
         fin = .TRUE.
      CASE(.TRUE.)
         fin=.FALSE.
         dzdt(1)=z(2)*sum(system_mass(1:2))/product(system_mass(1:2))
         dzdt(2)=-this%rovibrpot%getderiv(z(1))
   END SELECT
   RETURN
END SUBROUTINE TIME_DERIVS_INITSURFACEGROW
!#########################################################################################
!# SUBROUTINE: MMID_INITSURFACEGROW ######################################################
!#########################################################################################
!> @brief
!! Modified midpoint method. Given a vector of components @f$y_{i}(x)@f$, this subroutine
!! can advance @f$y_{i}(x+H)@f$ by a sequence of n substeps. It is used as an integrator in
!! the more powerful Burlisch-Stoer technique.
!
!> @details
!! - Adapted from Numerical recipes FORTRAN 77
!! - Adapted for the integration of the equations of motion of an atom
!
!> @param[in] this - We need information stored in this variable
!> @param[in] y(1:6) - array with positions and momenta
!> @param[in] dydx - derivatives of y
!> @param[in] xs - X at which y and dydx are meassured
!> @param[in] htot - Total step size
!> @param[in] nstep - substeps to be used
!> @param[out] yout - output array
!> @param[in,out] switch - Controls if there was a problem calculating the potential at
!!                         a given value
!
!> @see Fortran 77 numerical recipes
!-----------------------------------------------------------------------------------------
SUBROUTINE MMID_INITSURFACEGROW(this,y,dydx,xs,htot,nstep,yout,switch)
   IMPLICIT NONE
   ! I/O variables
   CLASS(InitSurfaceGrowDiatomic),INTENT(IN):: this
   INTEGER,INTENT(IN):: nstep
   REAL(KIND=8),DIMENSION(2), INTENT(IN) :: y,dydx
   REAL(KIND=8),DIMENSION(2), INTENT(OUT) :: yout
   REAL(KIND=8),INTENT(IN) :: xs,htot
   LOGICAL,INTENT(INOUT) :: switch
   ! Local variables
   REAL(KIND=8),DIMENSION(6) :: ym,yn
   INTEGER :: i,n ! counters
   INTEGER,PARAMETER :: nvar = 2
   REAL(KIND=8) :: h,h2,swap,x
   ! ROCK THE CASBAH !!! ---------------------
   h=htot/DFLOAT(nstep) ! Stepsize this trip.
   DO i=1,nvar
      ym(i)=y(i)
      yn(i)=y(i)+h*dydx(i) ! First step.
   END DO
   x=xs+h
   CALL this%TIME_DERIVS(yn,yout,switch) ! Will use yout for temporary storage of derivatives.
   SELECT CASE(switch)
      CASE(.TRUE.)
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   h2=2.D0*h
   DO n=2,nstep ! General step.
      DO i=1,nvar
         swap=ym(i)+h2*yout(i)
         ym(i)=yn(i)
         yn(i)=swap
      END DO
      x=x+h
      CALL this%TIME_DERIVS(yn,yout,switch)
      SELECT CASE(switch)
         CASE(.TRUE.)
            RETURN
         CASE(.FALSE.)
            ! do nothing
      END SELECT
   END DO
   DO i=1,nvar
      yout(i)=0.5D0*(ym(i)+yn(i)+h*yout(i))
   END DO
   RETURN
END SUBROUTINE MMID_INITSURFACEGROW
!##################################################################################################
!# SUBROUTINE: RZEXTR_INITSURFACEGROW #############################################################
!##################################################################################################
!> @brief
!! - A part of the Burlich-Stoer algorithm. Uses diagonal rational function extrapolation.
!
!>@details
!! Taken from Numerical recipes in Fortran 77
!> @see pzextr
!--------------------------------------------------------------------------------------------------
SUBROUTINE RZEXTR_INITSURFACEGROW(this,iest,xest,yest,yz,dy,nv)
   IMPLICIT NONE
   ! I/O variables
   CLASS(InitSurfaceGrowDiatomic),INTENT(IN):: this
   INTEGER,INTENT(IN) :: iest, nv
   REAL(KIND=8),INTENT(IN) :: xest
   REAL(KIND=8),DIMENSION(nv), INTENT(IN) :: yest
   REAL(KIND=8),DIMENSION(nv), INTENT(OUT) :: dy, yz
   ! Local variables
   INTEGER, PARAMETER :: IMAX = 13
   INTEGER, PARAMETER :: NMAX = 50
   INTEGER :: j,k
   REAL(KIND=8), DIMENSION(NMAX,IMAX) :: d
   REAL(KIND=8), DIMENSION(IMAX) :: fx, x
   REAL(KIND=8) :: b,b1,c,ddy,v,yy
   SAVE d,x
   ! HEY, HO!, LETS GO!! ------------------------
   x(iest)=xest
   !Save current independent variable.
   IF(iest.EQ.1) THEN
      DO j=1,nv
         yz(j)=yest(j)
         d(j,1)=yest(j)
         dy(j)=yest(j)
      END DO
   ELSE
      DO k=1,iest-1
         fx(k+1)=x(iest-k)/xest
      END DO
      DO j=1,nv
         ! Evaluate next diagonal in tableau.
         yy=yest(j)
         v=d(j,1)
         c=yy
         d(j,1)=yy
         DO k=2,iest
            b1=fx(k)*v
            b=b1-c
            IF(b.NE.0.) THEN
               b=(c-v)/b
               ddy=c*b
               c=b1*b
            ELSE
               ! Care needed to avoid division by 0.
               ddy=v
            END IF
            IF (k.NE.iest) v=d(j,k)
            d(j,k)=ddy
            yy=yy+ddy
         END DO
         dy(j)=ddy
         yz(j)=yy
      END DO
   END IF
   RETURN
END SUBROUTINE RZEXTR_INITSURFACEGROW
!############################################################################################
!# SUBROUTINE : PZEXTR_INITSURFACEGROW ######################################################
!############################################################################################
!!> @brief
!! Uses polynomial extrapolation to evaluate nv functions at x = 0 by fitting a polynomial to a
!! sequence of estimates with progressively smaller values x = xest, and corresponding function
!! vectors yest(1:nv).
!
!> @details
!! - Extrapolated function values are output as yz(1:nv), and their estimated error is output as dy(1:nv).
!! - Maximum expected value of iest is IMAX; of nv is NMAX.
!! - Taken from Numerical recipes
!
!> @see Numerical recipes in fortran 77
!--------------------------------------------------------------------------------------------
SUBROUTINE PZEXTR_INITSURFACEGROW(this,iest,xest,yest,yz,dy,nv)
   IMPLICIT NONE
   ! I/O variables
   CLASS(InitSurfaceGrowDiatomic),INTENT(IN):: this
   INTEGER,INTENT(IN) :: iest, nv
   REAL(KIND=8),INTENT(IN) :: xest
   REAL(KIND=8),DIMENSION(nv),INTENT(IN) :: yest
   REAL(KIND=8),DIMENSION(nv),INTENT(OUT) :: dy, yz
   ! Local Variables
   INTEGER,PARAMETER :: IMAX = 13
   INTEGER,PARAMETER :: NMAX = 50
   INTEGER :: j, k1 ! counters
   REAL(KIND=8),DIMENSION(NMAX) :: d
   REAL(KIND=8),DIMENSION(IMAX) :: x
   REAL(KIND=8),DIMENSION(NMAX,IMAX) :: qcol
   REAL(KIND=8):: delta,f1,f2,q
   ! ROCK THE CASBAH !!!! --------------
   SAVE qcol,x
   x(iest)=xest !Save current independent variable.
   DO j=1,nv
      dy(j)=yest(j)
      yz(j)=yest(j)
   END DO
   IF (iest.eq.1) THEN ! Store 1st estimate in 1st column.
      DO j=1,nv
         qcol(j,1)=yest(j)
      END DO
   ELSE
      DO j=1,nv
         d(j)=yest(j)
      END DO
      DO k1=1,iest-1
         delta=1./(x(iest-k1)-xest)
         f1=xest*delta
         f2=x(iest-k1)*delta
         DO j=1,nv ! Propagate tableau 1 diagonal more.
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
         END DO
      END DO
      DO j=1,nv
         qcol(j,iest)=dy(j)
      END DO
   END IF
   RETURN
END SUBROUTINE PZEXTR_INITSURFACEGROW
!###############################################################################################
!# SUBROUTINE : BSSTEP_INITSURFACEGROW #########################################################
!###############################################################################################
!> @brief
!! Bulirsch-Stoer step with monitoring of local truncation error to ensure accuracy and adjust
!! stepsize.
!
!> @details
!! - Input are the dependent variable vector y(1:nv) and its derivative dydx(1:nv)
!!   at the starting value of the independent variable x. Also input are the stepsize to be attempted
!!   htry, the required accuracy eps, and the vector yscal(1:nv) against which the
!!   error is scaled. On output, y and x are replaced by their new values, hdid is the stepsize
!!   that was actually accomplished, and hnext is the estimated next stepsize. derivs is the
!!   user-supplied subroutine that computes the right-hand side derivatives. Be sure to set htry
!!   on successive steps to the value of hnext returned from the previous step, as is the case
!!   if the routine is called by odeint.
!!
!! - Parameters: NMAX is the maximum value of nv; KMAXX is the maximum row number used
!!   in the extrapolation; IMAX is the next row number; SAFE1 and SAFE2 are safety factors;
!!   REDMAX is the maximum factor used when a stepsize is reduced, REDMIN the minimum;
!!   TINY prevents division by zero; 1/SCALMX is the maximum factor by which a stepsize can
!!   be increased.
!! - Adapted to integrate equations of motion for an atom
!
!> @todo
!! - Should be generalized (pending task)
!! - Should use dynamic memory (pending task)
!------------------------------------------------------------------------------------------------
SUBROUTINE BSSTEP_INITSURFACEGROW(this,y,dydx,x,htry,eps,yscal,hdid,hnext,switch)
   IMPLICIT NONE
   ! I/O variables
   CLASS(InitSurfaceGrowDiatomic),INTENT(IN) :: this
   REAL(KIND=8), INTENT(IN) :: eps     ! required accuracy
   REAL(KIND=8), INTENT(IN) ::  htry   ! step to try
   REAL(KIND=8), DIMENSION(2) :: yscal ! factors to scale error
   REAL(KIND=8), DIMENSION(2), INTENT(IN) :: dydx
   REAL(KIND=8), DIMENSION(2), INTENT(INOUT) :: y ! initial/final values for: X,Y,Z,Px,Py,Pz (in this order)
   REAL(KIND=8), INTENT(INOUT) :: x
   LOGICAL, INTENT(INOUT) :: switch ! .TRUE. if potential could not be calculated
   REAL(KIND=8), INTENT(OUT) :: hdid  ! step actually used
   REAL(KIND=8), INTENT(OUT) :: hnext ! guess of the next step
   ! Parameters for this routine
   INTEGER, PARAMETER :: nv = 2
   REAL(KIND=8),PARAMETER :: SAFE1 = 0.25D0
   REAL(KIND=8),PARAMETER :: SAFE2 = 0.7D0
   REAL(KIND=8),PARAMETER :: TINY = 1.D-30
   REAL(KIND=8),PARAMETER :: SCALMX = 0.1D0
   REAL(KIND=8),PARAMETER :: REDMIN = 0.7D0
   REAL(KIND=8),PARAMETER :: REDMAX = 1.D-5
   INTEGER,PARAMETER :: NMAX = 50
   INTEGER,PARAMETER :: KMAXX = 8
   INTEGER,PARAMETER :: IMAX = KMAXX+1
   CHARACTER(LEN=21), PARAMETER :: routinename = "BSSTEP_INITDIATOMIC: "
   ! Local variables
   INTEGER, DIMENSION(IMAX) :: nseq
   REAL(KIND=8), DIMENSION(KMAXX) :: err
   REAL(KIND=8), DIMENSION(NMAX) :: yerr, ysav, yseq
   REAL(KIND=8), DIMENSION(IMAX) :: a
   REAL(KIND=8), DIMENSION(KMAXX,KMAXX) :: alf
   INTEGER :: i,iq,k,kk,km,kmax,kopt
   REAL(KIND=8) :: eps1,epsold,errmax,fact,h,red,scale,work,wrkmin,xest, xnew
   LOGICAL :: first,reduct
   SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
   DATA first/.TRUE./,epsold/-1./
   DATA nseq /2,4,6,8,10,12,14,16,18/
   ! HEY, HO! LET'S GO !!! -----------------------------
   switch = .FALSE.
   IF (eps.NE.epsold) THEN !A new tolerance, so reinitialize.
      hnext = -1.D29  ! Impossible values.
      xnew  = -1.D29
      eps1 = SAFE1*eps
      a(1)=nseq(1)+1  ! Compute work coecients Ak.
      DO k=1,KMAXX
         a(k+1)=a(k)+nseq(k+1)
      END DO
      DO iq=2,KMAXX ! Compute alpha(k,q)
         DO k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.)*(2*k+1)))
         END DO
      END DO
      epsold = eps
      DO kopt=2,KMAXX-1 ! Determine optimal row number for convergence
         if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt)) goto 1
      END DO
1 kmax = kopt
   END IF
   h = htry
   DO i=1,nv ! Save the starting values.
      ysav(i) = y(i)
   END DO
   IF (h.ne.hnext.or.x.ne.xnew) THEN ! A new stepsize or a new integration: re-establish
      first=.true.              ! the order window.
      kopt=kmax
   END IF
   reduct=.false.
2 DO k=1,kmax             ! Evaluate the sequence of modified midpoint
      xnew=x+h  ! integrations.
      CALL this%MMID(ysav,dydx,x,h,nseq(k),yseq,switch)
      IF(switch) RETURN
      xest=(h/nseq(k))**2.D0 ! Squared, since error series is even.
      IF(this%extrapol.EQ."Rational") THEN
         CALL this%RATIONAL_EXTRAPOL(k,xest,yseq,y,yerr,nv) ! Perform extrapolation with Rational funcions
      ELSE IF (this%extrapol.EQ."Polinomi") THEN
         CALL this%POLINOM_EXTRAPOL(k,xest,yseq,y,yerr,nv) ! Perform extrapolation with Polinoms
      ELSE
         WRITE(0,*) "BSSTEP_ATOM ERR: Wrong keyword for extrapolation variable"
         CALL EXIT(1)
      END IF
      IF(k.NE.1) THEN ! Compute normalized error estimate
         errmax=TINY
         DO i=1,nv
            errmax=MAX(errmax,DABS(yerr(i)/yscal(i)))
         END DO
         errmax=errmax/eps ! Scale error relative to tolerance.
         km=k-1
         err(km)=(errmax/SAFE1)**(1./(2*km+1))
      END IF
      IF(k.NE.1.AND.(k.GE.kopt-1.OR.first)) THEN ! In order window.
         IF(errmax.lt.1.) GO TO 4  ! Converged.
         IF (k.eq.kmax.or.k.eq.kopt+1) THEN ! Check for possible stepsize reduction.
            red=SAFE2/err(km)
            GO TO 3
         ELSE IF (k.EQ.kopt) THEN
            IF(alf(kopt-1,kopt).LT.err(km)) THEN
               red=1./err(km)
               GO TO 3
            END IF
         ELSE IF(kopt.EQ.kmax) THEN
            IF(alf(km,kmax-1).LT.err(km)) THEN
               red=alf(km,kmax-1)*SAFE2/err(km)
               GO TO 3
            END IF
         ELSE IF (alf(km,kopt).LT.err(km)) THEN
            red=alf(km,kopt-1)/err(km)
            GO TO 3
         END IF
      END IF
   END DO
3 red=MIN(red,REDMIN)       ! Reduce stepsize by at least REDMIN and at
   red=MAX(red,REDMAX) ! most REDMAX.
   h=h*red
   reduct=.TRUE.
   GO TO 2 ! Try again.
4 x=xnew ! Successful step taken.
   hdid=h
   first=.FALSE.
   wrkmin=1.D35 ! Compute optimal row for convergence and
   DO kk=1,km   ! corresponding stepsize.
      fact= MAX(err(kk),SCALMX)
      work=fact*a(kk+1)
      IF(work.LT.wrkmin) THEN
         scale=fact
         wrkmin=work
         kopt=kk+1
      END IF
   END DO
   hnext=h/scale
   IF(kopt.GE.k.AND.kopt.NE.kmax.AND..NOT.reduct) THEN ! Check for possible order increase,
                         ! but not if stepsize
                         ! was just reduced.
      fact=MAX(scale/alf(kopt-1,kopt),SCALMX)
      IF(a(kopt+1)*fact.LE.wrkmin) THEN
         hnext=h/fact
         kopt=kopt+1
      END IF
   END IF
   RETURN
END SUBROUTINE BSSTEP_INITSURFACEGROW
!###########################################################
!# SUBROUTINE: SET_PERIOD_INITDIATOMIC
!###########################################################
!> @brief
!! Sets period of vibrational motion
!
!> @warning
!! - Invoking this subroutine is only meanful while doing
!!   @b quasiclassical dynamics, i.e., internal energy higher
!!   than 0.D0.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_PERIOD_INITSURFACEGROW(this)
   ! Initial declarations
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   ! I/O variables
   CLASS(InitSurfaceGrowDiatomic),INTENT(INOUT):: this
   ! Local variables
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: x,f
   REAL(KIND=8) :: init_r,init_pr
   REAL(KIND=8) :: erot
   REAL(KIND=8) :: turn1,turn2
   REAL(KIND=8) :: vturn1,vturn2
   REAL(KIND=8) :: t,dt,init_t,init_Eint,dt_did,dt_next
   REAL(KIND=8) :: time_turn2, time_cycle, time_finish
   INTEGER(KIND=4) :: cycles
   REAL(KIND=8) :: Eint ! internal energy, potential
   LOGICAL :: switch
   LOGICAL :: reached_turn2, onecycle_completed
   INTEGER(KIND=4) :: n_req_passedby
   REAL(KIND=8),DIMENSION(2) :: z,dzdt ! stores DOFS (r,Pr) and time derivatives
   REAL(KIND=8),DIMENSION(2) :: s ! scaling vector
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: pr_t,r_t,time
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: aux
   INTEGER(KIND=4) :: nold
   REAL(KIND=8) :: mu
   INTEGER(KIND=4) :: i
   TYPE(Vacuumpot):: fdummy
   CHARACTER(LEN=25),PARAMETER :: routinename="SET_PERIOD_INITDIATOMIC: "
   ! Run section
   mu = product(system_mass(1:2))/sum(system_mass(1:2))
   erot=dfloat(this%init_qn(2)*(this%init_qn(2)+1))/(2.D0*mu)
   ALLOCATE(x(this%vibrpot%n))
   ALLOCATE(f(this%vibrpot%n))
   DO i = 1, this%vibrpot%n
      x(i)=this%vibrpot%rpot%x(i)
      f(i)=this%vibrpot%rpot%f(i)+erot/(x(i)**2.D0)
   END DO
   CALL this%rovibrpot%INITIALIZE_DIRECT(x,f,0.D0)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Initialize effective potential, rovibrpot")
   CALL VERBOSE_WRITE(routinename,'Effective potential is: Veff(r) = V(r)+(L**2)/(2*mu*r**2)')
   CALL this%rovibrpot%PLOT(1000,"effectivepot.dat")
   CALL VERBOSE_WRITE(routinename,'Effective potential stored at effectivepot.dat')
#endif
   ! Find turning points
   DO i = 1, this%vibrpot%n
      x(i)=this%rovibrpot%rpot%x(i)
      f(i)=this%rovibrpot%rpot%f(i)-this%evirot%getvalue()
   END DO
   CALL fdummy%INITIALIZE_DIRECT(x,f,0.D0)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Initialize dummy potential to find turning points")
   CALL VERBOSE_WRITE(routinename,'Dummy potential is Vdummy(r)= Veff(r)-Erovibr')
   CALL fdummy%PLOT(1000,"dummypot.dat")
   CALL VERBOSE_WRITE(routinename,'Dummy potential stored at dummypot.dat')
#endif
   CALL fdummy%SET_ROOTS()
   turn1=fdummy%root(1)
   turn2=fdummy%root(2)
   vturn1=this%rovibrpot%getpot(turn1)
   vturn2=this%rovibrpot%getpot(turn2)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Turning points are: ",[turn1,turn2])
   CALL VERBOSE_WRITE(routinename,"Effective potential at turning points: ",[vturn1,vturn2])
#endif
   ! Start to integrate equations of motion
   t=0.D0
   dt=this%delta_t%getvalue()
   z(1)=turn1
   z(2)=0.D0 ! initial momentum (at a turning point should be zero)
   ! Allocate for first time
   ALLOCATE(r_t(1))
   ALLOCATE(pr_t(1))
   ALLOCATE(time(1))
   r_t(1)=z(1)
   pr_t(1)=z(2)
   time(1)=t
   ! Initial  control values
   Eint=this%evirot%getvalue()
   cycles=0
   reached_turn2=.FALSE.
   onecycle_completed=.FALSE.
   n_req_passedby=0
   DO
      switch=.FALSE.
      cycles=cycles+1
#ifdef DEBUG
      CALL DEBUG_WRITE(routinename,"Cycle: ",cycles)
#endif
      init_t=t
      init_r=z(1)
      init_pr=z(2)
      init_Eint=Eint
      ! Allocating r and p results. Doesn't cycle
      CALL this%TIME_DERIVS(z,dzdt,switch)
      SELECT CASE(switch)
         CASE(.TRUE.)
            SELECT CASE(cycles)
               CASE(1)
                  WRITE(0,*) "SET_PERIOD_INITDIATOMIC ERR: Initial position for integration is not adequate"
                  WRITE(0,*) "DOF'S: ", z
                  CALL EXIT(1)
               CASE DEFAULT
                  WRITE(0,*) "SET_PERIOD_INITDIATOMIC ERR: This error is quite uncommon, guess what is happening by your own."
                  WRITE(0,*) "DOF'S: ",z
                  CALL EXIT(1)
            END SELECT
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      FORALL (i=1:2) s(i) = DABS(z(i))+DABS(dt*dzdt(i)) ! Numerical Recipes.  Smart scaling
#ifdef DEBUG
      CALL DEBUG_WRITE(routinename,"Energy before integration:",init_Eint)
      CALL DEBUG_WRITE(routinename,"Position before integration:",z)
      CALL DEBUG_WRITE(routinename,"Time derivatives before integration:",dzdt)
#endif
      CALL this%BSSTEP(z,dzdt,t,dt,this%eps,s,dt_did,dt_next,switch)
#ifdef DEBUG
      CALL DEBUG_WRITE(routinename,"Position after integration:",z)
      CALL DEBUG_WRITE(routinename,"Time derivatives after integration:",dzdt)
#endif
      SELECT CASE(switch)
         CASE(.TRUE.)
            ! Problem detected in integration, reducing time-step
            ! reboot to previous step values
            dt = dt/2.D0
            t = init_t
            z(1)=init_r
            z(2)=init_pr
            Eint=init_Eint
#ifdef DEBUG
            CALL DEBUG_WRITE(routinename,"Error encountered while integration of eq. of motion")
            CALL DEBUG_WRITE(routinename,"Halving time-step to: ", dt)
#endif
            CYCLE
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      Eint=(z(2)**2.D0)/(2.D0*mu)+this%rovibrpot%getpot(z(1))
#ifdef DEBUG
      CALL DEBUG_WRITE(routinename,"Energy after integration:",Eint)
#endif
      ! Check conservation of energy. CAN cycle
      SELECT CASE (DABS(Eint-init_Eint) > this%eps*init_Eint)
         CASE(.TRUE.)
            ! Problems with energy conservation
            ! reboot to previous step values
            dt = dt/2.D0
            t = init_t
            z(1)=init_r
            z(2)=init_pr
            Eint=init_Eint
#ifdef DEBUG
            CALL DEBUG_WRITE(routinename,"Poor energy conservation. Cycling.")
#endif
            CYCLE
      CASE(.FALSE.)
         ! do nothing
      END SELECT
      ! Check time step. Doesn't cycle
      SELECT CASE(dt_next>this%delta_t%getvalue()) ! there's a maximum time-step defined
         CASE(.TRUE.)
#ifdef DEBUG
            CALL DEBUG_WRITE(routinename,"Predicted time-step too large: ",dt_next)
            CALL DEBUG_WRITE(routinename,"Recovering smaller time-step: ",this%delta_t%getvalue())
#endif
            dt_next=this%delta_t%getvalue()
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      ! Check if turningpoint2 reached. Doesn't cycle
      SELECT CASE(z(1)<init_r .AND. .NOT.reached_turn2)
         CASE(.TRUE.)
            reached_turn2=.TRUE.
            time_turn2=t
         CASE(.FALSE.)
      END SELECT
      ! Check if turnpoint1 reached. At this point one cycle has been made. Doesn't cycle
      SELECT CASE(z(1)> init_r .AND. reached_turn2 .AND. .NOT.onecycle_completed)
         CASE(.TRUE.)
            onecycle_completed=.TRUE.
            time_cycle=t
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      ! Check how many times req has been passed by
      SELECT CASE(onecycle_completed .AND. z(1)>= this%rovibrpot%getreq())
         CASE(.TRUE.)
            n_req_passedby=3
            time_finish=t
#ifdef DEBUG
            CALL DEBUG_WRITE(routinename,"Integration stopped at time: ",time_finish)
#endif
            EXIT
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      dt=dt_next
      nold=size(r_t)
      ALLOCATE(aux(nold))
      ! copy and expand r_t
      aux=r_t
      DEALLOCATE(r_t)
      ALLOCATE(r_t(nold+1))
      r_t(1:nold)=aux
      r_t(nold+1)=z(1)
      ! copy and expand pr_t
      aux=pr_t
      DEALLOCATE(pr_t)
      ALLOCATE(pr_t(nold+1))
      pr_t(1:nold)=aux
      pr_t(nold+1)=z(2)
      ! copy and expand time
      aux=time
      DEALLOCATE(time)
      ALLOCATE(time(nold+1))
      time(1:nold)=aux
      time(nold+1)=t
      DEALLOCATE(aux)
   END DO
   CALL this%pr_t%READ(time,pr_t)
   CALL this%pr_t%INTERPOL(0.D0,0,0.D0,0)
   CALL this%r_t%READ(time,r_t)
   CALL this%r_t%INTERPOL(0.D0,0,0.D0,0)
   CALL this%pr_t%SET_MINIMUM()
   CALL this%pr_t%SET_XROOT()
   this%period=maxval(this%pr_t%xroot)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'Radial motion performed inside Veff(r) potential')
   CALL VERBOSE_WRITE(routinename,"Minimums found in pr(t)",this%pr_t%xmin)
   CALL VERBOSE_WRITE(routinename,"Roots found in pr(t)",this%pr_t%xroot)
   CALL VERBOSE_WRITE(routinename,"Time specifications:",[time_turn2,time_cycle,time_finish])
   CALL VERBOSE_WRITE(routinename,"Period of vibrational motion :",this%period)
   CALL this%pr_t%PLOT(1000,"pr.dat")
   CALL this%r_t%PLOT(1000,"r.dat")
   CALL VERBOSE_WRITE(routinename,'Radial momentum as a function of time printed inside pr.dat')
   CALL VERBOSE_WRITE(routinename,'Radial generalized position as a function of time printed inside r.dat')
#endif
   RETURN
END SUBROUTINE SET_PERIOD_INITSURFACEGROW
END MODULE INITSURFACEGROW_MOD
