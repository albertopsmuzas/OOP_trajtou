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
!/////////////////////////////////////////////////////
! TYPE: INITSURFACEGROW
!> @brief
!! Sets initial conditions for atoms
!----------------------------------------------------
TYPE,EXTENDS(Inicond):: InitSurfaceGrowDiatomic
   CHARACTER(LEN=:),ALLOCATABLE::extrapol
   logical:: fixed_theta
   logical:: is_classic
   REAL(KIND=8):: eps
   LOGICAL:: control_vel,control_posX,control_posY,control_out
   REAL(KIND=8):: impact_x, impact_y
   INTEGER(KIND=4),DIMENSION(3):: init_qn
   type(VacuumPot):: vibrPot
   type(Energy):: E_norm, evirot
   type(Angle):: vz_angle, vpar_angle
   type(Length):: init_z ! initial Z value
   CONTAINS
      ! Initialization block
      procedure,public:: INITIALIZE => INITIALIZE_INITSURFACEGROW
      procedure,public:: generate_trajs => generate_trajs_INITSURFACEGROW
END TYPE InitSurfaceGrowDiatomic
!////////////////////////////////////////////////////
CONTAINS
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
!# SUBROUTINE: INITIALIZE_INITSURFACEGROW ################################################
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
   CHARACTER(LEN=*), PARAMETER :: routinename = "INITIALIZE_INITSURFACEGROW: "
   ! Lua variables
   TYPE(flu_State):: conf
   INTEGER(KIND=4):: inicond_table,trajlist_table,magnitude_table,out_table
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
         WRITE(0,*) "INITIALIZE_INITSURFACEGROW ERR: wrong kind of initial conditions"
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
         WRITE(0,*) "INITIALIZE_INITSURFACEGROW ERR: wrong incidence angle: ",this%vz_angle%getvalue()
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
!# SUBROUTINE: GENERATE_TRAJS_INITSURFACEGROW ##########################
!####################################################################
!> @brief
!! Creates trajectories of atoms. All data stored as atomic units.
!
!> @param[in] this - Initial conditions for atoms to be used
!> @param[in] thisPes - Dummy PES. We do not need it.
!
!> @warning
!! - Auxiliar cartesian coordinates.
!! - For atom/surface systems
!--------------------------------------------------------------------
SUBROUTINE GENERATE_TRAJS_INITSURFACEGROW(this,thispes)
   IMPLICIT NONE
   ! I/O variables
   CLASS(InitSurfaceGrowDiatomic),INTENT(INOUT):: this
   CLASS(PES),INTENT(IN):: thispes
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
   open(unit=runit,file='init.raw.dat',status='old',action='read')
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
      write(*,*) 0.5d0*(p(4)**2.d0)/mu+0.5d0*(p(5)**2.d0+(p(6)/dsin(r(5)))**2.d0)/(mu*r(4)**2.d0)+this%vibrPot%getPot(r(4))
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
END SUBROUTINE GENERATE_TRAJS_INITSURFACEGROW

END MODULE INITSURFACEGROW_MOD
