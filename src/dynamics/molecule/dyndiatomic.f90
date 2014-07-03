!##########################################################
! MODULE: DYNDIATOMIC_MOD
!> @brief
!! Provides tools to run dynamics on atoms
!> @todo
!! - Generalize the use of different PES, not only CRP3D
!##########################################################
MODULE DYNDIATOMIC_MOD
   USE DYNAMICS_MOD
   USE INITDIATOMIC_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
IMPLICIT NONE
!//////////////////////////////////////////////////////////
! TYPE: DYNDIATOMIC
!> @brief
!! Stores information for a dynamics job with atoms
!
!> @param extrapol - Method to extrapol stepsize in differential equations (Polinomi,Rational)
!> @param scaling - Method to scale differential equations (Equal,Smart)
!> @param eps - Precission of equations of motion. This is a factor.
!> @param zstop - Trajectories are stopped when Z is in between zstop and zstop+dzstop
!> @param dzstop - Tolerance value zstop
!> @param zscatt - Z value to consider a trajectory to be scattered
!> @param zads - Z for adsorbed atoms
!> @param zabs - Z for absorbed atoms
!----------------------------------------------------------
TYPE,EXTENDS(Dynamics) :: DYNDIATOMIC
   CHARACTER(LEN=8) :: extrapol
   CHARACTER(LEN=10) :: scaling
   REAL(KIND=8) :: eps
   TYPE(Length) :: zstop, dzstop
   TYPE(Length) :: zscatt,zads,zabs
   CONTAINS
      ! Initialization block
      PROCEDURE,PUBLIC :: INITIALIZE => INITIALIZE_DYNDIATOMIC
      ! Tools block
      PROCEDURE,PUBLIC :: RUN => RUN_DYNDIATOMIC
      ! Private block
      PROCEDURE,PRIVATE :: DO_DYNAMICS => DO_DYNAMICS_DYNDIATOMIC
      PROCEDURE,PRIVATE :: TIME_DERIVS => TIME_DERIVS_DYNDIATOMIC
      PROCEDURE,PRIVATE :: MMID => MMID_DYNDIATOMIC
      PROCEDURE,PRIVATE :: BSSTEP => BSSTEP_DYNDIATOMIC
      PROCEDURE,PRIVATE :: POLINOM_EXTRAPOL => PZEXTR_DYNDIATOMIC
      PROCEDURE,PRIVATE :: RATIONAL_EXTRAPOL => RZEXTR_DYNDIATOMIC
END TYPE DYNDIATOMIC
!//////////////////////////////////////////////////////////
CONTAINS
!#####################################################################
!# SUBROUTINE: INITIALIZE_DYNDIATOMIC ####################################
!#####################################################################
!> @brief 
!! Specific implementation of READ subroutine for atomic dynamics
!---------------------------------------------------------------------
SUBROUTINE INITIALIZE_DYNDIATOMIC(this,filename)
	IMPLICIT NONE
	! I/O variables
	CLASS(Dyndiatomic),INTENT(OUT) :: this
	CHARACTER(LEN=*),INTENT(IN) :: filename
   ! IMPORTANT: unit used to read info
   INTEGER(KIND=4),PARAMETER :: runit=108
	! Local variables
   CHARACTER(LEN=20) :: filenameinicond,filenamepes
	CHARACTER(LEN=18), PARAMETER :: routinename = "READ_DYNDIATOMIC: "
	CHARACTER(LEN=6) :: follow
   CHARACTER(LEN=10) :: units
   REAL(KIND=8) :: aux
	INTEGER :: i ! counters
	! YIPEE KI YAY -----------------------------
	this%filename = filename
	OPEN (runit,FILE=this%filename, STATUS="old")
   READ(runit,*) ! dummy line
	READ(runit,*) this%alias
	READ(runit,*) this%kind
	IF(this%kind.EQ."Diato") THEN
      READ(runit,*) filenamepes
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"File for PES: ",filenamepes)
#endif
      SELECT CASE(filenamepes)
         CASE("INcrp6d.inp")
            ALLOCATE(CRP6D::this%thispes)
            CALL this%thispes%INITIALIZE(filenamepes)
         CASE DEFAULT
            WRITE(0,*) "READ_DYNDIATOMIC ERR: PES file not implemented"
            WRITE(0,*) "Available files: INcrp6d.inp"
            CALL EXIT(1)
      END SELECT
      READ(runit,*) filenameinicond
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"File for Initial conditions: ",filenameinicond)
#endif
      SELECT CASE(filenameinicond)
         CASE("INinicond6d.inp")
            ALLOCATE(Initdiatomic::this%thisinicond)
            CALL this%thisinicond%INITIALIZE(filenameinicond)
            CALL this%thisinicond%GENERATE_TRAJS(this%thispes)
         CASE DEFAULT
            WRITE(0,*) "READ_DYNDIATOMIC ERR: Initial conditions file not implemented"
            WRITE(0,*) "Available files: INinicond6d.inp"
            CALL EXIT(1)
      END SELECT
		READ(runit,*) this%eps
#ifdef DEBUG
		CALL VERBOSE_WRITE(routinename,"Read: precision of integration (dimensionless factor): ",this%eps)
#endif
		READ(runit,*) this%scaling
		IF((this%scaling.NE."Equal").AND.(this%scaling.NE."Smart")) THEN
			WRITE(0,*) "INITIALIZE_DYNDIATOMIC ERR: Wrong error scaling keyword: "
			WRITE(0,*) "Only available: Equal and Smart"
			WRITE(0,*) "You have written: ", this%scaling
			CALL EXIT(1)
		END IF
		READ(runit,*) this%extrapol
		IF ((this%extrapol.NE."Polinomi").AND.(this%extrapol.NE."Rational")) THEN
			WRITE(0,*) "INITIALIZE_DYNDIATOMIC ERR: Wrong extrapolation keyword: "
			WRITE(0,*) "Only available: Polinomi and Rational"
			WRITE(0,*) "You have written: ", this%extrapol
			CALL EXIT(1)
		END IF
		
      READ(runit,*) aux,units
		CALL this%delta_t%READ(aux,units)
		CALL this%delta_t%TO_STD()
		
      READ(runit,*) aux,units
		CALL this%max_t%READ(aux,units)
		CALL this%max_t%TO_STD()
		
      READ(runit,*) aux,units
		CALL this%zstop%READ(aux,units)
		CALL this%zstop%TO_STD()
		
      READ(runit,*) aux,units
		CALL this%dzstop%READ(aux,units)
		CALL this%dzstop%TO_STD()
		
      READ(runit,*) aux,units
		CALL this%zscatt%READ(aux,units)
		CALL this%zscatt%TO_STD()
		
      READ(runit,*) aux,units
		CALL this%zads%READ(aux,units)
		CALL this%zads%TO_STD()
		
      READ(runit,*) aux,units
		CALL this%zabs%READ(aux,units)
		CALL this%zabs%TO_STD()
#ifdef DEBUG
		CALL VERBOSE_WRITE(routinename,"Initial time step in au: ",this%delta_t%getvalue())
		CALL VERBOSE_WRITE(routinename,"Maximum time in au: ",this%max_t%getvalue())
		CALL VERBOSE_WRITE(routinename,"ZSTOP in au: ",this%zstop%getvalue())
		CALL VERBOSE_WRITE(routinename,"DZSTOP in au: ",this%dzstop%getvalue())
		CALL VERBOSE_WRITE(routinename,"Z scattering benchmark in au: ",this%zscatt%getvalue())
		CALL VERBOSE_WRITE(routinename,"Z adsorption benchmark in au: ",this%zads%getvalue())
		CALL VERBOSE_WRITE(routinename,"Z absortion benchmark in au: ",this%zabs%getvalue())
#endif
		READ(runit,*) follow, this%nfollow
		IF (follow.NE."FOLLOW") THEN
			WRITE(0,*) "INITIALIZE_DYNDIATOMIC ERR: wrong FOLLOW keyword"
			CALL EXIT(1)
		ELSE IF (this%nfollow.NE.0) THEN
#ifdef DEBUG
			CALL VERBOSE_WRITE(routinename,"There are trajectories to follow")
#endif
			ALLOCATE(this%followtraj(this%nfollow))
			DO i=1, this%nfollow
				READ(runit,*) this%followtraj(i)
#ifdef DEBUG
            CALL VERBOSE_WRITE(routinename,"Follow trajectories: ",this%followtraj(i))
#endif
			END DO
		ELSE
#ifdef DEBUG
			CALL VERBOSE_WRITE(routinename,"No trajectories to follow")
#endif
		END IF
	END IF
	CLOSE(runit)
	RETURN
END SUBROUTINE INITIALIZE_DYNDIATOMIC
!###############################################################
!# SUBROUTINE: RUN_DYNDIATOMIC #####################################
!###############################################################
!> @brief
!! Launch trajectories as said in @b inicondat variable
!---------------------------------------------------------------
SUBROUTINE RUN_DYNDIATOMIC(this)
   IMPLICIT NONE
   ! I/O variables 
   CLASS(DYNDIATOMIC),INTENT(INOUT) :: this
   ! Local variables 
   INTEGER :: i ! counters
   CHARACTER(LEN=17),PARAMETER :: routinename = "RUN_DYNDIATOMIC: "
   ! HEY HO! LET'S GO !!! ------
   DO i=this%thisinicond%nstart, this%thisinicond%ntraj
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Start trajectory: ",i)
#endif
      CALL this%DO_DYNAMICS(i)
   END DO
   RETURN
END SUBROUTINE RUN_DYNDIATOMIC
!##############################################################
!# SUBROUTINE : DO_DYNAMICS_DYNDIATOMIC ###########################
!##############################################################
!> @brief
!! Integrates equation of motion following specifications stored
!! in @b this
!--------------------------------------------------------------
SUBROUTINE DO_DYNAMICS_DYNDIATOMIC(this,idtraj)
   IMPLICIT NONE
   ! I/O variables
   CLASS(DYNDIATOMIC),TARGET,INTENT(INOUT) :: this
   INTEGER(KIND=4),INTENT(IN) :: idtraj
   ! IMPORTANT: units used to write
   INTEGER(KIND=4),PARAMETER :: wunit1=799,wunit2=789,wunit3=724
   ! Local variables
   INTEGER :: i, cycles ! counters
   REAL(KIND=8) :: t,dt,E,init_t,init_E,v,dt_did,dt_next,zmin, angle
   REAL(KIND=8) :: Eint,Ecm
   REAL(KIND=8),DIMENSION(6) ::  dummy
   REAL(KIND=8),DIMENSION(6) :: atomiccoord
   REAL(KIND=8),DIMENSION(12) :: molec_dofs, s, dfdt
   REAL(KIND=8) :: ma,mb
   LOGICAL :: maxtime_reached
   LOGICAL :: switch, file_exists, in_list
   CHARACTER(LEN=27) :: filename_follow
   CHARACTER(LEN=21),PARAMETER :: routinename = "DO_DYNAMICS_DYNDIATOMIC: "
   CHARACTER(LEN=9),PARAMETER :: format_string = '(I10.10)'
   CHARACTER(LEN=10) :: x1
   INTEGER :: control
   CLASS(Dynobject),POINTER:: molecule
   REAL(KIND=8) :: masa
   REAL(KIND=8) :: mu
   ! HEY HO!, LET'S GO!!! -------------------------
   molecule => this%thisinicond%trajs(idtraj)
   ma=this%thispes%atomdat(1)%getmass()
   mb=this%thispes%atomdat(2)%getmass()
   masa = ma+mb
   mu = ma*mb/masa
   INQUIRE(FILE="OUTdynamics6d.MOLEC.out",EXIST=file_exists)
   SELECT CASE(file_exists)
      CASE(.TRUE.)
         OPEN(wunit1, FILE="OUTdynamics6d.MOLEC.out",STATUS="old",ACCESS="append")
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Previous output file found: OUTdynamics6d.MOLEC.out")
         CALL VERBOSE_WRITE(routinename,"Tajectories will be added to this file")
#endif
      CASE(.FALSE.)
         OPEN(wunit1,FILE="OUTdynamics6d.MOLEC.out",STATUS="new")
         WRITE(wunit1,*) "# DYNAMICS RESULTS MOLECULAR COORDINATES -------------------------------------"
         WRITE(wunit1,*) "# Format: id, status, ireb, ixyboun, Etot, Eint, &
            &t, X,Y,Z,R,THETA,PHI (a.u.), Px,Py,Pz,Pr,Ptheta,Pphi (a.u.)"
         WRITE(wunit1,*) "# ----------------------------------------------------------------------------"
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"New output file created: OUTdynamics6d.MOLEC.out")
         CALL VERBOSE_WRITE(routinename,"Header printed to that file")
         CALL VERBOSE_WRITE(routinename,"Tajectories will be added to this file")
#endif
   END SELECT
   INQUIRE(FILE="OUTturning6d.MOLEC.out",EXIST=file_exists)
   SELECT CASE(file_exists)
      CASE(.TRUE.)
         OPEN(wunit2, FILE="OUTturning6d.MOLEC.out",STATUS="old",ACCESS="append")
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Previous output file found: OUTturning6d.MOLEC.out")
         CALL VERBOSE_WRITE(routinename,"Tajectories will be added to this file")
#endif
      CASE(.FALSE.)
         OPEN(wunit2,FILE="OUTturning6d.MOLEC.out",STATUS="new")
         WRITE(wunit2,*) "# TURNING POINTS MOLECULAR COORDINATES ---------------------"
         WRITE(wunit2,*) "# Description: positions of scattered atoms at their lowest "
         WRITE(wunit2,*) "#              Z value reached during the dynamics.         "
         WRITE(wunit2,*) "#              XY values, projected into IWS cell           "
         WRITE(wunit2,*) "# Format: id, X,Y,Z,R,THETA,PHI (a.u.)                      "
         WRITE(wunit2,*) "# ----------------------------------------------------------"
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"New output file created: OUTturning6d.MOLEC.out")
         CALL VERBOSE_WRITE(routinename,"Header printed to that file")
         CALL VERBOSE_WRITE(routinename,"Tajectories will be added to this file")
#endif
   END SELECT
   in_list = .FALSE.
   SELECT CASE(this%nfollow)
      CASE(0)
         ! do nothing
      CASE DEFAULT
         DO i=1,this%nfollow
            SELECT CASE(this%followtraj(i)==idtraj)
               CASE(.TRUE.)
                  in_list=.TRUE.
                  EXIT
               CASE(.FALSE.)
                  ! do nothing
            END SELECT
         END DO
         SELECT CASE(in_list)
            CASE(.TRUE.)
               WRITE(x1,format_string) idtraj
               filename_follow = 'OUTtraj'//TRIM(x1)//'.MOLEC.out'
               OPEN(wunit3,FILE=filename_follow,STATUS="replace")
               WRITE(wunit3,*) "# TIME EVOLUTION FOR A TRAJECTORY --------------------------------"
               WRITE(wunit3,*) "# Format: t,dt,Etot,Enorm,Eint,Pot,X,Y,Z,R,THETA,PHI,&
                  &Px,Py,Pz,Pr,Ptheta,Pphi,Xa,Ya,Za,Xb,Yb,Zb"
               WRITE(wunit3,*) "# First and last position are not printed here                    "
               WRITE(wunit3,*) "# You can find them at OUTdynamics6d.out and initial conditions output file"
               WRITE(wunit3,*) "# ----------------------------------------------------------------"
#ifdef DEBUG
               CALL VERBOSE_WRITE(routinename,"Trajectory followed: ",idtraj)
               CALL VERBOSE_WRITE(routinename,"File created: ",filename_follow)
#endif
            CASE(.FALSE.)
               ! do nothing
         END SELECT
   END SELECT
   cycles=0 
   zmin = molecule%init_r(3)
   t=0.D0
   molecule%ireb=0
   molecule%ixyboun=0
   dt=this%delta_t%getvalue()
   ! Iterations
   DO
      switch = .FALSE.
      cycles = cycles+1
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Cycle: ",cycles)
#endif
      ! We cannot go beyond maximum time defined
      SELECT CASE(t+dt > this%max_t%getvalue())
         CASE(.TRUE.)
            dt=this%max_t%getvalue()-t
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      init_t = t
      ! Storing atomic DOF's in only one array
      molec_dofs(1:6)=molecule%r(1:6)
      molec_dofs(7:12)=molecule%p(1:6)
      SELECT CASE(dsin(molec_dofs(5))<=0.D0) ! when sinus(theta) <= 0
         CASE(.TRUE.)
            molec_dofs(6)=0.D0 ! phi cannot be defined
            molec_dofs(12)=0.D0 ! neither its conjugate momentum
         CASE(.FALSE.)
            ! Usual definition of coordinates
      END SELECT
      ! Initial values for the derivatives
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Molecular DOFS : ",molec_dofs)
#endif
      CALL this%TIME_DERIVS(molec_dofs,dfdt,switch)
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Time derivatives: ",dfdt)
#endif
      SELECT CASE(switch)
         CASE(.TRUE.)
            SELECT CASE(cycles)
               CASE(1)
                  WRITE(0,*) "DO_DYNAMICS_DYNDIATOMIC ERR: Initial position is not adequate"
                  WRITE(0,*) "molecule id: ", idtraj
                  WRITE(0,*) "Atom DOF'S: ", molec_dofs(:)
                  CALL EXIT(1)
               CASE DEFAULT
                  WRITE(0,*) "DO_DYNAMICS_DYNDIATOMIC ERR: This error is quite uncommon, guess what is happening by your own."
                  WRITE(0,*) "Atom id: ", idtraj
                  WRITE(0,*) "Atom DOF'S: ", molec_dofs(:)
                  CALL EXIT(1)
            END SELECT
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      ! Construct scaling array
      SELECT CASE(this%scaling)
         CASE("Smart")
            FORALL (i=1:12) s(i) = DABS(molec_dofs(i))+DABS(dt*dfdt(i)) ! Numerical Recipes.  
         CASE("Equal")
            s(1:12)=1.D0
         CASE DEFAULT
            WRITE(0,*) "DO_DYNAMICS ERR: incorrect scaling Keyword"
            CALL EXIT(1)
      END SELECT
      ! Energy before time-integration (a.u.)
      init_E = molecule%E
      ! Integrate and choose the correct time-step
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Energy before integration:",molecule%E)
#endif
      ! Call integrator
      CALL this%BSSTEP(molec_dofs,dfdt,t,dt,this%eps,s,dt_did,dt_next,switch)
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Atomic DOFs after integration: ",molec_dofs)
#endif
      SELECT CASE(switch)
         CASE(.TRUE.)
            ! Problem detected in integration, reducing time-step
            ! reboot to previous step values 
            dt = dt/2.D0 
            t = init_t
#ifdef DEBUG
            CALL VERBOSE_WRITE(routinename,"Error encountered while integration of eq. of motion")
            CALL VERBOSE_WRITE(routinename,"Halving time-step to: ", dt)
#endif
            CYCLE
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      ! At this point, molec_dofs contains the new values for position and momenta
      ! Let's obtain the potential value for this configuration
      CALL this%thispes%GET_V_AND_DERIVS(molec_dofs(1:6),v,dummy)
      Ecm = (molec_dofs(7)**2.D0+molec_dofs(8)**2.D0+molec_dofs(9)**2.D0)/(2.D0*masa)
      Eint= molec_dofs(10)**2.D0+(molec_dofs(11)/molec_dofs(4))**2.D0+(molec_dofs(12)/(molec_dofs(4)*dsin(molec_dofs(5))))**2.D0
      Eint= Eint/(2.D0*mu)
      E=Ecm+Eint+v
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Energy after integration:",E)
#endif
      SELECT CASE (DABS(E-molecule%E) > 100.D0*this%eps*molecule%E)
         CASE(.TRUE.)
            ! Problems with energy conservation
            ! reboot to previous step values
            dt = dt/2.D0
            t = init_t
#ifdef DEBUG
            CALL VERBOSE_WRITE(routinename,"Poor energy conservation. Cycling.")
#endif
            CYCLE
      CASE(.FALSE.)
         ! do nothing
      END SELECT
      ! Check  bouncing points in Z direction (Z- Turning point)
      SELECT CASE((molecule%p(3) < 0.D0).AND.(molec_dofs(9) > 0.D0))
         CASE(.TRUE.)
            molecule%ireb = molecule%ireb +1
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      ! Check bouncing points in XY (sign of Px or Py changes respect to previous integration step)
      SELECT CASE ((DSIGN(molec_dofs(7),molec_dofs(7)).NE.DSIGN(molec_dofs(7),molecule%p(1))).OR. &
                  (DSIGN(molec_dofs(8),molec_dofs(8)).NE.DSIGN(molec_dofs(8),molecule%p(2))))
         CASE(.TRUE.)
            molecule%ixyboun = molecule%ixyboun+1
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      ! Store last minimum Z value reached. This will be the turning-point
      SELECT CASE (molec_dofs(3) <= zmin)
         CASE(.TRUE.)
            zmin = molec_dofs(3)
            molecule%turning_point(1:6) = molec_dofs(1:6)
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      ! EXIT CHANNELS, AKA SWITCH SECTION -------------------------------------------
      SELECT CASE(t > this%max_t%getvalue())
         CASE(.TRUE.)
            maxtime_reached=.TRUE.
         CASE(.FALSE.)
            maxtime_reached=.FALSE.
      END SELECT
      SELECT CASE ((molec_dofs(3) <= this%zstop%getvalue()+this%dzstop%getvalue()).AND. &
                  (molec_dofs(3) > this%zstop%getvalue()-this%dzstop%getvalue()))
         CASE(.TRUE.)
            molecule%stat = "Stopped"
            molecule%r(1:6) = molec_dofs(1:6)
            molecule%p(1:6) = molec_dofs(7:12)
            molecule%E = E
            molecule%Eint = Eint
            molecule%Ecm = Ecm
            EXIT
        CASE(.FALSE.)
           ! do nothing, next switch
      END SELECT
      SELECT CASE ((molecule%r(3) >= this%zscatt%getvalue()).AND.(molecule%p(3) > 0.D0))
         CASE(.TRUE.)
            molecule%stat="Scattered"                                            
            molecule%r(1:6)=molec_dofs(1:6)
            molecule%p(1:6)=molec_dofs(7:12)
            molecule%E=E
            molecule%Eint=Eint
            molecule%Ecm=Ecm
            molecule%turning_point(1:2) = this%thispes%surf%project_iwscell(molecule%turning_point(1:2))
            WRITE(wunit2,'(I7,1X,6(F10.5,1X))') idtraj, molecule%turning_point(:)
            EXIT
         CASE(.FALSE.)
            ! do nothing next switch
      END SELECT
      SELECT CASE ((molec_dofs(3) <= this%zabs%getvalue()).AND.(molec_dofs(9) < 0.D0))
         CASE(.TRUE.)
            molecule%stat = "Absorbed"
            molecule%r(1:6) = molec_dofs(1:6)
            molecule%p(1:6) = molec_dofs(7:12)
            molecule%E = E
            molecule%Eint =Eint
            molecule%Ecm =Ecm
            EXIT
         CASE(.FALSE.)
            ! do nothing, next switch
      END SELECT
      SELECT CASE ((v < 0.D0).AND.(molec_dofs(3) <= this%zads%getvalue()).AND.maxtime_reached)
         CASE(.TRUE.)
            molecule%stat = "Adsorbed"
            molecule%r(1:6) = molec_dofs(1:6)
            molecule%p(1:6) = molec_dofs(7:12)
            molecule%E = E
            molecule%Eint = Eint
            molecule%Ecm = Ecm
            EXIT
         CASE(.FALSE.)
            ! do nothing, next switch
      END SELECT
      SELECT CASE (molec_dofs(3) <= this%zads%getvalue() .AND. maxtime_reached)
         CASE(.TRUE.)
            molecule%stat = "Trapped"
            molecule%r(1:6) = molec_dofs(1:6)
            molecule%p(1:6) = molec_dofs(7:12)
            molecule%E = E
            molecule%Eint = Eint
            molecule%Ecm = Ecm
            EXIT
         CASE(.FALSE.)
            ! do nothing, next switch
      END SELECT
      SELECT CASE (maxtime_reached)
         CASE(.TRUE.)
            molecule%stat = "Time-out"
            molecule%r(1:6) = molec_dofs(1:6)
            molecule%p(1:6) = molec_dofs(7:12)
            molecule%E = E
            molecule%Eint = Eint
            molecule%Ecm = Ecm
            EXIT
         CASE(.FALSE.)
            !do nothing next switch
      END SELECT
      SELECT CASE((cycles > 1000).AND.(dt < 1.D-9))
         CASE(.TRUE.)
            molecule%stat = "Patologic"
            molecule%r(1:6) = molec_dofs(1:6)
            molecule%p(1:6) = molec_dofs(7:12)
            molecule%E = E
            molecule%Eint = Eint
            molecule%Ecm = Ecm
            EXIT
         CASE(.FALSE.)
            ! do nothing, next switch
      END SELECT
      SELECT CASE(.NOT.maxtime_reached)
         CASE(.TRUE.)
            molecule%r(1:6) = molec_dofs(1:6)
            molecule%p(1:6) = molec_dofs(7:12)
            molecule%E = E
            molecule%Eint = Eint
            molecule%Ecm = Ecm
            dt = dt_next
            SELECT CASE((this%nfollow.NE.0).AND.(in_list))
               CASE(.TRUE.)
                  CALL FROM_MOLECULAR_TO_ATOMIC(ma,mb,molecule%r,atomiccoord)
                  WRITE(wunit3,*) t,dt_did,molecule%E,(molecule%p(3)**2.D0)/(2.D0*masa),&
                     molecule%Eint,v,molecule%r(:),molecule%p(:),atomiccoord
               CASE(.FALSE.)
                  ! do nothing
            END SELECT
            CYCLE
         CASE(.FALSE.)
            WRITE(0,*) "D=_DYNAMICS_DYNDIATOMIC ERR: Strange trajectory conditions. Switchs cannot classify this traj"
            CALL EXIT(1)
      END SELECT
   END DO
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Writing in dynamics.out")
#endif
   WRITE(wunit1,'(I7,1X,A10,1X,I5,1X,I5,1X,15(F15.5,1X))') idtraj,molecule%stat,molecule%ireb,&
      molecule%ixyboun,molecule%E,molecule%Eint,t,molecule%r,molecule%p
   CLOSE(wunit1)
   CLOSE(wunit2)
   SELECT CASE(this%nfollow/=0)
      CASE(.TRUE.)
         CLOSE(wunit3)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE DO_DYNAMICS_DYNDIATOMIC
!###############################################################
!# SUBROUTINE : TIME_DERIVS_DYNDIATOMIC ########################
!###############################################################
!> @brief
!! Gives dzdt at z and t values from Hamilton equations of motion
!
!> @param[in] this - Provides some data
!> @param[in] z - array of positions and momenta. z(1:3) -> positions, z(4:6) -> momenta 
!> @param[out] dzdt - time derivatives of position and momenta
!> @param[out] fin - controls errors
!--------------------------------------------------------------
SUBROUTINE TIME_DERIVS_DYNDIATOMIC(this,z,dzdt,fin)
   IMPLICIT NONE
   ! I/O variables
   CLASS(DYNDIATOMIC),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(12),INTENT(IN) :: z
   REAL(KIND=8),DIMENSION(12),INTENT(OUT) :: dzdt
   LOGICAL, INTENT(OUT) :: fin
   ! Local variables
   INTEGER :: i ! counters
   REAL(KIND=8) :: x,y,zeta,r,theta,phi
   REAL(KIND=8) :: px,py,pzeta,pr,ptheta,pphi
   REAL(KIND=8) :: mass
   REAL(KIND=8) :: mu
   REAL(KIND=8) :: v ! dummy variable
   CHARACTER(LEN=21),PARAMETER :: routinename = "TIME_DERIVS_DYNDIATOMIC: "
   ! ROCK THE CASBAH ! ---------------------
   SELECT CASE(this%thispes%is_allowed(z(1:6)))
      CASE(.FALSE.)
         fin = .TRUE.
      CASE(.TRUE.)
         fin=.FALSE.
         mass=this%thispes%atomdat(1)%getmass()+this%thispes%atomdat(2)%getmass()
         mu=this%thispes%atomdat(1)%getmass()*this%thispes%atomdat(2)%getmass()
         mu=mu/(this%thispes%atomdat(1)%getmass()+this%thispes%atomdat(2)%getmass())
         SELECT CASE(dsin(z(5))<=0.D0)
            CASE(.TRUE.)
               ! Set time derivatives of position
               dzdt(1)=z(7)/mass
               dzdt(2)=z(8)/mass
               dzdt(3)=z(9)/mass
               dzdt(4)=z(10)/mu
               dzdt(5)=z(11)/(mu*(z(4)**2.D0))
               dzdt(6)=0.D0
               ! Set time derivatives of momenta
               CALL this%thispes%GET_V_AND_DERIVS(z(1:6),v,dzdt(7:12))
               dzdt(7)=-dzdt(7)
               dzdt(8)=-dzdt(8)
               dzdt(9)=-dzdt(9)         
               dzdt(10)=-dzdt(10)+(z(11)**2.D0)/(mu*(z(4)**3.D0))
               dzdt(11)=-dzdt(11)
               dzdt(12)=0.D0 ! this property stands in the potential but, anyway
            CASE(.FALSE.) ! usual scheme
               ! Set time derivatives of position
               dzdt(1)=z(7)/mass
               dzdt(2)=z(8)/mass
               dzdt(3)=z(9)/mass
               dzdt(4)=z(10)/mu
               dzdt(5)=z(11)/(mu*(z(4)**2.D0))
               dzdt(6)=z(12)/(mu*(z(4)*dsin(z(5))**2.D0))
               ! Set time derivatives of momenta
               CALL this%thispes%GET_V_AND_DERIVS(z(1:6),v,dzdt(7:12))
               dzdt(7)=-dzdt(7)
               dzdt(8)=-dzdt(8)
               dzdt(9)=-dzdt(9)         
               dzdt(10)=-dzdt(10)+(z(11)**2.D0+(z(12)/dsin(z(5)))**2.D0)/(mu*(z(4)**3.D0))
               dzdt(11)=-dzdt(11)+((z(12)/z(4))**2.D0)*dcos(z(5))/(mu*dsin(z(5))**3.D0)
               dzdt(12)=-dzdt(12)
         END SELECT
   END SELECT
   RETURN
END SUBROUTINE TIME_DERIVS_DYNDIATOMIC
!#########################################################################################
!# SUBROUTINE: MMID_ATOM #################################################################
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
SUBROUTINE MMID_DYNDIATOMIC(this,y,dydx,xs,htot,nstep,yout,switch)
   IMPLICIT NONE
   ! I/O variables
   CLASS(DYNDIATOMIC),INTENT(IN) :: this
   INTEGER,INTENT(IN) :: nstep
   REAL(KIND=8),DIMENSION(12), INTENT(IN) :: y,dydx
   REAL(KIND=8),DIMENSION(12), INTENT(OUT) :: yout
   REAL(KIND=8),INTENT(IN) :: xs,htot
   LOGICAL,INTENT(INOUT) :: switch
   ! Local variables
   REAL(KIND=8),DIMENSION(12) :: ym,yn
   INTEGER :: i,n ! counters
   INTEGER,PARAMETER :: nvar = 12
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
END SUBROUTINE MMID_DYNDIATOMIC
!##################################################################################################
!# SUBROUTINE: RZEXTR_DYNDIATOMIC #####################################################################
!################################################################################################## 
!> @brief 
!! - A part of the Burlich-Stoer algorithm. Uses diagonal rational function extrapolation. 
!
!>@details
!! Taken from Numerical recipes in Fortran 77
!> @see pzextr
!--------------------------------------------------------------------------------------------------
SUBROUTINE RZEXTR_DYNDIATOMIC(this,iest,xest,yest,yz,dy,nv)
	IMPLICIT NONE
	! I/O variables
   CLASS(DYNDIATOMIC),INTENT(IN):: this
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
END SUBROUTINE RZEXTR_DYNDIATOMIC
!############################################################################################
!# SUBROUTINE : PZEXTR_DYNDIATOMIC ##############################################################
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
SUBROUTINE PZEXTR_DYNDIATOMIC(this,iest,xest,yest,yz,dy,nv)
   IMPLICIT NONE
   ! I/O variables
   CLASS(DYNDIATOMIC),INTENT(IN):: this
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
END SUBROUTINE PZEXTR_DYNDIATOMIC
!###############################################################################################
!# SUBROUTINE : BSSTEP_ATOM ####################################################################
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
SUBROUTINE BSSTEP_DYNDIATOMIC(this,y,dydx,x,htry,eps,yscal,hdid,hnext,switch)
	IMPLICIT NONE
	! I/O variables
	CLASS(DYNDIATOMIC),INTENT(IN) :: this
	REAL(KIND=8), INTENT(IN) :: eps     ! required accuracy
	REAL(KIND=8), INTENT(IN) ::  htry   ! step to try
	REAL(KIND=8), DIMENSION(12) :: yscal ! factors to scale error 
	REAL(KIND=8), DIMENSION(12), INTENT(IN) :: dydx 
	REAL(KIND=8), DIMENSION(12), INTENT(INOUT) :: y ! initial/final values for: X,Y,Z,Px,Py,Pz (in this order)
	REAL(KIND=8), INTENT(INOUT) :: x
	LOGICAL, INTENT(INOUT) :: switch ! .TRUE. if potential could not be calculated
	REAL(KIND=8), INTENT(OUT) :: hdid  ! step actually used
	REAL(KIND=8), INTENT(OUT) :: hnext ! guess of the next step
	! Parameters for this routine
	INTEGER, PARAMETER :: nv = 12
	REAL(KIND=8),PARAMETER :: SAFE1 = 0.25D0
	REAL(KIND=8),PARAMETER :: SAFE2 = 0.7D0
	REAL(KIND=8),PARAMETER :: TINY = 1.D-30 
	REAL(KIND=8),PARAMETER :: SCALMX = 0.1D0
	REAL(KIND=8),PARAMETER :: REDMIN = 0.7D0
	REAL(KIND=8),PARAMETER :: REDMAX = 1.D-5
	INTEGER,PARAMETER :: NMAX = 50
	INTEGER,PARAMETER :: KMAXX = 8
	INTEGER,PARAMETER :: IMAX = KMAXX+1
	CHARACTER(LEN=16), PARAMETER :: routinename = "BSSTEP_DYNDIATOMIC: "
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
END SUBROUTINE BSSTEP_DYNDIATOMIC

END MODULE DYNDIATOMIC_MOD
