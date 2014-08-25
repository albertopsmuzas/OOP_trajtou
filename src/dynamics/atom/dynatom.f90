!##########################################################
! MODULE: DYNATOM_MOD
!> @brief
!! Provides tools to run dynamics on atoms
!> @todo
!! - Generalize the use of different PES, not only CRP3D
!##########################################################
MODULE DYNATOM_MOD
   USE DYNAMICS_MOD
   USE INITATOM_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
IMPLICIT NONE
!//////////////////////////////////////////////////////////
! TYPE: Dynatom
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
TYPE,EXTENDS(Dynamics) :: Dynatom
   CHARACTER(LEN=8) :: extrapol
   CHARACTER(LEN=10) :: scaling
   REAL(KIND=8) :: eps
   TYPE(Length) :: zstop, dzstop
   TYPE(Length) :: zscatt,zads,zabs
   INTEGER(KIND=4),PRIVATE :: wusc=800 ! write unit for scattered trajs
   INTEGER(KIND=4),PRIVATE :: wupa=801 ! write unit for pathologic trajs
   INTEGER(KIND=4),PRIVATE :: wuto=802 ! write unit for timed out trajs
   INTEGER(KIND=4),PRIVATE :: wutr=803 ! write unit for trapped trajs
   INTEGER(KIND=4),PRIVATE :: wuad=804 ! write unit for adsorbed trajs
   INTEGER(KIND=4),PRIVATE :: wuab=805 ! write unit for absorbed trajs
   INTEGER(KIND=4),PRIVATE :: wust=806 ! write unit for stopped trajs
   INTEGER(KIND=4),PRIVATE :: wutp=807 ! write unit for turning points
   INTEGER(KIND=4),PRIVATE :: wufo=808 ! write unit for trajectory step by step

   CONTAINS
      ! Initialization block
      PROCEDURE,PUBLIC :: INITIALIZE => INITIALIZE_DYNATOM
      ! Tools block
      PROCEDURE,PUBLIC :: RUN => RUN_DYNATOM
      ! Private block
      PROCEDURE,PRIVATE :: DO_DYNAMICS => DO_DYNAMICS_DYNATOM
      PROCEDURE,PRIVATE :: TIME_DERIVS => TIME_DERIVS_DYNATOM
      PROCEDURE,PRIVATE :: MMID => MMID_DYNATOM
      PROCEDURE,PRIVATE :: BSSTEP => BSSTEP_DYNATOM
      PROCEDURE,PRIVATE :: POLINOM_EXTRAPOL => PZEXTR_DYNATOM
      PROCEDURE,PRIVATE :: RATIONAL_EXTRAPOL => RZEXTR_DYNATOM
END TYPE Dynatom
!//////////////////////////////////////////////////////////
CONTAINS
!#####################################################################
!# SUBROUTINE: INITIALIZE_DYNATOM ####################################
!#####################################################################
!> @brief 
!! Specific implementation of READ subroutine for atomic dynamics
!---------------------------------------------------------------------
SUBROUTINE INITIALIZE_DYNATOM(this,filename)
	IMPLICIT NONE
	! I/O variables
	CLASS(Dynatom),INTENT(OUT) :: this
	CHARACTER(LEN=*),INTENT(IN) :: filename
   ! IMPORTANT: unit used to read info
   INTEGER(KIND=4),PARAMETER :: runit=11
	! Local variables
   CHARACTER(LEN=20) :: filenameinicond,filenamepes
	CHARACTER(LEN=14), PARAMETER :: routinename = "READ_DYNATOM: "
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
	IF(this%kind.EQ."Atoms") THEN
      READ(runit,*) filenamepes
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"File for PES: ",filenamepes)
#endif
      SELECT CASE(filenamepes)
         CASE("INcrp3d.inp")
            ALLOCATE(CRP3D::this%thispes)
            CALL this%thispes%INITIALIZE(filenamepes)
         CASE DEFAULT
            WRITE(0,*) "READ_DYNATOM ERR: PES file not implemented"
            WRITE(0,*) "Available files: INcrp3d.inp"
            CALL EXIT(1)
      END SELECT
      READ(runit,*) filenameinicond
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"File for Initial conditions: ",filenameinicond)
#endif
      SELECT CASE(filenameinicond)
         CASE("INinicond3d.inp")
            ALLOCATE(Initatom::this%thisinicond)
            CALL this%thisinicond%INITIALIZE(filenameinicond)
            CALL this%thisinicond%GENERATE_TRAJS(this%thispes)
         CASE DEFAULT
            WRITE(0,*) "READ_DYNATOM ERR: Initial conditions file not implemented"
            WRITE(0,*) "Available files: INinicond3d.inp"
            CALL EXIT(1)
      END SELECT
		READ(runit,*) this%eps
#ifdef DEBUG
		CALL VERBOSE_WRITE(routinename,"Read: precision of integration (dimensionless factor): ",this%eps)
#endif
		READ(runit,*) this%scaling
		IF((this%scaling.NE."Equal").AND.(this%scaling.NE."Smart")) THEN
			WRITE(0,*) "INITIALIZE_DYNATOM ERR: Wrong error scaling keyword: "
			WRITE(0,*) "Only available: Equal and Smart"
			WRITE(0,*) "You have written: ", this%scaling
			CALL EXIT(1)
		END IF
		READ(runit,*) this%extrapol
		IF ((this%extrapol.NE."Polinomi").AND.(this%extrapol.NE."Rational")) THEN
			WRITE(0,*) "INITIALIZE_DYNATOM ERR: Wrong extrapolation keyword: "
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
			WRITE(0,*) "INITIALIZE_DYNATOM ERR: wrong FOLLOW keyword"
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
END SUBROUTINE INITIALIZE_DYNATOM
!###############################################################
!# SUBROUTINE: RUN_DYNATOM #####################################
!###############################################################
!> @brief
!! Launch trajectories as said in @b inicondat variable
!---------------------------------------------------------------
SUBROUTINE RUN_DYNATOM(this)
   IMPLICIT NONE
   ! I/O variables 
   CLASS(Dynatom),INTENT(INOUT) :: this
   ! Local variables 
   INTEGER :: i ! counters
   CHARACTER(LEN=20),PARAMETER :: routinename = "RUN_DYNAMICS_ATOMS: "
   ! HEY HO! LET'S GO !!! ------
   ! Check files for all traj status
   CALL FILE_TRAJSTATUS_DYNATOM(this%wusc,"OUTDYN3Dscattered.out","SCATTERED TRAJS")
   CALL FILE_TRAJSTATUS_DYNATOM(this%wupa,"OUTDYN3Dpatologic.out","PATOLOGIC TRAJS")
   CALL FILE_TRAJSTATUS_DYNATOM(this%wuto,"OUTDYN3Dtimeout.out","TIME-OUT TRAJS")
   CALL FILE_TRAJSTATUS_DYNATOM(this%wutr,"OUTDYN3Dtrapped.out","TRAPPED TRAJS")
   CALL FILE_TRAJSTATUS_DYNATOM(this%wuad,"OUTDYN3Dadsorbed.out","ADSORBED TRAJS")
   CALL FILE_TRAJSTATUS_DYNATOM(this%wuab,"OUTDYN3Dabsorbed.out","ABSORBED TRAJS")
   CALL FILE_TRAJSTATUS_DYNATOM(this%wust,"OUTDYN3Dstopped.out","STOPPED TRAJS")
   ! Check turning points file
   CALL FILE_TURNING_DYNATOM(this%wutp)
   ! Run trajectories one by one
   DO i=this%thisinicond%nstart,this%thisinicond%ntraj
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Start trajectory: ",i)
#endif
      CALL this%DO_DYNAMICS(i)
   END DO
   ! Closing opened units
   CLOSE(this%wusc)
   CLOSE(this%wupa)
   CLOSE(this%wuto)
   CLOSE(this%wutr)
   CLOSE(this%wuad)
   CLOSE(this%wuab)
   CLOSE(this%wust)
   CLOSE(this%wutp)
   RETURN
END SUBROUTINE RUN_DYNATOM
!##############################################################
!# SUBROUTINE : DO_DYNAMICS_DYNATOM ###########################
!##############################################################
!> @brief
!! Integrates equation of motion following specifications stored
!! in @b this
!--------------------------------------------------------------
SUBROUTINE DO_DYNAMICS_DYNATOM(this,idtraj)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Dynatom),TARGET,INTENT(INOUT) :: this
   INTEGER(KIND=4),INTENT(IN) :: idtraj
   ! Local variables
   INTEGER :: i, cycles ! counters
   REAL(KIND=8) :: t,dt,E,init_t,init_E,v,dt_did,dt_next,zmin,angle
   REAL(KIND=8),DIMENSION(3) :: r0, p0, dummy
   REAL(KIND=8),DIMENSION(6) :: atom_dofs, s, dfdt
   LOGICAL :: maxtime_reached
   LOGICAL :: switch, file_exists, in_list
   CHARACTER(LEN=21),PARAMETER :: routinename = "DO_DYNAMICS_DYNATOM: "
   INTEGER :: control
   CLASS(Dynobject),POINTER:: atomo
   REAL(KIND=8) :: masa
   ! Some Formats
10 FORMAT(I7,1X,A10,1X,I4,1X,I5,1X,8(F15.5,1X)) ! Format to print in status files
11 FORMAT(I7,1X,3(F10.5,1X)) ! Format to print in turning points file
   ! HEY HO!, LET'S GO!!! -------------------------
   atomo => this%thisinicond%trajs(idtraj)
   masa = this%thispes%atomdat(1)%getmass()
   ! Check if there are trajectories fo follow step by step
   in_list = .FALSE.
   SELECT CASE(this%nfollow)
      CASE(0) ! Doh!, there ain't trajs to follow, continue with normal dynamics
         ! do nothing
      CASE DEFAULT ! Ohh some trajs to follow
         DO i=1,this%nfollow
            SELECT CASE(this%followtraj(i)==idtraj)
               CASE(.TRUE.)
                  in_list=.TRUE.
                  EXIT
               CASE(.FALSE.)
                  ! do nothing
            END SELECT
         END DO
         ! Check if this traj is in list
         SELECT CASE(in_list)
            CASE(.TRUE.) ! Follow that traj!
               CALL FILE_FOLLOWTRAJ_DYNATOM(this%wufo,idtraj)
#ifdef DEBUG
               CALL VERBOSE_WRITE(routinename,"Trajectory followed: ",idtraj)
#endif
            CASE(.FALSE.) ! skip dynamics, there may be other trajs to follow
#ifdef DEBUG
               CALL VERBOSE_WRITE(routinename,"Trajectory skipped: ",idtraj) 
#endif
               RETURN
         END SELECT
   END SELECT
   cycles=0 
   zmin = atomo%init_r(3)
   t=0.D0
   atomo%ireb=0
   atomo%ixyboun=0
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
      atom_dofs(1:3)=atomo%r(1:3)
      atom_dofs(4:6)=atomo%p(1:3)
      ! Initial values for the derivatives
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Atomic DOFS : ",atom_dofs)
#endif
      CALL this%TIME_DERIVS(atom_dofs,dfdt,switch)
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Time derivatives: ",dfdt)
#endif
      SELECT CASE(switch)
         CASE(.TRUE.)
            SELECT CASE(cycles)
               CASE(1)
                  WRITE(0,*) "DO_DYNAMICS_DYNATOM ERR: Initial position is not adequate"
                  WRITE(0,*) "Atomo id: ", idtraj
                  WRITE(0,*) "Atom DOF'S: ", atom_dofs(:)
                  CALL EXIT(1)
               CASE DEFAULT
                  WRITE(0,*) "DO_DYNAMICS_DYNATOM ERR: This error is quite uncommon, guess what is happening by your own."
                  WRITE(0,*) "Atom id: ", idtraj
                  WRITE(0,*) "Atom DOF'S: ", atom_dofs(:)
                  CALL EXIT(1)
            END SELECT
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      ! Construct scaling array
      SELECT CASE(this%scaling)
         CASE("Smart")
            FORALL (i=1:6) s(i) = DABS(atom_dofs(i))+DABS(dt*dfdt(i)) ! Numerical Recipes.  
         CASE("Equal")
            s(1:6)=1.D0
         CASE DEFAULT
            WRITE(0,*) "DO_DYNAMICS ERR: incorrect scaling Keyword"
            CALL EXIT(1)
      END SELECT
      ! Energy before time-integration (a.u.)
      init_E = atomo%E
      ! Integrate and choose the correct time-step
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Energy before integration:",atomo%E)
#endif
      ! Call integrator
      CALL this%BSSTEP(atom_dofs,dfdt,t,dt,this%eps,s,dt_did,dt_next,switch)
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Atomic DOFs after integration: ",atom_dofs)
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
      ! At this point, atom_dofs contains the new values for position and momenta
      ! Let's obtain the potential value for this configuration
      CALL this%thispes%GET_V_AND_DERIVS(atom_dofs(1:3),v,dummy)
      E = (atom_dofs(4)**2.D0+atom_dofs(5)**2.D0+atom_dofs(6)**2.D0)/(2.D0*masa)+v
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Energy after integration:",E)
#endif
      SELECT CASE (DABS(E-atomo%E) > 100.D0*this%eps*atomo%E)
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
      ! Check  bouncing points in Z direction (Z-Turning point)
      SELECT CASE((atomo%p(3) < 0.D0).AND.(atom_dofs(6) > 0.D0))
         CASE(.TRUE.)
            atomo%ireb = atomo%ireb +1
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      ! Check bouncing points in XY (sign of Px or Py changes respect to previous integration step)
      SELECT CASE ((DSIGN(atom_dofs(4),atom_dofs(4)).NE.DSIGN(atom_dofs(4),atomo%p(1))).OR. &
                  (DSIGN(atom_dofs(5),atom_dofs(5)).NE.DSIGN(atom_dofs(5),atomo%p(2))))
         CASE(.TRUE.)
            atomo%ixyboun = atomo%ixyboun+1
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      ! Store last minimum Z value reached. This will be the turning-point
      SELECT CASE (atom_dofs(3) <= zmin)
         CASE(.TRUE.)
            zmin = atom_dofs(3)
            atomo%turning_point(1:3) = atom_dofs(1:3)
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
      SELECT CASE ((atom_dofs(3) <= this%zstop%getvalue()+this%dzstop%getvalue()).AND. &
                  (atom_dofs(3) > this%zstop%getvalue()-this%dzstop%getvalue()))
         CASE(.TRUE.)
            atomo%stat = "Stopped"
            atomo%r(1:3) = atom_dofs(1:3)
            atomo%p(1:3) = atom_dofs(4:6)
            atomo%E = E
            WRITE(this%wust,10) idtraj,atomo%stat,atomo%ireb,atomo%ixyboun,atomo%E,t,atomo%r,atomo%p
            EXIT
        CASE(.FALSE.)
           ! do nothing, next switch
      END SELECT
      SELECT CASE ((atomo%r(3) >= this%zscatt%getvalue()).AND.(atomo%p(3) > 0.D0))
         CASE(.TRUE.)
            atomo%stat="Scattered"                                            
            atomo%r(1:3)=atom_dofs(1:3)
            atomo%p(1:3)=atom_dofs(4:6)
            atomo%E=E
            atomo%turning_point(1:2) = this%thispes%surf%project_unitcell(atomo%turning_point(1:2))
            WRITE(this%wusc,10) idtraj,atomo%stat,atomo%ireb,atomo%ixyboun,atomo%E,t,atomo%r,atomo%p
            WRITE(this%wutp,11) idtraj,atomo%turning_point(:)
            EXIT
         CASE(.FALSE.)
            ! do nothing next switch
      END SELECT
      SELECT CASE ((atom_dofs(3) <= this%zabs%getvalue()).AND.(atom_dofs(6) < 0.D0))
         CASE(.TRUE.)
            atomo%stat = "Absorbed"
            atomo%r(1:3) = atom_dofs(1:3)
            atomo%p(1:3) = atom_dofs(4:6)
            atomo%E = E
            WRITE(this%wuab,10) idtraj,atomo%stat,atomo%ireb,atomo%ixyboun,atomo%E,t,atomo%r,atomo%p
            EXIT
         CASE(.FALSE.)
            ! do nothing, next switch
      END SELECT
      SELECT CASE ((v < 0.D0).AND.(atom_dofs(3) <= this%zads%getvalue()).AND.maxtime_reached)
         CASE(.TRUE.)
            atomo%stat = "Adsorbed"
            atomo%r(1:3) = atom_dofs(1:3)
            atomo%p(1:3) = atom_dofs(4:6)
            atomo%E = E
            WRITE(this%wuad,10) idtraj,atomo%stat,atomo%ireb,atomo%ixyboun,atomo%E,t,atomo%r,atomo%p
            EXIT
         CASE(.FALSE.)
            ! do nothing, next switch
      END SELECT
      SELECT CASE (atom_dofs(3) <= this%zads%getvalue() .AND. maxtime_reached)
         CASE(.TRUE.)
            atomo%stat = "Trapped"
            atomo%r(1:3) = atom_dofs(1:3)
            atomo%p(1:3) = atom_dofs(4:6)
            atomo%E = E
            WRITE(this%wutr,10) idtraj,atomo%stat,atomo%ireb,atomo%ixyboun,atomo%E,t,atomo%r,atomo%p
            EXIT
         CASE(.FALSE.)
            ! do nothing, next switch
      END SELECT
      SELECT CASE (maxtime_reached)
         CASE(.TRUE.)
            atomo%stat = "Time-out"
            atomo%r(1:3) = atom_dofs(1:3)
            atomo%p(1:3) = atom_dofs(4:6)
            atomo%E = E
            WRITE(this%wuto,10) idtraj,atomo%stat,atomo%ireb,atomo%ixyboun,atomo%E,t,atomo%r,atomo%p
            EXIT
         CASE(.FALSE.)
            !do nothing next switch
      END SELECT
      SELECT CASE((cycles > 1000).AND.(dt < 1.D-9))
         CASE(.TRUE.)
            atomo%stat = "Patologic"
            atomo%r(1:3) = atom_dofs(1:3)
            atomo%p(1:3) = atom_dofs(4:6)
            atomo%E = E
            WRITE(this%wupa,10) idtraj,atomo%stat,atomo%ireb,atomo%ixyboun,atomo%E,t,atomo%r,atomo%p
            EXIT
         CASE(.FALSE.)
            ! do nothing, next switch
      END SELECT
      SELECT CASE(.NOT.maxtime_reached)
         CASE(.TRUE.)
            atomo%r(1:3) = atom_dofs(1:3)
            atomo%p(1:3) = atom_dofs(4:6)
            atomo%E = E
            dt = dt_next
            SELECT CASE((this%nfollow.NE.0).AND.(in_list))
               CASE(.TRUE.)
                  WRITE(this%wufo,*) t,dt_did,atomo%E,(atomo%p(3)**2.D0)/(2.D0*masa),v,atomo%r(:),atomo%p(:) 
               CASE(.FALSE.)
                  ! do nothing
            END SELECT
            CYCLE
         CASE(.FALSE.)
            WRITE(0,*) "DO_DYNAMICS_DYNATOM ERR: Switches failed to classify this traj"
            CALL EXIT(1)
      END SELECT
   END DO
   SELECT CASE(this%nfollow/=0)
      CASE(.TRUE.)
         CLOSE(this%wufo)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE DO_DYNAMICS_DYNATOM
!###########################################################
!# SUBROUTINE: FILE_TRAJSTATUS_DYNATOM
!###########################################################
!> @brief
!! If file exists, open unit in append mode, else, create the file
!! and let it open
!
!> @param[in] wunit  - unit of the file to be opened
!> @param[in] filename  - name of the file
!> @param[in] title - some short title to append at the topmost line of the file
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Aug/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE FILE_TRAJSTATUS_DYNATOM(wunit,filename,title)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: wunit
   CHARACTER(LEN=*),INTENT(IN) :: filename
   CHARACTER(LEN=*),INTENT(IN) :: title
   ! Local variables
   LOGICAL :: file_exists
   CHARACTER(LEN=24),PARAMETER :: routinename="FILE_TRAJSTATUS_DYNATOM "
   ! Run section
   INQUIRE(FILE=filename,EXIST=file_exists)
   SELECT CASE(file_exists)
      CASE(.TRUE.)
         OPEN(wunit, FILE=filename,STATUS="old",ACCESS="append")
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Previous file found: ",filename)
         CALL VERBOSE_WRITE(routinename,"Appending info to this file")
#endif
      CASE(.FALSE.)
         OPEN(wunit,FILE=filename,STATUS="new")
         WRITE(wunit,*) "# ***** ",title," *****"
         WRITE(wunit,*) "# Format: id/status/ireb/ixyboun/Etot(a.u.)/t(a.u.)/X,Y,Z(a.u.)/Px,Py,Pz(a.u.)"
         WRITE(wunit,*) "# -----------------------------------------------------------"
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"New file created: ",filename)
         CALL VERBOSE_WRITE(routinename,"Header printed")
         CALL VERBOSE_WRITE(routinename,"Appending info to this file")
#endif
   END SELECT
   RETURN
END SUBROUTINE FILE_TRAJSTATUS_DYNATOM
!###########################################################
!# SUBROUTINE: FILE_TURNING_DYNATOM
!###########################################################
!> @brief
!! If  OUTDYN3Dturning.out exists, open unit in append mode, else, create the file
!! and let it open.
!
!> @param[in] wunit  - unit of the file to be opened
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Aug/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE FILE_TURNING_DYNATOM(wunit)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: wunit
   ! Local variables
   CHARACTER(LEN=19),PARAMETER :: filename="OUTDYN3Dturning.out"
   CHARACTER(LEN=34),PARAMETER :: title="TURNING POINTS FOR SCATTERED TRAJS"
   LOGICAL :: file_exists
   CHARACTER(LEN=21),PARAMETER :: routinename="FILE_TURNING_DYNATOM "
   ! Run section
   INQUIRE(FILE=filename,EXIST=file_exists)
   SELECT CASE(file_exists)
      CASE(.TRUE.)
         OPEN(wunit, FILE=filename,STATUS="old",ACCESS="append")
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Previous file found: ",filename)
         CALL VERBOSE_WRITE(routinename,"Appending info to this file")
#endif
      CASE(.FALSE.)
         OPEN(wunit,FILE=filename,STATUS="new")
         WRITE(wunit,*) "# ***** ",title," *****"
         WRITE(wunit,*) "# Description: positions of scattered atoms at their lowest  "
         WRITE(wunit,*) "#              Z value reached during the dynamics. X and Y  "
         WRITE(wunit,*) "#              values are projected in the unit cell.        "
         WRITE(wunit,*) "# Format: id/X,Y,Z(a.u.)                                   "
         WRITE(wunit,*) "# -----------------------------------------------------------"
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"New file created: ",filename)
         CALL VERBOSE_WRITE(routinename,"Header printed")
         CALL VERBOSE_WRITE(routinename,"Appending info to this file")
#endif
   END SELECT
   RETURN
END SUBROUTINE FILE_TURNING_DYNATOM
!###########################################################
!# SUBROUTINE: FILE_FOLLOWTRAJ_DYNATOM
!###########################################################
!> @brief
!! Open file OUTDYN3D'trajid'.out in replace mode. The unit is
!! let opened
!
!> @param[in] wunit  - integer(kind=4): unit of the file to be opened
!> @param[in] idtraj - integer(kind=4): id of the trajectory
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Aug/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE FILE_FOLLOWTRAJ_DYNATOM(wunit,idtraj)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: wunit,idtraj
   ! Local variables
   CHARACTER(LEN=10) :: idstring
   CHARACTER(LEN=26) :: filename
   CHARACTER(LEN=24),PARAMETER :: title="TIME EVOLUTION OF A TRAJ"
   CHARACTER(LEN=24),PARAMETER :: routinename="FILE_FOLLOWTRAJ_DYNATOM "
   ! Run section
   WRITE(idstring,'(I10.10)') idtraj
   filename='OUTDYN3Dtraj'//trim(idstring)//'.out'
   OPEN(wunit,FILE=filename,STATUS="replace",ACTION="write")
   WRITE(wunit,*) "# ***** ",title," *****"
   WRITE(wunit,*) "# Format: t(a.u.)/dt(a.u.)/Etot(a.u.)/Enorm(a.u.)/Pot(a.u.)/X,Y,Z(a.u.)/Px,Py,Pz(a.u.)"
   WRITE(wunit,*) "# ----------------------------------------------------------------"
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"New file created: ",filename)
   CALL VERBOSE_WRITE(routinename,"Header printed")
   CALL VERBOSE_WRITE(routinename,"Appending info to this file")
#endif
   RETURN
END SUBROUTINE FILE_FOLLOWTRAJ_DYNATOM
!###############################################################
!# SUBROUTINE : TIME_DERIVS_DYNATOM ########################
!###############################################################
!> @brief
!! Gives dzdt at z and t values from Hamilton equations of motion
!
!> @param[in] this - Provides some data
!> @param[in] z - array of positions and momenta. z(1:3) -> positions, z(4:6) -> momenta 
!> @param[out] dzdt - time derivatives of position and momenta
!> @param[out] fin - controls errors
!--------------------------------------------------------------
SUBROUTINE TIME_DERIVS_DYNATOM(this,z,dzdt,fin)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Dynatom),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: z
   REAL(KIND=8),DIMENSION(6),INTENT(OUT) :: dzdt
   LOGICAL, INTENT(OUT) :: fin
   ! Local variables
   INTEGER :: i ! counters
   REAL(KIND=8) :: mass
   REAL(KIND=8) :: v ! dummy variable
   CHARACTER(LEN=21),PARAMETER :: routinename = "TIME_DERIVS_DYNATOM: "
   ! ROCK THE CASBAH ! ---------------------
   SELECT CASE(this%thispes%is_allowed(z(1:3)))
      CASE(.FALSE.)
         fin = .TRUE.
      CASE(.TRUE.)
         fin=.FALSE.
         mass=this%thispes%atomdat(1)%getmass() ! there'd be only one atom defined in the PES
         DO i = 1, 3
            dzdt(i) = z(i+3)/mass ! Px/m,  etc...
         END DO
         CALL this%thispes%GET_V_AND_DERIVS(z(1:3),v,dzdt(4:6))
         FORALL (i=4:6) dzdt(i) = -dzdt(i) ! -dV/dx (minus sign comes from here)
   END SELECT
   RETURN
END SUBROUTINE TIME_DERIVS_DYNATOM
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
SUBROUTINE MMID_DYNATOM(this,y,dydx,xs,htot,nstep,yout,switch)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Dynatom),INTENT(IN) :: this
   INTEGER,INTENT(IN) :: nstep
   REAL(KIND=8),DIMENSION(6), INTENT(IN) :: y,dydx
   REAL(KIND=8),DIMENSION(6), INTENT(OUT) :: yout
   REAL(KIND=8),INTENT(IN) :: xs,htot
   LOGICAL,INTENT(INOUT) :: switch
   ! Local variables
   REAL(KIND=8),DIMENSION(6) :: ym,yn
   INTEGER :: i,n ! counters
   INTEGER,PARAMETER :: nvar = 6
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
END SUBROUTINE MMID_DYNATOM
!##################################################################################################
!# SUBROUTINE: RZEXTR_DYNATOM #####################################################################
!################################################################################################## 
!> @brief 
!! - A part of the Burlich-Stoer algorithm. Uses diagonal rational function extrapolation. 
!
!>@details
!! Taken from Numerical recipes in Fortran 77
!> @see pzextr
!--------------------------------------------------------------------------------------------------
SUBROUTINE RZEXTR_DYNATOM(this,iest,xest,yest,yz,dy,nv)
	IMPLICIT NONE
	! I/O variables
   CLASS(Dynatom),INTENT(IN):: this
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
END SUBROUTINE RZEXTR_DYNATOM
!############################################################################################
!# SUBROUTINE : PZEXTR_DYNATOM ##############################################################
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
SUBROUTINE PZEXTR_DYNATOM(this,iest,xest,yest,yz,dy,nv)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Dynatom),INTENT(IN):: this
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
END SUBROUTINE PZEXTR_DYNATOM
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
SUBROUTINE BSSTEP_DYNATOM(this,y,dydx,x,htry,eps,yscal,hdid,hnext,switch)
	IMPLICIT NONE
	! I/O variables
	CLASS(Dynatom),INTENT(IN) :: this
	REAL(KIND=8), INTENT(IN) :: eps     ! required accuracy
	REAL(KIND=8), INTENT(IN) ::  htry   ! step to try
	REAL(KIND=8), DIMENSION(6) :: yscal ! factors to scale error 
	REAL(KIND=8), DIMENSION(6), INTENT(IN) :: dydx 
	REAL(KIND=8), DIMENSION(6), INTENT(INOUT) :: y ! initial/final values for: X,Y,Z,Px,Py,Pz (in this order)
	REAL(KIND=8), INTENT(INOUT) :: x
	LOGICAL, INTENT(INOUT) :: switch ! .TRUE. if potential could not be calculated
	REAL(KIND=8), INTENT(OUT) :: hdid  ! step actually used
	REAL(KIND=8), INTENT(OUT) :: hnext ! guess of the next step
	! Parameters for this routine
	INTEGER, PARAMETER :: nv = 6
	REAL(KIND=8),PARAMETER :: SAFE1 = 0.25D0
	REAL(KIND=8),PARAMETER :: SAFE2 = 0.7D0
	REAL(KIND=8),PARAMETER :: TINY = 1.D-30 
	REAL(KIND=8),PARAMETER :: SCALMX = 0.1D0
	REAL(KIND=8),PARAMETER :: REDMIN = 0.7D0
	REAL(KIND=8),PARAMETER :: REDMAX = 1.D-5
	INTEGER,PARAMETER :: NMAX = 50
	INTEGER,PARAMETER :: KMAXX = 8
	INTEGER,PARAMETER :: IMAX = KMAXX+1
	CHARACTER(LEN=16), PARAMETER :: routinename = "BSSTEP_DYNATOM: "
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
END SUBROUTINE BSSTEP_DYNATOM
END MODULE DYNATOM_MOD
