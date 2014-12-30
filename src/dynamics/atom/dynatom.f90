!##########################################################
! MODULE: DYNATOM_MOD
!> @brief
!! Provides tools to run dynamics on atoms
!##########################################################
MODULE DYNATOM_MOD
   USE DYNAMICS_MOD
   USE INITATOM_MOD
   USE AOTUS_MODULE, ONLY: flu_State, OPEN_CONFIG_FILE, CLOSE_CONFIG, AOT_GET_VAL
   USE AOT_TABLE_MODULE, ONLY: AOT_TABLE_OPEN, AOT_TABLE_CLOSE, AOT_TABLE_LENGTH
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
   CHARACTER(LEN=20):: string
   CHARACTER(LEN=255),DIMENSION(:),ALLOCATABLE:: key
   CHARACTER(LEN=255),DIMENSION(:),ALLOCATABLE:: rawdata
   CHARACTER(LEN=20):: filenameinicond,filenamepes
   CHARACTER(LEN=6):: follow
   CHARACTER(LEN=10):: units
   INTEGER(KIND=4):: iErr,dyn_table,magnit_table
   TYPE(flu_State):: conf
   REAL(KIND=8):: aux
   INTEGER(KIND=4):: i ! counters
   ! Parameters 
   CHARACTER(LEN=20),PARAMETER :: routinename = "INITIALIZE_DYNATOM: "
   ! YIPEE KI YAY -----------------------------
   this%filename=filename
   CALL OPEN_CONFIG_FILE(L=conf,filename=this%filename,ErrCode=iErr)
   SELECT CASE (iErr)
      CASE(0)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Lua config file was opened successfuly: ",this%filename)
#endif
         ! do nothing
      CASE DEFAULT
         WRITE(0,*) 'INITIALIZE_DYNATOM: Fatal Error when opening the Lua config file'
         CALL EXIT(1)
   END SELECT
   CALL AOT_TABLE_OPEN(L=conf,thandle=dyn_table,key='dynamics')
   SELECT CASE(dyn_table)
      CASE(0)
         WRITE(0,*) "READ_DYNATOM ERR: empty dynamics table"
         CALL EXIT(1)
      CASE DEFAULT
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Dynamics table was opened successfuly")
#endif
         ! do nothing
   END SELECT
   ALLOCATE(key(aot_table_length(L=conf,thandle=dyn_table)))
   ALLOCATE(rawdata(size(key)))
   ! Detect kind of dynamics
   CALL AOT_GET_VAL(L=conf,val=string,thandle=dyn_table,ErrCode=iErr,key='kind')
   this%kind=trim(string)
   SELECT CASE(this%kind)
      CASE("Atom")
         ! do nothing
      CASE DEFAULT
         WRITE(0,*) "INITIALIZE_DYNATOM ERR: Expected Atom for dynamics.kind"
         CALL EXIT(1)
   END SELECT
   !
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
      CALL this%thisintegrator%SET_ERRSCALING(s(:))
      ! Energy before time-integration (a.u.)
      init_E = atomo%E
      ! Integrate and choose the correct time-step
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Energy before integration:",atomo%E)
#endif
      ! Call integrator
      CALL this%thisintegrator%INTEGRATE(atom_dofs,dfdt,t,TIME_DERIVS_DYNATOM,switch)
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
END MODULE DYNATOM_MOD
