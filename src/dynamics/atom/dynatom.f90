!##########################################################
! MODULE: DYNATOM_MOD
!> @brief
!! Provides tools to run dynamics on atoms
!> @todo
!! - Generalize the use of different PES, not only CRP
!##########################################################
MODULE DYNATOM_MOD
USE DYNAMICS_MOD
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
   REAL*8 :: eps
   TYPE(Length) :: zstop, dzstop
   TYPE(Length) :: zscatt,zads,zabs
   CONTAINS
      PROCEDURE,PUBLIC :: READ => READ_DYNATOM
      PROCEDURE,PUBLIC :: RUN => RUN_DYNAMICS_ATOMS
END TYPE Dynatom
!//////////////////////////////////////////////////////////
CONTAINS
!###############################################################
!# SUBROUTINE: READ_DYNATOM ####################################
!###############################################################
!> @brief 
!! Reads from file specifications for a dynamics job (atoms)
!
!---------------------------------------------------------------
SUBROUTINE READ_DYNATOM(dinamica,filename)
#ifdef DEBUG
   USE DEBUG_MOD
#endif
	IMPLICIT NONE
	! I/O variables
	CLASS(Dynatom),INTENT(OUT) :: dinamica
	CHARACTER(LEN=*),INTENT(IN) :: filename
	! Local variables
	CHARACTER(LEN=14), PARAMETER :: routinename = "READ_DYNATOM: "
	CHARACTER(LEN=6) :: follow
   CHARACTER(LEN=10) :: units
   REAL(KIND=8) :: aux
	INTEGER :: i ! counters
	! YIPEE KI YAY -----------------------------
	dinamica%filename = filename
	OPEN (11,FILE=dinamica%filename, STATUS="old")
	READ(11,*) dinamica%alias
	READ(11,*) dinamica%kind
	IF(dinamica%kind.EQ."Atoms") THEN
#ifdef DEBUG
		CALL VERBOSE_WRITE(routinename,"Read: precision of integration (dimensionless factor)")
#endif
		READ(11,*) dinamica%eps
		READ(11,*) dinamica%scaling
		IF((dinamica%scaling.NE."Equal").AND.(dinamica%scaling.NE."Smart")) THEN
			WRITE(0,*) "INITIALIZE_DYNATOM ERR: Wrong error scaling keyword: "
			WRITE(0,*) "Only available: Equal and Smart"
			WRITE(0,*) "You have written: ", dinamica%scaling
			CALL EXIT(1)
		END IF
		READ(11,*) dinamica%extrapol
		IF ((dinamica%extrapol.NE."Polinomi").AND.(dinamica%extrapol.NE."Rational")) THEN
			WRITE(0,*) "INITIALIZE_DYNATOM ERR: Wrong extrapolation keyword: "
			WRITE(0,*) "Only available: Polinomi and Rational"
			WRITE(0,*) "You have written: ", dinamica%extrapol
			CALL EXIT(1)
		END IF
		
      READ(11,*) aux,units
		CALL dinamica%delta_t%READ(aux,units)
		CALL dinamica%delta_t%TO_STD()
		
      READ(11,*) aux,units
		CALL dinamica%max_t%READ(aux,units)
		CALL dinamica%max_t%TO_STD()
		
      READ(11,*) aux,units
		CALL dinamica%zstop%READ(aux,units)
		CALL dinamica%zstop%TO_STD()
		
      READ(11,*) aux,units
		CALL dinamica%dzstop%READ(aux,units)
		CALL dinamica%dzstop%TO_STD()
		
      READ(11,*) aux,units
		CALL dinamica%zscatt%READ(aux,units)
		CALL dinamica%zscatt%TO_STD()
		
      READ(11,*) aux,units
		CALL dinamica%zads%READ(aux,units)
		CALL dinamica%zads%TO_STD()
		
      READ(11,*) aux,units
		CALL dinamica%zabs%READ(aux,units)
		CALL dinamica%zabs%TO_STD()
#ifdef DEBUG
		CALL VERBOSE_WRITE(routinename,"Initial time step in au: ",dinamica%delta_t%getvalue())
		CALL VERBOSE_WRITE(routinename,"Maximum time in au: ",dinamica%max_t%getvalue())
		CALL VERBOSE_WRITE(routinename,"ZSTOP in au: ",dinamica%zstop%getvalue())
		CALL VERBOSE_WRITE(routinename,"DZSTOP in au: ",dinamica%dzstop%getvalue())
		CALL VERBOSE_WRITE(routinename,"Z scattering benchmark in au: ",dinamica%zscatt%getvalue())
		CALL VERBOSE_WRITE(routinename,"Z adsorption benchmark in au: ",dinamica%zads%getvalue())
		CALL VERBOSE_WRITE(routinename,"Z absortion benchmark in au: ",dinamica%zabs%getvalue())
#endif
		READ(11,*) follow, dinamica%nfollow
		IF (follow.NE."FOLLOW") THEN
			WRITE(0,*) "INITIALIZE_DYNATOM ERR: wrong FOLLOW keyword"
			CALL EXIT(1)
		ELSE IF (dinamica%nfollow.NE.0) THEN
#ifdef DEBUG
			CALL VERBOSE_WRITE(routinename,"There are trajectories to follow")
#endif
			ALLOCATE(dinamica%follow_atom(1:dinamica%nfollow))
			DO i=1, dinamica%nfollow
				READ(11,*) dinamica%follow_atom(i)
			END DO
		ELSE
#ifdef DEBUG
			CALL VERBOSE_WRITE(routinename,"No trajectories to follow")
#endif
		END IF
	END IF
	CLOSE(11)
	RETURN
END SUBROUTINE READ_DYNATOM
!###############################################################
!# SUBROUTINE: RUN_DYNAMICS_ATOMS ##############################
!###############################################################
!> @brief
!! Launch trajectories as said in @b inicondat variable
!
!---------------------------------------------------------------
SUBROUTINE RUN_DYNAMICS_ATOMS(dinamica,inicondat,thispes,trajs)
   USE INITATOM_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   USE CRP_MOD
   IMPLICIT NONE
   ! I/O variables 
   CLASS(Dynatom),INTENT(IN) :: dinamica
   TYPE(Initatom),INTENT(IN) :: inicondat
   TYPE(CRP),INTENT(IN) :: thispes
   TYPE(Atom_trajs),INTENT(INOUT) :: trajs
   ! Local variables 
   INTEGER :: i ! counters
   CHARACTER(LEN=20), PARAMETER :: routinename = "RUN_DYNAMICS_ATOMS: "
   ! HEY HO! LET'S GO !!! ------
   DO i=inicondat%nstart, inicondat%ntraj
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Start trajectory: ", i)
#endif
      CALL DO_DYNAMICS_ATOM(dinamica,inicondat,thispes,trajs%atomo(i))
   END DO
   RETURN
END SUBROUTINE RUN_DYNAMICS_ATOMS
!##############################################################
!# SUBROUTINE : DO DYNAMICS ATOM ##############################
!##############################################################
!> @brief
!! Integrates equation of motion following specifications stored
!! in @b dinamica and @b inicondat
!--------------------------------------------------------------
SUBROUTINE DO_DYNAMICS_ATOM(dinamica,inicondat,thispes,atomo)
   USE INITATOM_MOD
   USE CRP_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
	IMPLICIT NONE
	! I/O variables
	TYPE(Dynatom),INTENT(IN) :: dinamica
	TYPE(Initatom),INTENT(IN) :: inicondat
   TYPE(CRP),INTENT(IN) :: thispes
	TYPE(Atom),INTENT(INOUT) :: atomo
	! Local variables
	INTEGER :: i, cycles ! counters
	REAL*8 :: t, dt, E, init_t, init_E, v, dt_did, dt_next, zmin, angle
	REAL*8, DIMENSION(3) :: r0, p0, dummy
	REAL*8, DIMENSION(6) :: atom_dofs, s, dfdt
	LOGICAL :: switch, file_exists, in_list
	CHARACTER(LEN=18), PARAMETER :: routinename = "DO_DYNAMICS_ATOM: "
	CHARACTER(LEN=18) :: filename_follow
	CHARACTER(LEN=9), PARAMETER :: format_string = '(I10.10)'
	CHARACTER(LEN=10) :: x1
	INTEGER :: control
	! HEY HO!, LET'S GO!!! -------------------------
	INQUIRE(FILE="dynamics.out",EXIST=file_exists)
	IF(file_exists) THEN 
		OPEN(11, FILE="dynamics.out",STATUS="old",ACCESS="append")
	ELSE
		OPEN(11,FILE="dynamics.out",STATUS="new")
		WRITE(11,*) "# DYNAMICS RESULTS -----------------------------------------------------------"
		WRITE(11,*) "# Format: id, status, ireb, ixyboun, Etot, t, X,Y,Z (a.u.), Px,Py,Pz (a.u.)"
		WRITE(11,*) "#------------------------------------------------------------------------------"
	END IF
	INQUIRE(FILE="turning.out",EXIST=file_exists)
	IF(file_exists) THEN
		OPEN(12, FILE="turning.out",STATUS="old",ACCESS="append")
	ELSE
		OPEN(12,FILE="turning.out",STATUS="new")
		WRITE(12,*) "# TURNING POINTS --------------------------------------------"
		WRITE(12,*) "# Description: positions of scattered atoms at their lowest"
		WRITE(12,*) "#              Z value reached during the dynamics."
		WRITE(12,*) "#               XY values, projected into IWS cell"
		WRITE(12,*) "# Format: id, X,Y,Z (a.u.)"
		WRITE(12,*) "#-------------------------------------------------------------"
	END IF
	in_list = .FALSE.
	IF(dinamica%nfollow.NE.0) THEN
		DO i=1,dinamica%nfollow
			IF (dinamica%follow_atom(i).EQ.atomo%id) in_list = .TRUE.
		END DO
		IF (in_list) THEN
			WRITE(x1,format_string) atomo%id
			filename_follow = 'traj'//TRIM(x1)//'.out'
			OPEN(13,FILE=filename_follow,STATUS="replace")
			WRITE(13,*) "# TIME EVOLUTION FOR A TRAJECTORY ---------------------"
			WRITE(13,*) "# Format: t, dt, Etot,Enorm,Pot(a.u) X,Y,Z(a.u.) Px,Py,Pz (a.u.)  "
			WRITE(13,*) "# First and last position are not printed here"
			WRITE(13,*) "# You can find them in dynamics.out and inicond.out    "
			WRITE(13,*) "#------------------------------------------------------"
		END IF
	END IF
	cycles = 0 
	zmin = atomo%init_r(3)
	t=0.D0
	atomo%ireb = 0
	atomo%ixyboun = 0
	dt = dinamica%delta_t%getvalue()
	! Iterations
	DO
		switch = .FALSE.
		cycles = cycles+1
		! We cannot go beyond maximum time defined
		IF (t+dt.GT.dinamica%max_t%getvalue()) THEN
			dt = dinamica%max_t%getvalue()-t
		END IF
		init_t = t
		! Storing atomic DOF's in only one array
		FORALL (i=1:3) 
			atom_dofs(i) = atomo%r(i)
			atom_dofs(i+3) = atomo%p(i)
		END FORALL 
      ! Initial values for the derivatives
		CALL ATOM_H_DERIVS(inicondat,thispes,t,atom_dofs,dfdt,switch)
		IF (switch) THEN
			IF (cycles.EQ.1) THEN ! first cycle
				WRITE(0,*) "DO_DYNAMICS_ATOMS ERR: Initial position is not adequate"
				WRITE(0,*) "Atomo id: ", atomo%id
				WRITE(0,*) "Atom DOF'S: ", atom_dofs(:)
				STOP 1
			ELSE ! other steps
				WRITE(0,*) "DO_DYNAMICS_ATOMS ERR: This error is quite uncommon, guess what is happening by your own."
				WRITE(0,*) "Atom id: ", atomo%id
				WRITE(0,*) "Atom DOF'S: ", atom_dofs(:)
				STOP 1
			END IF
		END IF
		! Construct scaling array
		IF(dinamica%scaling.EQ."Smart") FORALL (i=1:6) s(i) = DABS(atom_dofs(i))+DABS(dt*dfdt(i)) ! Numerical Recipes.  
		IF(dinamica%scaling.EQ."Equal") FORALL (i=1:6) s(i) = 1.D0 ! same scaling
		! Energy before time-integration (a.u.)
		init_E = atomo%E
		! Integrate and choose the correct time-step
		CALL BSSTEP_ATOM(inicondat,dinamica,thispes,atom_dofs,dfdt,t,dt,dinamica%eps,s,dt_did,dt_next,switch)
		IF(switch) THEN
			! Problem detected in integration, reducing time-step
			! reboot to previous step values 
			dt = dt/2.D0 
			t = init_t
			CYCLE
		END IF
		! At this point, atom_dofs contains the new values for position and momenta
		! Let's obtain the potential value for this configuration
		CALL thispes%GET_V_AND_DERIVS(atom_dofs(1:3),v,dummy) ! just to  obtain the potential, should work
		E = (atom_dofs(4)**2.D0+atom_dofs(5)**2.D0+atom_dofs(6)**2.D0)/(2.D0*inicondat%masss%getvalue())+v
		IF (DABS(E-atomo%E).GT.100.D0*dinamica%eps*atomo%E) THEN
			! Problems with energy conservation
			! reboot to previous step values
			dt = dt/2.D0
			t = init_t
			CYCLE
		END IF
		! Check  bouncing points in Z direction (Z- Turning point)
		IF ((atomo%p(3).LT.0.D0).AND.(atom_dofs(6).GT.0.D0)) THEN
			atomo%ireb = atomo%ireb +1
		END IF
		! Check bouncing points in XY (sign of Px or Py changes respect to previous integration step)
		IF((DSIGN(atom_dofs(4),atom_dofs(4)).NE.DSIGN(atom_dofs(4),atomo%p(1))).OR. &
		   (DSIGN(atom_dofs(5),atom_dofs(5)).NE.DSIGN(atom_dofs(5),atomo%p(2)))) THEN
			atomo%ixyboun = atomo%ixyboun+1
		END IF
		! Store last minimum Z value reached. This will be the turning-point
		IF (atom_dofs(3).LE.zmin) THEN
			zmin = atom_dofs(3)
			FORALL(i=1:3) atomo%turning_point(i) = atom_dofs(i)
		END IF
		! EXIT CHANNELS -------------------------------------------
		IF ((atom_dofs(3).LE.dinamica%zstop%getvalue()+dinamica%dzstop%getvalue()).AND. &
		   (atom_dofs(3).GT.dinamica%zstop%getvalue()-dinamica%dzstop%getvalue())) THEN
			atomo%stat = "Stopped"
			FORALL (i=1:3) 
				atomo%r(i) = atom_dofs(i)
				atomo%p(i) = atom_dofs(i+3)
			END FORALL
			atomo%E = E
			EXIT
		ELSE IF ((atomo%r(3).GE.dinamica%zscatt%getvalue()).AND.(atomo%p(3).GT.0.D0)) THEN ! Scattering: Z position equal or higher than initial Z
			atomo%stat = "Scattered"                                            ! and Pz pointing to the vacuum
			FORALL (i=1:3) 
				atomo%r(i) = atom_dofs(i)
				atomo%p(i) = atom_dofs(i+3)
			END FORALL
			atomo%E = E
         atomo%turning_point(1:2) = thispes%surf%project_iwscell(atomo%turning_point(1:2))
			WRITE(12,'(I7,1X,3(F10.5,1X))') atomo%id, atomo%turning_point(:)
			EXIT
		ELSE IF ((atom_dofs(3).LE.dinamica%zabs%getvalue()).AND.(atom_dofs(6).LT.0.D0)) THEN
			atomo%stat = "Absorbed"
			FORALL (i=1:3) 
				atomo%r(i) = atom_dofs(i)
				atomo%p(i) = atom_dofs(i+3)
			END FORALL
			atomo%E = E
			EXIT
		ELSE IF (t.GE.dinamica%max_t%getvalue()) THEN ! exit channels evaluated when time is out
			IF ((v.LT.0.D0).AND.(atom_dofs(3).LE.dinamica%zads%getvalue())) THEN ! Adsorption: interaction potential v lower than zero and Z lower than a certain value 
				atomo%stat = "Adsorbed"
				FORALL (i=1:3) 
					atomo%r(i) = atom_dofs(i)
					atomo%p(i) = atom_dofs(i+3)
				END FORALL
				atomo%E = E
				EXIT
			ELSE IF (atom_dofs(3).LE.dinamica%zads%getvalue()) THEN
				atomo%stat = "Trapped"
				FORALL (i=1:3) 
					atomo%r(i) = atom_dofs(i)
					atomo%p(i) = atom_dofs(i+3)
				END FORALL
				atomo%E = E
				EXIT
			ELSE
				atomo%stat = "Time-out"
				FORALL (i=1:3) 
					atomo%r(i) = atom_dofs(i)
					atomo%p(i) = atom_dofs(i+3)
				END FORALL
				atomo%E = E
				EXIT
			END IF
		ELSE IF ((cycles.GT.1000).AND.(dt.LT.1.D-9)) THEN ! Too many cycles and low delta-t
			atomo%stat = "Patologic"
			FORALL (i=1:3)
				atomo%r(i) = atom_dofs(i)
				atomo%p(i) = atom_dofs(i+3)
			END FORALL
			atomo%E = E
			EXIT
		ELSE ! max_time not reached, unclassificable traj => dyniamics is not finished, cycle
			FORALL (i=1:3) 
				atomo%r(i) = atom_dofs(i)
				atomo%p(i) = atom_dofs(i+3)
			END FORALL
			atomo%E = E
			dt = dt_next
			IF((dinamica%nfollow.NE.0).AND.(in_list)) WRITE(13,*) t, dt_did, atomo%E, &
				(atomo%p(3)**2.D0)/(2.D0*inicondat%masss%getvalue()), v, atomo%r(:), atomo%p(:) 
			CYCLE
		END IF
	END DO
#ifdef DEBUG
	CALL VERBOSE_WRITE(routinename,"Writing in dynamics.out")
#endif
	WRITE(11,'(I7,1X,A10,1X,I4,1X,I5,1X,8(F15.5,1X))') atomo%id, atomo%stat, atomo%ireb, atomo%ixyboun, atomo%E, t, atomo%r, atomo%p
	CLOSE(11)
	CLOSE(12)
	IF(dinamica%nfollow.NE.0) CLOSE(13)
	RETURN
END SUBROUTINE DO_DYNAMICS_ATOM
!###############################################################
!# SUBROUTINE : ATOM_H_DIFFERENTIAL_EQ #########################
!###############################################################
!> @brief
!! Gives dzdt at z and t values from Hamilton equations of motion
!
!> @param[in] inicondat - Initial conditions needed
!> @param[in] thispes - PES used
!> @param[in] z - array of positions and momenta. z(1:3) -> positions, z(4:6) -> momenta 
!> @param[out] dzdt - equations of motion for positions and momenta
!> @param[out] fin - controls errors
!--------------------------------------------------------------
SUBROUTINE ATOM_H_DERIVS(inicondat,thispes,t,z,dzdt,fin)
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   USE INITATOM_MOD
   USE CRP_MOD
	IMPLICIT NONE
	! I/O variables
	TYPE(Initatom), INTENT(IN) :: inicondat
   TYPE(CRP),INTENT(IN) :: thispes
	REAL*8, DIMENSION(6), INTENT(IN) :: z
	REAL*8, DIMENSION(6), INTENT(OUT) :: dzdt
	REAL*8, INTENT(IN) :: t ! time?
	LOGICAL, INTENT(OUT) :: fin
	! Local variables
	INTEGER :: i ! counters
	REAL*8 :: v ! really a dummy variable
	CHARACTER(LEN=15), PARAMETER :: routinename = "ATOM_H_DERIVS: "
	LOGICAL :: deb_mode, verb_mode
	! ROCK THE CASBAH ! ---------------------
	fin = .FALSE. ! initial value
   IF (.NOT.thispes%is_allowed(z(1:3))) THEN
		fin = .TRUE.
		RETURN
   END IF
   DO i = 1, 3
	   dzdt(i) = z(i+3)/inicondat%masss%getvalue() ! Px/m,  etc...
   END DO
#ifdef DEBUG
	deb_mode = get_debugmode()
	verb_mode = get_verbosemode()
   CALL SET_VERBOSE_MODE(.FALSE.) ! we don't want much output
	CALL SET_DEBUG_MODE(.FALSE.)
#endif
	CALL thispes%GET_V_AND_DERIVS(z(1:3),v,dzdt(4:6))
	FORALL (i=4:6) dzdt(i) = -dzdt(i) ! -dV/dx (minus sign comes from here)
#ifdef DEBUG
	CALL SET_VERBOSE_MODE(verb_mode)
	CALL SET_DEBUG_MODE(deb_mode)
#endif
	RETURN
END SUBROUTINE ATOM_H_DERIVS
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
!> @param[in] inicondat - We need information stored in this variable
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
SUBROUTINE MMID_ATOM(inicondat,thispes,y,dydx,xs,htot,nstep,yout,switch)
   USE CRP_MOD
   USE INITATOM_MOD
	IMPLICIT NONE
	! I/O variables
   TYPE(CRP),INTENT(IN) :: thispes
	TYPE(Initatom), INTENT(IN) :: inicondat
	INTEGER, INTENT(IN) :: nstep
	REAL*8, DIMENSION(6), INTENT(IN) :: y,dydx
	REAL*8, DIMENSION(6), INTENT(OUT) :: yout
	REAL*8, INTENT(IN) :: xs,htot
	LOGICAL,INTENT(INOUT) :: switch
	! Local variables
	REAL*8, DIMENSION(6) :: ym,yn
	INTEGER :: i,n ! counters
	INTEGER, PARAMETER :: nvar = 6
	REAL*8 :: h,h2,swap,x
	! ROCK THE CASBAH !!! ---------------------
	h=htot/DFLOAT(nstep) ! Stepsize this trip.
	DO i=1,nvar
		ym(i)=y(i)
		yn(i)=y(i)+h*dydx(i) ! First step.
	END DO
	x=xs+h
	CALL ATOM_H_DERIVS(inicondat,thispes,x,yn,yout,switch) ! Will use yout for temporary storage of derivatives.
	IF (switch) RETURN
	h2=2.D0*h
	DO n=2,nstep ! General step.
		DO i=1,nvar
			swap=ym(i)+h2*yout(i)
			ym(i)=yn(i)
			yn(i)=swap
		END DO
		x=x+h
		CALL ATOM_H_DERIVS(inicondat,thispes,x,yn,yout,switch)
		IF (switch) RETURN
	END DO
	DO i=1,nvar ! Last step.
		yout(i)=0.5D0*(ym(i)+yn(i)+h*yout(i))
	END DO 
	RETURN
END SUBROUTINE MMID_ATOM
!##################################################################################################
!# SUBROUTINE: RZEXTR #############################################################################
!################################################################################################## 
!> @brief 
!! - A part of the Burlich-Stoer algorithm. Uses diagonal rational function extrapolation. 
!
!>@details
!! Taken from Numerical recipes in Fortran 77
!> @see pzextr
!--------------------------------------------------------------------------------------------------
SUBROUTINE RZEXTR (iest,xest,yest,yz,dy,nv)
	IMPLICIT NONE
	! I/O variables
	INTEGER, INTENT(IN) :: iest, nv
	REAL*8, INTENT(IN) :: xest
	REAL*8, DIMENSION(nv), INTENT(IN) :: yest
	REAL*8, DIMENSION(nv), INTENT(OUT) :: dy, yz
	! Local variables 	
	INTEGER, PARAMETER :: IMAX = 13
	INTEGER, PARAMETER :: NMAX = 50
	INTEGER :: j,k
	REAL*8, DIMENSION(NMAX,IMAX) :: d
	REAL*8, DIMENSION(IMAX) :: fx, x
	REAL*8 :: b,b1,c,ddy,v,yy
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
END SUBROUTINE RZEXTR
!############################################################################################
!# SUBROUTINE : PZEXTR ######################################################################
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
SUBROUTINE PZEXTR(iest,xest,yest,yz,dy,nv)
	IMPLICIT NONE
	! I/O variables
	INTEGER, INTENT(IN) :: iest, nv
	REAL*8, INTENT(IN) :: xest
	REAL*8, DIMENSION(nv), INTENT(IN) :: yest
	REAL*8, DIMENSION(nv), INTENT(OUT) :: dy, yz
	! Local Variables
	INTEGER, PARAMETER :: IMAX = 13
	INTEGER, PARAMETER :: NMAX = 50 
	INTEGER :: j, k1 ! counters
	REAL*8, DIMENSION(NMAX) :: d
	REAL*8, DIMENSION(IMAX) :: x
	REAL*8, DIMENSION(NMAX,IMAX) :: qcol
	REAL*8 :: delta,f1,f2,q
	! ROCK THE CASBAH !!!! --------------
	SAVE qcol,x
	x(iest)=xest !Save current independent variable.
	DO j=1,nv
		dy(j)=yest(j)
		yz(j)=yest(j)
	END DO
	IF (iest.eq.1) THEN ! Store rst estimate in rst column.
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
END SUBROUTINE PZEXTR
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
SUBROUTINE BSSTEP_ATOM (inicondat,dinamica,thispes,y,dydx,x,htry,eps,yscal,hdid,hnext,switch)
   USE INITATOM_MOD
   USE CRP_MOD
	IMPLICIT NONE
	! Parameters for this routine
	INTEGER, PARAMETER :: nv = 6
	REAL*8, PARAMETER :: SAFE1 = 0.25D0
	REAL*8, PARAMETER :: SAFE2 = 0.7D0
	REAL*8, PARAMETER :: TINY = 1.D-30 
	REAL*8, PARAMETER :: SCALMX = 0.1D0
	REAL*8, PARAMETER :: REDMIN = 0.7D0
	REAL*8, PARAMETER :: REDMAX = 1.D-5
	INTEGER, PARAMETER :: NMAX = 50
	INTEGER, PARAMETER :: KMAXX = 8
	INTEGER, PARAMETER :: IMAX = KMAXX+1
	CHARACTER(LEN=13), PARAMETER :: routinename = "BSSTEP_ATOM: "
	! I/O variables
	TYPE(Initatom),INTENT(IN) :: inicondat
	TYPE(Dynatom),INTENT(IN) :: dinamica
   TYPE(CRP),INTENT(IN) :: thispes
	REAL*8, INTENT(IN) :: eps     ! required accuracy
	REAL*8, INTENT(IN) ::  htry   ! step to try
	REAL*8, DIMENSION(6) :: yscal ! factors to scale error 
	REAL*8, DIMENSION(6), INTENT(IN) :: dydx 
	REAL*8, DIMENSION(6), INTENT(INOUT) :: y ! initial/final values for: X,Y,Z,Px,Py,Pz (in this order)
	REAL*8, INTENT(INOUT) :: x
	LOGICAL, INTENT(INOUT) :: switch ! .TRUE. if potential could not be calculated
	REAL*8, INTENT(OUT) :: hdid  ! step actually used
	REAL*8, INTENT(OUT) :: hnext ! guess of the next step
	! Local variables
	INTEGER, DIMENSION(IMAX) :: nseq
	REAL*8, DIMENSION(KMAXX) :: err
	REAL*8, DIMENSION(NMAX) :: yerr, ysav, yseq
	REAL*8, DIMENSION(IMAX) :: a
	REAL*8, DIMENSION(KMAXX,KMAXX) :: alf
	INTEGER :: i,iq,k,kk,km,kmax,kopt
	REAL*8 :: eps1,epsold,errmax,fact,h,red,scale,work,wrkmin,xest, xnew
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
		CALL MMID_ATOM(inicondat,thispes,ysav,dydx,x,h,nseq(k),yseq,switch)
		IF(switch) RETURN
		xest=(h/nseq(k))**2.D0 ! Squared, since error series is even.
		IF(dinamica%extrapol.EQ."Rational") THEN
			CALL RZEXTR (k,xest,yseq,y,yerr,nv) ! Perform extrapolation with Rational funcions
		ELSE IF (dinamica%extrapol.EQ."Polinomi") THEN
			CALL PZEXTR(k,xest,yseq,y,yerr,nv) ! Perform extrapolation with Polinoms
		ELSE
			WRITE(0,*) "BSSTEP_ATOM ERR: Wrong keyword for extrapolation variable"
			STOP 1
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
END SUBROUTINE BSSTEP_ATOM
END MODULE DYNATOM_MOD
