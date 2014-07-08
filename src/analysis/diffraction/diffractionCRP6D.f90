MODULE DIFFRACTIONCRP6D_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
USE INITDIATOMIC_MOD
USE CONSTANTS_MOD
USE SURFACE_MOD
USE CRP6D_MOD
IMPLICIT NONE
!======================================================
! Peak derived data
!---------------------
TYPE :: PeakCRP6D
   PRIVATE
   INTEGER(KIND=4) :: id ! identification number
   INTEGER(KIND=4) :: order ! order of the peak
   INTEGER(KIND=4), DIMENSION(2) :: g ! coordinates in the reciprocal space lattice
	REAL(KIND=8) :: Psi ! azimuthal exit angle
	REAL(KIND=8) :: Phi ! deflection angle respect to incidence plane
	REAL(KIND=8) :: Theta_out ! deflection angle respect to surface plane
	REAL(KIND=8) :: Prob ! probability
END TYPE PeakCRP6D
!======================================================
! Allowed_peaksCRP6D derived data
!----------------------------
TYPE :: Allowed_peaksCRP6D
   PRIVATE
   TYPE(Surface):: surf
   TYPE(Initdiatomic):: inicond
   TYPE(CRP6D):: thispes
   REAL(KIND=8):: E
   REAL(KIND=8),DIMENSION(6):: conic
   TYPE(PeakCRP6D),DIMENSION(:),ALLOCATABLE:: peaks
   CONTAINS
      PROCEDURE,PUBLIC:: INITIALIZE => INITIALIZE_ALLOWEDPEAKSCRP6D
      PROCEDURE,PUBLIC:: SETUP => SETUP_ALLOWEDPEAKSCRP6D
      PROCEDURE,PUBLIC:: ASSIGN_PEAKS => ASSIGN_PEAKS_TO_TRAJS_ALLOWEDPEAKSCRP6D
      PROCEDURE,PUBLIC:: PRINT_XY_EXIT_ANGLES => PRINT_XY_EXIT_ANGLES_ALLOWEDPEAKSCRP6D
      PROCEDURE,PUBLIC:: evaluate_peak => evaluate_peak_ALLOWEDPEAKSCRP6D
END TYPE Allowed_peaksCRP6D
!=======================================================
CONTAINS
SUBROUTINE INITIALIZE_ALLOWEDPEAKSCRP6D(this,surfname,inicondname)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Allowed_peaksCRP6D),INTENT(OUT):: this
   CHARACTER(LEN=*),INTENT(IN) :: surfname,inicondname
   ! Local variables
   LOGICAL :: exists
   ! Run section
   CALL this%surf%INITIALIZE(surfname)
   CALL this%thispes%INITIALIZE("INcrp6d.inp")
   CALL this%inicond%INITIALIZE(inicondname)
   INQUIRE(FILE="OUTinicond6d.inp",EXIST=exists)
   SELECT CASE(exists)
      CASE(.TRUE.)
         CALL this%inicond%GENERATE_TRAJS_FROM_FILE("OUTinicond6d.inp")
      CASE(.FALSE.)
         CALL this%inicond%GENERATE_TRAJS(this%thispes)
   END SELECT
   RETURN
END SUBROUTINE 
!######################################################
!# SUBROUTINE : SETUP_ALLOWEDPEAKSCRP6D #####################
!######################################################
! Just for a square primitive cell
!------------------------------------------------------
SUBROUTINE SETUP_ALLOWEDPEAKSCRP6D(this)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Allowed_peaksCRP6D),INTENT(INOUT) :: this
   ! Local variables
	REAL(KIND=8), DIMENSION(2) :: kinit_par
	REAL(KIND=8), DIMENSION(3) :: p ! momentum
	LOGICAL, DIMENSION(:), ALLOCATABLE :: allowed
	LOGICAL :: need_cycle
	REAL(KIND=8) :: pinit_par
	REAL(KIND=8) :: E ! Total energy
	REAL(KIND=8) :: gamma ! angle between surface main vectors.
	REAL(KIND=8) :: Psi ! azimuthal exit angle
	REAL(KIND=8) :: Phi ! deflection angle respect to incidence plane
	REAL(KIND=8) :: Theta_out ! deflection angle respect to surface plane
   REAL(KIND=8) :: mass
	INTEGER(KIND=4) :: i, k ! Counters
	INTEGER(KIND=4) :: order
	INTEGER(KIND=4) :: count_peaks
	INTEGER(KIND=4), DIMENSION(2) :: g ! (n,m) vector
	CHARACTER(LEN=19), PARAMETER :: routinename = "SET_ALLOWED_PEAKS: "
	! Pointer definitions
	REAL(KIND=8) :: a, b ! axis longitude
	REAL(KIND=8) :: beta ! angle of incident beam projected on unit cell surface
	REAL(KIND=8) :: theta_in ! incidence angle measured from surface plane
	! Pointer adjudication
	a=this%surf%norm_s1
	b=this%surf%norm_s2
	beta=this%inicond%vpar_angle%getvalue()
	theta_in=this%inicond%vz_angle%getvalue()
	! FIRE IN THE HOLE>! ......................
#ifdef DEBUG 
	CALL VERBOSE_WRITE(routinename,"Starting job")
#endif
	E = (this%inicond%E_norm%getvalue())/(dsin(theta_in))**2.D0
	this%E = E
	gamma = DACOS(DOT_PRODUCT(this%surf%s1,this%surf%s2)/(a*b))
   mass=this%thispes%atomdat(1)%getmass()+this%thispes%atomdat(2)%getmass()
	pinit_par = DSQRT(2.D0*mass*(E - this%inicond%E_norm%getvalue()))
	kinit_par(1) = pinit_par*a*DCOS(beta)/(2.D0*PI)
	Kinit_par(2) = pinit_par*b*DCOS(gamma-beta)/(2.D0*PI)
	! Setting conic equation
	this%conic(1) = b**2.D0
	this%conic(2) = -2.D0*a*b*DCOS(gamma)
	this%conic(3) = a**2.D0
	this%conic(4) = 2.D0*b*(b*kinit_par(1)-a*kinit_par(2)*DCOS(gamma))
	this%conic(5) = 2.D0*a*(a*kinit_par(2)-b*kinit_par(1)*DCOS(gamma))
	this%conic(6) = -((a*b*DSIN(gamma)/(2.D0*PI))**2.D0)*2.D0*mass*this%inicond%E_norm%getvalue()
	! Debug messages:
#ifdef DEBUG
	CALL VERBOSE_WRITE(routinename,"Total energy: ", E)
	CALL VERBOSE_WRITE(routinename, "Angle between surface main axis (deg): ", gamma*180.D0/PI)
	CALL VERBOSE_WRITE(routinename, "Angle respect to surface main axis u1 (deg): ", beta*180.D0/PI)
	CALL VERBOSE_WRITE(routinename, "Incidence angle respect to surface plane (deg): ", theta_in*180.D0/PI)
	CALL VERBOSE_WRITE(routinename, "kinit_par 1: ", kinit_par(1))
	CALL VERBOSE_WRITE(routinename, "kinit_par 2: ", kinit_par(2))
	CALL VERBOSE_WRITE(routinename, "Conic A: ", this%conic(1))
	CALL VERBOSE_WRITE(routinename, "Conic B: ", this%conic(2))
	CALL VERBOSE_WRITE(routinename, "Conic C: ", this%conic(3))
	CALL VERBOSE_WRITE(routinename, "Conic D: ", this%conic(4))
	CALL VERBOSE_WRITE(routinename, "Conic E: ", this%conic(5))
	CALL VERBOSE_WRITE(routinename, "Conic F: ", this%conic(6))
#endif
	! operations
	OPEN(11,FILE="OUTallowedpeaks.out", STATUS="replace")
	WRITE(11,*) "# ----- ALLOWED PEAKS---------------------------------------"
	WRITE(11,*) "# Format: id, order, n, m, Psi(rad), Phi(rad), Theta_out(rad)   "
	WRITE(11,*) "#-----------------------------------------------------------"
	order = 0
	count_peaks = 0
	ALLOCATE(allowed(1))
	g(1) = 0
	g(2) = 0
	allowed(1) = this%evaluate_peak(g(1),g(2))
	IF(allowed(1)) THEN
		count_peaks = count_peaks + 1
		p(1) = pinit_par+(2.D0*PI/(a*b*DSIN(gamma)))*(DFLOAT(g(1))*b*DSIN(gamma-beta)+DFLOAT(g(2))*a*DSIN(beta))
		p(2) = (2.D0*PI/(a*b*DSIN(gamma)))*(DFLOAT(g(2))*a*DCOS(beta)+DFLOAT(g(1))*b*DCOS(beta+gamma))
		p(3) = DSQRT(2.D0*mass*E-p(1)**2.D0-p(2)**2.D0)
		Psi = DATAN(p(2)/p(1))
		Phi = DATAN(p(2)/p(3))
		Theta_out = DATAN(p(3)/(DSQRT(p(1)**2.D0+p(2)**2.D0)))
		WRITE(11,*) count_peaks, order, g, Psi, Phi, Theta_out
	END IF
	DEALLOCATE(allowed)
	DO
		order = order + 1
		i = 0
		ALLOCATE(allowed(8*order))
		!-----
		g(1) = order
		DO k = -order, order
			g(2) = k
			i = i + 1
			allowed(i) = this%evaluate_peak(g(1),g(2))
			IF(allowed(i)) THEN
				count_peaks = count_peaks + 1
				p(1) = pinit_par+(2.D0*PI/(a*b*DSIN(gamma)))*(DFLOAT(g(1))*b*DSIN(gamma-beta)+DFLOAT(g(2))*a*DSIN(beta))
				p(2) = (2.D0*PI/(a*b*DSIN(gamma)))*(DFLOAT(g(2))*a*DCOS(beta)-DFLOAT(g(1))*b*DCOS(gamma-beta))
				p(3) = DSQRT(2.D0*mass*E-p(1)**2.D0-p(2)**2.D0)
				Psi = DATAN(p(2)/p(1))
				Phi = DATAN(p(2)/p(3))
				Theta_out = DATAN(p(3)/(DSQRT(p(1)**2.D0+p(2)**2.D0)))
				WRITE(11,*) count_peaks, order, g, Psi, Phi, Theta_out
			END IF
		END DO
		!----------------------
		g(1) = -order
		DO k = -order, order
			g(2) = k
			i = i + 1
			allowed(i) = this%evaluate_peak(g(1),g(2))
			IF(allowed(i)) THEN
				count_peaks = count_peaks + 1
				p(1) = pinit_par+(2.D0*PI/(a*b*DSIN(gamma)))*(DFLOAT(g(1))*b*DSIN(gamma-beta)+DFLOAT(g(2))*a*DSIN(beta))
				p(2) = (2.D0*PI/(a*b*DSIN(gamma)))*(DFLOAT(g(2))*a*DCOS(beta)-DFLOAT(g(1))*b*DCOS(gamma-beta))
				p(3) = DSQRT(2.D0*mass*E-p(1)**2.D0-p(2)**2.D0)
				Psi = DATAN(p(2)/p(1))
				Phi = DATAN(p(2)/p(3))
				Theta_out = DATAN(p(3)/(DSQRT(p(1)**2.D0+p(2)**2.D0)))
				WRITE(11,*) count_peaks, order, g, Psi, Phi, Theta_out
			END IF
		END DO
		!----
		g(2) = order
		DO k = -order +1, order-1
			g(1) = k
			i = i + 1
			allowed(i) = this%evaluate_peak(g(1),g(2))
			IF(allowed(i)) THEN
				count_peaks = count_peaks + 1
				p(1) = pinit_par+(2.D0*PI/(a*b*DSIN(gamma)))*(DFLOAT(g(1))*b*DSIN(gamma-beta)+DFLOAT(g(2))*a*DSIN(beta))
				p(2) = (2.D0*PI/(a*b*DSIN(gamma)))*(DFLOAT(g(2))*a*DCOS(beta)-DFLOAT(g(1))*b*DCOS(gamma-beta))
				p(3) = DSQRT(2.D0*mass*E-p(1)**2.D0-p(2)**2.D0)
				Psi = DATAN(p(2)/p(1))
				Phi = DATAN(p(2)/p(3))
				Theta_out = DATAN(p(3)/(DSQRT(p(1)**2.D0+p(2)**2.D0)))
				WRITE(11,*) count_peaks, order, g, Psi, Phi, Theta_out
			END IF
		END DO
		!----
		g(2) = -order
		DO k = -order +1, order-1
			g(1) = k
			i = i + 1
			allowed(i) = this%evaluate_peak(g(1),g(2))
			IF(allowed(i)) THEN
				count_peaks = count_peaks + 1
				p(1) = pinit_par+(2.D0*PI/(a*b*DSIN(gamma)))*(DFLOAT(g(1))*b*DSIN(gamma-beta)+DFLOAT(g(2))*a*DSIN(beta))
				p(2) = (2.D0*PI/(a*b*DSIN(gamma)))*(DFLOAT(g(2))*a*DCOS(beta)-DFLOAT(g(1))*b*DCOS(gamma-beta))
				p(3) = DSQRT(2.D0*mass*E-p(1)**2.D0-p(2)**2.D0)
				Psi = DATAN(p(2)/p(1))
				Phi = DATAN(p(2)/p(3))
				Theta_out = DATAN(p(3)/(DSQRT(p(1)**2.D0+p(2)**2.D0)))
				WRITE(11,*) count_peaks, order, g, Psi, Phi, Theta_out
			END IF
		END DO
		DO i = 1, 8*order
			IF (allowed(i)) THEN
				need_cycle = .TRUE.
				EXIT
			END IF
			need_cycle = .FALSE.
		END DO
		DEALLOCATE(allowed)
		IF (.NOT.need_cycle) EXIT
	END DO
	REWIND (11)
	ALLOCATE(this%peaks(1:count_peaks))
	READ(11,*) ! dummy line
	READ(11,*) ! dummy line
	READ(11,*) ! dummy line
	DO i=1, count_peaks
		READ(11,*) this%peaks(i)%id, this%peaks(i)%order, this%peaks(i)%g(1), this%peaks(i)%g(2), &
			   this%peaks(i)%Psi, this%peaks(i)%Phi, this%peaks(i)%Theta_out
		this%peaks(i)%prob = 0.D0
	END DO
	CLOSE(11)
	RETURN
END SUBROUTINE SETUP_ALLOWEDPEAKSCRP6D
!####################################################################################
! SUBROUTINE: ASSIGN PEAKS TO TRAJS
!####################################################################################
! - GENERATE_TRAJS_ATOMS and SET_Allowed_peaksCRP6D should have been executed before
! - At the moment only works with C4v cells
!------------------------------------------------------------------------------------
SUBROUTINE ASSIGN_PEAKS_TO_TRAJS_ALLOWEDPEAKSCRP6D(this)
	IMPLICIT NONE
	! I/O variables
	CLASS(Allowed_peaksCRP6D),INTENT(INOUT):: this
	! Local variables
	INTEGER(KIND=4) :: lines
	INTEGER(KIND=4) :: id
	INTEGER(KIND=4) :: dummy_int
	INTEGER(KIND=4) :: i,j ! counters
	INTEGER(KIND=4), DIMENSION(2) :: g
	REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: peaks_prob
	REAL(KIND=8), DIMENSION(2,2) :: to_rec_space
	REAL(KIND=8) :: dummy_real
   REAL(KIND=8),DIMENSION(9) :: dummy
	REAL(KIND=8) :: gamma ! angle between unit cell surface vectors
	REAL(KIND=8), DIMENSION(2) :: p ! final momentum
	REAL(KIND=8), DIMENSION(2) :: dp ! variation of momentum
	REAL(KIND=8), DIMENSION(2) :: dk ! variation of momentum in rec. space coord
	CHARACTER(LEN=10) :: stat
	CHARACTER(LEN=23), PARAMETER :: routinename = "ASSIGN_PEAKS_TO_TRAJS: "
   INTEGER(KIND=4) :: ioerr
	! Pointer definitions
	REAL(KIND=8) :: beta ! angle between incident parallel momentum respect to u1 (surface vector)
	REAL(KIND=8) :: a, b ! length of surface main axis
	! Pointers assignation
	beta=this%inicond%vpar_angle%getvalue()
	a=this%surf%norm_s1
	b=this%surf%norm_s2
	! RUN SECTION -------------------------
#ifdef DEBUG
	CALL VERBOSE_WRITE(routinename, "Starting job")
#endif
	gamma = DACOS(DOT_PRODUCT(this%surf%s1,this%surf%s2)/(a*b))
	to_rec_space(1,1) = a/(2.D0*PI)
	to_rec_space(1,2) = 0.D0
	to_rec_space(2,1) = b*DCOS(gamma)/(2.D0*PI)
	to_rec_space(2,2) = b*DSIN(gamma)/(2.D0*PI)
	!---------
   OPEN(12,FILE="OUTmappingpeaks.out",STATUS="replace")
   WRITE(12,*) "#·-------- MAPPING TRAJECTORIES WITH DIFFRACTION PEAKS ----------"
   WRITE(12,*) "# Format: Trajectory id, ----> , Allowed peak id"
   WRITE(12,*) "#----------------------------------------------------------------"
   OPEN(11,FILE="OUTdynamics6d.MOLEC.out",STATUS="old")
   READ(11,*) ! dummy line
   READ(11,*) ! dummy line
   READ(11,*) ! dummy line
   i=0
   DO 
      i=i+1
      READ(11,*,IOSTAT=ioerr) id,stat,dummy_int,dummy_int,dummy(:),p(:)
      SELECT CASE(ioerr==0)
         CASE(.TRUE.)
            ! do nothing
         CASE(.FALSE.)
            WRITE(*,*) "ASSIGN_PEAKS_TO_TRAJS: EOF reached"
            EXIT
      END SELECT
      IF (stat.EQ."Scattered") THEN
         dp(1) = p(1)-this%inicond%trajs(id)%init_p(1)
         dp(2) = p(2)-this%inicond%trajs(id)%init_p(2)
         dk = MATMUL(to_rec_space,dp)
         g(1) = NINT(dk(1))
         g(2) = NINT(dk(2))
         DO j= 1,SIZE(this%peaks)
            IF ((this%peaks(j)%g(1).EQ.g(1)).AND.(this%peaks(j)%g(2).EQ.g(2))) THEN
               WRITE(12,*) id," ----> ",this%peaks(j)%id
               EXIT
            END IF
         END DO
      END IF
   END DO
         lines=i-1
	REWIND(12)
	READ(12,*) ! dummy line
	READ(12,*) ! dummy line
	READ(12,*) ! dummy line
	DO i=1,lines
		READ(12,*) dummy_int, stat, id
		this%peaks(id)%prob = this%peaks(id)%prob + 1.D0/lines
	END DO
	CLOSE(12)
	CLOSE(11)
	OPEN(13,FILE="OUTseenpeaks.out",STATUS = "replace") ! re-write allowed peaks file with probabilities printed
	WRITE(13,*) "# ----- ALLOWED PEAKS---------------------------------------"
	WRITE(13,*) "# Format: id, order, n, m, Psi(rad), Phi(rad), Theta_out(rad), Prob"
	WRITE(13,*) "#-----------------------------------------------------------"
	DO i=1, SIZE(this%peaks)
		IF (this%peaks(i)%prob.NE.0.D0) THEN
			WRITE(13,*) this%peaks(i)%id,this%peaks(i)%order,this%peaks(i)%g(1),this%peaks(i)%g(2), &
				    this%peaks(i)%Psi,this%peaks(i)%Phi,this%peaks(i)%Theta_out,this%peaks(i)%prob
		END IF
	END DO
	CLOSE(13)
	RETURN
END SUBROUTINE ASSIGN_PEAKS_TO_TRAJS_ALLOWEDPEAKSCRP6D
!####################################################################################
! FUNCTION: EVALUATE_PEAK ###########################################################
!####################################################################################
! - TRUE if A*(n^2) + B*nm + C*(m^2) + D*n + E*m < -F
! - FALSE otherwise
! - in_n and in_m are integer numbers
!------------------------------------------------------------------------------------
LOGICAL FUNCTION evaluate_peak_ALLOWEDPEAKSCRP6D(this,in_n, in_m)
	IMPLICIT NONE
	! I/O variables
	CLASS(Allowed_peaksCRP6D),TARGET,INTENT(IN) :: this
	INTEGER(KIND=4), INTENT(IN) :: in_n, in_m
	! Local variables
	REAL(KIND=8) :: left_term
	REAL(KIND=8), POINTER :: A, B, C, D, E, F
	REAL(KIND=8) :: n, m
	! HEY, HO ! LET'S GO!
	! Pointers & stuff
	A => this%conic(1)
	B => this%conic(2)
	C => this%conic(3)
	D => this%conic(4)
	E => this%conic(5)
	F => this%conic(6)
	n = DFLOAT(in_n)
	m = DFLOAT(in_m)
	! left-term of the equation (A, B, C, D, E)
	left_term = A*(n**2.D0) + B*n*m +C*(m**2.D0) + D*n + E*m 
	! Check value
   SELECT CASE(left_term<-F)
      CASE(.TRUE.)
         evaluate_peak_ALLOWEDPEAKSCRP6D = .TRUE.
      CASE(.FALSE.)
         evaluate_peak_ALLOWEDPEAKSCRP6D = .FALSE.
   END SELECT
	RETURN
END FUNCTION evaluate_peak_ALLOWEDPEAKSCRP6D
!###########################################################################3########
! SUBROUTINE: PRINT_XY_EXIT_ANGLES 
!###########################################################################3########
! - Reads data form "input_file" (should be dynamics-like output) and creates file "xy_exit_angle.out" with
!   information about exit angles in XY plane (taken from momenta information)
! - Only trajectories with "Scattered" status will be taken into account
!------------------------------------------------------------------------------------
SUBROUTINE PRINT_XY_EXIT_ANGLES_ALLOWEDPEAKSCRP6D(this,input_file)
	IMPLICIT NONE
	! I/O Variables
   CLASS(Allowed_peaksCRP6D),INTENT(IN):: this
	CHARACTER(LEN=*), INTENT(IN) :: input_file
	! Local variables
	INTEGER(KIND=4) :: lines
	INTEGER(KIND=4) :: dummy_int
	REAL(KIND=8) :: dummy_real
	INTEGER(KIND=4) :: i ! counters
	INTEGER(KIND=4) :: traj_id
	CHARACTER(LEN=10) :: stat
	REAL(KIND=8), DIMENSION(2) :: p
	REAL(KIND=8) :: angle
   INTEGER(KIND=4) :: ioerr
   ! RUN !! --------------------------
   OPEN(12,FILE="xy_exit_angles.out",STATUS="replace")
   WRITE(12,*) "# XY EXIT ANGLES ----------------------------------------------"
   WRITE(12,*) "# Format: traj id, angle (radians)"
   WRITE(12,*) "#--------------------------------------------------------------"
	OPEN(11,FILE=input_file,STATUS="old")
	READ(11,*) ! Dummy line
	READ(11,*) ! Dummy line
	READ(11,*) ! Dummy line
	DO i=1,lines
		READ(11,*,iostat=ioerr) traj_id, stat, dummy_int, dummy_int, dummy_real,&
         dummy_real, dummy_real, dummy_real, dummy_real, p(1), p(2)
      SELECT CASE (ioerr==0)
         CASE(.TRUE.)
            ! do nothing
         CASE(.FALSE.)
            EXIT ! break cycle
         END SELECT
		IF(stat.EQ."Scattered") THEN
			angle = DATAN(p(2)/p(1))
			WRITE(12,*) traj_id, angle
		END IF
	END DO
	CLOSE(11)
	CLOSE(12)
	RETURN
END SUBROUTINE PRINT_XY_EXIT_ANGLES_ALLOWEDPEAKSCRP6D
!
END MODULE DIFFRACTIONCRP6D_MOD
