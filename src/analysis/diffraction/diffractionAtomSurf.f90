MODULE DIFFRACTIONATOMSURF_MOD
use SYSTEM_MOD
use INITATOMSURF_MOD, only: InitAtomSurf
use CRP3D_MOD, only: Crp3d
#ifdef DEBUG
use DEBUG_MOD, only: verbose_write,debug_write
#endif
IMPLICIT NONE
!======================================================
! Peak derived data
!---------------------
TYPE :: PeakAtomSurf
   PRIVATE
   INTEGER(KIND=4) :: id ! identification number
   INTEGER(KIND=4) :: order ! order of the peak
   INTEGER(KIND=4), DIMENSION(2) :: g ! coordinates in the reciprocal space lattice
	REAL(KIND=8) :: Psi ! azimuthal exit angle
	REAL(KIND=8) :: Phi ! deflection angle respect to incidence plane
	REAL(KIND=8) :: Theta_out ! deflection angle respect to surface plane
	REAL(KIND=8) :: Prob ! probability
END TYPE PeakAtomSurf
!======================================================
! Allowed_peaksAtomSurf derived data
!----------------------------
TYPE :: Allowed_peaksAtomSurf
   PRIVATE
   TYPE(InitAtomSurf):: inicond
   TYPE(CRP3D):: thispes
   REAL(KIND=8):: E
   REAL(KIND=8),DIMENSION(6):: conic
   TYPE(PeakAtomSurf),DIMENSION(:),ALLOCATABLE:: peaks
   CONTAINS
      PROCEDURE,PUBLIC:: INITIALIZE => INITIALIZE_ALLOWEDPEAKSATOMSURF
      PROCEDURE,PUBLIC:: SETUP => SETUP_ALLOWEDPEAKSATOMSURF
      PROCEDURE,PUBLIC:: ASSIGN_PEAKS => ASSIGN_PEAKS_TO_TRAJS_ALLOWEDPEAKSATOMSURF
      PROCEDURE,PUBLIC:: PRINT_LABMOMENTA_AND_ANGLES => PRINT_LABMOMENTA_AND_ANGLES_ALLOWEDPEAKSATOMSURF
      PROCEDURE,PUBLIC:: evaluate_peak => evaluate_peak_ALLOWEDPEAKSATOMSURF
END TYPE Allowed_peaksAtomSurf
!=======================================================
CONTAINS
SUBROUTINE INITIALIZE_ALLOWEDPEAKSATOMSURF(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Allowed_peaksAtomSurf),INTENT(OUT):: this
   ! Local variables
   ! Run section
   CALL this%thispes%INITIALIZE()
   CALL this%inicond%INITIALIZE()
   CALL this%inicond%GENERATE_TRAJS(this%thispes)
   RETURN
END SUBROUTINE 
!######################################################
!# SUBROUTINE : SETUP_ALLOWEDPEAKSATOMSURF #####################
!######################################################
! Just for a square primitive cell
!------------------------------------------------------
SUBROUTINE SETUP_ALLOWEDPEAKSATOMSURF(this)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Allowed_peaksAtomSurf),INTENT(INOUT):: this
   ! Local variables
	REAL(KIND=8),DIMENSION(2):: kinit_par
	REAL(KIND=8),DIMENSION(3):: p ! momentum
	LOGICAL,DIMENSION(:),ALLOCATABLE:: allowed
	LOGICAL:: need_cycle
	REAL(KIND=8):: pinit_par
	REAL(KIND=8):: E ! Total energy
	REAL(KIND=8):: gamma ! angle between surface main vectors.
	REAL(KIND=8):: Psi ! azimuthal exit angle
	REAL(KIND=8):: Phi ! deflection angle respect to incidence plane
	REAL(KIND=8):: Theta_out ! deflection angle respect to surface plane
   REAL(KIND=8):: mass
   INTEGER(KIND=4):: i, k ! Counters
   INTEGER(KIND=4):: order, realorder
   INTEGER(KIND=4):: count_peaks
   INTEGER(KIND=4),DIMENSION(2):: g ! (n,m) vector
   CHARACTER(LEN=*),PARAMETER:: routinename = "SET_Allowed_peaksAtomSurf: "
   REAL(KIND=8):: a, b ! axis longitude
   REAL(KIND=8):: beta ! angle of incident beam projected on unit cell surface
   REAL(KIND=8):: theta_in ! incidence angle measured from surface plane
	! Pointer adjudication
	a=system_surface%norm_s1
	b=system_surface%norm_s2
	beta=this%inicond%vpar_angle%getvalue()
	theta_in=this%inicond%vz_angle%getvalue()
	! FIRE IN THE HOLE>! ......................
#ifdef DEBUG 
	CALL VERBOSE_WRITE(routinename,"Starting job")
#endif
	E = this%inicond%E_norm%getvalue()/(DSIN(theta_in)**2.D0)
	this%E = E
	gamma = DACOS(DOT_PRODUCT(system_surface%s1,system_surface%s2)/(a*b))
   mass=system_mass(1)
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
	OPEN(11,FILE="OUTANA3Dallowedpeaks.out",STATUS="replace")
	WRITE(11,*) "# ***** ALLOWED PEAKS *****"
	WRITE(11,*) "# Format: id/order,n,m/Azimuthal,Polar,Deflection(rad)"
	WRITE(11,*) "# -----------------------------------------------------------"
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
      ! Determine Azimuthal angle. Avoid indetermination
      SELECT CASE(g(1)==0 .AND. g(2)==0)
         CASE(.TRUE.)
            Psi = 0.D0
         CASE(.FALSE.)
            Psi = DATAN(p(2)/p(1))
      END SELECT
		Phi = DATAN(p(2)/p(3))
		Theta_out = DATAN(p(3)/(DSQRT(p(1)**2.D0+p(2)**2.D0)))
      CALL SET_REALORDER_CUAD(order,g,realorder)
		WRITE(11,*) count_peaks, realorder, g, Psi, Phi, Theta_out
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
            CALL SET_REALORDER_CUAD(order,g,realorder)
				WRITE(11,*) count_peaks,realorder, g, Psi, Phi, Theta_out
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
            CALL SET_REALORDER_CUAD(order,g,realorder)
				WRITE(11,*) count_peaks,realorder, g, Psi, Phi, Theta_out
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
            CALL SET_REALORDER_CUAD(order,g,realorder)
				WRITE(11,*) count_peaks,realorder, g, Psi, Phi, Theta_out
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
            CALL SET_REALORDER_CUAD(order,g,realorder)
				WRITE(11,*) count_peaks,realorder, g, Psi, Phi, Theta_out
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
END SUBROUTINE SETUP_ALLOWEDPEAKSATOMSURF
!####################################################################################
! SUBROUTINE: ASSIGN PEAKS TO TRAJS
!####################################################################################
! - GENERATE_TRAJS_ATOMS and SET_Allowed_peaksAtomSurf should have been executed before
! - At the moment only works with C4v cells
!------------------------------------------------------------------------------------
SUBROUTINE ASSIGN_PEAKS_TO_TRAJS_ALLOWEDPEAKSATOMSURF(this)
	IMPLICIT NONE
	! I/O variables
	CLASS(Allowed_peaksAtomSurf),INTENT(INOUT):: this
	! Local variables
	INTEGER(KIND=4) :: tottrajs ! total number of trajs
   INTEGER(KIND=4) :: totscatt ! total number of scattered trajs
	INTEGER(KIND=4) :: id
	INTEGER(KIND=4) :: dummy_int
	INTEGER(KIND=4) :: i,j ! counters
	INTEGER(KIND=4), DIMENSION(2) :: g
	REAL(KIND=8), DIMENSION(2,2) :: to_rec_space
	REAL(KIND=8),dimension(8):: dummy_real
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
	a=system_surface%norm_s1
	b=system_surface%norm_s2
	! RUN SECTION -------------------------
#ifdef DEBUG
	CALL VERBOSE_WRITE(routinename, "Starting job")
#endif
	gamma = DACOS(DOT_PRODUCT(system_surface%s1,system_surface%s2)/(a*b))
	to_rec_space(1,1) = a/(2.D0*PI)
	to_rec_space(1,2) = 0.D0
	to_rec_space(2,1) = b*DCOS(gamma)/(2.D0*PI)
	to_rec_space(2,2) = b*DSIN(gamma)/(2.D0*PI)
	!---------
   OPEN(12,FILE="OUTANA3Dmappingpeaks.out",STATUS="replace")
   WRITE(12,*) "# ***** MAPPING TRAJECTORIES WITH DIFFRACTION PEAKS *****"
   WRITE(12,*) "# Format: traj id/ ----> /peak id"
   WRITE(12,*) "# ----------------------------------------------------------------"
   OPEN(11,FILE="OUTDYN3Dscattered.out",STATUS="old")
   call skipHeaderFromFile("OUTDYN3Dscattered.out")
   i=0
   DO 
      i=i+1
      READ(11,*,IOSTAT=ioerr) id,stat,dummy_int,dummy_int,dummy_real(:),p
      SELECT CASE(ioerr==0)
         CASE(.TRUE.)
            ! do nothing
         CASE(.FALSE.)
            IF(ioerr/=-1) WRITE(*,*) routinename//'Unexpected error in scattered trajs file. Err Code: ',ioerr
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
   totscatt=i-1
   tottrajs=this%inicond%ntraj-this%inicond%nstart+1
   WRITE(*,*) "==========================================================="
   WRITE(*,*) "ASSIGN PEAKS TO TRAJS: total trajs: ",tottrajs
   WRITE(*,*) "ASSIGN PEAKS TO TRAJS: scattered trajs: ",totscatt
   WRITE(*,*) "ASSIGN PEAKS TO TRAJS: probability: ",totscatt/tottrajs
   WRITE(*,*) "==========================================================="
	REWIND(12)
	READ(12,*) ! dummy line
	READ(12,*) ! dummy line
	READ(12,*) ! dummy line
	DO i=1,totscatt
		READ(12,*) dummy_int,stat,id
		this%peaks(id)%prob = this%peaks(id)%prob + 1.D0/tottrajs
	END DO
	CLOSE(12)
	CLOSE(11)
	OPEN(13,FILE="OUTANA3Dseenpeaks.out",STATUS = "replace") ! re-write allowed peaks file with probabilities printed
	WRITE(13,*) "# ***** ALLOWED PEAKS *****"
	WRITE(13,*) "# Format: id/order,n,m/Azimuthal,Polar,Deflection(rad)/Prob"
	WRITE(13,*) "# -----------------------------------------------------------"
	DO i=1, SIZE(this%peaks)
		IF (this%peaks(i)%prob.NE.0.D0) THEN
			WRITE(13,*) this%peaks(i)%id,this%peaks(i)%order,this%peaks(i)%g(1),this%peaks(i)%g(2),&
				    this%peaks(i)%Psi,this%peaks(i)%Phi,this%peaks(i)%Theta_out,this%peaks(i)%prob
		END IF
	END DO
	CLOSE(13)
	RETURN
END SUBROUTINE ASSIGN_PEAKS_TO_TRAJS_ALLOWEDPEAKSATOMSURF
!####################################################################################
! FUNCTION: EVALUATE_PEAK ###########################################################
!####################################################################################
!> @brief
! - TRUE if A*(n^2) + B*nm + C*(m^2) + D*n + E*m < -F
! - FALSE otherwise
! - in_n and in_m are integer numbers
!------------------------------------------------------------------------------------
function evaluate_peak_ALLOWEDPEAKSATOMSURF(this,in_n, in_m) result(isAllowed)
   implicit none
   ! I/O variables
   class(Allowed_peaksAtomSurf),intent(in):: this
   integer(kind=4),intent(in):: in_n, in_m
   ! Function dummy variable
   logical:: isAllowed
   ! Local variables
   real(kind=8):: left_term
   real(kind=8):: n, m
   ! HEY, HO ! LET'S GO!
   n = DFLOAT(in_n)
   m = DFLOAT(in_m)
   ! left-term of the equation (A, B, C, D, E)
   left_term=this%conic(1)*(n**2.D0)+this%conic(2)*n*m+this%conic(3)*(m**2.D0)+this%conic(4)*n+this%conic(5)*m 
   ! Check value
   isAllowed=(left_term<this%conic(6))
   return
   end function evaluate_peak_ALLOWEDPEAKSATOMSURF
!###########################################################################3########
! SUBROUTINE: PRINT_XY_EXIT_ANGLES 
!###########################################################################3########
! - Reads data form "input_file" (should be dynamics-like output) and creates file "xy_exit_angle.out" with
!   information about exit angles in XY plane (taken from momenta information)
! - Only trajectories with "Scattered" status will be taken into account
!------------------------------------------------------------------------------------
SUBROUTINE PRINT_LABMOMENTA_AND_ANGLES_ALLOWEDPEAKSATOMSURF(this)
   ! Initial declarations
   implicit none
   ! I/O Variables
   class(Allowed_peaksAtomSurf),intent(in):: this
   ! Local variables
   integer(kind=4):: dummy_int
   real(kind=8),dimension(8):: dummy_real
   integer(kind=4):: traj_id
   character(len=10):: stat
   real(kind=8),dimension(3):: p,plab
   integer(kind=4) :: ioerr
   real(kind=8):: psi,Theta,thetaout,beta
   real(kind=8),dimension(2,2):: mtrx
   ! Opened units
   integer(kind=4),parameter:: ru=11
   integer(kind=4),parameter:: wu=12
   ! RUN !! --------------------------
   beta=this%inicond%vpar_angle%getvalue()
   mtrx(1,:)=[dcos(beta),dsin(beta)]
   mtrx(2,:)=[-dsin(beta),dcos(beta)]
   open(unit=wu,file="OUTANA3Dfinalpandangles.out",status="replace",action='write')
   write(wu,*) "# ***** FINAL MOMENTA AND EXIT ANGLES *****"
   write(wu,*) "# Format: id/Ppara,Pperp,Pnormal(a.u.)/Azimuthal,Polar,Deflection(rad)"
   write(wu,*) "# ---------------------------------------------------------------------"
   open(ru,file="OUTDYN3Dscattered.out",status="old",action='read')
   call skipHeaderFromFile("OUTDYN3Dscattered.out")
   ioErr=0
   do while(ioErr==0)
      read(ru,*,iostat=ioerr) traj_id,stat,dummy_int,dummy_int,dummy_real(:),p(:)
      select case( stat )
      case('Scattered')
         plab(1:2)=matmul(mtrx,p(1:2))
         plab(3)=p(3)
         psi = datan(plab(2)/plab(1))
         Theta = datan(plab(2)/plab(3))
         thetaout=datan(plab(3)/dsqrt(plab(1)**2.D0+plab(2)**2.D0))
         write(wu,*) traj_id,plab(:),psi,Theta,thetaout
      case default
         ! do nothing
      end select
   enddo
   close(ru)
   close(wu)
   return
END SUBROUTINE PRINT_LABMOMENTA_AND_ANGLES_ALLOWEDPEAKSATOMSURF
!###########################################################
!# SUBROUTINE: SET_REALORDER_CUAD
!###########################################################
!> @brief
!! Given environment order and peak labels, sets real diffraction
!! order.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Nov/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_REALORDER_CUAD(order,g,realorder)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: order
   INTEGER(KIND=4),DIMENSION(2),INTENT(IN) :: g
   INTEGER(KIND=4),INTENT(OUT) :: realorder
   ! Run section
   SELECT CASE(g(1)==0 .AND. g(2)==0)
      CASE(.TRUE.)
         realorder=order
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(g(1)==0 .OR. g(2)==0)
      CASE(.TRUE.)
         realorder=order*(order+1)/2
         RETURN   
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(abs(g(1))==abs(g(2)))
      CASE(.TRUE.)
         realorder=((order+1)*(order+2)/2)-1
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(abs(g(1))==order)
      CASE(.TRUE.)
         realorder=(order*(order+1)/2)+abs(g(2))
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(abs(g(2))==order)
      CASE(.TRUE.)
         realorder=(order*(order+1)/2)+abs(g(1))
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE SET_REALORDER_CUAD
!
END MODULE DIFFRACTIONATOMSURF_MOD
