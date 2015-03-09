MODULE DIFFRACTIONCRP6D_MOD
use SYSTEM_MOD
use INITDIATOMIC_MOD, only: InitDiatomic
use CRP6D_MOD, only: CRP6D
#ifdef DEBUG
use DEBUG_MOD, only: VERBOSE_WRITE, DEBUG_WRITE
#endif
IMPLICIT NONE
!////////////////////////////////////////////////////////
! TYPE: SubPeakCRP6D
!////////////////////////////////////////////////////////
!> @param robivrState(:) - integer(kind=4): V,J,mJ state
!> @param phiOut - real(kind=8): deflection angle respect to the perpendicular plane to the surface
!> @param thetaOut - real(kind=8): deflection angle respect to surface plane
!> @param prob - real(kind=8): Probability
!> @param dE - real(kind=8): internal energy exchange after data binning
!------------------------------------------------------
type:: SubPeakCRP6D
   private
	real(kind=8):: phiOut
	real(kind=8):: thetaOut
	integer(kind=4),dimension(3):: rovibrState
	real(kind=8):: prob
	real(kind=8):: dE
end type SubPeakCRP6D
!======================================================
! TYPE: PeakCRP6D
!------------------------------------------------------
!> @brief
!! All information that defines a fiffraction peak
!> @param id - integer(kind=4): Identification number
!> @param envOrder - integer(kind=4): environment order
!> @param diffOrder - integer(kind=4): diffraction order based on momentum exchange
!> @param dkxy - real(kind=8): momentum exchange (in XY plane)
!> @param g - integer(kind=4),dimension(2): diffraction state
!> @param psiOut - real(kind=8): azimuthal exit angle
!------------------------------------------------------
TYPE :: PeakCRP6D
   PRIVATE
   INTEGER(KIND=4):: id
   integer(kind=4):: envOrder
   integer(kind=4):: diffOrder
   real(kind=8):: dkx
   real(kind=8):: dky
   real(kind=8):: dkxy
   INTEGER(KIND=4),DIMENSION(2):: g
	REAL(KIND=8):: psiOut
	type(SubPeakCRP6D),dimension(:),allocatable:: subPeaks
END TYPE PeakCRP6D
!======================================================
! Allowed_peaksCRP6D derived data
!----------------------------
TYPE :: Allowed_peaksCRP6D
   PRIVATE
   TYPE(Initdiatomic):: inicond
   TYPE(CRP6D):: thispes
   REAL(KIND=8):: E
   REAL(KIND=8),DIMENSION(6):: conic
   TYPE(PeakCRP6D),DIMENSION(:),ALLOCATABLE:: peaks
   character(len=24):: fileNameAllowed='OUTANA6DallowedPeaks.out'
   character(len=21):: fileNameSeen='OUTANA6DseenPeaks.out'
   CONTAINS
      ! private tools section
      procedure,private:: getPeakId => getPeakId_ALLOWEDPEAKSCRP6D
      procedure,private:: createNewPeak => createNewPeak_ALLOWEDPEAKSCRP6D
      procedure,private:: addProbToSubPeak => addProbToSubPeak_ALLOWEDPEAKSCRP6D
      ! public tools section
      PROCEDURE,PUBLIC:: INITIALIZE => INITIALIZE_ALLOWEDPEAKSCRP6D
      PROCEDURE,PUBLIC:: SETUP => SETUP_ALLOWEDPEAKSCRP6D
      procedure,public:: printAllowedPeaks => printAllowedPeaks_ALLOWEDPEAKSCRP6D
      procedure,public:: printSeenPeaks => printSeenPeaks_ALLOWEDPEAKSCRP6D
      procedure,public:: sortByDiffOrder => sortByDiffOrder_ALLOWEDPEAKSCRP6D
      PROCEDURE,PUBLIC:: ASSIGN_PEAKS => ASSIGN_PEAKS_TO_TRAJS_ALLOWEDPEAKSCRP6D
      PROCEDURE,PUBLIC:: PRINT_LABMOMENTA_AND_ANGLES => PRINT_LABMOMENTA_AND_ANGLES_ALLOWEDPEAKSCRP6D
      PROCEDURE,PUBLIC:: evaluate_peak => evaluate_peak_ALLOWEDPEAKSCRP6D
END TYPE Allowed_peaksCRP6D
!=======================================================
CONTAINS
SUBROUTINE INITIALIZE_ALLOWEDPEAKSCRP6D(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Allowed_peaksCRP6D),INTENT(OUT):: this
   ! Local variables
   LOGICAL :: exists
   ! Run section
   CALL this%thispes%INITIALIZE()
   CALL this%inicond%INITIALIZE()
   CALL this%inicond%GENERATE_TRAJS(this%thispes)
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
   CLASS(Allowed_peaksCRP6D),INTENT(INOUT):: this
   ! Local variables
	REAL(KIND=8),DIMENSION(2):: kinit_par
	REAL(KIND=8),DIMENSION(3):: p ! momentum
	LOGICAL,DIMENSION(:),ALLOCATABLE:: allowed
	LOGICAL:: need_cycle
	REAL(KIND=8):: pinit_par
	REAL(KIND=8):: E ! Total energy
	REAL(KIND=8):: gama ! angle between surface main vectors.
	REAL(KIND=8):: psiOut ! azimuthal exit angle
	REAL(KIND=8):: phiOut ! deflection angle respect to incidence plane
	REAL(KIND=8):: thetaOut ! deflection angle respect to surface plane
   REAL(KIND=8):: mass
	INTEGER(KIND=4):: i,k ! Counters
	INTEGER(KIND=4):: order,diffOrder
	INTEGER(KIND=4):: count_peaks
	INTEGER(KIND=4),DIMENSION(2):: g ! (n,m) vector
   REAL(KIND=8):: a, b ! axis longitude
   REAL(KIND=8):: beta ! angle of incident beam projected on unit cell surface
   REAL(KIND=8):: theta_in ! incidence angle measured from surface plane
   ! Parameters
   character(len=*),parameter:: routinename = "SET_ALLOWED_PEAKS: "
   integer(kind=4),parameter:: wuAllowed=11
   character(len=*),parameter:: formatAllowed='(I5,1X,4(I5,1X),2(F10.5,1X))'
	! FIRE IN THE HOLE>! ......................
	a=system_surface%norm_s1
	b=system_surface%norm_s2
	beta=this%inicond%vpar_angle%getvalue()
	theta_in=this%inicond%vz_angle%getvalue()
#ifdef DEBUG 
	CALL VERBOSE_WRITE(routinename,"Starting job")
#endif
	E = (this%inicond%E_norm%getvalue())/(dsin(theta_in))**2.D0
	this%E = E
	gama = DACOS(DOT_PRODUCT(system_surface%s1,system_surface%s2)/(a*b))
   mass=sum(system_mass(1:2))
	pinit_par = DSQRT(2.D0*mass*(E - this%inicond%E_norm%getvalue()))
	kinit_par(1) = pinit_par*a*DCOS(beta)/(2.D0*PI)
	Kinit_par(2) = pinit_par*b*DCOS(gama-beta)/(2.D0*PI)
	! Setting conic equation
	this%conic(1) = b**2.D0
	this%conic(2) = -2.D0*a*b*DCOS(gama)
	this%conic(3) = a**2.D0
	this%conic(4) = 2.D0*b*(b*kinit_par(1)-a*kinit_par(2)*DCOS(gama))
	this%conic(5) = 2.D0*a*(a*kinit_par(2)-b*kinit_par(1)*DCOS(gama))
	this%conic(6) = -((a*b*DSIN(gama)/(2.D0*PI))**2.D0)*2.D0*mass*this%inicond%E_norm%getvalue()
	! Debug messages:
#ifdef DEBUG
	CALL VERBOSE_WRITE(routinename,"Total energy: ", E)
	CALL VERBOSE_WRITE(routinename, "Angle between surface main axis (deg): ", gama*180.D0/PI)
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
	order = 0
	count_peaks = 0
	ALLOCATE(allowed(1))
	g(1) = 0
	g(2) = 0
	allowed(1) = this%evaluate_peak(g(1),g(2))
	IF(allowed(1)) THEN
		count_peaks = count_peaks + 1
		call this%createNewPeak(order,g)
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
				call this%createNewPeak(order,g)
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
				call this%createNewPeak(order,g)
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
				call this%createNewPeak(order,g)
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
				call this%createNewPeak(order,g)
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
	RETURN
END SUBROUTINE SETUP_ALLOWEDPEAKSCRP6D
!####################################################################################
! SUBROUTINE: ASSIGN PEAKS TO TRAJS
!####################################################################################
! - GENERATE_TRAJS_ATOMS and SET_Allowed_peaksCRP6D should have been executed before
! - At the moment only works with C4v cells
!> @param[in] dJ - integer(kind=2): Variation of J
!> @param[in] Ed - real(kind=8): Dissociation energy from Morse potential fit of vacuum potential
!> @param[in] A  - real(kind=8): Parameter of Morse potential fit
!------------------------------------------------------------------------------------
SUBROUTINE ASSIGN_PEAKS_TO_TRAJS_ALLOWEDPEAKSCRP6D(this)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Allowed_peaksCRP6D),INTENT(INOUT):: this
   ! Local variables
   INTEGER(KIND=4):: totscatt
   INTEGER(KIND=4):: tottrajs
   INTEGER(KIND=4):: allowedScatt
   INTEGER(KIND=4):: id
   INTEGER(KIND=4):: dummy_int
   INTEGER(KIND=4):: i,j ! counters
   INTEGER(KIND=4),DIMENSION(2):: g
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: peaks_prob
	REAL(KIND=8),DIMENSION(2,2):: to_rec_space
	REAL(KIND=8):: dummy_real
   REAL(KIND=8),DIMENSION(2):: dummy
	REAL(KIND=8):: gama ! angle between unit cell surface vectors
	REAL(KIND=8),DIMENSION(6):: p,r ! final momentum and position
	REAL(KIND=8),DIMENSION(2):: dp ! variation of momentum
	REAL(KIND=8),DIMENSION(2):: dk ! variation of momentum in rec. space coord
	real(kind=8)::finalJ,finalV,L2,Etot,Ecm,Erot,Evibr,masa,mu
   real(kind=8):: morseEd,morseWidth
   integer(kind=4):: dJ
   CHARACTER(LEN=10):: stat
   INTEGER(KIND=4):: ioerr
   LOGICAL:: isAllowed
   integer(kind=4),dimension(3):: rovibrState
   integer(kind=4):: peakId
	REAL(KIND=8):: beta ! angle between incident parallel momentum respect to u1 (surface vector)
	REAL(KIND=8):: a,b ! length of surface main axis
   ! Some parameters
   character(len=*),parameter:: routinename = "ASSIGN_PEAKS_TO_TRAJS: "
   character(len=*),parameter:: formatMap='(I6," ---> ",I6,I5,I5,3(I4))' ! Id/--->/peak id/n/m/V/J/mJ
   character(len=*),parameter:: formatUnmap='(I6," ---> ",2(I5))'        ! Id/--->/n/m
   ! Read/write units
   INTEGER(KIND=4),PARAMETER:: rwuMap=12
   INTEGER(KIND=4),PARAMETER:: wuUnmap=13
   INTEGER(KIND=4),PARAMETER:: ruScatt=14
   ! Pointers assignation
   beta=this%inicond%vpar_angle%getvalue()
   a=system_surface%norm_s1
   b=system_surface%norm_s2
   masa=sum(system_mass(:))
   mu=product(system_mass(:))/masa
   ! RUN SECTION -------------------------
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename, "Starting job")
#endif
   gama = DACOS(DOT_PRODUCT(system_surface%s1,system_surface%s2)/(a*b))
   to_rec_space(1,1) = a/(2.D0*PI)
   to_rec_space(1,2) = 0.D0
   to_rec_space(2,1) = b*DCOS(gama)/(2.D0*PI)
   to_rec_space(2,2) = b*DSIN(gama)/(2.D0*PI)
   ! binning parameters (only Morse implemented)
   morseEd=system_binningParam(1)
   morseWidth=system_binningParam(3)
   dJ=system_binningdJ
   !---------
   OPEN(unit=rwuMap,file="OUTANA6Dmappingpeaks.out",status="replace",action='readwrite')
   WRITE(rwuMap,*) "# ***** MAPPING TRAJECTORIES WITH DIFFRACTION PEAKS *****"
   WRITE(rwuMap,*) "# Format: traj id/ ----> /peak id/n,m/V,J,mJ"
   WRITE(rwuMap,*) "# ----------------------------------------------------------------"
   OPEN(unit=wuUnmap,file='OUTANA6Dunmappedtrajs.out',status='replace',action='write')
   WRITE(wuUnmap,*) '# ***** LIST OF UNMAPPED TRAJS *****'
   WRITE(wuUnmap,*) "# Format: id/n,m/dkx,dky"
   WRITE(wuUnmap,*) "# ----------------------------------------------------------------"
   OPEN(unit=ruScatt,file="OUTDYN6Dscattered.out",status="old",action='read')
   call skipHeaderFromFile(unit=ruScatt)
   i=0
   allowedScatt=0
   DO 
      i=i+1
      READ(ruScatt,*,IOSTAT=ioerr) id,stat,dummy_int,dummy_int,Etot,dummy(:),r(:),p(:)
      SELECT CASE(ioerr==0)
         CASE(.TRUE.)
            ! do nothing
         CASE(.FALSE.)
            IF(ioerr/=-1) WRITE(*,*) routinename//'Unexpected error in scattered trajs file. Err Code: ',ioerr
            EXIT
      END SELECT
      SELECT CASE(stat)
         CASE('Scattered')
            dp(1) = p(1)-this%inicond%trajs(id)%init_p(1)
            dp(2) = p(2)-this%inicond%trajs(id)%init_p(2)
            dk = MATMUL(to_rec_space,dp)
            g(1) = NINT(dk(1))
            g(2) = NINT(dk(2))
            DO j= 1,SIZE(this%peaks)
               SELECT CASE( this%peaks(j)%g(1)==g(1) .and. this%peaks(j)%g(2)==g(2) )
               CASE(.true.)
                  select case( dsin(r(5))==0.d0 )
                  case(.true.)
                     L2=p(5)**2.d0
                  case(.false.)
                     L2=p(5)**2.d0+(p(6)/dsin(r(5)))**2.d0
                  end select
                  Ecm=0.5d0*dot_product(p(1:3),p(1:3))/masa
                  Erot=0.5d0*L2/(mu*r(4)**2.d0)
                  Evibr=Etot-Ecm-Erot
                  finalJ=(-1.d0+dsqrt(1.d0+4.d0*L2))*0.5d0
                  finalV=dsqrt(1.d0-Evibr/morseEd)+dsqrt(2.d0*morseEd/masa)/morseWidth-0.5d0
                  select case( discretizeJ(finalJ,dJ)==0 )
                  case(.true.)
                     rovibrState=[ nint(finalV),0,0 ]

                  case(.false.)
                     select case( nint(dabs(p(6)))>abs(discretizeJ(finalJ,dJ)) )
                     case(.true.) ! aboid weird quantum states
                        rovibrState=[ nint(finalV),discretizeJ(finalJ,dJ),sign(discretizeJ(finalJ,dJ),nint(p(6))) ]

                     case(.false.) ! generic case
                        rovibrState=[ nint(finalV),discretizeJ(finalJ,dJ),nint(p(6)) ]

                     end select

                  end select
                  call this%addProbToSubPeak(peakId=j, rovibrState=rovibrState)
                  isAllowed=.true.
                  EXIT

               CASE(.false.)
                  isAllowed=.false.

               END SELECT
         END DO
      CASE DEFAULT
            ! do nothing
      END SELECT
      SELECT CASE(isAllowed)
         CASE(.true.)
            allowedScatt=allowedScatt+1
         CASE(.false.)
            WRITE(wuUnmap,formatUnmap) id,g(:)
      END SELECT
   END DO
   totscatt=i-1
   tottrajs=this%inicond%ntraj-this%inicond%nstart+1
   WRITE(*,*) "==========================================================="
   WRITE(*,*) "ASSIGN PEAKS TO TRAJS: total trajs: ",tottrajs
   WRITE(*,*) "ASSIGN PEAKS TO TRAJS: scattered trajs: ",totscatt
   WRITE(*,*) "ASSIGN PEAKS TO TRAJS: allowed scattered trajs:",allowedScatt
   WRITE(*,*) "ASSIGN PEAKS TO TRAJS: probability: ",dfloat(allowedScatt)/dfloat(tottrajs)
   WRITE(*,*) "==========================================================="
   CLOSE(unit=rwuMap)
   CLOSE(unit=ruScatt)
   CLOSE(unit=wuUnmap)
   RETURN
   contains
   ! included function
	function discretizeJ(J,dJ) result(finalJ)
	   implicit none
	   ! I/O variables
	   real(kind=8),intent(in):: J
	   integer(kind=4),intent(in):: dJ
	   ! function dummy variable
	   integer(kind=4):: finalJ
	   ! Local variables
	   real(kind=8):: deltaJ,diff
	   integer(kind=4):: i ! counter
	   integer(kind=4):: firstI
	   ! Run section
	   deltaJ=dfloat(dJ)*0.5d0
	   firstI=mod(this%inicond%init_qn(2),2)
      do i=firstI,1000,dJ ! almost infinite loop
	      diff=J-dfloat(i)
	      select case( dabs(diff)<=deltaJ )
	      case(.true.)
	         exit
	      case(.false.)
	         cycle
	      end select
	   enddo
	   finalJ=i
	   return
	end function discretizeJ
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
   CLASS(Allowed_peaksCRP6D),INTENT(IN):: this
   INTEGER(KIND=4),INTENT(IN):: in_n, in_m
   ! Local variables
   REAL(KIND=8):: left_term
   REAL(KIND=8):: n,m
   ! HEY, HO ! LET'S GO!
   ! Pointers & stuff
   n = DFLOAT(in_n)
   m = DFLOAT(in_m)
   ! left-term of the equation (A, B, C, D, E)
   left_term = this%conic(1)*(n**2.D0)+this%conic(2)*n*m+this%conic(3)*(m**2.D0)+this%conic(4)*n+this%conic(5)*m 
   ! Check value
   SELECT CASE( left_term<-this%conic(6) )
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
subroutine PRINT_LABMOMENTA_AND_ANGLES_ALLOWEDPEAKSCRP6D(this)
   implicit none
	! I/O Variables
   class(Allowed_peaksCRP6D),intent(in):: this
	! Local variables
	integer(kind=4):: dummy_int
	real(kind=8),dimension(9):: dummy_real
	integer(kind=4):: i ! counters
	integer(kind=4):: traj_id
	character(len=10):: stat
	real(kind=8),dimension(3):: p
	real(kind=8),dimension(3):: plab
   integer(kind=4):: ioerr
   real(kind=8):: psi,Theta,thetaout,beta
   real(kind=8),dimension(2,2):: mtrx
   ! Open units
   integer(kind=4),parameter:: wuFinal=12
   integer(kind=4),parameter:: ruScatt=11
   character(len=*),parameter:: formatFinal='(I6,1X,3(F10.5,1X),3(F10.5,1X))'
   ! RUN !! --------------------------
   beta=this%inicond%vpar_angle%getvalue()
   mtrx(1,:)=[dcos(beta),dsin(beta)]
   mtrx(2,:)=[-dsin(beta),dcos(beta)]
   open(unit=wuFinal,file="OUTANA6Dfinalpandangles.out",status="replace",action='write')
   write(wuFinal,*) "# ***** FINAL MOMENTA AND EXIT ANGLES *****"
   write(wuFinal,*) "# Format: id/Px,Py,Pz(a.u.)/Azimuthal,Polar,Deflection(rad)"
   write(wuFinal,*) "# -----------------------------------------------------------------------"
   open(unit=ruScatt,file="OUTDYN6Dscattered.out",status="old",action='read')
   call skipHeaderFromFile(unit=ruScatt)
   ioErr=0
   do while( ioErr==0 ) 
      read(ruScatt,*,iostat=ioerr) traj_id,stat,dummy_int,dummy_int,dummy_real(:),p(:)
      select case( stat=="Scattered" .and. ioErr==0 )
      case(.true.)
         plab(1:2)=matmul(mtrx,p(1:2))
         plab(3)=p(3)
         psi = datan(plab(2)/plab(1))
         Theta = datan(plab(2)/plab(3))
         thetaout=datan(plab(3)/dsqrt(plab(1)**2.D0+plab(2)**2.D0))
         write(wuFinal,formatFinal) traj_id,plab(:),psi,Theta,thetaout

      case(.false.)
         ! do nothing

      end select
   enddo
   close(unit=wuFinal)
   close(unit=ruScatt)
   return
end subroutine PRINT_LABMOMENTA_AND_ANGLES_ALLOWEDPEAKSCRP6D
!###########################################################
!# FUNCTION: getDiffOrderC4
!###########################################################
!> @brief
!! Given environment order and peak labels, sets real diffraction
!! order.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Nov/2014
!> @version 1.0
!-----------------------------------------------------------
function getDiffOrderC4(envOrder,g) result(diffOrder)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: envOrder
   INTEGER(KIND=4),DIMENSION(2),INTENT(IN) :: g
   ! Function dummy variable
   INTEGER(KIND=4):: diffOrder
   ! Run section
   SELECT CASE(g(1)==0 .AND. g(2)==0)
      CASE(.TRUE.)
         diffOrder=envOrder
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(g(1)==0 .OR. g(2)==0)
      CASE(.TRUE.)
         diffOrder=envOrder*(envOrder+1)/2
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(abs(g(1))==abs(g(2)))
      CASE(.TRUE.)
         diffOrder=((envOrder+1)*(envOrder+2)/2)-1
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(abs(g(1))==envOrder)
      CASE(.TRUE.)
         diffOrder=(envOrder*(envOrder+1)/2)+abs(g(2))
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(abs(g(2))==envOrder)
      CASE(.TRUE.)
         diffOrder=(envOrder*(envOrder+1)/2)+abs(g(1))
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
end function getDiffOrderC4
!######################################################
! FUNCTION: getPeakId_ALLOWEDPEAKSCRP6D
!######################################################
!> @brief
!! - Gets allowed peak ID number given its diffraction numbers.
!! - If there is not an allowed peak with diffraction number g(:),
!!   this function exits with a negative integer
!------------------------------------------------------
function getPeakId_ALLOWEDPEAKSCRP6D(this,g) result(peakId)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Allowed_peaksCRP6D),intent(in):: this
   integer(kind=4),dimension(2):: g
   ! Function dummy variable
   integer(kind=4):: peakId
   ! Local variables
   integer(kind=4):: i ! counter
   ! Run section
   select case( .not.allocated(this%peaks) )
   case(.true.) ! stop! badness!
      write(0,*) 'getPeakId_ALLOWEDPEAKSCRP6D ERR: Peaks are not allocated'
      call exit(1)
   case(.false.) ! Initialize values
      peakId = -1
      i=1
   end select
   do while( peakId < 0 .and. i <= size(this%peaks) )
      if( this%peaks(i)%g(1) == g(1) .and. this%peaks(i)%g(2) == g(2) ) peakId=this%peaks(i)%id
      i=i+1
   enddo
   return
end function getPeakId_ALLOWEDPEAKSCRP6D
!############################################################
! SUBROUTINE: createNewPeak_PEAKCRP6D
!############################################################
!------------------------------------------------------------
subroutine createNewPeak_ALLOWEDPEAKSCRP6D(this,envOrder,g)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Allowed_peaksCRP6D),intent(inout):: this
   integer(kind=4),dimension(2),intent(in):: g
   integer(kind=4),intent(in):: envOrder
   ! Local variables
   type(PeakCRP6D),dimension(:),allocatable:: auxListPeaks
   integer(kind=4):: oldN
   integer(kind=4):: idNew
   real(kind=8):: pinit_par,beta,theta_in,gama,mass,a,b
   ! Run section ----------------------------------------------
   select case( .not.allocated(this%peaks) )
   case(.true.)  ! this is the first allocation
      allocate( this%peaks(1) )
      idNew=1

   case(.false.) ! add new peak to the list
      oldN=size( this%peaks )
      idNew=oldN+1
      call move_alloc( from=this%peaks,to=auxListPeaks )
      allocate( this%peaks(idNew) )
      this%peaks(1:oldN)=auxListPeaks(1:oldN)

   end select
   ! Some parameters
	a=system_surface%norm_s1
	b=system_surface%norm_s2
   mass=sum(system_mass(:))
	beta=this%inicond%vpar_angle%getvalue()
	theta_in=this%inicond%vz_angle%getvalue()
	gama = dacos(dot_product(system_surface%s1,system_surface%s2)/(a*b))
	pinit_par = this%inicond%E_norm%getvalue()/(dtan(theta_in)**2.d0)
   this%peaks(idNew)%id=idNew
   this%peaks(idNew)%g(:)=g(:)
   this%peaks(idNew)%envOrder=envOrder
   this%peaks(idNew)%diffOrder=getDiffOrderC4(envOrder,g)
   ! Change in momentum in laboratory coordinates (by definition both components are orthogonal)
   this%peaks(idNew)%dkx=(2.d0*pi/(a*b*DSIN(gama)))*(DFLOAT(g(1))*b*DSIN(gama-beta)+DFLOAT(g(2))*a*DSIN(beta))
   this%peaks(idNew)%dky=(2.d0*pi/(a*b*DSIN(gama)))*(DFLOAT(g(2))*a*DCOS(beta)-DFLOAT(g(1))*b*DCOS(gama-beta))
   this%peaks(idNew)%dkxy=norm2([this%peaks(idNew)%dkx,this%peaks(idNew)%dky])
   ! Get angles
   select case( g(1)==0 .and. g(2)==0 )
   case(.true.)
      this%peaks(idNew)%psiOut=0.d0
   case(.false.)
      this%peaks(idNew)%psiOut=datan(this%peaks(idNew)%dky/(pinit_par+this%peaks(idNew)%dkx))
   end select
   return
end subroutine createNewPeak_ALLOWEDPEAKSCRP6D
!############################################################
! SUBROUTINE: addProbToPeak_ALLOWEDPEAKSCRP6D
!############################################################
!> @brief
!! Adds probability to a given subpeak. If it does not exist, this
!! routine will initialize it. Subpeaks are classified by their quantum
!! state.
!------------------------------------------------------------
subroutine addProbToSubPeak_ALLOWEDPEAKSCRP6D(this,peakId,rovibrState)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Allowed_peaksCRP6D),intent(inout):: this
   integer(kind=4):: peakId
   integer(kind=4),dimension(3):: rovibrState
   ! Local variables
   real(kind=8):: initE,kz,masa,theta_in,pinit_par
   integer(kind=4):: N,newCol,oldCol
   type(subPeakCRP6D),dimension(:),allocatable:: auxListSubPeaks
   logical:: isNew
   integer(kind=4):: i ! counter
   ! Run section
   ! Some
   masa=sum(system_mass(:))
   initE=this%inicond%evirot%getValue()
	theta_in=this%inicond%vz_angle%getvalue()
   N=size( this%inicond%trajs(:) )
	pinit_par = this%inicond%E_norm%getvalue()/(dtan(theta_in)**2.d0)
   isNew=.true.
   i=1
   select case( .not.allocated(this%peaks(peakId)%subPeaks) )
   case(.true.)
      oldCol=0
      newCol=oldCol+1

      allocate( this%peaks(peakId)%subPeaks(newCol) )

   case(.false.)
      oldCol=size( this%peaks(peakId)%subPeaks )
      newCol=oldCol+1
      do while( i<= oldCol .and. isNew )
         if( all( this%peaks(peakId)%subPeaks(i)%rovibrState == rovibrState ) ) isNew=.false.
         i=i+1
      enddo
      select case( isNew )
      case(.true.)
         call move_alloc( from=this%peaks(peakId)%subPeaks, to=auxListSubPeaks )
         allocate( this%peaks(peakId)%subPeaks(neWcol) )
         this%peaks(peakId)%subPeaks(1:oldCol) = auxListSubPeaks(1:oldCol)

      case(.false.)
         ! do nothing

      end select

   end select

   select case( isNew )
   case(.true.)
      this%peaks(peakId)%subPeaks(newCol)%rovibrState(:)=rovibrState(:)
      this%peaks(peakId)%subPeaks(newCol)%dE=initE-evaluateEnergyRovibrState(rovibrState)
      kz=dsqrt( 2.d0*masa*this%inicond%E_norm%getValue()-2.d0*pinit_par*this%peaks(peakId)%dkx &
                -this%peaks(peakId)%dkxy**2.d0-2.d0*masa*this%peaks(peakId)%subPeaks(newCol)%dE )
      this%peaks(peakId)%subPeaks(newCol)%prob=1.d0/N
      this%peaks(peakId)%subPeaks(newCol)%phiOut=datan(this%peaks(peakId)%dky/kz)
      this%peaks(peakId)%subPeaks(newCol)%thetaOut=datan( kz/norm2([this%peaks(peakId)%dkx+pinit_par,this%peaks(peakId)%dky]) )

   case(.false.)
      this%peaks(peakId)%subPeaks(i-1)%prob=this%peaks(peakId)%subPeaks(i-1)%prob+1.d0/N

   end select

   return
end subroutine addProbToSubPeak_ALLOWEDPEAKSCRP6D
!#######################################################
! SUBROUTINE: printAllowedPeaks_ALLOWEDPEAKSCRP6D
!#######################################################
!> @brief
!! Prints all allowed peaks.
!-------------------------------------------------------
subroutine printAllowedPeaks_ALLOWEDPEAKSCRP6D(this)
   ! initial declarations
   implicit none
   ! I/O variables
   class(Allowed_PeaksCRP6D),intent(in):: this
   ! Local variables
   integer(kind=4),parameter:: wuAllowed=20
   integer(kind=4):: i ! counter
   character(len=*),parameter:: formatAllowed='(3(I10,1X),4(F10.5,1X))'
   ! Run section ---------------------------------------
   open(unit=wuAllowed,file=this%fileNameAllowed,status='replace',action='write')
   write(wuAllowed,*) '# ------ ALLOWED PEAKS ---------------------------------------------'
   write(wuAllowed,*) '# Format: diffOrder/n,m/dkx,dky,dkxy(au)/psiOut'
   write(wuAllowed,*) '# ------------------------------------------------------------------'
   do i=1,size( this%peaks(:) )
      write(wuAllowed,formatAllowed) this%peaks(i)%diffOrder,this%peaks(i)%g(:),this%peaks(i)%dkx,this%peaks(i)%dky,&
                                     this%peaks(i)%dkxy,this%peaks(i)%psiOut
   enddo
   close(unit=wuAllowed)
   return
end subroutine printAllowedPeaks_ALLOWEDPEAKSCRP6D
!########################################################
! SUBROUTINE: sortByDiffOrder_ALLOWEDCRP6D
!########################################################
! @brief
!! sort peaks by diffraction order and not by peak id
!--------------------------------------------------------
subroutine sortByDiffOrder_ALLOWEDPEAKSCRP6D(this)
   ! initial declarations
   implicit none
   ! I/O variables
   class(Allowed_PeaksCRP6D),intent(inout):: this
   ! Local variables
   integer(kind=4),dimension(:),allocatable:: intList
   integer(kind=4):: N
   integer(kind=4):: diffOrder
   integer(kind=4):: i,j,init_i ! counters
   logical:: stopLoop
   type(peakCRP6D),dimension(:),allocatable:: auxPeaksList
   ! Run section -----------------------------------------
   N=size( this%peaks )
   allocate( intList(N) )
   i=1
   diffOrder=0
   stopLoop=.false.
   do while( .not.stopLoop .and. i<=N )
      init_i=i
      do j=1,N
         select case( this%peaks(j)%diffOrder==diffOrder )
         case(.true.)
            intList(i)=j
            i=i+1
         case(.false.)
            ! do nothing
         end select
      enddo
      select case( init_i==i )
      case(.true.)
         stopLoop=.true.
         if( init_i/=N ) write(0,*) 'sortByDiffOrder WARNING: cycling stopped before final intList'
      case(.false.)
         diffOrder=diffOrder+1
      end select
   enddo
   allocate( auxPeaksList(N) )
   auxPeaksList(:)=this%peaks(:)
   do i=1,N
      this%peaks(i)=auxPeaksList(intList(i))
   enddo
   return
end subroutine sortByDiffOrder_ALLOWEDPEAKSCRP6D
!###########################################################
! SUBROUTINE: printSeenPeaks_ALLOWEDPEAKSCRP6D
!###########################################################
!> @brief
!! Prints Seen peaks file
!-----------------------------------------------------------
subroutine printSeenPeaks_ALLOWEDPEAKSCRP6D(this)
   ! initial declarations
   implicit none
   ! I/O variables
   class(Allowed_PeaksCRP6D),intent(in):: this
   ! Local variables
   integer(kind=4),parameter:: wuSeen=123
   character(len=*),parameter:: formatSeen='(6(I10,1X),5(F10.5,1X))'
   integer(kind=4):: i,j ! counters
   ! Run section --------------------------------------------
   open(unit=wuSeen,file=this%fileNameSeen,status='replace',action='write')
   write(wuSeen,*) '# --------------------------------------- SEEN PEAKS --------------------------------------------------------'
   write(wuSeen,*) '# Format: diffOrder/n,m,V,J,mJ/Azimuthal angle/Deflection from perp. surf./Deflection from surf./dE(au)/prob'
   write(wuSeen,*) '# -----------------------------------------------------------------------------------------------------------'
   do i=1,size( this%peaks )
      select case( allocated(this%peaks(i)%subPeaks) )
      case(.true.)
         do j =1,size( this%peaks(i)%subPeaks )
            write(wuSeen,formatSeen) this%peaks(i)%diffOrder,this%peaks(i)%g(:),this%peaks(i)%subPeaks(j)%rovibrState(:),&
                                     this%peaks(i)%psiOut,this%peaks(i)%subPeaks(j)%PhiOut,this%peaks(i)%subPeaks(j)%thetaOut,&
                                     this%peaks(i)%subPeaks(j)%dE,this%peaks(i)%subPeaks(j)%prob
         enddo
      case(.false.)
         ! do nothing. do not print anything
      end select
   enddo
   close(unit=wuSeen)
   return
end subroutine printSeenPeaks_ALLOWEDPEAKSCRP6D
END MODULE DIFFRACTIONCRP6D_MOD
