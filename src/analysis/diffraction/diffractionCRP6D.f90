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
   character(len=25):: fileNameUnmapped='OUTANA6DunmappedTrajs.out'
   CONTAINS
      ! private tools section
      procedure,private:: getPeakId => getPeakId_ALLOWEDPEAKSCRP6D
      procedure,private:: addNewPeak => addNewPeak_ALLOWEDPEAKSCRP6D
      procedure,private:: addProbToSubPeak => addProbToSubPeak_ALLOWEDPEAKSCRP6D
      ! public tools section
      PROCEDURE,PUBLIC:: INITIALIZE => INITIALIZE_ALLOWEDPEAKSCRP6D
      PROCEDURE,PUBLIC:: getEnvOrder => getEnvOrder_ALLOWEDPEAKSCRP6D
      procedure,public:: printAllowedPeaks => printAllowedPeaks_ALLOWEDPEAKSCRP6D
      procedure,public:: printSeenPeaks => printSeenPeaks_ALLOWEDPEAKSCRP6D
      procedure,public:: sortByDiffOrder => sortByDiffOrder_ALLOWEDPEAKSCRP6D
      PROCEDURE,PUBLIC:: assignTrajsToPeaks => AssignTrajsToPeaks_ALLOWEDPEAKSCRP6D
      PROCEDURE,PUBLIC:: PRINT_LABMOMENTA_AND_ANGLES => PRINT_LABMOMENTA_AND_ANGLES_ALLOWEDPEAKSCRP6D
      PROCEDURE,PUBLIC:: isAllowed => isAllowed_ALLOWEDPEAKSCRP6D
      procedure,public:: quantizeRovibrState => quantizeRovibrState_ALLOWEDPEAKSCRP6D
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
!######################################################################
!# FUNCTION: getEnvOrder_ALLOWEDPEAKSCRP6D ############################
!######################################################################
!> @brief
!! - Search for the environmental order of a given diffraction state.
!! - Environmental orders are general for all cells
!----------------------------------------------------------------------
function getEnvOrder_ALLOWEDPEAKSCRP6D(this,diffState) result(envOrder)
   implicit none
   ! I/O variables
   class(Allowed_peaksCRP6D),intent(inout):: this
   integer(kind=4),dimension(2),intent(in):: diffState
   ! Function dummy variable
   integer(kind=4):: envOrder
   ! Local variables
	INTEGER(KIND=4):: i,k ! Counters
	INTEGER(KIND=4):: order
	INTEGER(KIND=4),DIMENSION(2):: g ! (n,m) vector
   ! Parameters
   character(len=*),parameter:: routinename="SET_ALLOWED_PEAKS: "
	! FIRE IN THE HOLE ! -----------------------------------------------
	order = 0
	g(:) = 0
	select case( all(diffState(:)==g(:)) ) ! order zero
	case(.true.) ! I found you!
	   envOrder=order
	   return
	case(.false.) ! continue searching
	end select
	do
		order = order + 1
		i = 0
		!-----
		g(1) = order
		do k = -order, order
			g(2) = k
			i = i + 1
      	select case( all(diffState(:)==g(:)) )
	      case(.true.) ! I found you!
	         envOrder=order
	         return
	      case(.false.) ! continue searching
	      end select
		enddo
		!----------------------
		g(1) = -order
		do k = -order, order
			g(2) = k
			i = i + 1
      	select case( all(diffState(:)==g(:)) )
	      case(.true.) ! I found you!
	         envOrder=order
	         return
	      case(.false.) ! continue searching
	      end select
		enddo
		!----
		g(2) = order
		do k = -order +1, order-1
			g(1) = k
			i = i + 1
      	select case( all(diffState(:)==g(:)) )
	      case(.true.) ! I found you!
	         envOrder=order
	         return
	      case(.false.) ! continue searching
	      end select
		enddo
		!----
		g(2) = -order
		do k = -order +1, order-1
			g(1) = k
			i = i + 1
      	select case( all(diffState(:)==g(:)) )
	      case(.true.) ! I found you!
	         envOrder=order
	         return
	      case(.false.) ! continue searching
	      end select
		enddo
	enddo
end function getEnvOrder_ALLOWEDPEAKSCRP6D
!####################################################################################
! SUBROUTINE: assignTrajsToPeaks_ALLOWEDPEAKSCRP6D
!####################################################################################
! - At the moment only works with rectancular cells
!------------------------------------------------------------------------------------
subroutine assignTrajsToPeaks_ALLOWEDPEAKSCRP6D(this)
   implicit none
   ! I/O variables
   class(Allowed_peaksCRP6D),intent(inout):: this
   ! Local variables
   integer(kind=4):: totScatt,totTrajs,allowedScatt
   integer(kind=4),dimension(2):: g
	real(kind=8),dimension(2,2):: to_rec_space
	real(kind=8),dimension(6):: p,r ! final momentum and position
	real(kind=8),dimension(2):: dp ! variation of momentum
	real(kind=8),dimension(2):: dk ! variation of momentum in rec. space coord
	real(kind=8):: gama,a,b,Etot,dE
   character(len=10):: stat
   integer(kind=4):: ioErr
   integer(kind=4),dimension(3):: rovibrState
   integer(kind=4):: peakId
   integer(kind=4):: i ! counters
   integer(kind=4):: id
   ! Auxiliar variables
	real(kind=8),dimension(2):: auxReal
	integer(kind=4),dimension(2):: auxInt
   ! Some parameters
   character(len=*),parameter:: routinename = "assigtTrajsToPeaks_ALLOWEDPEAKSCVRP6D: "
   character(len=*),parameter:: formatUnmap='(I6," ---> ",5(I5))'        ! Id/--->/n/m/v,J,mJ
   ! Read/write units
   integer(kind=4),parameter:: wuUnmap=13
   integer(kind=4),parameter:: ruScatt=14
   ! RUN SECTION -------------------------
	a=system_surface%norm_s1
	b=system_surface%norm_s2
   gama=dacos(dot_product(system_surface%s1,system_surface%s2)/(a*b))
   to_rec_space(1,1) = a/(2.D0*PI)
   to_rec_space(1,2) = 0.D0
   to_rec_space(2,1) = b*DCOS(gama)/(2.D0*PI)
   to_rec_space(2,2) = b*DSIN(gama)/(2.D0*PI)
   ! binning parameters (only Morse implemented)
   !---------
   open(unit=wuUnmap,file=this%fileNameUnmapped,status='replace',action='write')
   write(wuUnmap,'("# *************** LIST OF UNMAPPED TRAJS ***************")')
   write(wuUnmap,'("# Format: id/n,m/v,J,mJ")')
   write(wuUnmap,'("# ----------------------------------------------------------------")')
   open(unit=ruScatt,file="OUTDYN6Dscattered.out",status="old",action='read')
   call skipHeaderFromFile(unit=ruScatt)
   i=0
   totScatt=0
   allowedScatt=0
   ioErr=0
   do while( ioErr == 0 )
      ! read from scattered file
      read(ruScatt,*,ioStat=ioErr) id,stat,auxInt(:),Etot,auxReal(:),r(:),p(:)
      select case( ioErr==0 .and. stat=='Scattered' )
      case(.true.) ! secure to operate
         dp(1) = p(1)-this%inicond%trajs(id)%init_p(1)
         dp(2) = p(2)-this%inicond%trajs(id)%init_p(2)
         dk = matmul(to_rec_space,dp)
         g(1) = nint(dk(1))
         g(2) = nint(dk(2))
         rovibrState=this%quantizeRovibrState( Etot=Etot,position=r(:),momenta=p(:) )
         totScatt=totScatt+1
         select case( this%isAllowed(diffState=g(:),rovibrState=rovibrState,diffEnergy=dE) )
         case(.true.) ! Allowed peak: add to list
            call this%addNewPeak(g, peakId=peakId)
            call this%addProbToSubPeak( peakId=peakId,rovibrState=rovibrState )
            allowedScatt=allowedScatt+1
         case(.false.)
            write(wuUnmap,formatUnmap) id,g(:),rovibrState(:)
         end select
      case(.false.) ! not secure to operate, next switch
         select case( ioErr )
         case(-1)
            ! do nothing, EOF reached, let it break the cycle
         case default
            write(0,*) routinename//'ERR: unexpected error encountered. Error Code: ',ioErr
            call exit(1)
         end select
      end select
   enddo
   totTrajs=this%inicond%ntraj-this%inicond%nstart+1
   write(*,'("===========================================================")')
   write(*,'("ASSIGN PEAKS TO TRAJS: total trajs: ",I10)')             totTrajs
   write(*,'("ASSIGN PEAKS TO TRAJS: scattered trajs: ",I10)')         totScatt
   write(*,'("ASSIGN PEAKS TO TRAJS: allowed scattered trajs: ",I10)') allowedScatt
   write(*,'("ASSIGN PEAKS TO TRAJS: probability: ",F10.5)')           dfloat(allowedScatt)/dfloat(totTrajs)
   write(*,'("===========================================================")')
   close(unit=ruScatt)
   close(unit=wuUnmap)
   call deleteIfEmptyFile( this%fileNameUnmapped )
   return
end subroutine assignTrajsToPeaks_ALLOWEDPEAKSCRP6D
!####################################################################################
! FUNCTION: EVALUATE_PEAK ###########################################################
!####################################################################################
! - TRUE if A*(n^2) + B*nm + C*(m^2) + D*n + E*m + F < 0
!------------------------------------------------------------------------------------
function isAllowed_ALLOWEDPEAKSCRP6D(this,diffState,rovibrState,diffEnergy) result(isAllowed)
   IMPLICIT NONE
   ! I/O variables
   class(Allowed_peaksCRP6D),intent(in):: this
   integer(kind=4),dimension(2),intent(in):: diffState
   integer(kind=4),dimension(3),intent(in):: rovibrState
   real(kind=8),optional,intent(out):: diffEnergy
   ! Function dummy variable
   logical:: isAllowed
   ! Local variables
   real(kind=8):: n,m,a,b,dE,beta,theta_in,masa,gama,pinit_par
   real(kind=8),dimension(2):: kinit_par
   real(kind=8),dimension(6):: C ! conic coefficients
   ! HEY, HO ! LET'S GO! ---------------------------------------
   ! Some parameters
	a=system_surface%norm_s1
	b=system_surface%norm_s2
	masa=sum( system_mass(:) )
   gama=dacos(dot_product(system_surface%s1,system_surface%s2)/(a*b))
	beta=this%inicond%vpar_angle%getvalue()
	theta_in=this%inicond%vz_angle%getvalue()
	pinit_par = (1.d0/dtan(theta_in))*dsqrt(2.d0*masa*this%inicond%E_norm%getvalue())
	kinit_par(1)=pinit_par*a*dcos(beta)/(2.d0*pi)
	kinit_par(2)=pinit_par*b*dcos(gama-beta)/(2.d0*pi)
   n=dfloat(diffState(1))
   m=dfloat(diffState(2))
   dE=evaluateEnergyRovibrState(rovibrState)-this%inicond%evirot%getValue()
   ! Set conic coefficients
	C(1) =  b**2.D0
	C(2) = -2.D0*a*b*dcos(gama)
	C(3) =  a**2.D0
	C(4) =  2.D0*b*(b*kinit_par(1)-a*kinit_par(2)*dcos(gama))
	C(5) =  2.D0*a*(a*kinit_par(2)-b*kinit_par(1)*dcos(gama))
   C(6) = -((a*b*dsin(gama)/(2.D0*pi))**2.D0)*2.D0*masa*(this%inicond%E_norm%getvalue()-dE)
	!C(6) = -((a*b*dsin(gama)/(2.D0*pi))**2.D0)*2.D0*masa*(this%inicond%E_norm%getvalue())
	! Allowed condition
   isAllowed=( C(1)*(n**2.D0)+C(2)*n*m+C(3)*(m**2.D0)+C(4)*n+C(5)*m+C(6) < 0.d0 )
	if( present(diffEnergy) ) diffEnergy=dE
   return
end function isAllowed_ALLOWEDPEAKSCRP6D
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
! SUBROUTINE: addNewPeak_PEAKCRP6D
!############################################################
!------------------------------------------------------------
subroutine addNewPeak_ALLOWEDPEAKSCRP6D(this,g,peakId)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Allowed_peaksCRP6D),intent(inout):: this
   integer(kind=4),dimension(2),intent(in):: g
   integer(kind=4),optional,intent(out):: peakId
   ! Local variables
   type(PeakCRP6D),dimension(:),allocatable:: auxListPeaks
   integer(kind=4):: oldN
   integer(kind=4):: idNew
   integer(kind=4):: i ! counters
   real(kind=8):: pinit_par,beta,theta_in,gama,masa,a,b
   logical:: alreadyExists
   ! Run section ----------------------------------------------
   alreadyExists=.false.
   select case( .not.allocated(this%peaks) )
   case(.true.)  ! this is the first allocation
      allocate( this%peaks(1) )
      idNew=1

   case(.false.) ! add new peak to the list
      oldN=size( this%peaks )
      do i=1,oldN
         select case( all(this%peaks(i)%g(:)==g(:)) )
         case(.true.)
            if( present(peakId) ) peakId=i
            return
         case(.false.)
            ! do nothing
         end select
      enddo
      idNew=oldN+1
      call move_alloc( from=this%peaks,to=auxListPeaks )
      allocate( this%peaks(idNew) )
      this%peaks(1:oldN)=auxListPeaks(1:oldN)
   end select
   ! Some parameters
	a=system_surface%norm_s1
	b=system_surface%norm_s2
   masa=sum( system_mass(:) )
	beta=this%inicond%vpar_angle%getvalue()
	theta_in=this%inicond%vz_angle%getvalue()
	gama = dacos(dot_product(system_surface%s1,system_surface%s2)/(a*b))
	pinit_par = (1.d0/dtan(theta_in))*dsqrt(2.d0*masa*this%inicond%E_norm%getvalue())
   this%peaks(idNew)%id=idNew
   this%peaks(idNew)%g(:)=g(:)
   this%peaks(idNew)%envOrder=this%getEnvOrder( diffState=g(:) )
   this%peaks(idNew)%diffOrder=getDiffOrderC4( this%peaks(idNew)%envOrder,g(:) )
   ! Change in momentum in laboratory coordinates (by definition both components are orthogonal)
   this%peaks(idNew)%dkx=(2.d0*pi/(a*b*DSIN(gama)))*(DFLOAT(g(1))*b*DSIN(gama-beta)+DFLOAT(g(2))*a*DSIN(beta))
   this%peaks(idNew)%dky=(2.d0*pi/(a*b*DSIN(gama)))*(DFLOAT(g(2))*a*DCOS(beta)-DFLOAT(g(1))*b*DCOS(gama-beta))
   this%peaks(idNew)%dkxy=norm2([this%peaks(idNew)%dkx,this%peaks(idNew)%dky])
   ! Get angles
   select case( all( g(:)==[0,0] ) )
   case(.true.)
      this%peaks(idNew)%psiOut=0.d0
   case(.false.)
      this%peaks(idNew)%psiOut=datan(this%peaks(idNew)%dky/(pinit_par+this%peaks(idNew)%dkx))
   end select
   if( present(peakId) ) peakId=idNew
   return
end subroutine addNewPeak_ALLOWEDPEAKSCRP6D
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
   integer(kind=4),intent(in):: peakId
   integer(kind=4),dimension(3),intent(in):: rovibrState
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
	pinit_par = (1.d0/dtan(theta_in))*dsqrt(2.d0*masa*this%inicond%E_norm%getvalue())
   i=1
   isNew=.true.
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
      this%peaks(peakId)%subPeaks(newCol)%dE=evaluateEnergyRovibrState(rovibrState)-initE
      kz=dsqrt( 2.d0*masa*this%inicond%E_norm%getValue()-2.d0*pinit_par*this%peaks(peakId)%dkx &
                -this%peaks(peakId)%dkxy**2.d0-2.d0*masa*this%peaks(peakId)%subPeaks(newCol)%dE )
      this%peaks(peakId)%subPeaks(newCol)%prob=1.d0/N
      this%peaks(peakId)%subPeaks(newCol)%phiOut=datan( this%peaks(peakId)%dky/kz )
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
   write(wuAllowed,'("# ------ ALLOWED PEAKS ---------------------------------------------")')
   write(wuAllowed,'("# Format: diffOrder/n,m/dkx,dky,dkxy(au)/psiOut")')
   write(wuAllowed,'("# ------------------------------------------------------------------")')
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
   integer(kind=4):: i,j,diffOrder ! counters
   type(peakCRP6D),dimension(:),allocatable:: auxPeaksList
   ! Run section -----------------------------------------
   N=size( this%peaks )
   allocate( intList(N) )
   i=1
   diffOrder=0
   do while( i<=N )
      do j=1,N
         select case( this%peaks(j)%diffOrder==diffOrder )
         case(.true.)
            intList(i)=j
            i=i+1
         case(.false.)
            ! do nothing
         end select
      enddo
      diffOrder=diffOrder+1
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
   character(len=*),parameter:: formatSeen='(6(I10,1X),6(F10.5,1X))'
   integer(kind=4):: i,j ! counters
   ! Run section --------------------------------------------
   open(unit=wuSeen,file=this%fileNameSeen,status='replace',action='write')
   write(wuSeen,'("# --------------------------------------- SEEN PEAKS ----------------------------------------------------------------")')
   write(wuSeen,'("# Format: diffOrder/n,m,V,J,mJ/Azimuthal angle/Deflection from perp. surf./Deflection from surf./dKxy(au)/dE(au)/prob")')
   write(wuSeen,'("# -------------------------------------------------------------------------------------------------------------------")')
   do i=1,size( this%peaks )
      select case( allocated(this%peaks(i)%subPeaks) )
      case(.true.)
         do j =1,size( this%peaks(i)%subPeaks )
            write(wuSeen,formatSeen) this%peaks(i)%diffOrder,this%peaks(i)%g(:),this%peaks(i)%subPeaks(j)%rovibrState(:),&
                                     this%peaks(i)%psiOut,this%peaks(i)%subPeaks(j)%PhiOut,this%peaks(i)%subPeaks(j)%thetaOut,&
                                     this%peaks(i)%dkxy,this%peaks(i)%subPeaks(j)%dE,this%peaks(i)%subPeaks(j)%prob
         enddo
      case(.false.)
         ! do nothing. do not print anything
      end select
   enddo
   close(unit=wuSeen)
   return
end subroutine printSeenPeaks_ALLOWEDPEAKSCRP6D
!###########################################################
! FUNCTION: quantizeRovibrState_ALLOWEDPEAKSCRP6D
!###########################################################
!-----------------------------------------------------------
function quantizeRovibrState_ALLOWEDPEAKSCRP6D(this,Etot,position,momenta) result(rovibrState)
   ! initial declarations
   implicit none
   ! I/O variables
   class(Allowed_peaksCRP6D),intent(in):: this
   real(kind=8),intent(in):: Etot
   real(kind=8),dimension(6),intent(in):: position
   real(kind=8),dimension(6),intent(in):: momenta
   ! Dummy function variable
   integer(kind=4),dimension(3):: rovibrState
   ! Local variables
   real(kind=8):: L2,mu,masa,Ecm,Evibr,Erot,finalJ,finalV
   ! Parameters
   character(len=*),parameter:: routinename='quantizeRovibrState_ALLOWEDPEAKSCRP6D: '
   ! Run section ------------------------------------------
   masa=sum( system_mass(:) )
   mu=product( system_mass(:) )/masa
   select case( dsin(position(5))==0.d0 )
   case(.true.)
      L2=momenta(5)**2.d0
   case(.false.)
      L2=momenta(5)**2.d0+(momenta(6)/dsin(position(5)))**2.d0
   end select
   Ecm=0.5d0*dot_product( momenta(1:3),momenta(1:3))/masa
   Erot=0.5d0*L2/(mu*position(4)**2.d0)
   Evibr=Etot-Ecm-Erot
   select case( system_binningScheme)
   case('Morse')
      finalJ=(-1.d0+dsqrt(1.d0+4.d0*L2))*0.5d0
      finalV=dsqrt(1.d0-Evibr/system_binningParam(1))+dsqrt(2.d0*system_binningParam(1)/masa)/system_binningParam(3)-0.5d0
   case default
      write(0,*) routinename//'ERR: bad binning Scheme'
      write(0,*) 'Implemented ones: Morse'
      write(0,*) 'Case sensitive'
      call exit(1)
   end select

   if( discretizeJ(finalJ)==0 ) then
      rovibrState=[ nint(finalV),0,0 ]
   else if( nint(dabs(momenta(6)))>abs(discretizeJ(finalJ)) ) then
      rovibrState=[ nint(finalV),discretizeJ(finalJ),sign(discretizeJ(finalJ),nint(momenta(6))) ]
   else
      rovibrState=[ nint(finalV),discretizeJ(finalJ),nint(momenta(6)) ]
   endif
   return ! ACTUAL END OF THIS FUNCTION

   contains
      ! ////////////////////////////////////////////
      !             INCLUDED FUNCTIONS
      !/////////////////////////////////////////////

      !#############################################
      ! FUNCTION: discretizeJ
      !#############################################
      !> @brief
      !! simple function to discretize J quantum number taking
      !! into account the selection rule dJ
      !---------------------------------------------
	   function discretizeJ(J) result(finalJ)
	      implicit none
	      ! I/O variables
	      real(kind=8),intent(in):: J
	      ! function dummy variable
	      integer(kind=4):: finalJ
	      ! Local variables
	      real(kind=8):: deltaJ,diff
	      integer(kind=4):: i ! counter
	      integer(kind=4):: firstI
	      ! Run section
	      deltaJ=dfloat(system_binningdJ)*0.5d0
	      firstI=mod(this%inicond%init_qn(2),2)
         do i=firstI,1000,system_binningdJ ! almost infinite loop
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

end function quantizeRovibrState_ALLOWEDPEAKSCRP6D
END MODULE DIFFRACTIONCRP6D_MOD
