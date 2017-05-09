MODULE DIFFRACTIONGROW_MOD
use SYSTEM_MOD
use INITSURFACEGROW_MOD
#ifdef DEBUG
use DEBUG_MOD, only: VERBOSE_WRITE, DEBUG_WRITE
#endif
IMPLICIT NONE
!////////////////////////////////////////////////////////
! TYPE: SubPeakGROW
!////////////////////////////////////////////////////////
!> @param robivrState(:) - integer(kind=4): V,J,mJ state
!> @param phiOut - real(kind=8): deflection angle respect to the perpendicular plane to the surface
!> @param thetaOut - real(kind=8): deflection angle respect to surface plane
!> @param prob - real(kind=8): Probability
!> @param dE - real(kind=8): internal energy exchange after data binning
!------------------------------------------------------
type:: SubPeakGROW
   private
   real(kind=8):: phiOut
   real(kind=8):: thetaOut
   integer(kind=4),dimension(3):: rovibrState
   real(kind=8):: prob
   real(kind=8):: deviation
   integer(kind=4):: ntrajs
   real(kind=8):: dE
end type SubPeakGROW
!======================================================
! TYPE: PeakGROW
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
TYPE :: PeakGROW
   PRIVATE
   INTEGER(KIND=4):: id
   integer(kind=4):: envOrder
   integer(kind=4):: diffOrder
   real(kind=8):: dkx
   real(kind=8):: dky
   real(kind=8):: dkxy
   INTEGER(KIND=4),DIMENSION(2):: g
	REAL(KIND=8):: psiOut
	type(SubPeakGROW),dimension(:),allocatable:: subPeaks
END TYPE PeakGROW
!======================================================
! Allowed_peaksGROW derived data
!----------------------------
TYPE :: Allowed_peaksGROW
   PRIVATE
   TYPE(InitSurfaceGrowDiatomic):: inicond
   class(Pes),allocatable:: thispes
   REAL(KIND=8):: E
   REAL(KIND=8),DIMENSION(6):: conic
   TYPE(PeakGROW),DIMENSION(:),ALLOCATABLE:: peaks
   integer(kind=4):: totAllowedScatt
   integer(kind=4):: totScatt
   character(len=24):: fileNameAllowed='OUTANA6DallowedPeaks.out'
   character(len=21):: fileNameSeen='OUTANA6DseenPeaks.out'
   character(len=25):: fileNameUnmapped='OUTANA6DunmappedTrajs.out'
   CONTAINS
      ! private tools section
      procedure,private:: getPeakId => getPeakId_ALLOWEDPEAKSGROW
      procedure,private:: addNewPeak => addNewPeak_ALLOWEDPEAKSGROW
      procedure,private:: addCountToSubPeak => addCountToSubPeak_ALLOWEDPEAKSGROW
      procedure,private:: doSubPeakStatistics => doSubPeakStatistics_ALLOWEDPEAKSGROW
      ! public tools section
      PROCEDURE,PUBLIC:: INITIALIZE => INITIALIZE_ALLOWEDPEAKSGROW
      PROCEDURE,PUBLIC:: getEnvOrder => getEnvOrder_ALLOWEDPEAKSGROW
      procedure,public:: printAllowedPeaks => printAllowedPeaks_ALLOWEDPEAKSGROW
      procedure,public:: printSeenPeaks => printSeenPeaks_ALLOWEDPEAKSGROW
      procedure,public:: sortByDiffOrder => sortByDiffOrder_ALLOWEDPEAKSGROW
      PROCEDURE,PUBLIC:: assignTrajsToPeaks => AssignTrajsToPeaks_ALLOWEDPEAKSGROW
      PROCEDURE,PUBLIC:: PRINT_LABMOMENTA_AND_ANGLES => PRINT_LABMOMENTA_AND_ANGLES_ALLOWEDPEAKSGROW
      PROCEDURE,PUBLIC:: isAllowed => isAllowed_ALLOWEDPEAKSGROW
      procedure,public:: quantizeRovibrState => quantizeRovibrState_ALLOWEDPEAKSGROW
      procedure,public:: quantizeDiffState => quantizeDiffState_ALLOWEDPEAKSGROW
END TYPE Allowed_peaksGROW
!=======================================================
CONTAINS
subroutine initialize_ALLOWEDPEAKSGROW(this)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(Allowed_peaksGROW),intent(out):: this
   ! Run section
   call this%inicond%initialize()
   call this%inicond%generate_trajs_from_file('start1.dat')
   return
end subroutine
!######################################################################
!# FUNCTION: getEnvOrder_ALLOWEDPEAKSGROW ############################
!######################################################################
!> @brief
!! - Search for the environmental order of a given diffraction state.
!! - Environmental orders are general for all cells
!----------------------------------------------------------------------
function getEnvOrder_ALLOWEDPEAKSGROW(this,diffState) result(envOrder)
   implicit none
   ! I/O variables
   class(Allowed_peaksGROW),intent(inout):: this
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
end function getEnvOrder_ALLOWEDPEAKSGROW
!####################################################################################
! SUBROUTINE: assignTrajsToPeaks_ALLOWEDPEAKSGROW
!####################################################################################
!------------------------------------------------------------------------------------
subroutine assignTrajsToPeaks_ALLOWEDPEAKSGROW(this)
   implicit none
   ! I/O variables
   class(Allowed_peaksGROW),intent(inout):: this
   ! Local variables
   integer(kind=4):: totScatt,totTrajs,allowedScatt
   integer(kind=4),dimension(2):: g
	real(kind=8),dimension(2,2):: to_rec_space
	real(kind=8),dimension(6):: p,r ! final momentum and position
	real(kind=8),dimension(12):: phaseSpaceVect
	real(kind=8),dimension(2):: dp ! variation of momentum
	real(kind=8),dimension(2):: dk ! variation of momentum in rec. space coord
	real(kind=8):: gama,a,b,Etot,dE
   integer(kind=4):: ioErr
   integer(kind=4),dimension(3):: rovibrState
   integer(kind=4):: peakId
   integer(kind=4):: i ! counters
   integer(kind=4):: id
   integer(kind=4):: stat
   real(kind=8):: mtot, mu
   ! Auxiliar variables
	real(kind=8),dimension(2):: auxReal
   integer(kind=4),dimension(2):: auxInt
   integer(kind=1):: controlSurfDyn
   character(len=4):: auxString
   character(len=1):: controlChar
   ! Some parameters
   character(len=*),parameter:: routinename = "assigtTrajsToPeaks_ALLOWEDPEAKSCVRP6D: "
   character(len=*),parameter:: formatUnmap='(I6," ---> ",5(I5))'        ! Id/--->/n/m/v,J,mJ
   ! Read/write units
   integer(kind=4),parameter:: wuUnmap=13
   integer(kind=4),parameter:: ruScatt=14
   integer(kind=4),parameter:: ruStat=9
   integer(kind=4),parameter:: ruSurfdyn=8
   ! RUN SECTION -------------------------
   a=system_surface%norm_s1
   b=system_surface%norm_s2
   mtot=sum( system_mass(:) )
   mu=product( system_mass(:) )/mtot
   gama=system_surface%angle
   ! binning parameters (only Morse implemented)
   !---------
   ioErr=0

   open(unit=ruStat,file='OUT_INTERPSUMMARY',status='old',action='read')
   do while( ioErr == 0 )
      read(ruStat,*,iostat=ioErr) id,auxReal(1),auxInt(:),stat,auxReal(:),auxReal(:),controlChar
      if (ioErr==0) then
         select case(controlChar)
         case('X')
            ! as expected, yeah
         case default
            write(0,*) routinename//'ERR: Reading OUT_INTERPSUMMARY: need X as final character'
            call exit(1)
         end select
         select case(stat)
         case(12)
            this%inicond%trajs(id)%stat='Scattered'
         case(-1,-2,-12)
            this%inicond%trajs(id)%stat='React'
         case default
            this%inicond%trajs(id)%stat='Unknown'
         end select
      endif
   enddo
   close(ruStat)
   
   ioErr=0
   open(unit=ruSurfdyn,file='IN_SURFDYN',status='old',action='read',iostat=ioErr)
   select case(ioErr)
   case(0)
      read(ruSurfdyn,*)
      read(ruSurfdyn,*) controlSurfDyn
   case default
      controlSurfDyn=0   
   end select
   close(unit=ruSurfdyn,iostat=ioErr)

   open(unit=wuUnmap,file=this%fileNameUnmapped,status='replace',action='write')
   write(wuUnmap,'("# *************** LIST OF UNMAPPED TRAJS ***************")')
   write(wuUnmap,'("# Format: id/n,m/v,J,mJ")')
   write(wuUnmap,'("# ----------------------------------------------------------------")')
   open(unit=ruScatt,file="OUT_FINALCV",status="old",action='read')
   i=0
   totScatt=0
   allowedScatt=0
   ioErr=0
   do while( ioErr == 0 )
      ! read from scattered file
      read(ruScatt,*,iostat=ioErr) auxString,auxString,auxString,auxString,id
      read(ruScatt,*,iostat=ioErr) r(1:3)
      read(ruScatt,*,iostat=ioErr) r(4:6)
      if(controlSurfDyn==1) read(ruScatt,*,iostat=ioErr) 
      read(ruScatt,*,iostat=ioErr)
      read(ruScatt,*,iostat=ioErr) p(1:3)
      read(ruScatt,*,iostat=ioErr) p(4:6)
      if(controlSurfDyn==1) read(ruScatt,*,iostat=ioErr) 
      ! correct momenta
      p(1:3)=system_mass(1)*p(1:3)/dsqrt(pmass2au)
      p(4:6)=system_mass(2)*p(4:6)/dsqrt(pmass2au)
      phaseSpaceVect(1:6)=r(:)
      phaseSpaceVect(7:12)=p(:)
      phaseSpaceVect(:)=from_atomic_to_molecular_phaseSpace(atomcoord=phaseSpaceVect)
      r(:)=phaseSpaceVect(1:6)
      p(:)=phaseSpaceVect(7:12)
      select case( ioErr==0 .and. this%inicond%trajs(id)%stat=='Scattered' )
      case(.true.) ! secure to operate
         dp(:) = p(1:2)-this%inicond%trajs(id)%init_p(1:2)
         g(:) = this%quantizeDiffState( dp )
         Etot=0.5d0*dot_product(p(1:3),p(1:3))/mtot+&                      ! Kinetic energy
              0.5d0*(p(4)**2.d0)/mu+&                                      ! Internal kinetic energy (1)
              0.5d0*(p(5)**2.d0+(p(6)/dsin(r(5)))**2.d0)/(mu*r(4)**2.d0)+& ! Internal kinetic energy (2)
              this%inicond%vibrPot%getPot( r(4) )                          ! Potential energy
         rovibrState=this%quantizeRovibrState( Etot=Etot,position=r(:),momenta=p(:) )
         totScatt=totScatt+1
         !write(*,*) '-----------------------------'
         !write(*,*) 'Id: ',id
         !write(*,*) 'Limit: ',this%inicond%trajs(id)%init_r(3)
         !write(*,*) 'initial momentum'
         !write(*,*) this%inicond%trajs(id)%init_p(:)
         !write(*,*) 'final momentum'
         !write(*,*) p(:)
         !write(*,*) 'Momentum variation in aux, recip and recip rounded coordinates'
         !write(*,*) dp(:)
         !write(*,*) g(:)
         !write(*,*) 'Ekin: ',0.5d0*dot_product(p(1:3),p(1:3))/mtot
         !write(*,*) 'Eint: ',Etot-0.5d0*dot_product(p(1:3),p(1:3))/mtot
         select case( this%isAllowed(diffState=g(:),rovibrState=rovibrState,diffEnergy=dE) )
         case(.true.) ! Allowed peak: add to list
            call this%addNewPeak(g, peakId=peakId)
            call this%addCountToSubPeak( peakId=peakId,rovibrState=rovibrState )
            allowedScatt=allowedScatt+1
         case(.false.)
            write(wuUnmap,formatUnmap) id,g(:),rovibrState(:)
         end select
      case(.false.) ! not secure to operate, next switch
         select case( ioErr )
         case(-1,5001)
            ! do nothing, EOF reached, let it break the cycle
         case(0)
            ! just ignore it
            write(wuUnmap,formatUnmap) id,g(:),rovibrState(:)
         case default
            write(0,*) routinename//'ERR: unexpected error encountered. Error Code: ',ioErr
            call exit(1)
         end select
      end select
   enddo
   totTrajs=this%inicond%ntraj-this%inicond%nstart+1
   this%totScatt=totScatt
   this%totAllowedScatt=allowedScatt
   write(*,'("===========================================================")')
   write(*,'("ASSIGN PEAKS TO TRAJS: total trajs initialized: ",I10)')             totTrajs
   write(*,'("ASSIGN PEAKS TO TRAJS: scattered trajs: ",I10)')         this%totScatt
   write(*,'("ASSIGN PEAKS TO TRAJS: allowed scattered trajs: ",I10)') this%totAllowedScatt
   write(*,'("ASSIGN PEAKS TO TRAJS: probability: ",F10.5)')           dfloat(allowedScatt)/dfloat(totTrajs)
   write(*,'("===========================================================")')
   close(unit=ruScatt)
   close(unit=wuUnmap)
   call deleteIfEmptyFile( this%fileNameUnmapped )
   call this%doSubPeakStatistics()
   return
end subroutine assignTrajsToPeaks_ALLOWEDPEAKSGROW
!####################################################################################
! SUBROUTINE: doSubPeakStatistics_ALLOWEDPEAKSGROW
!####################################################################################
!------------------------------------------------------------------------------------
subroutine doSubPeakStatistics_ALLOWEDPEAKSGROW(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Allowed_peaksGrow),intent(inout):: this
   ! Local variables
   integer(kind=4):: nPeaks
   integer(kind=4):: i,j ! counters
   ! Run section ------------------------------------
   nPeaks=size( this%peaks(:) )
   do i=1,nPeaks
      this%peaks(i)%subPeaks(:)%prob=dfloat(this%peaks(i)%subPeaks(:)%ntrajs)/dfloat(this%totAllowedScatt)
      this%peaks(i)%subPeaks(:)%deviation=dsqrt(this%peaks(i)%subPeaks(:)%prob*(1.d0-this%peaks(i)%subPeaks(:)%prob)/dfloat(this%totAllowedScatt-1))
   enddo
   return
end subroutine doSubPeakStatistics_ALLOWEDPEAKSGROW
!####################################################################################
! FUNCTION: EVALUATE_PEAK ###########################################################
!####################################################################################
! - TRUE if A*(n^2) + B*nm + C*(m^2) + D*n + E*m + F < 0
!------------------------------------------------------------------------------------
function isAllowed_ALLOWEDPEAKSGROW(this,diffState,rovibrState,diffEnergy) result(isAllowed)
   IMPLICIT NONE
   ! I/O variables
   class(Allowed_peaksGROW),intent(in):: this
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
end function isAllowed_ALLOWEDPEAKSGROW
!###########################################################################3########
! SUBROUTINE: PRINT_XY_EXIT_ANGLES 
!###########################################################################3########
! - Reads data form "input_file" (should be dynamics-like output) and creates file "xy_exit_angle.out" with
!   information about exit angles in XY plane (taken from momenta information)
! - Only trajectories with "Scattered" status will be taken into account
!------------------------------------------------------------------------------------
subroutine PRINT_LABMOMENTA_AND_ANGLES_ALLOWEDPEAKSGROW(this)
   implicit none
   ! I/O Variables
   class(Allowed_peaksGROW),intent(inout):: this
   ! Local variables
   integer(kind=4):: i ! counters
   character(len=10):: stat
	real(kind=8),dimension(6):: p,r
   real(kind=8):: dpz
	real(kind=8),dimension(2):: dp ! variation of momentum
	real(kind=8),dimension(2):: dk ! variation of momentum in rec. space coord
   integer(kind=4),dimension(2):: g
	real(kind=8),dimension(3):: plab
   integer(kind=4):: ioErr
   integer(kind=4),dimension(3):: rovibrState
   character(len=4):: auxString
   character(len=1):: controlChar
	real(kind=8),dimension(12):: phaseSpaceVect
   real(kind=8):: psi,Theta,thetaout,beta
	real(kind=8):: gama,a,b,Etot,dE
   real(kind=8),dimension(2,2):: mtrx
   real(kind=8):: mtot, mu
   real(kind=8),dimension(2):: auxReal
   integer(kind=4),dimension(2):: auxInt
   integer(kind=4):: id
   integer(kind=4):: trajStat
   integer(kind=1):: controlSurfDyn
   real(kind=8):: dL,dEvibr,dErot,dEkin,dEpar,dEnorm
   ! Open units
   integer(kind=4),parameter:: wuFinal=12
   integer(kind=4),parameter:: ruScatt=11
   integer(kind=4),parameter:: ruStat=9
   integer(kind=4),parameter:: ruSurfDyn=8
   character(len=*),parameter:: formatFinal='(7(I4,1X),14(F10.5,1X))'
   character(len=*),parameter:: routinename = "print_labmomenta_and_angles_ALLOWEDPEAKSGROW6D: "
   ! RUN !! --------------------------
   beta=this%inicond%vpar_angle%getvalue()
   mtrx(1,:)=[dcos(beta),dsin(beta)]
   mtrx(2,:)=[-dsin(beta),dcos(beta)]
   a=system_surface%norm_s1
   b=system_surface%norm_s2
   mtot=sum( system_mass(:) )
   mu=product( system_mass(:) )/mtot
   gama=system_surface%angle
   ioErr=0

   open(unit=ruStat,file='OUT_INTERPSUMMARY',status='old',action='read')
   do while( ioErr == 0 )
      read(ruStat,*,iostat=ioErr) id,auxReal(1),auxInt(:),trajStat,auxReal(:),auxReal(:),controlChar
      if (ioErr==0) then
         select case(controlChar)
         case('X')
            ! as expected, yeah
         case default
            write(0,*) routinename//'ERR: Reading OUT_INTERPSUMMARY: need X as final character'
            call exit(1)
         end select
         select case(trajStat)
         case(12)
            this%inicond%trajs(id)%stat='Scattered'
         case(-1,-2,-12)
            this%inicond%trajs(id)%stat='React'
         case default
            this%inicond%trajs(id)%stat='Unknown'
         end select
      endif
   enddo
   close(ruStat)

   ioErr=0
   open(unit=ruSurfdyn,file='IN_SURFDYN',status='old',action='read',iostat=ioErr)
   select case(ioErr)
   case(0)
      read(ruSurfdyn,*)
      read(ruSurfdyn,*) controlSurfDyn
   case default
      controlSurfDyn=0   
   end select
   close(unit=ruSurfdyn,iostat=ioErr)

   open(unit=wuFinal,file="OUTANA6Dfinalpandangles.out",status="replace",action='write')
   write(wuFinal,'("# *********************************************** FINAL MOMENTA AND EXIT ANGLES ***************************************************")')
   write(wuFinal,'("# Format: id/diffOrder,n,m,v,J,mJ/dpx,dpy/dEpar,dEnorm,dEkin,dEvibr,dErot/dL/Ppar,Pperp,Pnorm(a.u.)/Azimuthal,Polar,Deflection(rad)")')
   write(wuFinal,'("# ---------------------------------------------------------------------------------------------------------------------------------")')
   open(unit=ruScatt,file="OUT_FINALCV",status="old",action='read')
   ioErr=0
   do while( ioErr == 0 )
      ! read from scattered file
      read(ruScatt,*,iostat=ioErr) auxString,auxString,auxString,auxString,id
      read(ruScatt,*,iostat=ioErr) r(1:3)
      read(ruScatt,*,iostat=ioErr) r(4:6)
      if(controlSurfDyn==1) read(ruScatt,*,iostat=ioErr)
      read(ruScatt,*,iostat=ioErr)
      read(ruScatt,*,iostat=ioErr) p(1:3)
      read(ruScatt,*,iostat=ioErr) p(4:6)
      if(controlSurfDyn==1) read(ruScatt,*,iostat=ioErr)
      ! correct momenta
      p(1:3)=system_mass(1)*p(1:3)/dsqrt(pmass2au)
      p(4:6)=system_mass(2)*p(4:6)/dsqrt(pmass2au)
      phaseSpaceVect(1:6)=r(:)
      phaseSpaceVect(7:12)=p(:)
      phaseSpaceVect(:)=from_atomic_to_molecular_phaseSpace(atomcoord=phaseSpaceVect)
      r(:)=phaseSpaceVect(1:6)
      p(:)=phaseSpaceVect(7:12)

      select case( ioErr==0 .and. this%inicond%trajs(id)%stat=='Scattered' )
      case(.true.) ! secure to operate
         dp(:) = p(1:2)-this%inicond%trajs(id)%init_p(1:2)
         dL=dsqrt(p(5)**2.d0+(p(6)/dsin(r(5)))**2.d0)-dsqrt(&
            this%inicond%trajs(id)%init_p(5)**2.d0&
            +this%inicond%trajs(id)%init_p(6)**2.d0&
            *dsin(this%inicond%trajs(id)%init_r(5))**(-2.d0))
         g(:) = this%quantizeDiffState( dp )
         dEvibr=this%inicond%vibrPot%getPot( r(4) )+(0.5d0/mu)*p(4)**2.d0&
                -this%inicond%vibrPot%getPot( this%inicond%trajs(id)%init_r(4))-(0.5d0/mu)*this%inicond%trajs(id)%init_p(4)**2.d0
         dErot=(0.5d0/mu)*(p(5)**2.d0+(p(6)/dsin(r(5)))**2.d0)/r(4)**2.d0&
              -(0.5d0/mu)*(this%inicond%trajs(id)%init_p(5)**2.d0+(this%inicond%trajs(id)%init_p(6)/dsin(this%inicond%trajs(id)%init_r(5)))**2.d0)/this%inicond%trajs(id)%init_r(4)**2.d0
         dEkin=(0.5d0/mtot)*dot_product(p(1:3),p(1:3))-(0.5d0/mtot)*dot_product(this%inicond%trajs(id)%init_p(1:3),this%inicond%trajs(id)%init_p(1:3))
         dEpar=(0.5d0/mtot)*dot_product(p(1:2),p(1:2))-(0.5d0/mtot)*dot_product(this%inicond%trajs(id)%init_p(1:2),this%inicond%trajs(id)%init_p(1:2))
         dEnorm=(0.5d0/mtot)*p(3)**2.d0-(0.5d0/mtot)*this%inicond%trajs(id)%init_p(3)**2.d0
         Etot=0.5d0*dot_product(p(1:3),p(1:3))/mtot+&                      ! Kinetic energy
              0.5d0*(p(4)**2.d0)/mu+&                                      ! Internal kinetic energy (1)
              0.5d0*(p(5)**2.d0+(p(6)/dsin(r(5)))**2.d0)/(mu*r(4)**2.d0)+& ! Internal kinetic energy (2)
              this%inicond%vibrPot%getPot( r(4) )                          ! Potential energy
         rovibrState=this%quantizeRovibrState( Etot=Etot,position=r(:),momenta=p(:) )
         plab(1:2)=matmul(mtrx,p(1:2))
         plab(3)=p(3)
         psi = datan2(plab(2),plab(1))
         Theta = datan(plab(2)/plab(3))
         thetaout=datan(plab(3)/norm2(plab(1:2)))
         write(wuFinal,formatFinal) id,getDiffOrder(g=g),g(:),rovibrState(:),dp(1:2),dEpar,dEnorm,dEkin,dEvibr,dErot,dL,plab(:),psi,Theta,thetaout
         
      case(.false.) ! not secure to operate, next switch
         select case( ioErr )
         case(-1,0,5001)
            ! do nothing, EOF reached, let it break the cycle
         case default
            write(0,*) routinename//'ERR: unexpected error encountered while reading OUT_FINALCV. Error Code: ',ioErr
            call exit(1)
         end select
      end select
 
   enddo
   close(unit=wuFinal)
   close(unit=ruScatt)
   return
end subroutine PRINT_LABMOMENTA_AND_ANGLES_ALLOWEDPEAKSGROW
!###########################################################
!# FUNCTION: getDiffOrder
!###########################################################
!> @brief
!! Given environment order and peak labels, sets real diffraction
!! order.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Nov/2014
!> @version 1.0
!-----------------------------------------------------------
function getDiffOrder(envOrder,g) result(diffOrder)
   ! Initial declarations
   implicit none
   ! I/O variables
   integer(kind=4),optional,intent(in) :: envOrder
   integer(kind=4),dimension(2),intent(in) :: g
   ! Function dummy variable
   integer(kind=4):: diffOrder
   ! Run section
   if( system_surface%order==4 .and. present(envOrder) ) then
      if( all( g(:)==[0,0] ) ) then
         diffOrder=envOrder
      elseif(g(1)==0 .OR. g(2)==0) then
         diffOrder=envOrder*(envOrder+1)/2
      elseif(abs(g(1))==abs(g(2))) then
         diffOrder=((envOrder+1)*(envOrder+2)/2)-1
      elseif(abs(g(1))==envOrder) then
         diffOrder=(envOrder*(envOrder+1)/2)+abs(g(2))
      elseif(abs(g(2))==envOrder) then
         diffOrder=(envOrder*(envOrder+1)/2)+abs(g(1))
      else
         write(0,*) 'ERR getDiffOrder: Unclassificable difraction peak. Surface order: 4'
         call exit(1)
      endif
   elseif( system_surface%order==6 .and. dcos(system_surface%angle)>0.d0 ) then
      if( all( g(:)==[0,0] ) ) then
         diffOrder=0
      elseif( g(1)==0 .or. g(2)==0 ) then
         diffOrder=checkLoschianOrder( sum(g)**2 )
      elseif( g(1)==g(2) ) then
         diffOrder=checkLoschianOrder( g(1)**2 )
      elseif( g(1)==-g(2) ) then
         diffOrder=checkLoschianOrder( 3*g(1)**2 )
      elseif( sign( 1,g(1) )*sign( 1,g(2) ) > 0 ) then
         diffOrder=checkLoschianOrder( g(1)**2+g(2)**2-abs(product(g)) )
      elseif( sign( 1,g(1) )*sign( 1,g(2) ) < 0 ) then
         diffOrder=checkLoschianOrder( g(1)**2+g(2)**2+abs(product(g)) )
      else
         write(0,*) 'ERR getDiffOrder: Unclassificable difraction peak. Surface order: 6'
         call exit(1)
      endif
   elseif( system_surface%order==6 .and. dcos(system_surface%angle)<0.d0 ) then
      if( all( g(:)==[0,0] ) ) then
         diffOrder=0
      elseif( g(1)==0 .or. g(2)==0 ) then
         diffOrder=checkLoschianOrder( sum(g)**2 )
      elseif( g(1)==g(2) ) then
         diffOrder=checkLoschianOrder( 3*g(1)**2 )
      elseif( g(1)==-g(2) ) then
         diffOrder=checkLoschianOrder( g(1)**2 )
      elseif( sign( 1,g(1) )*sign( 1,g(2) ) > 0 ) then
         diffOrder=checkLoschianOrder( g(1)**2+g(2)**2+abs(product(g)) )
      elseif( sign( 1,g(1) )*sign( 1,g(2) ) < 0 ) then
         diffOrder=checkLoschianOrder( g(1)**2+g(2)**2-abs(product(g)) )
      else
         write(0,*) 'ERR getDiffOrder: Unclassificable difraction peak. Surface order: 6'
         call exit(1)
      endif
   else
      write(0,*) 'ERR getDiffOrder: wrong surface order or bar usage of this routine'
      call exit(1)
   endif
   return
end function getDiffOrder
!######################################################
! FUNCTION: getPeakId_ALLOWEDPEAKSGROW
!######################################################
!> @brief
!! - Gets allowed peak ID number given its diffraction numbers.
!! - If there is not an allowed peak with diffraction number g(:),
!!   this function exits with a negative integer
!------------------------------------------------------
function getPeakId_ALLOWEDPEAKSGROW(this,g) result(peakId)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Allowed_peaksGROW),intent(in):: this
   integer(kind=4),dimension(2):: g
   ! Function dummy variable
   integer(kind=4):: peakId
   ! Local variables
   integer(kind=4):: i ! counter
   ! Run section
   select case( .not.allocated(this%peaks) )
   case(.true.) ! stop! badness!
      write(0,*) 'getPeakId_ALLOWEDPEAKSGROW ERR: Peaks are not allocated'
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
end function getPeakId_ALLOWEDPEAKSGROW
!############################################################
! SUBROUTINE: addNewPeak_PEAKGROW
!############################################################
!------------------------------------------------------------
subroutine addNewPeak_ALLOWEDPEAKSGROW(this,g,peakId)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Allowed_peaksGROW),intent(inout):: this
   integer(kind=4),dimension(2),intent(in):: g
   integer(kind=4),optional,intent(out):: peakId
   ! Local variables
   type(PeakGROW),dimension(:),allocatable:: auxListPeaks
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
   this%peaks(idNew)%diffOrder=getDiffOrder( this%peaks(idNew)%envOrder,g(:) )
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
end subroutine addNewPeak_ALLOWEDPEAKSGROW
!############################################################
! SUBROUTINE: addCountToPeak_ALLOWEDPEAKSGROW
!############################################################
!> @brief
!! Adds probability to a given subpeak. If it does not exist, this
!! routine will initialize it. Subpeaks are classified by their quantum
!! state.
!------------------------------------------------------------
subroutine addCountToSubPeak_ALLOWEDPEAKSGROW(this,peakId,rovibrState)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Allowed_peaksGROW),intent(inout):: this
   integer(kind=4),intent(in):: peakId
   integer(kind=4),dimension(3),intent(in):: rovibrState
   ! Local variables
   real(kind=8):: initE,kz,masa,theta_in,pinit_par
   integer(kind=4):: newCol,oldCol
   type(subPeakGROW),dimension(:),allocatable:: auxListSubPeaks
   logical:: isNew
   integer(kind=4):: i ! counter
   ! Run section
   ! Some
   masa=sum(system_mass(:))
   initE=this%inicond%evirot%getValue()
	theta_in=this%inicond%vz_angle%getvalue()
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
      this%peaks(peakId)%subPeaks(newCol)%phiOut=datan( this%peaks(peakId)%dky/kz )
      this%peaks(peakId)%subPeaks(newCol)%thetaOut=datan( kz/norm2([this%peaks(peakId)%dkx+pinit_par,this%peaks(peakId)%dky]) )
      this%peaks(peakId)%subPeaks(newCol)%ntrajs=1
   case(.false.)
      this%peaks(peakId)%subPeaks(i-1)%ntrajs=this%peaks(peakId)%subPeaks(i-1)%ntrajs+1
   end select

   return
end subroutine addCountToSubPeak_ALLOWEDPEAKSGROW
!#######################################################
! SUBROUTINE: printAllowedPeaks_ALLOWEDPEAKSGROW
!#######################################################
!> @brief
!! Prints all allowed peaks.
!-------------------------------------------------------
subroutine printAllowedPeaks_ALLOWEDPEAKSGROW(this)
   ! initial declarations
   implicit none
   ! I/O variables
   class(Allowed_PeaksGROW),intent(in):: this
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
end subroutine printAllowedPeaks_ALLOWEDPEAKSGROW
!########################################################
! SUBROUTINE: sortByDiffOrder_ALLOWEDGROW
!########################################################
! @brief
!! sort peaks by diffraction order and not by peak id
!--------------------------------------------------------
subroutine sortByDiffOrder_ALLOWEDPEAKSGROW(this)
   ! initial declarations
   implicit none
   ! I/O variables
   class(Allowed_PeaksGROW),intent(inout):: this
   ! Local variables
   integer(kind=4),dimension(:),allocatable:: intList
   integer(kind=4):: N
   integer(kind=4):: i,j,diffOrder ! counters
   type(peakGROW),dimension(:),allocatable:: auxPeaksList
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
end subroutine sortByDiffOrder_ALLOWEDPEAKSGROW
!###########################################################
! SUBROUTINE: printSeenPeaks_ALLOWEDPEAKSGROW
!###########################################################
!> @brief
!! Prints Seen peaks file
!-----------------------------------------------------------
subroutine printSeenPeaks_ALLOWEDPEAKSGROW(this)
   ! initial declarations
   implicit none
   ! I/O variables
   class(Allowed_PeaksGROW),intent(in):: this
   ! Local variables
   integer(kind=4),parameter:: wuSeen=123
   character(len=*),parameter:: formatSeen='(6(I10,1X),6(F10.5,1X),F10.5)'
   integer(kind=4):: i,j ! counters
   ! Run section --------------------------------------------
   open(unit=wuSeen,file=this%fileNameSeen,status='replace',action='write')
   write(wuSeen,'("# --------------------------------------- SEEN PEAKS ----------------------------------------------------------------")')
   write(wuSeen,'("# Initialized trajectories: ",I15)') size( this%inicond%trajs(:) )
   write(wuSeen,'("# Total scattered trajectories: ",I15)') this%totScatt
   write(wuSeen,'("# Total physical scattered trajectories (probabilities based on this number): ",I15)')  this%totAllowedScatt
   write(wuSeen,'("# Format: diffOrder/n,m,V,J,mJ/Azimuthal angle/Deflection from perp. surf./Deflection from surf./dKxy(au)/dE(au)/prob/deviat")')
   write(wuSeen,'("# -------------------------------------------------------------------------------------------------------------------")')
   do i=1,size( this%peaks )
      select case( allocated(this%peaks(i)%subPeaks) )
      case(.true.)
         do j =1,size( this%peaks(i)%subPeaks )
            write(wuSeen,formatSeen) this%peaks(i)%diffOrder,this%peaks(i)%g(:),this%peaks(i)%subPeaks(j)%rovibrState(:),&
                                     this%peaks(i)%psiOut,this%peaks(i)%subPeaks(j)%PhiOut,this%peaks(i)%subPeaks(j)%thetaOut,&
                                     this%peaks(i)%dkxy,this%peaks(i)%subPeaks(j)%dE,this%peaks(i)%subPeaks(j)%prob,&
                                     this%peaks(i)%subPeaks(j)%deviation
         enddo
      case(.false.)
         ! do nothing. do not print anything
      end select
   enddo
   close(unit=wuSeen)
   return
end subroutine printSeenPeaks_ALLOWEDPEAKSGROW
!###########################################################
! FUNCTION: quantizeRovibrState_ALLOWEDPEAKSGROW
!###########################################################
!-----------------------------------------------------------
function quantizeRovibrState_ALLOWEDPEAKSGROW(this,Etot,position,momenta) result(rovibrState)
   ! initial declarations
   implicit none
   ! I/O variables
   class(Allowed_peaksGROW),intent(in):: this
   real(kind=8),intent(in):: Etot
   real(kind=8),dimension(6),intent(in):: position
   real(kind=8),dimension(6),intent(in):: momenta
   ! Dummy function variable
   integer(kind=4),dimension(3):: rovibrState
   ! Local variables
   real(kind=8):: L2,mu,masa,Ecm,Evibr,Erot,finalJ,finalV
   ! Parameters
   character(len=*),parameter:: routinename='quantizeRovibrState_ALLOWEDPEAKSGROW: '
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
   case( 'Morse' )
      finalJ=(-1.d0+dsqrt(1.d0+4.d0*L2))*0.5d0
      finalV=dsqrt(1.d0-Evibr/system_binningParam(1))+dsqrt(2.d0*system_binningParam(1)/masa)/system_binningParam(3)-0.5d0
   case( 'Harmonic' )
      finalJ=(-1.d0+dsqrt(1.d0+4.d0*L2))*0.5d0
      finalV=dsqrt(mu/system_binningParam(1))*Evibr-0.5d0
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

end function quantizeRovibrState_ALLOWEDPEAKSGROW
!###################################################################
! FUNCTION: quantizeDiffState_ALLOWEDPEAKSGROW
!###################################################################
function quantizeDiffState_ALLOWEDPEAKSGROW(this,p) result(g)
   implicit none
   ! I/O variables
   class(Allowed_peaksGROW),intent(in):: this
   real(kind=8),dimension(2),intent(in):: p
   ! Dummy output variables
   integer(kind=4),dimension(2):: g
   ! Local variables
   real(kind=8),dimension(2):: auxVect
   integer(kind=4):: i
   integer(kind=4),dimension(7,2):: c
   real(kind=8),dimension(7,2):: a
   real(kind=8),dimension(7):: dist
   integer(kind=4),dimension(1):: auxInt
   ! Run section ..................................
   auxVect(:)=system_surface%cart2recip( p )
   select case( system_surface%order )
   case(4)
      g(:)=nint( auxVect(:) )
   case(6)
      c(1,:)=nint( auxVect(:) )
      select case( dcos(system_surface%angle)>0.d0 )
      case(.true.)
         c(2,:)=c(1,:)+[1,0]
         c(3,:)=c(1,:)-[1,0]
         c(4,:)=c(1,:)+[0,1]
         c(5,:)=c(1,:)-[0,1]
         c(6,:)=c(1,:)+[1,1]
         c(7,:)=c(1,:)-[1,1]
      case(.false.)
         c(2,:)=c(1,:)+[1,0]
         c(3,:)=c(1,:)-[1,0]
         c(4,:)=c(1,:)+[0,1]
         c(5,:)=c(1,:)-[0,1]
         c(6,:)=c(1,:)+[1,-1]
         c(7,:)=c(1,:)-[1,-1]
      end select
      do i=1,7
         a(i,:)=system_surface%recip2cart( dfloat( c(i,:) ) )
         dist(i)=norm2( p(:)-a(i,:) )
      enddo
      auxInt(:)=minloc( array=dist(:) )
      g(:)=c(auxInt(1),:)
   end select
   return
end function quantizeDiffState_ALLOWEDPEAKSGROW
END MODULE DIFFRACTIONGROW_MOD
