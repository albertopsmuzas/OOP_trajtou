!##############################################################################
! MODULE: SYSTEM_MOD
!> @brief
!! Public module which constains interesting variables and subprobrams to
!! keep during runtime. Debug variables are separated inside module DEBUG_MOD
!##############################################################################
MODULE SYSTEM_MOD
! Initial declarations
use UNITS_MOD
use SURFACE_MOD
use MATHS_MOD, only: INV_MTRX
use AOTUS_MODULE, only: flu_State, OPEN_CONFIG_FILE, CLOSE_CONFIG, AOT_GET_VAL
use AOT_TABLE_MODULE, only: AOT_TABLE_OPEN, AOT_TABLE_CLOSE, AOT_TABLE_LENGTH, AOT_TABLE_GET_VAL
#ifdef DEBUG
   use DEBUG_MOD
#endif
#ifdef INTEL
   use IFPORT, only:getPid ! for ifort compatibility
#endif
IMPLICIT NONE
! Global variables
REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: system_mass
type(Surface):: system_surface
real(kind=8):: system_surfMass=0.d0
real(kind=8),dimension(3):: system_surfFreqs=[0.d0,0.d0,0.d0]
CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE:: system_atomsymbols
CHARACTER(LEN=:),ALLOCATABLE:: system_inputfile
CHARACTER(LEN=:),ALLOCATABLE:: system_pespath
character(len=:),allocatable:: system_binningScheme
real(kind=8),dimension(:),allocatable:: system_binningParam
integer(kind=4):: system_binningdJ
INTEGER(KIND=4):: system_natoms
integer(kind=4),dimension(:),allocatable:: system_iSeed

! Global parameters
character(len=*),parameter:: system_seedFilename='INseed.inp'
! Contains section
CONTAINS
!########################################################################
!# SUBROUTINE: generate_seed
!########################################################################
!> @brief
!! Generates a good seed and stores it in system_iSeed
!
!> @details
!! Taken from https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
!! Seems to be quite secure for this job. Feel free to implement a better seed
!! generation algorithm.
!
!> @param[out] iStat - integer(kind=1),optional: if given, it'll store
!!                     the status of the seed generation:
!!                     - 0) read from default seed input file
!!                     - 1) read from /dev/urandom
!!                     - 2) generated from system_clock
!------------------------------------------------------------------------
subroutine generate_seed(iStat)
   ! initial declarations
   implicit none
   ! I/O variables
   integer(kind=1),intent(out),optional:: iStat
   ! Local variables
   logical:: fileExists
   integer(kind=4):: seedLength
   integer(kind=4):: i ! counter
   integer(kind=1):: iErr
   integer(kind=8):: t
   integer(kind=4),dimension(8):: dt
   integer(kind=4):: pid
   ! R/W units and parameters
   integer(kind=1),parameter:: ruSeed=17
   integer(kind=1),parameter:: wuSeed=18
   character(len=*),parameter:: routinename='GENERATE_SEED: '
   ! Run section
   select case(allocated(system_iSeed))
      case(.true.)  ! bad news
         write(0,*) 'ERR: '//routinename//'System iseed was already allocated'
         call exit(1)
      case(.false.) ! allocate seed with proper length
         call random_seed(size=seedLength)
         allocate(system_iSeed(seedLength))
   end select
   inquire(file=system_seedFilename,exist=fileExists)
   select case(fileExists)
      case(.true.)
         open(unit=ruSeed,file=system_seedFilename,status='old',action='read')
         read(ruSeed,*,iostat=iErr) system_iSeed(:) ! we are assuming that this is a good seed. No checks
         close(ruSeed)
         select case(iErr)
            case(0)
               ! do nothing
            case default
               write(0,*) 'ERR: '//routinename//'Unexpected error while reading '//system_seedFilename//' file'
               write(0,*) 'IOSTAT error code: ',iErr
               call exit(1)
         end select
         if(present(iStat)) iStat=0
      case(.false.)
         ! try to get random iSeed from computer
         open(unit=ruSeed,file='/dev/urandom',access='stream',form='unformatted',action='read',status='old',iostat=iErr)
         select case(iErr)
            case(0)
               read(ruSeed) system_iSeed(:)
               close(ruSeed)
               if(present(iStat)) iStat=1

            case default
               call system_clock(t)
               select case(t==0)
                  case(.true.)
                     call date_and_time(values=dt)
                     t = (dt(1)-1970)*365_8*24*60*60*1000 &
                         +dt(2)*31_8*24*60*60*1000 &
                         +dt(3)*24_8*60*60*1000 &
                         +dt(5)*60*60*1000 &
                         +dt(6)*60*1000+dt(7)*1000+dt(8)
                     if(present(iStat)) iStat=2

                  case(.false.)
                     ! do nothing
               end select
               pid=getpid()
               t=ieor(t,int(pid,kind(t)))
               do i=1,seedLength
                  system_iSeed(i)=lcg(t)
               end do

         end select
               ! Store used seed
               open(unit=wuSeed,file=system_seedFilename,status='new',action='write')
               write(wuSeed,*) system_iSeed(:)
               close(wuSeed)
   end select
   return
   contains
   ! seeding function
   function lcg(iFeed) result(iSeed)
      implicit none
      integer(kind=8),intent(inout):: iFeed
      integer(kind=4):: iSeed
      select case(iFeed)
         case(0) ! badneeeessss
            iFeed=104729_8
         case default
            iFeed=mod(iFeed,4294967296_8)
      end select
      iFeed=mod(iFeed * 279470273_8, 4294967291_8)
      iSeed = int(mod(iFeed, int(huge(0),kind=8)), kind=4)
   end function lcg
end subroutine generate_seed
!###########################################################
!# SUBROUTINE: INITIALIZE_SYSTEM
!###########################################################
!> @brief
!! Load system parameters
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INITIALIZE_SYSTEM(filename)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CHARACTER(LEN=*),INTENT(IN):: filename
   ! Local variables
   TYPE(flu_State):: conf ! Lua file
   INTEGER(KIND=4):: ierr
   INTEGER(KIND=4):: sys_table,sym_table,magnitude_table,osciSurf_table,binning_table
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE:: subtable
   CHARACTER(LEN=*),PARAMETER:: routinename="INITIALIZE_SYSTEM: "
   TYPE(Mass):: masa
   REAL(KIND=8):: numero
   CHARACTER(LEN=10):: units
   INTEGER(KIND=4):: i ! counters
   CHARACTER(LEN=1024):: auxString
   real(kind=8):: auxReal
   LOGICAL:: auxbool
   ! Run section
   ALLOCATE(system_inputfile,source=filename)
   CALL OPEN_CONFIG_FILE(L=conf,filename=system_inputfile,ErrCode=ierr)
   SELECT CASE(ierr)
      CASE(0)
         ! do nothing
      CASE DEFAULT
         WRITE(0,*) "INITIALIZE_SYSTEM: error opening Lua config file"
         CALL EXIT(1)
   END SELECT
   ! Open table system
   CALL AOT_TABLE_OPEN(L=conf,thandle=sys_table,key='system')
#ifdef DEBUG
   SELECT CASE(get_debugmode())
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=sys_table,key='debugMode',val=auxbool)
         CALL SET_DEBUG_MODE(auxbool)
   END SELECT
   SELECT CASE(get_verbosemode())
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=sys_table,key='verboseMode',val=auxbool)
         CALL SET_VERBOSE_MODE(auxbool)
   END SELECT
#endif
   ! Get surface
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=sys_table,key='surface',val=auxString)
   CALL system_surface%INITIALIZE(trim(auxString))
   ! Get number of atoms, masses and atomic symbols
   CALL AOT_TABLE_OPEN(L=conf,parent=sys_table,thandle=sym_table,key='symbols')
   CALL AOT_TABLE_OPEN(L=conf,parent=sys_table,thandle=magnitude_table,key='masses')
   system_natoms=aot_table_length(L=conf,thandle=sym_table)
   SELECT CASE(system_natoms==aot_table_length(L=conf,thandle=magnitude_table) .AND. system_natoms/=0)
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         WRITE(0,*) "INITIALIZE_SYSTEM ERR: dimensions mismatch (or zero) between tables system.masses and system.symbols"
         CALL EXIT(1)
   END SELECT
   ALLOCATE(system_mass(system_natoms))
   ALLOCATE(system_atomsymbols(system_natoms))
   ALLOCATE(subtable(system_natoms))
   DO i = 1, system_natoms
      CALL AOT_TABLE_OPEN(L=conf,parent=magnitude_table,thandle=subtable(i),pos=i)
      CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=subtable(i),pos=1,val=numero)
      CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=subtable(i),pos=2,val=units)
      CALL masa%READ(numero,trim(units))
      CALL masa%TO_STD()
      system_mass(i)=masa%getvalue()
      CALL AOT_TABLE_CLOSE(L=conf,thandle=subtable(i))
      CALl AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=sym_table,pos=i,val=system_atomsymbols(i))
   END DO
   CALL AOT_TABLE_CLOSE(L=conf,thandle=magnitude_table)
   CALL AOT_TABLE_CLOSE(L=conf,thandle=sym_table)
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=sys_table,key='pesPath',val=auxString)
   system_pespath=trim(auxString)
   ! get Surface oscillator parameters
   call aot_table_open(L=conf,parent=sys_table,thandle=osciSurf_table,key='surfaceOscillator')
   call aot_get_val(L=conf,ErrCode=ierr,thandle=osciSurf_table,key='wx',val=system_surfFreqs(1))
   call aot_get_val(L=conf,ErrCode=ierr,thandle=osciSurf_table,key='wy',val=system_surfFreqs(2))
   call aot_get_val(L=conf,ErrCode=ierr,thandle=osciSurf_table,key='wz',val=system_surfFreqs(3))
   call aot_table_open(L=conf,parent=osciSurf_table,thandle=magnitude_table,key='mass')
   call aot_get_val(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=1,val=auxReal)
   call aot_get_val(L=conf,ErrCode=ierr,thandle=magnitude_table,pos=2,val=auxString)
   call masa%READ(auxReal,trim(auxString))
   call masa%to_std()
   system_surfMass=masa%getValue()
   call aot_table_close(L=conf,thandle=magnitude_table)
   call aot_table_close(L=conf,thandle=osciSurf_table)
   ! get binning parameters
   call aot_table_open(L=conf,parent=sys_table,thandle=binning_table,key='binning')
   call aot_get_val(L=conf,ErrCode=iErr,thandle=binning_table,key='kind',val=auxString)
   system_binningScheme=trim(auxString)
   if( system_binningScheme=='Morse' ) then
      allocate( system_binningParam(3) )
      call aot_get_val(L=conf,ErrCode=iErr,thandle=binning_table,key='dissociationEnergy',val=system_binningParam(1))
      call aot_get_val(L=conf,ErrCode=iErr,thandle=binning_table,key='equilibriumDistance',val=system_binningParam(2))
      call aot_get_val(L=conf,ErrCode=iErr,thandle=binning_table,key='width',val=system_binningParam(3))
      call aot_get_val(L=conf,ErrCode=iErr,thandle=binning_table,key='dJ',val=system_binningdJ)
#ifdef DEBUG
      call verbose_write( routinename,'Kind of data binning: '//system_binningScheme )
      call verbose_write( routinename,'Dissociation energy(au): ',system_binningParam(1) )
      call verbose_write( routinename,'Equilibrium Distance(au): ',system_binningParam(2) )
      call verbose_write( routinename,'Morse width parameter(au): ',system_binningParam(3) )
      call verbose_write( routinename,'dJ: ',system_binningdJ )
#endif
   else
      write(0,*) 'INITIALIZE SYSTEM ERR: wrong binning scheme type: '//system_binningScheme
      write(0,*) 'Implemented ones: Morse'
      write(0,*) 'Case sensitive'
      call exit(1)
   endif
   call aot_table_close(L=conf,thandle=binning_table)
   ! CLOSE LUA FILE /////////
   CALL CLOSE_CONFIG(conf)
   ! ////////////////////////
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'config file: ',system_inputfile)
   CALL VERBOSE_WRITE(routinename,'default surface input file: ',system_surface%getfilename())
   CALL VERBOSE_WRITE(routinename,'default PES path: ',system_pespath)
   CALL VERBOSE_WRITE(routinename,'Natoms: ',system_natoms)
   CALL VERBOSE_WRITE(routinename,'masses(au): ',system_mass(:))
   CALL VERBOSE_WRITE(routinename,'atomic symbols: ',system_atomsymbols(:))
   call verbose_write(routinename,'Surface mass (au): ',system_surfMass)
   call verbose_write(routinename,'Surface frequencies (au): ',system_surfFreqs(:))
#endif
   RETURN
END SUBROUTINE INITIALIZE_SYSTEM
!###########################################################
! FUNCTION: from_molecular_to_atomic
!###########################################################
!> @brief
!! Go from molecular coordinates x,y,z,r,theta,phi to
!! xa,ya,za,xb,yb,zb
!-----------------------------------------------------------
PURE FUNCTION from_molecular_to_atomic(molcoord) result(atomcoord)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   REAL(KIND=8),DIMENSION(6),INTENT(IN):: molcoord
   ! Dymmy function variable
   REAL(KIND=8),DIMENSION(6):: atomcoord
   ! Local
   REAL(KIND=8):: masa
   ! Run section
   masa=sum(system_mass(:))
   atomcoord(1)=molcoord(1)+(system_mass(2)/(masa))*molcoord(4)*dcos(molcoord(6))*dsin(molcoord(5))
   atomcoord(2)=molcoord(2)+(system_mass(2)/(masa))*molcoord(4)*dsin(molcoord(6))*dsin(molcoord(5))
   atomcoord(3)=molcoord(3)+(system_mass(2)/(masa))*molcoord(4)*dcos(molcoord(5))
   atomcoord(4)=molcoord(1)-(system_mass(1)/(masa))*molcoord(4)*dcos(molcoord(6))*dsin(molcoord(5))
   atomcoord(5)=molcoord(2)-(system_mass(1)/(masa))*molcoord(4)*dsin(molcoord(6))*dsin(molcoord(5))
   atomcoord(6)=molcoord(3)-(system_mass(1)/(masa))*molcoord(4)*dcos(molcoord(5))
   RETURN
END FUNCTION from_molecular_to_atomic

FUNCTION from_molecular_to_atomic_phaseSpace(molcoord) result(atomcoord)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   REAL(KIND=8),DIMENSION(12),INTENT(IN):: molcoord
   ! Dymmy function variable
   REAL(KIND=8),DIMENSION(12):: atomcoord
   ! Local
   REAL(KIND=8):: masa
   REAL(KIND=8):: nua,nub
   REAL(KIND=8),DIMENSION(6,6):: mtrx
   DATA mtrx(:,1)/1.d0,0.d0,0.d0,1.d0,0.d0,0.d0/
   DATA mtrx(:,2)/0.d0,1.d0,0.d0,0.d0,1.d0,0.d0/
   DATA mtrx(:,3)/0.d0,0.d0,1.d0,0.d0,0.d0,1.d0/
   DATA mtrx(:,6)/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/
   ! Run section
   masa=sum(system_mass(:))
   nua=system_mass(1)/masa
   nub=system_mass(2)/masa
   ! Transformation of position coordinates
   atomcoord(1)=molcoord(1)+(system_mass(2)/(masa))*molcoord(4)*dcos(molcoord(6))*dsin(molcoord(5))
   atomcoord(2)=molcoord(2)+(system_mass(2)/(masa))*molcoord(4)*dsin(molcoord(6))*dsin(molcoord(5))
   atomcoord(3)=molcoord(3)+(system_mass(2)/(masa))*molcoord(4)*dcos(molcoord(5))
   atomcoord(4)=molcoord(1)-(system_mass(1)/(masa))*molcoord(4)*dcos(molcoord(6))*dsin(molcoord(5))
   atomcoord(5)=molcoord(2)-(system_mass(1)/(masa))*molcoord(4)*dsin(molcoord(6))*dsin(molcoord(5))
   atomcoord(6)=molcoord(3)-(system_mass(1)/(masa))*molcoord(4)*dcos(molcoord(5))
   ! Transformation matrix for momenta
   mtrx(1,4)=dsin(molcoord(5))*dcos(molcoord(6))/nua
   mtrx(2,4)=dsin(molcoord(5))*dsin(molcoord(6))/nua
   mtrx(3,4)=dcos(molcoord(5))/nua
   mtrx(4,4)=-dsin(molcoord(5))*dcos(molcoord(6))/nub
   mtrx(5,4)=-dsin(molcoord(5))*dsin(molcoord(6))/nub
   mtrx(6,4)=-dcos(molcoord(5))/nub
   mtrx(1,5)=dcos(molcoord(5))*dcos(molcoord(6))/(nua*molcoord(4))
   mtrx(2,5)=dcos(molcoord(5))*dsin(molcoord(6))/(nua*molcoord(4))
   mtrx(3,5)=dsin(molcoord(5))/(nua*molcoord(4))
   mtrx(4,5)=-dcos(molcoord(5))*dcos(molcoord(6))/(nub*molcoord(4))
   mtrx(5,5)=-dcos(molcoord(5))*dsin(molcoord(6))/(nub*molcoord(4))
   mtrx(6,5)=-dsin(molcoord(5))/(nub*molcoord(4))
   SELECT CASE(dsin(molcoord(5)) /= 0.d0)
       CASE(.true.)
          mtrx(1,6)=-dsin(molcoord(6))/(nua*dsin(molcoord(5)))
          mtrx(2,6)=dcos(molcoord(6))/(nua*dsin(molcoord(5)))
          mtrx(3,6)=0.d0
          mtrx(4,6)=dsin(molcoord(6))/(nub*dsin(molcoord(5)))
          mtrx(5,6)=-dcos(molcoord(6))/(nub*dsin(molcoord(5)))
          mtrx(6,6)=0.d0
       CASE(.false.)
          ! do nothing
   END SELECT
   atomcoord(7:12)=matmul(mtrx,molcoord(7:12))
   RETURN
END FUNCTION from_molecular_to_atomic_phaseSpace
!###########################################################################
! FUNCTION: from_atomic_to_molecular
!###########################################################################
!> @brief
!! Goes from atomic coordinates xa,ya,za,xb,yb,zb to molecular
!! coordinates x,y,z,r,theta,phi.
!> @details
!! - We have enforced @f$\theta \in [0,\pi]@f$ and @f$\phi \in [0,2\pi)@f$
!--------------------------------------------------------------------------
PURE FUNCTION from_atomic_to_molecular(atomcoord) result(molcoord)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   REAL(KIND=8),DIMENSION(6),INTENT(IN):: atomcoord
   ! Dummy function variable
   REAL(KIND=8),DIMENSION(6):: molcoord
   ! Local variables
   REAL(KIND=8):: masa
   ! Run section
   masa=sum(system_mass(:))
   molcoord(1)=(1.D0/(masa))*(atomcoord(1)*system_mass(1)+atomcoord(4)*system_mass(2))
   molcoord(2)=(1.D0/(masa))*(atomcoord(2)*system_mass(1)+atomcoord(5)*system_mass(2))
   molcoord(3)=(1.D0/(masa))*(atomcoord(3)*system_mass(1)+atomcoord(6)*system_mass(2))
   molcoord(4)=dsqrt((atomcoord(1)-atomcoord(4))**2.D0+&
      (atomcoord(2)-atomcoord(5))**2.D0+(atomcoord(3)-atomcoord(6))**2.D0)
   molcoord(5)=dacos((atomcoord(3)-atomcoord(6))/molcoord(4))
   molcoord(6)=datan2((atomcoord(2)-atomcoord(5)),(atomcoord(1)-atomcoord(4)))
   SELECT CASE(molcoord(6)<0.d0)
      CASE(.true.)
         molcoord(6)=molcoord(6)+2.d0*pi
      CASE(.false.)
         ! do nothing
   END SELECT
END FUNCTION from_atomic_to_molecular

FUNCTION from_atomic_to_molecular_phaseSpace(atomcoord) result(molcoord)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   REAL(KIND=8),DIMENSION(12),INTENT(IN):: atomcoord
   ! Dummy function variable
   REAL(KIND=8),DIMENSION(12):: molcoord
   ! Local variables
   REAL(KIND=8):: masa,nua,nub
   REAL(KIND=8),DIMENSION(6,6):: mtrx, invMtrx
   DATA mtrx(:,1)/1.d0,0.d0,0.d0,1.d0,0.d0,0.d0/
   DATA mtrx(:,2)/0.d0,1.d0,0.d0,0.d0,1.d0,0.d0/
   DATA mtrx(:,3)/0.d0,0.d0,1.d0,0.d0,0.d0,1.d0/
   DATA mtrx(:,6)/1.d0,1.d0,0.d0,-1.d0,-1.d0,0.d0/
   ! Run section
   masa=sum(system_mass(:))
   nua=system_mass(1)/masa
   nub=system_mass(2)/masa
   !
   molcoord(1)=(1.D0/(masa))*(atomcoord(1)*system_mass(1)+atomcoord(4)*system_mass(2))
   molcoord(2)=(1.D0/(masa))*(atomcoord(2)*system_mass(1)+atomcoord(5)*system_mass(2))
   molcoord(3)=(1.D0/(masa))*(atomcoord(3)*system_mass(1)+atomcoord(6)*system_mass(2))
   molcoord(4)=dsqrt((atomcoord(1)-atomcoord(4))**2.D0+&
      (atomcoord(2)-atomcoord(5))**2.D0+(atomcoord(3)-atomcoord(6))**2.D0)
   molcoord(5)=dacos((atomcoord(3)-atomcoord(6))/molcoord(4))
   molcoord(6)=datan2((atomcoord(2)-atomcoord(5)),(atomcoord(1)-atomcoord(4)))
   SELECT CASE(molcoord(6)<0.d0)
      CASE(.true.)
         molcoord(6)=molcoord(6)+2.d0*pi
      CASE(.false.)
         ! do nothing
   END SELECT
   mtrx(1,4)=dsin(molcoord(5))*dcos(molcoord(6))/nua
   mtrx(2,4)=dsin(molcoord(5))*dsin(molcoord(6))/nua
   mtrx(3,4)=dcos(molcoord(5))/nua
   mtrx(4,4)=-dsin(molcoord(5))*dcos(molcoord(6))/nub
   mtrx(5,4)=-dsin(molcoord(5))*dsin(molcoord(6))/nub
   mtrx(6,4)=-dcos(molcoord(5))/nub

   mtrx(1,5)=dcos(molcoord(5))*dcos(molcoord(6))/(nua*molcoord(4))
   mtrx(2,5)=dcos(molcoord(5))*dsin(molcoord(6))/(nua*molcoord(4))
   mtrx(3,5)=dsin(molcoord(5))/(nua*molcoord(4))
   mtrx(4,5)=-dcos(molcoord(5))*dcos(molcoord(6))/(nub*molcoord(4))
   mtrx(5,5)=-dcos(molcoord(5))*dsin(molcoord(6))/(nub*molcoord(4))
   mtrx(6,5)=-dsin(molcoord(5))/(nub*molcoord(4))
   SELECT CASE(dsin(molcoord(5)) /= 0.d0)
       CASE(.true.)
          mtrx(1,6)=-dsin(molcoord(6))/(nua*dsin(molcoord(5)))
          mtrx(2,6)=dcos(molcoord(6))/(nua*dsin(molcoord(5)))
          mtrx(3,6)=0.d0
          mtrx(4,6)=dsin(molcoord(6))/(nub*dsin(molcoord(5)))
          mtrx(5,6)=-dcos(molcoord(6))/(nub*dsin(molcoord(5)))
          mtrx(6,6)=0.d0
          CALL INV_MTRX(6,mtrx,invMtrx)
          molcoord(7:12)=matmul(invMtrx,atomcoord(7:12))
       CASE(.false.)
         CALL INV_MTRX(6,mtrx,invMtrx)
         molcoord(7:12)=matmul(invMtrx,atomcoord(7:12))
         molcoord(12)=0.d0
   END SELECT
   RETURN
END FUNCTION from_atomic_to_molecular_phaseSpace
!###############################################################################################
!# FUNCTION : correctSphPoint ##################################################################
!###############################################################################################
!> @brief
!! Calculates the correct definition of a given point in phase space restricted to system DOFs.
!
!> @param[in] phaseSpacePoint(12) - real(kind=8) Point in phase space. Assumed in spherical coordinates
!------------------------------------------------------------------------------------------------
FUNCTION correctSphPoint(phaseSpacePoint) result(goodPoint)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   REAL(KIND=8),INTENT(IN),DIMENSION(12):: phaseSpacePoint
   ! Dummy function variable
   REAL(KIND=8),DIMENSION(12):: goodPoint
   CHARACTER(LEN=*),PARAMETER:: routinename='correctSphPoint: '
   ! Run section
   goodPoint(:)=phaseSpacePoint(:)
   ! Theta switches
   SELECT CASE(dsin(goodPoint(5))<0.d0)
      CASE(.true.) ! project into 0~pi range
         goodPoint(5)=dacos(dcos(goodPoint(5)))
      CASE(.false.)
         ! do nothing
   END SELECT
   SELECT CASE(dsin(goodPoint(5))==0.d0)
      CASE(.true.) ! phi angle and momenta cannot be deffined. Set them to zero
         goodPoint(6)=0.d0
         goodPoint(12)=0.d0
      CASE(.false.)
         ! do nothing
   END SELECT
   ! Phi Switch
   SELECT CASE(goodPoint(6)<0)
      CASE(.true.)
         goodPoint(6)=goodPoint(6)+2.d0*pi
      CASE(.false.)
         ! do nothing
   END SELECT
   SELECT CASE(goodPoint(6)>2.0*pi)
      CASE(.true.)
         goodPoint(6)=goodPoint(6)-2.d0*pi
      CASE(.false.)
         ! do nothing
   END SELECT
END FUNCTION correctSphPoint
!###############################################################################################
!# FUNCTION : normalDistRandom #################################################################
!###############################################################################################
!> @brief
!! Review of old gasdev subroutine to generate a normally distributed deviate with zero mean and
!! unit variance. Old file can be checked in FORTRAN'77 Numerical Recipes. Result will lie
!! between -1.d0 and 1.d0
!
!> @details
!! - Uses random_number intrinsic subroutine to generate uniform deviates
!! - Assumes that random_number was initialized at some point in the code prior to the call
!------------------------------------------------------------------------------------------------
function normalDistRandom() result(rndReal)
   ! Initial declarations
   implicit none
   ! I/O variables
   real(kind=8):: rndReal ! dummy function variable
   ! Local variables
   real(kind=8):: r2
   real(kind=8),dimension(2):: v
   real(kind=8):: rndStored
   logical:: useStored
   ! Initialization/Storage variables section
   save rndStored,useStored  ! conserve this values from call to call
   data useStored/.false./             ! don't use useStored by default
   ! Run section --------------------------
   select case(useStored)
      case(.true.)
         rndReal=rndStored
         useStored=.false. ! next call won't use stored value

      case(.false.)
      r2=3.d0 ! ensures bad initial r2 number
      do while (r2>=1.0 .or. r2==0.d0)
         call random_number(v)
         v(:)=2.d0*v(:)-1.d0
         r2=dot_product(v,v)
      end do
      rndReal=v(1)*dsqrt(-2.d0*dlog(r2)/r2)
      rndStored=v(2)*dsqrt(-2.d0*dlog(r2)/r2)
      useStored=.true. ! next call will use stored value

   end select
   return
end function normalDistRandom
!############################################################################
! SUBROUTINE: skipHeaderFromFile ############################################
!############################################################################
!> @brief
!! Routine that skips all header lines that start with a hash '#' symbol.
!! Useful to read typical OOPtrajtou output files
!----------------------------------------------------------------------------
subroutine skipHeaderFromFile(fileName,unit,skippedRows,errCode)
   ! initial declarations
   implicit none
   ! I/O variables
   character(len=*),optional,intent(in):: fileName
   integer(kind=4),optional,intent(in):: unit
   integer(kind=4),optional,intent(out):: skippedRows
   integer(kind=4),optional,intent(out):: errCode
   ! Local variables
   integer(kind=4):: ru,counter,i
   integer(kind=4),dimension(2):: iErr
   logical:: ready
   character(len=1024):: auxString
   character(len=1):: hash
   ! Run section -------------------------
   select case( present(fileName) .and. .not.present(unit) )
   case(.true.)
      inquire(file=fileName,number=ru,opened=ready)
      auxString=fileName
   case(.false.)
      ! do nothing
   end select

   select case( present(unit) .and. .not.present(fileName) )
   case(.true.)
      inquire(unit=unit,opened=ready,name=auxString)
      ru=unit
   case(.false.)
      write(0,*) "ERR SKIPHEADERFROMFILE: set unit or filename. Don't set both or neither of them"
      call exit(1)
   end select

   select case(ready)
   case(.true.)
      ! go on
   case(.false.)
      write(0,*) 'SKIPHEADERFROMFILE: '//trim(auxString)//' was not opened'
      call exit(1)
   end select
   iErr(:)=0
   counter=0
   do while( iErr(1)==0 .and. iErr(2)==0 )
      read(ru,*,iostat=iErr(1)) hash
      select case(hash)
      case('#')
         counter=counter+1
      case default
        iErr(2)=1
      end select
   enddo
   select case( iErr(1) /= 0 )
   case(.true.) ! exit with error code stored at errCode
      if( present(errCode) ) errCode=iErr(1)
      return
   case(.false.)
      ! do nothing
   end select
   rewind(unit=ru)
   if ( present(skippedRows) ) skippedRows=counter
   do i=1,counter
      read(ru,*)
   enddo
   return ! Now we've skipped all commented lines
end subroutine skipHeaderFromFile
!###########################################################
! SUBROUTINE: DELETEEMPTYFILE
!###########################################################
!> @brief
!! Deletes a file if it is empty or if it only has a header: lines
!! started by a hash symbol '#'
!-----------------------------------------------------------
subroutine deleteIfEmptyFile(fileName,errCode)
   ! initial declaration
   implicit none
   ! I/O variables
   character(len=*),intent(in):: fileName
   integer(kind=4),optional,intent(out):: errCode
   ! Local variables
   integer(kind=4):: iErr
   logical:: isOpened,exists
   ! Parameters
   integer(kind=4),parameter:: ru=623

   ! Run section -------------------------
   inquire(file=fileName,opened=isOpened,exist=exists)
   select case( .not.isOpened .and. exists )
   case(.true.)
      ! go on
   case(.false.)
      write(0,*) 'ERR DELETEEMPTYFILE '//fileName//' is still opened or does not exist'
      call exit(1)
   end select
   open(unit=ru,file=fileName,status='old',action='read')
   call skipHeaderFromFile(unit=ru,errCode=iErr)
   select case( iErr )
   case(-1)
      close(unit=ru,status='delete')
   case default
      close(unit=ru)
   end select
   if( present(errCode) ) errCode=iErr
   return
end subroutine deleteIfEmptyfile
!##########################################################
! FUNCTION: evaluateEnergyRovibrState
!##########################################################
!> @brief
!! - Gives vibration energy and average rotational energy for
!!   a given rovibrational state.
!
!> @details
!! @b Implemented @b methods:
!!    - Morse: vibrational energy is calculated with Morse corrections.
!> @warnings
!! - Rotational energy is always considered as : J(J+1)/2muReq
function evaluateEnergyRovibrState(rovibrState,eVibr,eRot) result(energy)
   ! initial declarations
   implicit none
   ! I/O variables
   integer(kind=4),dimension(3),intent(in):: rovibrState
   real(kind=8),optional,intent(out):: eVibr,eRot
   ! function dummy variable
   real(kind=8):: energy
   real(kind=8):: vibr,rot,ed,width,req,omega
   real(kind=8):: v,J,mu
   ! Run section
   if( system_binningScheme == 'Morse' ) then
      mu=product(system_mass(:))/sum(system_mass(:))
      ed=system_binningParam(1)
      req=system_binningParam(2)
      width=system_binningParam(3)
      omega=width*dsqrt(2.d0*ed/mu)
      v=dfloat(rovibrState(1))
      J=dfloat(rovibrState(2))
      vibr=omega*(v+0.5d0)-((omega**2.d0)/(4.d0*ed))*(v+0.5d0)**2.d0
      rot=j*(j+1)/(2.d0*mu*req**2.d0)
   elseif( system_binningScheme == 'Dong' )
      if( rovibrState(1) == 0 ) then
         vibr=0.01d0
         if( rovibrState(2) == 0) then
            rot=0.d0
         elseif( rovibrState(2) == 1) then
            rot=0.0005d0
         elseif( rovibrState(2) == 2) then
            rot=0.0016d0
         elseif( rovibrState(2) == 3) then
            rot=0.0031d0
         elseif( rovibrState(2) == 4) then
            rot=0.0052d0
         elseif( rovibrState(2) == 5) then
            rot=0.0077d0
         else
            write(0,*) 'evaluateEnergyRovibrState ERR: rotational number not implemented'
            write(0,*) 'implemented ones: up to J=5'
            call exit(1)
         endif
      elseif( rovibrState(1) == 1 ) then
         vibr=0.0288d0
         if( rovibrState(2) == 0) then
            rot=0.d0
         elseif( rovibrState(2) == 1) then
            rot=0.0005d0
         elseif( rovibrState(2) == 2) then
            rot=0.0015d0
         elseif( rovibrState(2) == 3) then
            rot=0.0028d0
         elseif( rovibrState(2) == 4) then
            rot=0.0049d0
         elseif( rovibrState(2) == 5) then
            rot=0.0072d0
         else
            write(0,*) 'evaluateEnergyRovibrState ERR: rotational number not implemented'
            write(0,*) 'implemented ones: up to J=5'
            call exit(1)
         endif
      elseif( rovibrState(1) == 2 ) then
         vibr=0.0464d0
         if( rovibrState(2) == 0) then
            rot=0.d0
         elseif( rovibrState(2) == 1) then
            rot=0.0004d0
         elseif( rovibrState(2) == 2) then
            rot=0.0013d0
         elseif( rovibrState(2) == 3) then
            rot=0.0027d0
         elseif( rovibrState(2) == 4) then
            rot=0.0045d0
         elseif( rovibrState(2) == 5) then
            rot=0.0067d0
         else
            write(0,*) 'evaluateEnergyRovibrState ERR: rotational number not implemented'
            write(0,*) 'implemented ones: up to J=5'
            call exit(1)
         endif
      elseif( rovibrState(1) == 3 ) then
         vibr=0.0628d0
         if( rovibrState(2) == 0) then
            rot=0.d0
         elseif( rovibrState(2) == 1) then
            rot=0.0004d0
         elseif( rovibrState(2) == 2) then
            rot=0.0013d0
         elseif( rovibrState(2) == 3) then
            rot=0.0026d0
         elseif( rovibrState(2) == 4) then
            rot=0.0043d0
         elseif( rovibrState(2) == 5) then
            rot=0.0064d0
         else
            write(0,*) 'evaluateEnergyRovibrState ERR: rotational number not implemented'
            write(0,*) 'implemented ones: up to J=5'
            call exit(1)
         endif
      elseif( rovibrState(1) == 4 ) then
         vibr=0.0783d0
         if( rovibrState(2) == 0) then
            rot=0.d0
         elseif( rovibrState(2) == 1) then
            rot=0.0004d0
         elseif( rovibrState(2) == 2) then
            rot=0.0012d0
         elseif( rovibrState(2) == 3) then
            rot=0.0024d0
         elseif( rovibrState(2) == 4) then
            rot=0.0040d0
         elseif( rovibrState(2) == 5) then
            rot=0.0060d0
         else
            write(0,*) 'evaluateEnergyRovibrState ERR: rotational number not implemented'
            write(0,*) 'implemented ones: up to J=5'
            call exit(1)
         endif
      elseif( rovibrState(1) == 5 ) then
         vibr=0.0929d0
         if( rovibrState(2) == 0) then
            rot=0.d0
         elseif( rovibrState(2) == 1) then
            rot=0.0004d0
         elseif( rovibrState(2) == 2) then
            rot=0.0012d0
         elseif( rovibrState(2) == 3) then
            rot=0.0023d0
         elseif( rovibrState(2) == 4) then
            rot=0.0038d0
         elseif( rovibrState(2) == 5) then
            rot=0.0057d0
         else
            write(0,*) 'evaluateEnergyRovibrState ERR: rotational number not implemented'
            write(0,*) 'implemented ones: up to J=5'
            call exit(1)
         endif
      else
         write(0,*) 'evaluateEnergyRovibrState ERR: vibrational number not implemented'
         write(0,*) 'implemented ones: up to V=5'
         call exit(1)
      endif
         
   else
      write(0,*) 'evaluateEnergyRovibrState ERR: wrong binning type: '//system_binningScheme
      write(0,*) 'implemented ones: Morse'
      write(0,*) 'case sensitive'
      call exit(1)
   endif
   energy=vibr+rot
   if( present(eVibr) ) eVibr=vibr
   if( present(eRot) )  eRot=rot
   return
end function
END MODULE SYSTEM_MOD
