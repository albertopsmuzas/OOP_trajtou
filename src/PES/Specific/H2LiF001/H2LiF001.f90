!#######################################################
! MODULE: PES_H2LiF001_MOD
!#######################################################
!> @brief
!! CRP6D PES explicit implementation
!#######################################################
module PES_H2LiF001_MOD
! Initial declarations
use PES_MOD, only: PES
use LiF001SURF_MOD, only: LiF001Surf,pi
use LOGISTIC_FUNCTION_MOD, only: Logistic_func
use XEXPONENTIAL_FUNCTION_MOD, only: Xexponential_func
use PES_HLIF001_NS_MOD, only: Pes_HLiF001_NS
use EXTRAPOL_TO_VACUUM_MOD, only: Vacuumpot
use FOURIER_P4MM_MOD, only: Fourierp4mm
use WYCKOFF_P4MM_MOD, only: WyckoffSitio, Wyckoffp4mm
implicit none
! Local module variable, used to simulate SYSTEM_MOD
type(LiF001Surf),private:: sysLiF001Surf
!/////////////////////////////////////////////////
! TYPE: PES_H2LiF001
!-------------------------------------------------
type,extends(PES):: PES_H2LiF001
   integer(kind=4):: nsites=4
   integer(kind=4):: natomic=1
   logical:: is_interpolated=.false.
   logical:: is_homonucl=.false.
   logical:: is_smooth=.false.
   logical:: is_shifted=.false.
   logical:: is_resized=.false.
   real(kind=8):: zvacuum
   integer(kind=4),dimension(2):: grid=[25,50]
   type(Wyckoffp4mm),dimension(:),allocatable:: wyckoffSite
   type(Pes_HLiF001_NS),dimension(:),allocatable:: atomicCrp
   type(VacuumPot):: farpot
   type(Logistic_func):: dumpfunc
   character(len=30):: extrapol2vac_flag='Xexponential'
   integer(kind=4),dimension(:,:),allocatable:: xyklist
   contains
      ! Initialization block
      procedure,public:: initialize                  => INITIALIZE_PES_H2LiF001
      ! Set block
      procedure,public:: set_smooth                  => SET_SMOOTH_PES_H2LiF001
      ! Get block
      procedure,public:: get_v_and_derivs            => GET_V_AND_DERIVS_PES_H2LiF001
      procedure,public:: get_v_and_derivs_pure       => GET_V_AND_DERIVS_PURE_PES_H2LiF001
      procedure,public:: get_atomicpot_and_derivs    => GET_ATOMICPOT_AND_DERIVS_PES_H2LiF001
      ! Tools block
      procedure,public:: smooth                      => SMOOTH_PES_H2LiF001
      procedure,public:: interpol                    => INTERPOL_PES_H2LiF001
      procedure,public:: extract_vacuumsurf          => EXTRACT_VACUUMSURF_PES_H2LiF001
      procedure,public:: add_vacuumsurf              => ADD_VACUUMSURF_PES_H2LiF001
      procedure,public:: interpol_new_rzgrid         => INTERPOL_NEW_RZGRID_PES_H2LiF001
      ! Plot toolk
      procedure,public:: plot1d_theta                => PLOT1D_THETA_PES_H2LiF001
      procedure,public:: plot1d_atomic_interac_theta => PLOT1D_ATOMIC_INTERAC_THETA_PES_H2LiF001
      procedure,public:: plot1d_phi                  => PLOT1D_PHI_PES_H2LiF001
      procedure,public:: plot1d_atomic_interac_phi   => PLOT1D_ATOMIC_INTERAC_PHI_PES_H2LiF001
      procedure,public:: plot1d_r                    => PLOT1D_R_PES_H2LiF001
      procedure,public:: plot1d_z                    => PLOT1D_Z_PES_H2LiF001
      procedure,public:: plot_xymap                  => PLOT_XYMAP_PES_H2LiF001
      procedure,public:: plot_rzmap                  => PLOT_RZMAP_PES_H2LiF001
      procedure,public:: plot_atomic_interac_rz      => PLOT_ATOMIC_INTERAC_RZ_PES_H2LiF001
      ! Enquire block
      procedure,public:: is_allowed                  => is_allowed_PES_H2LiF001
end type PES_H2LiF001

private initialize_PES_H2LiF001,set_smooth_PES_H2LiF001,get_v_and_derivs_PES_H2LiF001,get_v_and_derivs_pure_PES_H2LiF001,&
        get_atomicPot_and_derivs_PES_H2LiF001,smooth_PES_H2LiF001,interpol_PES_H2LiF001,extract_vacuumSurf_PES_H2LiF001,&
        add_vacuumSurf_PES_H2LiF001,interpol_new_RZgrid_PES_H2LiF001,plot1d_theta_PES_H2LiF001,&
        plot1d_atomic_interac_theta_PES_H2LiF001,plot1d_phi_PES_H2LiF001,plot1d_atomic_interac_phi_PES_H2LiF001,&
        plot1d_R_PES_H2LiF001,plot1d_Z_PES_H2LiF001,plot_XYmap_PES_H2LiF001,plot_RZmap_PES_H2LiF001,&
        plot_atomic_interac_RZ_PES_H2LiF001,is_allowed_PES_H2LiF001,from_molecular_to_atomic

contains
!###########################################################
!# SUBROUTINE: GET_ATOMICPOT_AND_DERIVS_PES_H2LiF001
!###########################################################
SUBROUTINE GET_ATOMICPOT_AND_DERIVS_PES_H2LiF001(this,molecx,atomicx,v,dvdu)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),INTENT(IN):: this
   REAL(KIND=8),DIMENSION(6),INTENT(IN):: molecx
   REAL(KIND=8),DIMENSION(2),INTENT(OUT):: v
   REAL(KIND=8),DIMENSION(6),INTENT(OUT):: dvdu
   REAL(KIND=8),DIMENSION(6),INTENT(OUT):: atomicx
   ! Local variables
   REAL(KIND=8):: vcorra,vcorrb
   ! Run section
   atomicx(:)=from_molecular_to_atomic(molecx)
   SELECT CASE(this%natomic)
      CASE(1)
         CALL this%atomiccrp(1)%GET_V_AND_DERIVS(atomicx(1:3),v(1),dvdu(1:3))
         CALL this%atomiccrp(1)%GET_V_AND_DERIVS(atomicx(4:6),v(2),dvdu(4:6))
      CASE(2)
         CALL this%atomiccrp(1)%GET_V_AND_DERIVS(atomicx(1:3),v(1),dvdu(1:3))
         CALL this%atomiccrp(2)%GET_V_AND_DERIVS(atomicx(4:6),v(2),dvdu(4:6))
      CASE DEFAULT
         WRITE(0,*) "GET_ATOMICPOT_AND_DERIVS_PES_H2LiF001 ERR: wrong number of atomic potentials"
         CALL EXIT(1)
   END SELECT
   vcorra=this%dumpfunc%getvalue(atomicx(3))
   vcorrb=this%dumpfunc%getvalue(atomicx(6))
   dvdu(1)=dvdu(1)*vcorra
   dvdu(2)=dvdu(2)*vcorra
   dvdu(3)=dvdu(3)*vcorra+v(1)*this%dumpfunc%getderiv(atomicx(3))
   dvdu(4)=dvdu(4)*vcorrb
   dvdu(5)=dvdu(5)*vcorrb
   dvdu(6)=dvdu(6)*vcorrb+v(2)*this%dumpfunc%getderiv(atomicx(6))
   v(1)=v(1)*vcorra
   v(2)=v(2)*vcorrb
   RETURN
END SUBROUTINE GET_ATOMICPOT_AND_DERIVS_PES_H2LiF001
!###########################################################
!# SUBROUTINE: SET_SMOOTH_PES_H2LiF001
!###########################################################
SUBROUTINE SET_SMOOTH_PES_H2LiF001(this,is_smooth)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),INTENT(INOUT):: this
   LOGICAL,INTENT(IN) :: is_smooth
   ! Run section
   this%is_smooth=is_smooth
   RETURN
END SUBROUTINE SET_SMOOTH_PES_H2LiF001
!###########################################################
!# SUBROUTINE: INTERPOL_NEW_RZGRID_PES_H2LiF001
!###########################################################
SUBROUTINE INTERPOL_NEW_RZGRID_PES_H2LiF001(this,nRpoints,nZpoints)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),INTENT(INOUT) :: this
   INTEGER(KIND=4),INTENT(IN) :: nrpoints,nzpoints
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   ! Run section
   DO i = 1, this%nsites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         CALL this%wyckoffsite(i)%zrcut(j)%interrz%INTERPOL_NEWGRID(nrpoints,nzpoints)
      END DO
   END DO
   RETURN
END SUBROUTINE INTERPOL_NEW_RZGRID_PES_H2LiF001
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_PURE_PES_H2LiF001
!###########################################################
SUBROUTINE GET_V_AND_DERIVS_PURE_PES_H2LiF001(this,x,v,dvdu)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),TARGET,INTENT(IN):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN):: x
   REAL(KIND=8),INTENT(OUT):: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT):: dvdu
   ! Local variables
   INTEGER(KIND=4):: i ! counters
   REAL(KIND=8):: ma,mb
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: f ! smooth function and derivs
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: xy
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: aux1
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: aux2
   REAL(KIND=8),DIMENSION(6):: atomicx
   REAL(KIND=8),DIMENSION(2):: atomic_v
   REAL(KIND=8),DIMENSION(6):: atomic_dvdu
   REAL(KIND=8),DIMENSION(3):: dvdu_atomicA,dvdu_atomicB
   TYPE(Fourierp4mm):: xyinterpol
   CHARACTER(LEN=*),PARAMETER:: routinename="GET_V_AND_DERIVS_PES_H2LiF001: "
   ! Run section
   SELECT CASE(this%is_smooth)
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         WRITE(0,*) "GET_V_AND_DERIVS_PES_H2LiF001 ERR: smooth the PES first (CALL thispes%SMOOTH())"
         CALl EXIT(1)
   END SELECT
   ALLOCATE(f(5,this%nsites))
   ALLOCATE(xy(this%nsites,2))
   DO i = 1, this%nsites 
      CALL this%wyckoffsite(i)%GET_V_AND_DERIVS(x(3:6),f(1,i),f(2:5,i))
      xy(i,1)=this%wyckoffsite(i)%x
      xy(i,2)=this%wyckoffsite(i)%y
   END DO
   ! f(1,:) smooth potential values
   ! f(2,:) smooth dvdz
   ! f(3,:) smooth dvdr
   ! f(4,:) smooth dvdtheta
   ! f(5,:) smooth dvdphi
   CALL xyinterpol%READ(xy,f,this%xyklist)
   CALL xyinterpol%INTERPOL(sysLiF001Surf)
   ALLOCATE(aux1(5))
   ALLOCATE(aux2(5,2))
   CALL xyinterpol%GET_F_AND_DERIVS(sysLiF001Surf,x(1:2),aux1,aux2)
   !--------------------------------------
   ! Results for the real potential
   !-------------------------------------
   CALL this%GET_ATOMICPOT_AND_DERIVS(x,atomicx,atomic_v,atomic_dvdu)
   dvdu_atomicA=atomic_dvdu(1:3)
   dvdu_atomicB=atomic_dvdu(4:6)
   v=aux1(1)+sum(atomic_v)
   ! It does not matter ma and mb values inside this subroutine as long ma=mb for D2 and H2. I let
   ! equations to depend on explicit values of Ma and Mb just in case someone want to
   ! implement an explicit HD PES in the future
   ma=1.d0
   mb=1.d0
   dvdu(1)=aux2(1,1)+dvdu_atomicA(1)+dvdu_atomicB(1)
   dvdu(2)=aux2(1,2)+dvdu_atomicA(2)+dvdu_atomicB(2)
   dvdu(3)=aux1(2)+dvdu_atomicA(3)+dvdu_atomicB(3)
   dvdu(4)=aux1(3)&
      +(mb/(ma+mb))*dvdu_atomicA(1)*dcos(x(6))*dsin(x(5))&
      +(mb/(ma+mb))*dvdu_atomicA(2)*dsin(x(6))*dsin(x(5))&
      +(mb/(ma+mb))*dvdu_atomicA(3)*dcos(x(5))&
      -(ma/(ma+mb))*dvdu_atomicB(1)*dcos(x(6))*dsin(x(5))&
      -(ma/(ma+mb))*dvdu_atomicB(2)*dsin(x(6))*dsin(x(5))&
      -(ma/(ma+mb))*dvdu_atomicB(3)*dcos(x(5))
   dvdu(5)=aux1(4)&
      +(mb/(ma+mb))*x(4)*dvdu_atomicA(1)*dcos(x(6))*dcos(x(5))&
      +(mb/(ma+mb))*x(4)*dvdu_atomicA(2)*dsin(x(6))*dcos(x(5))&
      -(mb/(ma+mb))*x(4)*dvdu_atomicA(3)*dsin(x(5))&
      -(ma/(ma+mb))*x(4)*dvdu_atomicB(1)*dcos(x(6))*dcos(x(5))&
      -(ma/(ma+mb))*x(4)*dvdu_atomicB(2)*dsin(x(6))*dcos(x(5))&
      +(ma/(ma+mb))*x(4)*dvdu_atomicB(3)*dsin(x(5))
   dvdu(6)=aux1(5)&
      -(mb/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicA(1)*dsin(x(6))&
      +(mb/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicA(2)*dcos(x(6))&
      +(ma/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicB(1)*dsin(x(6))&
      -(ma/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicB(2)*dcos(x(6))
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_PURE_PES_H2LiF001
!###########################################################
!# SUBROUTINE: GET_V AND DERIVS_PES_H2LiF001
!###########################################################
subroutine get_v_and_derivs_PES_H2LiF001(this,X,v,dvdu,errCode)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   class(PES_H2LiF001),target,intent(in):: this
   real(kind=8),dimension(:),intent(in) :: x
   real(kind=8),intent(out) :: v
   real(kind=8),dimension(:),intent(out) :: dvdu
   integer(kind=1),optional,intent(out):: errCode
   ! Local variables
   real(kind=8) :: zcrp, zvac ! last PES_H2LiF001 z value and Z infinity
   real(kind=8) :: vzcrp, vzvac ! potentials at zcrp and zvac
   real(kind=8),dimension(6) :: dvducrp ! derivatives at zcrp
   real(kind=8),dimension(6) :: dvduvac ! derivatives at vacuum
   real(kind=8) :: alpha,beta,gama ! parameters
   type(Xexponential_func):: extrapolfunc
   integer(kind=4) :: i !counter
   ! local parameter
   real(kind=8),parameter :: zero=0.D0 ! what we will consider zero (a.u.)
   real(kind=8),parameter :: dz=0.5D0 ! 0.25 Angstroems approx
   character(len=*),parameter:: routinename='GET_V_AND_DERIVS_PES_H2LiF001: '
   ! Run section
   zcrp=this%wyckoffsite(1)%zrcut(1)%getlastZ()
   zvac=this%zvacuum
   ! Check if we are in the pure PES_H2LiF001 region
   select case(x(3)<= zcrp) !easy
   case(.true.)
      call this%get_v_and_derivs_pure(x,v,dvdu)
      return
   case(.false.)
      ! do nothing, next switch
   end select
   ! check if we are in the extrapolation region
   select case(x(3)>zcrp .AND. x(3)<zvac)
   case(.true.) ! uff
      ! Set potential and derivs
      vzvac=this%farpot%getpot(x(4))
      CALL this%GET_V_AND_DERIVS_PURE([x(1),x(2),zcrp,x(4),x(5),x(6)],vzcrp,dvducrp)
      dvduvac(1:3)=zero
      dvduvac(4)=this%farpot%getderiv(x(4))
      dvduvac(5:6)=zero
      ! Extrapol potential
      beta=-1.D0/zvac
      alpha=(vzcrp-vzvac)/(zcrp*dexp(beta*zcrp)-zvac*dexp(beta*zvac))
      gama=vzvac-alpha*zvac*dexp(beta*zvac)
      CALL extrapolfunc%READ([alpha,beta])
      v=extrapolfunc%getvalue(x(3))+gama
      dvdu(3)=extrapolfunc%getderiv(x(3))
      ! extrapol derivatives
      do i = 1, 6
         select case(i)
         case(3)
            ! skip dvdz
         case default
            beta=-1.D0/zvac
            alpha=(dvducrp(i)-dvduvac(i))/(zcrp*dexp(beta*zcrp)-zvac*dexp(beta*zvac))
            gama=dvduvac(i)-alpha*zvac*dexp(beta*zvac)
            call extrapolfunc%READ([alpha,beta])
            dvdu(i)=extrapolfunc%getvalue(x(3))+gama
         end select
      end do
      return

   case(.false.)
      ! do nothing
   end select
   ! check if we are in the Vacuum region
   select case(x(3)>=zvac) !easy
      case(.true.)
         v=this%farpot%getpot(x(4))
         dvdu(1:3)=0.D0
         dvdu(4)=this%farpot%getderiv(x(4))
         dvdu(5:6)=0.D0
      case(.false.) ! this's the last switch!
          errCode=1_1
   end select
   return
end subroutine get_v_and_derivs_PES_H2LiF001
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_SMOOTH_PES_H2LiF001
!###########################################################
SUBROUTINE GET_V_AND_DERIVS_SMOOTH_PES_H2LiF001(this,x,v,dvdu)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),INTENT(IN):: this
   REAL(KIND=8),DIMENSION(6),INTENT(IN):: x
   REAL(KIND=8),INTENT(OUT):: v
   REAL(KIND=8),DIMENSION(6),INTENT(OUT):: dvdu
   ! Local variables
   INTEGER(KIND=4):: i ! counters
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: f ! smooth function and derivs
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: xy
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: aux1
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: aux2
   TYPE(Fourierp4mm):: xyinterpol
   CHARACTER(LEN=*),PARAMETER:: routinename="GET_V_AND_DERIVS_PES_H2LiF001: "
   ! Run section
   SELECT CASE(x(3) >= this%wyckoffsite(1)%zrcut(1)%getlastz()) 
      CASE(.TRUE.)
         v=this%farpot%getpot(x(4))
         dvdu(1)=0.D0
         dvdu(2)=0.D0
         dvdu(3)=0.D0
         dvdu(4)=this%farpot%getderiv(x(4))
         dvdu(5)=0.D0
         dvdu(6)=0.D0
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   ALLOCATE(f(5,this%nsites))
   ALLOCATE(xy(this%nsites,2))
   DO i = 1, this%nsites 
      CALL this%wyckoffsite(i)%GET_V_AND_DERIVS(x(3:6),f(1,i),f(2:5,i))
      xy(i,1)=this%wyckoffsite(i)%x
      xy(i,2)=this%wyckoffsite(i)%y
   END DO
   ! f(1,:) smooth potential values
   ! f(2,:) smooth dvdz
   ! f(3,:) smooth dvdr
   ! f(4,:) smooth dvdtheta
   ! f(5,:) smooth dvdphi
   CALL xyinterpol%READ(xy,f,this%xyklist)
   CALL xyinterpol%INTERPOL(sysLiF001Surf)
   ALLOCATE(aux1(5))
   ALLOCATE(aux2(5,2))
   CALL xyinterpol%GET_F_AND_DERIVS(sysLiF001Surf,x(1:2),aux1,aux2)
   v=aux1(1)
   dvdu(1)=aux2(1,1)
   dvdu(2)=aux2(1,2)
   dvdu(3)=aux1(2)
   dvdu(4)=aux1(3)
   dvdu(5)=aux1(4)
   dvdu(6)=aux1(5)
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_SMOOTH_PES_H2LiF001
!###########################################################
!# SUBROUTINE: INTERPOL_PES_H2LiF001
!###########################################################
SUBROUTINE INTERPOL_PES_H2LiF001(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),INTENT(INOUT)::this
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   ! Run section
   CALL this%EXTRACT_VACUUMSURF()   
   CALL this%SMOOTH()
   !CALL this%SMOOTH_EXTRA()
   DO i = 1, this%nsites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         CALL this%wyckoffsite(i)%zrcut(j)%INTERPOL()
      END DO
   END DO
   this%is_interpolated=.TRUE.
   RETURN
END SUBROUTINE INTERPOL_PES_H2LiF001
!###########################################################
!# SUBROUTINE: SMOOTH_PES_H2LiF001
!###########################################################
SUBROUTINE SMOOTH_PES_H2LiF001(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
    CLASS(PES_H2LiF001),INTENT(INOUT):: this
   ! Local variables
   REAL(KIND=8),DIMENSION(6):: molcoord,atomcoord,dummy
   INTEGER(KIND=4):: nr,nz
   INTEGER(KIND=4):: i,j,k,l ! counters
   CHARACTER(LEN=*),PARAMETER:: routinename="SMOOTH_PES_H2LiF001: "
   REAL(KIND=8):: newpot
   REAL(KIND=8),DIMENSION(2):: atomic_v
   ! Run section
   DO i = 1, this%nsites ! cycle wyckoff sites
      molcoord(1)=this%wyckoffsite(i)%x
      molcoord(2)=this%wyckoffsite(i)%y
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         molcoord(5)=this%wyckoffsite(i)%zrcut(j)%theta
         molcoord(6)=this%wyckoffsite(i)%zrcut(j)%phi
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizeR()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizeZ()
         DO k = 1, nr
            DO l = 1, nz
               molcoord(3)=this%wyckoffsite(i)%zrcut(j)%getgridvalueZ(l)
               molcoord(4)=this%wyckoffsite(i)%zrcut(j)%getgridvalueR(k)
               CALL this%GET_ATOMICPOT_AND_DERIVS(molcoord,atomcoord,atomic_v,dummy)
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)-sum(atomic_v)
               CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
            END DO
         END DO
      END DO
   END DO
   this%is_smooth=.TRUE.
   RETURN
END SUBROUTINE SMOOTH_PES_H2LiF001
!###########################################################
!# SUBROUTINE: EXTRACT_VACUUMSURF_PES_H2LiF001
!###########################################################
SUBROUTINE EXTRACT_VACUUMSURF_PES_H2LiF001(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
    CLASS(PES_H2LiF001),INTENT(INOUT) :: this
   ! Local variables
   INTEGER(KIND=4) :: nr,nz
   INTEGER(KIND=4) :: i,j,k,l ! counters
   REAL(KIND=8) :: newpot
   CHARACTER(LEN=*),PARAMETER :: routinename="EXTRACT_VACUUMSURF_PES_H2LiF001: "
   ! Run section
   DO i = 1, this%nsites ! cycle wyckoff sites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizeR()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizeZ()
         DO k = 1, nr
            DO l = 1, nz
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)
               newpot=newpot-this%farpot%getscalefactor()
               CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
            END DO
         END DO
      END DO
   END DO
   CALL this%farpot%SHIFTPOT()
   this%is_shifted=.TRUE.
   RETURN
END SUBROUTINE EXTRACT_VACUUMSURF_PES_H2LiF001
!###########################################################
!# SUBROUTINE: ADD_VACUUMSURF_PES_H2LiF001
!###########################################################
SUBROUTINE ADD_VACUUMSURF_PES_H2LiF001(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
    CLASS(PES_H2LiF001),INTENT(INOUT) :: this
   ! Local variables
   INTEGER(KIND=4) :: nr,nz
   INTEGER(KIND=4) :: i,j,k,l ! counters
   REAL(KIND=8) :: newpot
   CHARACTER(LEN=*),PARAMETER :: routinename="ADD_VACUUMSURF_PES_H2LiF001: "
   ! Run section
   DO i = 1, this%nsites ! cycle wyckoff sites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizeR()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizeZ()
         DO k = 1, nr
            DO l = 1, nz
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)
               newpot=newpot+this%farpot%getscalefactor()
               CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
            END DO
         END DO
      END DO
   END DO
   CALL this%farpot%SHIFTPOT_UNDO()
   this%is_shifted=.FALSE.
   RETURN
END SUBROUTINE ADD_VACUUMSURF_PES_H2LiF001
!###########################################################
!# SUBROUTINE: READ_PES_H2LiF001
!###########################################################
subroutine initialize_PES_H2LiF001(this,filename,tablename)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(Pes_H2LiF001),intent(out):: this
   character(len=*),optional,intent(in):: filename,tablename
   ! local variables
   real(kind=8),dimension(14,14):: gridPot1414
   real(kind=8),dimension(14,16):: gridPot1416
   ! Parameters
   character(len=*),parameter:: routinename="READ_PES_H2LiF001: "
   real(kind=8),dimension(14),parameter:: gridR14=[ 0.7558904532d0,0.9448630664d0,1.1338356797d0,1.2283219864d0,&
                                                    1.3228082930d0,1.4172945997d0,1.5117809063d0,1.6062672130d0,&
                                                    1.8897261329d0,2.3621576661d0,2.8345891993d0,3.3070207325d0,&
                                                    3.7794522658d0,4.3463701056d0 ]
   real(kind=8),dimension(14),parameter:: gridZ14=[ 0.4724315332d0,0.9448630664d0,1.8897261329d0,2.8345891993d0,&
                                                    3.4015070392d0,3.9684248791d0,4.7243153322d0,5.4802057854d0,&
                                                    6.0471236252d0,6.4250688518d0,6.8030140784d0,7.1809593050d0,&
                                                    7.3699319183d0,7.5589045315d0 ]
   real(kind=8),dimension(16),parameter:: gridZ16=[ 0.4724315332d0,0.9448630664d0,1.3983973383d0,1.8897261329d0,&
                                                    2.3054658821d0,2.8345891993d0,3.4015070392d0,3.9684248791d0,&
                                                    4.7243153322d0,5.4802057854d0,6.0471236252d0,6.4250688518d0,&
                                                    6.8030140784d0,7.1809593050d0,7.3699319183d0,7.5589045315d0 ]
   ! Run section -----------------------
   call this%set_pesType('PES_H2LiF001')
   call this%set_alias('H2_on_LiF001')
   call this%set_dimensions(6)
   this%is_homonucl=.true.
   call this%farPot%initialize_direct( x=[ 0.7558904532d0,0.9448630664d0,1.1338356797d0,1.2283219864d0,1.3228082930d0,&
                                           1.4172945997d0,1.5117809063d0,1.6062672130d0,1.7007535196d0,1.8897261329d0,&
                                           2.3621576661d0,2.8345891993d0,3.3070207325d0,3.7794522658d0,4.3463701056d0 ],&
                                       f=[ -0.0412144309d0,-0.1744799612d0,-0.2298169263d0,-0.2422003318d0,-0.2483344814d0,&
                                           -0.2500191955d0,-0.2485072327d0,-0.2446916237d0,-0.2392079371d0,-0.2250493114d0,&
                                           -0.1828462465d0,-0.1423982815d0,-0.1082070069d0,-0.0810442795d0,-0.0564546003d0 ],&
                                       surfaceEnergy=-7.09306998104d0 )
   call this%dumpFunc%read([3.d0,4.d0])
   this%zvacuum=13.2280829302d0
   allocate( this%wyckoffSite(4) )
   allocate( this%atomicCrp(1) )
   call this%atomicCrp(1)%initialize()
   allocate( this%xyKlist(4,2) )
   this%xyklist(:,1)=[0,1,1,2]
   this%xyklist(:,2)=[0,0,1,0]
   call sysLiF001Surf%initialize('dummyString')
   ! Create wyckoff sites//////////////////////////////////
   ! Wickoff Top Li -----------> 'a'
   this%wyckoffSite(1)%id='a'
   this%wyckoffSite(1)%myNumber=1
   this%wyckoffSite(1)%is_homonucl=.true.
   this%wyckoffSite(1)%x=0.d0
   this%wyckoffSite(1)%y=0.d0
   this%wyckoffSite(1)%n2dcuts=5
   allocate( this%wyckoffSite(1)%nPhiPoints(3) )
   this%wyckoffSite(1)%nPhiPoints(:)=[1,2,2]
   allocate( this%wyckoffSite(1)%zrCut(5) )
   this%wyckoffSite(1)%zrCut(:)%x=0.d0
   this%wyckoffSite(1)%zrCut(:)%y=0.d0
   ! Reading zrcuts
   this%wyckoffSite(1)%zrcut(1)%alias='Top_Li_0_0'
   this%wyckoffSite(1)%zrCut(1)%theta=0.d0
   this%wyckoffSite(1)%zrCut(1)%phi=0.d0
   gridPot1416(:,:)=reshape( [ -2.3486329404d0, 1.0844017627d0, 8.8266781296d0, 9.1272988309d0, 6.5171540893d0, 1.6638643184d0,&
                               -0.9649307703d0,-2.4497027752d0,-4.5344580512d0,-5.7406206931d0,-6.3350386134d0,-6.6717466581d0,&
                               -6.8772623685d0,-7.0181287846d0,-5.8769513827d0,-5.7776624228d0,-5.5109630136d0,-5.2525253008d0,&
                               -5.0118020103d0,-4.6926259863d0,-4.2606813963d0,-3.6820328918d0, 0.6399452089d0, 1.3764483020d0,&
                               -4.7204565840d0,-5.8608453541d0,-6.4126558964d0,-6.7673050127d0,-6.6586128413d0,-6.7073472242d0,&
                               -6.6589171257d0,-6.6112389184d0,-6.5509718623d0,-6.4789630291d0,-6.3952271187d0,-6.2990585441d0,&
                               -5.9175156555d0,-4.6496075711d0, 2.0317609625d0, 0.0824348596d0,-4.9097362666d0,-6.0776842388d0,&
                               -6.9680495110d0,-7.0663885008d0,-7.0788310874d0,-7.0666236965d0,-7.0459661657d0,-7.0185276494d0,&
                               -6.9853878427d0,-6.9472269758d0,-6.8055553864d0,-6.4702699762d0,-5.9208833637d0,-4.8334825967d0,&
                                0.5489994402d0,-1.2821154148d0,-7.0648608314d0,-7.1814932679d0,-7.2164459187d0,-7.2171338661d0,&
                               -7.2105138426d0,-7.1983207840d0,-7.1817347109d0,-7.1615556564d0,-7.0839759930d0,-6.9074087142d0,&
                               -6.6625938957d0,-6.3011745703d0,-5.6766355873d0,-3.6644919811d0,-7.1108893614d0,-7.2375201869d0,&
                               -7.2847423350d0,-7.2924879904d0,-7.2935801803d0,-7.2897972047d0,-7.2823650212d0,-7.2721420938d0,&
                               -7.2300049498d0,-7.1361504802d0,-7.0179293704d0,-6.8652991322d0,-6.6558110727d0,-6.2631022691d0,&
                               -7.1266893664d0,-7.2574589009d0,-7.3097704631d0,-7.3204310749d0,-7.3246987741d0,-7.3243661927d0,&
                               -7.3206776629d0,-7.3145136985d0,-7.2867102614d0,-7.2261021715d0,-7.1583944794d0,-7.0828874780d0,&
                               -6.9915426198d0,-6.8415160713d0,-7.1319312902d0,-7.2641987272d0,-7.3183562080d0,-7.3300766703d0,&
                               -7.3355008708d0,-7.3364269538d0,-7.3341054989d0,-7.3294243698d0,-7.3068566091d0,-7.2583361073d0,&
                               -7.2088160238d0,-7.1609661972d0,-7.1120951068d0,-7.0441860284d0,-7.1342082784d0,-7.2671530054d0,&
                               -7.3221252188d0,-7.3343068852d0,-7.3402286715d0,-7.3416912946d0,-7.3399471717d0,-7.3358863712d0,&
                               -7.3154735909d0,-7.2718069400d0,-7.2295404384d0,-7.1928925411d0,-7.1615273594d0,-7.1278002985d0,&
                               -7.1347778929d0,-7.2679750878d0,-7.3232067515d0,-7.3355295352d0,-7.3416001563d0,-7.3432193316d0,&
                               -7.3416376406d0,-7.3377473570d0,-7.3178946365d0,-7.2753620697d0,-7.2346438172d0,-7.2003467743d0,&
                               -7.1728178546d0,-7.1471344861d0,-7.1348653563d0,-7.2681114278d0,-7.3234217350d0,-7.3357867805d0,&
                               -7.3419011333d0,-7.3435644078d0,-7.3420293884d0,-7.3381850415d0,-7.3184723359d0,-7.2761819472d0,&
                               -7.2357139575d0,-7.2017344288d0,-7.1747170598d0,-7.1503221226d0,-7.1348403668d0,-7.2681033430d0,&
                               -7.3234437846d0,-7.3358253673d0,-7.3419547873d0,-7.3436338640d0,-7.3421146469d0,-7.3382868371d0,&
                               -7.3186222731d0,-7.2764006057d0,-7.2359888425d0,-7.2020449606d0,-7.1750764682d0,-7.1508513129d0,&
                               -7.1347775254d0,-7.2680559364d0,-7.3234151201d0,-7.3358058902d0,-7.3419456000d0,-7.3436345990d0,&
                               -7.3421253042d0,-7.3383085192d0,-7.3186737222d0,-7.2764921115d0,-7.2360924756d0,-7.2021552086d0,&
                               -7.1751793663d0,-7.1509336313d0,-7.1346952069d0,-7.2679905226d0,-7.3233588937d0,-7.3357551761d0,&
                               -7.3419007658d0,-7.3435956447d0,-7.3420925973d0,-7.3382820597d0,-7.3186663723d0,-7.2765101187d0,&
                               -7.2361288574d0,-7.2021702758d0,-7.1751576842d0,-7.1508505779d0,-7.1346639700d0,-7.2679545082d0,&
                               -7.3233258193d0,-7.3357235717d0,-7.3418706313d0,-7.3435677152d0,-7.3420668727d0,-7.3382578052d0,&
                               -7.3186487326d0,-7.2765012988d0,-7.2361237125d0,-7.2021577810d0,-7.1751290197d0,-7.1507881040d0,&
                               -7.1346239132d0,-7.2679181264d0,-7.3232905399d0,-7.3356901298d0,-7.3418379244d0,-7.3435364783d0,&
                               -7.3420374733d0,-7.3382295082d0,-7.3186248456d0,-7.2764843941d0,-7.2361093803d0,-7.2021383039d0,&
                               -7.1750937403d0,-7.1507212202d0 ],shape( gridPot1416 ) )
   call this%wyckoffSite(1)%zrCut(1)%interRZ%read( gridR14,gridZ16,gridPot1416 )
   ! Second cut2d
   this%wyckoffSite(1)%zrcut(2)%alias='Top_Li_45_0'
   this%wyckoffSite(1)%zrCut(2)%theta=0.785398163397d0
   this%wyckoffSite(1)%zrCut(2)%phi=0.d0
   gridPot1414(:,:)=reshape( [  -4.3930691146d0,-4.5492390478d0,-4.7432240368d0,-4.8587091594d0,&
                                -4.9825143292d0,-5.1107735173d0,-5.2401098281d0,-5.3680834739d0,&
                                -5.7370930050d0,-6.2247672104d0,-6.5372845760d0,-6.7355163171d0,&
                                -6.8625587350d0,-6.9560923830d0,-6.0979945325d0,-6.1737356271d0,&
                                -6.1832555399d0,-6.1793263020d0,-6.1744048323d0,-6.1707174050d0,&
                                -6.1697615551d0,-6.1724810052d0,-6.2061441222d0,-6.3326963041d0,&
                                -6.4939343368d0,-6.6527818258d0,-6.7889858506d0,-6.9054363779d0,&
                                -6.9945082903d0,-7.1105093734d0,-7.1472392217d0,-7.1501071390d0,&
                                -7.1467232611d0,-7.1389735634d0,-7.1281839615d0,-7.1153139802d0,&
                                -7.0705907862d0,-6.9974831482d0,-6.9453715023d0,-6.9210739508d0,&
                                -6.9218636939d0,-6.9449492526d0,-7.1160849811d0,-7.2459383549d0,&
                                -7.2974362870d0,-7.3077808546d0,-7.3118203405d0,-7.3113694263d0,&
                                -7.3076919213d0,-7.3016881840d0,-7.2755064946d0,-7.2230618998d0,&
                                -7.1741801521d0,-7.1342861869d0,-7.1044027380d0,-7.0809349863d0,&
                                -7.1288891810d0,-7.2609574367d0,-7.3149429307d0,-7.3266031241d0,&
                                -7.3319935152d0,-7.3329199657d0,-7.3306393025d0,-7.3260441669d0,&
                                -7.3040971022d0,-7.2583669767d0,-7.2151321304d0,-7.1791530708d0,&
                                -7.1508605002d0,-7.1258511143d0,-7.1329521864d0,-7.2657940154d0,&
                                -7.3206644331d0,-7.3328001629d0,-7.3386826274d0,-7.3401143811d0,&
                                -7.3383496785d0,-7.3342800582d0,-7.3139334267d0,-7.2708294079d0,&
                                -7.2299994374d0,-7.1960985523d0,-7.1694585988d0,-7.1456619406d0,&
                                -7.1346040686d0,-7.2677667192d0,-7.3229965453d0,-7.3353237390d0,&
                                -7.3414024449d0,-7.3430344825d0,-7.3414726362d0,-7.3376066071d0,&
                                -7.3178681770d0,-7.2756876688d0,-7.2355974622d0,-7.2022717039d0,&
                                -7.1761469760d0,-7.1529993109d0,-7.1349616395d0,-7.2682202058d0,&
                                -7.3235529301d0,-7.3359323079d0,-7.3420624628d0,-7.3437452144d0,&
                                -7.3422318772d0,-7.3384128873d0,-7.3187949949d0,-7.2767273072d0,&
                                -7.2365823441d0,-7.2030360899d0,-7.1766008302d0,-7.1531158063d0,&
                                -7.1349311376d0,-7.2682345381d0,-7.3235962943d0,-7.3359900043d0,&
                                -7.3421341240d0,-7.3438297379d0,-7.3423281604d0,-7.3385201953d0,&
                                -7.3189250876d0,-7.2768500499d0,-7.2366176235d0,-7.2028868876d0,&
                                -7.1761833578d0,-7.1523242258d0,-7.1348921833d0,-7.2681808841d0,&
                                -7.3235540326d0,-7.3359528875d0,-7.3421017846d0,-7.3438018084d0,&
                                -7.3423057434d0,-7.3385014532d0,-7.3189125928d0,-7.2768283678d0,&
                                -7.2365544146d0,-7.2027380529d0,-7.1759014905d0,-7.1518387672d0,&
                                -7.1348006775d0,-7.2681037105d0,-7.3234834739d0,-7.3358852687d0,&
                                -7.3420367383d0,-7.3437411720d0,-7.3422473119d0,-7.3384448592d0,&
                                -7.3188607762d0,-7.2767717739d0,-7.2364698912d0,-7.2025943630d0,&
                                -7.1756651924d0,-7.1514495919d0,-7.1347231364d0,-7.2680191870d0,&
                                -7.3234011554d0,-7.3358036852d0,-7.3419569923d0,-7.3436614260d0,&
                                -7.3421694034d0,-7.3383684206d0,-7.3187876451d0,-7.2766964377d0,&
                                -7.2363776504d0,-7.2024650054d0,-7.1754740959d0,-7.1511522898d0,&
                                -7.1346746273d0,-7.2679761903d0,-7.3233585262d0,-7.3357614235d0,&
                                -7.3419147305d0,-7.3436195317d0,-7.3421267741d0,-7.3383261589d0,&
                                -7.3187464858d0,-7.2766545435d0,-7.2363306112d0,-7.2024051040d0,&
                                -7.1753925124d0,-7.1510321195d0,-7.1346327331d0,-7.2679342961d0,&
                                -7.3233155295d0,-7.3357187943d0,-7.3418724688d0,-7.3435769025d0,&
                                -7.3420841449d0,-7.3382838972d0,-7.3187042241d0,-7.2766111793d0,&
                                -7.2362828371d0,-7.2023477750d0,-7.1753186462d0,-7.1509273840d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(1)%zrCut(2)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! third cut2d
   this%wyckoffSite(1)%zrcut(3)%alias='Top_Li_45_45'
   this%wyckoffSite(1)%zrCut(3)%theta=0.785398163397d0
   this%wyckoffSite(1)%zrCut(3)%phi=0.785398163397d0
   gridPot1414(:,:)=reshape( [ -4.3930823443d0,-4.5492658748d0,-4.7432806308d0,-4.8587903754d0,&
                               -4.9826212697d0,-5.1109150022d0,-5.2403027620d0,-5.3683429241d0,&
                               -5.7369761422d0,-6.2233861707d0,-6.5348999123d0,-6.7319990392d0,&
                               -6.8575891237d0,-6.9488557059d0,-6.0980114372d0,-6.1737826663d0,&
                               -6.1833514556d0,-6.1794549246d0,-6.1745764517d0,-6.1709390035d0,&
                               -6.1700426874d0,-6.1728290213d0,-6.2067321114d0,-6.3337142604d0,&
                               -6.4951518419d0,-6.6510373354d0,-6.7852300695d0,-6.8995230440d0,&
                               -6.9945185801d0,-7.1105380379d0,-7.1472983881d0,-7.1501868851d0,&
                               -7.1468305692d0,-7.1391135783d0,-7.1283636657d0,-7.1155414585d0,&
                               -7.0710111985d0,-6.9983882841d0,-6.9468738147d0,-6.9230687042d0,&
                               -6.9239664903d0,-6.9462244542d0,-7.1160890235d0,-7.2459512171d0,&
                               -7.2974627465d0,-7.3078176040d0,-7.3118703196d0,-7.3114355750d0,&
                               -7.3077775472d0,-7.3017980645d0,-7.2757185382d0,-7.2235749204d0,&
                               -7.1751823062d0,-7.1359446840d0,-7.1067690271d0,-7.0839734205d0,&
                               -7.1288910185d0,-7.2609651540d0,-7.3149572629d0,-7.3266240712d0,&
                               -7.3320225471d0,-7.3329581849d0,-7.3306889141d0,-7.3261084782d0,&
                               -7.3042231524d0,-7.2586852259d0,-7.2157932507d0,-7.1803220669d0,&
                               -7.1527163411d0,-7.1286506779d0,-7.1329525539d0,-7.2657984253d0,&
                               -7.3206728855d0,-7.3328119226d0,-7.3386991646d0,-7.3401364307d0,&
                               -7.3383783430d0,-7.3343168075d0,-7.3140047204d0,-7.2710116846d0,&
                               -7.2303930227d0,-7.1968364787d0,-7.1706915387d0,-7.1476894009d0,&
                               -7.1346073760d0,-7.2677689241d0,-7.3230009552d0,-7.3353296189d0,&
                               -7.3414105298d0,-7.3430451398d0,-7.3414862334d0,-7.3376242468d0,&
                               -7.3179016188d0,-7.2757707222d0,-7.2357767989d0,-7.2026167801d0,&
                               -7.1767503999d0,-7.1540587940d0,-7.1349664170d0,-7.2682216758d0,&
                               -7.3235551351d0,-7.3359352478d0,-7.3420668727d0,-7.3437510943d0,&
                               -7.3422395946d0,-7.3384224421d0,-7.3188130021d0,-7.2767673639d0,&
                               -7.2366624576d0,-7.2031856597d0,-7.1768617504d0,-7.1535861977d0,&
                               -7.1349458373d0,-7.2682352731d0,-7.3235981318d0,-7.3359925767d0,&
                               -7.3421374314d0,-7.3438341478d0,-7.3423340403d0,-7.3385271777d0,&
                               -7.3189375823d0,-7.2768754070d0,-7.2366639276d0,-7.2029670012d0,&
                               -7.1763160229d0,-7.1525579515d0,-7.1348844660d0,-7.2681816190d0,&
                               -7.3235555026d0,-7.3359547249d0,-7.3421047245d0,-7.3438058508d0,&
                               -7.3423108883d0,-7.3385062306d0,-7.3189225151d0,-7.2768482125d0,&
                               -7.2365889590d0,-7.2027924419d0,-7.1759863815d0,-7.1519791497d0,&
                               -7.1347981050d0,-7.2681044455d0,-7.3234849438d0,-7.3358867387d0,&
                               -7.3420396782d0,-7.3437437445d0,-7.3422513544d0,-7.3384503716d0,&
                               -7.3188699636d0,-7.2767883111d0,-7.2364963507d0,-7.2026329498d0,&
                               -7.1757184789d0,-7.1515297054d0,-7.1347227689d0,-7.2680199220d0,&
                               -7.3234022579d0,-7.3358055227d0,-7.3419595647d0,-7.3436647334d0,&
                               -7.3421734458d0,-7.3383735655d0,-7.3187957299d0,-7.2767107700d0,&
                               -7.2363993325d0,-7.2024933024d0,-7.1755090077d0,-7.1511956540d0,&
                               -7.1346757298d0,-7.2679769253d0,-7.3233596286d0,-7.3357628935d0,&
                               -7.3419173030d0,-7.3436224717d0,-7.3421315516d0,-7.3383320388d0,&
                               -7.3187545707d0,-7.2766681408d0,-7.2363504559d0,-7.2024300935d0,&
                               -7.1754211768d0,-7.1510629890d0,-7.1346378780d0,-7.2679350311d0,&
                               -7.3233166319d0,-7.3357206317d0,-7.3418746738d0,-7.3435798425d0,&
                               -7.3420889223d0,-7.3382894096d0,-7.3187119415d0,-7.2766240416d0,&
                               -7.2363015793d0,-7.2023705596d0,-7.1753421658d0,-7.1509483311d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(1)%zrCut(3)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Fourth cut2d
   this%wyckoffSite(1)%zrcut(4)%alias='Top_Li_90_0'
   this%wyckoffSite(1)%zrCut(4)%theta=1.57079632679d0
   this%wyckoffSite(1)%zrCut(4)%phi=0.d0
   gridPot1414(:,:)=reshape( [ -5.0184217276d0,-5.2816556754d0,-5.4847890415d0,-5.5759494181d0,&
                               -5.6628281308d0,-5.7462241099d0,-5.8264578096d0,-5.9036012584d0,&
                               -6.1161872860d0,-6.4275860022d0,-6.6500858953d0,-6.7907248286d0,&
                               -6.8786692718d0,-6.9416892199d0,-6.2526276089d0,-6.4004646326d0,&
                               -6.4832391806d0,-6.5138748882d0,-6.5408786600d0,-6.5656657125d0,&
                               -6.5890941425d0,-6.6116655781d0,-6.6763572531d0,-6.7753040767d0,&
                               -6.8581234589d0,-6.9217192690d0,-6.9668467056d0,-7.0012264345d0,&
                               -7.0178705715d0,-7.1470955318d0,-7.1996522120d0,-7.2112785961d0,&
                               -7.2171221063d0,-7.2189841947d0,-7.2181095607d0,-7.2153695310d0,&
                               -7.2013492959d0,-7.1744274751d0,-7.1527652177d0,-7.1371669666d0,&
                               -7.1257772481d0,-7.1153966662d0,-7.1208634958d0,-7.2533448139d0,&
                               -7.3079951032d0,-7.3201080483d0,-7.3260346120d0,-7.3275802887d0,&
                               -7.3260015377d0,-7.3221906326d0,-7.3030306368d0,-7.2630286286d0,&
                               -7.2260657897d0,-7.1959989616d0,-7.1726399879d0,-7.1519269656d0,&
                               -7.1309383234d0,-7.2640822318d0,-7.3193319026d0,-7.3316918032d0,&
                               -7.3378226932d0,-7.3395252894d0,-7.3380542139d0,-7.3343010053d0,&
                               -7.3150167968d0,-7.2739395033d0,-7.2351939546d0,-7.2031580977d0,&
                               -7.1780414037d0,-7.1557624927d0,-7.1339095064d0,-7.2672309140d0,&
                               -7.3226455893d0,-7.3350745786d0,-7.3412620625d0,-7.3430076555d0,&
                               -7.3415637745d0,-7.3378215907d0,-7.3184671910d0,-7.2769338384d0,&
                               -7.2373967092d0,-7.2044329318d0,-7.1784511587d0,-7.1553832397d0,&
                               -7.1349785442d0,-7.2683197965d0,-7.3237454966d0,-7.3361741184d0,&
                               -7.3423564574d0,-7.3440913931d0,-7.3426309749d0,-7.3388656390d0,&
                               -7.3194017263d0,-7.2775390998d0,-7.2374955648d0,-7.2039129288d0,&
                               -7.1773075197d0,-7.1536185371d0,-7.1351152517d0,-7.2684395993d0,&
                               -7.3238465572d0,-7.3362630518d0,-7.3424314261d0,-7.3441505595d0,&
                               -7.3426721341d0,-7.3388869536d0,-7.3193528497d0,-7.2773344060d0,&
                               -7.2370821349d0,-7.2032411511d0,-7.1763524047d0,-7.1523444380d0,&
                               -7.1350259509d0,-7.2683414786d0,-7.3237399842d0,-7.3361513338d0,&
                               -7.3423134607d0,-7.3440256118d0,-7.3425405715d0,-7.3387473062d0,&
                               -7.3191871102d0,-7.2771164825d0,-7.2368010026d0,-7.2028839477d0,&
                               -7.1759099429d0,-7.1517972405d0,-7.1349362825d0,-7.2682466653d0,&
                               -7.3236414960d0,-7.3360502732d0,-7.3422101951d0,-7.3439201412d0,&
                               -7.3424317935d0,-7.3386363232d0,-7.3190665725d0,-7.2769786726d0,&
                               -7.2366437155d0,-7.2027038760d0,-7.1757034117d0,-7.1515587374d0,&
                               -7.1348275045d0,-7.2681433997d0,-7.3235349229d0,-7.3359418627d0,&
                               -7.3420999471d0,-7.3438080558d0,-7.3423178706d0,-7.3385205628d0,&
                               -7.3189456672d0,-7.2768478450d0,-7.2365037006d0,-7.2025528363d0,&
                               -7.1755406121d0,-7.1513805031d0,-7.1347371012d0,-7.2680408691d0,&
                               -7.3234283499d0,-7.3358338197d0,-7.3419904341d0,-7.3436967053d0,&
                               -7.3422050502d0,-7.3384059049d0,-7.3188273343d0,-7.2767239997d0,&
                               -7.2363747104d0,-7.2024201712d0,-7.1754028022d0,-7.1512368133d0,&
                               -7.1346856521d0,-7.2679916250d0,-7.3233776358d0,-7.3357820031d0,&
                               -7.3419375151d0,-7.3436430513d0,-7.3421506612d0,-7.3383511485d0,&
                               -7.3187703729d0,-7.2766652008d0,-7.2363144415d0,-7.2023580649d0,&
                               -7.1753406958d0,-7.1511736044d0,-7.1346474328d0,-7.2679438509d0,&
                               -7.3233276567d0,-7.3357312890d0,-7.3418856986d0,-7.3435904998d0,&
                               -7.3420977422d0,-7.3382971269d0,-7.3187148814d0,-7.2766078719d0,&
                               -7.2362560101d0,-7.2022988984d0,-7.1752804269d0,-7.1511129681d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(1)%zrCut(4)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Fifth cut2d
   this%wyckoffSite(1)%zrcut(5)%alias='Top_Li_90_45'
   this%wyckoffSite(1)%zrCut(5)%theta=1.57079632679d0
   this%wyckoffSite(1)%zrCut(5)%phi=0.785398163397d0
   gridpot1414(:,:)=reshape( [ -5.0184786890d0,-5.2817479162d0,-5.4849353038d0,-5.5761136876d0,&
                               -5.6629935027d0,-5.7463681673d0,-5.8265401280d0,-5.9035718590d0,&
                               -6.1151924818d0,-6.4139604549d0,-6.6235473700d0,-6.7500183359d0,&
                               -6.8211885494d0,-6.8491095844d0,-6.2527040475d0,-6.4006380894d0,&
                               -6.4835713945d0,-6.5143088978d0,-6.5414254900d0,-6.5663378577d0,&
                               -6.5898974827d0,-6.6125960710d0,-6.6775887230d0,-6.7761860605d0,&
                               -6.8564260076d0,-6.9137641425d0,-6.9476279109d0,-6.9587846386d0,&
                               -7.0179179781d0,-7.1472057798d0,-7.1998734430d0,-7.2115773681d0,&
                               -7.2175142216d0,-7.2194894979d0,-7.2187482640d0,-7.2161600090d0,&
                               -7.2027163708d0,-7.1770506419d0,-7.1566264693d0,-7.1416062851d0,&
                               -7.1294477708d0,-7.1155524833d0,-7.1208844430d0,-7.2533918530d0,&
                               -7.3080946939d0,-7.3202432859d0,-7.3262128463d0,-7.3278125444d0,&
                               -7.3262984722d0,-7.3225621683d0,-7.3037049869d0,-7.2644835344d0,&
                               -7.2285959807d0,-7.1997238732d0,-7.1774140928d0,-7.1572916321d0,&
                               -7.1309511857d0,-7.2641105288d0,-7.3193892315d0,-7.3317697118d0,&
                               -7.3379248563d0,-7.3396586894d0,-7.3382250983d0,-7.3345141514d0,&
                               -7.3154041347d0,-7.2747895152d0,-7.2367223590d0,-7.2055085845d0,&
                               -7.1812091956d0,-7.1595711928d0,-7.1339168562d0,-7.2672470837d0,&
                               -7.3226782962d0,-7.3351190453d0,-7.3413201265d0,-7.3430833591d0,&
                               -7.3416607927d0,-7.3379417610d0,-7.3186843795d0,-7.2774068022d0,&
                               -7.2382533359d0,-7.2057687697d0,-7.1802794377d0,-7.1575665171d0,&
                               -7.1349825867d0,-7.2683286163d0,-7.3237627687d0,-7.3361976380d0,&
                               -7.3423862244d0,-7.3441299799d0,-7.3426805865d0,-7.3389262754d0,&
                               -7.3195082993d0,-7.2777621682d0,-7.2378854752d0,-7.2045034905d0,&
                               -7.1780866054d0,-7.1544086476d0,-7.1351178242d0,-7.2684443767d0,&
                               -7.3238572145d0,-7.3362773840d0,-7.3424494332d0,-7.3441737116d0,&
                               -7.3427019011d0,-7.3389229680d0,-7.3194142211d0,-7.2774560463d0,&
                               -7.2372805813d0,-7.2035167711d0,-7.1766680814d0,-7.1524800430d0,&
                               -7.1350156611d0,-7.2683458885d0,-7.3237488040d0,-7.3361627261d0,&
                               -7.3423277930d0,-7.3440443539d0,-7.3425637236d0,-7.3387752357d0,&
                               -7.3192334144d0,-7.2772046809d0,-7.2369377101d0,-7.2030588745d0,&
                               -7.1760749473d0,-7.1517145545d0,-7.1349355475d0,-7.2682503403d0,&
                               -7.3236492133d0,-7.3360605630d0,-7.3422230574d0,-7.3439363109d0,&
                               -7.3424527407d0,-7.3386609453d0,-7.3191073642d0,-7.2770540087d0,&
                               -7.2367572709d0,-7.2028416860d0,-7.1758143946d0,-7.1513963054d0,&
                               -7.1348341194d0,-7.2681467072d0,-7.3235422728d0,-7.3359514175d0,&
                               -7.3421117069d0,-7.3438231230d0,-7.3423373478d0,-7.3385429799d0,&
                               -7.3189827840d0,-7.2769158312d0,-7.2366029237d0,-7.2026671267d0,&
                               -7.1756166832d0,-7.1511680920d0,-7.1347374687d0,-7.2680438091d0,&
                               -7.3234353323d0,-7.3358430070d0,-7.3420014589d0,-7.3437106701d0,&
                               -7.3422234249d0,-7.3384275870d0,-7.3188618787d0,-7.2767868411d0,&
                               -7.2364643788d0,-7.2025193944d0,-7.1754571912d0,-7.1509931653d0,&
                               -7.1346878571d0,-7.2679949325d0,-7.3233842507d0,-7.3357908230d0,&
                               -7.3419485399d0,-7.3436570161d0,-7.3421686684d0,-7.3383720956d0,&
                               -7.3188045498d0,-7.2767262047d0,-7.2364008024d0,-7.2024525106d0,&
                               -7.1753870000d0,-7.1509185641d0,-7.1346459628d0,-7.2679471583d0,&
                               -7.3233342716d0,-7.3357401089d0,-7.3418967234d0,-7.3436044645d0,&
                               -7.3421150144d0,-7.3383180741d0,-7.3187483233d0,-7.2766670383d0,&
                               -7.2363394311d0,-7.2023893018d0,-7.1753212187d0,-7.1508494754d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(1)%zrCut(5)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Wickoff Top F -----------> 'b'
   this%wyckoffSite(2)%id='b'
   this%wyckoffSite(2)%myNumber=2
   this%wyckoffSite(2)%is_homonucl=.true.
   this%wyckoffSite(2)%x=2.7216780628885480d0
   this%wyckoffSite(2)%y=2.7216780628885480d0
   this%wyckoffSite(2)%n2dcuts=5
   allocate( this%wyckoffSite(2)%nPhiPoints(3) )
   this%wyckoffSite(2)%nPhiPoints(:)=[1,2,2]
   allocate( this%wyckoffSite(2)%zrCut(5) )
   this%wyckoffSite(2)%zrCut(:)%x=2.7216780628885480d0
   this%wyckoffSite(2)%zrCut(:)%y=2.7216780628885480d0
   ! reading zrcuts
   this%wyckoffSite(2)%zrcut(1)%alias='Top_F_0_0'
   this%wyckoffSite(2)%zrCut(1)%theta=0.d0
   this%wyckoffSite(2)%zrCut(1)%phi=0.d0
   gridPot1416(:,:)=reshape( [ 33.4715168470d0,36.7280777391d0,32.7111077495d0,27.9642199020d0,&
                               21.3856928713d0,13.9533148245d0, 9.0618268879d0, 5.5913042852d0,&
                               -0.5049412292d0,-4.6567383910d0,-6.2356988133d0,-6.8413948663d0,&
                               -7.0578498187d0,-7.1242317001d0,-1.9033033846d0, 0.0034275311d0,&
                                3.2801723154d0, 5.7512859620d0, 9.1276728621d0,13.9376318364d0,&
                               20.9998533972d0,27.4752624360d0,35.7943159572d0,13.6062627973d0,&
                               -0.6447695306d0,-4.7270460195d0,-6.2813546984d0,-6.9360302369d0,&
                               -5.9031036726d0,-5.7674850646d0,-5.4007425994d0,-5.1231609799d0,&
                               -4.7715419645d0,-4.3329174513d0,-3.7898079910d0,-3.1189144748d0,&
                                0.0390016364d0,16.3720911156d0,35.6009409801d0,11.3358053823d0,&
                               -1.1717378205d0,-5.3378967244d0,-6.7182110598d0,-6.8361484648d0,&
                               -6.8604724758d0,-6.8476532086d0,-6.8195539394d0,-6.7756550328d0,&
                               -6.7145070953d0,-6.6338603332d0,-6.2405863250d0,-4.7196967512d0,&
                               -0.6132108564d0,13.6361947599d0,35.7370124293d0, 8.7378711640d0,&
                               -6.9382990946d0,-7.0706065884d0,-7.1257691633d0,-7.1377939101d0,&
                               -7.1429116211d0,-7.1424430672d0,-7.1370585561d0,-7.1269477142d0,&
                               -7.0648972132d0,-6.7898784917d0,-6.0360156927d0,-4.0851696691d0,&
                                1.2479442124d0,28.4739811097d0,-7.0611300399d0,-7.1933041336d0,&
                               -7.2486306104d0,-7.2614649448d0,-7.2683367012d0,-7.2710050697d0,&
                               -7.2706780007d0,-7.2681731667d0,-7.2523099529d0,-7.2062387936d0,&
                               -7.1133137144d0,-6.8802009836d0,-6.2704458868d0,-4.1779815604d0,&
                               -7.1104149276d0,-7.2431556959d0,-7.2984219039d0,-7.3109839258d0,&
                               -7.3174554820d0,-7.3196405969d0,-7.3187920550d0,-7.3157947800d0,&
                               -7.2994578674d0,-7.2642814131d0,-7.2290340327d0,-7.1898500645d0,&
                               -7.1240533373d0,-6.9118171633d0,-7.1272883804d0,-7.2604888828d0,&
                               -7.3159524346d0,-7.3284905695d0,-7.3348522452d0,-7.3368418537d0,&
                               -7.3357165893d0,-7.3323687258d0,-7.3146584909d0,-7.2771436770d0,&
                               -7.2422424752d0,-7.2129385631d0,-7.1874631958d0,-7.1537780291d0,&
                               -7.1336412363d0,-7.2670148280d0,-7.3225349738d0,-7.3350418717d0,&
                               -7.3413252714d0,-7.3431884622d0,-7.3418853311d0,-7.3383096217d0,&
                               -7.3196236922d0,-7.2798127805d0,-7.2426720748d0,-7.2122925100d0,&
                               -7.1884785796d0,-7.1664249420d0,-7.1348840985d0,-7.2682492378d0,&
                               -7.3237212420d0,-7.3361833057d0,-7.3424038641d0,-7.3441884113d0,&
                               -7.3427845871d0,-7.3390865025d0,-7.3199011496d0,-7.2788267961d0,&
                               -7.2400467030d0,-7.2081508610d0,-7.1834244449d0,-7.1616776641d0,&
                               -7.1349855266d0,-7.2683315563d0,-7.3237682811d0,-7.3362064578d0,&
                               -7.3423979842d0,-7.3441465171d0,-7.3427019011d0,-7.3389549399d0,&
                               -7.3195733457d0,-7.2779903815d0,-7.2384620721d0,-7.2056768964d0,&
                               -7.1801070833d0,-7.1576940372d0,-7.1349381200d0,-7.2682712874d0,&
                               -7.3236947825d0,-7.3361211994d0,-7.3423005985d0,-7.3440333291d0,&
                               -7.3425718085d0,-7.3388053701d0,-7.3193488073d0,-7.2775648243d0,&
                               -7.2377241456d0,-7.2045200277d0,-7.1784779857d0,-7.1555412618d0,&
                               -7.1348473491d0,-7.2681797816d0,-7.3235911494d0,-7.3360109514d0,&
                               -7.3421826331d0,-7.3439058090d0,-7.3424328960d0,-7.3386543304d0,&
                               -7.3191496259d0,-7.2772366528d0,-7.2371945879d0,-7.2036975778d0,&
                               -7.1772895125d0,-7.1538857047d0,-7.1347558433d0,-7.2680798234d0,&
                               -7.3234820039d0,-7.3358970285d0,-7.3420624628d0,-7.3437790238d0,&
                               -7.3422991285d0,-7.3385124780d0,-7.3189794766d0,-7.2769889624d0,&
                               -7.2368168048d0,-7.2031268607d0,-7.1764615502d0,-7.1526825317d0,&
                               -7.1347099067d0,-7.2680298443d0,-7.3234283499d0,-7.3358411695d0,&
                               -7.3420047664d0,-7.3437187549d0,-7.3422362871d0,-7.3384463292d0,&
                               -7.3189026705d0,-7.2768838593d0,-7.2366698075d0,-7.2029122447d0,&
                               -7.1761499160d0,-7.1522239002d0,-7.1346647050d0,-7.2679820702d0,&
                               -7.3233761658d0,-7.3357878830d0,-7.3419485399d0,-7.3436606910d0,&
                               -7.3421760182d0,-7.3383834879d0,-7.3188310093d0,-7.2767927210d0,&
                               -7.2365433898d0,-7.2027325405d0,-7.1758937732d0,-7.1518435446d0 ],shape( gridPot1416 ) )
   call this%wyckoffSite(2)%zrCut(1)%interRZ%read( gridR14,gridZ16,gridPot1416 )
   ! Second cut2d
   this%wyckoffSite(2)%zrcut(2)%alias='Top_F_45_0'
   this%wyckoffSite(2)%zrCut(2)%theta=0.785398163397d0
   this%wyckoffSite(2)%zrCut(2)%phi=0.d0
   gridPot1414(:,:)=reshape( [   7.3744382510d0,   5.3167479737d0,   2.8013896828d0,   1.5772557258d0,&
                               0.4404425759d0,  -0.5879814523d0,  -1.5016393216d0,  -2.3030362705d0,&
                              -4.1194461324d0,  -5.7720787303d0,  -6.5066907626d0,  -6.8351525581d0,&
                              -6.9866299699d0,  -7.0681296839d0,  -3.9607838600d0,  -3.8982769324d0,&
                              -3.8325867783d0,  -3.8276432590d0,  -3.8507247753d0,  -3.9053151632d0,&
                              -3.9921137623d0,  -4.1093819622d0,  -4.6009152142d0,  -5.5587819707d0,&
                              -6.2990761838d0,  -6.7339651281d0,  -6.9516990012d0,  -7.0609808376d0,&
                              -6.7478251362d0,  -6.8910640842d0,  -6.9587118749d0,  -6.9776399825d0,&
                              -6.9906194767d0,  -6.9990527119d0,  -7.0041998224d0,  -7.0068615761d0,&
                              -7.0055080984d0,  -6.9933573015d0,  -6.9929074897d0,  -7.0134132458d0,&
                              -7.0455703755d0,  -7.0787094472d0,  -7.0617937327d0,  -7.1940872617d0,&
                              -7.2495136967d0,  -7.2623947028d0,  -7.2693355479d0,  -7.2721314365d0,&
                              -7.2720208211d0,  -7.2698724555d0,  -7.2566155039d0,  -7.2277540537d0,&
                              -7.2007793138d0,  -7.1777385893d0,  -7.1580736578d0,  -7.1380107311d0,&
                              -7.1098879423d0,  -7.2422928218d0,  -7.2970941508d0,  -7.3093669555d0,&
                              -7.3155092378d0,  -7.3173279619d0,  -7.3160773823d0,  -7.3126449953d0,&
                              -7.2948961736d0,  -7.2577554679d0,  -7.2234165308d0,  -7.1949799028d0,&
                              -7.1719230085d0,  -7.1496892992d0,  -7.1267279532d0,  -7.2596098389d0,&
                              -7.3146621658d0,  -7.3269548152d0,  -7.3330412384d0,  -7.3347261950d0,&
                              -7.3332643068d0,  -7.3295471126d0,  -7.3105451389d0,  -7.2704394976d0,&
                              -7.2329629030d0,  -7.2020243810d0,  -7.1774750966d0,  -7.1548195050d0,&
                              -7.1333101249d0,  -7.2665032774d0,  -7.3217944749d0,  -7.3341665028d0,&
                              -7.3403003327d0,  -7.3419981515d0,  -7.3405131113d0,  -7.3367356481d0,&
                              -7.3173297993d0,  -7.2759213945d0,  -7.2367550659d0,  -7.2042899769d0,&
                              -7.1787583831d0,  -7.1558448112d0,  -7.1347264439d0,  -7.2680044873d0,&
                              -7.3233695510d0,  -7.3357684059d0,  -7.3419206104d0,  -7.3436276166d0,&
                              -7.3421414739d0,  -7.3383522509d0,  -7.3188376242d0,  -7.2769896974d0,&
                              -7.2371093294d0,  -7.2038273029d0,  -7.1775842421d0,  -7.1541616921d0,&
                              -7.1348980632d0,  -7.2681955838d0,  -7.3235768172d0,  -7.3359797145d0,&
                              -7.3421352265d0,  -7.3438433351d0,  -7.3423549875d0,  -7.3385609871d0,&
                              -7.3190118160d0,  -7.2770249767d0,  -7.2368958159d0,  -7.2032698156d0,&
                              -7.1766493393d0,  -7.1528243841d0,  -7.1348789536d0,  -7.2681801491d0,&
                              -7.3235657924d0,  -7.3359697922d0,  -7.3421256717d0,  -7.3438334128d0,&
                              -7.3423454326d0,  -7.3385488598d0,  -7.3189871939d0,  -7.2769525805d0,&
                              -7.2367322814d0,  -7.2029684712d0,  -7.1761756405d0,  -7.1521342318d0,&
                              -7.1348186847d0,  -7.2681184102d0,  -7.3235047885d0,  -7.3359095233d0,&
                              -7.3420657703d0,  -7.3437738789d0,  -7.3422844287d0,  -7.3384871210d0,&
                              -7.3189184727d0,  -7.2768566648d0,  -7.2365819766d0,  -7.2027318055d0,&
                              -7.1758228470d0,  -7.1516179038d0,  -7.1347345287d0,  -7.2680371942d0,&
                              -7.3234232050d0,  -7.3358283073d0,  -7.3419838193d0,  -7.3436904579d0,&
                              -7.3422002728d0,  -7.3384029650d0,  -7.3188302743d0,  -7.2767533992d0,&
                              -7.2364478416d0,  -7.2025462214d0,  -7.1755633967d0,  -7.1512434282d0,&
                              -7.1346941045d0,  -7.2679941975d0,  -7.3233794733d0,  -7.3357842081d0,&
                              -7.3419397201d0,  -7.3436459913d0,  -7.3421550711d0,  -7.3383566609d0,&
                              -7.3187836027d0,  -7.2767015826d0,  -7.2363861027d0,  -7.2024675778d0,&
                              -7.1754597636d0,  -7.1510997383d0,  -7.1346485353d0,  -7.2679519358d0,&
                              -7.3233353741d0,  -7.3357401089d0,  -7.3418948859d0,  -7.3436007896d0,&
                              -7.3421095020d0,  -7.3383103567d0,  -7.3187358285d0,  -7.2766497661d0,&
                              -7.2363273038d0,  -7.2023962842d0,  -7.1753686253d0,  -7.1509755256d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(2)%zrCut(2)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Third cut2d
   this%wyckoffSite(2)%zrcut(3)%alias='Top_F_45_45'
   this%wyckoffSite(2)%zrCut(3)%theta=0.785398163397d0
   this%wyckoffSite(2)%zrCut(3)%phi=0.785398163397d0
   gridPot1414(:,:)=reshape( [ 7.3744812477d0, 5.3168402145d0, 2.8015803383d0, 1.5775181527d0,&
                               0.4407941934d0,-0.5875210935d0,-1.5010470327d0,-2.3022891935d0,&
                              -4.1180518630d0,-5.7690799853d0,-6.5014565562d0,-6.8269313665d0,&
                              -6.9744773355d0,-7.0498873187d0,-3.9607515206d0,-3.8982133561d0,&
                              -3.8324555832d0,-3.8274635548d0,-3.8504888446d0,-3.9050020589d0,&
                              -3.9917109897d0,-4.1088751890d0,-4.5999582617d0,-5.5565281346d0,&
                              -6.2944865605d0,-6.7255278505d0,-6.9374490828d0,-7.0368005165d0,&
                              -6.7478137439d0,-6.8910317448d0,-6.9586633658d0,-6.9776921665d0,&
                              -6.9905220910d0,-6.9989365840d0,-7.0040340830d0,-7.0066451225d0,&
                              -7.0050560817d0,-6.9920865098d0,-6.9898841227d0,-7.0070067359d0,&
                              -7.0331186016d0,-7.0534803003d0,-7.0617874853d0,-7.1940758694d0,&
                              -7.2494905447d0,-7.2623634658d0,-7.2692940211d0,-7.2720770475d0,&
                              -7.2719502624d0,-7.2697824197d0,-7.2564380046d0,-7.2272704325d0,&
                              -7.1995985580d0,-7.1750970478d0,-7.1525983758d0,-7.1259330653d0,&
                              -7.1098838999d0,-7.2422865744d0,-7.2970812885d0,-7.3093489483d0,&
                              -7.3154853507d0,-7.3172970924d0,-7.3160369581d0,-7.3125939138d0,&
                              -7.2947973180d0,-7.2575015301d0,-7.2228380964d0,-7.1937462279d0,&
                              -7.1694189095d0,-7.1441272888d0,-7.1267253807d0,-7.2596065315d0,&
                              -7.3146551834d0,-7.3269452603d0,-7.3330291111d0,-7.3347100253d0,&
                              -7.3332440947d0,-7.3295213880d0,-7.3104933223d0,-7.2703072000d0,&
                              -7.2326722159d0,-7.2014393317d0,-7.1763524047d0,-7.1524418237d0,&
                              -7.1333090224d0,-7.2665021749d0,-7.3217922699d0,-7.3341635628d0,&
                              -7.3402962902d0,-7.3419930066d0,-7.3405061289d0,-7.3367264608d0,&
                              -7.3173106897d0,-7.2758703129d0,-7.2366400405d0,-7.2040621311d0,&
                              -7.1783401758d0,-7.1550286087d0,-7.1347257089d0,-7.2680041198d0,&
                              -7.3233691835d0,-7.3357676709d0,-7.3419202429d0,-7.3436268816d0,&
                              -7.3421407389d0,-7.3383507810d0,-7.3188335817d0,-7.2769746301d0,&
                              -7.2370689052d0,-7.2037402070d0,-7.1774177677d0,-7.1538397680d0,&
                              -7.1348833635d0,-7.2681959513d0,-7.3235771847d0,-7.3359800820d0,&
                              -7.3421359615d0,-7.3438440701d0,-7.3423560899d0,-7.3385620896d0,&
                              -7.3190129184d0,-7.2770224043d0,-7.2368818511d0,-7.2032319638d0,&
                              -7.1765673883d0,-7.1526542348d0,-7.1348833635d0,-7.2681805166d0,&
                              -7.3235661599d0,-7.3359701597d0,-7.3421267741d0,-7.3438348828d0,&
                              -7.3423469026d0,-7.3385510648d0,-7.3189901339d0,-7.2769547855d0,&
                              -7.2367286064d0,-7.2029489940d0,-7.1761256614d0,-7.1520210439d0,&
                              -7.1348139073d0,-7.2681187777d0,-7.3235055235d0,-7.3359106258d0,&
                              -7.3420672402d0,-7.3437760839d0,-7.3422870012d0,-7.3384900609d0,&
                              -7.3189228826d0,-7.2768618097d0,-7.2365845491d0,-7.2027244556d0,&
                              -7.1757938150d0,-7.1515414652d0,-7.1347334262d0,-7.2680379292d0,&
                              -7.3234239400d0,-7.3358294097d0,-7.3419856567d0,-7.3436926629d0,&
                              -7.3422028452d0,-7.3384051700d0,-7.3188354192d0,-7.2767603816d0,&
                              -7.2364544565d0,-7.2025462214d0,-7.1755472270d0,-7.1511905091d0,&
                              -7.1346918995d0,-7.2679949325d0,-7.3233802083d0,-7.3357853106d0,&
                              -7.3419415575d0,-7.3436481962d0,-7.3421580111d0,-7.3383599683d0,&
                              -7.3187887476d0,-7.2767093000d0,-7.2363941876d0,-7.2024701503d0,&
                              -7.1754480038d0,-7.1510541691d0,-7.1346503728d0,-7.2679526707d0,&
                              -7.3233364766d0,-7.3357412114d0,-7.3418967234d0,-7.3436029945d0,&
                              -7.3421124419d0,-7.3383140316d0,-7.3187413409d0,-7.2766582184d0,&
                              -7.2363364911d0,-7.2024006941d0,-7.1753605404d0,-7.1509373063d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(2)%zrCut(3)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Fourth cut2d
   this%wyckoffSite(2)%zrcut(4)%alias='Top_F_90_0'
   this%wyckoffSite(2)%zrCut(4)%theta=1.57079632679d0
   this%wyckoffSite(2)%zrCut(4)%phi=0.d0
   gridPot1414(:,:)=reshape( [ 0.9741794146d0,-0.9041065359d0,-2.4291775949d0,-3.0667038997d0,&
                              -3.6277203507d0,-4.1181084570d0,-4.5440853224d0,-4.9120283881d0,&
                              -5.7275583925d0,-6.4391914392d0,-6.7431337173d0,-6.8884126204d0,&
                              -6.9689340673d0,-7.0192523461d0,-5.1308742958d0,-5.5135512685d0,&
                              -5.8165858383d0,-5.9475990208d0,-6.0670850424d0,-6.1758244588d0,&
                              -6.2743773297d0,-6.3632111065d0,-6.5765523327d0,-6.7934101019d0,&
                              -6.9072381973d0,-6.9740712555d0,-7.0171866665d0,-7.0473416930d0,&
                              -6.7654780421d0,-6.9156820898d0,-6.9907973434d0,-7.0138200608d0,&
                              -7.0309430415d0,-7.0438607969d0,-7.0537309307d0,-7.0613593557d0,&
                              -7.0757500240d0,-7.0863525719d0,-7.0913372504d0,-7.0948104291d0,&
                              -7.0976735690d0,-7.1001504736d0,-7.0621351339d0,-7.1942206617d0,&
                              -7.2490638850d0,-7.2614899344d0,-7.2678648399d0,-7.2699841735d0,&
                              -7.2690941048d0,-7.2660762502d0,-7.2498238610d0,-7.2158435973d0,&
                              -7.1855169515d0,-7.1621311508d0,-7.1452110264d0,-7.1316078961d0,&
                              -7.1093502997d0,-7.2414038556d0,-7.2957072313d0,-7.3076665642d0,&
                              -7.3134498056d0,-7.3148595097d0,-7.3131506661d0,-7.3092151808d0,&
                              -7.2897101088d0,-7.2492542465d0,-7.2121241981d0,-7.1824083261d0,&
                              -7.1600621638d0,-7.1413644745d0,-7.1261903106d0,-7.2587840816d0,&
                              -7.3134806750d0,-7.3255671606d0,-7.3314275755d0,-7.3328633717d0,&
                              -7.3311299060d0,-7.3271194521d0,-7.3071039321d0,-7.2649770779d0,&
                              -7.2254700831d0,-7.1930773902d0,-7.1680977713d0,-7.1466108082d0,&
                              -7.1329995931d0,-7.2660317835d0,-7.3211340895d0,-7.3333991769d0,&
                              -7.3394198188d0,-7.3409956299d0,-7.3393801296d0,-7.3354644889d0,&
                              -7.3155897188d0,-7.2732162766d0,-7.2329114540d0,-7.1993148532d0,&
                              -7.1729170778d0,-7.1497058364d0,-7.1345768741d0,-7.2677766415d0,&
                              -7.3230516693d0,-7.3354009126d0,-7.3415009331d0,-7.3431520803d0,&
                              -7.3416067712d0,-7.3377554419d0,-7.3180350189d0,-7.2757659448d0,&
                              -7.2353582241d0,-7.2014639538d0,-7.1746167341d0,-7.1507576021d0,&
                              -7.1348091298d0,-7.2680673286d0,-7.3233982154d0,-7.3357735508d0,&
                              -7.3419003983d0,-7.3435776375d0,-7.3420569504d0,-7.3382295082d0,&
                              -7.3185711915d0,-7.2763664288d0,-7.2359642205d0,-7.2020126212d0,&
                              -7.1750529486d0,-7.1510089675d0,-7.1348161122d0,-7.2680926857d0,&
                              -7.3234426821d0,-7.3358286748d0,-7.3419650771d0,-7.3436518711d0,&
                              -7.3421407389d0,-7.3383224840d0,-7.3186887894d0,-7.2765130586d0,&
                              -7.2361204051d0,-7.2021570460d0,-7.1751646665d0,-7.1510622540d0,&
                              -7.1347620907d0,-7.2680577738d0,-7.3234195300d0,-7.3358114026d0,&
                              -7.3419544198d0,-7.3436470937d0,-7.3421414739d0,-7.3383290989d0,&
                              -7.3187119415d0,-7.2765564228d0,-7.2361729566d0,-7.2022062901d0,&
                              -7.1751981084d0,-7.1510644589d0,-7.1347054967d0,-7.2679938300d0,&
                              -7.3233622011d0,-7.3357573811d0,-7.3419037057d0,-7.3435996871d0,&
                              -7.3420981097d0,-7.3382894096d0,-7.3186825420d0,-7.2765409881d0,&
                              -7.2361656068d0,-7.2021989403d0,-7.1751841437d0,-7.1510343245d0,&
                              -7.1346720549d0,-7.2679574482d0,-7.3233272892d0,-7.3357235717d0,&
                              -7.3418709988d0,-7.3435677152d0,-7.3420672402d0,-7.3382596426d0,&
                              -7.3186560825d0,-7.2765200410d0,-7.2361475996d0,-7.2021820356d0,&
                              -7.1751650340d0,-7.1510108049d0,-7.1346327331d0,-7.2679192289d0,&
                              -7.3232898049d0,-7.3356868224d0,-7.3418346170d0,-7.3435324358d0,&
                              -7.3420326959d0,-7.3382258333d0,-7.3186248456d0,-7.2764924790d0,&
                              -7.2361229775d0,-7.2021581485d0,-7.1751404120d0,-7.1509828754d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(2)%zrCut(4)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Fifth cut2d
   this%wyckoffSite(2)%zrcut(5)%alias='Top_F_90_45'
   this%wyckoffSite(2)%zrCut(5)%theta=1.57079632679d0
   this%wyckoffSite(2)%zrCut(5)%phi=0.785398163397d0
   gridPot1414(:,:)=reshape( [ 0.9743341660d0,-0.9037241960d0,-2.4283785544d0,-3.0655893662d0,&
                              -3.6261854049d0,-4.1160920215d0,-4.5414735478d0,-4.9086048210d0,&
                              -5.7208854500d0,-6.4226847447d0,-6.7098175139d0,-6.8282418475d0,&
                              -6.8621427326d0,-6.8124943940d0,-5.1307636803d0,-5.5132782211d0,&
                              -5.8160221036d0,-5.9468214051d0,-6.0660284993d0,-6.1744349668d0,&
                              -6.2725773477d0,-6.3609102312d0,-6.5721475585d0,-6.7828579006d0,&
                              -6.8859570305d0,-6.9347939440d0,-6.9476834024d0,-6.9162590542d0,&
                              -6.7654368829d0,-6.9156012413d0,-6.9906360139d0,-7.0136025048d0,&
                              -7.0306549268d0,-7.0434877912d0,-7.0532557619d0,-7.0607603416d0,&
                              -7.0746335795d0,-7.0836871433d0,-7.0857983920d0,-7.0843280515d0,&
                              -7.0792331250d0,-7.0668537473d0,-7.0621171267d0,-7.1941787675d0,&
                              -7.2489830365d0,-7.2613818914d0,-7.2677233550d0,-7.2698044693d0,&
                              -7.2688695664d0,-7.2657998953d0,-7.2493420774d0,-7.2148443831d0,&
                              -7.1837232169d0,-7.1591739326d0,-7.1405813464d0,-7.1241062563d0,&
                              -7.1093403774d0,-7.2413807035d0,-7.2956623971d0,-7.3076062953d0,&
                              -7.3133700595d0,-7.3147577141d0,-7.3130235134d0,-7.3090564237d0,&
                              -7.2894308140d0,-7.2486750771d0,-7.2111055068d0,-7.1807983381d0,&
                              -7.1576752951d0,-7.1376660224d0,-7.1261859006d0,-7.2587734242d0,&
                              -7.3134589929d0,-7.3255377612d0,-7.3313875188d0,-7.3328122901d0,&
                              -7.3310663297d0,-7.3270397061d0,-7.3069613447d0,-7.2646738959d0,&
                              -7.2249269281d0,-7.1922079012d0,-7.1667964777d0,-7.1445263865d0,&
                              -7.1329981231d0,-7.2660284761d0,-7.3211282096d0,-7.3333910920d0,&
                              -7.3394080591d0,-7.3409805627d0,-7.3393617549d0,-7.3354409694d0,&
                              -7.3155463546d0,-7.2731166860d0,-7.2327218275d0,-7.1989852118d0,&
                              -7.1723713503d0,-7.1486489258d0,-7.1345765066d0,-7.2677773765d0,&
                              -7.3230538743d0,-7.3354020151d0,-7.3415016681d0,-7.3431531828d0,&
                              -7.3416086087d0,-7.3377565444d0,-7.3180357539d0,-7.2757571250d0,&
                              -7.2353266197d0,-7.2013786953d0,-7.1744142453d0,-7.1501699804d0,&
                              -7.1348157447d0,-7.2680691661d0,-7.3234007879d0,-7.3357779607d0,&
                              -7.3419055432d0,-7.3435842524d0,-7.3420657703d0,-7.3382394305d0,&
                              -7.3185858913d0,-7.2763888459d0,-7.2359862700d0,-7.2020104162d0,&
                              -7.1749673227d0,-7.1505808378d0,-7.1348238296d0,-7.2680945231d0,&
                              -7.3234470920d0,-7.3358345546d0,-7.3419720595d0,-7.3436606910d0,&
                              -7.3421521312d0,-7.3383353462d0,-7.3187093690d0,-7.2765468680d0,&
                              -7.2361619318d0,-7.2021842405d0,-7.1751205674d0,-7.1506896158d0,&
                              -7.1347800979d0,-7.2680599788d0,-7.3234246750d0,-7.3358180175d0,&
                              -7.3419621372d0,-7.3436570161d0,-7.3421547036d0,-7.3383441661d0,&
                              -7.3187358285d0,-7.2765975821d0,-7.2362266106d0,-7.2022511243d0,&
                              -7.1751782638d0,-7.1507241602d0,-7.1347143166d0,-7.2679964024d0,&
                              -7.3233673460d0,-7.3357647309d0,-7.3419121581d0,-7.3436107119d0,&
                              -7.3421120744d0,-7.3383055793d0,-7.3187086340d0,-7.2765865573d0,&
                              -7.2362262431d0,-7.2022547993d0,-7.1751782638d0,-7.1507116654d0,&
                              -7.1346713199d0,-7.2679600206d0,-7.3233328016d0,-7.3357309216d0,&
                              -7.3418794512d0,-7.3435791075d0,-7.3420819400d0,-7.3382761798d0,&
                              -7.3186829095d0,-7.2765670801d0,-7.2362108084d0,-7.2022412020d0,&
                              -7.1751642990d0,-7.1506940257d0,-7.1346331006d0,-7.2679221688d0,&
                              -7.3232956848d0,-7.3356945397d0,-7.3418438043d0,-7.3435441956d0,&
                              -7.3420477631d0,-7.3382431054d0,-7.3186524076d0,-7.2765406206d0,&
                              -7.2361880238d0,-7.2022202549d0,-7.1751429844d0,-7.1506697712d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(2)%zrCut(5)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Wickoff Hollow -----------> 'c'
   this%wyckoffSite(3)%id='c'
   this%wyckoffSite(3)%myNumber=3
   this%wyckoffSite(3)%is_homonucl=.true.
   this%wyckoffSite(3)%x=2.7216780628885480d0
   this%wyckoffSite(3)%y=0.d0
   this%wyckoffSite(3)%n2dcuts=4
   allocate( this%wyckoffSite(3)%nPhiPoints(2) )
   this%wyckoffSite(3)%nPhiPoints(:)=[1,3]
   allocate( this%wyckoffSite(3)%zrCut(4) )
   this%wyckoffSite(3)%zrCut(:)%x=2.7216780628885480d0
   this%wyckoffSite(3)%zrCut(:)%y=0.d0
   ! Reading zrcuts
   this%wyckoffSite(3)%zrcut(1)%alias='Hollow_0_0'
   this%wyckoffSite(3)%zrCut(1)%theta=0.d0
   this%wyckoffSite(3)%zrCut(1)%phi=0.d0
   gridPot1414(:,:)=reshape( [ -6.9155832341d0,-7.0431757894d0,-7.0949357443d0,-7.1065338314d0,&
                               -7.1126007775d0,-7.1149446495d0,-7.1148035321d0,-7.1130384620d0,&
                               -7.1032602015d0,-7.0866983830d0,-7.0775919002d0,-7.0752053990d0,&
                               -7.0765959935d0,-7.0795432894d0,-6.9649287583d0,-7.0927587143d0,&
                               -7.1442276144d0,-7.1554471835d0,-7.1609845718d0,-7.1626478463d0,&
                               -7.1616850140d0,-7.1589677689d0,-7.1456549582d0,-7.1216381716d0,&
                               -7.1043222570d0,-7.0941989203d0,-7.0893777763d0,-7.0877354490d0,&
                               -7.0637296871d0,-7.1931130371d0,-7.2449909573d0,-7.2559113868d0,&
                               -7.2608001495d0,-7.2614708247d0,-7.2591846492d0,-7.2548324266d0,&
                               -7.2349852184d0,-7.1962117402d0,-7.1619216796d0,-7.1347734830d0,&
                               -7.1143239534d0,-7.0973924367d0,-7.1139211808d0,-7.2455660842d0,&
                               -7.2993351247d0,-7.3109945831d0,-7.3164617803d0,-7.3175466204d0,&
                               -7.3155074003d0,-7.3112382312d0,-7.2907258602d0,-7.2484402489d0,&
                               -7.2086991610d0,-7.1749118312d0,-7.1468114595d0,-7.1193957278d0,&
                               -7.1260984372d0,-7.2585944550d0,-7.3131826380d0,-7.3252162046d0,&
                               -7.3310288454d0,-7.3324242173d0,-7.3306598821d0,-7.3266292161d0,&
                               -7.3066232509d0,-7.2646639736d0,-7.2250246813d0,-7.1914078684d0,&
                               -7.1634706637d0,-7.1357473401d0,-7.1316064261d0,-7.2645588705d0,&
                               -7.3195983352d0,-7.3318457829d0,-7.3378583400d0,-7.3394396635d0,&
                               -7.3378440078d0,-7.3339640140d0,-7.3143060648d0,-7.2726654042d0,&
                               -7.2333421561d0,-7.2004088806d0,-7.1736645591d0,-7.1478499955d0,&
                               -7.1342582575d0,-7.2674587598d0,-7.3227451799d0,-7.3351083880d0,&
                               -7.3412289881d0,-7.3429087998d0,-7.3414009750d0,-7.3375966848d0,&
                               -7.3180864679d0,-7.2764487473d0,-7.2370380357d0,-7.2043381185d0,&
                               -7.1784996678d0,-7.1548683816d0,-7.1349043106d0,-7.2681768416d0,&
                               -7.3235360254d0,-7.3359301029d0,-7.3420778975d0,-7.3437812288d0,&
                               -7.3422895737d0,-7.3384952058d0,-7.3189710242d0,-7.2771392671d0,&
                               -7.2373522425d0,-7.2042657223d0,-7.1783008540d0,-7.1551998605d0,&
                               -7.1349465723d0,-7.2682411530d0,-7.3236165064d0,-7.3360186688d0,&
                               -7.3421719758d0,-7.3438771445d0,-7.3423869594d0,-7.3385903866d0,&
                               -7.3190386430d0,-7.2770812032d0,-7.2370453856d0,-7.2036148918d0,&
                               -7.1772990674d0,-7.1539459736d0,-7.1348936533d0,-7.2681988912d0,&
                               -7.3235804921d0,-7.3359852269d0,-7.3421392689d0,-7.3438448051d0,&
                               -7.3423538850d0,-7.3385562097d0,-7.3189916038d0,-7.2769742626d0,&
                               -7.2368226847d0,-7.2032099142d0,-7.1766677140d0,-7.1530525975d0,&
                               -7.1348230946d0,-7.2681268625d0,-7.3235117709d0,-7.3359172406d0,&
                               -7.3420727526d0,-7.3437782888d0,-7.3422866337d0,-7.3384882234d0,&
                               -7.3189151652d0,-7.2768588698d0,-7.2366253408d0,-7.2028777003d0,&
                               -7.1761510184d0,-7.1522830666d0,-7.1347323238d0,-7.2680423391d0,&
                               -7.3234283499d0,-7.3358341872d0,-7.3419896992d0,-7.3436952354d0,&
                               -7.3422032127d0,-7.3384040675d0,-7.3188258644d0,-7.2767453143d0,&
                               -7.2364584989d0,-7.2026171476d0,-7.1757548607d0,-7.1516767027d0,&
                               -7.1346893270d0,-7.2679989749d0,-7.3233846182d0,-7.3357900880d0,&
                               -7.3419456000d0,-7.3436507687d0,-7.3421587461d0,-7.3383592333d0,&
                               -7.3187795602d0,-7.2766901904d0,-7.2363850002d0,-7.2025102071d0,&
                               -7.1755961036d0,-7.1514326872d0,-7.1346470653d0,-7.2679563457d0,&
                               -7.3233405190d0,-7.3357463563d0,-7.3419007658d0,-7.3436059345d0,&
                               -7.3421139119d0,-7.3383132966d0,-7.3187321536d0,-7.2766369038d0,&
                               -7.2363170140d0,-7.2024153938d0,-7.1754593961d0,-7.1512235835d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(3)%zrCut(1)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! 2nd cut2d
   this%wyckoffSite(3)%zrcut(2)%alias='Hollow_90_0'
   this%wyckoffSite(3)%zrCut(2)%theta=1.57079632679d0
   this%wyckoffSite(3)%zrCut(2)%phi=0.d0
   gridPot1414(:,:)=reshape( [ -6.8925656617d0,-7.0062897566d0,-7.0400972985d0,-7.0410130916d0,&
                               -7.0351728889d0,-7.0243105233d0,-7.0095927859d0,-6.9918072149d0,&
                               -6.9244964155d0,-6.7748156782d0,-6.5731195782d0,-6.3010488876d0,&
                               -5.9367693621d0,-5.3591787448d0,-6.9498236831d0,-7.0685791281d0,&
                               -7.1083084563d0,-7.1125651307d0,-7.1103590687d0,-7.1034671002d0,&
                               -7.0931085678d0,-7.0801298086d0,-7.0305788557d0,-6.9265415155d0,&
                               -6.8019491751d0,-6.6550058950d0,-6.4870559656d0,-6.2760931557d0,&
                               -7.0595226244d0,-7.1864867662d0,-7.2353225772d0,-7.2444584595d0,&
                               -7.2473785609d0,-7.2458938882d0,-7.2412638407d0,-7.2343828970d0,&
                               -7.2058959224d0,-7.1498528337d0,-7.0959875100d0,-7.0482593236d0,&
                               -7.0073962787d0,-6.9684544886d0,-7.1126202547d0,-7.2435981578d0,&
                               -7.2965730454d0,-7.3077830596d0,-7.3127673706d0,-7.3133333102d0,&
                               -7.3107413803d0,-7.3058882644d0,-7.2834737483d0,-7.2379409667d0,&
                               -7.1959019434d0,-7.1618176791d0,-7.1358506057d0,-7.1139513152d0,&
                               -7.1253333163d0,-7.2574423637d0,-7.3115737525d0,-7.3233526463d0,&
                               -7.3288922396d0,-7.3299961893d0,-7.3279216899d0,-7.3235632199d0,&
                               -7.3024830719d0,-7.2586312044d0,-7.2176347594d0,-7.1841175372d0,&
                               -7.1584154265d0,-7.1365396556d0,-7.1311393422d0,-7.2638514460d0,&
                               -7.3186057359d0,-7.3306929565d0,-7.3365338943d0,-7.3379300012d0,&
                               -7.3361366341d0,-7.3320464342d0,-7.3116770181d0,-7.2686351057d0,&
                               -7.2279385353d0,-7.1942588810d0,-7.1680599195d0,-7.1453458964d0,&
                               -7.1340234293d0,-7.2671008214d0,-7.3222380392d0,-7.3345167238d0,&
                               -7.3405465532d0,-7.3421275091d0,-7.3405127438d0,-7.3365937957d0,&
                               -7.3166833787d0,-7.2741684516d0,-7.2336368857d0,-7.1997595201d0,&
                               -7.1730611352d0,-7.1495048176d0,-7.1347863453d0,-7.2680044873d0,&
                               -7.3232901724d0,-7.3356438256d0,-7.3417478886d0,-7.3434019757d0,&
                               -7.3418585041d0,-7.3380082772d0,-7.3182860168d0,-7.2759959956d0,&
                               -7.2355493206d0,-7.2015991913d0,-7.1746839854d0,-7.1507355525d0,&
                               -7.1348697662d0,-7.2681397248d0,-7.3234750215d0,-7.3358521943d0,&
                               -7.3419805118d0,-7.3436588535d0,-7.3421389014d0,-7.3383118267d0,&
                               -7.3186524076d0,-7.2764406624d0,-7.2360285318d0,-7.2020626003d0,&
                               -7.1750830830d0,-7.1510100699d0,-7.1348532290d0,-7.2681268625d0,&
                               -7.3234794314d0,-7.3358668941d0,-7.3420040314d0,-7.3436911929d0,&
                               -7.3421804282d0,-7.3383618058d0,-7.3187262737d0,-7.2765453980d0,&
                               -7.2361490696d0,-7.2021809331d0,-7.1751815712d0,-7.1510666639d0,&
                               -7.1347782604d0,-7.2680743110d0,-7.3234382722d0,-7.3358312472d0,&
                               -7.3419742644d0,-7.3436669384d0,-7.3421616860d0,-7.3383489435d0,&
                               -7.3187292136d0,-7.2765692851d0,-7.2361843489d0,-7.2022165800d0,&
                               -7.1752072958d0,-7.1510696039d0,-7.1347161541d0,-7.2680019148d0,&
                               -7.3233717559d0,-7.3357676709d0,-7.3419143630d0,-7.3436103444d0,&
                               -7.3421087670d0,-7.3382993319d0,-7.3186898919d0,-7.2765439280d0,&
                               -7.2361674442d0,-7.2022022477d0,-7.1751881861d0,-7.1510383669d0,&
                               -7.1346702174d0,-7.2679633281d0,-7.3233346391d0,-7.3357316565d0,&
                               -7.3418790837d0,-7.3435761675d0,-7.3420753251d0,-7.3382669925d0,&
                               -7.3186608599d0,-7.2765200410d0,-7.2361464971d0,-7.2021820356d0,&
                               -7.1751676065d0,-7.1510141124d0,-7.1346360405d0,-7.2679236388d0,&
                               -7.3232956848d0,-7.3356930697d0,-7.3418412319d0,-7.3435390507d0,&
                               -7.3420393107d0,-7.3382317131d0,-7.3186277855d0,-7.2764906415d0,&
                               -7.2361196701d0,-7.2021563111d0,-7.1751411470d0,-7.1509854479d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(3)%zrCut(2)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! third cut2d
   this%wyckoffSite(3)%zrcut(3)%alias='Hollow_90_45'
   this%wyckoffSite(3)%zrCut(3)%theta=1.57079632679d0
   this%wyckoffSite(3)%zrCut(3)%phi=0.785398163397d0
   gridPot1414(:,:)=reshape( [-6.9006196438d0,-7.0199344136d0,-7.0617540434d0,-7.0677937950d0,&
                              -7.0679308700d0,-7.0639924448d0,-7.0572294665d0,-7.0485217138d0,&
                              -7.0162065620d0,-6.9582091442d0,-6.9074149616d0,-6.8693941095d0,&
                              -6.8470615445d0,-6.8419978549d0,-6.9566623651d0,-7.0798291991d0,&
                              -7.1255971764d0,-7.1335912572d0,-7.1356613467d0,-7.1336261691d0,&
                              -7.1287407137d0,-7.1218876995d0,-7.0950257801d0,-7.0458228434d0,&
                              -7.0029690875d0,-6.9706870101d0,-6.9502775373d0,-6.9413033520d0,&
                              -7.0635521879d0,-7.1927352540d0,-7.2443246920d0,-7.2550529225d0,&
                              -7.2597149420d0,-7.2601269019d0,-7.2575489367d0,-7.2528725851d0,&
                              -7.2318997451d0,-7.1910605872d0,-7.1550694004d0,-7.1271204360d0,&
                              -7.1068954448d0,-7.0912597093d0,-7.1144422862d0,-7.2463558272d0,&
                              -7.3004295196d0,-7.3122484701d0,-7.3178788343d0,-7.3191297813d0,&
                              -7.3172588731d0,-7.3131587509d0,-7.2931858600d0,-7.2521210614d0,&
                              -7.2145650883d0,-7.1844163092d0,-7.1614744404d0,-7.1417330702d0,&
                              -7.1263681773d0,-7.2589986976d0,-7.3137335104d0,-7.3258413106d0,&
                              -7.3317259801d0,-7.3331922782d0,-7.3314955618d0,-7.3275270021d0,&
                              -7.3076871438d0,-7.2660402359d0,-7.2272009763d0,-7.1954870435d0,&
                              -7.1709803884d0,-7.1495140049d0,-7.1317052818d0,-7.2646988855d0,&
                              -7.3197762019d0,-7.3320379818d0,-7.3380600938d0,-7.3396458272d0,&
                              -7.3380483340d0,-7.3341576829d0,-7.3144141079d0,-7.2724522581d0,&
                              -7.2327949586d0,-7.1999818535d0,-7.1743253120d0,-7.1516046740d0,&
                              -7.1342685473d0,-7.2674672122d0,-7.3227426075d0,-7.3350955257d0,&
                              -7.3412014261d0,-7.3428617607d0,-7.3413289463d0,-7.3374919492d0,&
                              -7.3178351026d0,-7.2757380153d0,-7.2355901123d0,-7.2020188686d0,&
                              -7.1754821807d0,-7.1517399115d0,-7.1348962257d0,-7.2681614069d0,&
                              -7.3235069934d0,-7.3358926186d0,-7.3420290209d0,-7.3437172849d0,&
                              -7.3422094601d0,-7.3383941452d0,-7.3187791927d0,-7.2766615259d0,&
                              -7.2363614807d0,-7.2025102071d0,-7.1756056584d0,-7.1513841781d0,&
                              -7.1349359150d0,-7.2682235133d0,-7.3235907819d0,-7.3359863294d0,&
                              -7.3421322865d0,-7.3438297379d0,-7.3423296304d0,-7.3385216653d0,&
                              -7.3189225151d0,-7.2768074207d0,-7.2364739336d0,-7.2025506313d0,&
                              -7.1755409796d0,-7.1511677246d0,-7.1348829960d0,-7.2681827215d0,&
                              -7.3235577075d0,-7.3359569299d0,-7.3421061945d0,-7.3438069533d0,&
                              -7.3423101533d0,-7.3385047606d0,-7.3189122253d0,-7.2768004383d0,&
                              -7.2364599689d0,-7.2025164544d0,-7.1754748308d0,-7.1510523317d0,&
                              -7.1348164797d0,-7.2681117953d0,-7.3234911912d0,-7.3358926186d0,&
                              -7.3420440882d0,-7.3437466844d0,-7.3422517218d0,-7.3384485342d0,&
                              -7.3188607762d0,-7.2767530317d0,-7.2364103573d0,-7.2024572880d0,&
                              -7.1753991272d0,-7.1509501685d0,-7.1347304863d0,-7.2680272719d0,&
                              -7.3234085053d0,-7.3358106676d0,-7.3419636071d0,-7.3436669384d0,&
                              -7.3421730783d0,-7.3383709931d0,-7.3187861751d0,-7.2766817380d0,&
                              -7.2363394311d0,-7.2023830544d0,-7.1753171762d0,-7.1508542528d0,&
                              -7.1346882246d0,-7.2679842752d0,-7.3233655085d0,-7.3357680384d0,&
                              -7.3419206104d0,-7.3436243092d0,-7.3421308166d0,-7.3383287314d0,&
                              -7.3187446484d0,-7.2766416812d0,-7.2362990068d0,-7.2023422626d0,&
                              -7.1752741795d0,-7.1508075812d0,-7.1346444929d0,-7.2679412785d0,&
                              -7.3233217768d0,-7.3357243067d0,-7.3418772462d0,-7.3435805775d0,&
                              -7.3420870849d0,-7.3382853672d0,-7.3187016517d0,-7.2765990520d0,&
                              -7.2362567451d0,-7.2022996334d0,-7.1752300803d0,-7.1507605420d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(3)%zrCut(3)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! third cut2d
   this%wyckoffSite(3)%zrcut(4)%alias='Hollow_90_90'
   this%wyckoffSite(3)%zrCut(4)%theta=1.57079632679d0
   this%wyckoffSite(3)%zrCut(4)%phi=1.57079632679d0
   gridPot1414(:,:)=reshape( [-6.9064025177d0,-7.0281412729d0,-7.0723304993d0,-7.0794484761d0,&
                              -7.0805373586d0,-7.0773614819d0,-7.0711074817d0,-7.0625684085d0,&
                              -7.0274220886d0,-6.9381322527d0,-6.7570407645d0,-6.3143201714d0,&
                              -5.2611778488d0,-2.3857877032d0,-6.9621012652d0,-7.0876817950d0,&
                              -7.1359285143d0,-7.1451173156d0,-7.1483089945d0,-7.1472829534d0,&
                              -7.1432526548d0,-7.1370596585d0,-7.1104990836d0,-7.0518978743d0,&
                              -6.9716663796d0,-6.8353256474d0,-6.5833230284d0,-6.0721726190d0,&
                              -7.0672947392d0,-7.1982932220d0,-7.2518888056d0,-7.2636651270d0,&
                              -7.2693903044d0,-7.2708613799d0,-7.2693241556d0,-7.2656547354d0,&
                              -7.2473230694d0,-7.2088042640d0,-7.1716738481d0,-7.1390496345d0,&
                              -7.1113230035d0,-7.0848671641d0,-7.1162290384d0,-7.2490282381d0,&
                              -7.3041066571d0,-7.3164647202d0,-7.3226551441d0,-7.3244797481d0,&
                              -7.3231931542d0,-7.3196843286d0,-7.3014621756d0,-7.2628959636d0,&
                              -7.2267963662d0,-7.1967898071d0,-7.1728149147d0,-7.1509769956d0,&
                              -7.1273927485d0,-7.2605318795d0,-7.3158432891d0,-7.3282616212d0,&
                              -7.3344693172d0,-7.3362681967d0,-7.3349114116d0,-7.3312908681d0,&
                              -7.3125064504d0,-7.2725187744d0,-7.2349481016d0,-7.2039033740d0,&
                              -7.1794658076d0,-7.1576418532d0,-7.1322664440d0,-7.2655375051d0,&
                              -7.3209279258d0,-7.3333572826d0,-7.3395539539d0,-7.3413186565d0,&
                              -7.3399041749d0,-7.3362013129d0,-7.3170288224d0,-7.2759996705d0,&
                              -7.2371486512d0,-7.2049139804d0,-7.1795793630d0,-7.1571009031d0,&
                              -7.1345107253d0,-7.2678277230d0,-7.3232346809d0,-7.3356577904d0,&
                              -7.3418364545d0,-7.3435702876d0,-7.3421120744d0,-7.3383515159d0,&
                              -7.3189225151d0,-7.2771925536d0,-7.2373739246d0,-7.2040823432d0,&
                              -7.1777859959d0,-7.1544273897d0,-7.1349884666d0,-7.2683124466d0,&
                              -7.3237127897d0,-7.3361252418d0,-7.3422906761d0,-7.3440079721d0,&
                              -7.3425284443d0,-7.3387425288d0,-7.3192117323d0,-7.2772204831d0,&
                              -7.2370262760d0,-7.2032712856d0,-7.1764880098d0,-7.1526086656d0,&
                              -7.1349899365d0,-7.2683014218d0,-7.3236958850d0,-7.3361039272d0,&
                              -7.3422642166d0,-7.3439748977d0,-7.3424880200d0,-7.3386940197d0,&
                              -7.3191316188d0,-7.2770665034d0,-7.2367697657d0,-7.2028839477d0,&
                              -7.1759503671d0,-7.1518916863d0,-7.1349252577d0,-7.2682323331d0,&
                              -7.3236231213d0,-7.3360304286d0,-7.3421885130d0,-7.3438966217d0,&
                              -7.3424068040d0,-7.3386094962d0,-7.3190364380d0,-7.2769478031d0,&
                              -7.2366205634d0,-7.2026961587d0,-7.1757170089d0,-7.1516013666d0,&
                              -7.1348396318d0,-7.2681430323d0,-7.3235316155d0,-7.3359374528d0,&
                              -7.3420940672d0,-7.3438003384d0,-7.3423086833d0,-7.3385095380d0,&
                              -7.3189302325d0,-7.2768287353d0,-7.2364864284d0,-7.2025429140d0,&
                              -7.1755413471d0,-7.1513963054d0,-7.1347455535d0,-7.2680467490d0,&
                              -7.3234327598d0,-7.3358374946d0,-7.3419926391d0,-7.3436974403d0,&
                              -7.3422046827d0,-7.3384040675d0,-7.3188207195d0,-7.2767122399d0,&
                              -7.2363625831d0,-7.2024102489d0,-7.1753987597d0,-7.1512408557d0,&
                              -7.1346959419d0,-7.2679993424d0,-7.3233838832d0,-7.3357878830d0,&
                              -7.3419422925d0,-7.3436467262d0,-7.3421532337d0,-7.3383518834d0,&
                              -7.3187666980d0,-7.2766556460d0,-7.2363037842d0,-7.2023492450d0,&
                              -7.1753355509d0,-7.1511743394d0,-7.1346533127d0,-7.2679526707d0,&
                              -7.3233353741d0,-7.3357386389d0,-7.3418923134d0,-7.3435960122d0,&
                              -7.3421021521d0,-7.3383004344d0,-7.3187134114d0,-7.2766005220d0,&
                              -7.2362460878d0,-7.2022908136d0,-7.1752752820d0,-7.1511126006d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(3)%zrCut(4)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Wickoff Bridge -----------> 'f'
   this%wyckoffSite(4)%id='f'
   this%wyckoffSite(4)%myNumber=4
   this%wyckoffSite(4)%is_homonucl=.true.
   this%wyckoffSite(4)%x=1.36083903144d0
   this%wyckoffSite(4)%y=1.36083903144d0
   this%wyckoffSite(4)%n2dcuts=7
   allocate( this%wyckoffSite(4)%nPhiPoints(3) )
   this%wyckoffSite(4)%nPhiPoints(:)=[1,3,3]
   allocate( this%wyckoffSite(4)%zrCut(7) )
   this%wyckoffSite(4)%zrCut(:)%x=1.36083903144d0
   this%wyckoffSite(4)%zrCut(:)%y=1.36083903144d0
   ! First zrcuts
   this%wyckoffSite(4)%zrcut(1)%alias='Bridge_0_0'
   this%wyckoffSite(4)%zrCut(1)%theta=0.d0
   this%wyckoffSite(4)%zrCut(1)%phi=0.d0
   gridPot1414(:,:)=reshape( [-6.6746098227d0,-6.8113661897d0,-6.8747940554d0,-6.8930176784d0,&
                              -6.9061375550d0,-6.9158812712d0,-6.9234034906d0,-6.9294829315d0,&
                              -6.9436067997d0,-6.9648464398d0,-6.9878140332d0,-7.0109010619d0,&
                              -7.0315031012d0,-7.0497109220d0,-6.8048167250d0,-6.9354557993d0,&
                              -6.9912331904d0,-7.0050615941d0,-7.0134922568d0,-7.0183001711d0,&
                              -7.0206954921d0,-7.0215127971d0,-7.0195106939d0,-7.0150743153d0,&
                              -7.0165020266d0,-7.0234755786d0,-7.0339050371d0,-7.0478271516d0,&
                              -7.0202012137d0,-7.1482479906d0,-7.1989161230d0,-7.2092731854d0,&
                              -7.2136151182d0,-7.2137437409d0,-7.2109063254d0,-7.2059793434d0,&
                              -7.1841307670d0,-7.1406596224d0,-7.1002507992d0,-7.0671889011d0,&
                              -7.0430262197d0,-7.0263350436d0,-7.1061104792d0,-7.2372419945d0,&
                              -7.2905292513d0,-7.3019446943d0,-7.3071700808d0,-7.3080101704d0,&
                              -7.3057206874d0,-7.3011883931d0,-7.2797308295d0,-7.2348826878d0,&
                              -7.1904781104d0,-7.1494077994d0,-7.1114538311d0,-7.0698249302d0,&
                              -7.1236255751d0,-7.2559290265d0,-7.3103279504d0,-7.3222714811d0,&
                              -7.3279973935d0,-7.3293089769d0,-7.3274630583d0,-7.3233511763d0,&
                              -7.3030813509d0,-7.2603488678d0,-7.2190055093d0,-7.1820992642d0,&
                              -7.1488407573d0,-7.1112994839d0,-7.1308677647d0,-7.2637540603d0,&
                              -7.3187314186d0,-7.3309502018d0,-7.3369370344d0,-7.3384952058d0,&
                              -7.3368797055d0,-7.3329828070d0,-7.3132870061d0,-7.2715566771d0,&
                              -7.2318762255d0,-7.1979650505d0,-7.1692509651d0,-7.1391918544d0,&
                              -7.1341200800d0,-7.2673055151d0,-7.3225798080d0,-7.3349378711d0,&
                              -7.3410544288d0,-7.3427313005d0,-7.3412220058d0,-7.3374173481d0,&
                              -7.3179177885d0,-7.2763373968d0,-7.2369954065d0,-7.2042741747d0,&
                              -7.1781832561d0,-7.1537067354d0,-7.1348774836d0,-7.2681456047d0,&
                              -7.3235033185d0,-7.3358966610d0,-7.3420440882d0,-7.3437466844d0,&
                              -7.3422568668d0,-7.3384636014d0,-7.3189467697d0,-7.2771462495d0,&
                              -7.2374176563d0,-7.2044072072d0,-7.1785037102d0,-7.1553869146d0,&
                              -7.1349395899d0,-7.2682290257d0,-7.3236040117d0,-7.3360061740d0,&
                              -7.3421602160d0,-7.3438657522d0,-7.3423763021d0,-7.3385804642d0,&
                              -7.3190327631d0,-7.2770900230d0,-7.2370836049d0,-7.2037008852d0,&
                              -7.1774464322d0,-7.1541605896d0,-7.1348899784d0,-7.2681904389d0,&
                              -7.3235724072d0,-7.3359771420d0,-7.3421319190d0,-7.3438378227d0,&
                              -7.3423472701d0,-7.3385503298d0,-7.3189886639d0,-7.2769801425d0,&
                              -7.2368462043d0,-7.2032624658d0,-7.1767639972d0,-7.1532095171d0,&
                              -7.1348190522d0,-7.2681198802d0,-7.3235055235d0,-7.3359106258d0,&
                              -7.3420665052d0,-7.3437724089d0,-7.3422811213d0,-7.3384834460d0,&
                              -7.3189129603d0,-7.2768621772d0,-7.2366297507d0,-7.2029085698d0,&
                              -7.1762094499d0,-7.1523848622d0,-7.1347271788d0,-7.2680368267d0,&
                              -7.3234224700d0,-7.3358283073d0,-7.3419838193d0,-7.3436897230d0,&
                              -7.3421980678d0,-7.3383992901d0,-7.3188232919d0,-7.2767471518d0,&
                              -7.2364665837d0,-7.2026355223d0,-7.1757894051d0,-7.1517384416d0,&
                              -7.1346896945d0,-7.2679934625d0,-7.3233791058d0,-7.3357845756d0,&
                              -7.3419400876d0,-7.3436456238d0,-7.3421539686d0,-7.3383544559d0,&
                              -7.3187762528d0,-7.2766909253d0,-7.2363912476d0,-7.2025241718d0,&
                              -7.1756218282d0,-7.1514800938d0,-7.1346481678d0,-7.2679515683d0,&
                              -7.3233353741d0,-7.3357412114d0,-7.3418956209d0,-7.3436011571d0,&
                              -7.3421091345d0,-7.3383088867d0,-7.3187288462d0,-7.2766365363d0,&
                              -7.2363217914d0,-7.2024264186d0,-7.1754792408d0,-7.1512599654d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(4)%zrCut(1)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Second cuts2d
   this%wyckoffSite(4)%zrcut(2)%alias='Bridge_90_0'
   this%wyckoffSite(4)%zrCut(2)%theta=1.57079632679d0
   this%wyckoffSite(4)%zrCut(2)%phi=0.d0
   gridPot1414(:,:)=reshape( [-6.6458520056d0,-6.7681765451d0,-6.8148717079d0,-6.8238907273d0,&
                              -6.8272709303d0,-6.8267626871d0,-6.8235563085d0,-6.8184977638d0,&
                              -6.7975223514d0,-6.7619155625d0,-6.7443820919d0,-6.7518649895d0,&
                              -6.7767203957d0,-6.8098524850d0,-6.7978108336d0,-6.9250487579d0,&
                              -6.9767300692d0,-6.9882950819d0,-6.9942833845d0,-6.9964604145d0,&
                              -6.9960293449d0,-6.9938210779d0,-6.9816030297d0,-6.9576579043d0,&
                              -6.9401854375d0,-6.9320693490d0,-6.9309404097d0,-6.9331009026d0,&
                              -7.0248044342d0,-7.1553163559d0,-7.2088935649d0,-7.2208580428d0,&
                              -7.2269036743d0,-7.2288289714d0,-7.2278775314d0,-7.2249239881d0,&
                              -7.2094822891d0,-7.1778910990d0,-7.1495481818d0,-7.1266504121d0,&
                              -7.1087391584d0,-7.0929325386d0,-7.1077520715d0,-7.2397670406d0,&
                              -7.2940575540d0,-7.3060282793d0,-7.3118405526d0,-7.3132976634d0,&
                              -7.3116571735d0,-7.3078109891d0,-7.2886936225d0,-7.2491697230d0,&
                              -7.2129532629d0,-7.1837290968d0,-7.1613101709d0,-7.1419517287d0,&
                              -7.1242609710d0,-7.2568903888d0,-7.3116538660d0,-7.3237888608d0,&
                              -7.3297095446d0,-7.3312188394d0,-7.3295721021d0,-7.3256619739d0,&
                              -7.3060264418d0,-7.2647716491d0,-7.2263006178d0,-7.1948744322d0,&
                              -7.1706051777d0,-7.1495489168d0,-7.1310765009d0,-7.2640598147d0,&
                              -7.3191327212d0,-7.3313948686d0,-7.3374206555d0,-7.3390104313d0,&
                              -7.3374188181d0,-7.3335369868d0,-7.3138250162d0,-7.2719285803d0,&
                              -7.2323374295d0,-7.1995959856d0,-7.1740515295d0,-7.1516961799d0,&
                              -7.1341498469d0,-7.2673441019d0,-7.3226147198d0,-7.3349650656d0,&
                              -7.3410694961d0,-7.3427268906d0,-7.3411911363d0,-7.3373526693d0,&
                              -7.3176866353d0,-7.2755730109d0,-7.2354082032d0,-7.2018362244d0,&
                              -7.1753594380d0,-7.1519361529d0,-7.1348532290d0,-7.2681342124d0,&
                              -7.3234768590d0,-7.3358610142d0,-7.3419959465d0,-7.3436820056d0,&
                              -7.3421712408d0,-7.3383537209d0,-7.3187281112d0,-7.2765869247d0,&
                              -7.2362626250d0,-7.2024036340d0,-7.1755475945d0,-7.1516311335d0,&
                              -7.1349134979d0,-7.2682106510d0,-7.3235753472d0,-7.3359697922d0,&
                              -7.3421146469d0,-7.3438098932d0,-7.3423072133d0,-7.3384970433d0,&
                              -7.3188879707d0,-7.2767504592d0,-7.2363938201d0,-7.2024631679d0,&
                              -7.1755016579d0,-7.1514312172d0,-7.1348587414d0,-7.2681724317d0,&
                              -7.3235448452d0,-7.3359429652d0,-7.3420914948d0,-7.3437896811d0,&
                              -7.3422899411d0,-7.3384830785d0,-7.3188809884d0,-7.2767482543d0,&
                              -7.2363850002d0,-7.2024352384d0,-7.1754421239d0,-7.1513228067d0,&
                              -7.1348072924d0,-7.2681029755d0,-7.3234797989d0,-7.3358793888d0,&
                              -7.3420301234d0,-7.3437301472d0,-7.3422322447d0,-7.3384272196d0,&
                              -7.3188299068d0,-7.2767012152d0,-7.2363379611d0,-7.2023797470d0,&
                              -7.1753715652d0,-7.1512276260d0,-7.1347161541d0,-7.2680191870d0,&
                              -7.3233974805d0,-7.3357981728d0,-7.3419496424d0,-7.3436504012d0,&
                              -7.3421539686d0,-7.3383496785d0,-7.3187556732d0,-7.2766310239d0,&
                              -7.2362685049d0,-7.2023080858d0,-7.1752929217d0,-7.1511364876d0,&
                              -7.1346812422d0,-7.2679765578d0,-7.3233544837d0,-7.3357551761d0,&
                              -7.3419066457d0,-7.3436077720d0,-7.3421113394d0,-7.3383074168d0,&
                              -7.3187137789d0,-7.2765909672d0,-7.2362284481d0,-7.2022680290d0,&
                              -7.1752517624d0,-7.1510923884d0,-7.1346319981d0,-7.2679335611d0,&
                              -7.3233107520d0,-7.3357114444d0,-7.3418629140d0,-7.3435640403d0,&
                              -7.3420676077d0,-7.3382636851d0,-7.3186707822d0,-7.2765487055d0,&
                              -7.2361869214d0,-7.2022265023d0,-7.1752091332d0,-7.1510475543d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(4)%zrCut(2)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! third cut2d
   this%wyckoffSite(4)%zrcut(3)%alias='Bridge_90_54'
   this%wyckoffSite(4)%zrCut(3)%theta=1.57079632679d0
   this%wyckoffSite(4)%zrCut(3)%phi=0.942477796077d0
   gridPot1414(:,:)=reshape( [-6.6059602453d0,-6.6999852318d0,-6.7044785719d0,-6.6853839899d0,&
                              -6.6547934839d0,-6.6134240333d0,-6.5613171649d0,-6.4980344590d0,&
                              -6.2319002544d0,-5.4521493932d0,-4.1165565330d0,-2.3944015246d0,&
                              -1.3930856731d0,-2.5257612532d0,-6.7789617371d0,-6.8945967969d0,&
                              -6.9304090145d0,-6.9319638784d0,-6.9262486234d0,-6.9147740140d0,&
                              -6.8984782606d0,-6.8779552324d0,-6.7935721664d0,-6.5820665690d0,&
                              -6.3018904471d0,-6.0311177426d0,-5.8992990149d0,-6.0217481346d0,&
                              -7.0222896778d0,-7.1517156570d0,-7.2041786264d0,-7.2155856170d0,&
                              -7.2210711888d0,-7.2224312814d0,-7.2209069194d0,-7.2173642844d0,&
                              -7.2000171328d0,-7.1646561969d0,-7.1321749382d0,-7.1059245276d0,&
                              -7.0866524463d0,-7.0723352767d0,-7.1074573419d0,-7.2393683105d0,&
                              -7.2935853252d0,-7.3055398808d0,-7.3113495816d0,-7.3128217596d0,&
                              -7.3112143441d0,-7.3074199763d0,-7.2885734522d0,-7.2497860092d0,&
                              -7.2143835466d0,-7.1857572921d0,-7.1635621695d0,-7.1438546088d0,&
                              -7.1241768150d0,-7.2567805083d0,-7.3115259784d0,-7.3236598706d0,&
                              -7.3295849644d0,-7.3311049165d0,-7.3294761864d0,-7.3255917827d0,&
                              -7.3060859757d0,-7.2652097011d0,-7.2272142061d0,-7.1961889556d0,&
                              -7.1720979353d0,-7.1507866340d0,-7.1310562887d0,-7.2640300478d0,&
                              -7.3190989119d0,-7.3313625292d0,-7.3373912561d0,-7.3389861768d0,&
                              -7.3374022809d0,-7.3335307394d0,-7.3138694829d0,-7.2721244542d0,&
                              -7.2327438771d0,-7.2001990420d0,-7.1747291870d0,-7.1521529740d0,&
                              -7.1341483770d0,-7.2673407945d0,-7.3226121474d0,-7.3349646981d0,&
                              -7.3410709660d0,-7.3427313005d0,-7.3411992212d0,-7.3373647966d0,&
                              -7.3177178722d0,-7.2756542269d0,-7.2355552005d0,-7.2020409182d0,&
                              -7.1755622942d0,-7.1518813965d0,-7.1348561690d0,-7.2681356824d0,&
                              -7.3234801664d0,-7.3358661591d0,-7.3420032964d0,-7.3436919279d0,&
                              -7.3421837356d0,-7.3383695231d0,-7.3187564081d0,-7.2766435187d0,&
                              -7.2363497209d0,-7.2025076346d0,-7.1756163158d0,-7.1514231324d0,&
                              -7.1349134979d0,-7.2682128560d0,-7.3235797571d0,-7.3359760396d0,&
                              -7.3421230992d0,-7.3438205506d0,-7.3423204431d0,-7.3385135805d0,&
                              -7.3189155327d0,-7.2768022758d0,-7.2364695237d0,-7.2025476914d0,&
                              -7.1755428171d0,-7.1511890392d0,-7.1348804235d0,-7.2681750042d0,&
                              -7.3235492552d0,-7.3359492125d0,-7.3420999471d0,-7.3438003384d0,&
                              -7.3423035384d0,-7.3384992482d0,-7.3189078154d0,-7.2767982334d0,&
                              -7.2364573964d0,-7.2025142495d0,-7.1754755658d0,-7.1510696039d0,&
                              -7.1348072924d0,-7.2681055479d0,-7.3234842089d0,-7.3358860037d0,&
                              -7.3420385758d0,-7.3437408045d0,-7.3422458420d0,-7.3384433893d0,&
                              -7.3188563663d0,-7.2767500918d0,-7.2364074173d0,-7.2024547156d0,&
                              -7.1753994947d0,-7.1509663382d0,-7.1347220339d0,-7.2680221270d0,&
                              -7.3234022579d0,-7.3358047877d0,-7.3419580947d0,-7.3436614260d0,&
                              -7.3421675659d0,-7.3383662157d0,-7.3187821327d0,-7.2766791656d0,&
                              -7.2363364911d0,-7.2023797470d0,-7.1753164413d0,-7.1508696875d0,&
                              -7.1346794047d0,-7.2679794978d0,-7.3233592612d0,-7.3357621585d0,&
                              -7.3419154655d0,-7.3436187968d0,-7.3421249367d0,-7.3383235865d0,&
                              -7.3187406059d0,-7.2766387413d0,-7.2362956994d0,-7.2023389552d0,&
                              -7.1752734445d0,-7.1508226484d0,-7.1346367755d0,-7.2679361336d0,&
                              -7.3233158970d0,-7.3357184268d0,-7.3418717338d0,-7.3435750651d0,&
                              -7.3420812050d0,-7.3382802223d0,-7.3186972417d0,-7.2765961121d0,&
                              -7.2362538051d0,-7.2022963260d0,-7.1752293454d0,-7.1507756092d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(4)%zrCut(3)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Fourh cut2d
   this%wyckoffSite(4)%zrcut(4)%alias='Bridge_90_144'
   this%wyckoffSite(4)%zrCut(4)%theta=1.57079632679d0
   this%wyckoffSite(4)%zrCut(4)%phi=2.51327412287d0
   gridPot1414(:,:)=reshape( [-6.6769162103d0,-6.8143105457d0,-6.8780765051d0,-6.8963177678d0,&
                              -6.9093439337d0,-6.9188741362d0,-6.9260608343d0,-6.9316820111d0,&
                              -6.9436865458d0,-6.9589536855d0,-6.9728078137d0,-6.9831847207d0,&
                              -6.9870246577d0,-6.9796001915d0,-6.8134903007d0,-6.9478138625d0,&
                              -7.0070769271d0,-7.0225675027d0,-7.0325515595d0,-7.0387871850d0,&
                              -7.0424702024d0,-7.0444271040d0,-7.0449041102d0,-7.0413250934d0,&
                              -7.0393380574d0,-7.0382565247d0,-7.0354992229d0,-7.0268263821d0,&
                              -7.0271527161d0,-7.1585418442d0,-7.2128746193d0,-7.2251544064d0,&
                              -7.2314664705d0,-7.2336100587d0,-7.2328320755d0,-7.2300126672d0,&
                              -7.2148267435d0,-7.1838044329d0,-7.1570226271d0,-7.1367344270d0,&
                              -7.1216907232d0,-7.1080640733d0,-7.1080434937d0,-7.2401624634d0,&
                              -7.2945224330d0,-7.3065137379d0,-7.3123322586d0,-7.3137820195d0,&
                              -7.3121191125d0,-7.3082350763d0,-7.2889218358d0,-7.2489157852d0,&
                              -7.2123439590d0,-7.1831547049d0,-7.1611253218d0,-7.1422740203d0,&
                              -7.1243451269d0,-7.2570050467d0,-7.3117887361d0,-7.3239299782d0,&
                              -7.3298524995d0,-7.3313581193d0,-7.3297021947d0,-7.3257769993d0,&
                              -7.3060609862d0,-7.2645687929d0,-7.2258250816d0,-7.1942228667d0,&
                              -7.1699179654d0,-7.1488264250d0,-7.1311011229d0,-7.2640950941d0,&
                              -7.3191742480d0,-7.3314400703d0,-7.3374673272d0,-7.3390574705d0,&
                              -7.3374643872d0,-7.3335792485d0,-7.3138452283d0,-7.2718694138d0,&
                              -7.2321569903d0,-7.1992861887d0,-7.1736274423d0,-7.1510288121d0,&
                              -7.1341549919d0,-7.2673532892d0,-7.3226264796d0,-7.3349793978d0,&
                              -7.3410852983d0,-7.3427448978d0,-7.3412109810d0,-7.3373739839d0,&
                              -7.3177123598d0,-7.2755961629d0,-7.2354129806d0,-7.2017983726d0,&
                              -7.1752352252d0,-7.1515135357d0,-7.1348734412d0,-7.2681386223d0,&
                              -7.3234845764d0,-7.3358698340d0,-7.3420069713d0,-7.3436952354d0,&
                              -7.3421870430d0,-7.3383720956d0,-7.3187567756d0,-7.2766321264d0,&
                              -7.2363192189d0,-7.2024521431d0,-7.1755369372d0,-7.1513268491d0,&
                              -7.1349274627d0,-7.2682139585d0,-7.3235815946d0,-7.3359778770d0,&
                              -7.3421249367d0,-7.3438223880d0,-7.3423226480d0,-7.3385157854d0,&
                              -7.3189170027d0,-7.2767997033d0,-7.2364599689d0,-7.2025289492d0,&
                              -7.1755148876d0,-7.1511537598d0,-7.1348800560d0,-7.2681757392d0,&
                              -7.3235507251d0,-7.3359506825d0,-7.3421014171d0,-7.3438021759d0,&
                              -7.3423053759d0,-7.3385010857d0,-7.3189100203d0,-7.2767982334d0,&
                              -7.2364533540d0,-7.2025046947d0,-7.1754604986d0,-7.1510508617d0,&
                              -7.1348083949d0,-7.2681059154d0,-7.3234849438d0,-7.3358871062d0,&
                              -7.3420400457d0,-7.3437422745d0,-7.3422476794d0,-7.3384452267d0,&
                              -7.3188589388d0,-7.2767519292d0,-7.2364066823d0,-7.2024499382d0,&
                              -7.1753910424d0,-7.1509556809d0,-7.1347249739d0,-7.2680224945d0,&
                              -7.3234029929d0,-7.3358055227d0,-7.3419591972d0,-7.3436628959d0,&
                              -7.3421690359d0,-7.3383676856d0,-7.3187847051d0,-7.2766817380d0,&
                              -7.2363372261d0,-7.2023779095d0,-7.1753116638d0,-7.1508630726d0,&
                              -7.1346816097d0,-7.2679794978d0,-7.3233599961d0,-7.3357628935d0,&
                              -7.3419165680d0,-7.3436198992d0,-7.3421264067d0,-7.3383254239d0,&
                              -7.3187428109d0,-7.2766413138d0,-7.2362971693d0,-7.2023378527d0,&
                              -7.1752697696d0,-7.1508171360d0,-7.1346356730d0,-7.2679365010d0,&
                              -7.3233162644d0,-7.3357191618d0,-7.3418724688d0,-7.3435761675d0,&
                              -7.3420826750d0,-7.3382816922d0,-7.3186994467d0,-7.2765986845d0,&
                              -7.2362552751d0,-7.2022955910d0,-7.1752264054d0,-7.1507711993d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(4)%zrCut(4)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Fifth cut2d
   this%wyckoffSite(4)%zrcut(5)%alias='Bridge_135_54'
   this%wyckoffSite(4)%zrCut(5)%theta=2.35619449019d0
   this%wyckoffSite(4)%zrCut(5)%phi=0.942477796077d0
   gridPot1414(:,:)=reshape( [-6.6524169050d0,-6.7781863263d0,-6.8274892213d0,-6.8369484976d0,&
                              -6.8398568392d0,-6.8377022263d0,-6.8314265440d0,-6.8216280714d0,&
                              -6.7745113938d0,-6.6526102065d0,-6.5226774542d0,-6.4649644761d0,&
                              -6.5300482664d0,-6.7171592941d0,-6.8002509888d0,-6.9286968634d0,&
                              -6.9801139471d0,-6.9902630083d0,-6.9934620370d0,-6.9910229843d0,&
                              -6.9836290200d0,-6.9715513542d0,-6.9059141191d0,-6.6704281169d0,&
                              -6.2204734192d0,-5.5837439500d0,-5.0443605039d0,-5.1022491438d0,&
                              -7.0260785333d0,-7.1577601860d0,-7.2129686976d0,-7.2258658733d0,&
                              -7.2328570650d0,-7.2356783107d0,-7.2354901542d0,-7.2330624937d0,&
                              -7.2163794025d0,-7.1587663826d0,-7.0204808761d0,-6.6657627901d0,&
                              -5.7925348748d0,-3.1559945441d0,-7.1090416054d0,-7.2418514624d0,&
                              -7.2972157911d0,-7.3098443292d0,-7.3163971015d0,-7.3186825420d0,&
                              -7.3179538029d0,-7.3150983803d0,-7.2993468844d0,-7.2658902986d0,&
                              -7.2340874324d0,-7.2030827616d0,-7.1616262151d0,-7.0569255495d0,&
                              -7.1252370331d0,-7.2584173233d0,-7.3138867550d0,-7.3264414271d0,&
                              -7.3328310323d0,-7.3348618000d0,-7.3337901897d0,-7.3305099449d0,&
                              -7.3130830473d0,-7.2763318844d0,-7.2426015161d0,-7.2150950135d0,&
                              -7.1926900523d0,-7.1682153691d0,-7.1317181441d0,-7.2650527815d0,&
                              -7.3205622700d0,-7.3330790902d0,-7.3393856420d0,-7.3412852146d0,&
                              -7.3400331651d0,-7.3365217670d0,-7.3181232173d0,-7.2790572144d0,&
                              -7.2427992275d0,-7.2132965016d0,-7.1902866465d0,-7.1692652973d0,&
                              -7.1344710360d0,-7.2678383803d0,-7.3233210419d0,-7.3357941304d0,&
                              -7.3420312259d0,-7.3438345153d0,-7.3424571506d0,-7.3387899354d0,&
                              -7.3197291628d0,-7.2789653411d0,-7.2405703809d0,-7.2090174101d0,&
                              -7.1844909103d0,-7.1627716915d0,-7.1350016963d0,-7.2683624257d0,&
                              -7.3238017230d0,-7.3362417372d0,-7.3424369384d0,-7.3441876763d0,&
                              -7.3427474702d0,-7.3390056539d0,-7.3196435369d0,-7.2780888697d0,&
                              -7.2385649702d0,-7.2057195257d0,-7.1799865455d0,-7.1572320982d0,&
                              -7.1350141911d0,-7.2683385387d0,-7.3237572563d0,-7.3361811008d0,&
                              -7.3423590299d0,-7.3440895556d0,-7.3426247275d0,-7.3388557167d0,&
                              -7.3193859241d0,-7.2775519620d0,-7.2376069153d0,-7.2042058209d0,&
                              -7.1778561871d0,-7.1544046051d0,-7.1349351800d0,-7.2682591601d0,&
                              -7.3236683230d0,-7.3360870225d0,-7.3422576017d0,-7.3439800426d0,&
                              -7.3425063947d0,-7.3387263591d0,-7.3192172447d0,-7.2772807520d0,&
                              -7.2371773157d0,-7.2035535204d0,-7.1769301041d0,-7.1531242586d0,&
                              -7.1348436742d0,-7.2681628769d0,-7.3235650574d0,-7.3359793470d0,&
                              -7.3421451488d0,-7.3438624448d0,-7.3423829169d0,-7.3385962665d0,&
                              -7.3190625300d0,-7.2770639310d0,-7.2368631090d0,-7.2030989313d0,&
                              -7.1762906659d0,-7.1522272076d0,-7.1347481260d0,-7.2680621837d0,&
                              -7.3234588518d0,-7.3358705690d0,-7.3420330634d0,-7.3437463169d0,&
                              -7.3422627466d0,-7.3384716862d0,-7.3189225151d0,-7.2768864318d0,&
                              -7.2366279133d0,-7.2027784771d0,-7.1758537164d0,-7.1516120239d0,&
                              -7.1347029243d0,-7.2680133071d0,-7.3234074028d0,-7.3358176500d0,&
                              -7.3419786744d0,-7.3436904579d0,-7.3422050502d0,-7.3384125198d0,&
                              -7.3188578363d0,-7.2768085232d0,-7.2365308951d0,-7.2026527944d0,&
                              -7.1756883444d0,-7.1513841781d0,-7.1346599276d0,-7.2679662680d0,&
                              -7.3233570562d0,-7.3357662009d0,-7.3419261228d0,-7.3436364364d0,&
                              -7.3421495587d0,-7.3383555584d0,-7.3187960974d0,-7.2767368620d0,&
                              -7.2364441666d0,-7.2025447514d0,-7.1755494320d0,-7.1511963890d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(4)%zrCut(5)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Sixth cut2d
   this%wyckoffSite(4)%zrcut(6)%alias='Bridge_135_144'
   this%wyckoffSite(4)%zrCut(6)%theta=2.35619449019d0
   this%wyckoffSite(4)%zrCut(6)%phi=2.51327412287d0
   gridPot1414(:,:)=reshape( [-6.6740618902d0,-6.8098815170d0,-6.8718555793d0,-6.8891685540d0,&
                              -6.9012737818d0,-6.9099135482d0,-6.9162675066d0,-6.9211371597d0,&
                              -6.9316298271d0,-6.9481552637d0,-6.9695275689d0,-6.9941073552d0,&
                              -7.0183938819d0,-7.0433815857d0,-6.8076059988d0,-6.9389466178d0,&
                              -6.9949077555d0,-7.0086891200d0,-7.0169621281d0,-7.0215116946d0,&
                              -7.0235615720d0,-7.0239614047d0,-7.0204643389d0,-7.0137564845d0,&
                              -7.0144201773d0,-7.0223764063d0,-7.0347561515d0,-7.0511419407d0,&
                              -7.0229416109d0,-7.1522345575d0,-7.2042286055d0,-7.2152754527d0,&
                              -7.2203273825d0,-7.2211954016d0,-7.2191344994d0,-7.2150329072d0,&
                              -7.1960820151d0,-7.1593966335d0,-7.1283408811d0,-7.1059399623d0,&
                              -7.0915511314d0,-7.0825082249d0,-7.1067510199d0,-7.2382243040d0,&
                              -7.2918672942d0,-7.3034811836d0,-7.3089149388d0,-7.3099766268d0,&
                              -7.3079241770d0,-7.3036513329d0,-7.2831933510d0,-7.2413299895d0,&
                              -7.2027843570d0,-7.1713611114d0,-7.1469713191d0,-7.1254119598d0,&
                              -7.1237880072d0,-7.2561756144d0,-7.3106642067d0,-7.3226533066d0,&
                              -7.3284229507d0,-7.3297793683d0,-7.3279786514d0,-7.3239130735d0,&
                              -7.3038251572d0,-7.2618813147d0,-7.2227462231d0,-7.1905483016d0,&
                              -7.1652666032d0,-7.1424765091d0,-7.1308699697d0,-7.2637555303d0,&
                              -7.3187233338d0,-7.3309336646d0,-7.3369080024d0,-7.3384496366d0,&
                              -7.3368128217d0,-7.3328894637d0,-7.3130900297d0,-7.2712071910d0,&
                              -7.2318049318d0,-7.1992924361d0,-7.1738240512d0,-7.1510413069d0,&
                              -7.1340689984d0,-7.2672525961d0,-7.3224985920d0,-7.3348393829d0,&
                              -7.3409349935d0,-7.3425876107d0,-7.3410500189d0,-7.3372126544d0,&
                              -7.3175818997d0,-7.2756571668d0,-7.2358852094d0,-7.2028762303d0,&
                              -7.1770105852d0,-7.1540742287d0,-7.1348491866d0,-7.2681092229d0,&
                              -7.3234496645d0,-7.3358323497d0,-7.3419687520d0,-7.3436577510d0,&
                              -7.3421517637d0,-7.3383408586d0,-7.3187512632d0,-7.2767453143d0,&
                              -7.2366738499d0,-7.2031996244d0,-7.1768150787d0,-7.1533521045d0,&
                              -7.1349138654d0,-7.2682051386d0,-7.3235709373d0,-7.3359668522d0,&
                              -7.3421146469d0,-7.3438132007d0,-7.3423160332d0,-7.3385110080d0,&
                              -7.3189280275d0,-7.2768750395d0,-7.2366698075d0,-7.2029622238d0,&
                              -7.1762715562d0,-7.1524080143d0,-7.1348782186d0,-7.2681724317d0,&
                              -7.3235481527d0,-7.3359481101d0,-7.3420992121d0,-7.3438014409d0,&
                              -7.3423064783d0,-7.3385043931d0,-7.3189239851d0,-7.2768526224d0,&
                              -7.2365918989d0,-7.2027850920d0,-7.1759507346d0,-7.1518773540d0,&
                              -7.1348069249d0,-7.2681062829d0,-7.3234871488d0,-7.3358896786d0,&
                              -7.3420426182d0,-7.3437470519d0,-7.3422539268d0,-7.3384533116d0,&
                              -7.3188743735d0,-7.2767923535d0,-7.2364970857d0,-7.2026248650d0,&
                              -7.1756923869d0,-7.1514620866d0,-7.1347187265d0,-7.2680254344d0,&
                              -7.3234081378d0,-7.3358121376d0,-7.3419665471d0,-7.3436713483d0,&
                              -7.3421793257d0,-7.3383794454d0,-7.3188019773d0,-7.2767144449d0,&
                              -7.2363982300d0,-7.2024852175d0,-7.1754884281d0,-7.1511504524d0,&
                              -7.1346768323d0,-7.2679835402d0,-7.3233662435d0,-7.3357702433d0,&
                              -7.3419250203d0,-7.3436301890d0,-7.3421381664d0,-7.3383382862d0,&
                              -7.3187611856d0,-7.2766714482d0,-7.2363489859d0,-7.2024223762d0,&
                              -7.1754031697d0,-7.1510255047d0,-7.1346312631d0,-7.2679420134d0,&
                              -7.3233239818d0,-7.3357279816d0,-7.3418831261d0,-7.3435879273d0,&
                              -7.3420955372d0,-7.3382960245d0,-7.3187185563d0,-7.2766273490d0,&
                              -7.2363001093d0,-7.2023632098d0,-7.1753259961d0,-7.1509170942d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(4)%zrCut(6)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Sixth cut2d, THE LAST ONE (thank goodness (>_<)Uu )
   this%wyckoffSite(4)%zrcut(7)%alias='Bridge_135_234'
   this%wyckoffSite(4)%zrCut(7)%theta=2.35619449019d0
   this%wyckoffSite(4)%zrCut(7)%phi=4.08407044967d0
   gridPot1414(:,:)=reshape( [-6.6335189320d0,-6.7466719423d0,-6.7809414232d0,-6.7825547186d0,&
                              -6.7777637091d0,-6.7683389771d0,-6.7554906779d0,-6.7400603712d0,&
                              -6.6836281071d0,-6.5747909375d0,-6.4733613294d0,-6.4107511363d0,&
                              -6.4123273148d0,-6.4961907454d0,-6.7815885789d0,-6.8965412037d0,&
                              -6.9309933288d0,-6.9318922173d0,-6.9257462601d0,-6.9142606259d0,&
                              -6.8985885086d0,-6.8795170787d0,-6.8064597873d0,-6.6435941270d0,&
                              -6.4368394824d0,-6.1999474510d0,-5.9866477515d0,-5.9159890885d0,&
                              -7.0153969744d0,-7.1397228821d0,-7.1849440295d0,-7.1918088035d0,&
                              -7.1920825860d0,-7.1875241997d0,-7.1793342450d0,-7.1683454617d0,&
                              -7.1235502390d0,-7.0221548078d0,-6.8890998328d0,-6.7113896500d0,&
                              -6.4628594747d0,-5.9963936726d0,-7.1044207452d0,-7.2345192370d0,&
                              -7.2863780475d0,-7.2969210615d0,-7.3011509088d0,-7.3008638966d0,&
                              -7.2973072969d0,-7.2913623585d0,-7.2647613593d0,-7.2080847122d0,&
                              -7.1471863026d0,-7.0828007496d0,-7.0102976380d0,-6.9030274596d0,&
                              -7.1225451450d0,-7.2542451724d0,-7.3078775054d0,-7.3193675494d0,&
                              -7.3245870561d0,-7.3253378448d0,-7.3228727001d0,-7.3180798531d0,&
                              -7.2953919220d0,-7.2474564695d0,-7.1997073360d0,-7.1553549427d0,&
                              -7.1129807655d0,-7.0605369057d0,-7.1301992945d0,-7.2627265492d0,&
                              -7.3172614456d0,-7.3292255560d0,-7.3349330937d0,-7.3361855107d0,&
                              -7.3342366940d0,-7.3299774472d0,-7.3090130595d0,-7.2645845951d0,&
                              -7.2217330442d0,-7.1846504024d0,-7.1531088239d0,-7.1201814284d0,&
                              -7.1337981559d0,-7.2668093992d0,-7.3218749559d0,-7.3341146862d0,&
                              -7.3401026213d0,-7.3416402131d0,-7.3399798786d0,-7.3360127889d0,&
                              -7.3159465547d0,-7.2731299157d0,-7.2322268141d0,-7.1978283430d0,&
                              -7.1702781087d0,-7.1447968615d0,-7.1347216664d0,-7.2679206988d0,&
                              -7.3231843344d0,-7.3355269628d0,-7.3416185310d0,-7.3432604908d0,&
                              -7.3417048919d0,-7.3378421703d0,-7.3180813230d0,-7.2757380153d0,&
                              -7.2352597359d0,-7.2013033592d0,-7.1743491990d0,-7.1501343335d0,&
                              -7.1348539640d0,-7.2681062829d0,-7.3234323923d0,-7.3358073601d0,&
                              -7.3419316352d0,-7.3436063020d0,-7.3420830424d0,-7.3382522928d0,&
                              -7.3185829513d0,-7.2763645913d0,-7.2359660579d0,-7.2020346708d0,&
                              -7.1750838180d0,-7.1508696875d0,-7.1348388968d0,-7.2681081204d0,&
                              -7.3234581168d0,-7.3358444770d0,-7.3419808793d0,-7.3436676734d0,&
                              -7.3421572761d0,-7.3383382862d0,-7.3187020191d0,-7.2765270233d0,&
                              -7.2361483346d0,-7.2022066576d0,-7.1752179531d0,-7.1509325289d0,&
                              -7.1347789954d0,-7.2680643887d0,-7.3234287174d0,-7.3358224274d0,&
                              -7.3419658121d0,-7.3436599560d0,-7.3421565411d0,-7.3383449011d0,&
                              -7.3187310511d0,-7.2765839848d0,-7.2362163208d0,-7.2022632516d0,&
                              -7.1752392677d0,-7.1508836523d0,-7.1347099067d0,-7.2679986074d0,&
                              -7.3233702860d0,-7.3357684059d0,-7.3419165680d0,-7.3436147543d0,&
                              -7.3421161168d0,-7.3383088867d0,-7.3187086340d0,-7.2765799424d0,&
                              -7.2362199958d0,-7.2022588417d0,-7.1752083982d0,-7.1507969239d0,&
                              -7.1346632350d0,-7.2679614906d0,-7.3233357416d0,-7.3357349640d0,&
                              -7.3418845961d0,-7.3435842524d0,-7.3420878199d0,-7.3382813247d0,&
                              -7.3186858494d0,-7.2765634052d0,-7.2362063985d0,-7.2022426720d0,&
                              -7.1751826737d0,-7.1507498847d0,-7.1346400829d0,-7.2679243738d0,&
                              -7.3232993598d0,-7.3356993171d0,-7.3418500517d0,-7.3435508105d0,&
                              -7.3420540105d0,-7.3382500878d0,-7.3186579200d0,-7.2765402531d0,&
                              -7.2361854514d0,-7.2022202549d0,-7.1751529068d0,-7.1507024781d0 ],shape( gridPot1414  ) )
   call this%wyckoffSite(4)%zrCut(7)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   call this%interpol()
   return
end subroutine initialize_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT1D_PHI_PES_H2LiF001 #######################################
!#######################################################################
SUBROUTINE PLOT1D_PHI_PES_H2LiF001(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
   INTEGER,INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin 
   REAL(KIND=8),DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=*),PARAMETER :: routinename = "PLOT1D_PHI_PES_H2LiF001: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_PHI_PES_H2LiF001 ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:5)=x(1:5)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(6)=xmin
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(6),v,dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(6)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(6),v,dvdu(:)
   END DO
   ! Final value
   r(6) = xmax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(6),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_PHI_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT1D_ATOMIC_INTERAC_PHI_PES_H2LiF001 #######################################
!#######################################################################
SUBROUTINE PLOT1D_ATOMIC_INTERAC_PHI_PES_H2LiF001(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
   INTEGER,INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta
   REAL(KIND=8),DIMENSION(2) :: v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8),DIMENSION(6) :: r, dvdu, ratom
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=*),PARAMETER :: routinename = "PLOT1D_PHI_PES_H2LiF001: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_PHI_PES_H2LiF001 ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:5)=x(1:5)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(6)=xmin
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
   WRITE(11,*) r(6),sum(v),v(1),v(2),ratom(:),dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(6)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
      WRITE(11,*) r(6),sum(v),v(1),v(2),ratom(:),dvdu(:)  
   END DO
   ! Final value
   r(6) = xmax
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
   WRITE(11,*) r(6),sum(v),v(1),v(2),ratom(:),dvdu(:)  
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_ATOMIC_INTERAC_PHI_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT1D_THETA_PES_H2LiF001 #########################################
!#######################################################################
SUBROUTINE PLOT1D_THETA_PES_H2LiF001(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin 
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=*),PARAMETER :: routinename = "PLOT1D_THETA_PES_H2LiF001: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_THETA_PES_H2LiF001 ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:4)=x(1:4)
   r(6)=x(6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(5)=xmin
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(5),v,dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(5)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(5),v,dvdu(:)
   END DO
   ! Final value
   r(5) = xmax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(5),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_THETA_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT1D_ATOMIC_INTERAC_THETA_PES_H2LiF001 ########################
!#######################################################################
SUBROUTINE PLOT1D_ATOMIC_INTERAC_THETA_PES_H2LiF001(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8), DIMENSION(6) :: r,dvdu,ratom
   REAL(KIND=8),DIMENSION(2) :: v
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=*),PARAMETER :: routinename = "PLOT1D_THETA_PES_H2LiF001: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_THETA_PES_H2LiF001 ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:4)=x(1:4)
   r(6)=x(6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(5)=xmin
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
   WRITE(11,*) r(5),sum(v),v(:),ratom(:),dvdu(:)  

   ! cycle for inpoints
   DO i=1, inpoints
      r(5)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
      WRITE(11,*) r(5),sum(v),v(:),ratom(:),dvdu(:)
   END DO
   ! Final value
   r(5) = xmax
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
   WRITE(11,*) r(5),sum(v),v(:),ratom(:),dvdu(:)  
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_ATOMIC_INTERAC_THETA_PES_H2LiF001
!#######################################################################
!# SUBROUTINE: PLOT1D_R_PES_H2LiF001 ###################################
!#######################################################################
SUBROUTINE PLOT1D_R_PES_H2LiF001(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=*),PARAMETER :: routinename = "PLOT1D_R_PES_H2LiF001: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_R_PES_H2LiF001 ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = thispes%wyckoffsite(1)%zrcut(1)%getfirstr()
   xmax = thispes%wyckoffsite(1)%zrcut(1)%getlastr()
   r(1:3)=x(1:3)
   r(5:6)=x(5:6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=(xmax-xmin)/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(4)=xmin
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),v,dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(4)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(4),v,dvdu(:)
   END DO
   ! Final value
   r(4) = xmax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_R_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT1D_Z_PES_H2LiF001 #########################################
!#######################################################################
SUBROUTINE PLOT1D_Z_PES_H2LiF001(thispes,npoints,X,L,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   REAL(KIND=8),INTENT(IN) :: L
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=*),PARAMETER :: routinename = "PLOT1D_Z_PES_H2LiF001: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_Z_PES_H2LiF001 ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = thispes%wyckoffsite(1)%zrcut(1)%getfirstz()
   xmax = xmin+L
   r(1:2)=x(1:2)
   r(4:6)=x(4:6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=(xmax-xmin)/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(3)=xmin
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(3)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(3),v,dvdu(:)
   END DO
   ! Final value
   r(3) = xmax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_Z_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT_XYMAP_PES_H2LiF001
!#######################################################################
SUBROUTINE PLOT_XYMAP_PES_H2LiF001(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: init_point ! Initial position to start the scan (in a.u. and radians)
   INTEGER,INTENT(IN) :: nxpoints, nypoints ! number of points in XY plane
   CHARACTER(LEN=*),INTENT(IN) :: filename ! filename
   REAL(KIND=8),INTENT(IN) :: Lx ! Length of X axis 
   REAL(KIND=8),INTENT(IN) :: Ly ! Length of X axis 
   ! Local variables
   REAL(KIND=8) :: xmin, ymin, xmax, ymax
   REAL(KIND=8),DIMENSION(6) :: r,dvdu
   REAL(KIND=8) :: xdelta, ydelta
   INTEGER :: xinpoints, nxdelta
   INTEGER :: yinpoints, nydelta
   INTEGER :: i, j ! counters
   REAL(KIND=8) :: v ! potential
   ! GABBA, GABBA HEY! ---------
   xmin = init_point(1)
   ymin = init_point(2)
   xmax = init_point(1)+Lx
   ymax = init_point(2)+Ly
   ! For X, grid parameters
   xinpoints=nxpoints-2
   nxdelta=nxpoints-1
   xdelta=Lx/DFLOAT(nxdelta)
   ! For Y, grid parameters
   yinpoints=nypoints-2
   nydelta=nypoints-1
   ydelta=(ymax-ymin)/DFLOAT(nydelta)
   ! Let's go! 
   ! 1st XY point
   OPEN(11,file=filename,status="replace")
   r(1) = xmin
   r(2) = ymin
   r(3:6)=init_point(3:6)
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
   END DO
   r(2) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   ! inpoints in XY
   DO i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      r(3:6) = init_point(3:6)
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
      DO j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
         WRITE(11,*) r(1:2),v,dvdu(:)
      END DO
      r(2) = ymax
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
   END DO
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   r(3:6) = init_point(3:6)
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
   END DO
   r(2) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_XYMAP_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT_RZMAP_PES_H2LiF001
!#######################################################################
SUBROUTINE PLOT_RZMAP_PES_H2LiF001(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: init_point 
   INTEGER,INTENT(IN) :: nxpoints, nypoints 
   CHARACTER(LEN=*),INTENT(IN) :: filename 
   REAL(KIND=8),INTENT(IN) :: Lx
   REAL(KIND=8),INTENT(IN) :: Ly
   ! Local variables
   REAL(KIND=8) :: xmin, ymin, xmax, ymax
   REAL(KIND=8),DIMENSION(6) :: r,dvdu
   REAL(KIND=8) :: xdelta, ydelta
   INTEGER :: xinpoints, nxdelta
   INTEGER :: yinpoints, nydelta
   INTEGER :: i, j ! counters
   REAL(KIND=8) :: v ! potential
   ! GABBA, GABBA HEY! ---------
   xmin = init_point(4)
   ymin = init_point(3)
   xmax = init_point(4)+Lx
   ymax = init_point(3)+Ly
   ! For X, grid parameters
   xinpoints=nxpoints-2
   nxdelta=nxpoints-1
   xdelta=Lx/DFLOAT(nxdelta)
   ! For Y, grid parameters
   yinpoints=nypoints-2
   nydelta=nypoints-1
   ydelta=(ymax-ymin)/DFLOAT(nydelta)
   ! Let's go! 
   ! 1st XY point
   OPEN(11,file=filename,status="replace")
   r(4) = xmin
   r(3) = ymin
   r(1:2)=init_point(1:2)
   r(5:6)=init_point(5:6)
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),r(3),v,dvdu(:)
   DO i =1, yinpoints
      r(3) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(4),r(3),v,dvdu(:)
   END DO
   r(3) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),r(3),v,dvdu(:)
   ! inpoints in XY
   DO i = 1, xinpoints
      r(4) = xmin+DFLOAT(i)*xdelta
      r(3) = ymin
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(4),r(3),v,dvdu(:)
      DO j = 1, yinpoints
         r(3) = ymin + DFLOAT(j)*ydelta
         CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
         WRITE(11,*) r(4),r(3),v,dvdu(:)
      END DO
      r(3) = ymax
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(4),r(3),v,dvdu(:)
   END DO
   ! Last point in XY plane
   r(4) = xmax
   r(3) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),r(3),v,dvdu(:)
   DO i =1, yinpoints
      r(3) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(4),r(3),v,dvdu(:)
   END DO
   r(3) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),r(3),v,dvdu(:)
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_RZMAP_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT_ATOMIC_INTERAC_RZ_PES_H2LiF001
!#######################################################################
SUBROUTINE PLOT_ATOMIC_INTERAC_RZ_PES_H2LiF001(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: init_point 
   INTEGER,INTENT(IN) :: nxpoints, nypoints 
   CHARACTER(LEN=*),INTENT(IN) :: filename 
   REAL(KIND=8),INTENT(IN) :: Lx
   REAL(KIND=8),INTENT(IN) :: Ly
   ! Local variables
   REAL(KIND=8) :: xmin, ymin, xmax, ymax
   REAL(KIND=8),DIMENSION(6) :: r
   REAL(KIND=8) :: xdelta, ydelta
   INTEGER :: xinpoints, nxdelta
   INTEGER :: yinpoints, nydelta
   REAL(KIND=8),DIMENSION(2) :: v
   REAL(KIND=8),DIMENSION(6) :: atomicx
   REAL(KIND=8),DIMENSION(6) :: dvdu
   INTEGER :: i, j ! counters
   ! GABBA, GABBA HEY! ---------
   xmin = init_point(4)
   ymin = init_point(3)
   xmax = init_point(4)+Lx
   ymax = init_point(3)+Ly
   ! For X, grid parameters
   xinpoints=nxpoints-2
   nxdelta=nxpoints-1
   xdelta=Lx/DFLOAT(nxdelta)
   ! For Y, grid parameters
   yinpoints=nypoints-2
   nydelta=nypoints-1
   ydelta=(ymax-ymin)/DFLOAT(nydelta)
   ! Let's go! 
   ! 1st XY point
   OPEN(11,file=filename,status="replace")
   r(4) = xmin
   r(3) = ymin
   r(1:2)=init_point(1:2)
   r(5:6)=init_point(5:6)
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
   WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   DO i =1, yinpoints
      r(3) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
      WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   END DO
   r(3) = ymax
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
   WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   ! inpoints in XY
   DO i = 1, xinpoints
      r(4) = xmin+DFLOAT(i)*xdelta
      r(3) = ymin
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
      WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
      DO j = 1, yinpoints
         r(3) = ymin + DFLOAT(j)*ydelta
         CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
         WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
      END DO
      r(3) = ymax
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
      WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   END DO
   ! Last point in XY plane
   r(4) = xmax
   r(3) = ymax
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
   WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   DO i =1, yinpoints
      r(3) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
      WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   END DO
   r(3) = ymax
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
   WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_ATOMIC_INTERAC_RZ_PES_H2LiF001
!###########################################################
!# FUNCTION: is_allowed_PES_H2LiF001
!###########################################################
LOGICAL FUNCTION is_allowed_PES_H2LiF001(this,x)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8) :: zmin, rmin, rmax
   ! Run section
   zmin=this%wyckoffsite(1)%zrcut(1)%getfirstZ()
   rmin=this%wyckoffsite(1)%zrcut(1)%getfirstR()
   rmax=this%wyckoffsite(1)%zrcut(1)%getlastR()
   SELECT CASE(size(x)/=6)
      CASE(.TRUE.)
         WRITE(0,*) "is_allowed_PES_H2LiF001 ERR: wrong number of dimensions"
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(x(3)<zmin .OR. x(4)<rmin .OR. x(4)>rmax)
      CASE(.TRUE.)
         is_allowed_PES_H2LiF001=.FALSE.
      CASE(.FALSE.)
         is_allowed_PES_H2LiF001=.TRUE.
   END SELECT
   RETURN
END FUNCTION is_allowed_PES_H2LiF001
!###########################################################
! FUNCTION: from_molecular_to_atomic
!###########################################################
pure function from_molecular_to_atomic(molcoord) result(atomcoord)
   ! Initial declarations
   implicit none
   ! i/o variables
   real(kind=8),dimension(6),intent(in):: molcoord
   ! dymmy function variable
   real(kind=8),dimension(6):: atomcoord
   real(kind=8),dimension(2),parameter:: masa=[ 1.d0,1.d0 ]
   real(kind=8):: mTot
   ! run section
   mTot=sum(masa(:))
   atomcoord(1)=molcoord(1)+(masa(2)/(mTot))*molcoord(4)*dcos(molcoord(6))*dsin(molcoord(5))
   atomcoord(2)=molcoord(2)+(masa(2)/(mTot))*molcoord(4)*dsin(molcoord(6))*dsin(molcoord(5))
   atomcoord(3)=molcoord(3)+(masa(2)/(mTot))*molcoord(4)*dcos(molcoord(5))
   atomcoord(4)=molcoord(1)-(masa(1)/(mTot))*molcoord(4)*dcos(molcoord(6))*dsin(molcoord(5))
   atomcoord(5)=molcoord(2)-(masa(1)/(mTot))*molcoord(4)*dsin(molcoord(6))*dsin(molcoord(5))
   atomcoord(6)=molcoord(3)-(masa(1)/(mTot))*molcoord(4)*dcos(molcoord(5))
   return
end function from_molecular_to_atomic

end module PES_H2LiF001_MOD
