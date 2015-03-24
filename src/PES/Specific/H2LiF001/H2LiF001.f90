!#######################################################
! MODULE: PES_H2LiF001_MOD
!#######################################################
!> @brief
!! CRP6D PES explicit implementation
!#######################################################
module PES_H2LiF001_MOD
! Massive name spaces
use SYSTEM_MOD
use LOGISTIC_MOD, only: Logistic_func
! Selective name spaces
use UNITS_MOD, only: Length, pi
use PES_MOD, only: PES
use HLIF001_NS_MOD, only: HLi001_NS
use EXTRAPOL_TO_VACUUM_MOD, only: Vacuumpot
use FOURIER_P4MM_MOD, only: Fourierp4mm
use WYCKOFF_P4MM_MOD, only: WyckoffSitio, Wyckoffp4mm
implicit none
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
   type(Wyckoffp4mm),dimension(4):: wyckoffsite
   type(HLiF001_NS),dimension(1):: atomiccrp
   type(VacuumPot):: farpot
   type(Logistic_func):: dumpfunc
   character(len=30):: extrapol2vac_flag='Xexponential'
   integer(kind=4),dimension(4,2):: xyklist
   contains
      ! Initialization block
      procedure,public:: READ => READ_PES_H2LiF001
      procedure,public:: INITIALIZE => INITIALIZE_PES_H2LiF001
      ! Set block
      procedure,public:: SET_SMOOTH => SET_SMOOTH_PES_H2LiF001
      ! Get block
      procedure,public:: GET_V_AND_DERIVS => GET_V_AND_DERIVS_PES_H2LiF001
      procedure,public:: GET_V_AND_DERIVS_PURE => GET_V_AND_DERIVS_PURE_PES_H2LiF001
      procedure,public:: GET_ATOMICPOT_AND_DERIVS => GET_ATOMICPOT_AND_DERIVS_PES_H2LiF001
      ! Tools block
      procedure,public:: SMOOTH => SMOOTH_PES_H2LiF001
      procedure,public:: SMOOTH_EXTRA => SMOOTH_EXTRA_PES_H2LiF001
      procedure,public:: INTERPOL => INTERPOL_PES_H2LiF001
      procedure,public:: RAWINTERPOL => RAWINTERPOL_PES_H2LiF001
      procedure,public:: EXTRACT_VACUUMSURF => EXTRACT_VACUUMSURF_PES_H2LiF001
      procedure,public:: ADD_VACUUMSURF => ADD_VACUUMSURF_PES_H2LiF001
      procedure,public:: INTERPOL_NEW_RZGRID => INTERPOL_NEW_RZGRID_PES_H2LiF001
      ! Plot toolk
      procedure,public:: PLOT1D_THETA => PLOT1D_THETA_PES_H2LiF001
      procedure,public:: PLOT1D_THETA_SMOOTH => PLOT1D_THETA_SMOOTH_PES_H2LiF001
      procedure,public:: PLOT1D_ATOMIC_INTERAC_THETA => PLOT1D_ATOMIC_INTERAC_THETA_PES_H2LiF001
      procedure,public:: PLOT1D_PHI => PLOT1D_PHI_PES_H2LiF001
      procedure,public:: PLOT1D_PHI_SMOOTH => PLOT1D_PHI_SMOOTH_PES_H2LiF001
      procedure,public:: PLOT1D_ATOMIC_INTERAC_PHI => PLOT1D_ATOMIC_INTERAC_PHI_PES_H2LiF001
      procedure,public:: PLOT1D_R => PLOT1D_R_PES_H2LiF001
      procedure,public:: PLOT1D_R_SMOOTH => PLOT1D_R_SMOOTH_PES_H2LiF001
      procedure,public:: PLOT1D_Z => PLOT1D_Z_PES_H2LiF001
      procedure,public:: PLOT1D_Z_SMOOTH => PLOT1D_Z_SMOOTH_PES_H2LiF001
      procedure,public:: PLOT_XYMAP => PLOT_XYMAP_PES_H2LiF001
      procedure,public:: PLOT_RZMAP => PLOT_RZMAP_PES_H2LiF001
      procedure,public:: PLOT_ATOMIC_INTERAC_RZ => PLOT_ATOMIC_INTERAC_RZ_PES_H2LiF001
      ! Enquire block
      procedure,public:: is_allowed => is_allowed_PES_H2LiF001
end type PES_H2LiF001
contains
!###########################################################
!# SUBROUTINE: INITIALIZE_PES_H2LiF001
!###########################################################
SUBROUTINE INITIALIZE_PES_H2LiF001(this,filename,tablename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),INTENT(OUT)::this
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: filename,tablename
   ! Local variables
   CHARACTER(LEN=:),ALLOCATABLE:: auxstring
   ! Run section
   SELECT CASE(allocated(system_inputfile) .or. .not.present(filename))
      CASE(.TRUE.)
         auxstring=system_inputfile
      CASE(.FALSE.)
         auxstring=filename
   END SELECT
   SELECT CASE(present(tablename))
      CASE(.TRUE.)
         CALL this%READ(filename=auxstring,tablename=tablename)
      CASE(.FALSE.)
         CALL this%READ(filename=auxstring,tablename='pes')
   END SELECT
   CALL this%INTERPOL()
   SELECT CASE(this%is_resized)
      CASE(.TRUE.)
         CALL this%INTERPOL_NEW_RZGRID(this%grid(1),this%grid(2))
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE INITIALIZE_PES_H2LiF001
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
   CALL xyinterpol%INTERPOL(system_surface)
   ALLOCATE(aux1(5))
   ALLOCATE(aux2(5,2))
   CALL xyinterpol%GET_F_AND_DERIVS(system_surface,x(1:2),aux1,aux2)
   !--------------------------------------
   ! Results for the real potential
   !-------------------------------------
   CALL this%GET_ATOMICPOT_AND_DERIVS(x,atomicx,atomic_v,atomic_dvdu)
   dvdu_atomicA=atomic_dvdu(1:3)
   dvdu_atomicB=atomic_dvdu(4:6)
   v=aux1(1)+sum(atomic_v)

   ma=system_mass(1)
   mb=system_mass(2)
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
SUBROUTINE GET_V_AND_DERIVS_PES_H2LiF001(this,X,v,dvdu,errCode)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),TARGET,INTENT(IN):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   REAL(KIND=8),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: dvdu
   integer(kind=1),optional,intent(out):: errCode
   ! Local variables
   REAL(KIND=8) :: zcrp, zvac ! last PES_H2LiF001 z value and Z infinity
   REAL(KIND=8) :: vzcrp, vzvac ! potentials at zcrp and zvac
   REAL(KIND=8),DIMENSION(6) :: dvducrp ! derivatives at zcrp
   REAL(KIND=8),DIMENSION(6) :: dvduvac ! derivatives at vacuum
   REAL(KIND=8) :: alpha,beta,gama ! parameters
   CLASS(Function1d),ALLOCATABLE:: extrapolfunc
   INTEGER(KIND=4) :: i !counter
   ! Local Parameter
   REAL(KIND=8),PARAMETER :: zero=0.D0 ! what we will consider zero (a.u.)
   REAL(KIND=8),PARAMETER :: dz=0.5D0 ! 0.25 Angstroems approx
   character(len=*),parameter:: routinename='GET_V_AND_DERIVS_PES_H2LiF001: '
   ! Run section
   zcrp=this%wyckoffsite(1)%zrcut(1)%getlastZ()
   zvac=this%zvacuum
   ! Check if we are in the pure PES_H2LiF001 region
   SELECT CASE(x(3)<= zcrp) !easy
      CASE(.TRUE.)
         call this%get_v_and_derivs_pure(x,v,dvdu)
         RETURN
      CASE(.FALSE.)
         ! do nothing, next switch
   END SELECT
   ! Check if we are in the extrapolation region
   SELECT CASE(x(3)>zcrp .AND. x(3)<zvac)
      CASE(.TRUE.) ! uff
         ! Set potential and derivs
         vzvac=this%farpot%getpot(x(4))
         CALL this%GET_V_AND_DERIVS_PURE([x(1),x(2),zcrp,x(4),x(5),x(6)],vzcrp,dvducrp)
         dvduvac(1:3)=zero
         dvduvac(4)=this%farpot%getderiv(x(4))
         dvduvac(5:6)=zero
         ! Check kind of extrapolation
         SELECT CASE(this%extrapol2vac_flag)
            CASE("Xexponential")
               ALLOCATE(Xexponential_func::extrapolfunc)
               ! Extrapol potential
               beta=-1.D0/zvac
               alpha=(vzcrp-vzvac)/(zcrp*dexp(beta*zcrp)-zvac*dexp(beta*zvac))
               gama=vzvac-alpha*zvac*dexp(beta*zvac)
               CALL extrapolfunc%READ([alpha,beta])
               v=extrapolfunc%getvalue(x(3))+gama
               dvdu(3)=extrapolfunc%getderiv(x(3))
               ! Extrapol derivatives
               DO i = 1, 6
                  SELECT CASE(i)
                     CASE(3)
                        ! Skip dvdz
                     CASE DEFAULT
                        beta=-1.D0/zvac
                        alpha=(dvducrp(i)-dvduvac(i))/(zcrp*dexp(beta*zcrp)-zvac*dexp(beta*zvac))
                        gama=dvduvac(i)-alpha*zvac*dexp(beta*zvac)
                        CALL extrapolfunc%READ([alpha,beta])
                        dvdu(i)=extrapolfunc%getvalue(x(3))+gama
                  END SELECT
               END DO
               DEALLOCATE(extrapolfunc)
               RETURN

            CASE("Linear")
               ALLOCATE(Linear_func::extrapolfunc)
               ! Extrapol potential
               beta=(vzcrp*zvac-vzvac*zcrp)/(zvac-zcrp)
               alpha=(vzvac-beta)/zvac
               CALL extrapolfunc%READ([alpha,beta])
               v=extrapolfunc%getvalue(x(3))
               dvdu(3)=extrapolfunc%getderiv(x(3))
               ! Extrapol derivatives
               DO i = 1, 6
                  SELECT CASE(i)
                     CASE(3)
                        ! do nothing
                     CASE DEFAULT
                        beta=(dvducrp(i)*zvac-dvduvac(i)*zcrp)/(zvac-zcrp)
                        alpha=(dvduvac(i)-beta)/zvac
                        CALL extrapolfunc%READ([alpha,beta])
                        dvdu(i)=extrapolfunc%getvalue(x(3))
                  END SELECT
               END DO
               DEALLOCATE(extrapolfunc)
               RETURN

            CASE DEFAULT
               WRITE(0,*) "GET_V_AND_DERIVS_PES_H2LiF001 ERR: type of extrapolation function isn't implemented yet"
               WRITE(0,*) "Implemented ones: Linear, Xexponential, None"
               CALL EXIT(1)
         END SELECT
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   ! Check if we are in the Vacuum region
   SELECT CASE(x(3)>=zvac) !easy
      CASE(.TRUE.)
         v=this%farpot%getpot(x(4))
         dvdu(1:3)=0.D0
         dvdu(4)=this%farpot%getderiv(x(4))
         dvdu(5:6)=0.D0
      CASE(.FALSE.) ! this's the last switch!
          errCode=1_1
   END SELECT
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_PES_H2LiF001
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
   CALL xyinterpol%INTERPOL(system_surface)
   ALLOCATE(aux1(5))
   ALLOCATE(aux2(5,2))
   CALL xyinterpol%GET_F_AND_DERIVS(system_surface,x(1:2),aux1,aux2)
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
!# SUBROUTINE: RAWINTERPOL_PES_H2LiF001
!###########################################################
SUBROUTINE RAWINTERPOL_PES_H2LiF001(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),INTENT(INOUT)::this
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   ! Run section
   CALL this%EXTRACT_VACUUMSURF()   
   DO i = 1, this%nsites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         CALL this%wyckoffsite(i)%zrcut(j)%INTERPOL()
      END DO
   END DO
   this%is_interpolated=.TRUE.
   RETURN
END SUBROUTINE RAWINTERPOL_PES_H2LiF001
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
!# SUBROUTINE: ROUGH_PES_H2LiF001
!###########################################################
SUBROUTINE ROUGH_PES_H2LiF001(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
    CLASS(PES_H2LiF001),INTENT(INOUT):: this
   ! Local variables
   REAL(KIND=8),DIMENSION(6):: molcoord,atomcoord
   INTEGER(KIND=4):: nr,nz
   INTEGER(KIND=4):: i,j,k,l ! counters
   REAL(KIND=8):: newpot
   ! Parameters
   CHARACTER(LEN=*),PARAMETER:: routinename="ROUGH_PES_H2LiF001: "
   ! Run section
   DO i = 1, this%nsites ! cycle wyckoff sites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         molcoord(1)=this%wyckoffsite(i)%zrcut(j)%x
         molcoord(2)=this%wyckoffsite(i)%zrcut(j)%y
         molcoord(5)=this%wyckoffsite(i)%zrcut(j)%theta
         molcoord(6)=this%wyckoffsite(i)%zrcut(j)%phi
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizer()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizez()
         DO k = 1, nr
            DO l = 1, nz
               molcoord(3)=this%wyckoffsite(i)%zrcut(j)%getgridvalueZ(l)
               molcoord(4)=this%wyckoffsite(i)%zrcut(j)%getgridvalueR(k)
               atomcoord(:)=from_molecular_to_atomic(molcoord)
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)
               SELECT CASE(this%natomic)
                  CASE(1)
                     newpot=newpot+this%atomiccrp(1)%getpot(atomcoord(1:3))+&
                        this%atomiccrp(1)%getpot(atomcoord(4:6))
                     CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  CASE(2)
                     newpot=newpot+this%atomiccrp(1)%getpot(atomcoord(1:3))+&
                        this%atomiccrp(2)%getpot(atomcoord(4:6))
                     CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  CASE DEFAULT
                     WRITE(0,*) "SMOOTH_PES_H2LiF001 ERR: Something is wrong with number of atomic crp potentials"
                     CALL EXIT(1)
               END SELECT
            END DO
         END DO
      END DO
   END DO
   this%is_smooth=.FALSE.
   RETURN
END SUBROUTINE ROUGH_PES_H2LiF001
!###########################################################
!# SUBROUTINE: SMOOTH_EXTRA_PES_H2LiF001
!###########################################################
SUBROUTINE SMOOTH_EXTRA_PES_H2LiF001(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
    CLASS(PES_H2LiF001),INTENT(INOUT) :: this
   ! Local variables
   REAL(KIND=8),DIMENSION(6) :: molcoord,atomcoord
   INTEGER(KIND=4) :: nr,nz
   INTEGER(KIND=4) :: i,j,k,l ! counters
   CHARACTER(LEN=20),PARAMETER :: routinename="SMOOTH_EXTRA_PES_H2LiF001: "
   REAL(KIND=8) :: newpot
   ! Run section
   DO i = 1, this%nsites ! cycle wyckoff sites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         molcoord(1)=this%wyckoffsite(i)%zrcut(j)%x
         molcoord(2)=this%wyckoffsite(i)%zrcut(j)%y
         molcoord(5)=this%wyckoffsite(i)%zrcut(j)%theta
         molcoord(6)=this%wyckoffsite(i)%zrcut(j)%phi
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizer()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizez()
         DO k = 1, nr
            DO l = 1, nz
               molcoord(3)=this%wyckoffsite(i)%zrcut(j)%getgridvalueZ(l)
               molcoord(4)=this%wyckoffsite(i)%zrcut(j)%getgridvalueR(k)
               atomcoord(:)=from_molecular_to_atomic(molcoord)
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)
               SELECT CASE(this%natomic)
                  CASE(1)
                     newpot=newpot-this%atomiccrp(1)%getpot(atomcoord(1:3))&
                        -this%atomiccrp(1)%getpot(atomcoord(4:6))-this%farpot%getpot(molcoord(4))&
                        -0.8*dexp(-molcoord(4))
                     CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  CASE(2)
                     newpot=newpot-this%atomiccrp(1)%getpot(atomcoord(1:3))-&
                        this%atomiccrp(2)%getpot(atomcoord(4:6))-this%farpot%getpot(molcoord(4))
                     CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  CASE DEFAULT
                     WRITE(0,*) "SMOOTH_EXTRA_PES_H2LiF001 ERR: Something is wrong with number of atomic crp potentials"
                     CALL EXIT(1)
               END SELECT
            END DO
         END DO
      END DO
   END DO
   this%is_smooth=.TRUE.
   RETURN
END SUBROUTINE SMOOTH_EXTRA_PES_H2LiF001
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
   CHARACTER(LEN=26),PARAMETER :: routinename="EXTRACT_VACUUMSURF_PES_H2LiF001: "
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
   CHARACTER(LEN=22),PARAMETER :: routinename="ADD_VACUUMSURF_PES_H2LiF001: "
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
subroutine READ_PES_H2LiF001(this,filename,tablename)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(Pes_H2lif001),intent(out):: this
   character(len=*),intent(in):: filename,tablename
   ! local variables
   character(len=*),parameter:: routinename="READ_PES_H2LiF001: "
   real(kind=8),dimension(:),allocatable:: gridR,gridZ
   real(kind=8),dimension(:,:),allocatable:: gridPot
   ! Run section -----------------------
   call this%set_pesType('PES_H2LiF001')
   call this%set_alias('H2_on_LiF001')
   call this%set_dimensions(6)
   this%is_homonuclear=.true.
   call this%atomicCrp(1)%initialize()
   call this%farPot%initialize('Numerical')
   call this%dumpFunc%read([3.d0,4.d0])
   this%zvacuum=13.2280829302d0
   this%xyklist(:,1)=[0,1,1,2]
   this%xyklist(:,2)=[0,0,1,0]
   ! Create wyckoff sites//////////////////////////////////
   ! Wickoff Top Li -----------> 'a'
   this%wyckoffSite(1)%id='a'
   this%wyckoffSite(1)%myNumber=1
   this%wyckoffSite(1)%is_homonuclear=.true.
   this%wyckoffSite(1)%alias='Top Li'
   this%wyckoffSite(1)%x=0.d0
   this%wyckoffSite(1)%y=0.d0
   this%wyckoffSite(1)%n2dcuts=5
   allocate( this%wyckoffSite(1)%nPhiPoints(3) )
   this%wyckoffSite(1)%nPhiPoints(:)=[1,2,2]
   allocate( this%wyckoffSite(1)%zrCut(5) )
   ! Reading zrcuts
   this%wickoffSite(1)%zrcut(1)%alias='Top_Li_0_0'
   this%wickoffSite(1)%zrCut(1)%x=0.d0
   this%wickoffSite(1)%zrCut(1)%y=0.d0
   this%wickoffSite(1)%zrCut(1)%theta=0.d0
   this%wickoffSite(1)%zrCut(1)%phi=0.d0
   allocate( gridR(14) )
   allocate( gridZ(16) )
   allocate( gridPot(14,16) )
   !---------------------------------------------------------------------------------------------------------------------------
   gridR(:)=[ 0.7558904532d0,0.9448630664d0,1.1338356797d0,1.2283219864d0,1.3228082930d0,1.4172945997d0,1.5117809063d0,&
              1.6062672130d0,1.8897261329d0,2.3621576661d0,2.8345891993d0,3.3070207325d0,3.7794522658d0,4.3463701056d0 ]

   gridZ(:)=[ 0.4724315332d0,0.9448630664d0,1.3983973383d0,1.8897261329d0,2.3054658821d0,2.8345891993d0,&
              3.4015070392d0,3.9684248791d0,4.7243153322d0,5.4802057854d0,6.0471236252d0,6.4250688518d0,&
              6.8030140784d0,7.1809593050d0,7.3699319183d0,7.5589045315d0 ]

   gridPot(:,:)=reshape( [ -2.3486329404d0, 1.0844017627d0, 8.8266781296d0, 9.1272988309d0, 6.5171540893d0, 1.6638643184d0,&
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
                           -7.1750937403d0,-7.1507212202d0 ], shape( gridPot ) )
   !---------------------------------------------------------------------------------------------------------------------------
   call this%wyckoffSite(1)%zrCut(1)%interRZ%read( gridR,gridZ,gridPot )




   return
end subroutine READ_PES_H2LiF001
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
   CHARACTER(LEN=18),PARAMETER :: routinename = "PLOT1D_PHI_PES_H2LiF001: "
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
   CHARACTER(LEN=18),PARAMETER :: routinename = "PLOT1D_PHI_PES_H2LiF001: "
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
! SUBROUTINE: PLOT1D_PHI_SMOOTH_PES_H2LiF001 #######################################
!#######################################################################
SUBROUTINE PLOT1D_PHI_SMOOTH_PES_H2LiF001(thispes,npoints,X,filename)
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
   CHARACTER(LEN=25),PARAMETER :: routinename = "PLOT1D_PHI_SMOOTH_PES_H2LiF001: "
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
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(6),v,dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(6)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(6),v,dvdu(:)
   END DO
   ! Final value
   r(6) = xmax
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(6),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_PHI_SMOOTH_PES_H2LiF001
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
   CHARACTER(LEN=20),PARAMETER :: routinename = "PLOT1D_THETA_PES_H2LiF001: "
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
   CHARACTER(LEN=20),PARAMETER :: routinename = "PLOT1D_THETA_PES_H2LiF001: "
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
! SUBROUTINE: PLOT1D_THETA_SMOOTH_PES_H2LiF001 #########################################
!#######################################################################
SUBROUTINE PLOT1D_THETA_SMOOTH_PES_H2LiF001(thispes,npoints,X,filename)
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
   CHARACTER(LEN=27),PARAMETER :: routinename = "PLOT1D_THETA_SMOOTH_PES_H2LiF001: "
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
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(5),v,dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(5)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(5),v,dvdu(:)
   END DO
   ! Final value
   r(5) = xmax
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(5),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_THETA_SMOOTH_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT1D_R_PES_H2LiF001 #########################################
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
   CHARACTER(LEN=16),PARAMETER :: routinename = "PLOT1D_R_PES_H2LiF001: "
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
! SUBROUTINE: PLOT1D_R_SMOOTH_PES_H2LiF001 #########################################
!#######################################################################
SUBROUTINE PLOT1D_R_SMOOTH_PES_H2LiF001(thispes,npoints,X,filename)
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
   CHARACTER(LEN=23),PARAMETER :: routinename = "PLOT1D_R_SMOOTH_PES_H2LiF001: "
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
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(4),v,dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(4)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(4),v,dvdu(:)
   END DO
   ! Final value
   r(4) = xmax
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(4),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_R_SMOOTH_PES_H2LiF001
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
   CHARACTER(LEN=16),PARAMETER :: routinename = "PLOT1D_Z_PES_H2LiF001: "
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
! SUBROUTINE: PLOT1D_Z_SMOOTH_PES_H2LiF001 #########################################
!#######################################################################
SUBROUTINE PLOT1D_Z_SMOOTH_PES_H2LiF001(thispes,npoints,X,filename)
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
   CHARACTER(LEN=23),PARAMETER :: routinename = "PLOT1D_Z_SMOOTH_PES_H2LiF001: "
   ! HE HO ! LET'S GO ----------------------------
   SELECT CASE(npoints)
      CASE(: 1)
         WRITE(0,*) "PLOT1D_Z_PES_H2LiF001 ERR: Less than 2 points"
         CALL EXIT(1)
      CASE DEFAULT
         ! do nothing
   END SELECT
   !
   xmin = thispes%wyckoffsite(1)%zrcut(1)%getfirstz()
   xmax = thispes%wyckoffsite(1)%zrcut(1)%getlastz()
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
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(3)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(3),v,dvdu(:)
   END DO
   ! Final value
   r(3) = xmax
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_Z_SMOOTH_PES_H2LiF001
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

END MODULE PES_H2LiF001_MOD
