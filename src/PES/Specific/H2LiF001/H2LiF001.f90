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
   real(kind=8),dimension(14,14):: gridPot1414
   real(kind=8),dimension(14,16):: gridPot1416
   real(kind=8),dimension(14):: gridR14,gridZ14
   real(kind=8),dimension(16):: gridR16
   data gridR14(:)/0.7558904532d0,0.9448630664d0,1.1338356797d0,1.2283219864d0,1.3228082930d0,1.4172945997d0,1.5117809063d0&
                   1.6062672130d0,1.8897261329d0,2.3621576661d0,2.8345891993d0,3.3070207325d0,3.7794522658d0,4.3463701056d0/
   data gridZ14(:)/0.4724315332d0,0.9448630664d0,1.8897261329d0,2.8345891993d0,3.4015070392d0,3.9684248791d0,4.7243153322d0,&
                   5.4802057854d0,6.0471236252d0,6.4250688518d0,6.8030140784d0,7.1809593050d0,7.3699319183d0,7.5589045315d0/
   data gridZ16(:)/0.4724315332d0,0.9448630664d0,1.3983973383d0,1.8897261329d0,2.3054658821d0,2.8345891993d0,&
                   3.4015070392d0,3.9684248791d0,4.7243153322d0,5.4802057854d0,6.0471236252d0,6.4250688518d0,&
                   6.8030140784d0,7.1809593050d0,7.3699319183d0,7.5589045315d0/
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
   this%wickoffSite(1)%zrCut(:)%x=0.d0
   this%wickoffSite(1)%zrCut(:)%y=0.d0
   ! Reading zrcuts
   this%wickoffSite(1)%zrcut(1)%alias='Top_Li_0_0'
   this%wickoffSite(1)%zrCut(1)%theta=0.d0
   this%wickoffSite(1)%zrCut(1)%phi=0.d0
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
   this%wickoffSite(1)%zrcut(2)%alias='Top_Li_45_0'
   this%wickoffSite(1)%zrCut(2)%theta=0.785398163397d0
   this%wickoffSite(1)%zrCut(2)%phi=0.d0
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
   this%wickoffSite(1)%zrcut(3)%alias='Top_Li_45_45'
   this%wickoffSite(1)%zrCut(3)%theta=0.785398163397d0
   this%wickoffSite(1)%zrCut(3)%phi=0.785398163397d0
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
   this%wickoffSite(1)%zrcut(4)%alias='Top_Li_90_0'
   this%wickoffSite(1)%zrCut(4)%theta=1.57079632679d0
   this%wickoffSite(1)%zrCut(4)%phi=0.d0
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
   this%wickoffSite(1)%zrcut(5)%alias='Top_Li_90_45'
   this%wickoffSite(1)%zrCut(5)%theta=1.57079632679d0
   this%wickoffSite(1)%zrCut(5)%phi=0.785398163397d0
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
   this%wyckoffSite(2)%is_homonuclear=.true.
   this%wyckoffSite(2)%alias='Top F'
   this%wyckoffSite(2)%x=2.7216780628885480d0
   this%wyckoffSite(2)%y=2.7216780628885480d0
   this%wyckoffSite(2)%n2dcuts=5
   allocate( this%wyckoffSite(2)%nPhiPoints(3) )
   this%wyckoffSite(2)%nPhiPoints(:)=[1,2,2]
   allocate( this%wyckoffSite(2)%zrCut(5) )
   this%wyckoffSite(2)%zrCut(:)%x=2.7216780628885480d0
   this%wyckoffSite(2)%zrCut(:)%y=2.7216780628885480d0
   ! reading zrcuts
   this%wickoffSite(2)%zrcut(1)%alias='Top_F_0_0'
   this%wickoffSite(2)%zrCut(1)%theta=0.d0
   this%wickoffSite(2)%zrCut(1)%phi=0.d0
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
   this%wickoffSite(2)%zrcut(2)%alias='Top_F_45_0'
   this%wickoffSite(2)%zrCut(2)%theta=0.785398163397d0
   this%wickoffSite(2)%zrCut(2)%phi=0.d0
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
   this%wickoffSite(2)%zrcut(3)%alias='Top_F_45_45'
   this%wickoffSite(2)%zrCut(3)%theta=0.785398163397d0
   this%wickoffSite(2)%zrCut(3)%phi=0.785398163397d0
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
   this%wickoffSite(2)%zrcut(4)%alias='Top_F_90_0'
   this%wickoffSite(2)%zrCut(4)%theta=1.57079632679d0
   this%wickoffSite(2)%zrCut(4)%phi=0.d0
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
   call this%wyckofrSite(2)%zrCut(4)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Fifth cut2d
   this%wickoffSite(2)%zrcut(5)%alias='Top_F_90_45'
   this%wickoffSite(2)%zrCut(5)%theta=1.57079632679d0
   this%wickoffSite(2)%zrCut(5)%phi=0.785398163397d0
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
   call this%wyckofrSite(2)%zrCut(5)%interRZ%read( gridR14,gridZ14,gridPot1414 )









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
