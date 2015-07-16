!#######################################################
! MODULE: CRP6D_MOD
!#######################################################
!> @brief
!! Provides all structures and procedures to perform a complete
!! CRP-6D interpolation job. This method only works systematically with diatomic
!! molecules
!#######################################################
MODULE CRP6D_MOD
! Massive name spaces
use SYSTEM_MOD
use LINK_FUNCTION1D_MOD
! Selective name spaces
use UNITS_MOD, only: Length, pi
use PES_MOD, only: PES
use CRP3D_MOD, only: CRP3D
use EXTRAPOL_TO_VACUUM_MOD, only: Vacuumpot
use FOURIER_P4MM_MOD, only: Fourierp4mm
use FOURIER3D_P4MM_MOD, only: Fourier3d_p4mm
use WYCKOFF_P4MM_MOD, only: WyckoffSitio,Wyckoffp4mm,TermsInfo
use FOURIER3D_P4MM_MOD, only: Fourier3d_p4mm
use AOTUS_MODULE, only: flu_State, OPEN_CONFIG_FILE, CLOSE_CONFIG, AOT_GET_VAL
use AOT_TABLE_MODULE, only: AOT_TABLE_OPEN, AOT_TABLE_CLOSE, AOT_TABLE_LENGTH, AOT_TABLE_GET_VAL
#if DEBUG
use DEBUG_MOD, only: verbose_write, debug_write
#endif
IMPLICIT NONE
!/////////////////////////////////////////////////
! TYPE: CRP6D
!
!-------------------------------------------------
TYPE,EXTENDS(PES):: CRP6D
   INTEGER(KIND=4):: nsites
   INTEGER(KIND=4):: natomic
   LOGICAL:: is_interpolated=.FALSE.
   LOGICAL:: is_homonucl=.FALSE.
   LOGICAL:: is_smooth=.FALSE.
   LOGICAL:: is_shifted=.FALSE.
   LOGICAL:: is_resized=.FALSE.
   REAL(KIND=8):: zvacuum
   INTEGER(KIND=4),DIMENSION(2):: grid=[0,0]
   CLASS(Wyckoffsitio),DIMENSION(:),ALLOCATABLE:: wyckoffsite
   TYPE(CRP3D),DIMENSION(:),ALLOCATABLE:: atomiccrp
   TYPE(Vacuumpot):: farpot
   CLASS(Function1d),ALLOCATABLE:: dumpfunc
   CHARACTER(LEN=30):: extrapol2vac_flag
   INTEGER(KIND=4),DIMENSION(:,:),ALLOCATABLE:: kListXY
   character(len=1),dimension(:),allocatable:: parityListXY
   character(len=2),dimension(:),allocatable:: irrepListXY
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE:: kListAngle
   character(len=1),dimension(:),allocatable:: parityListAngle
   character(len=2),dimension(:),allocatable:: irrepListAngle
   integer(kind=4):: totalTerms=0

   CONTAINS
      ! Initialization block
      PROCEDURE,PUBLIC:: READ => READ_CRP6D
      PROCEDURE,PUBLIC:: INITIALIZE => INITIALIZE_CRP6D
      ! Set block
      PROCEDURE,PUBLIC:: SET_SMOOTH => SET_SMOOTH_CRP6D
      ! Get block
      PROCEDURE,PUBLIC:: GET_V_AND_DERIVS => GET_V_AND_DERIVS_CRP6D
      PROCEDURE,PUBLIC:: GET_V_AND_DERIVS_PURE => GET_V_AND_DERIVS_PURE_CRP6D
      PROCEDURE,PUBLIC:: GET_V_AND_DERIVS_SMOOTH => GET_V_AND_DERIVS_SMOOTH_CRP6D
      PROCEDURE,PUBLIC:: GET_ATOMICPOT_AND_DERIVS => GET_ATOMICPOT_AND_DERIVS_CRP6D
      ! Tools block
      PROCEDURE,PUBLIC:: SMOOTH => SMOOTH_CRP6D
      PROCEDURE,PUBLIC:: SMOOTH_EXTRA => SMOOTH_EXTRA_CRP6D
      PROCEDURE,PUBLIC:: INTERPOL => INTERPOL_CRP6D
      PROCEDURE,PUBLIC:: RAWINTERPOL => RAWINTERPOL_CRP6D
      PROCEDURE,PUBLIC:: EXTRACT_VACUUMSURF => EXTRACT_VACUUMSURF_CRP6D
      PROCEDURE,PUBLIC:: ADD_VACUUMSURF => ADD_VACUUMSURF_CRP6D
      PROCEDURE,PUBLIC:: INTERPOL_NEW_RZGRID => INTERPOL_NEW_RZGRID_CRP6D
      PROCEDURE,PUBLIC:: CHEAT_CARTWHEEL_ONTOP => CHEAT_CARTWHEEL_ONTOP_CRP6D
      ! Plot toolk
      PROCEDURE,PUBLIC:: PLOT1D_THETA => PLOT1D_THETA_CRP6D
      PROCEDURE,PUBLIC:: PLOT1D_THETA_SMOOTH => PLOT1D_THETA_SMOOTH_CRP6D
      PROCEDURE,PUBLIC:: PLOT1D_ATOMIC_INTERAC_THETA => PLOT1D_ATOMIC_INTERAC_THETA_CRP6D
      PROCEDURE,PUBLIC:: PLOT1D_PHI => PLOT1D_PHI_CRP6D
      PROCEDURE,PUBLIC:: PLOT1D_PHI_SMOOTH => PLOT1D_PHI_SMOOTH_CRP6D
      PROCEDURE,PUBLIC:: PLOT1D_ATOMIC_INTERAC_PHI => PLOT1D_ATOMIC_INTERAC_PHI_CRP6D
      PROCEDURE,PUBLIC:: PLOT1D_R => PLOT1D_R_CRP6D
      PROCEDURE,PUBLIC:: PLOT1D_R_SMOOTH => PLOT1D_R_SMOOTH_CRP6D
      PROCEDURE,PUBLIC:: PLOT1D_Z => PLOT1D_Z_CRP6D
      PROCEDURE,PUBLIC:: PLOT1D_Z_SMOOTH => PLOT1D_Z_SMOOTH_CRP6D
      PROCEDURE,PUBLIC:: PLOT_XYMAP => PLOT_XYMAP_CRP6D
      PROCEDURE,PUBLIC:: PLOT_XYMAP_SMOOTH => PLOT_XYMAP_SMOOTH_CRP6D
      PROCEDURE,PUBLIC:: PLOT_RZMAP => PLOT_RZMAP_CRP6D
      PROCEDURE,PUBLIC:: PLOT_ATOMIC_INTERAC_RZ => PLOT_ATOMIC_INTERAC_RZ_CRP6D
      ! Enquire block
      PROCEDURE,PUBLIC:: is_allowed => is_allowed_CRP6D
END TYPE CRP6D
CONTAINS
!###########################################################
!# SUBROUTINE: INITIALIZE_CRP6D 
!###########################################################
!> @brief
!! Specific implementation of initialize PES from file
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INITIALIZE_CRP6D(this,filename,tablename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(OUT)::this
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
END SUBROUTINE INITIALIZE_CRP6D
!###########################################################
!# SUBROUTINE: GET_ATOMICPOT_AND_DERIVS_CRP6D 
!###########################################################
!> @brief
!! Get atomic potentials and derivatives
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_ATOMICPOT_AND_DERIVS_CRP6D(this,molecx,atomicx,v,dvdu)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(IN):: this
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
         WRITE(0,*) "GET_ATOMICPOT_AND_DERIVS_CRP6D ERR: wrong number of atomic potentials"
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
END SUBROUTINE GET_ATOMICPOT_AND_DERIVS_CRP6D
!###########################################################
!# SUBROUTINE: SET_SMOOTH_CRP6D
!###########################################################
!> @brief
!! Common set function. Sets is_smooth atribute
!-----------------------------------------------------------
SUBROUTINE SET_SMOOTH_CRP6D(this,is_smooth)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(INOUT):: this
   LOGICAL,INTENT(IN) :: is_smooth
   ! Run section
   this%is_smooth=is_smooth
   RETURN
END SUBROUTINE SET_SMOOTH_CRP6D
!###########################################################
!# SUBROUTINE: CHEAT_CARTWHEEL_ONTOP_CRP6D 
!###########################################################
!> @brief
!! This subroutine changes data from a raw cut2d input by the
!! value interpolated from atomic potentials. Only data too close
!! to the plane @f$ \pi: z=\over{r}{2}+rumpling@f$
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE CHEAT_CARTWHEEL_ONTOP_CRP6D(this,wyckoff,cut2d,toptype,dmax)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O 
   CLASS(CRP6D),INTENT(INOUT) :: this
   INTEGER(KIND=4),INTENT(IN) :: wyckoff,cut2d,toptype
   REAL(KIND=8),INTENT(IN) :: dmax
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   REAL(KIND=8) :: r,z,rump,interac
   REAL(KIND=8) :: x,y ! x,y position
   REAL(KIND=8) :: m1, m2 ! masses
   REAL(KIND=8) :: dist ! distance to the plane
   INTEGER(KIND=4) :: nr,nz
   REAL(KIND=8),DIMENSION(3) :: pos1,pos2
   ! Run section
   nr=this%wyckoffsite(wyckoff)%zrcut(cut2d)%getgridsizeR()
   nz=this%wyckoffsite(wyckoff)%zrcut(cut2d)%getgridsizeZ()
   x=this%wyckoffsite(wyckoff)%x
   y=this%wyckoffsite(wyckoff)%y
   m1=system_mass(1)
   m2=system_mass(2)
   rump=this%atomiccrp(1)%getrumpling(toptype)
   DO i = 1, nr
      DO j = 1, nz
        r=this%wyckoffsite(wyckoff)%zrcut(cut2d)%getgridvalueR(i)
        z=this%wyckoffsite(wyckoff)%zrcut(cut2d)%getgridvalueZ(j)
        dist=(2.D0/sqrt(5.D0))*(z+rump-r/2.D0)
        pos1=(/x,y,z+(m2/(m1+m2))*r/)
        pos2=(/x,y,z-(m1/(m1+m2))*r/)
        SELECT CASE(dist <= dmax) 
           CASE(.TRUE.)
              SELECT CASE(this%is_homonucl)
                 CASE(.TRUE.)
                    interac=this%atomiccrp(1)%getpot(pos1)+this%atomiccrp(1)%getpot(pos2)!+this%farpot%getpot(r)
                 CASE(.FALSE.)
                    interac=this%atomiccrp(1)%getpot(pos1)+this%atomiccrp(2)%getpot(pos2)!+this%farpot%getpot(r)
              END SELECT
              CALL this%wyckoffsite(wyckoff)%zrcut(cut2d)%CHANGEPOT_AT_GRIDPOINT(i,j,interac)
           CASE(.FALSE.)
              ! do nothing
        END SELECT
      END DO
   END DO
   RETURN
END SUBROUTINE CHEAT_CARTWHEEL_ONTOP_CRP6D
!###########################################################
!# SUBROUTINE: INTERPOL_NEW_RZGRID_CRP6D 
!###########################################################
!> @brief
!! Creates a different R-Z grid using the ones stored in files
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INTERPOL_NEW_RZGRID_CRP6D(this,nRpoints,nZpoints)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(INOUT) :: this
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
END SUBROUTINE INTERPOL_NEW_RZGRID_CRP6D
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_PURE_CRP6D 
!###########################################################
!> @brief
!! Gets the potential and derivatives respect to all DOFs if we
!! are in the pure CRP6D region. Otherwise, there'd be errors and
!! unexpected behaviors.
!
!> @warning
!! - Assumed a Fourier interpolation in XY (symmetry adapted)
!! - Assumed a Fourier interpolation in theta (symmetry adapted)
!! - Assumed a Fourier interpolation in phi (symmetry adapted)
!! - Assumed a bicubic splines interpolation in R-Z (symmetry independent)
!! - Assumed that Z-R grids have the same size
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_PURE_CRP6D(this,x,v,dvdu)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),TARGET,INTENT(IN):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN):: x
   REAL(KIND=8),INTENT(OUT):: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT):: dvdu
   ! Local variables
   INTEGER(KIND=4):: i,j,h ! counters
   REAL(KIND=8):: ma,mb
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: f ! smooth function and derivs
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: xyList
   real(kind=8),dimension(:),allocatable:: phiList
   real(kind=8),dimension(:,:),allocatable:: completeGeomList
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: aux1
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: aux2
   REAL(KIND=8),DIMENSION(6):: atomicx
   REAL(KIND=8),DIMENSION(2):: atomic_v
   REAL(KIND=8),DIMENSION(6):: atomic_dvdu
   REAL(KIND=8),DIMENSION(3):: dvdu_atomicA,dvdu_atomicB
   TYPE(Fourier3d_p4mm):: fouInterpol
   CHARACTER(LEN=*),PARAMETER:: routinename="GET_V_AND_DERIVS_CRP6D: "
   ! Run section
   SELECT CASE(this%is_smooth)
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         WRITE(0,*) "GET_V_AND_DERIVS_CRP6D ERR: smooth the PES first (CALL thispes%SMOOTH())"
         CALl EXIT(1)
   END SELECT
   ! f(1,:) smooth potential values
   ! f(2,:) smooth dvdz
   ! f(3,:) smooth dvdr
   ! f(4,:) smooth dvdtheta
   ! f(5,:) smooth dvdphi
   allocate( f(5,this%totalTerms)      )
   allocate( xyList(this%totalTerms,2) )
   allocate( phiList(this%totalTerms)  )
   allocate( completeGeomList(this%totalTerms,3)  )
   h=0
   do i = 1, this%nsites
      do j=1,size(this%wyckoffsite(i)%phiList)
         h=h+1
         xyList(h,1)=this%wyckoffSite(i)%x
         xyList(h,2)=this%wyckoffSite(i)%y
         phiList(h)=this%wyckoffSite(i)%phiList(j)
         call this%wyckoffsite(i)%get_v_and_derivs( [x(3),x(4),x(5),phiList(h)],f(1,h),f(2:5,h) )
      enddo
   end do
   completeGeomList(:,1:2)=xyList(:,:)
   completeGeomList(:,3)=phiList(:)
   call fouInterpol%read( x=completeGeomList(:,:),f=f(1,:) )
   call fouInterpol%add_moreFuncs( f=f(2:4,1:this%totalTerms) )
   call fouInterpol%initializeTerms()
   call fouInterpol%setKlist( kListXY=this%kListXY,kListAngle=this%kListAngle )
   call fouInterpol%setParityList( parityListXY=this%parityListXY,parityListAngle=this%parityListAngle )
   call fouInterpol%setIrrepList( irrepListXY=this%irrepListXY,irrepListAngle=this%irrepListAngle )
   CALL fouInterpol%interpol()
   allocate(aux1(4))
   allocate(aux2(4,3))
   call fouInterpol%get_allFuncs_and_derivs( x=[x(1),x(2),x(6)],f=aux1,dfdx=aux2 )
#ifdef DEBUG
   !-------------------------------------
   ! Results for the smooth potential
   !-------------------------------------
   ! v=aux1(1)         ! v
   ! dvdu(1)=aux2(1,1) ! dvdx
   ! dvdu(2)=aux2(1,2) ! dvdy
   ! dvdu(3)=aux1(2)   ! dvdz
   ! dvdu(4)=aux1(3)   ! dvdr
   ! dvdu(5)=aux1(4)   ! dvdtheta
   ! dvdu(6)=aux2(1,3) ! dvdphi
   CALL DEBUG_WRITE(routinename,"Smooth V at each thetaCut:")
   CALL DEBUG_WRITE(routinename,f(1,:))
   CALL DEBUG_WRITE(routinename,"Smooth dv/dz at each thetaCut:")
   CALL DEBUG_WRITE(routinename,f(2,:))
   CALL DEBUG_WRITE(routinename,"Smooth dv/dr at each thetaCut:")
   CALL DEBUG_WRITE(routinename,f(3,:))
   CALL DEBUG_WRITE(routinename,"Smooth dv/dtheta at each thetaCut:")
   CALL DEBUG_WRITE(routinename,f(4,:))
   CALL DEBUG_WRITE(routinename,"Smooth interpolated values:")
   CALL DEBUG_WRITE(routinename,"v: ",aux1(1))   ! v
   CALL DEBUG_WRITE(routinename,"dvdx: ",aux2(1,1)) ! dvdx
   CALL DEBUG_WRITE(routinename,"dvdy: ",aux2(1,2)) ! dvdy
   CALL DEBUG_WRITE(routinename,"dvdz: ",aux1(2))   ! dvdz
   CALL DEBUG_WRITE(routinename,"dvdr: ",aux1(3))   ! dvdr
   CALL DEBUG_WRITE(routinename,"dvdtheta: ",aux1(4))   ! dvdtheta
   CALL DEBUG_WRITE(routinename,"dvdphi: ",aux2(1,3))   ! dvdphi
#endif
   !--------------------------------------
   ! Results for the real potential
   !-------------------------------------
   CALL this%GET_ATOMICPOT_AND_DERIVS(x,atomicx,atomic_v,atomic_dvdu)
   dvdu_atomicA=atomic_dvdu(1:3)
   dvdu_atomicB=atomic_dvdu(4:6)
#ifdef DEBUG
   CALL DEBUG_WRITE(routinename,"Contributions of the atomic potential: ")
   CALL DEBUG_WRITE(routinename, "Position Atom A: ",atomicx(1:3))
   CALL DEBUG_WRITE(routinename,"Va: ",atomic_v(1))
   CALL DEBUG_WRITE(routinename,"dVa/dxa; dVa/dya: ",dvdu_atomicA)
   CALL DEBUG_WRITE(routinename, "Position Atom B: ",atomicx(4:6))
   CALL DEBUG_WRITE(routinename,"Vb: ",atomic_v(2))
   CALL DEBUG_WRITE(routinename,"dVa/dxa; dVb/dyb: ",dvdu_atomicB)
#endif
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
   dvdu(6)=aux2(1,3)&
      -(mb/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicA(1)*dsin(x(6))&
      +(mb/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicA(2)*dcos(x(6))&
      +(ma/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicB(1)*dsin(x(6))&
      -(ma/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicB(2)*dcos(x(6))
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_PURE_CRP6D
!###########################################################
!# SUBROUTINE: GET_V AND DERIVS_CRP6D 
!###########################################################
!> @brief
!! Gets potential for any configuration. Discriminates between 3 regions:
!! - Pure CRP6D region: zmin to zcrp
!! - Extrapolation region zcrp to zvacuum
!! - Vacuum region: further than zvacuum
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_CRP6D(this,X,v,dvdu,errCode)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),TARGET,INTENT(IN):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   REAL(KIND=8),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: dvdu
   integer(kind=1),optional,intent(out):: errCode
   ! Local variables
   REAL(KIND=8) :: zcrp, zvac ! last CRP6D z value and Z infinity
   REAL(KIND=8) :: vzcrp, vzvac ! potentials at zcrp and zvac
   REAL(KIND=8),DIMENSION(6) :: dvducrp ! derivatives at zcrp
   REAL(KIND=8),DIMENSION(6) :: dvduvac ! derivatives at vacuum
   REAL(KIND=8) :: alpha,beta,gama ! parameters
   CLASS(Function1d),ALLOCATABLE:: extrapolfunc
   INTEGER(KIND=4) :: i !counter
   ! Local Parameter
   REAL(KIND=8),PARAMETER :: zero=0.D0 ! what we will consider zero (a.u.)
   REAL(KIND=8),PARAMETER :: dz=0.5D0 ! 0.25 Angstroems approx
   character(len=*),parameter:: routinename='GET_V_AND_DERIVS_CRP6D: '
   ! Run section
   zcrp=this%wyckoffsite(1)%zrcut(1)%getlastZ()
   zvac=this%zvacuum
   ! Check if we are in the pure CRP6D region
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
               WRITE(0,*) "GET_V_AND_DERIVS_CRP6D ERR: type of extrapolation function isn't implemented yet"
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
#ifdef DEBUG
         call debug_write(routinename,'Unclassificable point')
         call debug_write(routinename,'Asking potential at Z: ',x(3))
#endif
          errCode=1_1
   END SELECT
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_CRP6D
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_SMOOTH_CRP6D 
!###########################################################
!> @brief
!! Gets the SMOOTH potential and the derivatives respect to all DOFs
!
!> @details
!! This routine doesn't check PES smoothness, so it can be used to get
!! interpolated values for the raw PES
!> @warning
!! - Assumed a 2D Fourier interpolation in XY (symmetry adapted)
!! - Assumed a 1D Fourier interpolation in theta (symmetry adapted)
!! - Assumed a 1D Fourier interpolation in phi (symmetry adapted)
!! - Assumed a bicubic splines interpolation in R-Z (symmetry independent)
!! - Assumed that Z-R grids have the same size. This can be ensured using
!!   using subroutine INTERPOL_NEWGRID
!
!> @see
!! interpol_newgrid_crp6d
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_SMOOTH_CRP6D(this,x,v,dvdu)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),TARGET,INTENT(IN):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN):: x
   REAL(KIND=8),INTENT(OUT):: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT):: dvdu
   ! Local variables
   INTEGER(KIND=4):: i,j,h ! counters
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: f ! smooth function and derivs
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: xyList
   real(kind=8),dimension(:),allocatable:: phiList
   real(kind=8),dimension(:,:),allocatable:: completeGeomList
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: aux1
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: aux2
   TYPE(Fourier3d_p4mm):: fouInterpol
   CHARACTER(LEN=*),PARAMETER:: routinename="GET_V_AND_DERIVS_CRP6D: "
   ! Run section
   ! f(1,:) smooth potential values
   ! f(2,:) smooth dvdz
   ! f(3,:) smooth dvdr
   ! f(4,:) smooth dvdtheta
   ! f(5,:) smooth dvdphi
   allocate( f(5,this%totalTerms)      )
   allocate( xyList(this%totalTerms,2) )
   allocate( phiList(this%totalTerms)  )
   allocate( completeGeomList(this%totalTerms,3)  )
   h=0
   do i = 1, this%nsites
      do j=1,size(this%wyckoffsite(i)%phiList)
         h=h+1
         xyList(h,1)=this%wyckoffSite(i)%x
         xyList(h,2)=this%wyckoffSite(i)%y
         phiList(h)=this%wyckoffSite(i)%phiList(j)
         call this%wyckoffsite(i)%get_v_and_derivs( [x(3),x(4),x(5),phiList(h)],f(1,h),f(2:5,h) )
      enddo
   end do
   completeGeomList(:,1:2)=xyList(:,:)
   completeGeomList(:,3)=phiList(:)
   call fouInterpol%read( x=completeGeomList(:,:),f=f(1,:) )
   call fouInterpol%add_moreFuncs( f=f(2:4,1:this%totalTerms) )
   call fouInterpol%initializeTerms()
   call fouInterpol%setKlist( kListXY=this%kListXY,kListAngle=this%kListAngle )
   call fouInterpol%setParityList( parityListXY=this%parityListXY,parityListAngle=this%parityListAngle )
   call fouInterpol%setIrrepList( irrepListXY=this%irrepListXY,irrepListAngle=this%irrepListAngle )
   CALL fouInterpol%interpol()
   allocate(aux1(4))
   allocate(aux2(4,3))
   call fouInterpol%get_allFuncs_and_derivs( x=[x(1),x(2),x(6)],f=aux1,dfdx=aux2 )
#ifdef DEBUG
   !-------------------------------------
   ! Results for the smooth potential
   !-------------------------------------
   ! v=aux1(1)         ! v
   ! dvdu(1)=aux2(1,1) ! dvdx
   ! dvdu(2)=aux2(1,2) ! dvdy
   ! dvdu(3)=aux1(2)   ! dvdz
   ! dvdu(4)=aux1(3)   ! dvdr
   ! dvdu(5)=aux1(4)   ! dvdtheta
   ! dvdu(6)=aux2(1,3) ! dvdphi
   CALL DEBUG_WRITE(routinename,"Smooth V at each thetaCut:")
   CALL DEBUG_WRITE(routinename,f(1,:))
   CALL DEBUG_WRITE(routinename,"Smooth dv/dz at each thetaCut:")
   CALL DEBUG_WRITE(routinename,f(2,:))
   CALL DEBUG_WRITE(routinename,"Smooth dv/dr at each thetaCut:")
   CALL DEBUG_WRITE(routinename,f(3,:))
   CALL DEBUG_WRITE(routinename,"Smooth dv/dtheta at each thetaCut:")
   CALL DEBUG_WRITE(routinename,f(4,:))
   CALL DEBUG_WRITE(routinename,"Smooth interpolated values:")
   CALL DEBUG_WRITE(routinename,"v: ",aux1(1))   ! v
   CALL DEBUG_WRITE(routinename,"dvdx: ",aux2(1,1)) ! dvdx
   CALL DEBUG_WRITE(routinename,"dvdy: ",aux2(1,2)) ! dvdy
   CALL DEBUG_WRITE(routinename,"dvdz: ",aux1(2))   ! dvdz
   CALL DEBUG_WRITE(routinename,"dvdr: ",aux1(3))   ! dvdr
   CALL DEBUG_WRITE(routinename,"dvdtheta: ",aux1(4))   ! dvdtheta
   CALL DEBUG_WRITE(routinename,"dvdphi: ",aux2(1,3))   ! dvdphi
#endif
   v=aux1(1)
   dvdu(1)=aux2(1,1)
   dvdu(2)=aux2(1,2)
   dvdu(3)=aux1(2)
   dvdu(4)=aux1(3)
   dvdu(5)=aux1(4)
   dvdu(6)=aux2(1,3)
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_SMOOTH_CRP6D
!###########################################################
!# SUBROUTINE: INTERPOL_CRP6D 
!###########################################################
!> @brief
!! Interpolates CRP6D smooth 2dcuts in R and Z.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INTERPOL_CRP6D(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(INOUT)::this
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
END SUBROUTINE INTERPOL_CRP6D
!###########################################################
!# SUBROUTINE: RAWINTERPOL_CRP6D 
!###########################################################
!> @brief
!! Interpolates CRP6D 2dcuts in R and Z. There is not smoothing
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE RAWINTERPOL_CRP6D(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(INOUT)::this
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
END SUBROUTINE RAWINTERPOL_CRP6D
!###########################################################
!# SUBROUTINE: SMOOTH_CRP6D
!###########################################################
!> @brief
!! Smooths a RZ-2dcut of the potential using atomic potentials
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 25/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SMOOTH_CRP6D(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
    CLASS(CRP6D),INTENT(INOUT):: this
   ! Local variables
   REAL(KIND=8),DIMENSION(6):: molcoord,atomcoord,dummy
   INTEGER(KIND=4):: nr,nz
   INTEGER(KIND=4):: i,j,k,l ! counters
   CHARACTER(LEN=*),PARAMETER:: routinename="SMOOTH_CRP6D: "
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
END SUBROUTINE SMOOTH_CRP6D
!###########################################################
!# SUBROUTINE: ROUGH_CRP6D
!###########################################################
!> @brief
!! ROUGHs a RZ-2dcut of the potential using atomic potentials
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE ROUGH_CRP6D(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
    CLASS(CRP6D),INTENT(INOUT):: this
   ! Local variables
   REAL(KIND=8),DIMENSION(6):: molcoord,atomcoord
   INTEGER(KIND=4):: nr,nz
   INTEGER(KIND=4):: i,j,k,l ! counters
   REAL(KIND=8):: newpot
   ! Parameters
   CHARACTER(LEN=*),PARAMETER:: routinename="ROUGH_CRP6D: "
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
#ifdef DEBUG
               CALL DEBUG_WRITE(routinename,"Molecular coords:")
               CALL DEBUG_WRITE(routinename,molcoord)
               CALL DEBUG_WRITE(routinename,"Atomic coords:")
               CALL DEBUG_WRITE(routinename,atomcoord)
#endif
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
                     WRITE(0,*) "SMOOTH_CRP6D ERR: Something is wrong with number of atomic crp potentials"
                     CALL EXIT(1)
               END SELECT
            END DO
         END DO
      END DO
   END DO
   this%is_smooth=.FALSE.
   RETURN
END SUBROUTINE ROUGH_CRP6D
!###########################################################
!# SUBROUTINE: SMOOTH_EXTRA_CRP6D
!###########################################################
!> @brief
!! Smooths a RZ-2dcut of the potential using atomic potentials as well
!! as the vacuum potential
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SMOOTH_EXTRA_CRP6D(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
    CLASS(CRP6D),INTENT(INOUT) :: this
   ! Local variables
   REAL(KIND=8),DIMENSION(6) :: molcoord,atomcoord
   INTEGER(KIND=4) :: nr,nz
   INTEGER(KIND=4) :: i,j,k,l ! counters
   CHARACTER(LEN=20),PARAMETER :: routinename="SMOOTH_EXTRA_CRP6D: "
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
#ifdef DEBUG
               CALL DEBUG_WRITE(routinename,"Molecular coords:")
               CALL DEBUG_WRITE(routinename,molcoord)
               CALL DEBUG_WRITE(routinename,"Atomic coords:")
               CALL DEBUG_WRITE(routinename,atomcoord)
#endif
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
                     WRITE(0,*) "SMOOTH_EXTRA_CRP6D ERR: Something is wrong with number of atomic crp potentials"
                     CALL EXIT(1)
               END SELECT
            END DO
         END DO
      END DO
   END DO
   this%is_smooth=.TRUE.
   RETURN
END SUBROUTINE SMOOTH_EXTRA_CRP6D
!###########################################################
!# SUBROUTINE: EXTRACT_VACUUMSURF_CRP6D
!###########################################################
!> @brief
!! Extracts energy at equilibrium distance in the vacuum and surface energy
!! to an Rz-2dcut of the potential. Shifts the vacuum potential as well.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 25/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE EXTRACT_VACUUMSURF_CRP6D(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
    CLASS(CRP6D),INTENT(INOUT) :: this
   ! Local variables
   INTEGER(KIND=4) :: nr,nz
   INTEGER(KIND=4) :: i,j,k,l ! counters
   REAL(KIND=8) :: newpot
   CHARACTER(LEN=26),PARAMETER :: routinename="EXTRACT_VACUUMSURF_CRP6D: "
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
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Potential shifted: ",this%farpot%getscalefactor())
#endif
   this%is_shifted=.TRUE.
   RETURN
END SUBROUTINE EXTRACT_VACUUMSURF_CRP6D
!###########################################################
!# SUBROUTINE: ADD_VACUUMSURF_CRP6D
!###########################################################
!> @brief
!! Adds energy at equilibrium distance in the vacuum and surface energy
!! to an Rz-2dcut of the potential. Shifts the vacuum potential as well.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE ADD_VACUUMSURF_CRP6D(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
    CLASS(CRP6D),INTENT(INOUT) :: this
   ! Local variables
   INTEGER(KIND=4) :: nr,nz
   INTEGER(KIND=4) :: i,j,k,l ! counters
   REAL(KIND=8) :: newpot
   CHARACTER(LEN=22),PARAMETER :: routinename="ADD_VACUUMSURF_CRP6D: "
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
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Potential re-shifted: ",this%farpot%getscalefactor())
#endif
   this%is_shifted=.FALSE.
   RETURN
END SUBROUTINE ADD_VACUUMSURF_CRP6D
!###########################################################
!# SUBROUTINE: READ_CRP6D 
!###########################################################
!> @brief
!! Sets up a CRP6D Object
!-----------------------------------------------------------
SUBROUTINE READ_CRP6D(this,filename,tablename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(OUT):: this
   CHARACTER(LEN=*),INTENT(IN):: filename,tablename
   ! Local variables
   INTEGER(KIND=4):: i,j,k,n ! counters
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: param_damp
   ! Wyckoff variables
   CHARACTER(LEN=1),DIMENSION(:),ALLOCATABLE:: wyckoff_letters
   CHARACTER(LEN=1024),DIMENSION(:),ALLOCATABLE:: cuts2d_files
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE:: thetablocks_data
   type(TermsInfo),dimension(:),allocatable:: phiTerms
   type(TermsInfo):: thetaTerms
   INTEGER(KIND=4):: n2dcuts
   INTEGER(KIND=4):: nthetablocks
   real(kind=8),dimension(:),allocatable:: phiList

   ! Lua specifications
   TYPE(flu_State):: conf
   INTEGER(KIND=4):: ierr
   INTEGER(KIND=4):: pes_table,crp3d_table,vacfunc_table,dampfunc_table,param_table,extrapol_table
   INTEGER(KIND=4):: resize_table,wyckoff_table,inwyckoff_table,cut2d_table,files_table,kpoints_table
   integer(kind=4):: phi_table,term_table,theta_table,fourier_table,magnitude_table,phiList_table
   INTEGER(KIND=4):: inkpoints_table
   ! Auxiliar, dummy variables
   INTEGER(KIND=4):: auxInt,auxInt2
   REAL(KIND=8):: auxReal
   CHARACTER(LEN=1024):: auxString,auxString2
   TYPE(Length):: len
   type(Angle):: angl
   ! Parameters
   CHARACTER(LEN=*),PARAMETER:: routinename="READ_CRP6D: "
   ! Run section -----------------------
   ! Open Lua config file
   CALL OPEN_CONFIG_FILE(L=conf,ErrCode=ierr,filename=filename)
   ! Open PES table
   CALL AOT_TABLE_OPEN(L=conf,thandle=pes_table,key=tablename)
   ! get pes.kind
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='kind',val=auxstring)
   CALL this%SET_PESTYPE(trim(auxstring))
   SELECT CASE(trim(auxstring))
      CASE('CRP6D')
         ! do nothing
      CASE DEFAULT
         WRITE(0,*) 'READ CRP6D ERR: wrong kind of PES. Expected: CRP6D. Encountered: '//trim(auxstring)
   END SELECT
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'PES type: '//trim(auxstring))
#endif
   ! get pes.name
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='name',val=auxstring)
   CALL this%SET_ALIAS(trim(auxstring))
#ifdef DEBUG
   CALl VERBOSE_WRITE(routinename,'PES name: '//trim(auxstring))
#endif
   ! get pes.dimensions
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='dimensions',val=auxint)
   CALL this%SET_DIMENSIONS(auxint)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'PES dimensions: ',auxint)
#endif
   ! get crp3d subpes
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=crp3d_table,key='crp3dPes')
   this%natomic=aot_table_length(L=conf,thandle=crp3d_table)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Atomic potentials found: ",this%natomic)
#endif
   SELECT CASE(this%natomic)
      CASE(1)
         this%is_homonucl=.TRUE.
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,'This PES is for homonuclear projectiles')
#endif         
      CASE(2)
         this%is_homonucl=.FALSE.
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,'This PES is for heteronuclear projectiles')
#endif         
      CASE DEFAULT 
         WRITE(0,*) "READ_CRP6D ERR: Wrong number of atomic potentials. Allowed number: 1 or 2."
         CALL EXIT(1)
   END SELECT
   ALLOCATE(this%atomiccrp(this%natomic))
   DO i = 1, this%natomic
      CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=crp3d_table,pos=i,val=auxstring)
#ifdef DEBUG
      CALl VERBOSE_WRITE(routinename,'Atomic potential keyword: '//trim(auxstring))
#endif
      CALL this%atomiccrp(i)%INITIALIZE(filename=filename,tablename=trim(auxstring))
   END DO
   CALL AOT_TABLE_CLOSE(L=conf,thandle=crp3d_table)
   ! get pes.vacuumFunction
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=vacfunc_table,key='vacuumFunction')
   CALL AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=vacfunc_table,key='kind',val=auxstring)
   SELECT CASE(trim(auxstring))
      CASE('Numerical')
         CALL AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=vacfunc_table,key='source',val=auxstring)
         CALL this%farpot%INITIALIZE(trim(auxstring))
      CASE DEFAULT
         WRITE(0,*) 'READ_CRP6D ERR: wrong kind of vacuuum function: '//trim(auxstring)
         WRITE(0,*) 'Implemented ones: Numerical'
         WRITE(0,*) 'Case sensitive'
         CALL EXIT(1)
   END SELECT
   CALL AOT_TABLE_CLOSE(L=conf,thandle=vacfunc_table)
   ! get pes.dampFunction
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=dampfunc_table,key='dampFunction')
   CALL AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=dampfunc_table,key='kind',val=auxstring)
   SELECT CASE(trim(auxstring))
      CASE("Logistic")
         ALLOCATE(Logistic_func::this%dumpfunc)
         ALLOCATE(param_damp(2))
         CALL AOT_TABLE_OPEN(L=conf,parent=dampfunc_table,thandle=param_table,key='param')
         auxint=aot_table_length(L=conf,thandle=param_table)
         SELECT CASE(auxint)
            CASE(2)
               ! do nothing
            CASE DEFAULT
               WRITE(0,*) 'READ_CRP6D ERR: incorrect number of parameters in damp function: ',auxint
               CALL EXIT(1)
         END SELECT
         CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=param_table,pos=1,val=param_damp(1))
         CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=param_table,pos=2,val=param_damp(2))
         CALL AOT_TABLE_CLOSE(L=conf,thandle=param_table)
         CALL this%dumpfunc%READ(param_damp)
      CASE("None") 
         ALLOCATE(One_func::this%dumpfunc)
      CASE DEFAULT
         WRITE(0,*) "READ_CRP6D ERR: Keyword for dumping function needed"
         WRITE(*,*) "Currently implemented: Logistic, None"
         CALL EXIT(1)
   END SELECT
   CALL AOT_TABLE_CLOSE(L=conf,thandle=dampfunc_table)
   ! get pes.extrapolFunction
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=extrapol_table,key='extrapolFunction')
   CALL AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=extrapol_table,key='kind',val=auxstring)
   SELECT CASE(trim(auxstring))
      CASE("Xexponential","Linear")
         this%extrapol2vac_flag=trim(auxstring)
         CALL AOT_TABLE_OPEN(L=conf,parent=extrapol_table,thandle=param_table,key='upToZ')
         auxint=aot_table_length(L=conf,thandle=param_table)
         CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=param_table,pos=1,val=auxreal)
         CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=param_table,pos=2,val=auxstring)
         CALL len%READ(auxreal,trim(auxstring))
         CALL len%TO_STD()
         this%zvacuum=len%getvalue()
         CALl AOT_TABLE_CLOSE(L=conf,thandle=param_table)
      CASE("None")
         this%zvacuum=0.D0
      CASE DEFAULT
         WRITE(0,*) "READ_CRP6D ERR: Keyword for extrapolation function needed"
         WRITE(0,*) "Currently implemented: None, Xexponential, Linear"
         WRITE(0,*) "Case sensitive."
         CALL EXIT(1)
   END SELECT
   CALL AOT_TABLE_CLOSE(L=conf,thandle=extrapol_table)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Z vacuum: ",this%zvacuum)
#endif
   ! get pes.resize
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=resize_table,key='resize')
   auxint=aot_table_length(L=conf,thandle=resize_table)
   SELECT CASE(auxint)
      CASE(0)
         ! do nothing, resize is not required
      CASE DEFAULT
         CALL AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=resize_table,key='r',val=auxint)
         this%grid(1)=auxint
         CALL AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=resize_table,key='z',val=auxint)
         this%grid(2)=auxint
         if( this%grid(1)/=0 .and. this%grid(2)/=0 ) this%is_resized=.true.
   END SELECT
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,'New grid (R,Z):',this%grid(:))
         call verbose_write(routinename,'Is this PES going to be resized?: ',this%is_resized)
#endif
   CALL AOT_TABLE_CLOSE(L=conf,thandle=resize_table)
   ! get wyckoff sites
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=wyckoff_table,key='wyckoffSite')
   this%nsites=aot_table_length(L=conf,thandle=wyckoff_table)
   ALLOCATE(wyckoff_letters(this%nsites))
   SELECT CASE(this%nsites)
      CASE(0)
         WRITE(0,*) "CRP3D_READ ERR: there aren't Wyckoff sites"
         CALL EXIT(1)
      CASE DEFAULT
         ! do nothing
   END SELECT
   ! Allocate with the correct type (symmetry) all wyckoff sites
   SELECT CASE(system_surface%getsymmlabel())
      CASE("p4mm")
         ALLOCATE(Wyckoffp4mm::this%wyckoffsite(this%nsites))
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Allocated p4mm Wyckoff sites")
#endif
      CASE DEFAULT
         WRITE(0,*) "READ_CRP6D ERR: surface symmetry is not implemented yet"
         WRITE(0,*) "Good luck!"
         CALL EXIT(1)
   END SELECT
   DO i = 1, this%nsites
      CALL AOT_TABLE_OPEN(L=conf,parent=wyckoff_table,thandle=inwyckoff_table,pos=i)
      call aot_table_open( L=conf,parent=inWyckoff_table,thandle=phiList_table,key='phiList' )
      allocate( phiList(aot_table_length(L=conf,thandle=phiList_table)) )
      do k=1,size(phiList)
         call aot_table_open( L=conf,parent=phiList_table,thandle=magnitude_table,pos=k )
         call aot_get_val( L=conf,errCode=iErr,thandle=magnitude_table,pos=1,val=auxReal )
         call aot_get_val( L=conf,errCode=iErr,thandle=magnitude_table,pos=2,val=auxString )
         call angl%read( auxReal,trim(auxString) )
         call angl%to_std()
         phiList(k)=angl%getValue()
         call aot_table_close( L=conf,thandle=magnitude_table )
      enddo
      call aot_table_close( L=conf,thandle=phiList_table )
      CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=inwyckoff_table,key='kind',val=wyckoff_letters(i))
      CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=inwyckoff_table,key='n2dcuts',val=n2dcuts)
      ALLOCATE(cuts2d_files(n2dcuts))
      nthetablocks=aot_table_length(L=conf,thandle=inwyckoff_table)-7
      ALLOCATE(thetablocks_data(nthetablocks))
      call aot_table_open(L=conf,parent=inwyckoff_table,thandle=theta_table,key='thetaFourierTerms')
      do k=1,nThetaBlocks
         call aot_table_open( L=conf,parent=theta_table,thandle=term_table,pos=k )
         call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=1,val=auxString  )
         call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=2,val=auxString2 )
         call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=3,val=auxInt     )
         call thetaTerms%addTerm( irrep=trim(auxString),parity=trim(auxString2),kpoint=auxInt )
         call aot_table_close( L=conf,thandle=term_table )
      enddo
      call aot_table_close(L=conf,thandle=theta_table)
#ifdef DEBUG
      CALL VERBOSE_WRITE( routinename,'Wyckoff Site number: ',i)
      CALL VERBOSE_WRITE( routinename,'Wyckoff Letter: '//trim(wyckoff_letters(i)))
      CALL VERBOSE_WRITE( routinename,'Wyckoff number of cut2ds: ',n2dcuts)
      CALl VERBOSE_WRITE( routinename,'Wyckoff number of theta blocks: ',nthetablocks)
      call verbose_write( routinename,'Wyckoff theta terms irreps: ',thetaTerms%irrepList(:) )
      call verbose_write( routinename,'Wyckoff theta terms parities: ',thetaTerms%parityList(:) )
      call verbose_write( routinename,'Wyckoff theta terms kpoints: ',thetaTerms%kpointList(:) )
#endif
      allocate( phiTerms(nThetaBlocks) )
      DO j = 1, nThetaBlocks
         CALL AOT_TABLE_OPEN(L=conf,parent=inwyckoff_table,thandle=cut2d_table,pos=j)
         CALL AOT_TABLE_OPEN(L=conf,parent=cut2d_table,thandle=files_table,key='files')
         thetablocks_data(j)=aot_table_length(L=conf,thandle=files_table)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,'Wyckoff thetablock: ',j)
         CALL VERBOSE_WRITE(routinename,'Wyckoff cut2d inside: ',thetablocks_data(j))
#endif
         SELECT CASE(j)
            CASE(1)
               auxint=1
               auxint2=thetablocks_data(1)
            CASE DEFAULT
               auxint=sum(thetablocks_data(1:j-1))+1
               auxint2=sum(thetablocks_data(1:j-1))+thetablocks_data(j)
         END SELECT
         n=0
         DO k = auxint, auxint2
            n=n+1
            CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=files_table,pos=n,val=cuts2d_files(k))
#ifdef DEBUG
            CALL VERBOSE_WRITE(routinename,'Wyckoff pos array: ',k)
            CALL VERBOSE_WRITE(routinename,'Wyckoff cut2d filename: '//trim(cuts2d_files(k)))
#endif
         END DO
         CALL AOT_TABLE_CLOSE(L=conf,thandle=files_table)
         call aot_table_open( L=conf,parent=cut2d_table,thandle=phi_table,key='phiFourierTerms' )
         do k=1,thetaBlocks_data(j)
            call aot_table_open( L=conf,parent=phi_table,thandle=term_table,pos=k )
            call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=1,val=auxString )
            call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=2,val=auxString2 )
            call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=3,val=auxInt )
            call phiTerms(j)%addTerm( irrep=trim(auxString),parity=trim(auxString2),kpoint=auxInt )
            call aot_table_close( L=conf,thandle=term_table )
         enddo
#ifdef DEBUG
         call verbose_write( routinename,'Wyckoff phi terms irreps: ',  phiTerms(j)%irrepList(:)  )
         call verbose_write( routinename,'Wyckoff phi terms parities: ',phiTerms(j)%parityList(:) )
         call verbose_write( routinename,'Wyckoff phi terms kpoints: ', phiTerms(j)%kpointList(:) )
#endif
         call aot_table_close(L=conf,thandle=phi_table)
         CALL AOT_TABLE_CLOSE(L=conf,thandle=cut2d_table)
      END DO
      CALL this%wyckoffsite(i)%INITIALIZE(mynumber=i,letter=wyckoff_letters(i),nphipoints=thetablocks_data(:),&
                                          filenames=cuts2d_files(:),phiTerms=phiTerms,thetaTerms=thetaTerms,  &
                                          phiList=phiList(:) )
      CALL AOT_TABLE_CLOSE(L=conf,thandle=inwyckoff_table)
      deallocate(cuts2d_files)
      deallocate(thetablocks_data)
      deallocate(phiTerms)
      deallocate(phiList)
      call thetaTerms%reboot()
   END DO
   CALL AOT_TABLE_CLOSE(L=conf,thandle=wyckoff_table)
   ! Get kpoints for Fourier interpolation
   call aot_table_open( L=conf,parent=pes_table,thandle=fourier_table,key='fourierTerms' )
   this%totalTerms=aot_table_length( L=conf,thandle=fourier_table )
   allocate( this%kListXY(this%totalTerms,2)    )
   allocate( this%irrepListXY(this%totalTerms)  )
   allocate( this%parityListXY(this%totalTerms) )
   allocate( this%kListAngle(this%totalTerms)    )
   allocate( this%irrepListAngle(this%totalTerms)  )
   allocate( this%parityListAngle(this%totalTerms) )
   do i=1,this%totalTerms
      call aot_table_open( L=conf,parent=fourier_table,thandle=term_table,pos=i )
      call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=1,val=auxString )
      call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=2,val=auxString2 )
      call aot_table_open( L=conf,parent=term_table,thandle=kpoints_table,pos=3 )
      call aot_get_val( L=conf,errCode=iErr,thandle=kpoints_table,pos=1,val=auxInt )
      call aot_get_val( L=conf,errCode=iErr,thandle=kpoints_table,pos=2,val=auxInt2 )
      call aot_table_close( L=conf,thandle=kpoints_table )
      this%kListXY(i,1:2)=[auxInt,auxInt2]
      this%irrepListXY(i)=trim(auxString)
      this%parityListXY(i)=trim(auxString2)
      call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=4,val=auxString )
      call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=5,val=auxString2 )
      call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=6,val=auxInt )
      this%kListAngle(i)=auxInt
      this%irrepListAngle(i)=trim(auxString)
      this%parityListAngle(i)=trim(auxString2)
      call aot_table_close( L=conf,thandle=term_table )
   enddo
   call aot_table_close( L=conf,thandle=fourier_table )
   ! ENDDDDDD!!!!!!
   CALL CLOSE_CONFIG(conf)
   RETURN
END SUBROUTINE READ_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_PHI_CRP6D #######################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0. Initial PHI value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 09/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_PHI_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER,INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin 
   REAL(KIND=8),DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=18),PARAMETER :: routinename = "PLOT1D_PHI_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_PHI_CRP6D ERR: Less than 2 points"
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
END SUBROUTINE PLOT1D_PHI_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_ATOMIC_INTERAC_PHI_CRP6D #######################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the atomic interaction
!! PES
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0. Initial PHI value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_ATOMIC_INTERAC_PHI_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
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
   CHARACTER(LEN=18),PARAMETER :: routinename = "PLOT1D_PHI_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_PHI_CRP6D ERR: Less than 2 points"
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
END SUBROUTINE PLOT1D_ATOMIC_INTERAC_PHI_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_PHI_SMOOTH_CRP6D #######################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the smooth PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0. Initial PHI value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 09/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_PHI_SMOOTH_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER,INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin 
   REAL(KIND=8),DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=25),PARAMETER :: routinename = "PLOT1D_PHI_SMOOTH_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_PHI_CRP6D ERR: Less than 2 points"
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
END SUBROUTINE PLOT1D_PHI_SMOOTH_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_THETA_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial THETA value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_THETA_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin 
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=20),PARAMETER :: routinename = "PLOT1D_THETA_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_THETA_CRP6D ERR: Less than 2 points"
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
END SUBROUTINE PLOT1D_THETA_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_ATOMIC_INTERAC_THETA_CRP6D ########################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the atomic PES 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial THETA value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_ATOMIC_INTERAC_THETA_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
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
   CHARACTER(LEN=20),PARAMETER :: routinename = "PLOT1D_THETA_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_THETA_CRP6D ERR: Less than 2 points"
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
END SUBROUTINE PLOT1D_ATOMIC_INTERAC_THETA_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_THETA_SMOOTH_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the smooth PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial THETA value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_THETA_SMOOTH_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=27),PARAMETER :: routinename = "PLOT1D_THETA_SMOOTH_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_THETA_CRP6D ERR: Less than 2 points"
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
END SUBROUTINE PLOT1D_THETA_SMOOTH_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_R_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial R value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_R_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=16),PARAMETER :: routinename = "PLOT1D_R_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_R_CRP6D ERR: Less than 2 points"
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
END SUBROUTINE PLOT1D_R_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_R_SMOOTH_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the smooth PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial R value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_R_SMOOTH_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=23),PARAMETER :: routinename = "PLOT1D_R_SMOOTH_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_R_CRP6D ERR: Less than 2 points"
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
END SUBROUTINE PLOT1D_R_SMOOTH_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_Z_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!> @param[in] L - Length to plot from minimim Z
!
!> @warning
!! - The graph starts always at 0,0, Initial Z value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_Z_CRP6D(thispes,npoints,X,L,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
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
   CHARACTER(LEN=16),PARAMETER :: routinename = "PLOT1D_Z_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_Z_CRP6D ERR: Less than 2 points"
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
END SUBROUTINE PLOT1D_Z_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_Z_SMOOTH_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the smooth PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial Z value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_Z_SMOOTH_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=23),PARAMETER :: routinename = "PLOT1D_Z_SMOOTH_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   SELECT CASE(npoints)
      CASE(: 1)
         WRITE(0,*) "PLOT1D_Z_CRP6D ERR: Less than 2 points"
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
END SUBROUTINE PLOT1D_Z_SMOOTH_CRP6D
!#######################################################################
! SUBROUTINE: PLOT_XYMAP_CRP6D
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (X,Y) of the PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the file to print the output
!> @param[in] init_point - Initial position to start the scan (a.u. and radians). 6 DIMENSIONAL
!> @param[in] nxpoints - Number of points in X axis (auxiliar cartesian coordinates)
!> @param[in] nypoints - Number of points in Y axis (auxiliar cartesian coordinates)
!> @param[in] Lx - Length of X axis (a.u.)
!> @param[in] Ly - Length of Y axis (a.u.)
!
!> @warning
!! - Z,R,THETA,PHI parameters are taken from @b init_point
!
!> @author A.S. Muzas
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT_XYMAP_CRP6D(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(CRP6D),INTENT(IN) :: thispes
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
END SUBROUTINE PLOT_XYMAP_CRP6D
!######################################################################################
! SUBROUTINE: PLOT_XYMAP_SMOOTH_CRP6D
!######################################################################################
!> @brief
!! Same as plot_xymap_crp6d but calling get_v_and_derivs_smooth, i.e. potential
!! and derivatives are not corrected. Thus, we get smooth corrugationless potential.
!-------------------------------------------------------------------------------------
SUBROUTINE PLOT_XYMAP_SMOOTH_CRP6D(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(CRP6D),INTENT(IN) :: thispes
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
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
   END DO
   r(2) = ymax
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   ! inpoints in XY
   DO i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      r(3:6) = init_point(3:6)
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
      DO j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
         WRITE(11,*) r(1:2),v,dvdu(:)
      END DO
      r(2) = ymax
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
   END DO
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   r(3:6) = init_point(3:6)
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
   END DO
   r(2) = ymax
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_XYMAP_SMOOTH_CRP6D
!#######################################################################
! SUBROUTINE: PLOT_RZMAP_CRP6D
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (R,Z) of the PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the file to print the output
!> @param[in] init_point - Initial position to start the scan (a.u. and radians). 6 DIMENSIONAL
!> @param[in] nxpoints - Number of points in R axis
!> @param[in] nypoints - Number of points in Z axis
!> @param[in] Lx - Length of R axis (a.u.)
!> @param[in] Ly - Length of Z axis (a.u.)
!
!> @warning
!! - X,Y,THETA,PHI parameters are taken from @b init_point
!
!> @author A.S. Muzas
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT_RZMAP_CRP6D(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(CRP6D),INTENT(IN) :: thispes
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
END SUBROUTINE PLOT_RZMAP_CRP6D
!#######################################################################
! SUBROUTINE: PLOT_ATOMIC_INTERAC_RZ_CRP6D
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (R,Z) of the sub of atomic potentials. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the file to print the output
!> @param[in] init_point - Initial position to start the scan (a.u. and radians). 6 DIMENSIONAL
!> @param[in] nxpoints - Number of points in R axis
!> @param[in] nypoints - Number of points in Z axis
!> @param[in] Lx - Length of R axis (a.u.)
!> @param[in] Ly - Length of Z axis (a.u.)
!
!> @warning
!! - X,Y,THETA,PHI parameters are taken from @b init_point
!
!> @author A.S. Muzas
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT_ATOMIC_INTERAC_RZ_CRP6D(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(CRP6D),INTENT(IN) :: thispes
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
END SUBROUTINE PLOT_ATOMIC_INTERAC_RZ_CRP6D
!###########################################################
!# FUNCTION: is_allowed_CRP6D 
!###########################################################
!> @brief 
!! Enquires if the CRP6D potential can be calculated at @b X
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
LOGICAL FUNCTION is_allowed_CRP6D(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8) :: zmin, rmin, rmax
   ! Run section
   zmin=this%wyckoffsite(1)%zrcut(1)%getfirstZ()
   rmin=this%wyckoffsite(1)%zrcut(1)%getfirstR()
   rmax=this%wyckoffsite(1)%zrcut(1)%getlastR()
   SELECT CASE(size(x)/=6)
      CASE(.TRUE.)
         WRITE(0,*) "is_allowed_CRP6D ERR: wrong number of dimensions"
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(x(3)<zmin .OR. x(4)<rmin .OR. x(4)>rmax)
      CASE(.TRUE.)
         is_allowed_CRP6D=.FALSE.
      CASE(.FALSE.)
         is_allowed_CRP6D=.TRUE.
   END SELECT
   RETURN
END FUNCTION is_allowed_CRP6D

END MODULE CRP6D_MOD
