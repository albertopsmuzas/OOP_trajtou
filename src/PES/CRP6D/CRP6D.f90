!#######################################################
! MODULE: CRP6D_MOD
!#######################################################
!> @brief
!! Provides all structures and procedures to perform a complete
!! CRP-6D interpolation job. This method only works systematically with diatomic
!! molecules
!#######################################################
module CRP6D_MOD
! Massive name spaces
use SYSTEM_MOD
use LINK_FUNCTION1D_MOD
! Selective name spaces
use UNITS_MOD,              only: Length, pi
use PES_MOD,                only: PES
use CRP3D_MOD,              only: CRP3D
use EXTRAPOL_TO_VACUUM_MOD, only: Vacuumpot
use FOURIER_P4MM_MOD,       only: Fourierp4mm
use FOURIER3D_P4MM_MOD,     only: Fourier3d_p4mm
use WYCKOFF_P4MM_MOD,       only: WyckoffSitio,Wyckoffp4mm,TermsInfo
use FOURIER3D_P4MM_MOD,     only: Fourier3d_p4mm
use AOTUS_MODULE,           only: flu_State, OPEN_CONFIG_FILE, CLOSE_CONFIG, AOT_GET_VAL
use AOT_TABLE_MODULE,       only: AOT_TABLE_OPEN, AOT_TABLE_CLOSE, AOT_TABLE_LENGTH, AOT_TABLE_GET_VAL
#if DEBUG
use DEBUG_MOD, only: verbose_write, debug_write
#endif
implicit none
!/////////////////////////////////////////////////
! TYPE: CRP6D
!
!-------------------------------------------------
type,extends(PES):: CRP6D
   integer(kind=4):: nsites
   integer(kind=4):: natomic
   logical:: is_interpolated=.false.
   logical:: is_homonucl=.false.
   logical:: is_smooth=.false.
   logical:: is_shifted=.false.
   logical:: is_resized=.false.
   real(kind=8):: zvacuum
   integer(kind=4),dimension(2):: grid=[0,0]
   class(Wyckoffsitio),dimension(:),allocatable:: wyckoffsite
   type(CRP3D),dimension(:),allocatable:: atomiccrp
   type(Vacuumpot):: farpot
   class(Function1d),allocatable:: dumpfunc
   character(len=30):: extrapol2vac_flag
   integer(kind=4),dimension(:,:),allocatable:: kListXY
   character(len=1),dimension(:),allocatable:: parityListXY
   character(len=2),dimension(:),allocatable:: irrepListXY
   integer(kind=4),dimension(:),allocatable:: kListAngle
   character(len=1),dimension(:),allocatable:: parityListAngle
   character(len=2),dimension(:),allocatable:: irrepListAngle
   integer(kind=4):: totalTerms=0

   contains
      ! Initialization block
      procedure,public:: READ => READ_CRP6D
      procedure,public:: INITIALIZE => INITIALIZE_CRP6D
      ! Set block
      procedure,public:: SET_SMOOTH => SET_SMOOTH_CRP6D
      ! Get block
      procedure,public:: GET_V_AND_DERIVS => GET_V_AND_DERIVS_CRP6D
      procedure,public:: GET_V_AND_DERIVS_PURE => GET_V_AND_DERIVS_PURE_CRP6D
      procedure,public:: GET_V_AND_DERIVS_SMOOTH => GET_V_AND_DERIVS_SMOOTH_CRP6D
      procedure,public:: GET_ATOMICPOT_AND_DERIVS => GET_ATOMICPOT_AND_DERIVS_CRP6D
      ! Tools block
      procedure,public:: SMOOTH => SMOOTH_CRP6D
      procedure,public:: SMOOTH_EXTRA => SMOOTH_EXTRA_CRP6D
      procedure,public:: INTERPOL => INTERPOL_CRP6D
      procedure,public:: RAWINTERPOL => RAWINTERPOL_CRP6D
      procedure,public:: EXTRACT_VACUUMSURF => EXTRACT_VACUUMSURF_CRP6D
      procedure,public:: ADD_VACUUMSURF => ADD_VACUUMSURF_CRP6D
      procedure,public:: INTERPOL_NEW_RZGRID => INTERPOL_NEW_RZGRID_CRP6D
      procedure,public:: CHEAT_CARTWHEEL_ONTOP => CHEAT_CARTWHEEL_ONTOP_CRP6D
      ! Plot toolk
      procedure,public:: PLOT1D_THETA => PLOT1D_THETA_CRP6D
      procedure,public:: PLOT1D_THETA_SMOOTH => PLOT1D_THETA_SMOOTH_CRP6D
      procedure,public:: PLOT1D_ATOMIC_INTERAC_THETA => PLOT1D_ATOMIC_INTERAC_THETA_CRP6D
      procedure,public:: PLOT1D_PHI => PLOT1D_PHI_CRP6D
      procedure,public:: PLOT1D_PHI_SMOOTH => PLOT1D_PHI_SMOOTH_CRP6D
      procedure,public:: PLOT1D_ATOMIC_INTERAC_PHI => PLOT1D_ATOMIC_INTERAC_PHI_CRP6D
      procedure,public:: PLOT1D_R => PLOT1D_R_CRP6D
      procedure,public:: PLOT1D_R_SMOOTH => PLOT1D_R_SMOOTH_CRP6D
      procedure,public:: PLOT1D_Z => PLOT1D_Z_CRP6D
      procedure,public:: PLOT1D_Z_SMOOTH => PLOT1D_Z_SMOOTH_CRP6D
      procedure,public:: PLOT_XYMAP => PLOT_XYMAP_CRP6D
      procedure,public:: PLOT_XYMAP_SMOOTH => PLOT_XYMAP_SMOOTH_CRP6D
      procedure,public:: PLOT_RZMAP => PLOT_RZMAP_CRP6D
      procedure,public:: PLOT_ATOMIC_INTERAC_RZ => PLOT_ATOMIC_INTERAC_RZ_CRP6D
      ! Enquire block
      procedure,public:: is_allowed => is_allowed_CRP6D
end type CRP6D
contains
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
subroutine INITIALIZE_CRP6D(this,filename,tablename)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP6D),intent(out)::this
   character(len=*),optional,intent(in):: filename,tablename
   ! Local variables
   character(len=:),allocatable:: auxstring
   ! Run section
   select case(allocated(system_inputfile) .or. .not.present(filename))
      case(.true.)
         auxstring=system_inputfile
      case(.false.)
         auxstring=filename
   end select
   select case(present(tablename))
      case(.true.)
         call this%READ(filename=auxstring,tablename=tablename)
      case(.false.)
         call this%READ(filename=auxstring,tablename='pes')
   end select
   call this%INTERPOL()
   select case(this%is_resized)
      case(.true.)
         call this%INTERPOL_NEW_RZGRID(this%grid(1),this%grid(2))
      case(.false.)
         ! do nothing
   end select
   return
end subroutine INITIALIZE_CRP6D
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
subroutine GET_ATOMICPOT_AND_DERIVS_CRP6D(this,molecx,atomicx,v,dvdu)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP6D),intent(in):: this
   real(kind=8),dimension(6),intent(in):: molecx
   real(kind=8),dimension(2),intent(out):: v
   real(kind=8),dimension(6),intent(out):: dvdu
   real(kind=8),dimension(6),intent(out):: atomicx
   ! Local variables
   real(kind=8):: vcorra,vcorrb
   ! Run section
   atomicx(:)=from_molecular_to_atomic(molecx)
   select case(this%natomic)
      case(1)
         call this%atomiccrp(1)%GET_V_AND_DERIVS(atomicx(1:3),v(1),dvdu(1:3))
         call this%atomiccrp(1)%GET_V_AND_DERIVS(atomicx(4:6),v(2),dvdu(4:6))
      case(2)
         call this%atomiccrp(1)%GET_V_AND_DERIVS(atomicx(1:3),v(1),dvdu(1:3))
         call this%atomiccrp(2)%GET_V_AND_DERIVS(atomicx(4:6),v(2),dvdu(4:6))
      case default
         write(0,*) "GET_ATOMICPOT_AND_DERIVS_CRP6D ERR: wrong number of atomic potentials"
         call EXIT(1)
   end select
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
   return
end subroutine GET_ATOMICPOT_AND_DERIVS_CRP6D
!###########################################################
!# SUBROUTINE: SET_SMOOTH_CRP6D
!###########################################################
!> @brief
!! Common set function. Sets is_smooth atribute
!-----------------------------------------------------------
subroutine SET_SMOOTH_CRP6D(this,is_smooth)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP6D),intent(inout):: this
   logical,intent(in) :: is_smooth
   ! Run section
   this%is_smooth=is_smooth
   return
end subroutine SET_SMOOTH_CRP6D
!###########################################################
!# SUBROUTINE: CHEAT_CARTWHEEL_ONTOP_CRP6D 
!###########################################################
!> @brief
!! This subroutine changes data from a raw cut2d input by the
!! value interpolated from atomic potentials. Only data close
!! to the plane @f$ \pi: z=\over{r}{2}+rumpling@f$
!
!> param[in] wyckoff - Wyckoff ID site
!> param[in] cut2d   - Cut2d ID (depends on the chosen Wyckoff ID)
!> param[in] topTye  - ID of top site chosen. This value is used
!!                     to correctly apply a given rumpling value.
!!                     This ID depends on which order the repulsive
!!                     crp3d potentials were declared in the input.
!> param[in] dmax    - Points are only modified if they are closer
!!                     to the plane @f$ \pi: z=\over{r}{2}+rumpling@f$
!!                     than this value. Should be given in bohr radius.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jul/2016
!> @version 2.0
!-----------------------------------------------------------
subroutine CHEAT_CARTWHEEL_ONTOP_CRP6D(this,wyckoff,cut2d,toptype,dmax)
   ! Initial declarations   
   implicit none
   ! I/O 
   class(CRP6D),intent(inout) :: this
   integer(kind=4),intent(in) :: wyckoff,cut2d,toptype
   real(kind=8),intent(in) :: dmax
   ! Local variables
   integer(kind=4) :: i,j ! counters
   real(kind=8) :: r,z,rump,interac
   real(kind=8) :: x,y ! x,y position
   real(kind=8) :: m1, m2 ! masses
   real(kind=8) :: dist ! distance to the plane
   integer(kind=4) :: nr,nz
   real(kind=8),dimension(3) :: pos1,pos2
   ! Run section
   nr=this%wyckoffsite(wyckoff)%zrcut(cut2d)%getgridsizeR()
   nz=this%wyckoffsite(wyckoff)%zrcut(cut2d)%getgridsizeZ()
   x=this%wyckoffsite(wyckoff)%x
   y=this%wyckoffsite(wyckoff)%y
   m1=system_mass(1)
   m2=system_mass(2)
   rump=this%atomiccrp(1)%getrumpling(toptype)
   do i = 1, nr
      do j = 1, nz
        r=this%wyckoffsite(wyckoff)%zrcut(cut2d)%getgridvalueR(i)
        z=this%wyckoffsite(wyckoff)%zrcut(cut2d)%getgridvalueZ(j)
        dist=(2.D0/sqrt(5.D0))*(z+rump-r/2.D0)
        pos1=(/x,y,z+(m2/(m1+m2))*r/)
        pos2=(/x,y,z-(m1/(m1+m2))*r/)
        select case(dist <= dmax) 
           case(.true.)
              select case(this%is_homonucl)
                 case(.true.)
                    interac=this%atomiccrp(1)%getpot(pos1)+this%atomiccrp(1)%getpot(pos2)!+this%farpot%getpot(r)
                 case(.false.)
                    interac=this%atomiccrp(1)%getpot(pos1)+this%atomiccrp(2)%getpot(pos2)!+this%farpot%getpot(r)
              end select
              call this%wyckoffsite(wyckoff)%zrcut(cut2d)%CHANGEPOT_AT_GRIDPOINT(i,j,interac)
           case(.false.)
              ! do nothing
        end select
      end do
   end do
   return
end subroutine CHEAT_CARTWHEEL_ONTOP_CRP6D
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
subroutine INTERPOL_NEW_RZGRID_CRP6D(this,nRpoints,nZpoints)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP6D),intent(inout) :: this
   integer(kind=4),intent(in) :: nrpoints,nzpoints
   ! Local variables
   integer(kind=4) :: i,j ! counters
   ! Run section
   do i = 1, this%nsites
      do j = 1, this%wyckoffsite(i)%n2dcuts
         call this%wyckoffsite(i)%zrcut(j)%interrz%INTERPOL_NEWGRID(nrpoints,nzpoints)
      end do
   end do
   return
end subroutine INTERPOL_NEW_RZGRID_CRP6D
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
subroutine GET_V_AND_DERIVS_PURE_CRP6D(this,x,v,dvdu)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP6D),target,intent(in):: this
   real(kind=8),dimension(:),intent(in):: x
   real(kind=8),intent(out):: v
   real(kind=8),dimension(:),intent(out):: dvdu
   ! Local variables
   integer(kind=4):: i,j,h ! counters
   real(kind=8):: ma,mb
   real(kind=8),dimension(:,:),allocatable:: f ! smooth function and derivs
   real(kind=8),dimension(:,:),allocatable:: xyList
   real(kind=8),dimension(:),allocatable:: phiList
   real(kind=8),dimension(:,:),allocatable:: completeGeomList
   real(kind=8),dimension(:),allocatable:: aux1
   real(kind=8),dimension(:,:),allocatable:: aux2
   real(kind=8),dimension(6):: atomicx
   real(kind=8),dimension(2):: atomic_v
   real(kind=8),dimension(6):: atomic_dvdu
   real(kind=8),dimension(3):: dvdu_atomicA,dvdu_atomicB
   type(Fourier3d_p4mm):: fouInterpol
   character(len=*),parameter:: routinename="GET_V_AND_DERIVS_CRP6D: "
   ! Run section
   select case(this%is_smooth)
      case(.true.)
         ! do nothing
      case(.false.)
         write(0,*) "GET_V_AND_DERIVS_CRP6D ERR: smooth the PES first (CALL thispes%SMOOTH())"
         call EXIT(1)
   end select
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
   call fouInterpol%interpol()
   allocate(aux1(4))
   allocate(aux2(4,3))
   call fouInterpol%get_allFuncs_and_derivs( x=[x(1),x(2),x(6)],f=aux1,dfdx=aux2 )
   call fouInterpol%cleanAll()
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
   call DEBUG_WRITE(routinename,"Smooth V at each thetaCut:")
   call DEBUG_WRITE(routinename,f(1,:))
   call DEBUG_WRITE(routinename,"Smooth dv/dz at each thetaCut:")
   call DEBUG_WRITE(routinename,f(2,:))
   call DEBUG_WRITE(routinename,"Smooth dv/dr at each thetaCut:")
   call DEBUG_WRITE(routinename,f(3,:))
   call DEBUG_WRITE(routinename,"Smooth dv/dtheta at each thetaCut:")
   call DEBUG_WRITE(routinename,f(4,:))
   call DEBUG_WRITE(routinename,"Smooth interpolated values:")
   call DEBUG_WRITE(routinename,"v: ",aux1(1))   ! v
   call DEBUG_WRITE(routinename,"dvdx: ",aux2(1,1)) ! dvdx
   call DEBUG_WRITE(routinename,"dvdy: ",aux2(1,2)) ! dvdy
   call DEBUG_WRITE(routinename,"dvdz: ",aux1(2))   ! dvdz
   call DEBUG_WRITE(routinename,"dvdr: ",aux1(3))   ! dvdr
   call DEBUG_WRITE(routinename,"dvdtheta: ",aux1(4))   ! dvdtheta
   call DEBUG_WRITE(routinename,"dvdphi: ",aux2(1,3))   ! dvdphi
#endif
   !--------------------------------------
   ! Results for the real potential
   !-------------------------------------
   call this%GET_ATOMICPOT_AND_DERIVS(x,atomicx,atomic_v,atomic_dvdu)
   dvdu_atomicA=atomic_dvdu(1:3)
   dvdu_atomicB=atomic_dvdu(4:6)
#ifdef DEBUG
   call DEBUG_WRITE(routinename,"Contributions of the atomic potential: ")
   call DEBUG_WRITE(routinename, "Position Atom A: ",atomicx(1:3))
   call DEBUG_WRITE(routinename,"Va: ",atomic_v(1))
   call DEBUG_WRITE(routinename,"dVa/dxa; dVa/dya: ",dvdu_atomicA)
   call DEBUG_WRITE(routinename, "Position Atom B: ",atomicx(4:6))
   call DEBUG_WRITE(routinename,"Vb: ",atomic_v(2))
   call DEBUG_WRITE(routinename,"dVa/dxa; dVb/dyb: ",dvdu_atomicB)
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
   return
end subroutine GET_V_AND_DERIVS_PURE_CRP6D
!###########################################################
!# SUBROUTINE: GET_V AND DERIVS_CRP6D 
!###########################################################
!> @brief
!! Gets potential for any configuration. Discriminates between 3 regions:
!! - Pure CRP6D region: zmin to zcrp
!! - Extrapolation region zcrp to zvacuum
!! - Vacuum region: further than zvacuum
!
!> @param[in]  this     CRP6D potential used
!> @param[in]  x        Geometry array (x,y,z,r,theta,phi)
!> @param[out] v        Potential energy in hartrees
!> @param[out] dvdu     Derivatives array (dvdx,dvdy,dvdz,dvdr,dvdtheta,dvdphi) @b x geometry
!> @param[out] errCode  Optional error code to know the status of energy evaluation

!> @details
!! ERROR CODE GENERATION
!! errCode ==  0  => Everithing is OK
!! errCode == -1  => r > rCRPmax. Then, we used r=rCRPmax
!! errCode == -2  => r < rCRPmin. Then, we used r=rCRPmin
!! errCode == -3  => z < zCRPmin. Then, we used z=zCRPmin
!! errCode == -10 => The geometry was unclassificable (beware of this problem >_< )
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014, Jan/2014
!> @version 1.2
!-----------------------------------------------------------
subroutine GET_V_AND_DERIVS_CRP6D(this,X,v,dvdu,errCode)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP6D),target,intent(in):: this
   real(kind=8),dimension(:),intent(in):: x
   real(kind=8),intent(out):: v
   real(kind=8),dimension(:),intent(out):: dvdu
   integer(kind=1),optional,intent(out):: errCode
   ! Local variables
   real(kind=8),dimension(6):: geom             !< @var auxiliar geometry that may change respect to @b x
   real(kind=8):: zCrpMin                       !< @var Minimum z value at which the potential is defined
   real(kind=8):: zCrpMax                       !< @var Maximum z value at which the potential is defined
   real(kind=8):: zVac                          !< @var Z infinity
   real(kind=8):: rCrpMin                       !< @var Minimum r value at which the potential is defined
   real(kind=8):: rCrpMax                       !< @var Maximum r value at which the potential is defined
   real(kind=8):: vzCrpMax                      !< @var Potential at (x,y,zCrpMax,r,theta,phi)
   real(kind=8):: vzvac                         !< @var Potential at (x,y,zVac,r,theta,phi)
   real(kind=8),dimension(6):: dvducrp          !< @var Derivatives at zcrp
   real(kind=8),dimension(6):: dvduvac          !< @var Derivatives at vacuum
   real(kind=8):: alpha,beta,gama               !< @var Internal parameters
   class(Function1d),allocatable:: extrapolfunc !< @var Function used to extrapolate potential between zcrp and zvac
   integer:: error                              !< @var Internal error code that may be inherit by @b errCode if present
   integer(kind=4):: i !counter
   ! Local Parameter
   real(kind=8),parameter:: zero=0.D0 ! what we will consider zero (a.u.)
   real(kind=8),parameter:: dz=0.5D0  ! 0.25 Angstroems approx
   character(len=*),parameter:: routinename='GET_V_AND_DERIVS_CRP6D: '
   ! Run section
   zCrpMin=this%wyckoffSite(1)%zrcut(1)%getFirstZ()
   zCrpMax=this%wyckoffSite(1)%zrcut(1)%getLastZ()
   zvac=this%zvacuum
   rCrpMin=this%wyckoffSite(1)%zrcut(1)%getFirstR()
   rCrpMax=this%wyckoffSite(1)%zrcut(1)%getLastR()
   geom(:)=x(1:6)
   ! Check if the energy evaluation  is outside of the boundaries defined by PES files. Choose an error code when appropriate.
   error=0
   select case( geom(4)>rCrpMax )
   case(.true.)
      geom(4)=rCrpMax
      error=-1
   case(.false.)
      ! do nothing
   end select

   select case( geom(4)<rCrpMin )
   case(.true.)
      geom(4)=rCrpMin
      error=-2
   case(.false.)
      ! do nothing
   end select

   select case( geom(3)<zCrpMin )
   case(.true.)
      geom(3)=zCrpMin
      error=-3
   case(.false.)
      ! do nothing
   end select

   ! *************************************************************************
   ! SWITCH 1: Check if we are inside Z's range
   select case( geom(3) >= zCrpMin .and. geom(3)<= zCrpMax ) !easy
   case(.true.)
      call this%get_v_and_derivs_pure(geom,v,dvdu)
      return
   case(.false.)
      ! do nothing, next switch
   end select

   ! *************************************************************************
   ! SWITCH 2: Check if we are inside the extrapolation region
   select case(geom(3)>zCrpMax .and. geom(3)<zvac)
      case(.true.) ! uff
         ! Set potential and derivs
         vzvac=this%farpot%getpot( geom(4) )
         call this%GET_V_AND_DERIVS_PURE([geom(1),geom(2),zCrpMax,geom(4),geom(5),geom(6)],vzCrpMax,dvducrp)
         dvduvac(1:3)=zero
         dvduvac(4)=this%farpot%getderiv(geom(4))
         dvduvac(5:6)=zero
         ! Check kind of extrapolation
         select case(this%extrapol2vac_flag)
            case("Xexponential")
               allocate(Xexponential_func::extrapolfunc)
               ! Extrapol potential
               beta=-1.D0/zvac
               alpha=(vzCrpMax-vzvac)/(zCrpMax*dexp(beta*zCrpMax)-zvac*dexp(beta*zvac))
               gama=vzvac-alpha*zvac*dexp(beta*zvac)
               call extrapolfunc%READ([alpha,beta])
               v=extrapolfunc%getvalue(geom(3))+gama
               dvdu(3)=extrapolfunc%getderiv(geom(3))
               ! Extrapol derivatives
               do i = 1, 6
                  select case(i)
                     case(3)
                        ! Skip dvdz
                     case default
                        beta=-1.D0/zvac
                        alpha=(dvducrp(i)-dvduvac(i))/(zCrpMax*dexp(beta*zCrpMax)-zvac*dexp(beta*zvac))
                        gama=dvduvac(i)-alpha*zvac*dexp(beta*zvac)
                        call extrapolfunc%READ([alpha,beta])
                        dvdu(i)=extrapolfunc%getvalue(geom(3))+gama
                  end select
               end do
               deallocate(extrapolfunc)
               return

            case("Linear")
               allocate(Linear_func::extrapolfunc)
               ! Extrapol potential
               beta=(vzCrpMax*zvac-vzvac*zCrpMax)/(zvac-zCrpMax)
               alpha=(vzvac-beta)/zvac
               call extrapolfunc%READ([alpha,beta])
               v=extrapolfunc%getvalue(geom(3))
               dvdu(3)=extrapolfunc%getderiv(geom(3))
               ! Extrapol derivatives
               do i = 1, 6
                  select case(i)
                     case(3)
                        ! do nothing
                     case default
                        beta=(dvducrp(i)*zvac-dvduvac(i)*zcrpmax)/(zvac-zcrpmax)
                        alpha=(dvduvac(i)-beta)/zvac
                        call extrapolfunc%READ([alpha,beta])
                        dvdu(i)=extrapolfunc%getvalue(geom(3))
                  end select
               end do
               deallocate(extrapolfunc)
               return

            case default
               write(0,*) "GET_V_AND_DERIVS_CRP6D ERR: type of extrapolation function isn't implemented yet"
               write(0,*) "Implemented ones: Linear, Xexponential, None"
               call EXIT(1)
         end select
      case(.false.)
         ! do nothing
   end select
   ! *************************************************************************
   ! SWITCH 3: Check if we are inside the Vacuum region
   select case(geom(3)>=zvac) !easy
      case(.true.)
         v=this%farpot%getpot(geom(4))
         dvdu(1:3)=0.d0
         dvdu(4)=this%farpot%getderiv(geom(4))
         dvdu(5:6)=0.d0
      case(.false.) ! this's the last switch!
         error=-10
#ifdef DEBUG
         call debug_write(routinename,'Unclassificable point')
         call debug_write(routinename,'Asking potential at Z: ',geom(3))
#endif
   end select
   ! Set error Code if present
   select case(present(errCode))
   case(.true.)
      errCode = error
   case(.false.)
      ! do nothing
   end select
   return
end subroutine GET_V_AND_DERIVS_CRP6D
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
subroutine GET_V_AND_DERIVS_SMOOTH_CRP6D(this,x,v,dvdu)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP6D),target,intent(in):: this
   real(kind=8),dimension(:),intent(in):: x
   real(kind=8),intent(out):: v
   real(kind=8),dimension(:),intent(out):: dvdu
   ! Local variables
   integer(kind=4):: i,j,h ! counters
   real(kind=8),dimension(:,:),allocatable:: f ! smooth function and derivs
   real(kind=8),dimension(:,:),allocatable:: xyList
   real(kind=8),dimension(:),allocatable:: phiList
   real(kind=8),dimension(:,:),allocatable:: completeGeomList
   real(kind=8),dimension(:),allocatable:: aux1
   real(kind=8),dimension(:,:),allocatable:: aux2
   type(Fourier3d_p4mm):: fouInterpol
   character(len=*),parameter:: routinename="GET_V_AND_DERIVS_CRP6D: "
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
   call fouInterpol%interpol()
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
   call DEBUG_WRITE(routinename,"Smooth V at each thetaCut:")
   call DEBUG_WRITE(routinename,f(1,:))
   call DEBUG_WRITE(routinename,"Smooth dv/dz at each thetaCut:")
   call DEBUG_WRITE(routinename,f(2,:))
   call DEBUG_WRITE(routinename,"Smooth dv/dr at each thetaCut:")
   call DEBUG_WRITE(routinename,f(3,:))
   call DEBUG_WRITE(routinename,"Smooth dv/dtheta at each thetaCut:")
   call DEBUG_WRITE(routinename,f(4,:))
   call DEBUG_WRITE(routinename,"Smooth interpolated values:")
   call DEBUG_WRITE(routinename,"v: ",aux1(1))   ! v
   call DEBUG_WRITE(routinename,"dvdx: ",aux2(1,1)) ! dvdx
   call DEBUG_WRITE(routinename,"dvdy: ",aux2(1,2)) ! dvdy
   call DEBUG_WRITE(routinename,"dvdz: ",aux1(2))   ! dvdz
   call DEBUG_WRITE(routinename,"dvdr: ",aux1(3))   ! dvdr
   call DEBUG_WRITE(routinename,"dvdtheta: ",aux1(4))   ! dvdtheta
   call DEBUG_WRITE(routinename,"dvdphi: ",aux2(1,3))   ! dvdphi
#endif
   v=aux1(1)
   dvdu(1)=aux2(1,1)
   dvdu(2)=aux2(1,2)
   dvdu(3)=aux1(2)
   dvdu(4)=aux1(3)
   dvdu(5)=aux1(4)
   dvdu(6)=aux2(1,3)
   return
end subroutine GET_V_AND_DERIVS_SMOOTH_CRP6D
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
subroutine INTERPOL_CRP6D(this)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP6D),intent(inout)::this
   ! Local variables
   integer(kind=4) :: i,j ! counters
   ! Run section
   call this%EXTRACT_VACUUMSURF()   
   call this%SMOOTH()
   !CALL this%SMOOTH_EXTRA()
   do i = 1, this%nsites
      do j = 1, this%wyckoffsite(i)%n2dcuts
         call this%wyckoffsite(i)%zrcut(j)%INTERPOL()
      end do
   end do
   this%is_interpolated=.true.
   return
end subroutine INTERPOL_CRP6D
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
subroutine RAWINTERPOL_CRP6D(this)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP6D),intent(inout)::this
   ! Local variables
   integer(kind=4) :: i,j ! counters
   ! Run section
   call this%EXTRACT_VACUUMSURF()   
   do i = 1, this%nsites
      do j = 1, this%wyckoffsite(i)%n2dcuts
         call this%wyckoffsite(i)%zrcut(j)%INTERPOL()
      end do
   end do
   this%is_interpolated=.true.
   return
end subroutine RAWINTERPOL_CRP6D
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
subroutine SMOOTH_CRP6D(this)
   ! Initial declarations   
   implicit none
   ! I/O variables
    class(CRP6D),intent(inout):: this
   ! Local variables
   real(kind=8),dimension(6):: molcoord,atomcoord,dummy
   integer(kind=4):: nr,nz
   integer(kind=4):: i,j,k,l ! counters
   character(len=*),parameter:: routinename="SMOOTH_CRP6D: "
   real(kind=8):: newpot
   real(kind=8),dimension(2):: atomic_v
   ! Run section
   do i = 1, this%nsites ! cycle wyckoff sites
      molcoord(1)=this%wyckoffsite(i)%x
      molcoord(2)=this%wyckoffsite(i)%y
      do j = 1, this%wyckoffsite(i)%n2dcuts
         molcoord(5)=this%wyckoffsite(i)%zrcut(j)%theta
         molcoord(6)=this%wyckoffsite(i)%zrcut(j)%phi
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizeR()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizeZ()
         do k = 1, nr
            do l = 1, nz
               molcoord(3)=this%wyckoffsite(i)%zrcut(j)%getgridvalueZ(l)
               molcoord(4)=this%wyckoffsite(i)%zrcut(j)%getgridvalueR(k)
               call this%GET_ATOMICPOT_AND_DERIVS(molcoord,atomcoord,atomic_v,dummy)
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)-sum(atomic_v)
               call this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
            end do
         end do
      end do
   end do
   this%is_smooth=.true.
   return
end subroutine SMOOTH_CRP6D
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
subroutine ROUGH_CRP6D(this)
   ! Initial declarations   
   implicit none
   ! I/O variables
    class(CRP6D),intent(inout):: this
   ! Local variables
   real(kind=8),dimension(6):: molcoord,atomcoord
   integer(kind=4):: nr,nz
   integer(kind=4):: i,j,k,l ! counters
   real(kind=8):: newpot
   ! Parameters
   character(len=*),parameter:: routinename="ROUGH_CRP6D: "
   ! Run section
   do i = 1, this%nsites ! cycle wyckoff sites
      do j = 1, this%wyckoffsite(i)%n2dcuts
         molcoord(1)=this%wyckoffsite(i)%zrcut(j)%x
         molcoord(2)=this%wyckoffsite(i)%zrcut(j)%y
         molcoord(5)=this%wyckoffsite(i)%zrcut(j)%theta
         molcoord(6)=this%wyckoffsite(i)%zrcut(j)%phi
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizer()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizez()
         do k = 1, nr
            do l = 1, nz
               molcoord(3)=this%wyckoffsite(i)%zrcut(j)%getgridvalueZ(l)
               molcoord(4)=this%wyckoffsite(i)%zrcut(j)%getgridvalueR(k)
               atomcoord(:)=from_molecular_to_atomic(molcoord)
#ifdef DEBUG
               call DEBUG_WRITE(routinename,"Molecular coords:")
               call DEBUG_WRITE(routinename,molcoord)
               call DEBUG_WRITE(routinename,"Atomic coords:")
               call DEBUG_WRITE(routinename,atomcoord)
#endif
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)
               select case(this%natomic)
                  case(1)
                     newpot=newpot+this%atomiccrp(1)%getpot(atomcoord(1:3))+&
                        this%atomiccrp(1)%getpot(atomcoord(4:6))
                     call this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  case(2)
                     newpot=newpot+this%atomiccrp(1)%getpot(atomcoord(1:3))+&
                        this%atomiccrp(2)%getpot(atomcoord(4:6))
                     call this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  case default
                     write(0,*) "SMOOTH_CRP6D ERR: Something is wrong with number of atomic crp potentials"
                     call EXIT(1)
               end select
            end do
         end do
      end do
   end do
   this%is_smooth=.false.
   return
end subroutine ROUGH_CRP6D
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
subroutine SMOOTH_EXTRA_CRP6D(this)
   ! Initial declarations   
   implicit none
   ! I/O variables
    class(CRP6D),intent(inout) :: this
   ! Local variables
   real(kind=8),dimension(6) :: molcoord,atomcoord
   integer(kind=4) :: nr,nz
   integer(kind=4) :: i,j,k,l ! counters
   character(len=20),parameter :: routinename="SMOOTH_EXTRA_CRP6D: "
   real(kind=8) :: newpot
   ! Run section
   do i = 1, this%nsites ! cycle wyckoff sites
      do j = 1, this%wyckoffsite(i)%n2dcuts
         molcoord(1)=this%wyckoffsite(i)%zrcut(j)%x
         molcoord(2)=this%wyckoffsite(i)%zrcut(j)%y
         molcoord(5)=this%wyckoffsite(i)%zrcut(j)%theta
         molcoord(6)=this%wyckoffsite(i)%zrcut(j)%phi
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizer()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizez()
         do k = 1, nr
            do l = 1, nz
               molcoord(3)=this%wyckoffsite(i)%zrcut(j)%getgridvalueZ(l)
               molcoord(4)=this%wyckoffsite(i)%zrcut(j)%getgridvalueR(k)
               atomcoord(:)=from_molecular_to_atomic(molcoord)
#ifdef DEBUG
               call DEBUG_WRITE(routinename,"Molecular coords:")
               call DEBUG_WRITE(routinename,molcoord)
               call DEBUG_WRITE(routinename,"Atomic coords:")
               call DEBUG_WRITE(routinename,atomcoord)
#endif
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)
               select case(this%natomic)
                  case(1)
                     newpot=newpot-this%atomiccrp(1)%getpot(atomcoord(1:3))&
                        -this%atomiccrp(1)%getpot(atomcoord(4:6))-this%farpot%getpot(molcoord(4))&
                        -0.8*dexp(-molcoord(4))
                     call this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  case(2)
                     newpot=newpot-this%atomiccrp(1)%getpot(atomcoord(1:3))-&
                        this%atomiccrp(2)%getpot(atomcoord(4:6))-this%farpot%getpot(molcoord(4))
                     call this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  case default
                     write(0,*) "SMOOTH_EXTRA_CRP6D ERR: Something is wrong with number of atomic crp potentials"
                     call EXIT(1)
               end select
            end do
         end do
      end do
   end do
   this%is_smooth=.true.
   return
end subroutine SMOOTH_EXTRA_CRP6D
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
subroutine EXTRACT_VACUUMSURF_CRP6D(this)
   ! Initial declarations   
   implicit none
   ! I/O variables
    class(CRP6D),intent(inout) :: this
   ! Local variables
   integer(kind=4) :: nr,nz
   integer(kind=4) :: i,j,k,l ! counters
   real(kind=8) :: newpot
   character(len=26),parameter :: routinename="EXTRACT_VACUUMSURF_CRP6D: "
   ! Run section
   do i = 1, this%nsites ! cycle wyckoff sites
      do j = 1, this%wyckoffsite(i)%n2dcuts
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizeR()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizeZ()
         do k = 1, nr
            do l = 1, nz
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)
               newpot=newpot-this%farpot%getscalefactor()
               call this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
            end do
         end do
      end do
   end do
   call this%farpot%SHIFTPOT()
#ifdef DEBUG
   call VERBOSE_WRITE(routinename,"Potential shifted: ",this%farpot%getscalefactor())
#endif
   this%is_shifted=.true.
   return
end subroutine EXTRACT_VACUUMSURF_CRP6D
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
subroutine ADD_VACUUMSURF_CRP6D(this)
   ! Initial declarations   
   implicit none
   ! I/O variables
    class(CRP6D),intent(inout) :: this
   ! Local variables
   integer(kind=4) :: nr,nz
   integer(kind=4) :: i,j,k,l ! counters
   real(kind=8) :: newpot
   character(len=22),parameter :: routinename="ADD_VACUUMSURF_CRP6D: "
   ! Run section
   do i = 1, this%nsites ! cycle wyckoff sites
      do j = 1, this%wyckoffsite(i)%n2dcuts
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizeR()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizeZ()
         do k = 1, nr
            do l = 1, nz
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)
               newpot=newpot+this%farpot%getscalefactor()
               call this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
            end do
         end do
      end do
   end do
   call this%farpot%SHIFTPOT_UNDO()
#ifdef DEBUG
   call VERBOSE_WRITE(routinename,"Potential re-shifted: ",this%farpot%getscalefactor())
#endif
   this%is_shifted=.false.
   return
end subroutine ADD_VACUUMSURF_CRP6D
!###########################################################
!# SUBROUTINE: READ_CRP6D 
!###########################################################
!> @brief
!! Sets up a CRP6D Object
!-----------------------------------------------------------
subroutine READ_CRP6D(this,filename,tablename)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP6D),intent(out):: this
   character(len=*),intent(in):: filename,tablename
   ! Local variables
   integer(kind=4):: i,j,k,n ! counters
   real(kind=8),dimension(:),allocatable:: param_damp
   ! Wyckoff variables
   character(len=1),dimension(:),allocatable:: wyckoff_letters
   character(len=1024),dimension(:),allocatable:: cuts2d_files
   integer(kind=4),dimension(:),allocatable:: thetablocks_data
   type(TermsInfo),dimension(:),allocatable:: phiTerms
   type(TermsInfo):: thetaTerms
   integer(kind=4):: n2dcuts
   integer(kind=4):: nthetablocks
   real(kind=8),dimension(:),allocatable:: phiList

   ! Lua specifications
   type(flu_State):: conf
   integer(kind=4):: ierr
   integer(kind=4):: pes_table,crp3d_table,vacfunc_table,dampfunc_table,param_table,extrapol_table
   integer(kind=4):: resize_table,wyckoff_table,inwyckoff_table,cut2d_table,files_table,kpoints_table
   integer(kind=4):: phi_table,term_table,theta_table,fourier_table,magnitude_table,phiList_table
   integer(kind=4):: inkpoints_table
   ! Auxiliar, dummy variables
   integer(kind=4):: auxInt,auxInt2
   real(kind=8):: auxReal
   character(len=1024):: auxString,auxString2
   type(Length):: len
   type(Angle):: angl
   ! Parameters
   character(len=*),parameter:: routinename="READ_CRP6D: "
   ! Run section -----------------------
   ! Open Lua config file
   call OPEN_CONFIG_FILE(L=conf,ErrCode=ierr,filename=filename)
   ! Open PES table
   call AOT_TABLE_OPEN(L=conf,thandle=pes_table,key=tablename)
   ! get pes.kind
   call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='kind',val=auxstring)
   call this%SET_PESTYPE(trim(auxstring))
   select case(trim(auxstring))
      case('CRP6D')
         ! do nothing
      case default
         write(0,*) 'READ CRP6D ERR: wrong kind of PES. Expected: CRP6D. Encountered: '//trim(auxstring)
   end select
#ifdef DEBUG
   call VERBOSE_WRITE(routinename,'PES type: '//trim(auxstring))
#endif
   ! get pes.name
   call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='name',val=auxstring)
   call this%SET_ALIAS(trim(auxstring))
#ifdef DEBUG
   call VERBOSE_WRITE(routinename,'PES name: '//trim(auxstring))
#endif
   ! get pes.dimensions
   call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='dimensions',val=auxint)
   call this%SET_DIMENSIONS(auxint)
#ifdef DEBUG
   call VERBOSE_WRITE(routinename,'PES dimensions: ',auxint)
#endif
   ! get crp3d subpes
   call AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=crp3d_table,key='crp3dPes')
   this%natomic=aot_table_length(L=conf,thandle=crp3d_table)
#ifdef DEBUG
   call VERBOSE_WRITE(routinename,"Atomic potentials found: ",this%natomic)
#endif
   select case(this%natomic)
      case(1)
         this%is_homonucl=.true.
#ifdef DEBUG
         call VERBOSE_WRITE(routinename,'This PES is for homonuclear projectiles')
#endif         
      case(2)
         this%is_homonucl=.false.
#ifdef DEBUG
         call VERBOSE_WRITE(routinename,'This PES is for heteronuclear projectiles')
#endif         
      case default 
         write(0,*) "READ_CRP6D ERR: Wrong number of atomic potentials. Allowed number: 1 or 2."
         call EXIT(1)
   end select
   allocate(this%atomiccrp(this%natomic))
   do i = 1, this%natomic
      call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=crp3d_table,pos=i,val=auxstring)
#ifdef DEBUG
      call VERBOSE_WRITE(routinename,'Atomic potential keyword: '//trim(auxstring))
#endif
      call this%atomiccrp(i)%INITIALIZE(filename=filename,tablename=trim(auxstring))
   end do
   call AOT_TABLE_CLOSE(L=conf,thandle=crp3d_table)
   ! get pes.vacuumFunction
   call AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=vacfunc_table,key='vacuumFunction')
   call AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=vacfunc_table,key='kind',val=auxstring)
   select case(trim(auxstring))
      case('Numerical')
         call AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=vacfunc_table,key='source',val=auxstring)
         call this%farpot%INITIALIZE(trim(auxstring))
      case default
         write(0,*) 'READ_CRP6D ERR: wrong kind of vacuuum function: '//trim(auxstring)
         write(0,*) 'Implemented ones: Numerical'
         write(0,*) 'Case sensitive'
         call EXIT(1)
   end select
   call AOT_TABLE_CLOSE(L=conf,thandle=vacfunc_table)
   ! get pes.dampFunction
   call AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=dampfunc_table,key='dampFunction')
   call AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=dampfunc_table,key='kind',val=auxstring)
   select case(trim(auxstring))
      case("Logistic")
         allocate(Logistic_func::this%dumpfunc)
         allocate(param_damp(2))
         call AOT_TABLE_OPEN(L=conf,parent=dampfunc_table,thandle=param_table,key='param')
         auxint=aot_table_length(L=conf,thandle=param_table)
         select case(auxint)
            case(2)
               ! do nothing
            case default
               write(0,*) 'READ_CRP6D ERR: incorrect number of parameters in damp function: ',auxint
               call EXIT(1)
         end select
         call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=param_table,pos=1,val=param_damp(1))
         call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=param_table,pos=2,val=param_damp(2))
         call AOT_TABLE_CLOSE(L=conf,thandle=param_table)
         call this%dumpfunc%READ(param_damp)
      case("fullCRP") 
         allocate(One_func::this%dumpfunc)
      case("fullRaw") 
         allocate(Zero_func::this%dumpfunc)
      case default
         write(0,*) "READ_CRP6D ERR: Keyword for dumping function needed"
         write(*,*) "Currently implemented: Logistic, fullCRP, fullRaw"
         call EXIT(1)
   end select
   call AOT_TABLE_CLOSE(L=conf,thandle=dampfunc_table)
   ! get pes.extrapolFunction
   call AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=extrapol_table,key='extrapolFunction')
   call AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=extrapol_table,key='kind',val=auxstring)
   select case(trim(auxstring))
      case("Xexponential","Linear")
         this%extrapol2vac_flag=trim(auxstring)
         call AOT_TABLE_OPEN(L=conf,parent=extrapol_table,thandle=param_table,key='upToZ')
         auxint=aot_table_length(L=conf,thandle=param_table)
         call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=param_table,pos=1,val=auxreal)
         call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=param_table,pos=2,val=auxstring)
         call len%READ(auxreal,trim(auxstring))
         call len%TO_STD()
         this%zvacuum=len%getvalue()
         call AOT_TABLE_CLOSE(L=conf,thandle=param_table)
      case("None")
         this%zvacuum=0.D0
      case default
         write(0,*) "READ_CRP6D ERR: Keyword for extrapolation function needed"
         write(0,*) "Currently implemented: None, Xexponential, Linear"
         write(0,*) "Case sensitive."
         call EXIT(1)
   end select
   call AOT_TABLE_CLOSE(L=conf,thandle=extrapol_table)
#ifdef DEBUG
   call VERBOSE_WRITE(routinename,"Z vacuum: ",this%zvacuum)
#endif
   ! get pes.resize
   call AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=resize_table,key='resize')
   auxint=aot_table_length(L=conf,thandle=resize_table)
   select case(auxint)
      case(0)
         ! do nothing, resize is not required
      case default
         call AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=resize_table,key='r',val=auxint)
         this%grid(1)=auxint
         call AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=resize_table,key='z',val=auxint)
         this%grid(2)=auxint
         if( this%grid(1)/=0 .and. this%grid(2)/=0 ) this%is_resized=.true.
   end select
#ifdef DEBUG
         call VERBOSE_WRITE(routinename,'New grid (R,Z):',this%grid(:))
         call verbose_write(routinename,'Is this PES going to be resized?: ',this%is_resized)
#endif
   call AOT_TABLE_CLOSE(L=conf,thandle=resize_table)
   ! get wyckoff sites
   call AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=wyckoff_table,key='wyckoffSite')
   this%nsites=aot_table_length(L=conf,thandle=wyckoff_table)
   allocate(wyckoff_letters(this%nsites))
   select case(this%nsites)
      case(0)
         write(0,*) "CRP3D_READ ERR: there aren't Wyckoff sites"
         call EXIT(1)
      case default
         ! do nothing
   end select
   ! Allocate with the correct type (symmetry) all wyckoff sites
   select case(system_surface%getsymmlabel())
      case("p4mm")
         allocate(Wyckoffp4mm::this%wyckoffsite(this%nsites))
#ifdef DEBUG
         call VERBOSE_WRITE(routinename,"Allocated p4mm Wyckoff sites")
#endif
      case default
         write(0,*) "READ_CRP6D ERR: surface symmetry is not implemented yet"
         write(0,*) "Good luck!"
         call EXIT(1)
   end select
   do i = 1, this%nsites
      call AOT_TABLE_OPEN(L=conf,parent=wyckoff_table,thandle=inwyckoff_table,pos=i)
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
      call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=inwyckoff_table,key='kind',val=wyckoff_letters(i))
      call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=inwyckoff_table,key='n2dcuts',val=n2dcuts)
      allocate(cuts2d_files(n2dcuts))
      nthetablocks=aot_table_length(L=conf,thandle=inwyckoff_table)-7
      allocate(thetablocks_data(nthetablocks))
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
      call VERBOSE_WRITE( routinename,'Wyckoff Site number: ',i)
      call VERBOSE_WRITE( routinename,'Wyckoff Letter: '//trim(wyckoff_letters(i)))
      call VERBOSE_WRITE( routinename,'Wyckoff number of cut2ds: ',n2dcuts)
      call VERBOSE_WRITE( routinename,'Wyckoff number of theta blocks: ',nthetablocks)
      call verbose_write( routinename,'Wyckoff theta terms irreps: ',thetaTerms%irrepList(:) )
      call verbose_write( routinename,'Wyckoff theta terms parities: ',thetaTerms%parityList(:) )
      call verbose_write( routinename,'Wyckoff theta terms kpoints: ',thetaTerms%kpointList(:) )
#endif
      allocate( phiTerms(nThetaBlocks) )
      do j = 1, nThetaBlocks
         call AOT_TABLE_OPEN(L=conf,parent=inwyckoff_table,thandle=cut2d_table,pos=j)
         call AOT_TABLE_OPEN(L=conf,parent=cut2d_table,thandle=files_table,key='files')
         thetablocks_data(j)=aot_table_length(L=conf,thandle=files_table)
#ifdef DEBUG
         call VERBOSE_WRITE(routinename,'Wyckoff thetablock: ',j)
         call VERBOSE_WRITE(routinename,'Wyckoff cut2d inside: ',thetablocks_data(j))
#endif
         select case(j)
            case(1)
               auxint=1
               auxint2=thetablocks_data(1)
            case default
               auxint=sum(thetablocks_data(1:j-1))+1
               auxint2=sum(thetablocks_data(1:j-1))+thetablocks_data(j)
         end select
         n=0
         do k = auxint, auxint2
            n=n+1
            call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=files_table,pos=n,val=cuts2d_files(k))
#ifdef DEBUG
            call VERBOSE_WRITE(routinename,'Wyckoff pos array: ',k)
            call VERBOSE_WRITE(routinename,'Wyckoff cut2d filename: '//trim(cuts2d_files(k)))
#endif
         end do
         call AOT_TABLE_CLOSE(L=conf,thandle=files_table)
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
         call AOT_TABLE_CLOSE(L=conf,thandle=cut2d_table)
      end do
      call this%wyckoffsite(i)%INITIALIZE(mynumber=i,letter=wyckoff_letters(i),nphipoints=thetablocks_data(:),&
                                          filenames=cuts2d_files(:),phiTerms=phiTerms,thetaTerms=thetaTerms,  &
                                          phiList=phiList(:) )
      call AOT_TABLE_CLOSE(L=conf,thandle=inwyckoff_table)
      deallocate(cuts2d_files)
      deallocate(thetablocks_data)
      deallocate(phiTerms)
      deallocate(phiList)
      call thetaTerms%reboot()
   end do
   call AOT_TABLE_CLOSE(L=conf,thandle=wyckoff_table)
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
   call CLOSE_CONFIG(conf)
   return
end subroutine READ_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_PHI_CRP6D #######################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values.
!
!> @warning
!! - The graph starts always at 0,0. Initial PHI value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 09/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
subroutine PLOT1D_PHI_CRP6D(thispes,npoints,X,filename)
   implicit none
   ! I/O variables -------------------------------
   class(CRP6D),intent(in) :: thispes
   integer,intent(in) :: npoints
   character(len=*),intent(in) :: filename
   real(kind=8),dimension(6),intent(in) :: X
   ! Local variables -----------------------------
   integer :: inpoints, ndelta
   real(kind=8) :: delta,v
   real(kind=8) :: xmax, xmin 
   real(kind=8),dimension(6) :: r, dvdu
   integer(kind=4) :: i ! Counter
   character(len=18),parameter :: routinename = "PLOT1D_PHI_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   if (npoints.lt.2) then
      write(0,*) "PLOT1D_PHI_CRP6D ERR: Less than 2 points"
      call EXIT(1)
   end if
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:5)=x(1:5)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   open(11,file=filename,status="replace")
   ! Initial value
   r(6)=xmin
   call thispes%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(6),v,dvdu(:)  
   ! cycle for inpoints
   do i=1, inpoints
      r(6)=xmin+DFLOAT(i)*delta
      call thispes%GET_V_AND_DERIVS(r,v,dvdu)
      write(11,*) r(6),v,dvdu(:)
   end do
   ! Final value
   r(6) = xmax
   call thispes%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(6),v,dvdu(:)
   write(*,*) routinename, "file created ",filename
   close(11)
end subroutine PLOT1D_PHI_CRP6D
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
subroutine PLOT1D_ATOMIC_INTERAC_PHI_CRP6D(thispes,npoints,X,filename)
   implicit none
   ! I/O variables -------------------------------
   class(CRP6D),intent(in) :: thispes
   integer,intent(in) :: npoints
   character(len=*),intent(in) :: filename
   real(kind=8),dimension(6),intent(in) :: X
   ! Local variables -----------------------------
   integer :: inpoints, ndelta
   real(kind=8) :: delta
   real(kind=8),dimension(2) :: v
   real(kind=8) :: xmax, xmin
   real(kind=8),dimension(6) :: r, dvdu, ratom
   integer(kind=4) :: i ! Counter
   character(len=18),parameter :: routinename = "PLOT1D_PHI_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   if (npoints.lt.2) then
      write(0,*) "PLOT1D_PHI_CRP6D ERR: Less than 2 points"
      call EXIT(1)
   end if
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:5)=x(1:5)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   open(11,file=filename,status="replace")
   ! Initial value
   r(6)=xmin
   call thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
   write(11,*) r(6),sum(v),v(1),v(2),ratom(:),dvdu(:)  
   ! cycle for inpoints
   do i=1, inpoints
      r(6)=xmin+DFLOAT(i)*delta
      call thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
      write(11,*) r(6),sum(v),v(1),v(2),ratom(:),dvdu(:)  
   end do
   ! Final value
   r(6) = xmax
   call thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
   write(11,*) r(6),sum(v),v(1),v(2),ratom(:),dvdu(:)  
   write(*,*) routinename, "file created ",filename
   close(11)
end subroutine PLOT1D_ATOMIC_INTERAC_PHI_CRP6D
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
subroutine PLOT1D_PHI_SMOOTH_CRP6D(thispes,npoints,X,filename)
   implicit none
   ! I/O variables -------------------------------
   class(CRP6D),intent(in) :: thispes
   integer,intent(in) :: npoints
   character(len=*),intent(in) :: filename
   real(kind=8),dimension(6),intent(in) :: X
   ! Local variables -----------------------------
   integer :: inpoints, ndelta
   real(kind=8) :: delta,v
   real(kind=8) :: xmax, xmin 
   real(kind=8),dimension(6) :: r, dvdu
   integer(kind=4) :: i ! Counter
   character(len=25),parameter :: routinename = "PLOT1D_PHI_SMOOTH_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   if (npoints.lt.2) then
      write(0,*) "PLOT1D_PHI_CRP6D ERR: Less than 2 points"
      call EXIT(1)
   end if
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:5)=x(1:5)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   open(11,file=filename,status="replace")
   ! Initial value
   r(6)=xmin
   call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(6),v,dvdu(:)  
   ! cycle for inpoints
   do i=1, inpoints
      r(6)=xmin+DFLOAT(i)*delta
      call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      write(11,*) r(6),v,dvdu(:)
   end do
   ! Final value
   r(6) = xmax
   call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(6),v,dvdu(:)
   write(*,*) routinename, "file created ",filename
   close(11)
end subroutine PLOT1D_PHI_SMOOTH_CRP6D
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
subroutine PLOT1D_THETA_CRP6D(thispes,npoints,X,filename)
   implicit none
   ! I/O variables -------------------------------
   class(CRP6D),intent(in) :: thispes
   integer, intent(in) :: npoints
   character(len=*),intent(in) :: filename
   real(kind=8),dimension(6),intent(in) :: X
   ! Local variables -----------------------------
   integer :: inpoints, ndelta
   real(kind=8) :: delta,v
   real(kind=8) :: xmax, xmin 
   real(kind=8), dimension(6) :: r, dvdu
   integer(kind=4) :: i ! Counter
   character(len=20),parameter :: routinename = "PLOT1D_THETA_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   if (npoints.lt.2) then
      write(0,*) "PLOT1D_THETA_CRP6D ERR: Less than 2 points"
      call EXIT(1)
   end if
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
   open(11,file=filename,status="replace")
   ! Initial value
   r(5)=xmin
   call thispes%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(5),v,dvdu(:)  
   ! cycle for inpoints
   do i=1, inpoints
      r(5)=xmin+DFLOAT(i)*delta
      call thispes%GET_V_AND_DERIVS(r,v,dvdu)
      write(11,*) r(5),v,dvdu(:)
   end do
   ! Final value
   r(5) = xmax
   call thispes%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(5),v,dvdu(:)
   write(*,*) routinename, "file created ",filename
   close(11)
end subroutine PLOT1D_THETA_CRP6D
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
subroutine PLOT1D_ATOMIC_INTERAC_THETA_CRP6D(thispes,npoints,X,filename)
   implicit none
   ! I/O variables -------------------------------
   class(CRP6D),intent(in) :: thispes
   integer, intent(in) :: npoints
   character(len=*),intent(in) :: filename
   real(kind=8),dimension(6),intent(in) :: X
   ! Local variables -----------------------------
   integer :: inpoints, ndelta
   real(kind=8) :: delta
   real(kind=8) :: xmax, xmin
   real(kind=8), dimension(6) :: r,dvdu,ratom
   real(kind=8),dimension(2) :: v
   integer(kind=4) :: i ! Counter
   character(len=20),parameter :: routinename = "PLOT1D_THETA_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   if (npoints.lt.2) then
      write(0,*) "PLOT1D_THETA_CRP6D ERR: Less than 2 points"
      call EXIT(1)
   end if
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
   open(11,file=filename,status="replace")
   ! Initial value
   r(5)=xmin
   call thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
   write(11,*) r(5),sum(v),v(:),ratom(:),dvdu(:)  

   ! cycle for inpoints
   do i=1, inpoints
      r(5)=xmin+DFLOAT(i)*delta
      call thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
      write(11,*) r(5),sum(v),v(:),ratom(:),dvdu(:)
   end do
   ! Final value
   r(5) = xmax
   call thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
   write(11,*) r(5),sum(v),v(:),ratom(:),dvdu(:)  
   write(*,*) routinename, "file created ",filename
   close(11)
end subroutine PLOT1D_ATOMIC_INTERAC_THETA_CRP6D
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
subroutine PLOT1D_THETA_SMOOTH_CRP6D(thispes,npoints,X,filename)
   implicit none
   ! I/O variables -------------------------------
   class(CRP6D),intent(in) :: thispes
   integer, intent(in) :: npoints
   character(len=*),intent(in) :: filename
   real(kind=8),dimension(6),intent(in) :: X
   ! Local variables -----------------------------
   integer :: inpoints, ndelta
   real(kind=8) :: delta,v
   real(kind=8) :: xmax, xmin
   real(kind=8), dimension(6) :: r, dvdu
   integer(kind=4) :: i ! Counter
   character(len=27),parameter :: routinename = "PLOT1D_THETA_SMOOTH_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   if (npoints.lt.2) then
      write(0,*) "PLOT1D_THETA_CRP6D ERR: Less than 2 points"
      call EXIT(1)
   end if
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
   open(11,file=filename,status="replace")
   ! Initial value
   r(5)=xmin
   call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(5),v,dvdu(:)  
   ! cycle for inpoints
   do i=1, inpoints
      r(5)=xmin+DFLOAT(i)*delta
      call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      write(11,*) r(5),v,dvdu(:)
   end do
   ! Final value
   r(5) = xmax
   call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(5),v,dvdu(:)
   write(*,*) routinename, "file created ",filename
   close(11)
end subroutine PLOT1D_THETA_SMOOTH_CRP6D
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
subroutine PLOT1D_R_CRP6D(thispes,npoints,X,filename)
   implicit none
   ! I/O variables -------------------------------
   class(CRP6D),intent(in) :: thispes
   integer, intent(in) :: npoints
   character(len=*),intent(in) :: filename
   real(kind=8),dimension(6),intent(in) :: X
   ! Local variables -----------------------------
   integer :: inpoints, ndelta
   real(kind=8) :: delta,v
   real(kind=8) :: xmax, xmin
   real(kind=8), dimension(6) :: r, dvdu
   integer(kind=4) :: i ! Counter
   character(len=16),parameter :: routinename = "PLOT1D_R_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   if (npoints.lt.2) then
      write(0,*) "PLOT1D_R_CRP6D ERR: Less than 2 points"
      call EXIT(1)
   end if
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
   open(11,file=filename,status="replace")
   ! Initial value
   r(4)=xmin
   call thispes%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(4),v,dvdu(:)  
   ! cycle for inpoints
   do i=1, inpoints
      r(4)=xmin+DFLOAT(i)*delta
      call thispes%GET_V_AND_DERIVS(r,v,dvdu)
      write(11,*) r(4),v,dvdu(:)
   end do
   ! Final value
   r(4) = xmax
   call thispes%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(4),v,dvdu(:)
   write(*,*) routinename, "file created ",filename
   close(11)
end subroutine PLOT1D_R_CRP6D
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
subroutine PLOT1D_R_SMOOTH_CRP6D(thispes,npoints,X,filename)
   implicit none
   ! I/O variables -------------------------------
   class(CRP6D),intent(in) :: thispes
   integer, intent(in) :: npoints
   character(len=*),intent(in) :: filename
   real(kind=8),dimension(6),intent(in) :: X
   ! Local variables -----------------------------
   integer :: inpoints, ndelta
   real(kind=8) :: delta,v
   real(kind=8) :: xmax, xmin
   real(kind=8), dimension(6) :: r, dvdu
   integer(kind=4) :: i ! Counter
   character(len=23),parameter :: routinename = "PLOT1D_R_SMOOTH_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   if (npoints.lt.2) then
      write(0,*) "PLOT1D_R_CRP6D ERR: Less than 2 points"
      call EXIT(1)
   end if
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
   open(11,file=filename,status="replace")
   ! Initial value
   r(4)=xmin
   call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(4),v,dvdu(:)  
   ! cycle for inpoints
   do i=1, inpoints
      r(4)=xmin+DFLOAT(i)*delta
      call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      write(11,*) r(4),v,dvdu(:)
   end do
   ! Final value
   r(4) = xmax
   call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(4),v,dvdu(:)
   write(*,*) routinename, "file created ",filename
   close(11)
end subroutine PLOT1D_R_SMOOTH_CRP6D
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
subroutine PLOT1D_Z_CRP6D(thispes,npoints,X,L,filename)
   implicit none
   ! I/O variables -------------------------------
   class(CRP6D),intent(in) :: thispes
   integer, intent(in) :: npoints
   character(len=*),intent(in) :: filename
   real(kind=8),dimension(6),intent(in) :: X
   real(kind=8),intent(in) :: L
   ! Local variables -----------------------------
   integer :: inpoints, ndelta
   real(kind=8) :: delta,v
   real(kind=8) :: xmax, xmin
   real(kind=8), dimension(6) :: r, dvdu
   integer(kind=4) :: i ! Counter
   character(len=16),parameter :: routinename = "PLOT1D_Z_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   if (npoints.lt.2) then
      write(0,*) "PLOT1D_Z_CRP6D ERR: Less than 2 points"
      call EXIT(1)
   end if
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
   open(11,file=filename,status="replace")
   ! Initial value
   r(3)=xmin
   call thispes%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(3),v,dvdu(:)  
   ! cycle for inpoints
   do i=1, inpoints
      r(3)=xmin+DFLOAT(i)*delta
      call thispes%GET_V_AND_DERIVS(r,v,dvdu)
      write(11,*) r(3),v,dvdu(:)
   end do
   ! Final value
   r(3) = xmax
   call thispes%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(3),v,dvdu(:)
   write(*,*) routinename, "file created ",filename
   close(11)
end subroutine PLOT1D_Z_CRP6D
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
subroutine PLOT1D_Z_SMOOTH_CRP6D(thispes,npoints,X,filename)
   implicit none
   ! I/O variables -------------------------------
   class(CRP6D),intent(in) :: thispes
   integer, intent(in) :: npoints
   character(len=*),intent(in) :: filename
   real(kind=8),dimension(6),intent(in) :: X
   ! Local variables -----------------------------
   integer :: inpoints, ndelta
   real(kind=8) :: delta,v
   real(kind=8) :: xmax, xmin
   real(kind=8), dimension(6) :: r, dvdu
   integer(kind=4) :: i ! Counter
   character(len=23),parameter :: routinename = "PLOT1D_Z_SMOOTH_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   select case(npoints)
      case(: 1)
         write(0,*) "PLOT1D_Z_CRP6D ERR: Less than 2 points"
         call EXIT(1)
      case default
         ! do nothing
   end select
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
   open(11,file=filename,status="replace")
   ! Initial value
   r(3)=xmin
   call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(3),v,dvdu(:)  
   ! cycle for inpoints
   do i=1, inpoints
      r(3)=xmin+DFLOAT(i)*delta
      call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      write(11,*) r(3),v,dvdu(:)
   end do
   ! Final value
   r(3) = xmax
   call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(3),v,dvdu(:)
   write(*,*) routinename, "file created ",filename
   close(11)
end subroutine PLOT1D_Z_SMOOTH_CRP6D
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
subroutine PLOT_XYMAP_CRP6D(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   implicit none
   class(CRP6D),intent(in) :: thispes
   real(kind=8),dimension(6),intent(in) :: init_point ! Initial position to start the scan (in a.u. and radians)
   integer,intent(in) :: nxpoints, nypoints ! number of points in XY plane
   character(len=*),intent(in) :: filename ! filename
   real(kind=8),intent(in) :: Lx ! Length of X axis 
   real(kind=8),intent(in) :: Ly ! Length of X axis 
   ! Local variables
   real(kind=8) :: xmin, ymin, xmax, ymax
   real(kind=8),dimension(6) :: r,dvdu
   real(kind=8) :: xdelta, ydelta
   integer :: xinpoints, nxdelta
   integer :: yinpoints, nydelta
   integer :: i, j ! counters
   real(kind=8) :: v ! potential
   ! Parameters section:
   character(len=*),parameter:: routineName='PLOT_XYMAP_CRP6D: '
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
   open(11,file=filename,status="replace")
   r(1) = xmin
   r(2) = ymin
   r(3:6)=init_point(3:6)
   call thispes%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(1:2),v,dvdu(:)
   do i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      call thispes%GET_V_AND_DERIVS(r,v,dvdu)
      write(11,*) r(1:2),v,dvdu(:)
   end do
   r(2) = ymax
   call thispes%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(1:2),v,dvdu(:)
   ! inpoints in XY
   do i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      r(3:6) = init_point(3:6)
      call thispes%GET_V_AND_DERIVS(r,v,dvdu)
      write(11,*) r(1:2),v,dvdu(:)
      do j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         call thispes%GET_V_AND_DERIVS(r,v,dvdu)
         write(11,*) r(1:2),v,dvdu(:)
      end do
      r(2) = ymax
      call thispes%GET_V_AND_DERIVS(r,v,dvdu)
      write(11,*) r(1:2),v,dvdu(:)
   end do
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   r(3:6) = init_point(3:6)
   call thispes%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(1:2),v,dvdu(:)
   do i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      call thispes%GET_V_AND_DERIVS(r,v,dvdu)
      write(11,*) r(1:2),v,dvdu(:)
   end do
   r(2) = ymax
   call thispes%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(1:2),v,dvdu(:)
   close(11)
   write(*,*) routineName//'Graph created: '//fileName
   return
end subroutine PLOT_XYMAP_CRP6D
!######################################################################################
! SUBROUTINE: PLOT_XYMAP_SMOOTH_CRP6D
!######################################################################################
!> @brief
!! Same as plot_xymap_crp6d but calling get_v_and_derivs_smooth, i.e. potential
!! and derivatives are not corrected. Thus, we get smooth corrugationless potential.
!-------------------------------------------------------------------------------------
subroutine PLOT_XYMAP_SMOOTH_CRP6D(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   implicit none
   class(CRP6D),intent(in) :: thispes
   real(kind=8),dimension(6),intent(in) :: init_point ! Initial position to start the scan (in a.u. and radians)
   integer,intent(in) :: nxpoints, nypoints ! number of points in XY plane
   character(len=*),intent(in) :: filename ! filename
   real(kind=8),intent(in) :: Lx ! Length of X axis
   real(kind=8),intent(in) :: Ly ! Length of X axis
   ! Local variables
   real(kind=8) :: xmin, ymin, xmax, ymax
   real(kind=8),dimension(6) :: r,dvdu
   real(kind=8) :: xdelta, ydelta
   integer :: xinpoints, nxdelta
   integer :: yinpoints, nydelta
   integer :: i, j ! counters
   real(kind=8) :: v ! potential
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
   open(11,file=filename,status="replace")
   r(1) = xmin
   r(2) = ymin
   r(3:6)=init_point(3:6)
   call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(1:2),v,dvdu(:)
   do i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      write(11,*) r(1:2),v,dvdu(:)
   end do
   r(2) = ymax
   call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(1:2),v,dvdu(:)
   ! inpoints in XY
   do i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      r(3:6) = init_point(3:6)
      call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      write(11,*) r(1:2),v,dvdu(:)
      do j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
         write(11,*) r(1:2),v,dvdu(:)
      end do
      r(2) = ymax
      call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      write(11,*) r(1:2),v,dvdu(:)
   end do
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   r(3:6) = init_point(3:6)
   call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(1:2),v,dvdu(:)
   do i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      write(11,*) r(1:2),v,dvdu(:)
   end do
   r(2) = ymax
   call thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(1:2),v,dvdu(:)
   close(11)
   return
end subroutine PLOT_XYMAP_SMOOTH_CRP6D
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
subroutine PLOT_RZMAP_CRP6D(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   implicit none
   class(CRP6D),intent(in) :: thispes
   real(kind=8),dimension(6),intent(in) :: init_point 
   integer,intent(in) :: nxpoints, nypoints 
   character(len=*),intent(in) :: filename 
   real(kind=8),intent(in) :: Lx
   real(kind=8),intent(in) :: Ly
   ! Local variables
   real(kind=8) :: xmin, ymin, xmax, ymax
   real(kind=8),dimension(6) :: r,dvdu
   real(kind=8) :: xdelta, ydelta
   integer :: xinpoints, nxdelta
   integer :: yinpoints, nydelta
   integer :: i, j ! counters
   real(kind=8) :: v ! potential
   ! Parameter section
   character(len=*),parameter:: routineName='PLOT_RZMAP_CRP6D: '
   integer(kind=4),parameter:: wunit=11 ! write unit
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
   open(unit=wunit,file=filename,status="replace")
   write(wunit,*) '# Initial geometry: ',init_point(:)
   write(wunit,*) '# X length (bohr): ',Lx
   write(wunit,*) '# Y length (bohr): ',Ly
   write(wunit,*) '# Grid: ',nxpoints,nypoints
   r(4) = xmin
   r(3) = ymin
   r(1:2)=init_point(1:2)
   r(5:6)=init_point(5:6)
   call thispes%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(4),r(3),v,dvdu(:)
   do i =1, yinpoints
      r(3) = ymin + DFLOAT(i)*ydelta
      call thispes%GET_V_AND_DERIVS(r,v,dvdu)
      write(wunit,*) r(4),r(3),v,dvdu(:)
   end do
   r(3) = ymax
   call thispes%GET_V_AND_DERIVS(r,v,dvdu)
   write(wunit,*) r(4),r(3),v,dvdu(:)
   ! inpoints in XY
   do i = 1, xinpoints
      r(4) = xmin+DFLOAT(i)*xdelta
      r(3) = ymin
      call thispes%GET_V_AND_DERIVS(r,v,dvdu)
      write(wunit,*) r(4),r(3),v,dvdu(:)
      do j = 1, yinpoints
         r(3) = ymin + DFLOAT(j)*ydelta
         call thispes%GET_V_AND_DERIVS(r,v,dvdu)
         write(wunit,*) r(4),r(3),v,dvdu(:)
      end do
      r(3) = ymax
      call thispes%GET_V_AND_DERIVS(r,v,dvdu)
      write(wunit,*) r(4),r(3),v,dvdu(:)
   end do
   ! Last point in XY plane
   r(4) = xmax
   r(3) = ymax
   call thispes%GET_V_AND_DERIVS(r,v,dvdu)
   write(wunit,*) r(4),r(3),v,dvdu(:)
   do i =1, yinpoints
      r(3) = ymin + DFLOAT(i)*ydelta
      call thispes%GET_V_AND_DERIVS(r,v,dvdu)
      write(wunit,*) r(4),r(3),v,dvdu(:)
   end do
   r(3) = ymax
   call thispes%GET_V_AND_DERIVS(r,v,dvdu)
   write(wunit,*) r(4),r(3),v,dvdu(:)
   close(wunit)
   write(*,*) routineName//'Graph created: '//fileName
   return
end subroutine PLOT_RZMAP_CRP6D
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
subroutine PLOT_ATOMIC_INTERAC_RZ_CRP6D(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   implicit none
   class(CRP6D),intent(in) :: thispes
   real(kind=8),dimension(6),intent(in) :: init_point 
   integer,intent(in) :: nxpoints, nypoints 
   character(len=*),intent(in) :: filename 
   real(kind=8),intent(in) :: Lx
   real(kind=8),intent(in) :: Ly
   ! Local variables
   real(kind=8) :: xmin, ymin, xmax, ymax
   real(kind=8),dimension(6) :: r
   real(kind=8) :: xdelta, ydelta
   integer :: xinpoints, nxdelta
   integer :: yinpoints, nydelta
   real(kind=8),dimension(2) :: v
   real(kind=8),dimension(6) :: atomicx
   real(kind=8),dimension(6) :: dvdu
   integer :: i, j ! counters
   ! parameters
   character(len=*),parameter:: routineName='PLOT_ATOMIC_INTERAC_RZ_CRP6D: '
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
   open(11,file=filename,status="replace")
   r(4) = xmin
   r(3) = ymin
   r(1:2)=init_point(1:2)
   r(5:6)=init_point(5:6)
   call thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
   write(11,*) r(4),r(3),sum(v),v(1),v(2)
   do i =1, yinpoints
      r(3) = ymin + DFLOAT(i)*ydelta
      call thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
      write(11,*) r(4),r(3),sum(v),v(1),v(2)
   end do
   r(3) = ymax
   call thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
   write(11,*) r(4),r(3),sum(v),v(1),v(2)
   ! inpoints in XY
   do i = 1, xinpoints
      r(4) = xmin+DFLOAT(i)*xdelta
      r(3) = ymin
      call thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
      write(11,*) r(4),r(3),sum(v),v(1),v(2)
      do j = 1, yinpoints
         r(3) = ymin + DFLOAT(j)*ydelta
         call thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
         write(11,*) r(4),r(3),sum(v),v(1),v(2)
      end do
      r(3) = ymax
      call thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
      write(11,*) r(4),r(3),sum(v),v(1),v(2)
   end do
   ! Last point in XY plane
   r(4) = xmax
   r(3) = ymax
   call thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
   write(11,*) r(4),r(3),sum(v),v(1),v(2)
   do i =1, yinpoints
      r(3) = ymin + DFLOAT(i)*ydelta
      call thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
      write(11,*) r(4),r(3),sum(v),v(1),v(2)
   end do
   r(3) = ymax
   call thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
   write(11,*) r(4),r(3),sum(v),v(1),v(2)
   close(11)
   write(*,*) routineName//'Graph created: '//fileName
   return
end subroutine PLOT_ATOMIC_INTERAC_RZ_CRP6D
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
logical function is_allowed_CRP6D(this,x) 
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP6D),intent(in) :: this
   real(kind=8),dimension(:),intent(in) :: x
   ! Local variables
   real(kind=8) :: zmin, rmin, rmax
   ! Run section
   zmin=this%wyckoffsite(1)%zrcut(1)%getfirstZ()
   rmin=this%wyckoffsite(1)%zrcut(1)%getfirstR()
   rmax=this%wyckoffsite(1)%zrcut(1)%getlastR()
   select case(size(x)/=6)
      case(.true.)
         write(0,*) "is_allowed_CRP6D ERR: wrong number of dimensions"
         call EXIT(1)
      case(.false.)
         ! do nothing
   end select
   select case(x(3)<zmin .or. x(4)<rmin .or. x(4)>rmax)
      case(.true.)
         is_allowed_CRP6D=.false.
      case(.false.)
         is_allowed_CRP6D=.true.
   end select
   return
end function is_allowed_CRP6D

end module CRP6D_MOD
