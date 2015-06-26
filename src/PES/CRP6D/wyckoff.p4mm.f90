!#########################################################
! MODULE: WYCKOFF_P4MM
!> @brief
!! Provides tools to interpolate through wyckoff sites belonging to
!! p4mm wallpaper symmetry
!##########################################################
module WYCKOFF_P4MM_MOD
! Initial declarations
use WYCKOFF_GENERIC_MOD
use FOURIER1D_2_MOD
use FOURIER1D_4MM_MOD
use FOURIER1D_M_MOD
use FOURIER1D_M45_MOD
use FOURIER1D_MM2_MOD
use FOURIER1D_M45M135_2_MOD
use FOURIER1D_E_MOD
#ifdef DEBUG
use DEBUG_MOD
use UNITS_MOD
#endif
implicit none
!/////////////////////////////////////////////////////////////////
! TYPE: Wyckoffp4mm
!> @brief
!! Subclass of Wyckoffsitio for generic p4mm symmetry.
!! in p4mm symmetry does not matter
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 20/03/2014 
!> @version 1.0
!----------------------------------------------------------------
type,extends(Wyckoffsitio) :: Wyckoffp4mm
   contains
      procedure,public :: GET_V_AND_DERIVS => GET_V_AND_DERIVS_WYCKOFFP4MM
end type Wyckoffp4mm
!/////////////////////////////////////////////////////////////////
contains
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_WYCKOFFP4MM 
!###########################################################
!> @brief
!! Symmetry adapted interpolation for p4mm wallpaper group.
!
!> @param[in] x - Array which stands for Z,r,theta,phi
!> @param[out] v - potential at X
!> @param[out] dvdu - Array with derivatives: dvdz, dvdr, dvdtheta dvdphi
!
!> @warning
!! - Only stands for a,b,c, and f p4mm Wyckoff sites.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.2
!-----------------------------------------------------------
subroutine GET_V_AND_DERIVS_WYCKOFFP4MM(this,x,v,dvdu)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(Wyckoffp4mm),intent(in) :: this
   real(kind=8),dimension(4),intent(in) ::x 
   real(kind=8),intent(out) :: v
   real(kind=8),dimension(4),intent(out) :: dvdu
   ! Local variables
   integer(kind=4) :: i,j,h ! counters
   real(kind=8) :: z,r,theta,phi
   class(Fourier1d),dimension(:),allocatable :: phicut
   class(Fourier1d),allocatable :: thetacut
   real(kind=8) :: aux_theta
   real(kind=8),dimension(:,:),allocatable :: aux
   real(kind=8),dimension(:),allocatable :: f
   real(kind=8),dimension(:),allocatable :: dfdz
   real(kind=8),dimension(:),allocatable :: dfdr
   real(kind=8),dimension(:),allocatable :: dfdphi
   real(kind=8),dimension(:),allocatable :: philist
   real(kind=8),dimension(:),allocatable :: thetalist
   character(len=2) :: theta_irrep, phi_irrep
   character(len=30),parameter :: routinename="GET_V_AND_DERIVS_WYCKOFFP4MM: "
#ifdef DEBUG
   type(Angle),dimension(:),allocatable :: beta
   real(kind=8),dimension(:),allocatable :: philistdeg
   character(len=18) :: filename
#endif
   ! Run section
   z=x(1)
   r=x(2)
   theta=x(3)
   phi=x(4)
   ! CONDITIONS: edit this part to include/edit symmetries
   select case( this%id )
      case("a" : "b")
         allocate( Fourier1d_4mm::phiCut(this%nphicuts) )
         allocate( Fourier1d_mm2::thetaCut )
         
      case("c")
         allocate( Fourier1d_mm2::phiCut(this%nphicuts) )
         allocate( Fourier1d_mm2::thetaCut )

      case("f")
         allocate( Fourier1d_m45::phiCut(this%nphicuts) )
         allocate( Fourier1d_2::thetaCut )

      case default
         write(0,*) "GET_V_AND_DERIVS_WYCKOFFP4MM: Unexpected error with Wyckoff id"
         call EXIT(1)
   end select
   ! PHI INTERPOLATION ----------------------------------------------------
   h=0 ! initialize h
   do i = 1, this%nphicuts ! loop over specific phi cuts (each one for a different theta value) 
      allocate(f(this%nphipoints(i)))
      allocate(dfdr(this%nphipoints(i)))
      allocate(dfdz(this%nphipoints(i)))
      allocate(philist(this%nphipoints(i)))
      allocate(aux(2,this%nphipoints(i)))
      do j = 1, this%nphipoints(i) ! loop over number of zrcuts inside
         h=h+1 ! numbering of zrcuts
         f(j)=this%zrcut(h)%interrz%getvalue((/r,z/)) ! storing potential at this site
         dfdr(j)=this%zrcut(h)%interrz%getderivx((/r,z/)) ! storing d/dr at this site
         dfdz(j)=this%zrcut(h)%interrz%getderivy((/r,z/)) ! storing d/dz at this site
         philist(j)=this%zrcut(h)%phi
         aux_theta=this%zrcut(h)%theta
      end do
      aux(1,:)=dfdz(:)
      aux(2,:)=dfdr(:)
      call phiCut(i)%read(philist,f)
      call phiCut(i)%add_morefuncs(aux)
      call phiCut(i)%setKlist( this%phiTerms(i)%kpointList(:) )
      call phiCut(i)%setParityList( this%phiTerms(i)%parityList(:) )
      call phiCut(i)%setIrrepList( this%phiTerms(i)%irrepList(:) )
      call phiCut(i)%initializeTerms()
      call phiCut(i)%interpol()
#ifdef DEBUG
      call DEBUG_WRITE(routinename,"NEW PHICUT")
      call DEBUG_WRITE(routinename,"For theta: ",this%zrcut(h)%theta)
      call DEBUG_WRITE(routinename,"Is homonuclear: ",this%is_homonucl)
      allocate(beta(size(philist)))
      allocate(philistdeg(size(philist)))
      do j = 1, size(philist)
         call beta(j)%READ(philist(j),"rad")
         call beta(j)%TO_DEG()
         philistdeg(j)=beta(j)%getvalue()
      end do
      call DEBUG_WRITE(routinename,"At Phi: (deg)",philistdeg)
      call DEBUG_WRITE(routinename,"At Phi: (rad)",philist)
      call DEBUG_WRITE(routinename,"Klist: ",phicut(i)%getklist())
      call DEBUG_WRITE(routinename,"f: ",f)
      call DEBUG_WRITE(routinename,"dfdz: ",aux(1,:))
      call DEBUG_WRITE(routinename,"dfdr: ",aux(2,:))
      select case(get_debugmode())
         case(.true.)
            write(filename,'(I1,A1,I1,A1,A14)') this%mynumber,"-",i,"-","wyckoffphi.raw" 
            call phicut(i)%PLOTDATA(filename)
            write(filename,'(I1,A1,I1,A1,A14)') this%mynumber,"-",i,"-","wyckoffphi.cyc" 
            call phicut(i)%PLOTCYCLIC_ALL(300,filename)
         case(.false.)
            ! do nothing
      end select
      deallocate(beta)
      deallocate(philistdeg)
#endif
      deallocate(philist)
      deallocate(f)
      deallocate(dfdr)
      deallocate(dfdz)
      deallocate(aux)
   end do
   ! THETA INTERPOLATION --------------------------
   allocate(thetalist(this%nphicuts))
   allocate(f(this%nphicuts))
   allocate(dfdz(this%nphicuts))
   allocate(dfdr(this%nphicuts))
   allocate(dfdphi(this%nphicuts))
   h=0 ! reboot h
   do i = 1, this%nphicuts
      do j = 1, this%nphipoints(i)
         h=h+1 
         thetalist(i)=this%zrcut(h)%theta 
      end do
      allocate(aux(2,3))
      call phicut(i)%GET_ALLFUNCS_AND_DERIVS(phi,aux(1,:),aux(2,:))
      f(i)=aux(1,1)
      dfdz(i)=aux(1,2)
      dfdr(i)=aux(1,3)
      dfdphi(i)=aux(2,1)
      deallocate(aux)
   end do
   allocate(aux(3,this%nphicuts))
   aux(1,:)=dfdz(:)
   aux(2,:)=dfdr(:)
   aux(3,:)=dfdphi(:)
   call thetaCut%read(thetalist,f)
   call thetaCut%add_morefuncs(aux)
   call thetaCut%setKlist     ( this%thetaTerms%kpointList(:)    )
   call thetaCut%setParityList( this%thetaTerms%parityList(:) )
   call thetaCut%setIrrepList ( this%thetaTerms%irrepList(:)  )
   call thetaCut%initializeTerms()
   call thetaCut%interpol()
#ifdef DEBUG
   call DEBUG_WRITE(routinename,"NEW THETACUT")
   allocate(beta(size(thetalist)))
   allocate(philistdeg(size(thetalist)))
   do j = 1, size(thetalist)
      call beta(j)%READ(thetalist(j),"rad")
      call beta(j)%TO_DEG()
      philistdeg(j)=beta(j)%getvalue()
   end do
   call DEBUG_WRITE(routinename,"At Theta: (deg) ",philistdeg)
   call DEBUG_WRITE(routinename,"At Theta: (rad) ",thetalist)
   call DEBUG_WRITE(routinename,"Klist:          ",thetacut%getklist())
   call DEBUG_WRITE(routinename,"f:              ",f)
   call DEBUG_WRITE(routinename,"dfdz:           ",dfdz)
   call DEBUG_WRITE(routinename,"dfdr:           ",dfdr)
   call DEBUG_WRITE(routinename,"dfdphi:         ",dfdphi)
   select case(get_debugmode())
      case(.true.)
         write(filename,'(I1,A1,A16)') this%mynumber,"-","wyckofftheta.raw" 
         call thetacut%PLOTDATA(filename)
         write(filename,'(I1,A1,A16)') this%mynumber,"-","wyckofftheta.cyc" 
         call thetacut%PLOTCYCLIC(300,filename)
      case(.false.)
         ! do nothing
   end select
   deallocate(beta)
   deallocate(philistdeg)
#endif
   deallocate(f)
   deallocate(dfdr)
   deallocate(dfdz)
   deallocate(dfdphi)
   deallocate(aux)
   deallocate(thetalist)
   allocate(aux(2,4))
   call thetacut%GET_ALLFUNCS_AND_DERIVS(theta,aux(1,:),aux(2,:))
   v=aux(1,1) ! value of the potential
   dvdu(1)=aux(1,2) ! dvdz
   dvdu(2)=aux(1,3) ! dvdr
   dvdu(3)=aux(2,1) ! dvdtheta
   dvdu(4)=aux(1,4) ! dvdphi
   deallocate(aux)
   return
end subroutine GET_V_AND_DERIVS_WYCKOFFP4MM
end module WYCKOFF_P4MM_MOD
