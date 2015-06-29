!##################################################################################
! MODULE: FOURIER2D_MOD
!> @brief
!! Provides tools to perform 2D interpolations with fourier series
!##################################################################################
module FOURIER2D_MOD
use SURFACE_MOD, only: Surface
implicit none
!////////////////////////////////////////////////////////////////
! TYPE: Interpol2d
!
!> @brief
!! Generic 2D interpolation type variable
!
!> @param n - Number of data points
!> @param xy(:,:) - Matrix that collects @f$(x_{i},y_{i})@f$ pairs
!> @param x(:) - Grid in X. Only if input has grid structure
!> @param y(:) - Grid in Y. Only if input has grid structure
!> @param fgrid(:,:) - Function evaluated in a grid
!> @param f - Array that stores couples @f$F(x_{i},y_{i})@f$. Non grid input.
!> @param dfdz - Array that stores couples @f$\frac{\partial F(x_{i},y_{i})}{\partial z}@f$
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Feb/2014, Jul/2014
!> @version 2.0 
!---------------------------------------------------------------
type,abstract :: Fourier2d
   ! public atributes
   integer(kind=4),public :: n
   integer(kind=4),public :: nfunc
   real(kind=8),dimension(:,:),allocatable,public :: xy
   real(kind=8),dimension(:,:),allocatable,public :: f
   integer(kind=4),dimension(:,:),allocatable,public:: kList
   character(len=1),dimension(:),allocatable,public:: parityList
   character(len=2),dimension(:),allocatable,public:: irrepList
contains
   ! initialize block
   procedure,public,non_overridable :: read => read_FOURIER2D
   ! tools block
   procedure(interpol_FOURIER2D),public,deferred :: interpol  ! deferred !!!! take a look to interface
   procedure(get_f_and_derivs_FOURIER2D),public,deferred :: get_f_and_derivs ! deferred !!!! take a look to interface
end type Fourier2d

abstract interface
   !###########################################################
   !# SUBROUTINE: INTERPOL_FOURIER2D 
   !###########################################################
   !-----------------------------------------------------------
   subroutine interpol_FOURIER2D(this,surf,filename)
      import fourier2d
      import surface
      class(fourier2d),intent(inout) :: this
      class(surface),intent(in) :: surf
      character(len=*),intent(in),optional :: filename
   end subroutine interpol_FOURIER2D
   !###########################################################
   !# SUBROUTINE: GET_F_AND_DERIVS 
   !###########################################################
   !-----------------------------------------------------------
   subroutine get_f_and_derivs_FOURIER2D(this,surf,r,v,dvdu)
      import fourier2d
      import surface
      class(fourier2d),intent(in):: this
      class(surface),intent(in):: surf
      real(kind=8),dimension(2),intent(in) :: r
      real(kind=8),dimension(:),intent(out) :: v
      real(kind=8),dimension(:,:),intent(out) :: dvdu
   end subroutine get_f_and_derivs_FOURIER2D
end interface
!////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: READ_FOURIER2D 
!###########################################################
!> @brief
!! Reads xy, f and klist values from arguments
!
!> @param[out] this - Interpol2d class to be read
!> @param[in] xy - Matrix that collects @f$(x_{i},y_{i})@f$ pairs
!> @param[in] f - Values of @f$F(x_{i},y_{i})@f$. It is a matrix so that
!!                the user can provide in each row a different function that
!!                will have the same @f$T^{-1}@f$ matrix during the
!!                interpolation
!> @param[in] klist - Kpoints to be used in the expansion. There should be as many
!!                    of them as numbers of evaluations of f
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Mar/2014
!> @version 1.0 
!-----------------------------------------------------------
subroutine read_FOURIER2D(this,xy,f,kList,irrepList,parityList)
   ! initial declarations
   implicit none
   ! i/o variables
   class(fourier2d),intent(inout):: this
   real(kind=8),dimension(:,:),intent(in):: xy
   real(kind=8),dimension(:,:),intent(in):: f
   integer(kind=4),dimension(:,:),intent(in):: kList
   character(len=1),dimension(:):: parityList
   character(len=2),dimension(:):: irrepList
   ! local variables
   integer(kind=4) :: ndata,nfunc
   ! run section
   ndata=size(f(1,:))
   nfunc=size(f(:,1))
   select case(size(xy(:,1)) == ndata .and. size(xy(1,:))==2)
      case(.false.)
         write(0,*) "READ_FOURIER2D ERR: dimensions mismatch in arrays xy or f"
         call exit(1)
      case(.true.)
         ! do nothing
   end select
   this%n=ndata
   allocate( this%xy(ndata,2),       source=xy(:,:)       )
   allocate( this%f(nfunc,ndata),    source=f(:,:)        )
   allocate( this%kList(ndata,2),    source=kList(:,:)    )
   allocate( this%parityList(ndata), source=parityList(:) )
   allocate( this%irrepList(ndata),  source=irrepList(:)  )
   this%nfunc = nfunc
   return
end subroutine read_fourier2d

end module FOURIER2D_MOD
