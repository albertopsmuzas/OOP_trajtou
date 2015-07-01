!########################################################
! MODULE : FOURIER3D_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes Interpol1d_mod in its scope
!########################################################
module FOURIER3D_MOD
use INTERPOL3D_MOD, only: Interpol3d
use MATHS_MOD, only: inv_mtrx
use FOURIER1D_MOD, only: Fourier1d
use FOURIER2D_MOD, only: Fourier2d
#ifdef DEBUG
use DEBUG_MOD, only: verbose_write, debug_write
#endif
implicit none
!/////////////////////////////////////////////////////////////////
! TYPE: Termcalculator
!> @brief
!! Abstract class to calculate terms of the series avoiding the use of
!! unnecessary switches
!----------------------------------------------------------------
type,abstract:: TermCalculator3d
   ! public atributes
   class(Fourier1d),allocatable,public:: angleFourier
   class(Fourier2d),allocatable,public:: xyFourier
   contains
      procedure(getvalue_termcalculator_example),public,deferred:: getValue
      procedure(getvalue_termcalculator_example),public,deferred:: getDeriv1
      procedure(getvalue_termcalculator_example),public,deferred:: getDeriv2
      procedure(getvalue_termcalculator_example),public,deferred:: getDeriv3
end type TermCalculator3d
!
abstract interface
   !###########################################################
   !# FUNCTION: getvalue_termcalculator_example 
   !###########################################################
   !> @brief 
   !! Just an example that child objects should override
   !-----------------------------------------------------------
   function getvalue_termcalculator_example(this,k,parityXY,irrepXY,l,parityAngle,irrepAngle,x) result(answer)
      import TermCalculator3d
      class(TermCalculator3d),intent(in):: this
      integer(kind=4),dimension(2),intent(in):: k
      integer(kind=4),intent(in):: l
      character(len=1),intent(in):: parityXY
      character(len=1),intent(in):: parityAngle
      character(len=2),intent(in):: irrepXY
      character(len=2),intent(in):: irrepAngle
      real(kind=8),dimension(3),intent(in):: x
      real(kind=8):: answer
   end function getvalue_termcalculator_example
   !-------------------------------------------------------------
end interface
!
!/////////////////////////////////////////////////////////////////////////////
! TYPE: FOURIER3D
!> @brief
!! Class to store all information needed for a 3D REAL combined fourier interpolation
!----------------------------------------------------------------------------
type,abstract,extends(Interpol3d):: Fourier3d
   ! public atributes
   class(TermCalculator3d),allocatable,public:: term
   integer(kind=4),dimension(:,:),allocatable,public:: kListXY
   character(len=1),dimension(:),allocatable,public:: parityListXY
   character(len=2),dimension(:),allocatable,public:: irrepListXY
   integer(kind=4),dimension(:),allocatable,public:: kListAngle
   character(len=1),dimension(:),allocatable,public:: parityListAngle
   character(len=2),dimension(:),allocatable,public:: irrepListAngle
   ! private atributes
   real(kind=8),dimension(:),allocatable,private:: coeff
   real(kind=8),dimension(:,:),allocatable,private:: extracoeff
   real(kind=8),dimension(:,:),allocatable,private:: extrafuncs
   contains
      ! initialize block
      procedure(initializeTerms_FOURIER3D),public,deferred:: initializeTerms
      ! get block
      procedure,public,non_overridable:: getValue  => getvalue_FOURIER3D
      procedure,public,non_overridable:: getDeriv1 => getDeriv1_FOURIER3D
      procedure,public,non_overridable:: getDeriv2 => getDeriv2_FOURIER3D
      procedure,public,non_overridable:: getDeriv3 => getDeriv3_FOURIER3D
      ! set block
      procedure,public,non_overridable:: setKlist => setKlist_FOURIER3D
      procedure,public,non_overridable:: setParityList => setParityList_FOURIER3D
      procedure,public,non_overridable:: setIrrepList => setIrrepList_FOURIER3D
      ! tools
      procedure,public,non_overridable:: interpol => interpol_FOURIER3D
      procedure,public,non_overridable:: add_morefuncs => add_more_funcs_FOURIER3D
      procedure,public,non_overridable:: get_allfuncs_and_derivs => get_allfunc_and_derivs_FOURIER3D
end type Fourier3d
abstract interface
   !###########################################################
   !# SUBROUTINE: SET_IRREP_FOURIER3D
   !###########################################################
   !> @brief
   !! Sets irrep for this fourier series. Should be overriden by
   !! child non-abstract classes
   !-----------------------------------------------------------
   subroutine initializeTerms_FOURIER3D(this)
      import Fourier3d
      class(Fourier3d),intent(inout):: this
   end subroutine initializeTerms_FOURIER3D
end interface
!/////////////////////////////////////////////////////////////////////////////
contains
!###################################################################
!# SUBROUTINE: setKlist_FOURIER3D
!###################################################################
!> @brief
!! Common set subroutine. Sets Klist atribute of a FOURIER3D object
!-------------------------------------------------------------------
subroutine setKlist_FOURIER3D(this,kListXY,kListAngle)
   implicit none
   ! I/O variables
   class(Fourier3d),intent(inout):: this
   integer(kind=4),dimension(:,:),intent(in):: kListXY
   integer(kind=4),dimension(:),intent(in):: kListAngle
   ! Local variables
   integer(kind=4):: n
   ! Run section
   n=size(kListXY(:,1))
   allocate( this%kListXY(n,2),  source=kListXY(:,1:2)  )
   allocate( this%kListAngle(n), source=kListAngle(:) )
   return
end subroutine setKlist_FOURIER3D
!###################################################################
!# SUBROUTINE: setParityList_FOURIER3D
!###################################################################
!> @brief
!! Common set subroutine. Sets parityList atribute of a FOURIER3D object
!-------------------------------------------------------------------
subroutine setParityList_FOURIER3D(this,parityListXY,parityListAngle)
   implicit none
   ! I/O variables
   class(Fourier3d),intent(inout):: this
   character(len=1),dimension(:),intent(in):: parityListXY
   character(len=1),dimension(:),intent(in):: parityListAngle
   ! Local variables
   integer(kind=4):: n
   ! Run section
   n=size(parityListXY)
   allocate( this%parityListXY(n),    source=parityListXY(:)    )
   allocate( this%parityListAngle(n), source=parityListAngle(:) )
   return
end subroutine setParityList_FOURIER3D
!###################################################################
!# SUBROUTINE: setIrrepList_FOURIER3D
!###################################################################
!> @brief
!! Common set subroutine. Sets IrrepList atribute of a FOURIER3D object
!-------------------------------------------------------------------
subroutine setIrrepList_FOURIER3D(this,irrepListXY,irrepListAngle)
   implicit none
   ! I/O variables
   class(Fourier3d),intent(inout):: this
   character(len=2),dimension(:),intent(in):: irrepListXY
   character(len=2),dimension(:),intent(in):: irrepListAngle
   ! Local variables
   integer(kind=4):: n
   ! Run section
   n=size( irrepListXY )
   allocate( this%irrepListXY(n),    source=irrepListXY(:)    )
   allocate( this%irrepListAngle(n), source=irrepListAngle(:) )
   return
end subroutine setIrrepList_FOURIER3D
!###########################################################
!# SUBROUTINE: GET_ALLFUNC_AND_DERIVS_FOURIER3D
!###########################################################
!> @brief
!! Get value of the potential and derivs for an specific point x
!! for the main function and extra ones
!-----------------------------------------------------------
subroutine GET_ALLFUNC_AND_DERIVS_FOURIER3D(this,x,f,dfdx)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(Fourier3d),intent(in) :: this
   real(kind=8),dimension(3),intent(in) :: x
   real(kind=8),dimension(:),intent(out) :: f
   real(kind=8),dimension(:,:),intent(out):: dfdx
   ! Local variables
   integer(kind=4) :: nfuncs
   integer(kind=4) :: i ! counters
   real(kind=8),dimension(:),allocatable :: terms
   real(kind=8),dimension(:),allocatable :: terms_dx
   real(kind=8),dimension(:),allocatable :: terms_dy
   real(kind=8),dimension(:),allocatable :: terms_dz
   ! Check section
   select case(allocated(this%extracoeff))
      case(.false.)
         write(0,*) "GET_ALLFUNCS_AND_DERIVS ERR: extra coefficients are not allocated"
         write(0,*) "GET_ALLFUNCS_AND_DERIVS ERR: did you use ADD_MOREFUNCS and INTERPOL before this?"
         call EXIT(1)
      case(.true.)
         ! do nothing
   end select
   nfuncs=size(this%extrafuncs(:,1))+1
   select case( size(f)/=nfuncs .or. size(dfdx(:,1))/=nfuncs .or. size(dfdx(1,:))/=3 )
      case(.true.)
         write(0,*) "GET_ALLFUNCS_AND DERIVS ERR: size mismatch of output arguments"
         write(0,*) "nfuncs: ",nfuncs
         write(0,*) "size f: ", size(f)
         write(0,*) "size dfdx: ",size(dfdx(:,1)),size(dfdx(1,:))
         call EXIT(1)
      case(.false.)
         ! do nothing
   end select
   ! Run section
   allocate(terms(this%n))
   allocate(terms_dx(this%n))
   allocate(terms_dy(this%n))
   allocate(terms_dz(this%n))
   do i = 1, this%n
      terms(i)=this%term%getValue( k=this%kListXY(i,1:2), parityXY=this%parityListXY(i),       irrepXY=this%irrepListXY(i),&
                                   l=this%kListAngle(i),  parityAngle=this%parityListAngle(i), irrepAngle=this%irrepListAngle(i),&
                                   x=x )
      terms_dx(i)=this%term%getDeriv1( k=this%kListXY(i,1:2), parityXY=this%parityListXY(i),       irrepXY=this%irrepListXY(i),&
                                       l=this%kListAngle(i),  parityAngle=this%parityListAngle(i), irrepAngle=this%irrepListAngle(i),&
                                       x=x )
      terms_dy(i)=this%term%getDeriv2( k=this%kListXY(i,1:2), parityXY=this%parityListXY(i),       irrepXY=this%irrepListXY(i),&
                                       l=this%kListAngle(i),  parityAngle=this%parityListAngle(i), irrepAngle=this%irrepListAngle(i),&
                                       x=x )
      terms_dz(i)=this%term%getDeriv3( k=this%kListXY(i,1:2), parityXY=this%parityListXY(i),       irrepXY=this%irrepListXY(i),&
                                       l=this%kListAngle(i),  parityAngle=this%parityListAngle(i), irrepAngle=this%irrepListAngle(i),&
                                       x=x )
   end do
   f(1)=dot_product(terms,this%coeff)
   dfdx(1,1)=dot_product(terms_dx,this%coeff)
   dfdx(1,2)=dot_product(terms_dy,this%coeff)
   dfdx(1,3)=dot_product(terms_dz,this%coeff)
   do i = 2, nfuncs
      f(i)=dot_product(terms,this%extracoeff(:,i-1))
      dfdx(i,1)=dot_product(terms_dx,this%extracoeff(:,i-1))
      dfdx(i,1)=dot_product(terms_dy,this%extracoeff(:,i-1))
      dfdx(i,1)=dot_product(terms_dz,this%extracoeff(:,i-1))
   end do
   return
end subroutine GET_ALLFUNC_AND_DERIVS_FOURIER3D
!###########################################################
!# SUBROUTINE: ADD_MORE_FUNCS_FOURIER3D
!###########################################################
!> @brief
!! Adds a new set of functions to interpolate at the same time
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine ADD_MORE_FUNCS_FOURIER3D(this,f)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(FOURIER3D),intent(inout) :: this
   real(kind=8),dimension(:,:),intent(in) :: f
   ! Local variables
   integer(kind=4) :: nfuncs, ndata
   ! Run section
   nfuncs=size(f(:,1)) ! number of rows
   ndata=size(f(1,:)) ! number of columns
   select case(ndata == this%n)
      case(.false.)
         write(0,*) "ADD_MORE_FUNCS_FOURIER3D ERR: size mismatch between extra functions and the original one"
         call EXIT(1)
      case(.true.)
         ! donothing
   end select
   allocate(this%extrafuncs(nfuncs,ndata))
   this%extrafuncs=f
   return
end subroutine ADD_MORE_FUNCS_FOURIER3D
!###########################################################
!# SUBROUTINE: INTERPOL_FOURIER3D
!###########################################################
!> @brief
!! Performs a generic FOURIER3D interpolation
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Mar/2014 
!> @version 1.0
!-----------------------------------------------------------
subroutine INTERPOL_FOURIER3D(this)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(FOURIER3D),intent(inout) :: this
   ! Local variables
   real(kind=8),dimension(:,:),allocatable :: terms,inv_terms
   integer(kind=4) :: i,j ! counters
   ! Run section
   allocate(this%coeff(this%n))
   allocate(terms(this%n,this%n))
   allocate(inv_terms(this%n,this%n))
   do i = 1, this%n ! loop over eq for different points
      do j = 1, this%n ! loop over coefficients
         terms(i,j)=this%term%getValue( k=this%kListXY(j,1:2), parityXY=this%parityListXY(j),       irrepXY=this%irrepListXY(j),&
                                        l=this%kListAngle(j),  parityAngle=this%parityListAngle(j), irrepAngle=this%irrepListAngle(j),&
                                        x=this%x(i,:) )
      end do
   end do
   call INV_MTRX(this%n,terms,inv_terms)
   this%coeff=matmul(inv_terms,this%f)
   deallocate(terms)
   ! Check if there are extra functions to be interpolated
   select case(allocated(this%extrafuncs))
      case(.true.)
         allocate(this%extracoeff(this%n,size(this%extrafuncs(:,1))))
         do i = 1, size(this%extrafuncs(:,1))
            this%extracoeff(:,i)=matmul(inv_terms,this%extrafuncs(i,:))
         end do
      case(.false.)
         ! do nothing
   end select
   return
end subroutine INTERPOL_FOURIER3D
!###########################################################
!# FUNCTION: getvalue_FOURIER3D
!###########################################################
!> @brief 
!! Get's fourier 3D interpolation value for a given point
!! inside the correct range
!-----------------------------------------------------------
function getValue_FOURIER3D(this,x) result(answer)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(Fourier3d),target,intent(in) :: this
   real(kind=8),dimension(3),intent(in) :: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Local variables
   integer(kind=4):: i ! counters
   real(kind=8),dimension(:),allocatable:: terms
   ! Run section
   allocate(terms(this%n))
   do i = 1, this%n
      terms(i)=this%term%getValue( k=this%kListXY(i,1:2), parityXY=this%parityListXY(i),       irrepXY=this%irrepListXY(i),&
                                   l=this%kListAngle(i),  parityAngle=this%parityListAngle(i), irrepAngle=this%irrepListAngle(i),&
                                   x=x )
   end do
   answer=dot_product(terms,this%coeff)
   deallocate(terms)
   return
end function getValue_FOURIER3D
!###########################################################
!# FUNCTION: getDeriv1_FOURIER3D
!###########################################################
!> @brief 
!! Get's derivative value at a given point X using the interpolation
!-----------------------------------------------------------
function getDeriv1_FOURIER3D(this,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier3d),target,intent(in) :: this
   real(kind=8),dimension(3),intent(in) :: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Local variables
   integer(kind=4):: i ! counters
   real(kind=8),dimension(:),allocatable:: terms
   ! Run section
   allocate(terms(this%n))
   do i = 1, this%n
      terms(i)=this%term%getDeriv1( k=this%kListXY(i,1:2), parityXY=this%parityListXY(i),       irrepXY=this%irrepListXY(i),&
                                    l=this%kListAngle(i),  parityAngle=this%parityListAngle(i), irrepAngle=this%irrepListAngle(i),&
                                    x=x )
   end do
   answer=dot_product(terms,this%coeff)
   deallocate(terms)
   return
end function getDeriv1_FOURIER3D
!###########################################################
!# FUNCTION: getDeriv2_FOURIER3D
!###########################################################
function getDeriv2_FOURIER3D(this,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier3d),target,intent(in) :: this
   real(kind=8),dimension(3),intent(in) :: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Local variables
   integer(kind=4):: i ! counters
   real(kind=8),dimension(:),allocatable:: terms
   ! Run section
   allocate(terms(this%n))
   do i = 1, this%n
      terms(i)=this%term%getDeriv1( k=this%kListXY(i,1:2), parityXY=this%parityListXY(i),  irrepXY=this%irrepListXY(i),&
                                    l=this%kListAngle(i),  parityAngle=this%parityListAngle(i), irrepAngle=this%irrepListAngle(i),&
                                    x=x )
   end do
   answer=dot_product(terms,this%coeff)
   deallocate(terms)
   return
end function getDeriv2_FOURIER3D
!###########################################################
!# FUNCTION: getDeriv3_FOURIER3D
!###########################################################
function getDeriv3_FOURIER3D(this,x) result(answer)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(Fourier3d),target,intent(in) :: this
   real(kind=8),dimension(3),intent(in) :: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Local variables
   integer(kind=4):: i ! counters
   real(kind=8),dimension(:),allocatable:: terms
   ! Run section
   allocate(terms(this%n))
   do i = 1, this%n
      terms(i)=this%term%getDeriv1( k=this%kListXY(i,1:2), parityXY=this%parityListXY(i),       irrepXY=this%irrepListXY(i),&
                                    l=this%kListAngle(i),  parityAngle=this%parityListAngle(i), irrepAngle=this%irrepListAngle(i),&
                                    x=x )
   end do
   answer=dot_product(terms,this%coeff)
   deallocate(terms)
   return
end function getDeriv3_FOURIER3D

end module FOURIER3D_MOD
