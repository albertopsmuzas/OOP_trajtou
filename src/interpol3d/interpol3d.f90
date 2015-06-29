!###################################################################################################
! MODULE: INTERPOL3D
!
!> @brief
!! Module that manages different interpolation schemes for one variable
!
!> @details
!! All types and subprograms intended to create interpolations in 1D should be placed inside this module
!##################################################################################################
module INTERPOL3D_MOD
! Initial declarations
#ifdef DEBUG
   use DEBUG_MOD
#endif
implicit none
!//////////////////////////////////////////////////////////////////////
! TYPE: interpol3d
!> @brief
!! Generic type of three dimensional interpolations
!-----------------------------------------------------------------------
type,abstract:: Interpol3d
   ! public atributes
   integer(kind=4),public:: n
   real(kind=8),dimension(:,:),allocatable,public:: x
   real(kind=8),dimension(:),allocatable,public:: f
   contains
      ! Initialization block
      procedure,non_overridable,public:: read => read_INTERPOL3D
      ! Get functions block
      procedure(getValue_INTERPOL3D),public,deferred:: getValue  ! child types, override this
      procedure(getValue_INTERPOL3D),public,deferred:: getDeriv1 ! child types, override this
      procedure(getValue_INTERPOL3D),public,deferred:: getDeriv2 ! child types, override this
      procedure(getValue_INTERPOL3D),public,deferred:: getDeriv3 ! child types, override this
end type Interpol3d
!//////////////////////////////////////////////////////////////////////
abstract interface
   !###########################################################
   !# FUNCTION: getvalue_interpol1d
   !###########################################################
   !> @brief
   !! Dummy function. Override it!!
   !-----------------------------------------------------------
   function getValue_INTERPOL3D(this,x,shift) result(answer)
      import Interpol3d
      class(Interpol3d),target,intent(in):: this
      real(kind=8),dimension(:,3),intent(in):: x
      real(kind=8):: answer
   end function getValue_INTERPOL3D
end interface
contains
!######################################################################
! SUBROUTINE: READ_INTERPOL3D #########################################
!######################################################################
!> @brief
!! Read main parameters for a 3D interpolation from arguments
!----------------------------------------------------------------------
subroutine read_INTERPOL3D(this,x,f)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Interpol3d),intent(inout):: this
   real(kind=8),dimension(:,:),intent(in)::x
   real(kind=8), dimension(:),intent(in):: f
   ! Run section -------
   if ( size(x(:,1))/=size(f) ) then
      write(0,*) "READ_INTERPOL3D: dimensions mismatch between x and f"
      call exit(1)
   end if
   this%n=size(x)
   allocate( this%x(this%n,3), source=x(:,:) )
   allocate( this%f(this%n),   source=f(:)   )
   return
end subroutine READ_INTERPOL3D
end module INTERPOL3D_MOD
