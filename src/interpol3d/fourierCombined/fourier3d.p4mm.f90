!########################################################
! MODULE : FOURIER3D_P4MM_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes FOURIER1D_mod in its scope
!########################################################
module FOURIER3D_P4MM_MOD
use FOURIER3D_MOD, only: TermCalculator3d,Fourier3d
use FOURIER1D_4MM_MOD, only: Fourier1d_4mm
use FOURIER_P4MM_MOD, only: Fourierp4mm
implicit none
!/////////////////////////////////////////////////////////////////
! TYPE: term_A1
!> @brief
!! Child class of abstract termcalculator. Structure to avoid unnecessary switches
!----------------------------------------------------------------
type,extends(TermCalculator3d) :: Term3d_p4mm
private
   ! some atributes
contains
   procedure,public:: getValue  => termFou3d_p4mm
   procedure,public:: getDeriv1 => termFou3d_dx_p4mm
   procedure,public:: getDeriv2 => termFou3d_dy_p4mm
   procedure,public:: getDeriv3 => termFou3d_dz_p4mm
end type term3d_p4mm
!/////////////////////////////////////////////////
! TYPE: FOURIER3D_P4MM
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
type,extends(Fourier3d):: Fourier3d_p4mm
   contains
      ! Set block
      procedure,public:: initializeTerms => initializeTerms_FOURIER3D_P4MM
end type Fourier3d_p4mm
!//////////////////////////////////////////////////
contains
!###########################################################
!# SUBROUTINE: SET_IRREP_FOURIER3D_P4MM
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!-----------------------------------------------------------
subroutine initializeTerms_FOURIER3D_P4MM(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier3d_p4mm),intent(inout):: this
   ! Run section
   allocate( Term3d_p4mm::   this%term )
   allocate( Fourierp4mm::   this%term%xyFourier )
   allocate( Fourier1d_4mm:: this%term%angleFourier )
   return
end subroutine initializeTerms_FOURIER3D_P4MM
!###########################################################
!# FUNCTION: termfou1d_4mm
!###########################################################
!-----------------------------------------------------------
function termFou3d_p4mm(this,k,parityXY,irrepXY,l,parityAngle,irrepAngle,x) result(answer)
   ! Initial declarations 
   implicit none
   ! I/O variables
   class(Term3d_p4mm),intent(in):: this
   integer(kind=4),dimension(2),intent(in):: k
   integer(kind=4),intent(in):: l
   character(len=1),intent(in):: parityXY
   character(len=1),intent(in):: parityAngle
   character(len=2),intent(in):: irrepXY
   character(len=2),intent(in):: irrepAngle
   real(kind=8),dimension(3),intent(in):: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou1d_4mm: '
   ! Run section
   answer=this%xyFourier%term%getValue( x=x(1:2),k=k,irrep=irrepXY,parity=parityXY )*&
          this%angleFourier%term%getValue( x=x(3),kpoint=l,irrep=irrepAngle,parity=parityAngle )
   return
end function termFou3d_p4mm
!###########################################################
!# FUNCTION: termfou1d_dx_4mm_A1
!###########################################################
!-----------------------------------------------------------
function termFou3d_dx_p4mm(this,k,parityXY,irrepXY,l,parityAngle,irrepAngle,x) result(answer)
   ! Initial declarations 
   implicit none
   ! I/O variables
   class(Term3d_p4mm),intent(in):: this
   integer(kind=4),dimension(2),intent(in):: k
   integer(kind=4),intent(in):: l
   character(len=1),intent(in):: parityXY
   character(len=1),intent(in):: parityAngle
   character(len=2),intent(in):: irrepXY
   character(len=2),intent(in):: irrepAngle
   real(kind=8),dimension(3),intent(in):: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou1d_4mm: '
   ! Run section
   answer=this%xyFourier%term%getDeriv1( x=x(1:2),k=k,irrep=irrepXY,parity=parityXY )*&
          this%angleFourier%term%getValue( x=x(3),kpoint=l,irrep=irrepAngle,parity=parityAngle )
   return
end function termFou3d_dx_p4mm
!###########################################################
!# FUNCTION: termfou1d_dy_4mm_A1
!###########################################################
!-----------------------------------------------------------
function termFou3d_dy_p4mm(this,k,parityXY,irrepXY,l,parityAngle,irrepAngle,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term3d_p4mm),intent(in):: this
   integer(kind=4),dimension(2),intent(in):: k
   integer(kind=4),intent(in):: l
   character(len=1),intent(in):: parityXY
   character(len=1),intent(in):: parityAngle
   character(len=2),intent(in):: irrepXY
   character(len=2),intent(in):: irrepAngle
   real(kind=8),dimension(3),intent(in):: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou1d_4mm: '
   ! Run section
   answer=this%xyFourier%term%getDeriv2( x=x(1:2),k=k,irrep=irrepXY,parity=parityXY )*&
          this%angleFourier%term%getValue( x=x(3),kpoint=l,irrep=irrepAngle,parity=parityAngle )
   return
end function termFou3d_dy_p4mm
!###########################################################
!# FUNCTION: termfou1d_dz_4mm_A1
!###########################################################
!-----------------------------------------------------------
function termFou3d_dz_p4mm(this,k,parityXY,irrepXY,l,parityAngle,irrepAngle,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term3d_p4mm),intent(in):: this
   integer(kind=4),dimension(2),intent(in):: k
   integer(kind=4),intent(in):: l
   character(len=1),intent(in):: parityXY
   character(len=1),intent(in):: parityAngle
   character(len=2),intent(in):: irrepXY
   character(len=2),intent(in):: irrepAngle
   real(kind=8),dimension(3),intent(in):: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou1d_4mm: '
   ! Run section
   answer=this%xyFourier%term%getValue( x=x(1:2),k=k,irrep=irrepXY,parity=parityXY )*&
          this%angleFourier%term%getDeriv( x=x(3),kpoint=l,irrep=irrepAngle,parity=parityAngle )
   return
end function termFou3d_dz_p4mm

end module FOURIER3D_P4MM_MOD
