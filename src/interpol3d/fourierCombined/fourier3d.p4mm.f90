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
type,extends(TermCalculator3d) :: term3d_p4mm
private
   ! some atributes
contains
   procedure,public:: getValue => termFou3d_p4mm
   procedure,public:: getDeriv => termFou3d_dx_p4mm
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
   allocate(Term_4mm::this%term)
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
   real(kind=8),intent(in):: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou1d_4mm: '
   ! Run section
   answer=this%xyFourier%getValue()





   select case( irrep )
   case('A1')
      ! kpoint check
      select case( mod(kpoint,4)/=0 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('+')
         answer=dcos(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep "'//irrep//'" not implemented yet'
      write(0,*) 'Implemented ones: A1'
      call exit(1)
   end select

   return
end function termFou3d_p4mm
!###########################################################
!# FUNCTION: termfou1d_dx_4mm_A1
!###########################################################
!-----------------------------------------------------------
function termFou3d_dx_p4mm(this,kpoint,parity,irrep,x) result(realNum)
   ! Initial declarations 
   implicit none
   ! I/O variables
   class(Term3d_p4mm),intent(in) :: this
   integer(kind=4),intent(in) :: kpoint
   real(kind=8),intent(in) :: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: realNum
   ! Parameters
   character(len=*),parameter:: routineName='termfou1d_dx_4mm: '
    select case( irrep )
   case('A1')
      select case( parity )
      case('+')
         realNum=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep "'//irrep//'" not implemented yet'
      write(0,*) 'Implemented ones: A1'
      call exit(1)
   end select

   return
end function termFou3d_dx_p4mm
end module FOURIER3D_P4MM_MOD
