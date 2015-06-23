!########################################################
! MODULE : FOURIER1D_E_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes FOURIER1D_mod in its scope
!########################################################
MODULE FOURIER1D_E_MOD
USE FOURIER1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: term_A
!> @brief
!! Child class of abstract termcalculator. Strcture to avoid unnecessary switches
!----------------------------------------------------------------
TYPE,EXTENDS(Termcalculator) :: term_E
PRIVATE
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_E
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_E
END TYPE term_E
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_E
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
TYPE,EXTENDS(FOURIER1D):: Fourier1d_E
   CONTAINS
      ! Set block
      PROCEDURE,PUBLIC :: initializeTerms => initializeTerms_FOURIER1D_E
END TYPE FOURIER1D_E
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: initializeTerms_FOURIER1D_E
!###########################################################
!> @brief
!! Sets Terms for this fourier series
!-----------------------------------------------------------
subroutine initializeTerms_FOURIER1D_E(this)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(Fourier1d_E),intent(inout):: this
   ! Run section
   allocate(Term_E::this%term)
   return
end subroutine initializeTerms_FOURIER1D_E
!###########################################################
!# FUNCTION: termfou1d_E
!###########################################################
!-----------------------------------------------------------
function termfou1d_E(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations 
   implicit none
   ! I/O variables
   class(Term_E),intent(in) :: this
   integer(kind=4),intent(in) :: kpoint
   real(kind=8),intent(in) :: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou_E: '
   ! Run section
   select case( irrep )
   case('A')
      ! check parity
      select case( parity )
      case('+')
            answer=dcos(dfloat(kpoint)*x)
      case('-')
            answer=dsin(dfloat(kpoint)*x)
      case('o')
            answer=dsin(dfloat(kpoint)*x)+dcos(dfloat(kpoint)*x)
      case default
         write(0,*) 'ERR '//routinename//'bad parity symbol'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep not implemented'
      write(0,*) 'Implemented ones: A'
      call exit(1)
   end select

   return
end function termfou1d_E
!###########################################################
!# FUNCTION: termfou1d_dx_E
!###########################################################
!-----------------------------------------------------------
function termfou1d_dx_E(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations 
   implicit none
   ! I/O variables
   class(Term_E),intent(in) :: this
   integer(kind=4),intent(in) :: kpoint
   real(kind=8),intent(in) :: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou_dx_E: '
   ! Run section
   select case( irrep )
   case('A')
      ! check parity
      select case( parity )
      case('+')
            answer=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
      case('-')
            answer=dfloat(kpoint)*dcos(dfloat(kpoint)*x)
      case('o')
            answer=dfloat(kpoint)*( dcos(dfloat(kpoint)*x)-dsin(dfloat(kpoint)*x) )
      case default
         write(0,*) 'ERR '//routinename//'bad parity symbol'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep not implemented'
      write(0,*) 'Implemented ones: A'
      call exit(1)
   end select

   return
end function termfou1d_dx_E
END MODULE FOURIER1D_E_MOD
