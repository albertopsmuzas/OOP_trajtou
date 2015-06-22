!########################################################
! MODULE : FOURIER1D_MM2_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes FOURIER1D_mod in its scope
!########################################################
MODULE FOURIER1D_MM2_MOD
USE FOURIER1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: term_A1
!> @brief
!! Child class of abstract termcalculator. Strcture to avoid unnecessary switches
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Termcalculator) :: term_mm2
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_MM2
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_MM2
END TYPE term_mm2
!///////////////////////////////////////////////////////////////////////////////
! TYPE: FOURIER1D_MM2
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!-------------------------------------------------------------------------------
TYPE,EXTENDS(FOURIER1D):: Fourier1d_MM2
   CONTAINS
      PROCEDURE,PUBLIC :: initializeTerms => initializeTerms_FOURIER1D_MM2
END TYPE FOURIER1D_MM2
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: initializeTerms_FOURIER1D_MM2
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE initializeTerms_FOURIER1D_MM2(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d_MM2),INTENT(INOUT):: this
   ! Run section
   ALLOCATE(Term_mm2::this%term)
   RETURN
END SUBROUTINE initializeTerms_FOURIER1D_MM2
!###########################################################
!# FUNCTION: termfou1d_MM2_A1
!###########################################################
!-----------------------------------------------------------
function termfou1d_MM2(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations 
   implicit none
   ! I/O variables
   class(Term_mm2),intent(in):: this
   integer(kind=4),intent(in):: kpoint
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   real(kind=8),intent(in):: x
   ! Dummy out variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou1d_MM2: '
   ! Run section
   select case( irrep )
   case('A1')
      ! check KPOINT
      select case( mod(kpoint,2)==1 )
      case(.true.)
         write(0,*) 'ERR '//routinename//'bad irrep'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! check parity
      select case( parity )
      case('+')
         answer=dcos(dfloat(kpoint)*x)
      case default
         write(0,*) 'ERR '//routinename//'bad parity'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep not implemented'
      call exit(1)
   end select
   return
end function termfou1d_MM2
!###########################################################
!# FUNCTION: termfou1d_dx_MM2_A1
!###########################################################
!-----------------------------------------------------------
function termfou1d_dx_MM2(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term_mm2),intent(in):: this
   integer(kind=4),intent(in):: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy out variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou1d_MM2: '
   ! Run section
   select case( irrep )
   case('A1')
      ! check KPOINT
      select case( mod(kpoint,2)==1 )
      case(.true.)
         write(0,*) 'ERR '//routinename//'bad irrep'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! check parity
      select case( parity )
      case('+')
         answer=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
      case default
         write(0,*) 'ERR '//routinename//'bad parity'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep not implemented'
      call exit(1)
   end select
   return
end function termfou1d_dx_MM2

END MODULE FOURIER1D_MM2_MOD
