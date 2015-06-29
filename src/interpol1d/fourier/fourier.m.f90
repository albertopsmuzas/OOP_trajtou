!########################################################
! MODULE : FOURIER1D_M_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes FOURIER1D_mod in its scope
!########################################################
MODULE FOURIER1D_M_MOD
USE FOURIER1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: term_Ap
!> @brief
!! Child class of abstract termcalculator. Strcture to avoid unnecessary switches
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Termcalculator) :: term_M
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_M
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_M
END TYPE term_M
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_M
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
TYPE,EXTENDS(FOURIER1D):: Fourier1d_M
   CONTAINS
      ! Set block
      PROCEDURE,PUBLIC :: initializeTerms => initializeTerms_FOURIER1D_M
END TYPE FOURIER1D_M
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_IRREP_FOURIER1D_M 
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine initializeTerms_FOURIER1D_M(this)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(Fourier1d_M),intent(inout):: this
   ! Run section
   allocate(term_M::this%term)
   return
end subroutine initializeTerms_FOURIER1D_M
!###########################################################
!# FUNCTION: termfou1d_M
!###########################################################
!-----------------------------------------------------------
function termfou1d_M(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations 
   implicit none
   ! I/O variables
   class(term_M),intent(in):: this
   integer(kind=4),intent(in):: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou_M: '
   ! Run section
   select case( irrep )
   case('Ap')
      select case( parity )
      case('+')
         answer=dcos(dfloat(kpoint)*(x+this%getShift()))
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
end function termfou1d_M
!###########################################################
!# FUNCTION: termfou1d_dx_M
!###########################################################
!-----------------------------------------------------------
function termfou1d_dx_M(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations 
   implicit none
   ! I/O variables
   class(term_M),intent(in):: this
   integer(kind=4),intent(in):: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou_M: '
   ! Run section
   select case( irrep )
   case('Ap')
      select case( parity )
      case('+')
         answer=-dsin(dfloat(kpoint)*(x+this%getshift()))*dfloat(kpoint)
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
end function termfou1d_dx_M

END MODULE FOURIER1D_M_MOD
