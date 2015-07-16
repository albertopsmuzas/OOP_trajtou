!########################################################
! MODULE : FOURIER1D_4MM_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes FOURIER1D_mod in its scope
!########################################################
MODULE FOURIER1D_4MM_MOD
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
TYPE,EXTENDS(Termcalculator) :: term_4mm
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_4mm
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_4mm
END TYPE term_4mm
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_4MM
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
type,extends(Fourier1d):: Fourier1d_4mm
   contains
      ! Set block
      procedure,public:: initializeTerms => initializeTerms_FOURIER1D_4MM
end type Fourier1d_4mm
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_IRREP_FOURIER1D_4MM
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!-----------------------------------------------------------
subroutine initializeTerms_FOURIER1D_4MM(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier1d_4mm),intent(inout):: this
   ! Run section
   allocate(Term_4mm::this%term)
   return
end subroutine initializeTerms_FOURIER1D_4MM
!###########################################################
!# FUNCTION: termfou1d_4mm
!###########################################################
!-----------------------------------------------------------
function termfou1d_4mm(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations 
   implicit none
   ! I/O variables
   class(Term_4mm),intent(in):: this
   integer(kind=4),intent(in):: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou1d_4mm: '
   ! Run section
   select case( irrep )
   case('A1')
      ! kpoint check
      select case( mod(kpoint,4)/=0 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
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

   case('A2')
      ! kpoint check
      select case( mod(kpoint,4)/=0 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('-')
         answer=dsin(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case('B1')
      ! kpoint check
      select case( mod(kpoint,4)/=2 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
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

   case('B2')
      ! kpoint check
      select case( mod(kpoint,4)/=2 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('-')
         answer=dsin(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case('E','EE')
      ! kpoint check
      select case( mod(kpoint,2)==0 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('+')
         answer=dcos(dfloat(kpoint)*x)
      case('-')
         answer=dsin(dfloat(kpoint)*x)
      case('o')
         answer=dsin(dfloat(kpoint)*x)+dcos(dfloat(kpoint)*x)
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
end function termfou1d_4mm
!###########################################################
!# FUNCTION: termfou1d_dx_4mm_A1
!###########################################################
!-----------------------------------------------------------
function termfou1d_dx_4mm(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations 
   implicit none
   ! I/O variables
   class(Term_4mm),intent(in) :: this
   integer(kind=4),intent(in) :: kpoint
   real(kind=8),intent(in) :: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routineName='termfou1d_dx_4mm: '
    select case( irrep )
   case('A1')
      select case( parity )
      case('+')
         answer=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case('A2')
      ! kpoint check
      select case( mod(kpoint,4)/=0 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('-')
         answer=dfloat(kpoint)*dcos(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case('B1')
      ! kpoint check
      select case( mod(kpoint,4)/=2 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('+')
         answer=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case('B2')
      ! kpoint check
      select case( mod(kpoint,4)/=2 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('-')
         answer=dfloat(kpoint)*dcos(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case('E','EE')
      ! kpoint check
      select case( mod(kpoint,2)==0 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('+')
         answer=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
      case('-')
         answer=dfloat(kpoint)*dcos(dfloat(kpoint)*x)
      case('o')
         answer=dfloat(kpoint)*( dcos(dfloat(kpoint)*x)-dsin(dfloat(kpoint)*x) )
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
end function termfou1d_dx_4mm
END MODULE FOURIER1D_4MM_MOD
