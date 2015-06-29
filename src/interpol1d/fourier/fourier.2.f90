!########################################################
! MODULE : FOURIER1D_2_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes FOURIER1D_mod in its scope
!########################################################
MODULE FOURIER1D_2_MOD
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
TYPE,EXTENDS(Termcalculator) :: term_2
PRIVATE
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_2
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_2
END TYPE term_2
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_2
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
TYPE,EXTENDS(Fourier1d):: Fourier1d_2
   CONTAINS
      ! Set block
      PROCEDURE,PUBLIC :: initializeTerms => initializeTerms_FOURIER1D_2
END TYPE FOURIER1D_2
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: initializeTerm_FOURIER1D_2
!###########################################################
!> @brief
!! Sets Term for this fourier series
!-----------------------------------------------------------
subroutine initializeTerms_FOURIER1D_2(this)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(Fourier1d_2),intent(inout):: this
   ! Run section
   allocate( Term_2::this%term )
   return
end subroutine InitializeTerms_FOURIER1D_2
!###########################################################
!# FUNCTION: termfou1d_2
!###########################################################
!-----------------------------------------------------------
function termfou1d_2(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations 
   implicit none
   ! I/O variables
   class(Term_2),intent(in):: this
   integer(kind=4),intent(in):: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Local
   character(len=1):: auxChar
   ! Parameters
   character(len=*),parameter:: routinename='termfou_2: '
   ! Run section
   auxChar=trim(irrep)
   select case( auxChar )
   case('A')
      ! check parity
      select case( parity )
      case('+')
         select case( mod(kpoint,2)==0 )
         case(.true.) ! is par
            answer=dcos(dfloat(kpoint)*x)
         case(.false.) ! is odd
            write(0,*) 'ERR '//routinename//'bad combination of parity and Kpoint'
            call exit(1)
         end select
      case('-')
         select case( mod(kpoint,2)==0 )
         case(.true.) ! is par
            answer=dsin(dfloat(kpoint)*x)
         case(.false.) ! is odd
            write(0,*) 'ERR '//routinename//'bad combination of parity and Kpoint'
            call exit(1)
         end select
      case('o')
         select case( mod(kpoint,2)==0 )
         case(.true.) ! is par
            answer=dsin(dfloat(kpoint)*x)+dcos(dfloat(kpoint)*x)
         case(.false.) ! is odd
            write(0,*) 'ERR '//routinename//'bad combination of parity and Kpoint'
            call exit(1)
         end select
      case default
         write(0,*) 'ERR '//routinename//'bad parity symbol: "'//parity//'"'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep not implemented: "'//auxChar//'"'
      write(0,*) 'Implemented ones: A'
      call exit(1)
   end select

   return
end function termfou1d_2
!###########################################################
!# FUNCTION: termfou1d_dx_2
!###########################################################
!-----------------------------------------------------------
function termfou1d_dx_2(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations 
   implicit none
   ! I/O variables
   class(Term_2),intent(in):: this
   integer(kind=4),intent(in):: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Local
   character(len=1):: auxChar
   ! Parameters
   character(len=*),parameter:: routinename='termfou_dx_2: '
   ! Run section
   auxChar=trim(irrep)
   select case( auxChar )
   case('A')
      ! check parity
      select case( parity )
      case('+')
         select case( mod(kpoint,2)==0 )
         case(.true.) ! kpoint is par
            answer=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
         case(.false.) ! kpoint is odd
            write(0,*) 'ERR '//routinename//'bad combination of parity and Kpoint'
            call exit(1)
         end select
      case('-')
         select case( mod(kpoint,2)==0 )
         case(.true.) ! is par
            answer=dfloat(kpoint)*dcos(dfloat(kpoint)*x)
         case(.false.) ! is odd
            write(0,*) 'ERR '//routinename//'bad combination of parity and Kpoint'
            call exit(1)
         end select
      case('o')
         select case( mod(kpoint,2)==0 )
         case(.true.) ! is par
            answer=dfloat(kpoint)*( dcos(dfloat(kpoint)*x)-dsin(dfloat(kpoint)*x) )
         case(.false.) ! is odd
            write(0,*) 'ERR '//routinename//'bad combination of parity and Kpoint'
            call exit(1)
         end select
      case default
         write(0,*) 'ERR '//routinename//'bad parity symbol: "'//parity//'"'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep not implemented: "'//auxChar//'"'
      write(0,*) 'Implemented ones: A'
      call exit(1)
   end select

   return
end function termfou1d_dx_2
END MODULE FOURIER1D_2_MOD
