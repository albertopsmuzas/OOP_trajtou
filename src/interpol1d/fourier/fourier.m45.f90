!########################################################
! MODULE : FOURIER1D_M45_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes FOURIER1D_mod in its scope
!########################################################
MODULE FOURIER1D_M45_MOD
USE FOURIER1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: term_Ap
!> @brief
!! Child class of abstract termcalculator. Strcture to avoid unnecessary switches
!----------------------------------------------------------------
TYPE,EXTENDS(Termcalculator):: term_M45
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_M45
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_M45
END TYPE term_M45
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_M45
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
TYPE,EXTENDS(FOURIER1D):: Fourier1d_M45
   CONTAINS
      ! Set block
      PROCEDURE,PUBLIC :: initializeTerms => initializeTerms_FOURIER1D_M45
END TYPE FOURIER1D_M45
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_IRREP_FOURIER1D_M45 
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine initializeTerms_FOURIER1D_M45(this)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(Fourier1d_M45),intent(inout):: this
   ! Run section
   allocate(term_M45::this%term)
   return
end subroutine initializeTerms_FOURIER1D_M45
!###########################################################
!# FUNCTION: termfou1d_M45_Ap
!###########################################################
!-----------------------------------------------------------
function termfou1d_M45(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations 
   implicit none
   ! I/O variables
   class(term_M45),intent(in) :: this
   integer(kind=4),intent(in) :: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou_M45: '
   ! Run section
   select case( irrep )
   case('Ap')
      select case( parity )
      case('+')
         select case( mod(kpoint,4) )
         case(0)
            answer=dcos(dfloat(kpoint)*x)
         case default
            write(0,*) "termfou1d_m45 ERR: incorrect Kpoint-parity combination"
            call exit(1)
         end select

      case('-')
         select case( mod(kpoint,4) )
         case(2)
            answer=dsin(dfloat(kpoint)*x)
         case default
            write(0,*) "termfou1d_m45 ERR: incorrect point-parity combination"
            call exit(1)
         end select

      case('o')
         select case( mod(kpoint,4) )
         case(1)
            answer=dcos(dfloat(kpoint)*x)+dsin(dfloat(kpoint)*x)
         case(3)
            answer=dcos(dfloat(kpoint)*x)-dsin(dfloat(kpoint)*x)
         case default
            write(0,*) "termfou1d_m45 ERR: incorrect Kpoint-parity combination"
            call exit(1)
         end select

      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep "'//irrep//'" not implemented yet'
      write(0,*) 'Implemented ones: Ap'
      call exit(1)
   end select
   return
end function termfou1d_M45
!###########################################################
!# FUNCTION: termfou1d_dx_M45_Ap
!###########################################################
!-----------------------------------------------------------
function termfou1d_dx_M45(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations 
   implicit none
   ! I/O variables
   class(term_M45),intent(in) :: this
   integer(kind=4),intent(in) :: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou_dx_M45: '
   ! Run section
   select case( irrep )
   case('Ap')
      select case( parity )
      case('+')
         select case( mod(kpoint,4) )
         case(0)
            answer=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
         case default
            write(0,*) "termfou1d_m45 ERR: incorrect Kpoint-parity combination"
            call exit(1)
         end select

      case('-')
         select case( mod(kpoint,4) )
         case(2)
            answer=dfloat(kpoint)*dcos(dfloat(kpoint)*x)
         case default
            write(0,*) "termfou1d_m45 ERR: incorrect Kpoint-parity combination"
            call exit(1)
         end select

      case('o')
         select case( mod(kpoint,4) )
         case(1)
            answer=dfloat(kpoint)*( dcos(dfloat(kpoint)*x)-dsin(dfloat(kpoint)*x) )
         case(3)
            answer=-dfloat(kpoint)*( dcos(dfloat(kpoint)*x)+dsin(dfloat(kpoint)*x) )
         case default
            write(0,*) "termfou1d_m45 ERR: incorrect Kpoint-parity combination"
            call exit(1)
         end select

      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep "'//irrep//'" not implemented yet'
      write(0,*) 'Implemented ones: Ap'
      call exit(1)
   end select
   return
end function termfou1d_dx_M45
END MODULE FOURIER1D_M45_MOD
