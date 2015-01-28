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
TYPE,EXTENDS(Termcalculator) :: term_A1
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_MM2_A1
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_MM2_A1
END TYPE term_A1
!///////////////////////////////////////////////////////////////////////////////
! TYPE: FOURIER1D_MM2
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!-------------------------------------------------------------------------------
TYPE,EXTENDS(FOURIER1D):: Fourier1d_MM2
   CONTAINS
      PROCEDURE,PUBLIC :: SET_IRREP => SET_IRREP_FOURIER1D_MM2
END TYPE FOURIER1D_MM2
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_IRREP_FOURIER1D_MM2 
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_IRREP_FOURIER1D_MM2(this,irrep)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d_MM2),INTENT(INOUT):: this
   CHARACTER(LEN=2),INTENT(IN) :: irrep
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   ! Run section
   SELECT CASE(irrep)
      CASE("A1")
         ALLOCATE(Term_A1::this%term)
         this%irrep=irrep
         ALLOCATE(this%klist(this%n))
         FORALL(i=1:this%n) this%klist(i)=(i-1)*2
         
      CASE DEFAULT
         WRITE(0,*) "SET_IRREP_FOURIER1D_MM2 ERR: irrep used is not implemented or does not exist"
         WRITE(0,*) "List of irreps implemented: A1"
         CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE SET_IRREP_FOURIER1D_MM2
!###########################################################
!# FUNCTION: termfou1d_MM2_A1
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_MM2_A1(this,kpoint,x)
   ! Initial declarations 
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_A1),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   termfou1d_MM2_A1=dcos(dfloat(kpoint)*x)
   RETURN
END FUNCTION termfou1d_MM2_A1
!###########################################################
!# FUNCTION: termfou1d_dx_MM2_A1
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_dx_MM2_A1(this,kpoint,x)
   ! Initial declarations 
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_A1),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   termfou1d_dx_MM2_A1=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
   RETURN
END FUNCTION termfou1d_dx_MM2_A1
END MODULE FOURIER1D_MM2_MOD
