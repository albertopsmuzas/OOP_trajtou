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
TYPE,EXTENDS(Termcalculator) :: term_Ap
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_M45_Ap
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_M45_Ap
END TYPE term_Ap
!/////////////////////////////////////////////////////////////////
! TYPE: term_expanded_m45m1352_A1
!> @brief
!! Child class of abstract termcalculator. Strcture to avoid unnecessary switches
!----------------------------------------------------------------
TYPE,EXTENDS(Termcalculator) :: term_expanded_m45m1352_A1
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_expanded_m45m1352_A1
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_expanded_m45m1352_A1
END TYPE term_expanded_m45m1352_A1
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_M45
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
TYPE,EXTENDS(FOURIER1D):: Fourier1d_M45
   CONTAINS
      ! Set block
      PROCEDURE,PUBLIC :: SET_IRREP => SET_IRREP_FOURIER1D_M45
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
SUBROUTINE SET_IRREP_FOURIER1D_M45(this,irrep)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d_M45),INTENT(INOUT):: this
   CHARACTER(LEN=2),INTENT(IN) :: irrep
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   ! Run section
   SELECT CASE(irrep)
      CASE("Ap")
         ALLOCATE(term_Ap::this%term)
         this%irrep=irrep
         ALLOCATE(this%klist(this%n))
         FORALL(i=1:this%n) this%klist(i)=i-1
      CASE("A1") ! expanded symmetry
         ALLOCATE(term_expanded_m45m1352_A1::this%term)
         this%irrep=irrep
         ALLOCATE(this%klist(this%n))
         FORALL(i=1:this%n) this%klist(i)=(i-1)*2
      CASE DEFAULT
         WRITE(0,*) "SET_IRREP_FOURIER1D_M45 ERR: irrep used is not implemented or does not exist"
         WRITE(0,*) "List of irreps implemented: Ap, A1"
         WRITE(0,*) "Be careful with A1, it is an expanded symmetry flag. In this case we're using symmetry m45m1352_A1"
         CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE SET_IRREP_FOURIER1D_M45
!###########################################################
!# FUNCTION: termfou1d_expanded_m45m1352_A1 
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_expanded_m45m1352_A1(this,kpoint,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_expanded_m45m1352_A1),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   SELECT CASE(mod(kpoint,4))
      CASE(0)
         termfou1d_expanded_m45m1352_A1=dcos(dfloat(kpoint)*x)
      CASE(2)
         termfou1d_expanded_m45m1352_A1=dsin(dfloat(kpoint)*x)
      CASE DEFAULT
          WRITE(0,*) "termfou1d_expanded_m45m1352_A1 ERR: Unclassificable kpoint"
          CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfou1d_expanded_m45m1352_A1
!###########################################################
!# FUNCTION: termfou1d_dx_expanded_m45m1352_A1 
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_dx_expanded_m45m1352_A1(this,kpoint,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_expanded_m45m1352_A1),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   SELECT CASE(mod(kpoint,4))
      CASE(0)
         termfou1d_dx_expanded_m45m1352_A1=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
      CASE(2)
         termfou1d_dx_expanded_m45m1352_A1=dfloat(kpoint)*dcos(dfloat(kpoint)*x)
      CASE DEFAULT
          WRITE(0,*) "termfou1d_dx_expanded_m45m1352_A1 ERR: Unclassificable kpoint"
          CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfou1d_dx_expanded_m45m1352_A1
!###########################################################
!# FUNCTION: termfou1d_M45_Ap
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_M45_Ap(this,kpoint,x)
   ! Initial declarations 
   IMPLICIT NONE
   ! I/O variables
   CLASS(term_Ap),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   SELECT CASE(mod(kpoint,4))
      CASE(0)
         termfou1d_M45_Ap=dcos(dfloat(kpoint)*x)
      CASE(1)   
         termfou1d_M45_Ap=dcos(dfloat(kpoint)*x)+dsin(dfloat(kpoint)*x)
      CASE(2)   
         termfou1d_M45_Ap=dsin(dfloat(kpoint)*x)
      CASE(3)   
         termfou1d_M45_Ap=dcos(dfloat(kpoint)*x)-dsin(dfloat(kpoint)*x)
      CASE DEFAULT
          WRITE(0,*) "termfou1d_m45_Ap ERR: Unclassificable kpoint"
          CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfou1d_M45_Ap
!###########################################################
!# FUNCTION: termfou1d_dx_M45_Ap
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_dx_M45_Ap(this,kpoint,x)
   ! Initial declarations 
   IMPLICIT NONE
   ! I/O variables
   CLASS(term_Ap),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   SELECT CASE(mod(kpoint,4))
      CASE(0)
         termfou1d_dx_M45_Ap=-dsin(dfloat(kpoint)*x)*dfloat(kpoint)
      CASE(1)   
         termfou1d_dx_M45_Ap=(-dsin(dfloat(kpoint)*x)+dcos(dfloat(kpoint)*x))*dfloat(kpoint)
      CASE(2)   
         termfou1d_dx_M45_Ap=-dcos(dfloat(kpoint)*x)*dfloat(kpoint)
      CASE(3)   
         termfou1d_dx_M45_Ap=(-dsin(dfloat(kpoint)*x)-cos(dfloat(kpoint)*x))*dfloat(kpoint)
      CASE DEFAULT
          WRITE(0,*) "termfou1d_m45_Ap ERR: Unclassificable kpoint"
          CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfou1d_dx_M45_Ap
END MODULE FOURIER1D_M45_MOD
