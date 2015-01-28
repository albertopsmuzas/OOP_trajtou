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
TYPE,EXTENDS(Termcalculator) :: term_A
PRIVATE
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_2_A
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_2_A
END TYPE term_A
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_2
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
TYPE,EXTENDS(Fourier1d):: Fourier1d_2
   CONTAINS
      ! Set block
      PROCEDURE,PUBLIC :: SET_IRREP => SET_IRREP_FOURIER1D_2
END TYPE FOURIER1D_2
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_IRREP_FOURIER1D_2 
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_IRREP_FOURIER1D_2(this,irrep)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d_2),INTENT(INOUT):: this
   CHARACTER(LEN=2),INTENT(IN) :: irrep
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   INTEGER(KIND=4) :: npar
   ! Run section
   ALLOCATE(this%klist(this%n))
   SELECT CASE(irrep)
      CASE("A")
         ALLOCATE(Term_A::this%term)
         this%irrep=irrep
         SELECT CASE(mod(this%n,2) == 0) 
            CASE(.TRUE.) ! case is even
               CALL this%SET_AVERAGE_LASTKPOINT(.TRUE.)
               this%klist(1)=0
               npar=(this%n-2)/2
                DO i = 1, npar 
                   this%klist(i+1)=2*i
                   this%klist(i+1+npar)=-2*i
                END DO
                this%klist(this%n)=npar+1
                CALL this%SET_LASTKPOINT(npar+1)
            CASE(.FALSE.) ! case is odd
               this%klist(1)=0
               npar=(this%n-1)/2
                DO i = 1, npar 
                   this%klist(i+1)=2*i
                   this%klist(i+1+npar)=-2*i
                END DO
                CALL this%SET_LASTKPOINT(npar)
         END SELECT
      CASE DEFAULT
         WRITE(0,*) "SET_IRREP_FOURIER1D_2 ERR: irrep used is not implemented or does not exist"
         WRITE(0,*) "List of irreps implemented: A"
         CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE SET_IRREP_FOURIER1D_2
!###########################################################
!# FUNCTION: termfou1d_2_A1
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_2_A(this,kpoint,x)
   ! Initial declarations 
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_A),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   SELECT CASE(this%getaveragelast() .AND. kpoint == this%getlastkpoint())
      CASE(.TRUE.)
         termfou1d_2_A=dcos(dfloat(kpoint)*x)+dsin(dfloat(kpoint)*x)
      CASE(.FALSE.)
         SELECT CASE(kpoint)
            CASE(0)
               termfou1d_2_A=1.D0
            CASE(: -1)
               termfou1d_2_A=dsin(dfloat(-kpoint)*x)
            CASE(1 :)
               termfou1d_2_A=dcos(dfloat(kpoint)*x)
            CASE DEFAULT
               WRITE(0,*) "Termfou1d_2_A ERR: Something went really wrong with kpoints of this interpolation"
               CALL EXIT(1)
         END SELECT
   END SELECT
   RETURN
END FUNCTION termfou1d_2_A
!###########################################################
!# FUNCTION: termfou1d_dx_2_A1
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_dx_2_A(this,kpoint,x)
   ! Initial declarations 
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_A),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   SELECT CASE(this%getaveragelast() .AND. kpoint == this%getlastkpoint())
      CASE(.TRUE.)
         termfou1d_dx_2_A=(-dsin(dfloat(kpoint)*x)+dcos(dfloat(kpoint)*x))*dfloat(kpoint)
      CASE(.FALSE.)
         SELECT CASE(kpoint)
            CASE(0)
               termfou1d_dx_2_A=0.D0
            CASE(: -1)
               termfou1d_dx_2_A=dfloat(-kpoint)*dcos(dfloat(-kpoint)*x)
            CASE(1 :)
               termfou1d_dx_2_A=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
            CASE DEFAULT
               WRITE(0,*) "Termfou1d_2_A ERR: Something went really wrong with kpoints of this interpolation"
               CALL EXIT(1)
         END SELECT
   END SELECT
   RETURN
END FUNCTION termfou1d_dx_2_A
END MODULE FOURIER1D_2_MOD
