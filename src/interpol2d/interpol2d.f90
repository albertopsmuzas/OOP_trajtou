!##################################################################################
! MODULE: INTERPOL2D_MOD
!> @brief
!! Provides tools to perform 2D interpolations
!##################################################################################
MODULE INTERPOL2D_MOD
IMPLICIT NONE
!////////////////////////////////////////////////////////////////
! TYPE: Interpol2d
!
!> @brief
!! Generic 2D interpolation type variable
!
!> @param n - Number of data points
!> @param xy - Matrix that collects @f$(x_{i},y_{i})@f$ pairs
!> @param f - Array that stores couples @f$F(x_{i},y_{i})@f$
!> @param dfdz - Array that stores couples @f$\frac{\partial F(x_{i},y_{i})}{\partial z}@f$
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 06/Feb/2014
!> @version 1.0 
!
!> @warning
!! - We are assuming that there is a third variable Z in which a derivative
!!   can be defined and should be given
!
!> @todo 
!! - Generalize to "n" extra variables, apart from x and y, not only Z
!---------------------------------------------------------------
TYPE Interpol2d
   INTEGER(KIND=4) :: n
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: xy
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: f, dfdz
CONTAINS
   PROCEDURE,PUBLIC :: READ => READ_INTERPOL2D
END TYPE Interpol2d
!////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: READ_INTERPOL2D 
!###########################################################
!> @brief
!! Reads xy and f values from arguments
!
!> @param[out] this - Interpol2d class to be read
!> @param[in] xy - Matrix that collects @f$(x_{i},y_{i})@f$ pairs
!> @param[in] f - Values of @f$F(x_{i},y_{i})@f$
!> @param[in] dfdz - Values of @f$\frac{\partial F(x_{i},y_{i}}{\partial z})@f$
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 06/Feb/2014
!> @version 1.0 
!
!> @see interpol2d
!-----------------------------------------------------------
SUBROUTINE READ_INTERPOL2D(this,xy,f,dfdz)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Interpol2d),INTENT(OUT) :: this
   REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: xy
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: f
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: dfdz
   ! Local variables
   INTEGER(KIND=4) :: n
   ! Run section
   n=size(f)
   IF ((size(xy(:,1))/=n).OR.(size(dfdz)/=n)) THEN
      WRITE(0,*) "READ_INTERPOL2D ERR: dimensions mismatch between arrays xy, f and dfdz"
      CALL EXIT(1)
   END IF
   IF (size(xy(1,:))/=2) THEN
      WRITE(0,*) "READ_INTERPOL2D ERR: xy should have 2 colummns"
      CALL EXIT(1)
   END IF
   this%n = n
   ALLOCATE(this%xy(n,2))
   ALLOCATE(this%f(n))
   ALLOCATE(this%dfdz(n))
   this%f = f
   this%xy = xy
   this%dfdz = dfdz
   RETURN
END SUBROUTINE READ_INTERPOL2D  
END MODULE INTERPOL2D_MOD
