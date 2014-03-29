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
!> @param xy(:,:) - Matrix that collects @f$(x_{i},y_{i})@f$ pairs
!> @param x(:) - Grid in X. Only if input has grid structure
!> @param y(:) - Grid in Y. Only if input has grid structure
!> @param fgrid(:,:) - Function evaluated in a grid
!> @param f - Array that stores couples @f$F(x_{i},y_{i})@f$. Non grid input.
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
TYPE,ABSTRACT :: Interpol2d
   INTEGER(KIND=4) :: n
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: x,y
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: fgrid
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: xy
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: f
CONTAINS
   PROCEDURE,PUBLIC :: READ => READ_INTERPOL2D
   PROCEDURE,PUBLIC :: READGRID => READGRID_INTERPOL2D
   PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_INTERPOL2D
END TYPE Interpol2d
!////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: INTERPOL_INTERPOL2D 
!###########################################################
!> @brief
!! Dummy subroutine to be overriden
!-----------------------------------------------------------
SUBROUTINE INTERPOL_INTERPOL2D(this,surf,filename)
   ! Initial declarations   
   USE SURFACE_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(Interpol2d),INTENT(INOUT) :: this
   TYPE(Surface),INTENT(IN),OPTIONAL :: surf
   CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
   ! Run section
   WRITE(0,*) "INTERPOL_INTERPOL2D ERR: You didn't allocate this variable with the proper type"
   CALL EXIT(1)
   RETURN
END SUBROUTINE INTERPOL_INTERPOL2D
!##########################################################
! SUBROUTINE: READGRID_INTERPOL2D
!> @brief
!! Reads input defined in a grid
!
!> @param[out] this - interpol 2D object to be set up
!> @param[in] x(:) - X grid
!> @param[in] y(:) - Y grid
!> @param[in] f(:,:) - stores F falues for each point in the grid
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 17/Feb/2014
!> @version 1.0 
!----------------------------------------------------------
SUBROUTINE READGRID_INTERPOL2D(this,x,y,f)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Interpol2d),INTENT(OUT) :: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x,y
   REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: f
   ! Local variables
   INTEGER(KIND=4) :: nx,ny
   INTEGER(KIND=4) :: i,j ! counters
   ! Run section ------------------------------------
   nx=size(x)
   ny=size(y)
   IF ((nx/=size(f(:,1))).OR.(ny/=size(f(1,:)))) THEN
      WRITE(0,*) "READGRID_INTERPOL2D: array mismatch x, y, fgrif"
      CALL EXIT(1)
   END IF
   ALLOCATE(this%x(nx))
   ALLOCATE(this%y(ny))
   ALLOCATE(this%fgrid(nx,ny))
   this%x = x
   this%y = y
   this%fgrid = f
   RETURN
END SUBROUTINE READGRID_INTERPOL2D
!###########################################################
!# SUBROUTINE: READ_INTERPOL2D 
!###########################################################
!> @brief
!! Reads xy and f values from arguments
!
!> @param[out] this - Interpol2d class to be read
!> @param[in] xy - Matrix that collects @f$(x_{i},y_{i})@f$ pairs
!> @param[in] f - Values of @f$F(x_{i},y_{i})@f$i. It is a matrix so that
!!                the user can provide in each row a different function that
!!                will have the same @f$T^{-1}@f$ matrix during the
!!                interpolation
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Mar/2014
!> @version 1.1 
!
!> @see interpol2d
!-----------------------------------------------------------
SUBROUTINE READ_INTERPOL2D(this,xy,f)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Interpol2d),INTENT(OUT) :: this
   REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: xy
   REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: f
   ! Local variables
   INTEGER(KIND=4) :: n,q
   ! Run section
   n=size(f(1,:))
   q=size(f(:,1))
   SELECT CASE(size(xy(:,1)) == n .OR. size(xy(1,:))==2)
      CASE(.FALSE.)
         WRITE(0,*) "READ_INTERPOL2D ERR: dimensions mismatch in arrays xy or f"
         CALL EXIT(1)
      CASE(.TRUE.)
         ! do nothing
   END SELECT
   this%n = n
   ALLOCATE(this%xy(n,2))
   ALLOCATE(this%f(n,q))
   this%f = f
   this%xy = xy
   RETURN
END SUBROUTINE READ_INTERPOL2D  
END MODULE INTERPOL2D_MOD
