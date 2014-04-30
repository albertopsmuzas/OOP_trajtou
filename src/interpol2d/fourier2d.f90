!##################################################################################
! MODULE: FOURIER2D_MOD
!> @brief
!! Provides tools to perform 2D interpolations with fourier series
!##################################################################################
MODULE FOURIER2D_MOD
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
!---------------------------------------------------------------
TYPE :: Fourier2d
   INTEGER(KIND=4) :: n
   INTEGER(KIND=4) :: nfunc
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: xy
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: f
   INTEGER(KIND=4),DIMENSION(:,:),ALLOCATABLE :: klist
CONTAINS
   PROCEDURE,PUBLIC :: READ => READ_FOURIER2D
   PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_FOURIER2D
END TYPE Fourier2d
!////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: INTERPOL_FOURIER2D 
!###########################################################
!> @brief
!! Dummy subroutine to be overriden
!-----------------------------------------------------------
SUBROUTINE INTERPOL_FOURIER2D(this,surf,filename)
   ! Initial declarations   
   USE SURFACE_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier2d),INTENT(INOUT) :: this
   TYPE(Surface),INTENT(IN) :: surf
   CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
   ! Run section
   WRITE(0,*) "INTERPOL_FOURIER2D ERR: You didn't allocate this variable with the proper type"
   CALL EXIT(1)
   RETURN
END SUBROUTINE INTERPOL_FOURIER2D
!###########################################################
!# SUBROUTINE: READ_FOURIER2D 
!###########################################################
!> @brief
!! Reads xy, f and klist values from arguments
!
!> @param[out] this - Interpol2d class to be read
!> @param[in] xy - Matrix that collects @f$(x_{i},y_{i})@f$ pairs
!> @param[in] f - Values of @f$F(x_{i},y_{i})@f$. It is a matrix so that
!!                the user can provide in each row a different function that
!!                will have the same @f$T^{-1}@f$ matrix during the
!!                interpolation
!> @param[in] klist - Kpoints to be used in the expansion. There should be as many
!!                    of them as numbers of evaluations of f
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Mar/2014
!> @version 1.0 
!-----------------------------------------------------------
SUBROUTINE READ_FOURIER2D(this,xy,f,klist)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier2d),INTENT(INOUT) :: this
   REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: xy
   REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: f
   INTEGER(KIND=4),DIMENSION(:,:),INTENT(IN) :: klist
   ! Local variables
   INTEGER(KIND=4) :: ndata,nfunc
   ! Run section
   ndata=size(f(1,:))
   nfunc=size(f(:,1))
   SELECT CASE(size(xy(:,1)) == ndata .AND. size(xy(1,:))==2)
      CASE(.FALSE.)
         WRITE(0,*) "READ_FOURIER2D ERR: dimensions mismatch in arrays xy or f"
         CALL EXIT(1)
      CASE(.TRUE.)
         ! do nothing
   END SELECT
   this%n=ndata
   ALLOCATE(this%xy(ndata,2))
   ALLOCATE(this%f(nfunc,ndata))
   ALLOCATE(this%klist(ndata,2))
   this%f = f
   this%nfunc = nfunc
   this%xy = xy
   this%klist=klist
   RETURN
END SUBROUTINE READ_FOURIER2D  
END MODULE FOURIER2D_MOD
