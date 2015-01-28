!##################################################################################
! MODULE: FOURIER2D_MOD
!> @brief
!! Provides tools to perform 2D interpolations with fourier series
!##################################################################################
MODULE FOURIER2D_MOD
use SURFACE_MOD, only: Surface
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
!> @date Feb/2014, Jul/2014
!> @version 2.0 
!---------------------------------------------------------------
TYPE,ABSTRACT :: Fourier2d
   INTEGER(KIND=4) :: n
   INTEGER(KIND=4) :: nfunc
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: xy
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: f
   INTEGER(KIND=4),DIMENSION(:,:),ALLOCATABLE :: klist
CONTAINS
   PROCEDURE,PUBLIC,NON_OVERRIDABLE :: READ => READ_FOURIER2D
   PROCEDURE(INTERPOL_FOURIER2D),PUBLIC,DEFERRED :: INTERPOL  ! DEFERRED !!!! take a look to interface
   PROCEDURE(GET_F_AND_DERIVS_FOURIER2D),PUBLIC,DEFERRED :: GET_F_AND_DERIVS ! DEFERRED !!!! take a look to interface
END TYPE Fourier2d

ABSTRACT INTERFACE
   !###########################################################
   !# SUBROUTINE: INTERPOL_FOURIER2D 
   !###########################################################
   !-----------------------------------------------------------
   SUBROUTINE INTERPOL_FOURIER2D(this,surf,filename)
      IMPORT Fourier2d
      IMPORT Surface
      CLASS(Fourier2d),INTENT(INOUT) :: this
      TYPE(Surface),INTENT(IN) :: surf
      CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
   END SUBROUTINE INTERPOL_FOURIER2D
   !###########################################################
   !# SUBROUTINE: GET_F_AND_DERIVS 
   !###########################################################
   !-----------------------------------------------------------
   SUBROUTINE GET_F_AND_DERIVS_FOURIER2D(this,surf,r,v,dvdu)
      IMPORT Fourier2d
      IMPORT Surface
      CLASS(Fourier2d),INTENT(IN):: this
      TYPE(Surface),INTENT(IN):: surf
      REAL(KIND=8),DIMENSION(2),INTENT(IN) :: r
      REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: v
      REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: dvdu
   END SUBROUTINE GET_F_AND_DERIVS_FOURIER2D
END INTERFACE
!////////////////////////////////////////////////////////////////
CONTAINS
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
