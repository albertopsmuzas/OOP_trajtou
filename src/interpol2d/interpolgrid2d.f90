!##################################################################################
! MODULE: INTERPOLGRID2D_MOD
!> @brief
!! Provides tools to perform 2D interpolations based on a reactangular/square grid
!##################################################################################
MODULE INTERPOLGRID2D_MOD
IMPLICIT NONE
!////////////////////////////////////////////////////////////////
! TYPE: Interpolgrid2d
!
!> @brief
!! Generic 2D interpolation type variable
!
!> @param n - Number of data points
!> @param x(:) - Grid in X. Only if input has grid structure
!> @param y(:) - Grid in Y. Only if input has grid structure
!> @param fgrid(:,:) - Function evaluated in a grid
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Mar/2014
!> @version 1.0 
!---------------------------------------------------------------
TYPE,ABSTRACT :: Interpolgrid2d
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: x,y
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: fgrid
CONTAINS
   PROCEDURE,PUBLIC :: READ => READ_INTERPOLGRID2D
   PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_INTERPOLGRID2D
   PROCEDURE,PUBLIC :: PLOTDATA => PLOTDATA_INTERPOLGRID2D
END TYPE Interpolgrid2d
!////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: INTERPOL_INTERPOLGRID2D 
!###########################################################
!> @brief
!! Dummy subroutine to be overriden
!-----------------------------------------------------------
SUBROUTINE INTERPOL_INTERPOLGRID2D(this,filename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Interpolgrid2d),INTENT(INOUT) :: this
   CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
   ! Run section
   WRITE(0,*) "INTERPOL_INTERPOLGRID2D ERR: You didn't allocate this variable with the proper type"
   CALL EXIT(1)
   RETURN
END SUBROUTINE INTERPOL_INTERPOLGRID2D
!##########################################################
! SUBROUTINE: READGRID_INTERPOL2D
!> @brief
!! Reads input defined in a grid
!
!> @param[out] this - interpol 2D object to be set up
!> @param[in] x(:) - X grid
!> @param[in] y(:) - Y grid
!> @param[in] f(:,:) - stores F values for each point in the grid
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 17/Feb/2014
!> @version 1.0 
!----------------------------------------------------------
SUBROUTINE READ_INTERPOLGRID2D(this,x,y,f)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Interpolgrid2d),INTENT(OUT) :: this
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
END SUBROUTINE READ_INTERPOLGRID2D
!###########################################################
!# SUBROUTINE: PLOTDATA_INTERPOLGRID2D 
!###########################################################
!> @brief
!! Plot data stored in the grid, not interpolated values
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE PLOTDATA_INTERPOLGRID2D(this,filename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Interpolgrid2d),INTENT(IN) :: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables
   INTEGER(KIND=4) :: nx,ny
   INTEGER(KIND=4) :: i,j ! counters
   ! Parameter section
   character(len=*),parameter:: routineName='PLOTDATA_INTERPOLGRID2D: '
   ! Run section
   nx=size(this%x)
   ny=size(this%y)
   OPEN (10,FILE=filename,STATUS="replace",ACTION="write")
   DO i = 1, nx
      DO j = 1, ny
         WRITE(10,*) this%x(i),this%y(j),this%fgrid(i,j)
      END DO
   END DO
   write(*,*) routineName//'Graph created: '//fileName
   RETURN
END SUBROUTINE PLOTDATA_INTERPOLGRID2D
END MODULE INTERPOLGRID2D_MOD
