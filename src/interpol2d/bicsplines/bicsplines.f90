!##########################################################
! MODULE: BICUBICSPLINES
!
!> @brief 
!! Provides tools to perform bicubic splines interpolations on
!! 2D functions
!##########################################################
MODULE BICSPLINES_MOD
USE INTERPOL2D_MOD
USE CUBICSPLINES_MOD
IMPLICIT NONE
!//////////////////////////////////////////////////////////
! TYPE: BICSPLINES
!> @brief
!! Class to store data for a bicubic splines interpolation 
!
!> @param xstring - Interpolation in 1D for each string of control
!!                  points in x
!> @param ystring - Interpolation in 1D for each string of control
!!                  points in y
!----------------------------------------------------------
TYPE,EXTENDS(Interpol2d) :: Bicsplines
   PRIVATE
   TYPE(Csplines),DIMENSION(:),ALLOCATABLE :: xcsplines
   TYPE(Csplines),DIMENSION(:),ALLOCATABLE :: ycsplines
   REAL(KIND=8),DIMENSION(:,:,:,:),ALLOCATABLE :: coeff
   CONTAINS
      PROCEDURE,PUBLIC :: SET_COEFF => SET_COEFF_BICSPLINES
      PROCEDURE,PUBLIC :: getvalue => getvalue_bicsplines
      PROCEDURE,PUBLIC :: getderivx => getderivx_bicsplines
      PROCEDURE,PUBLIC :: getderivy => getderivy_bicsplines
      PROCEDURE,PUBLIC :: PLOT_XYMAP => PLOT_XYMAP_BICSPLINES
      PROCEDURE,PUBLIC :: PLOT_SPLINES => PLOT_SPLINES_BICSPLINES
END TYPE
!//////////////////////////////////////////////////////////
CONTAINS
!##########################################################
! SUBROUTINE: SEC_COEFF_BICSPLINES
!> @brief
!! Sets coefficients
!
!> @param[in,out] this - bicubic splines object to set coefficients
!----------------------------------------------------------
SUBROUTINE SET_COEFF_BICSPLINES(this,filename)
   ! Initial declarations
   USE MATHS_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(Bicsplines),INTENT(INOUT) :: this
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: filename
   ! Local variables
   INTEGER(KIND=4) :: nx,ny
   INTEGER(KIND=4) :: i,j,k ! counters
   REAL(KIND=8) :: x0,x1,y0,y1
   REAL(KIND=8),DIMENSION(4,4) :: xmtrx,inv_xmtrx,ymtrx,inv_ymtrx,smtrx
   ! Run section -------------------------------
   nx=size(this%x)
   ny=size(this%y)
   ALLOCATE(this%xcsplines(nx))
   ALLOCATE(this%ycsplines(ny))
   ALLOCATE(this%coeff(nx,ny,4,4))
   DO i = 1, ny
      CALL this%xcsplines(i)%READ(this%x,this%fgrid(:,i))
      CALL this%xcsplines(i)%INTERPOL(0.D0,0,0.D0,0) ! last and initial 2 splines are equal
   END DO
   DO i = 1, nx
      CALL this%ycsplines(i)%READ(this%y,this%fgrid(i,:))
      CALL this%ycsplines(i)%INTERPOL(0.D0,0,0.D0,0) ! last and initial 2 splines are equal
   END DO
   !
   IF(present(filename)) OPEN (521,FILE=filename,STATUS="replace",ACTION="write")
   DO i = 1, nx-1
      DO j = 1, ny-1
         x0=this%xcsplines(j)%x(i)
         x1=this%xcsplines(j)%x(i+1)

         y0=this%ycsplines(i)%x(j)
         y1=this%ycsplines(i)%x(j+1)

         xmtrx(1,:)=(/1.D0,x0,x0**2.D0,x0**3.D0/)
         xmtrx(2,:)=(/1.D0,x1,x1**2.D0,x1**3.D0/)
         xmtrx(3,:)=(/0.D0,1.D0,2.D0*x0,3.D0*x0**2.D0/)
         xmtrx(4,:)=(/0.D0,1.D0,2.D0*x1,3.D0*x1**2.D0/)
         ymtrx(:,1)=(/1.D0,y0,y0**2.D0,y0**3.D0/)
         ymtrx(:,2)=(/1.D0,y1,y1**2.D0,y1**3.D0/)
         ymtrx(:,3)=(/0.D0,1.D0,2.D0*y0,3.D0*y0**2.D0/)
         ymtrx(:,4)=(/0.D0,1.D0,2.D0*y1,3.D0*y1**2.D0/)

         smtrx(1,1)=this%fgrid(i,i)
         smtrx(1,2)=this%fgrid(i,j+1)
         smtrx(2,1)=this%fgrid(i+1,j)
         smtrx(2,2)=this%fgrid(i+1,j+1)

         smtrx(1,3)=this%ycsplines(i)%getderiv(y0)
         smtrx(1,4)=this%ycsplines(i)%getderiv(y1)
         smtrx(2,3)=this%ycsplines(i+1)%getderiv(y0)
         smtrx(2,4)=this%ycsplines(i+1)%getderiv(y1)

         smtrx(3,1)=this%xcsplines(j)%getderiv(x0)
         smtrx(3,2)=this%xcsplines(j+1)%getderiv(x0)
         smtrx(4,1)=this%xcsplines(j)%getderiv(x1)
         smtrx(4,4)=this%xcsplines(j+1)%getderiv(x1)

         smtrx(3,3)=d2fdxdy_finitdiff(i,j,this%x,this%y,this%fgrid)
         smtrx(3,4)=d2fdxdy_finitdiff(i,j+1,this%x,this%y,this%fgrid)
         smtrx(4,3)=d2fdxdy_finitdiff(i+1,j,this%x,this%y,this%fgrid)
         smtrx(4,4)=d2fdxdy_finitdiff(i+1,j+1,this%x,this%y,this%fgrid)
         
         CALL INV_MTRX(4,xmtrx,inv_xmtrx)
         CALL INV_MTRX(4,ymtrx,inv_ymtrx)
         this%coeff(i,j,:,:)=matmul(matmul(inv_xmtrx,smtrx),inv_ymtrx)
         IF (present(filename)) THEN
            WRITE(521,*) "BICUBIC SPLINE :",i,j
            WRITE(521,*) "=========================="
            WRITE(521,*) "X matrix:"
            DO k = 1, 4
               WRITE(521,*) xmtrx(k,:)
            END DO
            WRITE(521,*) "inv_X matrix:"
            DO k = 1, 4
               WRITE(521,*) inv_xmtrx(k,:)
            END DO
            WRITE(521,*) "Y matrix:"
            DO k = 1, 4
               WRITE(521,*) ymtrx(k,:)
            END DO
            WRITE(521,*) "inv_Y matrix:"
            DO k = 1, 4
               WRITE(521,*) inv_ymtrx(k,:)
            END DO
            WRITE(521,*) "S matrix:"
            DO k = 1, 4
               WRITE(521,*) smtrx(k,:)
            END DO
            WRITE(521,*) "Coeff. matrix:"
            DO k = 1, 4
               WRITE(521,*) this%coeff(i,j,k,:)
            END DO
         END IF
      END DO
   END DO
   IF(present(filename)) CLOSE(521)
   RETURN
END SUBROUTINE SET_COEFF_BICSPLINES
!###########################################################
!# FUNCTION: d2fdxdy_finitdiff
!###########################################################
!> @brief
!! Calculates cross-term second derivatives at grid points.
!! Uses finite differences. The grid can have different steps.
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION d2fdxdy_finitdiff(i,j,x,y,f) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: i,j
   REAL(KIND=8),DIMENSION(:) :: x,y
   REAL(KIND=8),DIMENSION(:,:) :: f
   ! Local variables
   REAL(KIND=8) :: hx0,hx1
   REAL(KIND=8) :: hy0,hy1
   INTEGER(KIND=4) :: nx,ny
   ! Run section
   nx=size(x)
   ny=size(y)
   IF (nx/=size(f(:,1))) THEN
      WRITE(0,*) "d2fdxdy: size mismatch of arrays x and f"
      CALL EXIT(1)
   ELSE IF (ny/=size(f(1,:))) THEN
      WRITE(0,*) "d2fdxdy: size mismatch of arrays y and f"
      CALL EXIT(1)
   END IF
   ! initialize variables
   hx0=0.D0
   hx1=0.D0
   hy0=0.D0
   hy1=0.D0
   !
   ! Check if we are in a corner, edge or bulk point of the grid 
   IF ( i/=1 .AND. i/=nx .AND. j/=1 .AND. j/=ny) THEN ! we are in the 2D bulk
      hx0=x(i+1)-x(i)
      hx1=x(i)-x(i-1)
      hy0=y(j+1)-y(j)
      hy1=y(j)-y(j-1)
      d2fdxdy_finitdiff=((1.D0/hx1)-(1.D0/hx0))*((1.D0/hy1)-(1.D0/hy0))*f(j,i)
      d2fdxdy_finitdiff=d2fdxdy_finitdiff+((1.D0/hx1)-(1.D0/hx0))*((f(i,j+1)/hy0)-(f(i,j-1)/hy1))
      d2fdxdy_finitdiff=d2fdxdy_finitdiff+((1.D0/hy1)-(1.D0/hy0))*((f(i+1,j)/hx0)-(f(i-1,j)/hx1))
      d2fdxdy_finitdiff=d2fdxdy_finitdiff+(f(i+1,j+1)/(hx0*hy0))-(f(i+1,j-1)/(hx0*hy1))-(f(i-1,j+1)/(hx1*hy0))
      d2fdxdy_finitdiff=d2fdxdy_finitdiff+(f(i-1,j-1)/(hx1*hy1))
      d2fdxdy_finitdiff=0.25D0*d2fdxdy_finitdiff
      RETURN
   ELSE IF (i==1 .AND. j==1) THEN ! corner ++
      hx0=x(i+1)-x(i)
      hy0=y(j+1)-y(j)
      d2fdxdy_finitdiff=(1.D0/(hx0*hy0))*(f(i,j)-f(i,j+1)-f(i+1,j)+f(i+1,j+1))
      RETURN
   ELSE IF  (i==1 .AND. j==ny) THEN ! corner +-
      hx0=x(i+1)-x(i)
      hy1=y(j)-y(j-1)
      d2fdxdy_finitdiff=(1.D0/(hx0*hy1))*(f(i+1,j)-f(i+1,j-1)+f(i,j-1)-f(i,j))
     RETURN
   ELSE IF (i==nx .AND. j==1) THEN ! corner -+
      hx1=x(i)-x(i-1)
      hy0=y(j+1)-y(j)
      d2fdxdy_finitdiff=(1.D0/(hx1*hy0))*(f(i,j+1)-f(i-1,j+1)+f(i-1,j)-f(i,j))
     RETURN
   ELSE IF (i==nx .AND. j==ny) THEN ! corner --
      hx1=x(i)-x(i-1)
      hy1=y(j)-y(j-1)
      d2fdxdy_finitdiff=(1.D0/(hx1*hy1))*(f(i,j)-f(i,j-1)-f(i-1,j)+f(i-1,j-1))
      RETURN
   ELSE IF (i==1) THEN ! left edge
      hx0=x(i+1)-x(i)
      hy0=y(j+1)-y(j)
      hy1=y(j)-y(j-1)
      d2fdxdy_finitdiff=((1.D0/hy1)-(1.D0/hy0))*f(i+1,j)+(1.D0/hy0)*f(i+1,j+1)-(1.D0/hy1)*f(i+1,j-1)
      d2fdxdy_finitdiff=d2fdxdy_finitdiff-((1.D0/hy1)-(1.D0/hy0))*f(i,j)-(1.D0/hy0)*f(i,j+1)+(1.D0/hy1)*f(i,j-1)
      d2fdxdy_finitdiff=d2fdxdy_finitdiff*(1.D0/(2.D0*hx0))
      RETURN
   ELSE IF (i==nx) THEN ! right edge
      hx1=x(i)-x(i-1)
      hy0=y(j+1)-y(j)
      hy1=y(j)-y(j-1)
      d2fdxdy_finitdiff=((1.D0/hy1)-(1.D0/hy0))*f(i,j)+(1.D0/hy0)*f(i,j+1)-(1.D0/hy1)*f(i,j-1)
      d2fdxdy_finitdiff=d2fdxdy_finitdiff-((1.D0/hy1)-(1.D0/hy0))*f(i-1,j)-(1.D0/hy0)*f(i-1,j+1)+(1.D0/hy1)*f(i-1,j-1)
      d2fdxdy_finitdiff=d2fdxdy_finitdiff*(1.D0/(2.D0*hx1))
      RETURN
   ELSE IF (j==1) THEN ! down edge
      hx0=x(i+1)-x(i)
      hx1=x(i)-x(i-1)
      hy0=y(j+1)-y(j)
      d2fdxdy_finitdiff=((1.D0/hx1)-(1.D0/hx0))*f(i,j+1)+(1.D0/hx0)*f(i+1,j+1)-(1.D0/hx1)*f(i-1,j+1)
      d2fdxdy_finitdiff=d2fdxdy_finitdiff-((1.D0/hx1)-(1.D0/hx0))*f(i,j)-(1.D0/hx0)*f(i+1,j)+(1.D0/hx1)*f(i-1,j)
      d2fdxdy_finitdiff=d2fdxdy_finitdiff*(1.D0/(2.D0*hy0))
      RETURN
   ELSE IF (j==ny) THEN ! upper edge
      hx0=y(i+1)-y(i)
      hx1=y(i)-y(i-1)
      hy1=x(j)-x(j-1)
      d2fdxdy_finitdiff=((1.D0/hx1)-(1.D0/hx0))*f(i,j)+(1.D0/hx0)*f(i+1,j)-(1.D0/hx1)*f(i-1,j)
      d2fdxdy_finitdiff=d2fdxdy_finitdiff-((1.D0/hx1)-(1.D0/hx0))*f(i,j-1)-(1.D0/hx0)*f(i+1,j-1)+(1.D0/hx1)*f(i-1,j-1)
      d2fdxdy_finitdiff=d2fdxdy_finitdiff*(1.D0/(2.D0*hy1))
      RETURN
   END IF      
   WRITE(0,*) "d2fdxdy_finitdiff ERR: If this message was printed, something is wrong"
   CALL EXIT(1)
   !
END FUNCTION d2fdxdy_finitdiff
!#############################################################
! SUBROUTINE: getvalue_bicsplines
!#############################################################
!> @brief
!! Gives the value of the interpolation for a 2D point
!-------------------------------------------------------------
REAL(KIND=8) FUNCTION getvalue_bicsplines(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Bicsplines),TARGET,INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),POINTER:: x1,x2
   REAL(KIND=8),DIMENSION(:,:),POINTER :: coeff
   REAL(KIND=8),DIMENSION(4) :: vecy,vecx
   INTEGER(KIND=4) :: nx,ny
   INTEGER(KIND=4) :: i,j
   ! Run section
   nx=size(this%x)
   ny=size(this%y)
   DO i = 1, nx-1
      x1 => this%x(i)
      x2 => this%x(i+1)
      IF((x(1).LE.x2).AND.(x(1).GE.x1)) EXIT
   END DO
   DO j = 1, ny-1
      x1 => this%y(j)
      x2 => this%y(j+1)
      IF((x(2).LE.x2).AND.(x(2).GE.x1)) EXIT
   END DO
   ! Now, i and j have the correct value
   coeff => this%coeff(i,j,:,:)
   vecx=(/1.D0,x(1),x(1)**2.D0,x(1)**3.D0/)
   vecy=(/1.D0,x(2),x(2)**2.D0,x(2)**3.D0/)
   getvalue_bicsplines=dot_product(vecx,matmul(coeff,vecy))
   RETURN
END FUNCTION getvalue_bicsplines
!#############################################################
! SUBROUTINE: getderivx_bicsplines
!#############################################################
!> @brief
!! Gives the x derivative of the interpolation for a 2D point
!-------------------------------------------------------------
REAL(KIND=8) FUNCTION getderivx_bicsplines(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Bicsplines),TARGET,INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),POINTER :: x1,x2
   REAL(KIND=8),DIMENSION(:,:),POINTER :: coeff
   REAL(KIND=8),DIMENSION(4) :: vecy,vecx
   INTEGER(KIND=4) :: nx,ny
   INTEGER(KIND=4) :: i,j
   ! Run section
   nx=size(this%x)
   ny=size(this%y)
   DO i = 1, nx-1
      x1 => this%x(i)
      x2 => this%x(i+1)
      IF((x(1).LE.x2).AND.(x(1).GE.x1)) EXIT
   END DO
   DO j = 1, ny-1
      x1 => this%y(j)
      x2 => this%y(j+1)
      IF((x(2).LE.x2).AND.(x(2).GE.x1)) EXIT
   END DO
   ! Now, i and j have the correct value
   coeff => this%coeff(i,j,:,:)
   vecx=(/0.D0,1.D0,2.D0*x(1),3.D0*x(1)**2.D0/)
   vecy=(/1.D0,x(2),x(2)**2.D0,x(2)**3.D0/)
   getderivx_bicsplines=dot_product(vecx,matmul(coeff,vecy))
   RETURN
END FUNCTION getderivx_bicsplines
!#############################################################
! SUBROUTINE: getderivy_bicsplines
!#############################################################
!> @brief
!! Gives the y derivative of the interpolation for a 2D point
!-------------------------------------------------------------
REAL(KIND=8) FUNCTION getderivy_bicsplines(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Bicsplines),TARGET,INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),POINTER :: x1,x2
   REAL(KIND=8),DIMENSION(:,:),POINTER :: coeff
   REAL(KIND=8),DIMENSION(4) :: vecy,vecx
   INTEGER(KIND=4) :: nx,ny
   INTEGER(KIND=4) :: i,j
   ! Run section
   nx=size(this%x)
   ny=size(this%y)
   DO i = 1, nx-1
      x1 => this%x(i)
      x2 => this%x(i+1)
      IF((x(1).LE.x2).AND.(x(1).GE.x1)) EXIT
   END DO
   DO j = 1, ny-1
      x1 => this%y(j)
      x2 => this%y(j+1)
      IF((x(2).LE.x2).AND.(x(2).GE.x1)) EXIT
   END DO
   ! Now, i and j have the correct value
   coeff => this%coeff(i,j,:,:)
   vecx=(/1.D0,x(1),x(1)**2.D0,x(1)**3.D0/)
   vecy=(/0.D0,1.D0,2.D0*x(2),3.D0*x(2)**2.D0/)
   getderivy_bicsplines=dot_product(vecx,matmul(coeff,vecy))
   RETURN
END FUNCTION getderivy_bicsplines
!#############################################################
! SUBROUTINE: getderivxy_bicsplines
!#############################################################
!> @brief
!! Gives the xy derivative of the interpolation for a 2D point
!-------------------------------------------------------------
REAL(KIND=8) FUNCTION getderivxy_bicsplines(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Bicsplines),TARGET,INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),POINTER :: x1,x2
   REAL(KIND=8),DIMENSION(:,:),POINTER :: coeff
   REAL(KIND=8),DIMENSION(4) :: vecy,vecx
   INTEGER(KIND=4) :: nx,ny
   INTEGER(KIND=4) :: i,j
   ! Run section
   nx=size(this%x)
   ny=size(this%y)
   DO i = 1, nx-1
      x1 => this%x(i)
      x2 => this%x(i+1)
      IF((x(1).LE.x2).AND.(x(1).GE.x1)) EXIT
   END DO
   DO j = 1, ny-1
      x1 => this%y(j)
      x2 => this%y(j+1)
      IF((x(2).LE.x2).AND.(x(2).GE.x1)) EXIT
   END DO
   ! Now, i and j have the correct value
   coeff => this%coeff(i,j,:,:)
   vecx=(/0.D0,1.D0,2.D0*x(1),3.D0*x(1)**2.D0/)
   vecy=(/0.D0,1.D0,2.D0*x(2),3.D0*x(2)**2.D0/)
   getderivxy_bicsplines=dot_product(vecx,matmul(coeff,vecy))
   RETURN
END FUNCTION getderivxy_bicsplines
!###############################################################
! SUBROUTINE: PLOT_XYMAP_BICSPLINES
!##############################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (X,Y)
!
!> @param[in] this - Interpolation 2D object
!> @param[in] filename - Name of the file to print the output
!> @param[in] init_xy - Initial position to start the scan 
!> @param[in] nxpoints - Number of points in X axis 
!> @param[in] nypoints - Number of points in Y axis
!> @param[in] Lx - Length of X axis
!> @param[in] Ly - Length of Y axis
!
!> @author A.S. Muzas
!> @date 17/Feb/2014
!> @version 1.0
!---------------------------------------------------------------
SUBROUTINE PLOT_XYMAP_BICSPLINES(this,filename,init_xy,nxpoints,nypoints,Lx,Ly)
   IMPLICIT NONE
   CLASS(Bicsplines),INTENT(IN) :: this
   REAL*8,DIMENSION(2),INTENT(IN) :: init_xy ! Initial position to start the scan (in a.u.)
   INTEGER,INTENT(IN) :: nxpoints, nypoints ! number of points in XY plane
   CHARACTER(LEN=*),INTENT(IN) :: filename ! filename
   REAL*8,INTENT(IN) :: Lx ! Length of X axis 
   REAL*8,INTENT(IN) :: Ly ! Length of X axis 
   ! Local variables
   REAL*8 :: xmin, ymin, xmax, ymax
   REAL*8, DIMENSION(2) :: r
   REAL*8 :: xdelta, ydelta
   INTEGER :: xinpoints, nxdelta
   INTEGER :: yinpoints, nydelta
   INTEGER :: i, j ! counters
   ! GABBA, GABBA HEY! ---------
   xmin = init_xy(1)
   ymin = init_xy(2)
   xmax = init_xy(1)+Lx
   ymax = init_xy(2)+Ly
   ! For X, grid parameters
   xinpoints=nxpoints-2
   nxdelta=nxpoints-1
   xdelta=Lx/DFLOAT(nxdelta)
   ! For Y, grid parameters
   yinpoints=nypoints-2
   nydelta=nypoints-1
   ydelta=(ymax-ymin)/DFLOAT(nydelta)
   ! Let's go! 
   ! 1st XY point
   OPEN(11,file=filename,status="replace")
   r(1) = xmin
   r(2) = ymin
   WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r)
   END DO
   r(2) = ymax
   WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r)
   ! inpoints in XY
   DO i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r)
      DO j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r)
      END DO
      r(2) = ymax
      WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r)
   END DO
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r)
   END DO
   r(2) = ymax
   WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r)
   CLOSE(11)
   WRITE(*,*) "PLOT_XYMAP_BICSPLINES: Graph created: ",filename
   RETURN
END SUBROUTINE PLOT_XYMAP_BICSPLINES
!###########################################################
!# SUBROUTINE: PLOT_SPLINES_BICSPLINES 
!###########################################################
!> @brief
!! Plot internal splines
!-----------------------------------------------------------
SUBROUTINE PLOT_SPLINES_BICSPLINES(this,npoints)
   ! Initial declarations   
      
   IMPLICIT NONE
   ! I/O variables
   CLASS(Bicsplines),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: npoints
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   INTEGER(KIND=4) :: nx,ny
   CHARACTER(LEN=100) :: filename
   ! Run section
   nx=size(this%xcsplines)
   ny=size(this%ycsplines)
   DO i = 1, ny
      WRITE(filename,'(I4,A12)') i,"-xspline.dat"
      filename=adjustl(filename)
      CALL this%xcsplines(i)%PLOT_INTERPOL(npoints,filename)
   END DO
   DO i = 1, nx
      WRITE(filename,'(I4,A12)') i,"-yspline.dat"
      filename=adjustl(filename)
      CALL this%ycsplines(i)%PLOT_INTERPOL(npoints,filename)
   END DO
   RETURN
END SUBROUTINE PLOT_SPLINES_BICSPLINES
END MODULE BICSPLINES_MOD
