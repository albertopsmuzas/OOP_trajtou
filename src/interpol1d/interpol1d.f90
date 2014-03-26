!################################################################################################### 
! MODULE: INTERPOL1D
!
!> @brief  
!! Module that manages different interpolation schemes for one variable
!
!> @details 
!! All types and subprograms intended to create interpolations in 1D should be placed inside this module
!##################################################################################################
MODULE INTERPOL1D_MOD
! Initial declarations
IMPLICIT NONE
!//////////////////////////////////////////////////////////////////////
! TYPE: interpol1d
!> @brief
!! Generic type of one dimensional interpolations
!
!> @details
!! All kinds of interpolations should be declared as extensions of this type
!
!> @param n - Number of data
!> @param x - A set of @f$x_{i}@f$ data
!> @param f - A set of @f$f(x_{i})@f$ data
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 28/Jan/2014
!> @version 1.0
!-----------------------------------------------------------------------
TYPE :: Interpol1d
   INTEGER(KIND=4) :: n
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: x
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: f
   CONTAINS
      PROCEDURE,PUBLIC :: READ => READ_INTERPOL1D
      PROCEDURE,PUBLIC :: PLOTDATA => PLOT_DATA_INTERPOL1D
      PROCEDURE,PUBLIC :: PLOT_INTERPOL => PLOT_INTERPOL_INTERPOL1D
      PROCEDURE,PUBLIC :: getvalue => getvalue_interpol1d ! child types, override this
      PROCEDURE,PUBLIC :: getderiv => getderiv_interpol1d ! child types, override this
END TYPE Interpol1d
!//////////////////////////////////////////////////////////////////////
CONTAINS
!######################################################################
! SUBROUTINE: READ_INTERPOL1D #########################################
!######################################################################
!> @brief
!! Read main parameters for a 1D interpolation from arguments
!> @see interpol1d
!----------------------------------------------------------------------
SUBROUTINE READ_INTERPOL1D(this,x,f)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Interpol1d),INTENT(INOUT) :: this
   REAL(KIND=8), DIMENSION(:),INTENT(IN) :: x,f
   ! Run section -------
   IF (size(x)/=size(f)) THEN
      WRITE(0,*) "READ_INTERPOL1D: dimensions mismatch between x and f"
      CALL EXIT(1)
   END IF
   this%n=size(x)
   ALLOCATE(this%x(this%n))
   ALLOCATE(this%f(this%n))
   this%x = x
   this%f = f
   RETURN
END SUBROUTINE READ_INTERPOL1D
!######################################################################
! SUBROUTINE: READ_INTERPOL1D #########################################
!######################################################################
!> @brief
!! Plots in a file the data stored in an inteprol1d type variable
!> @see interpol1d
!----------------------------------------------------------------------
SUBROUTINE PLOT_DATA_INTERPOL1D(this,filename)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Interpol1d),INTENT(IN) :: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables
   INTEGER(KIND=4) :: i ! counter
   ! Run section ---------------------
   OPEN (UNIT=11,FILE=filename,STATUS="replace",ACTION="write")
   DO i = 1, this%n
      WRITE(11,*) this%x(i),this%f(i)
   END DO
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_DATA_INTERPOL1D
!######################################################################
! SUBROUTINE: PLOT_INTERPOL_INTERPOL1D ################################
!######################################################################
!> @brief
!! Creates a data file called @b filename with the interpolation graphic of
!! this interpol1d type variable. The number of points in that graphic is defined
!! by @b npoints. Cannot be less than two. It also plots the first derivative
!----------------------------------------------------------------------
SUBROUTINE PLOT_INTERPOL_INTERPOL1D(this,npoints,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   INTEGER,INTENT(IN) :: npoints
   CLASS(Interpol1d),INTENT(IN),TARGET :: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL*8 :: delta, interval, x
   INTEGER :: i ! Counter
   CHARACTER(LEN=26), PARAMETER :: routinename = "PLOT_INTERPOL_INTERPOL1D: " 
   ! Pointers ------------------------------------
   REAL*8, POINTER :: xmin, xmax
   INTEGER, POINTER :: n
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT_INTERPOL_INTERPOL1D ERR: Less than 2 points"
      STOP
   END IF
   !
   n => this%n
   xmin => this%x(1)
   xmax => this%x(n)
   !
   interval=xmax-xmin
   inpoints=npoints-2
   ndelta=npoints-1
   delta=interval/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   WRITE(11,*) xmin,this%getvalue(xmin),this%getderiv(xmin)
   DO i=1, inpoints
      x=xmin+(dfloat(i)*delta)
      WRITE(11,*) x,this%getvalue(x),this%getderiv(x)
   END DO
   WRITE(11,*) xmax,this%getvalue(xmax),this%getderiv(xmax)
   WRITE(*,*) routinename,filename," file created"
   CLOSE(11)
END SUBROUTINE PLOT_INTERPOL_INTERPOL1D
!###########################################################
!# FUNCTION: getvalue_interpol1d 
!###########################################################
!> @brief
!! Dummy function. Override it!!
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getvalue_interpol1d(this,x,shift) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Interpol1d),TARGET,INTENT(IN) :: this
   REAL(KIND=8),INTENT(IN) :: x
   REAL(KIND=8),OPTIONAL,INTENT(IN) :: shift
   ! Local variables
   REAL(KIND=8) :: dummy
   ! Run section
   dummy=shift+x+this%x(1) ! just to avoid warnings during compilation
   WRITE(0,*) "GET_VALUE_INTERPOL1D: If this routine was invoked, something is wrong"
   CALL EXIT(1)
   getvalue_interpol1d=dummy
   RETURN
END FUNCTION getvalue_interpol1d
!###########################################################
!# FUNCTION: getderiv_interpol1d 
!###########################################################
!> @brief
!! Dummy function. Override it!!
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_interpol1d(this,x,shift) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Interpol1d),TARGET,INTENT(IN) :: this
   REAL(KIND=8),INTENT(IN) :: x
   REAL(KIND=8),OPTIONAL,INTENT(IN) :: shift
   ! Local variables
   REAL(KIND=8) :: dummy
   ! Run section
   dummy=shift+x+this%x(1) ! Just to avoid warnings during compilation
   WRITE(0,*) "GET_DERIV_INTERPOL1D: If this routine was invoked, something is wrong"
   CALL EXIT(1)
   getderiv_interpol1d=dummy
   RETURN
END FUNCTION getderiv_interpol1d
END MODULE INTERPOL1D_MOD
