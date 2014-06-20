!#########################################################
! MODULE: FUNCTION1D_MOD
!> @brief
!! Types and procedures to manage a function 
!##########################################################
MODULE FUNCTION1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: Function1d
!> @brief
!! General type for functions
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!----------------------------------------------------------------
TYPE,ABSTRACT :: Function1d
   PRIVATE
   CHARACTER(LEN=10) :: id
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: param
CONTAINS
   ! Initialization block
   PROCEDURE,PUBLIC,NON_OVERRIDABLE :: READ => READ_FUNCTION1D
   ! Set block
   PROCEDURE,PUBLIC,NON_OVERRIDABLE :: SET_ID => SET_ID_FUNCTION1D
   ! Get block
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: getparam => getparam_function1d
   PROCEDURE(getvalue_function1d),PUBLIC,DEFERRED :: getvalue
   PROCEDURE(getvalue_function1d),PUBLIC,DEFERRED :: getderiv
   ! Plot block
   PROCEDURE,PUBLIC,NON_OVERRIDABLE :: PLOT => PLOT_FUNCTION1D
END TYPE Function1d

ABSTRACT INTERFACE
   !###########################################################
   !# SUBROUTINE: getvalue_function1d 
   !###########################################################
   !-----------------------------------------------------------
   REAL(KIND=8) FUNCTION getvalue_function1d(this,x)
      IMPORT Function1d
      CLASS(Function1d),INTENT(IN)::this
      REAL(KIND=8),INTENT(IN) :: x
   END FUNCTION getvalue_function1d
END INTERFACE

!//////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: READ_FUNCTION1D
!###########################################################
!> @brief
!! Read parameters.
!-----------------------------------------------------------
SUBROUTINE READ_FUNCTION1D(this,param)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Function1d),INTENT(INOUT):: this
   REAL(KIND=8),DIMENSION(:) :: param
   ! Run section
   ALLOCATE(this%param(size(param)))
   this%param=param
   RETURN
END SUBROUTINE READ_FUNCTION1D
!###########################################################
!# FUNCTION: getparam_function1d 
!###########################################################
!> @brief 
!! Common get function. Gets atribute "param"
!-----------------------------------------------------------
FUNCTION getparam_function1d(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Function1d),INTENT(IN):: this
   ! Local variables
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: getparam_function1d
   ! Run section
   ALLOCATE(getparam_function1d(size(this%param)))
   getparam_function1d=this%param
   RETURN
END FUNCTION getparam_function1d
!###########################################################
!# SUBROUTINE: SET_ID_FUNCTION1D 
!###########################################################
!> @brief
!! Common set function. Sets atribute "id"
!-----------------------------------------------------------
SUBROUTINE SET_ID_FUNCTION1D(this,id)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Function1d),INTENT(INOUT):: this
   CHARACTER(LEN=*),INTENT(IN) :: id
   ! Run section
   this%id=id
   RETURN
END SUBROUTINE SET_ID_FUNCTION1D
!###########################################################
!# SUBROUTINE: PLOT_FUNCTION1D 
!###########################################################
! - Common plot function
!-----------------------------------------------------------
SUBROUTINE PLOT_FUNCTION1D(this,npoints,xmin,xmax,filename)
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(Function1d),INTENT(IN) :: this
   INTEGER,INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables -----------------------------
   INTEGER :: inpoints,ndelta,nfunc
   REAL(KIND=8) :: delta, interval, x
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: f,dfdx
   INTEGER :: i ! Counter
   CHARACTER(LEN=17), PARAMETER :: routinename = "PLOT_FUNCTION1D: " 
   ! Pointers ------------------------------------
   REAL(KIND=8):: xmin, xmax
   INTEGER(KIND=4):: n
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT_FUNCTION1D ERR: Less than 2 points"
      STOP
   END IF
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
   WRITE(11,*) xmax,this%getvalue(xmax),this%getderiv(x)
#ifdef DEBUG   
   CALL VERBOSE_WRITE(routinename,filename," file created")
#endif
   CLOSE(11)
END SUBROUTINE PLOT_FUNCTION1D
END MODULE FUNCTION1D_MOD 
