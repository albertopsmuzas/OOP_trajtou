!########################################################
! MODULE : FOURIER1D
!
!> @brief
!! Provides tools to perform 1D periodical interpolations
!
!> @warning
!! - Includes Interpol1d_mod in its scope
!########################################################
MODULE FOURIER1D_MOD
USE INTERPOL1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////
! TYPE: FOURIER1D
!> @brief
!! Class to store all information needed for a 1D REAL fourier inteprolation
!------------------------------------------------
TYPE,EXTENDS(Interpol1d):: Fourier1d
   PRIVATE
   LOGICAL :: par
   INTEGER(KIND=4) :: order
   REAL(KIND=8) :: period
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: coeffpar,coeffodd
   CONTAINS
      PROCEDURE,PUBLIC :: ispar => ispar_fourier1d
      PROCEDURE,PUBLIC :: getperiod => getperiod_fourier1d
      PROCEDURE,PUBLIC :: getorder => getorder_fourier1d
      PROCEDURE,PUBLIC :: SET_COEFF => SET_COEFF_FOURIER1D
      PROCEDURE,PUBLIC :: getvalue => getvalue_fourier1d
      PROCEDURE,PUBLIC :: getderiv => getderiv_fourier1d
END TYPE Fourier1d
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# FUNCTION: getorder_fourier1d 
!###########################################################
!> @brief
!! Simple get function
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 03/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
INTEGER(KIND=4) FUNCTION getorder_fourier1d(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d),INTENT(IN) :: this
   ! Run section
   getorder_fourier1d=this%order
   RETURN
END FUNCTION getorder_fourier1d
!###########################################################
!# FUNCTION: getperiod_fourier1d 
!###########################################################
!> @brief 
!! Simple get function
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 03/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getperiod_fourier1d(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d),INTENT(IN) :: this
   ! Run section ---------------
   getperiod_fourier1d=this%period
   RETURN
END FUNCTION getperiod_fourier1d  
!##################################################
! FUNCTION: ispar_fourier1d
!> @brief
!! Basic enquire function. If the function is par,
!! sinus elements vanish from expansion
!--------------------------------------------------
LOGICAL FUNCTION ispar_fourier1d(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d),INTENT(IN) :: this
   ! Run section ----------------
   ispar_fourier1d = this%par
   RETURN
END FUNCTION ispar_fourier1d
!#################################################
! SUBROUTINE: SET_COEFF_FURIER1D
!> @brief
!! Performs interpolation
!-------------------------------------------------
SUBROUTINE SET_COEFF_FOURIER1D(this,period,ispar)
   ! Initial declarations
   USE MATHS_MOD
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d),INTENT(INOUT) :: this
   REAL(KIND=8),INTENT(IN) :: period
   LOGICAL,INTENT(IN),OPTIONAL :: ispar
   ! Local variables
   INTEGER(KIND=4) :: q ! max order of the expansion
   INTEGER(KIND=4) :: i,l ! counters
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: m ,inv_m
   REAL(KIND=8), DIMENSION(:),ALLOCATABLE :: aux
   ! Run section -----------------
   q=((this%n-1)/2) ! if this%n-1 is not a par number, this'd be an error
   this%order=q
   ALLOCATE(m(this%n,this%n))
   ALLOCATE(inv_m(this%n,this%n))
   SELECT CASE(ispar)
      CASE(.TRUE.) ! is par, only par coefficients
         this%par=.TRUE.
         ALLOCATE(this%coeffpar(this%n))
         FORALL(i=1:this%n,l=1:this%n)
            m(i,l)=dcos((2.D0*pi/period)*dfloat(l-1)*this%f(i))
         END FORALL
         CALL INV_MATRIX(this%n,m,inv_m)
         this%coeffpar=matmul(inv_m,this%f)
         RETURN
      CASE(.FALSE.)
         this%par=.FALSE. 
         ALLOCATE(this%coeffpar(q+1))
         ALLOCATE(this%coeffodd(q))
         ALLOCATE(aux(this%n))
         FORALL(i=1:this%n,l=1:q+1)
            m(i,l)=dcos((2.D0*pi/period)*dfloat(l-1)*this%f(i))
         END FORALL
         FORALL(i=1:this%n,l=q+2:this%n)
            m(i,l)=dsin((2.D0*pi/period)*dfloat(l-1)*this%f(i))
         END FORALL
         CALL INV_MATRIX(this%n,m,inv_m)
         aux=matmul(inv_m,this%f)
         this%coeffpar=aux(1:q+1)
         this%coeffodd=aux(q+2:this%n)
         RETURN
      CASE DEFAULT
         WRITE(0,*) "SET_COEFF_FOURIER1D ERR: Wrong value for ispar variable"
         CALL EXIT(1)
   END SELECT
END SUBROUTINE SET_COEFF_FOURIER1D
!###########################################################
!# FUNCTION: getvalue_fourier1d 
!###########################################################
!> @brief 
!! Get's fourier 1D interpolation value for a given point
!! inside the correct range
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 03/03/2014
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getvalue_fourier1d(this,x,shift)
   ! Initial declarations   
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d),TARGET,INTENT(IN) :: this
   REAL(KIND=8),INTENT(IN) :: x
   REAL(KIND=8),INTENT(IN),OPTIONAL :: shift
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: m
   ! Run section
   IF (x <= this%f(1) .OR. x >= this%f(this%n)) THEN
      WRITE(0,*) "getvalue_fourier1d ERR: requested f(x) value is out of range"
      CALL EXIT(1)
   END IF
   ALLOCATE(m(this%n))
   SELECT CASE(this%par)
      CASE(.TRUE.)
         FORALL(i=1:this%n) m(i)=dcos((2.D0*pi/this%period)*dfloat(i-1)*x)
         getvalue_fourier1d=dot_product(this%coeffpar,m)
      CASE(.FALSE.)
         FORALL(i=1:this%order+1) m(i)=dcos((2.D0*pi/this%period)*dfloat(i-1)*x)
         FORALL(i=this%order+2:this%n) m(i)=dsin((2.D0*pi/this%period)*dfloat(i-this%order-1)*x)
         getvalue_fourier1d=dot_product(m(1:(this%order+1)),this%coeffpar)
         getvalue_fourier1d=getvalue_fourier1d+dot_product(m(this%order+2:this%n),this%coeffodd)
      CASE DEFAULT
         ! do nothing
   END SELECT
   RETURN
END FUNCTION getvalue_fourier1d
!###########################################################
!# FUNCTION: getderiv_fourier1d 
!###########################################################
!> @brief 
!! Get's derivative value at a given point X using the interpolation
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 03/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_fourier1d(this,x,shift) 
   ! Initial declarations   
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d),TARGET,INTENT(IN) :: this
   REAL(KIND=8),INTENT(IN) :: x
   REAL(KIND=8),INTENT(IN),OPTIONAL :: shift
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: m
   ! Run section
   IF (x <= this%f(1) .OR. x >= this%f(this%n)) THEN
      WRITE(0,*) "getvalue_fourier1d ERR: requested f(x) value is out of range"
      CALL EXIT(1)
   END IF
   ALLOCATE(m(this%n))
   SELECT CASE(this%par)
      CASE(.TRUE.)
         FORALL(i=1:this%n) m(i)=-dsin((2.D0*pi/this%period)*dfloat(i-1)*x)*(2.D0*pi/this%period)*dfloat(i-1)
         getderiv_fourier1d=dot_product(this%coeffpar,m)
      CASE(.FALSE.)
         FORALL(i=1:this%order+1)
            m(i)=-dsin((2.D0*pi/this%period)*dfloat(i-1)*x)*(2.D0*pi/this%period)*dfloat(i-1)
         END FORALL
         FORALL(i=this%order+2:this%n) 
            m(i)=dcos((2.D0*pi/this%period)*dfloat(i-this%order-1)*x)*(2.D0*pi/this%period)*dfloat(i-this%order-1)
         END FORALL
         getderiv_fourier1d=dot_product(m(1:this%order+1),this%coeffpar)
         getderiv_fourier1d=getderiv_fourier1d+dot_product(m(this%order+2:this%n),this%coeffodd)
      CASE DEFAULT
         ! do nothing
   END SELECT
   RETURN
END FUNCTION getderiv_fourier1d
END MODULE FOURIER1D_MOD
