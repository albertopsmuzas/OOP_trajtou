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
   LOGICAL :: even=.FALSE.
   REAL(KIND=8) :: period
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE :: klist
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: coeff
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: extracoeff
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: extrafuncs
   CONTAINS
      ! get block
      PROCEDURE,PUBLIC :: getperiod => getperiod_fourier1d
      PROCEDURE,PUBLIC :: getvalue => getvalue_fourier1d
      PROCEDURE,PUBLIC :: getderiv => getderiv_fourier1d
      ! Initialization block
      PROCEDURE,PUBLIC :: READ_EXTRA => READ_EXTRA_DETAILS_FOURIER1D
      ! Enquire block
      PROCEDURE,PUBLIC :: is_even => is_even_fourier1d
      ! Tools
      PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_FOURIER1D
      PROCEDURE,PUBLIC :: ADD_MOREFUNCS => ADD_MORE_FUNCS_FOURIER1D
      PROCEDURE,PUBLIC :: GET_ALLFUNCS_AND_DERIVS => GET_ALLFUNC_AND_DERIVS_FOURIER1D
END TYPE Fourier1d
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: GET_ALLFUNC_AND_DERIVS_FOURIER1D
!###########################################################
!> @brief
!! Get value of the potential and derivs for an specific point x
!! for the main function and extra ones
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_ALLFUNC_AND_DERIVS_FOURIER1D(this,x,f,dfdx)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d),INTENT(IN) :: this
   REAL(KIND=8),INTENT(IN) :: x
   REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: f,dfdx
   ! Local variables
   INTEGER(KIND=4) :: nfuncs
   INTEGER(KIND=4) :: i ! counters
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: terms
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: terms_dx
   ! Check section
   SELECT CASE(allocated(this%extracoeff))
      CASE(.FALSE.)
         WRITE(0,*) "GET_ALLFUNCS_AND DERIVS ERR: extra coefficients are not allocated"
         WRITE(0,*) "GET_ALLFUNCS_AND DERIVS ERR: did you use ADD_MOREFUNCS and INTERPOL before this?"
         CALL EXIT(1)
      CASE(.TRUE.)
         ! do nothing
   END SELECT
   nfuncs=size(this%extrafuncs(:,1))+1
   SELECT CASE(size(f)/=nfuncs .OR. size(dfdx)/=nfuncs)
      CASE(.TRUE.)
         WRITE(0,*) "GET_ALLFUNCS_AND DERIVS ERR: size mismatch of output arguments"
         WRITE(0,*) "nfuncs: ",nfuncs
         WRITE(0,*) "size f: ", size(f)
         WRITE(0,*) "size dfdx: ",size(dfdx)
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   ! Run section
   ALLOCATE(terms(this%n))
   ALLOCATE(terms_dx(this%n))
   DO i = 1, this%n
      terms(i)=termfou1d(this%klist(i),this%period,x)
      terms_dx(i)=termfou1d_dx(this%klist(i),this%period,x)
   END DO
   f(1)=dot_product(terms,this%coeff)
   dfdx(1)=dot_product(terms_dx,this%coeff)
   DO i = 2, nfuncs
      f(i)=dot_product(terms,this%extracoeff(:,i-1))
      dfdx(i)=dot_product(terms_dx,this%extracoeff(:,i-1))
   END DO
   RETURN
END SUBROUTINE GET_ALLFUNC_AND_DERIVS_FOURIER1D
!###########################################################
!# SUBROUTINE: ADD_MORE_FUNCS_FOURIER1D
!###########################################################
!> @brief
!! brief description
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE ADD_MORE_FUNCS_FOURIER1D(this,f)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d),INTENT(INOUT) :: this
   REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: f
   ! Local variables
   INTEGER(KIND=4) :: nfuncs, ndata
   ! Run section
   nfuncs=size(f(:,1))
   ndata=size(f(1,:))
   SELECT CASE(ndata == this%n)
      CASE(.FALSE.)
         WRITE(0,*) "ADD_MORE_FUNCS_FOURIER1D ERR: size mismatch between extra functions and the original one"
         CALL EXIT(1)
      CASE(.TRUE.)
         ! donothing
   END SELECT
   ALLOCATE(this%extrafuncs(nfuncs,ndata))
   this%extrafuncs=f
   RETURN
END SUBROUTINE ADD_MORE_FUNCS_FOURIER1D
!###########################################################
!# SUBROUTINE: READ_EXTRA_DETAILS_FOURIER1D 
!###########################################################
!> @brief
!! Reads specific data for fourier interpolations
!-----------------------------------------------------------
SUBROUTINE READ_EXTRA_DETAILS_FOURIER1D(this,period,klist,is_even)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d),INTENT(INOUT):: this
   REAL(KIND=8),INTENT(IN) :: period
   INTEGER(KIND=4),DIMENSION(:) :: klist
   LOGICAL, INTENT(IN),OPTIONAL:: is_even
   ! Local variables
   INTEGER(KIND=4) :: i  ! counters
   ! Run section
   SELECT CASE(size(klist) == this%n)
      CASE(.FALSE.)
         WRITE(0,*) "READ_EXTRA_DETAILS_FOURIER1D ERR: mismatch between number of points and size of kpoints array"
         CALL EXIT(1)
      CASE(.TRUE.)
         ! do nothing
   END SELECT
   this%period=period
   ALLOCATE(this%klist(size(klist)))
   this%klist=klist
   SELECT CASE(present(is_even))
      CASE(.TRUE.)
         this%even=is_even
         DO i = 1, size(klist)
            SELECT CASE(klist(i)<0 .AND. is_even)
               CASE(.TRUE.)
                  WRITE(0,*) "READ_EXTRA_DETAILS_FOURIER1D ERR: Kpoints incompatible with selected parity"
                  CALL EXIT(1)
               CASE(.FALSE.)
                  ! do nothing
            END SELECT
         END DO
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE READ_EXTRA_DETAILS_FOURIER1D
!###########################################################
!# FUNCTION: termfou1d 
!###########################################################
!> @brief 
!! Value of the term of order kpoint of a one dimensional fourier
!! expansion. Positive values stand for even terms of the series.
!! Negative numbers stand for odd termns in the expansion
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d(kpoint,period,x)
   ! Initial declarations 
   USE CONSTANTS_MOD  
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x,period
   ! Run section
   SELECT CASE(kpoint)
      CASE(0)
         termfou1d=1.D0
      CASE(1 :) 
         termfou1d=dcos((2.D0*PI/period)*dfloat(kpoint)*x)
      CASE(: -1)
         termfou1d=dsin((2.D0*PI/period)*dfloat(-kpoint)*x)
      CASE DEFAULT
         WRITE(0,*) "Termfou1d ERR: This message was not supposed to be printed ever. Check the code"
         CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfou1d
!###########################################################
!# FUNCTION: termfou1d_dx
!###########################################################
!> @brief 
!! Value of the d(term)/dx of order kpoint of a one dimensional fourier
!! expansion. Positive values stand for even terms of the series.
!! Negative numbers stand for odd termns in the expansion
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_dx(kpoint,period,x)
   ! Initial declarations 
   USE CONSTANTS_MOD  
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x,period
   ! Run section
   SELECT CASE(kpoint)
      CASE(0)
         termfou1d_dx=0.D0
      CASE(1 :) 
         termfou1d_dx=-(2.D0*PI/period)*dfloat(kpoint)*dsin((2.D0*PI/period)*dfloat(kpoint)*x)
      CASE(: -1)
         termfou1d_dx=(2.D0*PI/period)*dfloat(-kpoint)*dcos((2.D0*PI/period)*dfloat(-kpoint)*x)
      CASE DEFAULT
         WRITE(0,*) "Termfou1d ERR: This message was not supposed to be printed ever. Check the code"
         CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfou1d_dx
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
! FUNCTION: is_even_fourier1d
!> @brief
!! Basic enquire function. If the function is even,
!! sinus elements vanish from expansion
!--------------------------------------------------
LOGICAL FUNCTION is_even_fourier1d(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d),INTENT(IN) :: this
   ! Run section ----------------
   is_even_fourier1d = this%even
   RETURN
END FUNCTION is_even_fourier1d
!###########################################################
!# SUBROUTINE: INTERPOL_FOURIER1D 
!###########################################################
!> @brief
!! Performs a generic fourier1d interpolation
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Mar/2014 
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INTERPOL_FOURIER1D(this)
   ! Initial declarations   
   USE MATHS_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d),INTENT(INOUT) :: this 
   ! Local variables
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: terms,inv_terms
   INTEGER(KIND=4) :: i,j ! counters
   ! Run section
   ALLOCATE(this%coeff(this%n))
   ALLOCATE(terms(this%n,this%n))
   ALLOCATE(inv_terms(this%n,this%n))
   SELECT CASE(allocated(this%klist))
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         WRITE(0,*) "INTERPOL_FOURIER ERR: Klist was not initialized before"
         WRITE(0,*) "INTERPOL_FOURIER_ERR: Use READ_EXTRA subroutine"
         CALL EXIT(1)
   END SELECT
   DO i = 1, this%n ! loop over eq for different points
      DO j = 1, this%n ! loop over coefficients
         terms(i,j)=termfou1d(this%klist(j),this%period,this%x(i))
      END DO
   END DO
   CALL INV_MTRX(this%n,terms,inv_terms)
   this%coeff=matmul(inv_terms,this%f)
   DEALLOCATE(terms)
   ! Check if there are extra functions to be interpolated
   SELECT CASE(allocated(this%extrafuncs))
      CASE(.TRUE.)
         ALLOCATE(this%extracoeff(this%n,size(this%extrafuncs(:,1))))
         DO i = 1, size(this%extrafuncs(:,1))
            this%extracoeff(:,i)=matmul(inv_terms,this%extrafuncs(i,:))
         END DO
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE INTERPOL_FOURIER1D
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
   REAL(KIND=8) :: r
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: terms
   ! Run section
   SELECT CASE(present(shift))
      CASE(.TRUE.)
         r=x+shift
      CASE(.FALSE.)
         r=x
   END SELECT
   ALLOCATE(terms(this%n))
   DO i = 1, this%n
      terms(i)=termfou1d(this%klist(i),this%period,r)
   END DO
   getvalue_fourier1d=dot_product(terms,this%coeff)
   DEALLOCATE(terms)
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
   REAL(KIND=8) :: r
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: terms
   ! Run section
   SELECT CASE(present(shift))
      CASE(.TRUE.)
         r=x+shift
      CASE(.FALSE.)
         r=x
   END SELECT
   ALLOCATE(terms(this%n))
   DO i = 1, this%n
      terms(i)=termfou1d_dx(this%klist(i),this%period,r)
   END DO
   getderiv_fourier1d=dot_product(terms,this%coeff)
   DEALLOCATE(terms)
   RETURN
END FUNCTION getderiv_fourier1d

END MODULE FOURIER1D_MOD
