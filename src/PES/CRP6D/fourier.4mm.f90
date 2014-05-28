!########################################################
! MODULE : FOURIER1D_4MM_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes Interpol1d_mod in its scope
!########################################################
MODULE FOURIER1D_4MM_MOD
USE INTERPOL1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: Termcalculator
!> @brief
!! Abstract class to calculate terms of the series avoiding the use of
!! unnecessary switches
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014 
!> @version 1.0
!----------------------------------------------------------------
TYPE,ABSTRACT :: Termcalculator
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE(getvalue_termcalculator_example),PUBLIC,DEFERRED :: getvalue 
   PROCEDURE(getvalue_termcalculator_example),PUBLIC,DEFERRED :: getderiv 
END TYPE Termcalculator
!/////////////////////////////////////////////////////////////////
! TYPE: term_A1
!> @brief
!! Child class of abstract termcalculator. Strcture to avoid unnecessary switches
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Termcalculator) :: term_A1
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_4mm_A1
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_4mm_A1
END TYPE term_A1
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_4MM
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
TYPE,EXTENDS(Interpol1d):: Fourier1d_4mm
   PRIVATE
   CHARACTER(LEN=2) :: irrep
   REAL(KIND=8) :: period=2.D0*dacos(-1.D0) ! 2pi
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE :: klist
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: coeff
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: extracoeff
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: extrafuncs
   CLASS(Termcalculator),ALLOCATABLE :: term
   CONTAINS
      ! get block
      PROCEDURE,PUBLIC :: getvalue => getvalue_FOURIER1D_4MM
      PROCEDURE,PUBLIC :: getderiv => getderiv_FOURIER1D_4MM
      ! Tools
      PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_FOURIER1D_4MM
      PROCEDURE,PUBLIC :: ADD_MOREFUNCS => ADD_MORE_FUNCS_FOURIER1D_4MM
      PROCEDURE,PUBLIC :: GET_ALLFUNCS_AND_DERIVS => GET_ALLFUNC_AND_DERIVS_FOURIER1D_4MM
      ! Plotting tools
      PROCEDURE,PUBLIC :: PLOTCYCLIC => PLOTCYCLIC_INTERPOL_FOURIER1D_4MM
      PROCEDURE,PUBLIC :: PLOTCYCLIC_ALL => PLOTCYCLIC_ALL_INTERPOL_FOURIER1D_4MM
END TYPE FOURIER1D_4MM
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# FUNCTION: getvalue_termcalculator_example 
!###########################################################
!> @brief 
!! Just an example that child objects should override
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getvalue_termcalculator_example(this,kpoint,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Termcalculator),INTENT(IN):: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   WRITE(0,*) "getvalue_termcalculator_example ERR: This is just an example."
   CALL EXIT(1)
   RETURN
END FUNCTION getvalue_termcalculator_example
!###########################################################
!# FUNCTION: getderiv_termcalculator_example 
!###########################################################
!> @brief 
!! Just an example that child objects should override
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_termcalculator_example(this,kpoint,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Termcalculator),INTENT(IN):: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   WRITE(0,*) "getderiv_termcalculator_example ERR: This is just an example."
   CALL EXIT(1)
   RETURN
END FUNCTION getderiv_termcalculator_example
!###########################################################
!# FUNCTION: termfou1d_4mm_A1
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_4mm_A1(this,kpoint,x)
   ! Initial declarations 
   USE CONSTANTS_MOD  
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_A1),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   termfou1d_4mm_A1=dcos(dfloat(kpoint)*x)
   RETURN
END FUNCTION termfou1d_4mm_A1
!###########################################################
!# FUNCTION: termfou1d_dx_4mm_A1
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_dx_4mm_A1(this,kpoint,x)
   ! Initial declarations 
   USE CONSTANTS_MOD  
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_A1),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   termfou1d_dx_4mm_A1=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
   RETURN
END FUNCTION termfou1d_dx_4mm_A1
!###########################################################
!# SUBROUTINE: GET_ALLFUNC_AND_DERIVS_FOURIER1D_4MM
!###########################################################
!> @brief
!! Get value of the potential and derivs for an specific point x
!! for the main function and extra ones
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_ALLFUNC_AND_DERIVS_FOURIER1D_4MM(this,x,f,dfdx)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(FOURIER1D_4MM),INTENT(IN) :: this
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
         WRITE(0,*) "GET_ALLFUNCS_AND_DERIVS ERR: extra coefficients are not allocated"
         WRITE(0,*) "GET_ALLFUNCS_AND_DERIVS ERR: did you use ADD_MOREFUNCS and INTERPOL before this?"
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
      terms(i)=this%term%getvalue(this%klist(i),x)
      terms_dx(i)=this%term%getderiv(this%klist(i),x)
   END DO
   f(1)=dot_product(terms,this%coeff)
   dfdx(1)=dot_product(terms_dx,this%coeff)
   DO i = 2, nfuncs
      f(i)=dot_product(terms,this%extracoeff(:,i-1))
      dfdx(i)=dot_product(terms_dx,this%extracoeff(:,i-1))
   END DO
   RETURN
END SUBROUTINE GET_ALLFUNC_AND_DERIVS_FOURIER1D_4MM
!###########################################################
!# SUBROUTINE: ADD_MORE_FUNCS_FOURIER1D_4MM
!###########################################################
!> @brief
!! Adds a new set of functions to interpolate at the same time
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE ADD_MORE_FUNCS_FOURIER1D_4MM(this,f)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(FOURIER1D_4MM),INTENT(INOUT) :: this
   REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: f
   ! Local variables
   INTEGER(KIND=4) :: nfuncs, ndata
   ! Run section
   nfuncs=size(f(:,1)) ! number of rows
   ndata=size(f(1,:)) ! number of columns
   SELECT CASE(ndata == this%n)
      CASE(.FALSE.)
         WRITE(0,*) "ADD_MORE_FUNCS_FOURIER1D_4MM ERR: size mismatch between extra functions and the original one"
         CALL EXIT(1)
      CASE(.TRUE.)
         ! donothing
   END SELECT
   ALLOCATE(this%extrafuncs(nfuncs,ndata))
   this%extrafuncs=f
   RETURN
END SUBROUTINE ADD_MORE_FUNCS_FOURIER1D_4MM
!###########################################################
!# SUBROUTINE: INTERPOL_FOURIER1D_4MM 
!###########################################################
!> @brief
!! Performs a generic FOURIER1D_4MM interpolation
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Mar/2014 
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INTERPOL_FOURIER1D_4MM(this,irrep)
   ! Initial declarations   
   USE MATHS_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(FOURIER1D_4MM),INTENT(INOUT) :: this
   CHARACTER(LEN=2),INTENT(IN) :: irrep
   ! Local variables
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: terms,inv_terms
   INTEGER(KIND=4) :: i,j ! counters
   ! Run section
   ALLOCATE(this%coeff(this%n))
   ALLOCATE(terms(this%n,this%n))
   ALLOCATE(inv_terms(this%n,this%n))
   ALLOCATE(this%klist(this%n))
   this%irrep=irrep
   SELECT CASE(this%irrep)
      CASE("A1")
         FORALL(i=1:this%n) this%klist(i)=(i-1)*4
         ALLOCATE(Term_A1::this%term)
         DO i = 1, this%n ! loop over eq for different points
            DO j = 1, this%n ! loop over coefficients
               terms(i,j)=this%term%getvalue(this%klist(j),this%x(i))
            END DO
         END DO
      CASE DEFAULT
         WRITE(0,*) "INTERPOL_FOURIER1D_4MM ERR: irrep used is not implemented or does not exist"
         WRITE(0,*) "List of irreps implemented: A1"
         CALL EXIT(1)
   END SELECT
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
END SUBROUTINE INTERPOL_FOURIER1D_4MM
!###########################################################
!# FUNCTION: getvalue_FOURIER1D_4MM 
!###########################################################
!> @brief 
!! Get's fourier 1D interpolation value for a given point
!! inside the correct range
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 03/03/2014
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getvalue_FOURIER1D_4MM(this,x,shift)
   ! Initial declarations   
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(FOURIER1D_4MM),TARGET,INTENT(IN) :: this
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
      terms(i)=this%term%getvalue(this%klist(i),r)
   END DO
   getvalue_FOURIER1D_4MM=dot_product(terms,this%coeff)
   DEALLOCATE(terms)
   RETURN
END FUNCTION getvalue_FOURIER1D_4MM
!###########################################################
!# FUNCTION: getderiv_FOURIER1D_4MM 
!###########################################################
!> @brief 
!! Get's derivative value at a given point X using the interpolation
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 03/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_FOURIER1D_4MM(this,x,shift)
   ! Initial declarations   
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(FOURIER1D_4MM),TARGET,INTENT(IN) :: this
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
      terms(i)=this%term%getvalue(this%klist(i),r)
   END DO
   getderiv_FOURIER1D_4MM=dot_product(terms,this%coeff)
   DEALLOCATE(terms)
   RETURN
END FUNCTION getderiv_FOURIER1D_4MM
!######################################################################
! SUBROUTINE: PLOT_INTERPOL_FOURIER1D_4MM ################################
!######################################################################
!> @brief
!! Creates a data file called @b filename with the interpolation graphic of
!! this FOURIER1D_4MM type variable. The number of points in that graphic is defined
!! by @b npoints. Cannot be less than two. It also plots the first derivative.
!! The graphic goes from 0 to @f$2\pi@f$
!----------------------------------------------------------------------
SUBROUTINE PLOTCYCLIC_INTERPOL_FOURIER1D_4MM(this,npoints,filename)
   USE CONSTANTS_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   ! I/O variables -------------------------------
   INTEGER,INTENT(IN) :: npoints
   CLASS(FOURIER1D_4MM),INTENT(IN) :: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables -----------------------------
   INTEGER :: inpoints,ndelta
   REAL(KIND=8) :: delta, interval, x
   INTEGER :: i ! Counter
   CHARACTER(LEN=35), PARAMETER :: routinename = "PLOTCYCLIC_INTERPOL_FOURIER1D_4MM: " 
   ! Pointers ------------------------------------
   REAL(KIND=8):: xmin, xmax
   INTEGER(KIND=4):: n
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOTCYCLIC_INTERPOL_FOURIER1D_4MM ERR: Less than 2 points"
      STOP
   END IF
   !
   n= this%n
   xmin=0
   xmax=2.D0*PI
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
#ifdef DEBUG   
   CALL DEBUG_WRITE(routinename,filename," file created")
#endif
   CLOSE(11)
END SUBROUTINE PLOTCYCLIC_INTERPOL_FOURIER1D_4MM
!######################################################################
! SUBROUTINE: PLOTCYCLIC_ALL_INTERPOL_FOURIER1D_4MM ################################
!######################################################################
!> @brief
!! Creates a data file called @b filename with the interpolation graphic of
!! this FOURIER1D_4MM type variable. The number of points in that graphic is defined
!! by @b npoints. Cannot be less than two. It also plots the first derivative.
!! The graphic goes from 0 to @f$2\pi@f$
!----------------------------------------------------------------------
SUBROUTINE PLOTCYCLIC_ALL_INTERPOL_FOURIER1D_4MM(this,npoints,filename)
   USE CONSTANTS_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   ! I/O variables -------------------------------
   INTEGER,INTENT(IN) :: npoints
   CLASS(FOURIER1D_4MM),INTENT(IN) :: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables -----------------------------
   INTEGER :: inpoints,ndelta,nfunc
   REAL(KIND=8) :: delta, interval, x
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: f,dfdx
   INTEGER :: i ! Counter
   CHARACTER(LEN=35), PARAMETER :: routinename = "PLOTCYCLIC_INTERPOL_FOURIER1D_4MM: " 
   ! Pointers ------------------------------------
   REAL(KIND=8):: xmin, xmax
   INTEGER(KIND=4):: n
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOTCYCLIC_INTERPOL_FOURIER1D_4MM ERR: Less than 2 points"
      STOP
   END IF
   !
   n= this%n
   xmin=0
   xmax=2.D0*PI
   nfunc=size(this%extrafuncs(:,1))+1
   ALLOCATE(f(nfunc))
   ALLOCATE(dfdx(nfunc))
   !
   interval=xmax-xmin
   inpoints=npoints-2
   ndelta=npoints-1
   delta=interval/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   CALL this%GET_ALLFUNCS_AND_DERIVS(xmin,f,dfdx)
   WRITE(11,*) xmin,f(:),dfdx(:)
   DO i=1, inpoints
      x=xmin+(dfloat(i)*delta)
      CALL this%GET_ALLFUNCS_AND_DERIVS(x,f,dfdx) 
      WRITE(11,*) x,f(:),dfdx(:)
   END DO
   CALL this%GET_ALLFUNCS_AND_DERIVS(xmax,f,dfdx) 
   WRITE(11,*) xmax,f(:),dfdx(:)
#ifdef DEBUG   
   CALL DEBUG_WRITE(routinename,filename," file created")
#endif
   CLOSE(11)
END SUBROUTINE PLOTCYCLIC_ALL_INTERPOL_FOURIER1D_4MM

END MODULE FOURIER1D_4MM_MOD
