!########################################################
! MODULE : FOURIER1D_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes Interpol1d_mod in its scope
!########################################################
MODULE FOURIER1D_MOD
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
   INTEGER(KIND=4) :: lastkpoint
   LOGICAL :: average_last=.FALSE.
   REAL(KIND=8) :: shift=0.D0
   CONTAINS
      PROCEDURE,PUBLIC,NON_OVERRIDABLE :: getshift => getshift_TERMCALCULATOR
      PROCEDURE,PUBLIC,NON_OVERRIDABLE :: getlastkpoint => getlastkpoint_TERMCALCULATOR
      PROCEDURE,PUBLIC,NON_OVERRIDABLE :: getaveragelast => getaveragelast_TERMCALCULATOR
      PROCEDURE(getvalue_termcalculator_example),PUBLIC,DEFERRED :: getvalue 
      PROCEDURE(getvalue_termcalculator_example),PUBLIC,DEFERRED :: getderiv 
END TYPE Termcalculator
!
ABSTRACT INTERFACE
   !###########################################################
   !# FUNCTION: getvalue_termcalculator_example 
   !###########################################################
   !> @brief 
   !! Just an example that child objects should override
   !-----------------------------------------------------------
   REAL(KIND=8) FUNCTION getvalue_termcalculator_example(this,kpoint,x) 
      IMPORT Termcalculator
      CLASS(Termcalculator),INTENT(IN):: this
      INTEGER(KIND=4),INTENT(IN) :: kpoint
      REAL(KIND=8),INTENT(IN) :: x
   END FUNCTION getvalue_termcalculator_example
   !-------------------------------------------------------------
END INTERFACE
!/////////////////////////////////////////////////////////////////////////////
! TYPE: FOURIER1D
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!----------------------------------------------------------------------------
TYPE,ABSTRACT,EXTENDS(Interpol1d):: FOURIER1D
   PRIVATE
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE,PUBLIC :: klist
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: coeff
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: extracoeff
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: extrafuncs
   CLASS(Termcalculator),ALLOCATABLE,PUBLIC :: term
   CONTAINS
      ! get block
      PROCEDURE,PUBLIC,NON_OVERRIDABLE :: getvalue => getvalue_FOURIER1D
      PROCEDURE,PUBLIC,NON_OVERRIDABLE :: getderiv => getderiv_FOURIER1D
      PROCEDURE,PUBLIC,NON_OVERRIDABLE :: getklist => getklist_FOURIER1D
      ! Set block
      PROCEDURE(SET_IRREP_FOURIER1D),PUBLIC,DEFERRED :: SET_IRREP 
      PROCEDURE,PUBLIC,NON_OVERRIDABLE :: SET_AVERAGE_LASTKPOINT => SET_AVERAGE_LASTKPOINT_FOURIER1D
      PROCEDURE,PUBLIC,NON_OVERRIDABLE :: SET_LASTKPOINT => SET_LASTKPOINT_FOURIER1D
      PROCEDURE,PUBLIC,NON_OVERRIDABLE :: SET_SHIFT => SET_SHIFT_FOURIER1D
      ! Tools
      PROCEDURE,PUBLIC,NON_OVERRIDABLE :: INTERPOL => INTERPOL_FOURIER1D
      PROCEDURE,PUBLIC,NON_OVERRIDABLE :: ADD_MOREFUNCS => ADD_MORE_FUNCS_FOURIER1D
      PROCEDURE,PUBLIC,NON_OVERRIDABLE :: GET_ALLFUNCS_AND_DERIVS => GET_ALLFUNC_AND_DERIVS_FOURIER1D
      ! Plotting tools
      PROCEDURE,PUBLIC,NON_OVERRIDABLE :: PLOTCYCLIC => PLOTCYCLIC_INTERPOL_FOURIER1D
      PROCEDURE,PUBLIC,NON_OVERRIDABLE :: PLOTCYCLIC_ALL => PLOTCYCLIC_ALL_INTERPOL_FOURIER1D
END TYPE FOURIER1D
ABSTRACT INTERFACE
   !###########################################################
   !# SUBROUTINE: SET_IRREP_FOURIER1D 
   !###########################################################
   !> @brief
   !! Sets irrep for this fourier series. Should be overriden by
   !! child non-abstract classes
   !-----------------------------------------------------------
   SUBROUTINE SET_IRREP_FOURIER1D(this,irrep)
      IMPORT Fourier1d
      CLASS(FOURIER1D),INTENT(INOUT):: this
      CHARACTER(LEN=2),INTENT(IN) :: irrep
   END SUBROUTINE SET_IRREP_FOURIER1D
END INTERFACE
!/////////////////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# FUNCTION: getklist_FOURIER1D 
!###########################################################
!> @brief 
!! Common get function. Gets Klist atribute
!-----------------------------------------------------------
FUNCTION getklist_FOURIER1D(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d),INTENT(IN):: this
   INTEGER(KIND=4),ALLOCATABLE,DIMENSION(:):: getklist_FOURIER1D
   ! Run section
   ALLOCATE(getklist_FOURIER1D(size(this%klist)))
   getklist_FOURIER1D=this%klist
   RETURN
END FUNCTION getklist_FOURIER1D
!###########################################################
!# FUNCTION: getshift_TERMCALCULATOR 
!###########################################################
!> @brief 
!! Common get function. Gets shift atribute
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getshift_TERMCALCULATOR(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Termcalculator),INTENT(IN):: this
   ! Run section
   getshift_TERMCALCULATOR=this%shift
   RETURN
END FUNCTION getshift_TERMCALCULATOR
!###########################################################
!# SUBROUTINE: SET_SHIFT_FOURIER1D 
!###########################################################
!> @brief
!! Common set subroutine. Sets shift atribute
!-----------------------------------------------------------
SUBROUTINE SET_SHIFT_FOURIER1D(this,shift)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d),INTENT(INOUT):: this
   REAL(KIND=8),INTENT(IN):: shift
   ! Run section
   this%term%shift=shift
   RETURN
END SUBROUTINE SET_SHIFT_FOURIER1D
!###########################################################
!# FUNCTION: getaveragelast_TERMCALCULATOR 
!###########################################################
!> @brief 
!! Common get function. Gets average_last atribute
!-----------------------------------------------------------
LOGICAL FUNCTION getaveragelast_TERMCALCULATOR(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Termcalculator),INTENT(IN):: this
   ! Run section
   getaveragelast_TERMCALCULATOR=this%average_last
   RETURN
END FUNCTION getaveragelast_TERMCALCULATOR
!###########################################################
!# FUNCTION: getlastkpoint_TERMCALCULATOR 
!###########################################################
!> @brief
!! Common get functions. Gets lastkpoint atribute.
!-----------------------------------------------------------
INTEGER(KIND=4) FUNCTION getlastkpoint_TERMCALCULATOR(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Termcalculator),INTENT(IN):: this
   ! Run section
   getlastkpoint_TERMCALCULATOR=this%lastkpoint
   RETURN
END FUNCTION getlastkpoint_TERMCALCULATOR
!###########################################################
!# SUBROUTINE: SET_LASTKPOINT_FOURIER1D 
!###########################################################
!> @brief
!! Common set function. Sets term%lastkpoint atribute
!-----------------------------------------------------------
SUBROUTINE SET_LASTKPOINT_FOURIER1D(this,lastkpoint)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d),INTENT(INOUT):: this
   INTEGER(KIND=4),INTENT(IN) :: lastkpoint
   ! Run section
   this%term%lastkpoint=lastkpoint
   RETURN
END SUBROUTINE SET_LASTKPOINT_FOURIER1D
!###########################################################
!# SUBROUTINE: SET_AVERAGE_LASTKPOINT_FOURIER1D 
!###########################################################
!> @brief
!! Common set function. Sets average_last_kpoint private atribute
!-----------------------------------------------------------
SUBROUTINE SET_AVERAGE_LASTKPOINT_FOURIER1D(this,bool)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d),INTENT(INOUT):: this
   LOGICAL,INTENT(IN) :: bool
   ! Run section
   this%term%average_last=bool
   RETURN
END SUBROUTINE SET_AVERAGE_LASTKPOINT_FOURIER1D
!###########################################################
!# SUBROUTINE: GET_ALLFUNC_AND_DERIVS_FOURIER1D
!###########################################################
!> @brief
!! Get value of the potential and derivs for an specific point x
!! for the main function and extra ones
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
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
END SUBROUTINE GET_ALLFUNC_AND_DERIVS_FOURIER1D
!###########################################################
!# SUBROUTINE: ADD_MORE_FUNCS_FOURIER1D
!###########################################################
!> @brief
!! Adds a new set of functions to interpolate at the same time
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE ADD_MORE_FUNCS_FOURIER1D(this,f)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(FOURIER1D),INTENT(INOUT) :: this
   REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: f
   ! Local variables
   INTEGER(KIND=4) :: nfuncs, ndata
   ! Run section
   nfuncs=size(f(:,1)) ! number of rows
   ndata=size(f(1,:)) ! number of columns
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
!# SUBROUTINE: INTERPOL_FOURIER1D 
!###########################################################
!> @brief
!! Performs a generic FOURIER1D interpolation
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
   CLASS(FOURIER1D),INTENT(INOUT) :: this
   ! Local variables
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: terms,inv_terms
   INTEGER(KIND=4) :: i,j ! counters
   ! Run section
   ALLOCATE(this%coeff(this%n))
   ALLOCATE(terms(this%n,this%n))
   ALLOCATE(inv_terms(this%n,this%n))
   DO i = 1, this%n ! loop over eq for different points
      DO j = 1, this%n ! loop over coefficients
         terms(i,j)=this%term%getvalue(this%klist(j),this%x(i))
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
!# FUNCTION: getvalue_FOURIER1D 
!###########################################################
!> @brief 
!! Get's fourier 1D interpolation value for a given point
!! inside the correct range
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 03/03/2014
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getvalue_FOURIER1D(this,x,shift)
   ! Initial declarations   
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(FOURIER1D),TARGET,INTENT(IN) :: this
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
   getvalue_FOURIER1D=dot_product(terms,this%coeff)
   DEALLOCATE(terms)
   RETURN
END FUNCTION getvalue_FOURIER1D
!###########################################################
!# FUNCTION: getderiv_FOURIER1D 
!###########################################################
!> @brief 
!! Get's derivative value at a given point X using the interpolation
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 03/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_FOURIER1D(this,x,shift)
   ! Initial declarations   
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(FOURIER1D),TARGET,INTENT(IN) :: this
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
   getderiv_FOURIER1D=dot_product(terms,this%coeff)
   DEALLOCATE(terms)
   RETURN
END FUNCTION getderiv_FOURIER1D
!######################################################################
! SUBROUTINE: PLOT_INTERPOL_FOURIER1D ################################
!######################################################################
!> @brief
!! Creates a data file called @b filename with the interpolation graphic of
!! this FOURIER1D type variable. The number of points in that graphic is defined
!! by @b npoints. Cannot be less than two. It also plots the first derivative.
!! The graphic goes from 0 to @f$2\pi@f$
!----------------------------------------------------------------------
SUBROUTINE PLOTCYCLIC_INTERPOL_FOURIER1D(this,npoints,filename)
   USE CONSTANTS_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   ! I/O variables -------------------------------
   INTEGER,INTENT(IN) :: npoints
   CLASS(FOURIER1D),INTENT(IN) :: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables -----------------------------
   INTEGER :: inpoints,ndelta
   REAL(KIND=8) :: delta, interval, x
   INTEGER :: i ! Counter
   CHARACTER(LEN=35), PARAMETER :: routinename = "PLOTCYCLIC_INTERPOL_FOURIER1D: " 
   ! Pointers ------------------------------------
   REAL(KIND=8):: xmin, xmax
   INTEGER(KIND=4):: n
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOTCYCLIC_INTERPOL_FOURIER1D ERR: Less than 2 points"
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
END SUBROUTINE PLOTCYCLIC_INTERPOL_FOURIER1D
!######################################################################
! SUBROUTINE: PLOTCYCLIC_ALL_INTERPOL_FOURIER1D ################################
!######################################################################
!> @brief
!! Creates a data file called @b filename with the interpolation graphic of
!! this FOURIER1D type variable. The number of points in that graphic is defined
!! by @b npoints. Cannot be less than two. It also plots the first derivative.
!! The graphic goes from 0 to @f$2\pi@f$
!----------------------------------------------------------------------
SUBROUTINE PLOTCYCLIC_ALL_INTERPOL_FOURIER1D(this,npoints,filename)
   USE CONSTANTS_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   ! I/O variables -------------------------------
   INTEGER,INTENT(IN) :: npoints
   CLASS(FOURIER1D),INTENT(IN) :: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables -----------------------------
   INTEGER :: inpoints,ndelta,nfunc
   REAL(KIND=8) :: delta, interval, x
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: f,dfdx
   INTEGER :: i ! Counter
   CHARACTER(LEN=35), PARAMETER :: routinename = "PLOTCYCLIC_INTERPOL_FOURIER1D: " 
   ! Pointers ------------------------------------
   REAL(KIND=8):: xmin, xmax
   INTEGER(KIND=4):: n
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOTCYCLIC_INTERPOL_FOURIER1D ERR: Less than 2 points"
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
END SUBROUTINE PLOTCYCLIC_ALL_INTERPOL_FOURIER1D

END MODULE FOURIER1D_MOD
