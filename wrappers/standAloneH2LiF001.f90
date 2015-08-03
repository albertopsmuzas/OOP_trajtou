!#########################################################
! MODULE MATHS
!> @brief
!! Contains some mathematical operations
!##########################################################
MODULE MATHS_MOD
! Initial declarations
implicit none
CONTAINS
!####################################################################
! SUBROUTINE: ORDER #################################################
!####################################################################
!> @brief
!! This subroutine orders an array using the insertion ordering algorithm while
!! a second array is ordered following the same permutations as the prior one.
!
!> @details
!! Let's define some notation:
!! - @b arregl1 has elements @f$v_{i}@f$
!! - @b arregl2 has elements @f$w_{i}@f$
!! - If @f$<@f$ and @f$\spadesuit@f$ are order relationships, it holds:
!!   @f$w_{i}\spadesuit w_{i+1} \Leftrightarrow v_{i}<v_{i+1}@f$
!
!> @param[in,out] arreg1 - This array is ordered so that @f$v_{i}<v_{i+1}@f$
!> @param[in,out] arreg2 - This array is ordered so that @f$w_{i}\spadesuit w_{i+1}@f$
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 27/Jan/2014
!> @version 1.0
! -----------------------------------------------------
SUBROUTINE ORDER(ARREG1,ARREG2)
   IMPLICIT NONE
   ! I/O variables ---------------------------
   REAL*8, DIMENSION(:), INTENT(INOUT) :: ARREG1, ARREG2
   ! Local variables -------------------------
   INTEGER :: NELEM
   LOGICAL :: CLAVE
   INTEGER :: I, K, POS
   REAL*8 :: AUX1, AUX2
   ! FIRE IN THE HOLE! ------------------------
   nelem=size(arreg1)
   IF (nelem/=size(arreg2)) THEN
      WRITE(0,*) "ORDER ERR: dimension mismatch in arreg1 and arreg2"
      WRITE(0,*) "arreg1: ", nelem
      WRITE(0,*) "arreg2: ", size(arreg2)
      CALL EXIT(1)
   END IF
   !
   SELECT CASE (nelem)
      CASE(: 1)
         WRITE(0,*) "ORDER ERR: Less than 2 elements inside arrays."
         CALL EXIT(1)
      CASE DEFAULT
         ! do nothing
   END SELECT
   !
   DO I=2,NELEM
      K = I
      AUX1 = ARREG1(K)
      AUX2 = ARREG2(K)
      CLAVE = .FALSE.
      DO WHILE((K.GT.1).AND.(.NOT.CLAVE))
         IF (ARREG1(K-1).GT.AUX1) THEN
            ARREG1(K) = ARREG1(K-1)
            ARREG2(K) = ARREG2(K-1)
            K = K-1
         ELSE
            CLAVE = .TRUE.
         ENDIF
      END DO
   POS = K
   ARREG1(POS) = AUX1
   ARREG2(POS) = AUX2
   ENDDO
END SUBROUTINE ORDER
!###########################################################
!# SUBROUTINE: ORDER_VECT
!###########################################################
!> @brief
!! This subroutine orders an order 1 array from low to high values
!> @param[in,out] arreg1 - Vector to order
!!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 31/Jan/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE ORDER_VECT(ARREG1)
   IMPLICIT NONE
   ! I/O variables ---------------------------
   REAL(KIND=8),DIMENSION(:),INTENT(INOUT) :: ARREG1
   ! Local variables -------------------------
   INTEGER(KIND=4) :: NELEM
   LOGICAL :: CLAVE
   INTEGER :: I, K, POS
   REAL*8 :: AUX1
   ! FIRE IN THE HOLE! ------------------------
   nelem=size(arreg1)
   IF (NELEM.LT.2)  THEN
      WRITE(0,*) "ORDER ERR: Less than 2 elements inside arrays."
      STOP
   END IF
   DO I=2,NELEM
      K = I
      AUX1 = ARREG1(K)
      CLAVE = .FALSE.
      DO WHILE((K.GT.1).AND.(.NOT.CLAVE))
         IF (ARREG1(K-1).GT.AUX1) THEN
            ARREG1(K) = ARREG1(K-1)
            K = K-1
         ELSE
            CLAVE = .TRUE.
         ENDIF
      END DO
   POS = K
   ARREG1(POS) = AUX1
   ENDDO
   RETURN
END SUBROUTINE ORDER_VECT
!####################################################################
! SUBROUTINE: SYMMETRIZE ############################################
!####################################################################
!> @brief
!! This subroutine takes  couples of values @f$(x_{i},F(x_{i}))@f$ and
!! makes them symmetric respect to a given value @f$x_{0}@f$. See details.
!
!> @details
!> - Actually, this routine projects points around a given value @f$x_{0}@f$ so
!! that the final distribution of points is symmetric respect that value.
!! - This procedure is only done if @f$F(x_{i})>f_{0}@f$, where @f$f_{0}@f$ is a
!! threshold value given to the routine.
!
!> @param[in] n - Number of initial @f$(x_{i},F(x_{i}))@f$ couples
!> @param[in,out] x,v - Couples @f$(x_{i},F(x_{i}))@f$
!> @param{in} zero - @f$x_{0}@f$ value
!> @param[in] vtop - @f$f_{0}@f$ threshold value
!> @param[out] nsym - Number of final points after symmetrization procedure
!> @param[out] xsym, vsym - Collection of @f$(x_{i},F(x_{i}))@f$ couples after
!!                          symmetrization procedure
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 03/Feb/2014
!> @version 1.0
! -----------------------------------------------------
SUBROUTINE SYMMETRIZE(n,x,v,zero,vtop,nsym,xsym,vsym)
   IMPLICIT NONE
   ! I/O variables ------------------------
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(OUT) :: nsym
   REAL*8, INTENT(IN) :: vtop, zero
   REAL*8, INTENT(INOUT), DIMENSION(n) :: x, v
   REAL*8, INTENT(OUT), DIMENSION(:), ALLOCATABLE :: xsym, vsym
   ! Internal variables -------------------
   INTEGER :: i,j, k ! Counters
   INTEGER :: npoints, nnonred
   REAL*8, DIMENSION(:), ALLOCATABLE :: auxx, auxv
   REAL*8, DIMENSION(:), ALLOCATABLE :: nonredx, nonredv
   INTEGER, DIMENSION(:), ALLOCATABLE :: redundant
   ! HEY, HO! LET'S GO! -------------------
   CALL ORDER(x,v) ! Order the initial array
   npoints=0
   ! Checking number of points that satisfies F(x) > vtop
   DO i=1,n
      IF (v(i).GE.vtop) THEN
         WRITE(*,*) "SYMMETRIZE: Pair: ",i, x(i), v(i)
         npoints=npoints+1
      END IF
   END DO
   WRITE(*,*) "SYMMETRIZE: Symmetrize subroutine was invoked. Take care of what you are doing."
   WRITE(*,*) "SYMMETRIZE: Found ", npoints, "points to symmetrize. F(X) > ", vtop
   ! Check if some of these points are redundant
   ALLOCATE(auxx(1:npoints))
   ALLOCATE(auxv(1:npoints))
   ALLOCATE(redundant(1:npoints))
   DO i=1, npoints
      redundant(i)=0
   END DO
   k=0
   DO i=1,n
      IF (v(i).GE.vtop) THEN
         k=k+1
         auxv(k)=v(i)
         auxx(k)=x(i)
      END IF
   END DO
   k=0
   DO i=1, npoints
      IF (redundant(i).EQ.0) THEN
         DO j=i+1, npoints
            IF (DABS(auxx(i)-zero).EQ.DABS(auxx(j)-zero)) THEN
               k=k+1
               redundant(i)=k
               redundant(j)=k
            END IF
         END DO
      END IF
   END DO
   WRITE(*,*) "SYMMETRIZE: Redundance mapping:"
   DO i=1, npoints
      WRITE(*,*) i, auxx(i), auxv(i), redundant(i)
   END DO

!debug   WRITE(*,*) "SYMMETRIZE: ", k, "redundant pairs found."
   ! Store non-redundant points
   nnonred=npoints-k*2
   ALLOCATE(nonredx(1:nnonred))
   ALLOCATE(nonredv(1:nnonred))
   k=0
   DO i=1, npoints
      IF (redundant(i).EQ.0) THEN
         k=k+1
         nonredx(k)=auxx(i)
         nonredv(k)=auxv(i)
      END IF
   END DO
   ! Adding the new points
   nsym=n+nnonred
   ALLOCATE(xsym(1:nsym))
   ALLOCATE(vsym(1:nsym))
   k=0
   DO i=1,nnonred
      IF (nonredx(i)-zero.LE.0) THEN
         k=k+1
         xsym(k)=nonredx(i)+2.D0*DABS(nonredx(i)-zero)
         vsym(k)=nonredv(i)
      ELSE IF (nonredx(i)-zero.GT.0) THEN
         k=k+1
         xsym(k)=nonredx(i)-2.D0*DABS(nonredx(i)-zero)
         vsym(k)=nonredv(i)
      END IF
   END DO
   DO i=nnonred+1,nsym
      xsym(i)=x(i-nnonred)
      vsym(i)=v(i-nnonred)
   END DO
   ! The vector is not ordered
   CALL ORDER(xsym,vsym)
   WRITE(*,*) "SYMMETRIZE: New points added."
   RETURN
END SUBROUTINE SYMMETRIZE
!###########################################################################
! SUBROUTINE: TRIDIA #######################################################
!###########################################################################
!> @brief
!! Solves a tridiagonal system. Thomas algorithm.
!! Conserve the system matrix.
!
!> @param[in] n dimensions of the arrays. number of equations
!> @param[in] a(n) - subdiagonal
!> @param[in] b(n) - diagonal
!> @param[in] c(n) - superdiagonal
!> @param[in] d(n) - right part of the equation
!> @param[out] x(n) - answer
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 11/Feb/2014
!> @version 1.0
! -------------------------------------------------------------------
SUBROUTINE TRIDIA(n,a,b,c,d,x)
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: n
        REAL*8,DIMENSION(n),INTENT(IN) :: a,b,c,d
        REAL*8,DIMENSION(n),INTENT(OUT) :: x
        REAL*8,DIMENSION(n) :: cp,dp
        REAL*8 :: m
        INTEGER :: i
! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
   DO i = 2,n
      m = b(i)-cp(i-1)*a(i)
      cp(i) = c(i)/m
      dp(i) = (d(i)-dp(i-1)*a(i))/m
   END DO
! initialize x
   x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
   DO i = n-1, 1, -1
      x(i) = dp(i)-cp(i)*x(i+1)
   END DO
END SUBROUTINE TRIDIA
!#####################################################################
! SUBROUTINE : INV_MTRX ##############################################
!#####################################################################
!> @brief
!! Invert a square matrix with Gauss method. @b mtrx is not destroyed during  the procedure
!
!> @param[in] n - Order of the matrix
!> @param[in] mtrx(n,n) - Initial matrix
!> @param[out] i_mtrx(n,n) - Inverse matrix
!
!> @warning
!! - This algorithm is inefficient for large matrices or sparse matrices. Use
!!   it wisely
!---------------------------------------------------------------------
SUBROUTINE INV_MTRX (n, mtrx, i_mtrx)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   REAL(8), DIMENSION(n,n), INTENT(IN):: mtrx
   REAL(8), DIMENSION(n,n), INTENT(OUT):: i_mtrx
   ! Local Variables ------------------------
   REAL*8, DIMENSION(n,n) :: B, A
   REAL*8, DIMENSION(n) :: temp
   INTEGER, DIMENSION(n) :: ipvt
   INTEGER, DIMENSION(1) :: imax
   REAL(8) :: c, d
   INTEGER :: i, j, k, m
   ! HEY, HO! LET'S GO! ---------------------
   A = mtrx
   B = A
   ipvt = (/ (i, i = 1, n) /)
   DO k = 1,n
      imax = MAXLOC(ABS(b(k:n,k)))
      m = k-1+imax(1)
      IF (m /= k) THEN
         ipvt( (/m,k/) ) = ipvt( (/k,m/) )
         B((/m,k/),:) = B((/k,m/),:)
      END IF
      d = 1/B(k,k)
      temp = B(:,k)
      DO j = 1, n
         c = B(k,j)*d
         B(:,j) = B(:,j)-temp*c
         B(k,j) = c
      END DO
      B(:,k) = temp*(-d)
      B(k,k) = d
   END DO
   A(:,ipvt) = B
   i_mtrx = A
END SUBROUTINE INV_MTRX
!###################################################################
! FUNCTION: cartesianPeriodizer
!###################################################################
!> @brief
!! Given a function f(x), and the cartesian periodizer function g(x,T),
!! it stands that f( g(x,T) ) is the f(x)'s segment between 0<=x<=T repeated
!! indefinitely
!-------------------------------------------------------------------
function cartesianPeriodizer(x,T) result(y)
   ! initial declarations
   implicit none
   ! I/O variables
   real(kind=8),intent(in):: x, T
   ! Dummy function variable
   real(kind=8):: y
   ! Local variables
   real(kind=8),parameter:: pi=dacos(-1.d0)
   ! Run section .....................................
   y=T/2.d0-(T/pi)*datan( 1.d0/dtan(x*pi/T) )
   return
end function cartesianPeriodizer
!###################################################################
! FUNCTION: polarPeriodizer
!###################################################################
!> @brief
!! Defined as cartesian periodizer function g(x,T) with some changes:
!! x=x-x0 and T=2pi/N. Designed to be used with figures with polar
!! symmetry and whole number of lobes.
!> @details
!! - Extracted from E. Chicurel-Uziel/Computer Aided Geometric Design/21(2004) 23-42
!-------------------------------------------------------------------
function polarPeriodizer(theta,N,theta0) result(y)
   ! initial declarations
   implicit none
   ! I/O variables
   real(kind=8),intent(in):: theta,theta0
   integer(kind=4),intent(in):: N
   ! dummy function out variable
   real(kind=8):: y
   ! Local variables
   real(kind=8),parameter:: pi=dacos(-1.d0)
   ! Run section ..............................
   y=pi/dfloat(N)-(2.d0/dfloat(N))*datan( 1.d0/dtan(dfloat(N)*(theta-theta0)/2.d0) )
   return
end function polarPeriodizer
!####################################################################
! FUNCTION: radialPolygonEquation
!####################################################################
!> @brief
!!
!--------------------------------------------------------------------
function radialPolygonEquation(theta,r,N,theta0) result(y)
   ! initial declarations
   implicit none
   ! I/O variables
   real(kind=8),intent(in):: theta,r,theta0
   integer(kind=4),intent(in):: N
   ! Dummy output variable
   real(kind=8):: y
   ! Local variable
   real(kind=8):: beta
   real(kind=8),parameter:: pi=dacos(-1.d0)
   ! Run section ......................................
   beta=pi*(dfloat(N-2)/dfloat(2*N))
   y=r*dtan(beta)/(dsin(polarPeriodizer(theta,N,theta0))+dtan(beta)*dcos(polarPeriodizer(theta,N,theta0)))
   return
end function radialPolygonEquation
!####################################################################
! FUNCTION: parametricPolygonEquation
!####################################################################
!> @brief
!!
!--------------------------------------------------------------------
subroutine parametricPolygonEquation(theta,r,N,theta0,x0,x)
   ! initial declarations
   implicit none
   ! I/O variables
   real(kind=8),intent(in):: theta,theta0,r
   real(kind=8),dimension(2),intent(in):: x0
   integer(kind=4),intent(in):: N
   real(kind=8),dimension(2),intent(out):: x
   ! Run section ...................................
   x(1)=x0(1)+radialPolygonEquation(theta,r,N,theta0)*dcos(theta)
   x(2)=x0(2)+radialPolygonEquation(theta,r,N,theta0)*dsin(theta)
   return
end subroutine parametricPolygonEquation
!##################################################################
! FUNCTION: checkLoschianOrder
!##################################################################
!> @brief
!! Given an integer number, check the order of the nearest Loschian
!! number. It will be used to give diffraction order for hexagonal
!! lattices.
!> @details
!! - It uses an array of first 64 Loschian numbers. If diffraction is
!!   too high, it may be insufficient and a longer series should be given
!!   and compiled. You can detect this unconviniency if you get too many
!!   peaks with order 63 and delta K changes within this order
!------------------------------------------------------------------
function checkLoschianOrder(num) result(order)
   implicit none
   ! I/O variables
   integer(kind=4),intent(in):: num
   ! Dummy function output variable
   integer(kind=4):: order
   ! Local variables
   integer(kind=4),dimension(1):: arrorder
   integer(kind=4),dimension(0:63),parameter:: loschianNum=[   0,  1,  3,  4,  7,  9, 12, 13, 16, 19, 21, 25, 27, 28, 31, 36,&
                                                              37, 39, 43, 48, 49, 52, 57, 61, 63, 64, 67, 73, 75, 76, 79, 81,&
                                                              84, 91, 93, 97,100,103,108,109,111,112,117,121,124,127,129,133,&
                                                             139,144,147,148,151,156,157,163,169,171,172,175,181,183,189,192]
   integer(kind=4),dimension(0:63):: intDistance
   ! Run section ...............................................................................................................
   intDistance(:)=( loschianNum(:) -num )**2
   arrorder(:)=minloc( array=intDistance(:) )
   order=arrorder(1)-1
   return
end function checkLoschianOrder
END MODULE MATHS_MOD
!#########################################################
! MODULE: UNITS_MOD
!> @brief
!! Provides tools to manage change of units
!
!> @details
!! - Routines inside the program (not in this module) will
!!   always use a.u. for internal calculations.
!! - Magnitudes read from input files that could be deffined with
!!   different units, should be declared as a "Quantity" subtype.
!
!> @todo
!! - Write subroutines to go from a.u. to other units. This way
!!   output units in oautput files could be managed.
!##########################################################
MODULE UNITS_MOD
IMPLICIT NONE
!//////////////////////////////////////////////////////////////////////
! TYPE: Quantity
! --------------
!> @brief
!! A general physical quantity class
!
!> @section* sec DETAILS
!> @param units - String (maximum 10 characters) with units shortcut
!> @param mag - Number, magnitude of the quantity
!
!> @author A.P. Muzas - alberto.muzas@uam.es
!> @date 21/Jan/2014
!> @version 1.0
!> @see length, mass, energy, angle, time
!----------------------------------------------------------------------
TYPE :: Quantity
   PRIVATE
   CHARACTER(LEN=10) :: units
   REAL(KIND=8) :: mag
CONTAINS
   PROCEDURE,PUBLIC :: READ => READ_QUANTITY_FROM_ARGUMENTS
   PROCEDURE,PUBLIC :: getvalue => get_magnitude_quantity
   PROCEDURE,PUBLIC :: getunits => get_units_quantity
END TYPE Quantity
!//////////////////////////////////////////////////////////////////////
! SUBTYPE: Length
! ---------------
!> @brief
!! Length quantity subtype
!
!> - Supported units: @a au, @a angst
!
!> @author A.P. Muzas - alberto.muzas@uam.es
!> @date 06/Nov/2013
!> @version 1.0
!> @see  quantity
!----------------------------------------------------------------------
TYPE, EXTENDS(Quantity) :: Length
CONTAINS
   PROCEDURE,PUBLIC :: TO_STD => LENGTH_AU
   PROCEDURE,PUBLIC :: TO_ANGST => TO_ANGST_LENGTH
END TYPE Length
!//////////////////////////////////////////////////////////////////////
! SUBTYPE: Temperature
! --------------------
!> @brief
!! Temperature quantity subtype
!
!> - Supported units: @a Celsius, @a Kelvin, @a Farenheit
!
!> @author A.P. Muzas - alberto.muzas@uam.es
!> @date Feb/2015
!> @version 1.0
!> @see  quantity
!----------------------------------------------------------------------
TYPE,EXTENDS(Quantity) :: Temperature
CONTAINS
   PROCEDURE,PUBLIC:: TO_STD => TO_KELVIN_UNITS
END TYPE Temperature
!//////////////////////////////////////////////////////////////////////
! SUBTYPE: Energy
! ---------------
!> @brief
!! Stands for an energy quantity
!
!> - Supported units: @a au, @a ev, @a kjmol, @a kcalmol
!
!> @author A.P. Muzas - alberto.muzas@uam.es
!> @date 06/Nov/2013
!> @version 1.0
!> @see quantity
!----------------------------------------------------------------------
TYPE, EXTENDS(Quantity) :: Energy
CONTAINS
   PROCEDURE,PUBLIC :: TO_STD => ENERGY_AU
END TYPE Energy
!//////////////////////////////////////////////////////////////////////
! SUBTYPE: Mass
! -------------
!> @brief
!! Stands for a mass quantity
!
!> - Supported units: @a au, @a hmass, @a dmass, @a pmass (hydrogen mass, deuterium mass, proton mass)
!
!> @author A.P. Muzas - alberto.muzas@uam.es
!> @date 06/Nov/2013
!> @version 1.0
!> @see quantity
!----------------------------------------------------------------------
TYPE, EXTENDS(Quantity) :: Mass
CONTAINS
   PROCEDURE,PUBLIC :: TO_STD => MASS_AU
END TYPE Mass
!//////////////////////////////////////////////////////////////////////
! SUBTYPE: Time
! -------------
!> @brief
!! Stands for a time quantity
!
!> - Supported units: @a au, @a ps, @a fs
!
!> @author A.P. Muzas - alberto.muzas@uam.es
!> @date 06/Nov/2013
!> @version 1.0
!> @see quantity
!----------------------------------------------------------------------
TYPE, EXTENDS(Quantity) :: Time
CONTAINS
   PROCEDURE,PUBLIC :: TO_STD => TIME_AU
END TYPE Time
!//////////////////////////////////////////////////////////////////////
! SUBTYPE: Angle
! --------------
!> @brief
!! Stands for an angle quantity
!
!! - Supported units: @a rad, @a deg
!
!> @author A.P. Muzas - alberto.muzas@uam.es
!> @date 06/Nov/2013
!> @version 1.0
!> @see quantity
!----------------------------------------------------------------------
TYPE, EXTENDS(Quantity) :: angle
CONTAINS
   PROCEDURE,PUBLIC :: TO_STD => TO_RAD
   PROCEDURE,PUBLIC :: TO_DEG
END TYPE
!//////////////////////////////////////////////////////////////////////
! Conversion factors:
REAL(KIND=8),PARAMETER:: au2ev = 27.21138386D0
REAL(KIND=8),PARAMETER:: au2kcalmol = 627.503D0
REAL(KIND=8),PARAMETER:: au2kjmol = 2.6255D3
REAL(KIND=8),PARAMETER:: au2angst = 0.52917720859D0
REAL(KIND=8),PARAMETER:: au2fs = 0.02418884326505D0
REAL(KIND=8),PARAMETER:: au2ps = 2.418884326505D-5
REAL(KIND=8),PARAMETER:: hmass2au = 1837.15264409D0
REAL(KIND=8),PARAMETER:: dmass2au = 3671.482934845d0
REAL(KIND=8),PARAMETER:: pmass2au = 1836.15267247d0
real(kind=8),parameter:: dalton2au=1822.88848325d0
real(kind=8),parameter:: kelvinParam=273.15d0
real(kind=8),parameter:: fahrenheitParam=459.67d0
real(kind=8),parameter:: boltzmann=8.617332478d-5/au2ev
REAL(KIND=8),PARAMETER:: pi=dacos(-1.D0)
!//////////////////////////////////////////////////////////////////////

! MODULE CONTAINS:
CONTAINS
!###########################################################
!# SUBROUTINE: READ_QUANTITY_FROM_ARGUMENTS ################
!###########################################################
!> @brief
!! Reads a quantity from arguments
!
!> \param[out] quant - A quantity type. It should be a Quantity subtype, in fact
!> \param[in] mag - Magnitude to be stored in quantity%mag
!> \param[in] units - Units key to be stored in quantity%units.
!!                    It should be a supported unit.
!
!> @author A.P. Muzas - alberto.muzas@uam.es
!> @date 06/Nov/2013
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE READ_QUANTITY_FROM_ARGUMENTS(quant,mag,units)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Quantity),INTENT(OUT) :: quant
   REAL(KIND=8),INTENT(IN) :: mag
   CHARACTER(LEN=*),INTENT(IN) :: units
   ! Run section
   ! Check variable type
   SELECT TYPE (quant)
   TYPE IS (Quantity)
      WRITE(0,*) "READ_QUANT_FROM_ARGUMENTS ERR: Pure Quantities're not allowed. &
      & Declare it with a specific sub-type instead (length, time, etc.)"
      CALL EXIT(1)
   END SELECT
   quant%mag=mag
   quant%units=units
   RETURN
END SUBROUTINE READ_QUANTITY_FROM_ARGUMENTS
!###########################################################
!# FUNCTION: get_magnitude_quantity ########################
!###########################################################
! - Typical get function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION get_magnitude_quantity(quant)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Quantity),INTENT(IN) :: quant
   ! Run section
   get_magnitude_quantity=quant%mag
   RETURN
END FUNCTION get_magnitude_quantity
!###########################################################
!# FUNCTION: get_units_quantity ############################
!###########################################################
! - Typical get function
!-----------------------------------------------------------
CHARACTER(LEN=10) FUNCTION get_units_quantity(quant)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Quantity),INTENT(IN) :: quant
   ! Run section
   get_units_quantity=quant%units
   RETURN
END FUNCTION get_units_quantity
!###########################################################
!# SUBROUTINE: GO TO RAD ####################################
!############################################################
!> @brief
!! Manages angle units. Changes them to radians
!
!> @param[in,out] this - Angle subtype variable. Can be omitted.
!
!> @author A.P. Muzas - alberto.muzas@uam.es
!> @date 06/Nov/2013
!> @version 1.0
!------------------------------------------------------------
SUBROUTINE TO_RAD(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(angle), INTENT(INOUT) :: this
   ! GO, GO, GO !!!------
   IF (this%units.EQ."rad") THEN
      RETURN
   ELSE IF (this%units.EQ."deg") THEN
      this%mag = this%mag*(pi/180.D0)
   ELSE
      WRITE(0,*) "TO_RAD ERR: incorrect kind"
      WRITE(0,*) "Supported ones: deg, rad"
      CALL EXIT(1)
   END IF
   this%units = "rad"
   RETURN
END SUBROUTINE TO_RAD
!############################################################
!# SUBROUTINE: TO DEG #######################################
!############################################################
!> @brief
!! Manages angle units. Changes to degrees
!
!> @param[in,out] this - Angle subtype variable.
!
!> @author A.P. Muzas - alberto.muzas@uam.es
!> @date 06/Nov/2013
!> @version 1.0
!------------------------------------------------------------
SUBROUTINE TO_DEG(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(angle), INTENT(INOUT) :: this
   ! GO, GO, GO !!!------
   SELECT CASE(this%units)
      CASE("deg")
         RETURN
      CASE("rad")
         this%mag = this%mag*180.D0/pi
         this%units = "deg"
      CASE DEFAULT
         WRITE(0,*) "TO_DEG ERR: incorrect kind"
         WRITE(0,*) "Supported ones: rad, deg"
         CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE TO_DEG
!############################################################
!# SUBROUTINE: LENGTH_AU ####################################
!############################################################
!> @brief
!! Manages length units. Changes them to au
!
!> @param[in,out] this - Length subtype variable. Can be omitted.
!
!> @author A.P. Muzas - alberto.muzas@uam.es
!> @date 06/Nov/2013
!> @version 1.0
!------------------------------------------------------------
SUBROUTINE LENGTH_AU(this)
   IMPLICIT NONE
   ! I/O variables
   CLASS(length), INTENT(INOUT) :: this
   ! Run section
   SELECT CASE(this%units)
      CASE("angst")
         this%mag = this%mag / au2angst
         this%units = "au"
      CASE("au")
         ! do nothing
      CASE DEFAULT
         WRITE(0,*) "TO_ANGST ERR: incorrect units"
         WRITE(0,*) "Supported ones: angst, au"
         CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE LENGTH_AU
!###########################################################
!# SUBROUTINE: TO_ANGST_LENGTH
!###########################################################
!> @brief
!! simple units change function. From au to angstroem units
!-----------------------------------------------------------
SUBROUTINE TO_ANGST_LENGTH(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Length),INTENT(INOUT):: this
   ! Run section
   SELECT CASE(this%units)
      CASE('angst','Angst')
         ! do nothing
      CASE('au','bohr','bohrs','Bohr','Bohrs')
         this%mag = this%mag*au2angst
         this%units = "angst"
      CASE DEFAULT
         WRITE(0,*) "TO_ANGST ERR: incorrect units"
         WRITE(0,*) "Supported ones: angst, au"
         CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE TO_ANGST_LENGTH
!###########################################################
!# SUBROUTINE: TO_CELSIUS_UNITS
!###########################################################
!> @brief
!! simple units change function.
!-----------------------------------------------------------
subroutine TO_KELVIN_UNITS(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Temperature),intent(inout):: this
   ! Run section
   select case(this%units)
      case( 'Kelvin','kelvin','K','k','Kelvins','kelvins' )
         ! do nothing
      case( 'Celsius','celsius','C','c' )
         this%mag=this%mag+kelvinParam
      case( 'Fahrenheit','fahrenheit','F','f' )
         this%mag=(this%mag+fahrenheitParam)*5.d0/9.d0
      case default
         write(0,*) 'TO_KELVIN ERR: incorrect units'
         write(0,*) 'Supported ones: Kelvin, Celsius, Fahrenheit'
         call exit(1)
   end select
   this%units='Kelvin'
   return
end subroutine TO_KELVIN_UNITS
!############################################################
!# SUBROUTINE: MASS_AU ######################################
!############################################################
!> @brief
!! Manage mass units. Changes them to au
!
!> @param[in,out] this - Mass subtype variable. Can be omitted.
!
!> @author A.P. Muzas - alberto.muzas@uam.es
!> @date 06/Nov/2013
!> @version 1.0
!------------------------------------------------------------
subroutine mass_au(this)
   ! Initial declarations
   implicit none
   ! I/O variables
    class(Mass),intent(inout) :: this
    ! Run section
    select case( this%units )
    case( 'hmass','Hmass' )
       this%mag = this%mag*hmass2au
    case( 'dmass','Dmass' )
       this%mag = this%mag*dmass2au
    case( "pmass" )
       this%mag = this%mag*pmass2au
    case( 'Da','da','dalton','Dalton','Daltons','daltons' )
       this%mag = this%mag*dalton2au
    case( "au" )
       return
    case default
       write(0,*) "MASS_AU ERR: incorrect units"
       write(0,*) "Supported ones: hmass, dmass, pmass, au, Da"
       call exit(1)
    end select
    this%units = "au"
    return
end subroutine mass_au
!############################################################
!# SUBROUTINE: ENERGY_AU ####################################
!############################################################
!> @brief
!! Manages energy units. Changes them to au
!
!> @param[in,out] this - Energy subtype variable. Can be omitted.
!
!> @author A.P. Muzas - alberto.muzas@uam.es
!> @date 06/Nov/2013
!> @version 1.0
!------------------------------------------------------------
SUBROUTINE ENERGY_AU(this)
        IMPLICIT NONE
        ! I/O variables
        CLASS(energy), INTENT(INOUT) :: this
        ! GO, GO, GO !!!
        IF (this%units.EQ."ev") THEN
                this%mag = this%mag/au2ev
        ELSE IF (this%units.EQ."kcalmol") THEN
                this%mag = this%mag/au2kcalmol
        ELSE IF (this%units.EQ."kjmol") THEN
                this%mag = this%mag/au2kjmol
        ELSE IF (this%units.EQ."au") THEN
                RETURN
        ELSE
                WRITE(0,*) "ENERGY_AU ERR: incorrect units"
                WRITE(0,*) "Supported ones: ev, kcalmol, kjmol, au"
                CALL EXIT(1)
        END IF
        this%units = "au"
        RETURN
END SUBROUTINE ENERGY_AU
!############################################################
!# SUBROUTINE: TIME_AU ######################################
!############################################################
!> @brief
!! Manages time units. Changes them to au
!
!> @param[in,out] this - Time subtype variable. Can be omitted.
!
!> @author A.P. Muzas - alberto.muzas@uam.es
!> @date 06/Nov/2013
!> @version 1.0
!------------------------------------------------------------
SUBROUTINE TIME_AU(this)
   ! Initial declarations
        IMPLICIT NONE
        ! I/O variables
        CLASS(time), INTENT(INOUT) :: this
        ! GO, GO, GO !!!
        IF (this%units.EQ."ps") THEN
                this%mag = this%mag/au2ps
        ELSE IF (this%units.EQ."fs") THEN
                this%mag = this%mag/au2fs
        ELSE IF (this%units.EQ."au") THEN
                RETURN
        ELSE
                WRITE(0,*) "TIME_AU ERR: incorrect units"
                WRITE(0,*) "Supported ones: fs, ps, au"
                CALL EXIT(1)
        END IF
        this%units = "au"
        RETURN
END SUBROUTINE TIME_AU
END MODULE UNITS_MOD
!#########################################################
! MODULE SURFACE_MOD
!
!> @brief
!! Should contain everything related with periodic 2D surfaces
!##########################################################
MODULE SURFACE_MOD
use UNITS_MOD, only: Length, pi
use MATHS_MOD, only: INV_MTRX
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////////////////////////
! TYPE : Atom_list
!> @brief
!! Auxiliary type data that contains info about a list of atoms of the same type
!
!> @param n - Number of atoms in this list
!> @param alias - Its periodic table symbol
!> @param atom - Matrix of positions of each atom in this list
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Feb/2014; Jun/2014
!> @version 2.0
!-------------------------------------------------------------------------------------
TYPE,PRIVATE :: Atom_list
   INTEGER(KIND=4) :: n ! number of atoms in this list
   CHARACTER(LEN=2) :: alias ! atom name, periodic table
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: atom
END TYPE
!/////////////////////////////////////////////////////////////////////////////////////
! TYPE: Surface
!> @brief
!! Class that contains all data needed from a 2D periodic surface
!
!> @param alias - Human-fiendly name for this surface
!> @param filename - Input file that contains all surface information
!> @param symmlabel - Name for the symmetry group of this surface
!> @param units -  Units in which distances are stored
!> @param initialized - Controls if this class was initialized or not
!> @param norm_s1, norm_s2 - Norms of surface vectors s1 & s2
!> @param surf2cart_mtrx - Transformation matrix from surface coordinates to auxiliary cartesian
!> @param cart2surf_mtrx - Transformation matrix from auxiliary cartesian coordinates to surface
!> @param surfunit2cart_mtrx - Transformation matrix from unit surface corrdinates to auxiliary cartesians
!> @param cart2surfunit_mtrx - Transformation matrix from auxiliary cartesians to unit surface coordinates
!> @param recip2cart_mtrx - Transformation matrix from reciprocal lattice coordinates to auxiliary cartesians
!> @param cart2recip_mtrx - Transformation matrix from auxiliary cartesian coordinates to reciprocal lattice
!> @param metricsurf_mtrx - Metric matrix for surface coordinates
!> @param diff_atoms - Number of different types of atoms in the surface
!> @param atomtype - Array of lists of atoms
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 03/Feb/2014
!> @version 1.0
!
!> @see atom_list
!-------------------------------------------------------------------------------------
Type Surface
   character(len=30):: alias
   character(len=:),allocatable:: filename
   integer(kind=4):: order
   logical:: initialized=.false.
   real(kind=8):: angle
   real(kind=8),public,dimension(2) :: s1,s2
   real(kind=8),dimension(2,2):: surf2cart_mtrx
   real(kind=8),dimension(2,2):: cart2surf_mtrx
   real(kind=8),dimension(2,2):: surfunit2cart_mtrx
   real(kind=8),dimension(2,2):: cart2surfunit_mtrx
   real(kind=8),dimension(2,2):: recip2cart_mtrx
   real(kind=8),dimension(2,2):: cart2recip_mtrx
   integer(kind=4):: diff_atoms
   type(atom_list),dimension(:),allocatable,public:: atomType
   real(kind=8),dimension(2,2),public:: metricSurf_mtrx
   character(len=10),public:: units
   real(kind=8),public:: norm_s1,norm_s2
   character(len=4):: symmLabel
contains
   ! Initiallize
   PROCEDURE,PUBLIC:: INITIALIZE => INITIALIZE_SURFACE
   ! Operations block
   PROCEDURE,PUBLIC:: surf2cart => surf2cart_SURFACE
   PROCEDURE,PUBLIC:: cart2surf => cart2surf_SURFACE
   PROCEDURE,PUBLIC:: surfunit2cart => surfunit2cart_SURFACE
   PROCEDURE,PUBLIC:: cart2surfunit => cart2surfunit_SURFACE
   PROCEDURE,PUBLIC:: recip2cart => recip2cart_SURFACE
   PROCEDURE,PUBLIC:: cart2recip => cart2recip_SURFACE
   PROCEDURE,PUBLIC:: project_unitcell => project_unitcell_SURFACE
   PROCEDURE,PUBLIC:: project_iwscell => project_iwscell_SURFACE
   ! Get block
   PROCEDURE,PUBLIC:: getsymmlabel => getsymmlabel_SURFACE
   PROCEDURE,PUBLIC:: getfilename => getfilename_SURFACE
   ! Tools block
   PROCEDURE,PUBLIC:: PRINT_PATTERN => PRINT_PATTERN_SURFACE
   PROCEDURE,PUBLIC:: MOVE_PATTERN => MOVE_PATTERN_SURFACE
   ! Enquire block
   PROCEDURE,PUBLIC:: is_initialized => is_initialized_SURFACE
end type
! MODULE CONTAINS
contains
!###########################################################
!# SUBROUTINE: MOVE_PATTERN_SURFACE
!###########################################################
!> @brief
!! Moves surface pattern and projects it uppon the unit cell
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE MOVE_PATTERN_SURFACE(this,dr)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(SUrface),INTENT(INOUT):: this
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: dr
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   INTEGER(KIND=4) :: natoms
   ! Run section
   DO i = 1, this%diff_atoms
      DO j = 1, this%atomtype(i)%n
         this%atomtype(i)%atom(j,1:2)=this%atomtype(i)%atom(j,1:2)+dr
         this%atomtype(i)%atom(j,1:2)=this%project_unitcell(this%atomtype(i)%atom(j,1:2))
      END DO
   END DO
   RETURN
END SUBROUTINE MOVE_PATTERN_SURFACE
!###########################################################
!# SUBROUTINE: PRINT_PATTERN_SURFACE
!###########################################################
!> @brief
!! Prints the pattern defined in surface to a specific unit.
!! It can be defined, as well, the order of the pattern.
!
!> @param[in] this - Surface class object
!> @param[in] wunit - Unit to print output
!> @param[in] order - Order of the pattern.
!> @param[in] format_out - string: XYZ available
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE PRINT_PATTERN_SURFACE(this,wunit,order,format_out)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface),INTENT(IN):: this
   INTEGER(KIND=4),INTENT(IN) :: wunit
   INTEGER(KIND=4),INTENT(IN) :: order
   CHARACTER(LEN=3),INTENT(IN) :: format_out
   ! Local variables
   INTEGER(KIND=4) :: i,j,n,k ! counters
   REAL(KIND=8),DIMENSION(2) :: aux
   REAL(KIND=8),DIMENSION(2) :: xypos
   REAL(KIND=8),DIMENSION(3) :: P
   ! Run section
   SELECT CASE(format_out)
      CASE("XYZ")
         DO i = 1, this%diff_atoms
            SELECT CASE(order)
               CASE(0)
                  DO j = 1, this%atomtype(i)%n
                     WRITE(wunit,*) this%atomtype(i)%alias,this%atomtype(i)%atom(j,:)
                  END DO
               CASE(1 :)
                  DO j = 1, this%atomtype(i)%n
                     WRITE(wunit,*) this%atomtype(i)%alias,this%atomtype(i)%atom(j,:)
                     xypos(1:2)=this%atomtype(i)%atom(j,1:2)
                     P(3)=this%atomtype(i)%atom(j,3)
                     DO n = 1, order
                        DO k = -n, n
                           aux(1)=dfloat(n)
                           aux(2)=dfloat(k)
                           aux=this%surf2cart(aux)
                           P(1:2)=xypos(1:2)+aux(1:2)
                           WRITE(wunit,*) this%atomtype(i)%alias,P
                           aux(1)=dfloat(-n)
                           aux(2)=dfloat(k)
                           aux=this%surf2cart(aux)
                           P(1:2)=xypos(1:2)+aux(1:2)
                           WRITE(wunit,*) this%atomtype(i)%alias,P
                        END DO
                        DO k = -n+1, n-1
                           aux(1)=dfloat(k)
                           aux(2)=dfloat(n)
                           aux=this%surf2cart(aux)
                           P(1:2) =xypos(1:2)+aux(1:2)
                           WRITE(wunit,*) this%atomtype(i)%alias,P
                           aux(1)=dfloat(k)
                           aux(2)=dfloat(-n)
                           aux=this%surf2cart(aux)
                           P(1:2)=xypos(1:2)+aux(1:2)
                           WRITE(wunit,*) this%atomtype(i)%alias,P
                        END DO
                     END DO
                  END DO
            END SELECT
         END DO
      CASE DEFAULT
         WRITE(0,*) "PRINT_PATTERN_SURFACE ERR: Wrong format specifier"
         WRITE(0,*) "Implemented ones: XYZ"
         CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE PRINT_PATTERN_SURFACE
!###########################################################
!# FUNCTION: getfilename
!###########################################################
!> @brief
!! Typical enquire function
!-----------------------------------------------------------
PURE FUNCTION getfilename_SURFACE(this) result(filename)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface),INTENT(IN):: this
   CHARACTER(LEN=:),ALLOCATABLE:: filename
   ! Run section
   filename=this%filename
   RETURN
END FUNCTION getfilename_SURFACE
!###########################################################
!# FUNCTION: getsymmlabel
!###########################################################
!> @brief
!! typical enquire function
!-----------------------------------------------------------
PURE FUNCTION getsymmlabel_SURFACE(this) result(symmlabel)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface),INTENT(IN):: this
   CHARACTER(LEN=:),ALLOCATABLE:: symmlabel
   ! Run section
   symmlabel=this%symmlabel
   RETURN
END FUNCTION getsymmlabel_SURFACE
!###########################################################
!# FUNCTION: is_initialized
!###########################################################
! - Check if surface type is already initialized
!-----------------------------------------------------------
PURE FUNCTION is_initialized_SURFACE(surf) result(bool)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface),INTENT(IN):: surf
   LOGICAL:: bool
   ! Run section
   bool=surf%initialized
   RETURN
END FUNCTION is_initialized_SURFACE

!###############################################################################
!# SUBROUTINE: INITIALIZE ######################################################
!###############################################################################
!> @brief
!! Initializes surface from file @b filename
!-------------------------------------------------------------------------------
SUBROUTINE INITIALIZE_SURFACE(surf,filename)
   IMPLICIT NONE
   ! I/O Variables -----------------------------------------------
   CLASS(Surface),INTENT(INOUT):: surf
   CHARACTER(LEN=*),INTENT(IN):: filename
   ! Local variables ---------------------------------------------
   INTEGER:: i,j ! Counters
   INTEGER:: control
   CHARACTER(LEN=*),PARAMETER:: routinename = "INITIALIZE_SURFACE: "
   TYPE(Length):: len
   REAL(KIND=8),DIMENSION(2,2):: aux_r
   ! Run section --------------------------------------------
   surf%filename=filename
   IF (.NOT.surf%is_initialized()) THEN
      surf%initialized=.FALSE.
      OPEN(10,FILE=surf%filename,STATUS="old")
      READ(10,*) ! dummy line
      READ(10,*) surf%alias
      READ(10,*) surf%units
      ! Read surface basis vectors in cartesian coordinates.
      ! They should be written horizontally
      ! Store results in a.u.
      ! Set surface coordinates to cartesian change matrix and metadata
      DO i=1,2
         READ(10,*)  aux_r(i,:)
         DO j = 1, 2
            CALL len%READ(aux_r(i,j),surf%units)
            CALL len%TO_STD()
            aux_r(i,j)=len%getvalue()
         END DO
      END DO
      surf%s1=aux_r(1,:)
      surf%s2=aux_r(2,:)
      surf%norm_s1=norm2( surf%s1 )
      surf%norm_s2=norm2( surf%s2 )
      surf%angle=dacos(dot_product(surf%s1,surf%s2)/(surf%norm_s1*surf%norm_s2))
      READ(10,*) surf%symmlabel
      select case( trim(surf%symmlabel) )
      case('p4mm')
        surf%order=4
      case('p3')
        surf%order=6
      case default
         write(0,*) routinename//'ERR surface not implemented'
         write(0,*) 'Implemented ones: p4mm'
         call exit(1)
      end select
      aux_r=transpose(aux_r)
      surf%surf2cart_mtrx=aux_r
      ! Read number of different atom types in the surface
      READ(10,*) surf%diff_atoms
      ! Set the the atomic basis, i.e., coordinates for non equivalent  atoms inside
      ! the unit cell
      ALLOCATE(surf%atomtype(1:surf%diff_atoms))
      DO i=1, surf%diff_atoms
         READ(10,*) control, surf%atomtype(i)%alias, surf%atomtype(i)%n
         IF (control.NE.i) THEN
            WRITE(0,*) "INITIALIZE_SURF ERR: Atom type definitions are not in the correct order"
            CALL EXIT(1)
         END IF
         ALLOCATE(surf%atomtype(i)%atom(surf%atomtype(i)%n,3))
      END DO
      DO i=1, surf%diff_atoms
            DO j=1, surf%atomtype(i)%n
               READ(10,*) control,surf%atomtype(i)%atom(j,1),surf%atomtype(i)%atom(j,2),surf%atomtype(i)%atom(j,3)
               CALL len%READ(surf%atomtype(i)%atom(j,1),surf%units)
               CALL len%TO_STD()
               surf%atomtype(i)%atom(j,1)=len%getvalue()
               CALL len%READ(surf%atomtype(i)%atom(j,2),surf%units)
               CALL len%TO_STD()
               surf%atomtype(i)%atom(j,2)=len%getvalue()
               CALL len%READ(surf%atomtype(i)%atom(j,3),surf%units)
               CALL len%TO_STD()
               surf%atomtype(i)%atom(j,3)=len%getvalue()
               IF (control.NE.i) THEN
                  WRITE(0,*) "INITIALIZE_SURF ERR: Atom type definitions are not in the correct order"
                  CALL EXIT(1)
               END IF
            END DO
      END DO
      surf%units="au"
      CLOSE(10)
      ! Set metric matrix associated with surface coordinates.
      surf%metricsurf_mtrx=MATMUL(TRANSPOSE(surf%surf2cart_mtrx),surf%surf2cart_mtrx)
      ! Set Matrix: from auxiliar cartesian coordinates to surface coordinates
      CALL INV_MTRX(2,surf%surf2cart_mtrx,surf%cart2surf_mtrx)
      ! Set Matrix: from normalized surface coordinates to auxiliar cartesian coordinates
      FORALL(i=1:2) surf%surfunit2cart_mtrx(i,1)=surf%surf2cart_mtrx(i,1)/surf%norm_s1
      FORALL(i=1:2) surf%surfunit2cart_mtrx(i,2)=surf%surf2cart_mtrx(i,2)/surf%norm_s2
      ! Set Matrix: from auxiliar cartesian coordinates to normalized surface coordinates
      CALL INV_MTRX(2,surf%surfunit2cart_mtrx,surf%cart2surfunit_mtrx)
      ! Set Matrix: form reciprocal space to auxiliar cartesian coordinates
      surf%recip2cart_mtrx=TRANSPOSE(surf%cart2surf_mtrx)
      !FORALL(i=1:2,j=1:2) surf%recip_mtrx(i,j)=2*PI*surf%recip_mtrx(i,j)
      surf%recip2cart_mtrx=2.D0*PI*surf%recip2cart_mtrx
      ! Set Matrix: from auxiliar cartesian coordinates to reciprocal space
      CALL INV_MTRX(2,surf%recip2cart_mtrx,surf%cart2recip_mtrx)
      surf%initialized=.TRUE.
   ELSE
      WRITE(0,*) "INITIALIZE_SURF ERR: surface already initialized"
      CALL EXIT(1)
   END IF
   RETURN
END SUBROUTINE INITIALIZE_SURFACE
!###########################################################
!# FUNCTION: cart2surf
!###########################################################
! - Goes from auxiliar cartesian coordinates to surface coordinates
!-----------------------------------------------------------
FUNCTION cart2surf_SURFACE(surf,r)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface), INTENT(IN) :: surf
   REAL(KIND=8), DIMENSION(2), INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: cart2surf_SURFACE
   ! Run section
   cart2surf_SURFACE=matmul(surf%cart2surf_mtrx,r)
   RETURN
END FUNCTION cart2surf_SURFACE
!###########################################################
!# FUNCTION: surf2cart
!###########################################################
! - Goes from surface to auxiliar cartesian coordinates
!-----------------------------------------------------------
FUNCTION surf2cart_SURFACE(surf,r)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface),INTENT(IN) :: surf
   REAL(KIND=8), DIMENSION(2),INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: surf2cart_SURFACE
   ! Run section
   surf2cart_SURFACE=matmul(surf%surf2cart_mtrx,r)
   RETURN
END FUNCTION surf2cart_SURFACE
!###########################################################
!# FUNCTION: surfunit2cart
!###########################################################
! - Goes from unit surface to auxiliar cartesian coordinates
!-----------------------------------------------------------
FUNCTION surfunit2cart_SURFACE(surf,r)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface),INTENT(IN) :: surf
   REAL(KIND=8), DIMENSION(2),INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: surfunit2cart_SURFACE
   ! Run section
   surfunit2cart_SURFACE=matmul(surf%surfunit2cart_mtrx,r)
   RETURN
END FUNCTION surfunit2cart_SURFACE
!###########################################################
!# FUNCTION: cart2surfunit
!###########################################################
! - Goes from auxiliar cartesian coordinates to surface coordinates
!-----------------------------------------------------------
FUNCTION cart2surfunit_SURFACE(surf,r)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface), INTENT(IN) :: surf
   REAL(KIND=8), DIMENSION(2), INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: cart2surfunit_SURFACE
   ! Run section
   cart2surfunit_SURFACE=matmul(surf%cart2surfunit_mtrx,r)
   RETURN
END FUNCTION cart2surfunit_SURFACE
!###########################################################
!# FUNCTION: cart2recip
!###########################################################
! - Goes from auxiliar cartesian to reciprocal space coordinates
!-----------------------------------------------------------
FUNCTION cart2recip_SURFACE(surf,r)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface),INTENT(IN) :: surf
   REAL(KIND=8), DIMENSION(2),INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: cart2recip_SURFACE
   ! Run section
   cart2recip_SURFACE=matmul(surf%cart2recip_mtrx,r)
   RETURN
END FUNCTION cart2recip_SURFACE
!###########################################################
!# FUNCTION: recip2cart
!###########################################################
! - Goes from auxiliar cartesian to surface coordinates
!-----------------------------------------------------------
FUNCTION recip2cart_SURFACE(surf,r)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface), INTENT(IN) :: surf
   REAL(KIND=8), DIMENSION(2), INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: recip2cart_SURFACE
   ! Run section
   recip2cart_SURFACE=matmul(surf%recip2cart_mtrx,r)
   RETURN
END FUNCTION recip2cart_SURFACE
!################################################################
!# SUBROUTINE: PROJECT_UNITCELL #################################
!################################################################
!> @brief
!! Projects 2D point into the C4v unit cell.
!
!> @warning
!! - Input/output in cartesian coordinates (r)
!! - Special care should be taken to project into the correct quadrant
!!   of the cell.
!----------------------------------------------------------------
function project_unitcell_SURFACE(surf,r) result(newR)
   implicit none
   ! I/O variables
   class(Surface),intent(in):: surf
   real(kind=8),dimension(2),intent(in):: r
   ! Dummy function output variable
   real(kind=8),dimension(2):: newR
   ! Local variables
   real(kind=8),dimension(2):: center
   ! Parameters
   character(len=*),parameter:: routinename='project_unitcell_SURFACE: '
   ! HEY, HO! LET'S GO !!! ----------------------
   newR=surf%cart2surf( r )
   if( surf%symmLabel=='p4mm' ) then
      if( newR(1)>=0.d0 ) then
         center(1)=dfloat( int( newR(1) ) )
      else
         center(1)=dfloat( int( newR(1) ) )-1.d0
      endif
      if( newR(2)>=0.d0 ) then
         center(2)=dfloat( int( newR(2) ) )
      else
         center(2)=dfloat( int( newR(2) ) )-1.d0
      endif
   else
      write(0,*) 'ERR '//routinename//'wallpaper symmetry not implemented'
      write(0,*) 'Implemented ones: p4mm'
      call exit(1)
   endif
   newR(:)=newR(:)-center(:)
   newR(:)=surf%surf2cart( newR(:) )
   return
end function project_unitcell_SURFACE
!################################################################
! SUBROUTINE: project_iwscell ###################################
!################################################################
!> @brief
!! Projects R into Irreducible WS cell. Needs information from
!! surface main vectors.
!
!> @param[in] surf - Surface specifications
!> @param[in] x - 2D point
!----------------------------------------------------------------
function project_iwscell_SURFACE(surf,x) result(r)
   implicit none
   ! I/O variables
   class(Surface),intent(in):: surf
   real(kind=8),dimension(:),intent(in):: x
   ! dummy variable
   real(kind=8),dimension(:),allocatable:: r
   ! Local variables
   real(kind=8):: auxReal
   ! Parameters
   character(len=*),parameter:: routineName='PROJECT_IWCELL_SURFACE: '
   ! HEY, HO! LET'S GO! ------------------
   ! Go to surface coordinates
   allocate(r(size(x)),source=x)
   r(1:2)=surf%project_unitCell( r(1:2) )
   r(1:2)=surf%cart2surf( r(1:2) )
   ! ----------------------------------------------------------
   if( surf%symmLabel=='p4mm' ) then
      if( r(1)>0.5d0 .and. r(2)<=1.d0-r(1) ) then ! sector II
         r(1)=1.d0-r(1)
         if( size(r)==6 ) r(6)=pi-r(6)
      elseif( r(1)>0.5d0 .and. r(2)<0.5d0 ) then ! sector III
         auxReal=r(1)
         r(1)=r(2)
         r(2)=1.d0-auxReal
         if( size(r)==6 ) r(6)=r(6)-pi/2.d0
      elseif( r(1)>0.5d0 .and. r(2)<r(1) ) then ! sector IV
         auxReal=r(1)
         r(1)=1.d0-r(2)
         r(2)=1.d0-auxReal
         if( size(r)==6 ) r(6)=3.d0*pi/2.d0-r(6)
      elseif( r(1)>0.5d0 ) then ! sector V
         r(1)=1.d0-r(1)
         r(2)=1.d0-r(2)
         if( size(r)==6 ) r(6)=r(6)-pi
      elseif( r(1)<=0.5d0 .and. r(2)>1.d0-r(1) ) then ! sector VI
         r(2)=1.d0-r(2)
         if( size(r)==6 ) r(6)=-r(6)
      elseif( r(1)<=0.5d0 .and. r(2)>0.5d0 ) then ! sector VII
         auxReal=r(1)
         r(1)=1.d0-r(2)
         r(2)=auxReal
         if( size(r)==6 ) r(6)=r(6)-3.d0*pi/2.d0
      elseif( r(1)<=0.5d0 .and. r(2)>r(1) ) then ! sector VIII
         auxReal=r(1)
         r(1)=r(2)
         r(2)=auxReal
         if( size(r)==6 ) r(6)=pi/2.d0-r(6)
      elseif( r(1)<=0.5d0 .and. r(2)<=r(1) ) then ! sector I
         ! do nothing
      else
         write(0,*) routinename//'sector selector had a weird problem. Rewrite this switch'
         call exit(1)
      endif
   ! --------------------------------------------------------
   else ! default case
      write(0,*) routinename//'surface is not implemented'
      call exit(1)
   endif
   ! Go to cartesian coordinates
   r(1:2)=surf%surf2cart( r(1:2) )
   return
end function project_iwscell_SURFACE
END MODULE SURFACE_MOD
!#########################################################
! MODULE SURFACE_MOD
!> @brief
!! Specific implementation of module SURFACE_MOD
!##########################################################
module LiF001SURF_MOD
use SURFACE_MOD
implicit none
!/////////////////////////////////////////////////////////////////////////////////////
! TYPE: LiF001SURF
!> @brief
!! Specific implementation for LiF001 extending Surface Type
!-------------------------------------------------------------------------------------
type,extends(Surface):: LiF001Surf
	contains
	   procedure,public:: initialize => initialize_LiF001
end type LiF001Surf

contains
!###############################################################################
!# SUBROUTINE: INITIALIZE_LiF001 ###############################################
!###############################################################################
!> @brief
!! Specific implementation of initialize
!-------------------------------------------------------------------------------
subroutine initialize_LiF001(surf,filename)
   implicit none
   ! I/O Variables -----------------------------------------------
   class(LiF001Surf),intent(inout):: surf
   character(len=*),intent(in):: filename
   ! Run section --------------------------------------------
   surf%alias='LiF001'
   surf%diff_atoms=2
   surf%symmLabel='p4mm'
   surf%s1=[5.4433561257770959d0,0.d0]
   surf%s2=[0.d0,5.4433561257770959d0]
   allocate( surf%atomType(2) )
   surf%units='au'
   surf%norm_s1=5.4433561257770959d0
   surf%norm_s2=5.4433561257770959d0
   ! Set basis vectors
   surf%surf2cart_mtrx(:,1)=surf%s1(:)
   surf%surf2cart_mtrx(:,2)=surf%s2(:)
   ! Set Lithium specifications
   surf%atomtype(1)%alias='Li'
   surf%atomtype(1)%n=1
   allocate(surf%atomtype(1)%atom(1,3))
   surf%atomtype(1)%atom(1,:)=[0.d0,0.d0,-0.123593758269d0]
   ! Set Fluorinne specifications
   surf%atomtype(2)%alias='F'
   surf%atomtype(2)%n=1
   allocate(surf%atomtype(2)%atom(1,3))
   surf%atomtype(2)%atom(1,:)=[2.72167806289d0,2.72167806289d0,0.d0]
   ! Set conversion matrices
   surf%metricsurf_mtrx=matmul(transpose(surf%surf2cart_mtrx),surf%surf2cart_mtrx)
   call inv_mtrx(2,surf%surf2cart_mtrx,surf%cart2surf_mtrx)
   surf%surfunit2cart_mtrx(:,1)=surf%surf2cart_mtrx(:,1)/surf%norm_s1
   surf%surfunit2cart_mtrx(:,2)=surf%surf2cart_mtrx(:,2)/surf%norm_s2
   call inv_mtrx(2,surf%surfunit2cart_mtrx,surf%cart2surfunit_mtrx)
   surf%recip2cart_mtrx=transpose(surf%cart2surf_mtrx)
   surf%recip2cart_mtrx=2.D0*pi*surf%recip2cart_mtrx
   call inv_mtrx(2,surf%recip2cart_mtrx,surf%cart2recip_mtrx)
   surf%initialized=.true.
   return
end subroutine initialize_LiF001
end module LiF001SURF_MOD

MODULE PES_MOD
IMPLICIT NONE
!////////////////////////////////////////////////////////////////////////////////////////
! TYPE: PES
!> @brief
!! Type for a generic Potential energy surface
!
!> @param alias - An alias for the PES
!> @param dimensions - Number of coordinates which this PES depends on
!> @param r - last point where the potential was calculated
!> @param dvdu - last calculated derivatives
!> @param v - last value of the potential (at r, obviously)
!> @param initialized - Controls status of the PES
!> @param atomdat - Contains array of atoms (with some info) related to this PES
!
!> @details
!! In order to add a new PES type, one should create a fortran module that declares an
!! extension of this abstrac type inside a module, using only the tools defined here. Take care of the deferred
!! type-bounded subtroutines, because an specific implementation should be given in this new file.
!! Later, load the new module in link_PES.f90 file. Only generic PES should be included in the repository.
!------------------------------------------------------------------------------------------
TYPE,ABSTRACT :: PES
PRIVATE
   CHARACTER(LEN=:),ALLOCATABLE:: alias
   CHARACTER(LEN=:),ALLOCATABLE:: pestype
   INTEGER:: dimensions
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: r
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: dvdu
   REAL(KIND=8):: v
   LOGICAL:: initialized = .FALSE.
   CONTAINS
      ! Initialization block
      PROCEDURE(INITIALIZE_PES),DEFERRED,PUBLIC:: INITIALIZE     ! DEFERRED, see INTERFACE !!!!!!!!!!!
      ! Set block
      PROCEDURE,NON_OVERRIDABLE,PUBLIC:: SET_ALIAS => SET_ALIAS_PES
      PROCEDURE,NON_OVERRIDABLE,PUBLIC:: SET_PESTYPE => SET_PESTYPE_PES
      PROCEDURE,NON_OVERRIDABLE,PUBLIC:: SET_DIMENSIONS => SET_DIMENSIONS_PES
      procedure,non_overridable,public:: set_lastPot => set_lastPot_PES
      procedure,non_overridable,public:: set_lastDeriv => set_lastDeriv_PES
      ! Get block
      PROCEDURE,NON_OVERRIDABLE,PUBLIC:: getAlias => getalias_PES
      PROCEDURE,NON_OVERRIDABLE,PUBLIC:: getDimensions => getdimensions_PES
      PROCEDURE,NON_OVERRIDABLE,PUBLIC:: getLastPotValue => getLastPotValue_PES
      procedure,non_overridable,public:: getLastDerivValue => getLastDerivValue_PES
      PROCEDURE(GET_V_AND_DERIVS_PES),DEFERRED:: GET_V_AND_DERIVS ! DEFERRED, see INTERFACE !!!!!!!!!!
      ! Enquire block
      PROCEDURE,NON_OVERRIDABLE,PUBLIC:: is_initialized => is_initialized_PES
      PROCEDURE(is_allowed_PES),DEFERRED,PUBLIC:: is_allowed      ! DEFERRED, see INTERFACE !!!!!!!!!
END TYPE PES
ABSTRACT INTERFACE
   !###########################################################
   !# SUBROUTINE: INITIALIZE_PES
   !###########################################################
   !-----------------------------------------------------------
   SUBROUTINE INITIALIZE_PES(this,filename,tablename)
      IMPORT PES
      CLASS(PES),INTENT(OUT):: this
      CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: filename,tablename
   END SUBROUTINE INITIALIZE_PES
   !###########################################################
   !# SUBROUTINE: GET_V_AND_DERIVS_PES
   !###########################################################
   !-----------------------------------------------------------
   SUBROUTINE GET_V_AND_DERIVS_PES(this,x,v,dvdu,errCode)
      IMPORT PES
      CLASS(PES),TARGET,INTENT(in):: this
      REAL(KIND=8),DIMENSION(:),INTENT(IN):: x
      REAL(KIND=8),INTENT(OUT):: v
      REAL(KIND=8),DIMENSION(:),INTENT(OUT):: dvdu
      integer(kind=1),intent(out),optional:: errCode
   END SUBROUTINE GET_V_AND_DERIVS_PES
   !###########################################################
   !# FUNCTION: is_allowed_PES
   !###########################################################
   !> @brief
   !! Enquires if the potential can be calculated at @b X
   !-----------------------------------------------------------
   LOGICAL FUNCTION is_allowed_PES(this,x)
      IMPORT PES
      CLASS(PES),INTENT(IN):: this
      REAL(KIND=8),DIMENSION(:),INTENT(IN):: x
   END FUNCTION is_allowed_PES
END INTERFACE

private set_lastPot_PES, set_lastDeriv_PES, set_pesType_PES, set_alias_PES,&
        set_dimensions_PES, is_initialized_PES, getAlias_PES, getDimensions_PES,&
        getLastPotValue_PES, getLastDerivValue_PES

!////////////////////////////////////////////////////////////////////////////////////////
CONTAINS
!###############################################################
! SUBROUTINE: set_lastPot ######################################
!###############################################################
!> @brief
!! Standard set subroutine. Sets v atribute
!---------------------------------------------------------------
subroutine set_lastPot_PES(this,lastPot)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Pes),intent(inout):: this
   real(kind=8),intent(in):: lastPot
   ! Run section
   this%v=lastPot
   return
end subroutine set_lastPot_PES
!###############################################################
! SUBROUTINE: set_lastDeriv ####################################
!###############################################################
!> @brief
!! Standard set subroutine. Sets dvdu atribute
!---------------------------------------------------------------
subroutine set_lastDeriv_PES(this,lastDeriv)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Pes),intent(inout):: this
   real(kind=8),dimension(:),intent(in):: lastDeriv
   ! Parameters
   character(len=*),parameter:: routinename='SET_LASTDERIV_PES: '
   ! Run section
   select case( size(lastDeriv) /= size(this%dvdu) )
   case(.true.)
      this%dvdu(:)=lastDeriv(:)
   case(.false.)
      write(0,*) 'ERR: '//routinename//'mismatch between derivatives'
      call exit(1)
   end select
   return
end subroutine set_lastDeriv_PES
!###############################################################
! SUBROUTINE: SET_ALIAS ########################################
!###############################################################
!> @brief
!! Standard set subroutine. Sets alias atribute
!> @details
!! - If no argument is given a default string will be loaded
!---------------------------------------------------------------
SUBROUTINE SET_PESTYPE_PES(this,pestype)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES),INTENT(INOUT):: this
   CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: pestype
   ! Run section
   SELECT CASE(present(pestype))
      CASE(.true.)
         ALLOCATE(this%pestype,source=trim(pestype))
      CASE(.false.)
         this%pestype="NoPesType"
   END SELECT
   RETURN
END SUBROUTINE SET_PESTYPE_PES
!###############################################################
! SUBROUTINE: SET_ALIAS ########################################
!###############################################################
!> @brief
!! Sets an alias for PES.
!> @details
!! - If no argument is given a default string will be loaded
!---------------------------------------------------------------
SUBROUTINE SET_ALIAS_PES(this,alias)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES),INTENT(INOUT) :: this
   CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: alias
   ! Run section
   IF (PRESENT(alias)) THEN
      this%alias=trim(alias)
   ELSE
   this%alias="NoAlias"
   END IF
   RETURN
END SUBROUTINE SET_ALIAS_PES
!###############################################################
! SUBROUTINE: SET_DIMENSIONS ###################################
!###############################################################
!> @brief
!! Sets number of dimensions and allocates some arrays
!> @details
!! - If no argument is given, default value will be 1
!---------------------------------------------------------------
SUBROUTINE SET_DIMENSIONS_PES(this,dimensions)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES), INTENT(INOUT) :: this
   INTEGER,INTENT(IN),OPTIONAL :: dimensions
   ! Run section
   IF (PRESENT(dimensions)) THEN
      this%dimensions=dimensions
   ELSE
      this%dimensions=1
   END IF
   ! Allocate arrays:
   ALLOCATE(this%r(1:this%dimensions))
   ALLOCATE(this%dvdu(1:this%dimensions))
   RETURN
END SUBROUTINE SET_DIMENSIONS_PES
!#########################################################
! FUNCTION: is_initiallized ##############################
!#########################################################
!> @brief
!! Enquires whether the PES is initiallized or not
!---------------------------------------------------------
LOGICAL FUNCTION is_initialized_PES(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES), INTENT(IN) :: this
   ! Run section
   IF(this%initialized) THEN
      is_initialized_PES=.TRUE.
   ELSE
      is_initialized_PES=.FALSE.
   END IF
   RETURN
END FUNCTION is_initialized_PES
!##########################################################
! FUNCTION: get_alias #####################################
!##########################################################
!> @brief
!! Common get function. Gets alias from PES derived type
!----------------------------------------------------------
CHARACTER(LEN=30) FUNCTION getAlias_PES(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O VAriables
   CLASS(PES), INTENT(IN) :: this
   ! Run section
   getalias_PES=this%alias
   RETURN
END FUNCTION getAlias_PES
!##########################################################
! FUNCTION: get_dimensions ################################
!##########################################################
!> @brief
!! Common get function. Gets dimensions from PES derived type
!----------------------------------------------------------
pure function getDimensions_PES(this) result(dimension)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(PES),intent(in) :: this
   ! Dummy function variable
   integer(kind=4):: dimension
   ! Run section
   dimension=this%dimensions
   return
end function getDimensions_PES
!##########################################################
! FUNCTION: getLastPotValue ###############################
!##########################################################
!> @brief
!! Common get function. Gets last value calculated of the PES
!----------------------------------------------------------
function getLastPotValue_PES(this) result(lastPot)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES), INTENT(IN) :: this
   ! Dummy funciton variable
   real(kind=8):: lastPot
   ! Run part
   lastPot=this%v
end function getLastPotValue_PES
!##########################################################
! FUNCTION: getLastDerivValue #############################
!##########################################################
!> @brief
!! Common get function. Gets last derivatives calculated of the PES
!----------------------------------------------------------
function getLastDerivValue_PES(this) result(lastDeriv)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES),INTENT(IN):: this
   ! Dummy funciton variable
   real(kind=8),dimension(:),allocatable:: lastDeriv
   ! Run part
   allocate(lastDeriv(this%dimensions))
   lastDeriv=this%dvdu
end function getLastDerivValue_PES
END MODULE PES_MOD

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
TYPE,ABSTRACT :: Interpol1d
   INTEGER(KIND=4) :: n
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: x
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: f
   CONTAINS
      ! Initialization block
      PROCEDURE,NON_OVERRIDABLE,PUBLIC :: READ => READ_INTERPOL1D
      ! Get functions block
      PROCEDURE(getvalue_interpol1d),PUBLIC,DEFERRED :: getvalue ! child types, override this
      PROCEDURE(getvalue_interpol1d),PUBLIC,DEFERRED :: getderiv ! child types, override this
      ! Plot tools block
      PROCEDURE,NON_OVERRIDABLE,PUBLIC :: PLOTDATA => PLOT_DATA_INTERPOL1D
      PROCEDURE,NON_OVERRIDABLE,PUBLIC :: PLOT => PLOT_INTERPOL_INTERPOL1D
END TYPE Interpol1d
!//////////////////////////////////////////////////////////////////////
ABSTRACT INTERFACE
   !###########################################################
   !# FUNCTION: getvalue_interpol1d
   !###########################################################
   !> @brief
   !! Dummy function. Override it!!
   !-----------------------------------------------------------
   REAL(KIND=8) FUNCTION getvalue_interpol1d(this,x,shift)
      IMPORT Interpol1d
      CLASS(Interpol1d),TARGET,INTENT(IN) :: this
      REAL(KIND=8),INTENT(IN) :: x
      REAL(KIND=8),OPTIONAL,INTENT(IN) :: shift
   END FUNCTION getvalue_interpol1d
END INTERFACE
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
SUBROUTINE PLOT_INTERPOL_INTERPOL1D(this,npoints,filename,shift)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   INTEGER,INTENT(IN) :: npoints
   CLASS(Interpol1d),INTENT(IN),TARGET :: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),OPTIONAL,INTENT(IN) :: shift
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
   SELECT CASE(present(shift))
      CASE(.TRUE.)
         ! body
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   WRITE(11,*) xmin,this%getvalue(xmin),this%getderiv(xmin)
   DO i=1, inpoints
      x=xmin+(dfloat(i)*delta)
      WRITE(11,*) x,this%getvalue(x),this%getderiv(x)
   END DO
   WRITE(11,*) xmax,this%getvalue(xmax),this%getderiv(xmax)
   CLOSE(11)
END SUBROUTINE PLOT_INTERPOL_INTERPOL1D
END MODULE INTERPOL1D_MOD

!###################################################################################################
! MODULE: CUBICSPLINES_MOD
!
!> @brief
!! Module that manages cubic splines interpolations
!
!##################################################################################################
MODULE CUBICSPLINES_MOD
! Initial declarations
use INTERPOL1D_MOD
use MATHS_MOD
IMPLICIT NONE
!//////////////////////////////////////////////////////////////////////
! TYPE: Cubic Splines
!> @brief
!! Extends type Interpol1d. It stores all information and actions that can be extracted
!! from a cubic splines interpolation
!
!> @param dv2vdz - Set of @f$\frac{d^{2}F(x_{i})}{dx^{2}}@f$, calculated by DSPLIN
!> @param coeff - Matrix of cubic splines coefficients. Coeff(i,j) stands for the
!!                j'th coefficient of the i'th cubic spline
!> @see dsplin
!------------------------------------------------------------------------
TYPE,EXTENDS(Interpol1d) :: Csplines
   PRIVATE
   LOGICAL :: is_initialized = .FALSE.
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: d2fdx
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: coeff
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE,PUBLIC :: xmin
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE,PUBLIC :: xroot
   CONTAINS
      ! Public procedures
      PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_CUBIC_SPLINES
      PROCEDURE,PUBLIC :: REINTERPOL => REINTERPOL_CSPLINES
      PROCEDURE,PUBLIC :: getvalue => get_csplines_value
      PROCEDURE,PUBLIC :: getderiv => get_csplines_dfdx_value
      PROCEDURE,PUBLIC :: getderiv2 => get_csplines_df2dx2_value
      PROCEDURE,PUBLIC :: GET_V_AND_DERIVS => GET_V_AND_DERIVS_CSPLINES
      PROCEDURE,PUBLIC :: SET_MINIMUM => SET_MINIMUM_CSPLINES
      PROCEDURE,PUBLIC :: SET_XROOT => SET_XROOT_CSPLINES
      ! Private procedures
      PROCEDURE,PRIVATE :: SET_SECOND_DERIVS => DSPLIN
      PROCEDURE,PRIVATE :: SET_COEFF => SET_CUBIC_SPLINES_COEFF
END TYPE Csplines
!//////////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: REINTERPOL_CSPLINES
!###########################################################
!> @brief
!! Interpolates an already interpolated CSplines job
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!
!> @see interpol_cubic_splines
!-----------------------------------------------------------
SUBROUTINE REINTERPOL_CSPLINES(this,dz1,id1,dz2,id2)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Csplines),INTENT(INOUT) :: this
   REAL(KIND=8),INTENT(IN) :: dz1,dz2
   INTEGER(KIND=4),INTENT(IN) :: id1,id2
   ! Local variables
   ! Run section
   SELECT CASE(this%is_initialized)
      CASE(.FALSE.)
         WRITE(0,*) "REINTERPOL_CSPLINES ERR: use INTERPOL intead of REINTERPOL, this csplines &
                     job was not prior initialized"
         CALL EXIT(1)
      CASE(.TRUE.)
         ! do nothing
   END SELECT
   DEALLOCATE(this%d2fdx)
   DEALLOCATE(this%coeff)
   this%is_initialized=.FALSE.
   CALL this%INTERPOL(dz1,id1,dz2,id2)
   SELECT CASE(ALLOCATED(this%xmin))
      CASE(.TRUE.)
         DEALLOCATE(this%xmin)
         CALL this%SET_MINIMUM()
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE REINTERPOL_CSPLINES
!###########################################################
!# SUBROUTINE: SET_MINIMUM_CSPLINES
!###########################################################
!> @brief
!! Locates minimums in the interpolation. Analytical solution
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 26/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_MINIMUM_CSPLINES(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Csplines),TARGET,INTENT(INOUT)::this
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   REAL(KIND=8),DIMENSION(:),POINTER :: coeff
   REAL(KIND=8),POINTER :: x1,x2
   REAL(KIND=8) :: s
   REAL(KIND=8) :: discriminant
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: aux
   CHARACTER(LEN=22),PARAMETER :: routinename="SET_MINIMUM_CSPLINES: "
   ! Run section
   DO i = 1, this%n-1 ! number of splines
      coeff => this%coeff(i,:)
      x1 => this%x(i)
      x2 => this%x(i+1)
      discriminant=coeff(2)**2.D0-3.D0*coeff(1)*coeff(3)
      SELECT CASE( discriminant > 0.D0)
         CASE(.TRUE.)
            s=(-coeff(2)+dsqrt(discriminant))/(3.D0*coeff(1))
         CASE(.FALSE.)
            CYCLE
      END SELECT
      SELECT CASE(s >= 0.D0 .AND. s <= x2-x1)
         CASE(.TRUE.)
            SELECT CASE(ALLOCATED(this%xmin))
               CASE(.TRUE.)
                  ALLOCATE(aux(size(this%xmin)))
                  aux=this%xmin
                  DEALLOCATE(this%xmin)
                  ALLOCATE(this%xmin(size(aux)+1))
                  this%xmin(size(aux)+1)=s+x1
                  DEALLOCATE(aux)
               CASE(.FALSE.)
                  ALLOCATE(this%xmin(1))
                  this%xmin=s+x1
            END SELECT
         CASE(.FALSE.)
            CYCLE
      END SELECT
   END DO
   RETURN
END SUBROUTINE SET_MINIMUM_CSPLINES
!###########################################################
!# SUBROUTINE: SET_XROOT_CSPLINES
!###########################################################
!> @brief
!! Locates roots in the interpolation. Analytical solution
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 26/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_XROOT_CSPLINES(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Csplines),TARGET,INTENT(INOUT)::this
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   REAL(KIND=8),DIMENSION(:),POINTER :: coeff
   REAL(KIND=8),POINTER :: x1,x2
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: s
   REAL(KIND=8) :: discriminant
   REAL(KIND=8) :: a1,a2,a3,theta
   REAL(KIND=8):: param_q, param_r
   REAL(KIND=8) :: param_s, param_t
   COMPLEX(KIND=8):: param_sc, param_tc
   INTEGER(KIND=4) :: nold
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: aux
   CHARACTER(LEN=22),PARAMETER :: routinename="SET_MINIMUM_CSPLINES: "
   ! Run section
   DO i = 1, this%n-1 ! number of splines
      coeff => this%coeff(i,:)
      x1 => this%x(i)
      x2 => this%x(i+1)
      a1=coeff(2)/coeff(1)
      a2=coeff(3)/coeff(1)
      a3=coeff(4)/coeff(1)
      param_q=(3.D0*a2-a1**2.D0)/9.D0
      param_r=(9.D0*a1*a2-27.D0*a3-2.D0*a1**3.D0)/54.D0
      discriminant=param_q**3.D0+param_r**2.D0
      SELECT CASE(discriminant > 0.D0) ! one real solution
         CASE(.TRUE.)
            param_s=param_r+dsqrt(discriminant)
            SELECT CASE(param_s >= 0.D0)
               CASE(.TRUE.)
                  param_s=(param_s)**(1.D0/3.D0)
               CASE(.FALSE.)
                  param_s=-(dabs(param_s))**(1.D0/3.D0)
            END SELECT
            param_t=param_r-dsqrt(discriminant)
            SELECT CASE(param_t >= 0.D0)
               CASE(.TRUE.)
                  param_t=(param_t)**(1.D0/3.D0)
               CASE(.FALSE.)
                  param_t=-(dabs(param_t))**(1.D0/3.D0)
            END SELECT
            ALLOCATE(s(1))
            s(1)=param_s+param_t-a1/3.D0
         CASE(.FALSE.) ! i.e. discriminant <= 0.D0
            param_sc=dcmplx(param_r,0.D0)+sqrt(dcmplx(discriminant,0.D0))
            param_sc=param_sc**(1.D0/3.D0)
            param_tc=dcmplx(param_r,0.D0)-sqrt(dcmplx(discriminant,0.D0))
            param_tc=param_tc**(1.D0/3.D0)
            ALLOCATE(s(3))
            s(1)=real(param_sc+param_tc-dcmplx(a1/3.D0,0.D0))
            s(2)=real(-0.5D0*(param_sc+param_tc)+dcmplx(0.D0,0.5D0*dsqrt(3.D0))*(param_sc-param_tc)-dcmplx(a1/3.D0,0.D0))
            s(3)=real(-0.5D0*(param_sc+param_tc)-dcmplx(0.D0,0.5D0*dsqrt(3.D0))*(param_sc-param_tc)-dcmplx(a1/3.D0,0.D0))
      END SELECT
      DO j = 1, size(s) ! loop over roots
         SELECT CASE(s(j) >= 0.D0 .AND. s(j) <= x2-x1)
            CASE(.TRUE.)
               SELECT CASE(allocated(this%xroot))
                  CASE(.TRUE.)
                     nold=size(this%xroot)
                     ALLOCATE(aux(nold))
                     aux=this%xroot
                     DEALLOCATE(this%xroot)
                     ALLOCATE(this%xroot(nold+1))
                     this%xroot(1:nold)=aux
                     this%xroot(nold+1)=s(j)+x1
                     DEALLOCATE(aux)
                  CASE(.FALSE.)
                     ALLOCATE(this%xroot(1))
                     this%xroot(1)=s(j)+x1
             END SELECT
          CASE(.FALSE.)
            CYCLE
         END SELECT
      END DO
      SELECT CASE(allocated(s))
         CASE(.TRUE.)
            DEALLOCATE(s)
         CASE(.FALSE.)
            ! do nothing
      END SELECT
   END DO
   RETURN
END SUBROUTINE SET_XROOT_CSPLINES

!######################################################################
! SUBROUTINE: DSPLIN ##################################################
!######################################################################
!> @brief
!> For cubic splines interpolation.
!
!> @detailed
!! Given a set of couples @f$(x_{i},F(x_{i}))@f$ and a set of two conditions (read parameters section),
!! this module calculates the set of second derivatives at nodes @f$d^{2}S_{i}(x_{i})\over{dx^{2}}@f$, where:
!! - @f$F(x)=S_{i}(x)@f$ if @f$x\in[x_{i},x_{i+1}]@f$
!! - @f$S_{i}(x)=A_{i}(x-x_{i})^{3}+B_{i}(x-x_{i})^{2}+C_{i}(x-x_{i})+D_{i}@f$.
!!
!! @b Constrains:
!! - @f$S_{i}(x_{i+1})=S_{i+1}(x_{i+1})@f$
!! - @f$\frac{dS_{i}(x_{i+1})}{dx}=\frac{dS_{i+1}(x_{i+1})}{dx}@f$
!! - @f$\frac{d^{2}S_{i}(x_{i+1})}{dx^{2}}=\frac{d^{2}S_{i+1}(x_{i+1})}{dx^{2}}@f$
!! - There should be two extra conditions to get a  system of compatible
!!   linear equations. @b id1 and @b id2 define the kind of conditions that
!!   are supported by this dspline version. @b id1 is the code for @b cond1 (first node condition)
!!   and @b id2 for @b cond2 (last node condition).
!!
!> @param[in,out] cubicspl - Csplines subtype variable. Contains input and output data
!> @param[in] cond1 - if @b id1=1, @b cond1@f$=\frac{d^{2}S(x_{1})}{dx^{2}}@f$; if @b id1=0, @f$S_{1}(x)=S_{2}(x)@f$
!> @param[in] id1 - Integer parameter that controls the meaning of @b cond1
!> @param[in] cond2 - if @b id2=1, @b cond2@f$=\frac{d^{2}S(x_{N})}{dx^{2}}@f$; if @b id2=0, @f$S_{N-1}(x)=S_{N}(x)@f$
!> @param[in] id2 - Integer parameter that controls the meaning of @b cond2
!
!> @author A.P Muzas - alberto.muzas@uam.es
!> @date    20/Jan/2014
!> @version 1.1
!
!> @see debug_mod
!----------------------------------------------------------------------
SUBROUTINE DSPLIN(cubicspl,cond1,id1,cond2,id2)
        IMPLICIT NONE
        ! I/O Variables
        CLASS(Csplines),TARGET,INTENT(INOUT) :: cubicspl
        REAL(KIND=8),INTENT(IN) :: cond1
        INTEGER,INTENT(IN) :: id1
        REAL(KIND=8),INTENT(IN) :: cond2
        INTEGER,INTENT(IN) :: id2
        ! Internal variables
        REAL(KIND=8),DIMENSION(cubicspl%n),TARGET :: diag, indep
        REAL(KIND=8),DIMENSION(cubicspl%n),TARGET :: supdiag, subdiag
        REAL(KIND=8),DIMENSION(cubicspl%n-1) :: h ! space between nodes
        REAL(KIND=8),DIMENSION(cubicspl%n-2) :: sigma
        REAL(KIND=8),DIMENSION(cubicspl%n-2) :: delta
        INTEGER :: i ! Counter
        CHARACTER(LEN=8), PARAMETER :: routinename = "DSPLIN: "
        !Pointers
        INTEGER(KIND=4),POINTER :: n
        REAL(KIND=8),DIMENSION(:),POINTER :: x,y,d2sdx
        REAL(KIND=8),DIMENSION(:),POINTER :: pointdiag, pointindep, pointsupdiag, pointsubdiag
        REAL*8,DIMENSION(:),POINTER :: pointd2sdx
        ! HEY, HO! LET'S GO!
        !=============================
        n => cubicspl%n
        ALLOCATE(cubicspl%d2fdx(1:n))
        x => cubicspl%x(1:n)
        y => cubicspl%f(1:n)
        d2sdx => cubicspl%d2fdx(1:n)
!----------------------------------------------- auxiliar vectors
        ! Get inter-space vector h
        DO i=1,n-1
                h(i)=x(i+1)-x(i)
        END DO
        ! Get sigma and delta auxiliary vectors
        DO i=1, n-2
                sigma(i)=h(i)+h(i+1)
                delta(i)=(1.D0/h(i))+(1.D0/h(i+1))
        END DO
!--------------------------------------------------------
        ! Set diagonal values
        DO i=1, n-2
                diag(i+1)=2.D0*sigma(i)
        END DO
        ! Set supdiag values
        DO i=2,n-1
                supdiag(i)=h(i)
        END DO
        supdiag(n)=0.D0 ! assumption needed by TRIDIA
        ! Set subdiag values
        subdiag(1)=0.D0 ! assumption needed by TRIDIA
        DO i=1,n-2
                subdiag(i+1)=h(i)
        END DO
        ! Set independent terms
        DO i=2, n-1
	        indep(i)=6.D0*((1.D0/h(i-1))*y(i-1)-delta(i-1)*y(i)+(1.D0/h(i))*y(i+1))
        END DO
!------------------- Set diag, supdiag and subdiag parts that depend on id1 and id2
        IF (id1.EQ.1) THEN
                diag(1)=2.D0
                supdiag(1)=1.D0
                indep(1)=(((y(2)-y(1))/h(1))-cond1)*6D0/h(1)
        ELSE IF (id1.EQ.0) THEN
                diag(1)=0.D0 ! Not needed
                diag(2)=2.D0*(sigma(1)+h(1))
                subdiag(2)=0.D0 ! Dummy for tridia
                supdiag(1)=0.D0 ! Not needed
                supdiag(2)=h(2)-h(1)
                ! indep(2) has its usual value
        END IF
!
        IF (id2.EQ.1) THEN
                diag(n)=2.D0
                subdiag(n)=1.D0
                indep(n)=(((y(n-1)-y(n))/h(n-1))+cond2)*6.D0/h(n-1)
        ELSE IF (id2.EQ.0) THEN
                diag(n)=0.D0 ! Not needed
                diag(n-1)=2.D0*(sigma(n-2)+h(n-1))
                supdiag(n-1)=0.D0 ! Dummy for tridia
                subdiag(n)=0.D0 ! Not needed
                subdiag(n-1)=h(n-2)-h(n-1)
                ! indep(n-1) has its usual value
        END IF
!------------------ Set the correct number of dimensions depending on id1 and id2
        IF ((id1.EQ.0).AND.(id2.EQ.0)) THEN
                pointdiag => diag(2:n-1)
                pointsupdiag => supdiag(2:n-1)
                pointsubdiag => subdiag(2:n-1)
                pointindep => indep(2:n-1)
                pointd2sdx => d2sdx(2:n-1)
                CALL TRIDIA(n-2,pointsubdiag,pointdiag,pointsupdiag,pointindep,pointd2sdx)
                d2sdx(1)=2*d2sdx(2)-d2sdx(3)
                d2sdx(n)=2*d2sdx(n-1)-d2sdx(n-2)
        ELSE IF (id1.EQ.0) THEN
                pointdiag => diag(2:n)
                pointsupdiag => supdiag(2:n)
                pointsubdiag => subdiag(2:n)
                pointindep => indep(2:n)
                pointd2sdx => d2sdx(2:n)
                CALL TRIDIA(n-1,pointsubdiag,pointdiag,pointsupdiag,pointindep,pointd2sdx)
                d2sdx(1)=2*d2sdx(1)-d2sdx(2)
        ELSE IF (id2.EQ.0) THEN
                pointdiag => diag(1:n-1)
                pointsupdiag => supdiag(1:n-1)
                pointsubdiag => subdiag(1:n-1)
                pointindep => indep(1:n-1)
                pointd2sdx => d2sdx(1:n-1)
                CALL TRIDIA(n-1,pointsubdiag,pointdiag,pointsupdiag,pointindep,pointd2sdx)
                d2sdx(n)=2*d2sdx(n-1)-d2sdx(n-2)
        ELSE
                CALL TRIDIA(n,subdiag,diag,supdiag,indep,d2sdx)
        END IF
END SUBROUTINE DSPLIN
!######################################################################
! SUBROUTINE: SET_CUBIC_SPLINES_COEFF #################################
!######################################################################
!> @brief
!! Calculates cubic splines coefficients a,b,c,d for a given
!! set of Z(i), F(Z(i)), F''(Z(i))
!
!> @warning
!! - DSPLIN routine should have been executed before, i.e. second derivatives
!!   should be available
!
!> @author A.P Muzas - alberto.muzas@uam.es
!> @date 20/Jan/2014
!> @version 1.0
!
!> @see dsplin
!----------------------------------------------------------------------
SUBROUTINE SET_CUBIC_SPLINES_COEFF(this)
        IMPLICIT NONE
        ! I/O variable ==============================================
        CLASS(Csplines),TARGET,INTENT(INOUT) :: this
        ! Internal Variables
        INTEGER :: i ! counter
        REAL(KIND=8) :: h ! step between nodes
        REAL(KIND=8) :: x1,x2,y1,y2 ! values of nodes and the function defined there
        CHARACTER(LEN=13), PARAMETER :: routinename = "SET_SPLINES: "
        ! Pointers
        INTEGER(KIND=4),POINTER :: n ! Number of nodes
        REAL(KIND=8),DIMENSION(:),POINTER :: z
        REAL(KIND=8),DIMENSION(:),POINTER :: v
        REAL(KIND=8), POINTER :: a,b,c,d ! pointers for coefficients
        REAL(KIND=8), POINTER :: m1,m2 ! pointers for second derivatives
        ! Run section ========================================
        n => this%n
        z => this%x(1:n)
        v => this%f(1:n)
        !
        SELECT CASE(this%is_initialized)
           CASE(.TRUE.)
               WRITE(0,*) "SET_SPLINES ERR: this splines were initialized before"
               CALL EXIT(1)
           CASE(.FALSE.)
              ! do nothing
        END SELECT
        ALLOCATE (this%coeff(n -1,4)) ! There are N-1 Splines and 4 coefficients for each one of them
        DO i=1, n-1
                a => this%coeff(i,1)
                b => this%coeff(i,2)
                c => this%coeff(i,3)
                d => this%coeff(i,4)
                m1 => this%d2fdx(i)
                m2 => this%d2fdx(i+1)
                y1 = v(i)
                y2 = v(i+1)
                x1 = z(i)
                x2 = z(i+1)
!
                h=x2-x1
                a=(m2-m1)/(6.D0*h)
                b=m1/2.D0
                c=((y2-y1)/h)-((m2+2.D0*m1)/6.D0)*h
                d=y1
        END DO
        this%is_initialized=.TRUE.
        RETURN
END SUBROUTINE SET_CUBIC_SPLINES_COEFF
!######################################################################
! SUBROUTINE: INTERPOL_CUBIC_SPLINES ##################################
!######################################################################
!> @brief
!! This subroutine coordinates a complete cubic spline interpolation.
!
!> @details
!! - Contidions @b id1 and @b id2 controls the meaning of @b dz1 and @b dz2. See dsplin
!!   routine
!! - This subroutine allocates csplines%d2vdz(:) with the same size as z(:) or v(:).
!! - Obviously, z(:) and v(:) must have the same size.
!
!> @param[in,out] this - Csplines subtype variable that contains all the cubic splines interpolation data
!> @param[in] dz1 - maybe, second derivative at first node
!> @param[in] id1 - Integer control variable that controls the meaning of @ dz1
!> @param[in] dz2 - maybe, second derivative at last node
!> @param[in] id2 - Integer control variable that controls the meaning of @ dz2
!> @param[in] filename - Output file in which cubic splines coefficients will be printed (optional). If it is not
!!                       given or set to "None", no output file will be printed.
!
!> @author A.P. Muzas - alberto.muzas@uam.es
!> @date 03/Feb/2014
!> @version 1.1
!
!> @see interpol_cubic_splines
!----------------------------------------------------------------------
SUBROUTINE INTERPOL_CUBIC_SPLINES(this,dz1,id1,dz2,id2)
        ! Initial declarations
        IMPLICIT NONE
        ! I/O variables ==================================
        CLASS(Csplines),TARGET,INTENT(INOUT) :: this
        REAL(KIND=8),INTENT(IN) :: dz1,dz2
        INTEGER(KIND=4),INTENT(IN) :: id1,id2
        ! Local variables
        CHARACTER(LEN=24),PARAMETER :: routinename="INTERPOL_CUBIC_SPLINES: "
        ! Pointers
        INTEGER(KIND=4),POINTER :: n
        REAL(KIND=8),DIMENSION(:),POINTER :: z, v
        ! Run section =====================================
        ! Allocate pointers
        n => this%n
        z => this%x(1:n)
        v => this%f(1:n)
        CALL this%SET_SECOND_DERIVS(dz1,id1,dz2,id2)
        CALL this%SET_COEFF()
        RETURN
END SUBROUTINE INTERPOL_CUBIC_SPLINES
!###########################################################
!# FUNCTION: get_csplines_value ############################
!###########################################################
!> @brief
!! Gets @f$F(r)@f$ using cubic interpolation
!
!> @param[in] this - Csplines subtype variable that contains all the information needed
!> @param[in] x - Point to evaluate @f$F(x)@f$
!> @param[in] shift - Optional parameter that if it exists and shift=@f$s@f$, we get @f$F(x+s)@f$ instead of @f$F(x)@f$
!
!> @warning
!! - This function doesn't work if we want to evaluate our fuction outside the range
!!   in which the interpolation was defined. Though, there is a confidence value of 1e-6 units around the
!!   extremes of this range.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 28/Jan/2014
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION get_csplines_value(this,x,shift)
   ! I/O Variables
   CLASS(Csplines),INTENT(IN),TARGET :: this
   REAL(KIND=8),INTENT(IN) :: x
   REAL(KIND=8),OPTIONAL,INTENT(IN) :: shift
   ! Local variables
   INTEGER,POINTER :: n
   REAL(KIND=8),DIMENSION(:),POINTER :: z,v
   REAL(KIND=8):: z1, z2, a, b, c, d
   INTEGER :: i ! Counter
   REAL(KIND=8) :: r
   ! MAY THE FORCE BE WITH YOU --------------------------------------
   n => this%n
   z => this%x(1:n)
   v => this%f(1:n)
   !
   SELECT CASE(present(shift))
      CASE(.TRUE.)
         r=x+shift
      CASE(.FALSE.)
         r=x
   END SELECT
   SELECT CASE(r<this%x(1) .OR. r>this%x(n))
      CASE(.TRUE.)
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: r outside interpolation interval"
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: r = ", r
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: range from ", z(1), "to", z(n)
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   !
   DO i=1,n-1
      z1 = z(i)
      z2 = z(i+1)
      a = this%coeff(i,1)
      b = this%coeff(i,2)
      c = this%coeff(i,3)
      d = this%coeff(i,4)
      SELECT CASE((r<=z2).AND.(r>=z1))
         CASE(.TRUE.)
            EXIT
         CASE(.FALSE.)
            ! do nothing
      END SELECT
   END DO
   ! Now, counter i has the correct label
   get_csplines_value=a*((r-z1)**3.D0)+b*((r-z1)**2.D0)+c*(r-z1)+d
   RETURN
END FUNCTION get_csplines_value
!#################################################################
! FUNCTION: get_csplines_dfdx_value
!> @brief
!! Gets @f$\frac{dF(x)}{dx}@f$ using cubic interpolation
!
!> @param[in] this - Csplines subtype variable that contains all the information needed
!> @param[in] x - Point to evaluate @f$\frac{dF(x)}{dx}@f$
!> @param[in] shift - Optional parameter that if it exists and shift=@f$s@f$,
!!                    we get @f$\frac{dF(x+s)}{dx}@f$ instead of @f$\frac{dF(x)}{dx}@f$
!
!> @warning
!! - This function doesn't work if we want to evaluate our fuction outside the range
!!   in which the interpolation was defined. Though, there is a confidence value of 1e-6 units around the
!!   extremes of this range.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 30/Jan/2014
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION get_csplines_dfdx_value(this,x,shift)
   ! Initial declaration
   IMPLICIT NONE
   ! I/O variables
   CLASS(Csplines),INTENT(IN),TARGET :: this
   REAL(KIND=8),INTENT(IN) :: x
   REAL(KIND=8),OPTIONAL,INTENT(IN) :: shift
   ! Local variables
   INTEGER :: i ! Counter
   REAL(KIND=8) :: r
   ! Pointers
   INTEGER,POINTER :: n
   REAL(KIND=8):: z1,z2,a,b,c,d
   REAL(KIND=8),DIMENSION(:),POINTER :: z,v
   ! Run section -----------------------------------
   n => this%n
   z => this%x(1:n)
   v => this%f(1:n)
   SELECT CASE(present(shift))
      CASE(.TRUE.)
         r=x+shift
      CASE(.FALSE.)
         r=x
   END SELECT
   SELECT CASE(r<this%x(1) .OR. r>this%x(n))
      CASE(.TRUE.)
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: r outside interpolation interval"
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: r = ", r
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: range from ", z(1), "to", z(n)
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   !
   DO i=1,n-1
      z1 = this%x(i)
      z2 = this%x(i+1)
      a = this%coeff(i,1)
      b = this%coeff(i,2)
      c = this%coeff(i,3)
      SELECT CASE((r<=z2).AND.(r>=z1))
         CASE(.TRUE.)
            EXIT
         CASE(.FALSE.)
            ! do nothing
      END SELECT
   END DO
   ! Now, counter i has the correct label
   get_csplines_dfdx_value=3.D0*(a*((r-z1)**2.D0))+2.D0*b*(r-z1)+c
   RETURN
END FUNCTION get_csplines_dfdx_value

REAL(KIND=8) FUNCTION get_csplines_df2dx2_value(this,x,shift)
   ! Initial declaration
   IMPLICIT NONE
   ! I/O variables
   CLASS(Csplines),INTENT(IN),TARGET :: this
   REAL(KIND=8),INTENT(IN) :: x
   REAL(KIND=8),OPTIONAL,INTENT(IN) :: shift
   ! Local variables
   INTEGER :: i ! Counter
   REAL(KIND=8) :: r
   ! Pointers
   INTEGER,POINTER :: n
   REAL(KIND=8):: z1,z2,a,b,c,d
   REAL(KIND=8),DIMENSION(:),POINTER :: z,v
   ! Run section -----------------------------------
   n => this%n
   z => this%x(1:n)
   v => this%f(1:n)
   SELECT CASE(present(shift))
      CASE(.TRUE.)
         r=x+shift
      CASE(.FALSE.)
         r=x
   END SELECT
   SELECT CASE(r<this%x(1) .OR. r>this%x(n))
      CASE(.TRUE.)
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: r outside interpolation interval"
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: r = ", r
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: range from ", z(1), "to", z(n)
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   !
   DO i=1,n-1
      z1 = this%x(i)
      z2 = this%x(i+1)
      a = this%coeff(i,1)
      b = this%coeff(i,2)
      SELECT CASE((r<=z2).AND.(r>=z1))
         CASE(.TRUE.)
            EXIT
         CASE(.FALSE.)
            ! do nothing
      END SELECT
   END DO
   ! Now, counter i has the correct label
   get_csplines_df2dx2_value=6.D0*(a*(r-z1))+2.D0*b
   RETURN
END FUNCTION get_csplines_df2dx2_value
!###############################################################
! SUBROUTINE: GET_V_AND_DERIVS_CSPLINE
!> @brief
!! Computes F(X) and F'(X) at the same time. Better time performance
!> @details
!! Some guidelines: let's call @f$f(x)=f(r)@f$ where @f$r=x+\delta@f$
!! f(x) is the function shifted and f(r) the old one. Take this into
!! account to use correctly the shift option.
!
!> @param[in]  x  - point to evaluate the function (shifted function)
!> @param[out] pot   - f(x)
!> @param[out] deriv - @f$\over{df(x)}{dx}@f$
!> @param[in,optional] shift - Turns on the shifted option
!---------------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_CSPLINES(this,x,pot,deriv,shift)
   IMPLICIT NONE
   CLASS(Csplines),TARGET,INTENT(IN) :: this
   REAL(KIND=8),INTENT(IN) :: x
   REAL(KIND=8),INTENT(OUT) :: pot,deriv
   REAL(KIND=8),OPTIONAL,INTENT(IN) :: shift
   ! Local variables
   INTEGER :: i ! Counter
   REAL(KIND=8) ::  r
   ! Pointers
   INTEGER,POINTER :: n
   REAL(KIND=8) :: z1,z2,a,b,c,d
   REAL(KIND=8),DIMENSION(:),POINTER :: z,v
   ! Run section -----------------------------------
   n => this%n
   z => this%x(1:n)
   v => this%f(1:n)
   SELECT CASE(present(shift))
      CASE(.TRUE.)
         r=x+shift
      CASE(.FALSE.)
         r=x
   END SELECT
   SELECT CASE(r<this%x(1).OR. r>this%x(n))
      CASE(.TRUE.)
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: r outside interpolation interval"
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: r = ", r
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: range from ", z(1), "to", z(n)
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   !
   DO i=1,n-1
      z1=this%x(i)
      z2=this%x(i+1)
      a=this%coeff(i,1)
      b=this%coeff(i,2)
      c=this%coeff(i,3)
      d=this%coeff(i,4)
      SELECT CASE((r<=z2).AND.(r>=z1))
         CASE(.TRUE.)
            EXIT
         CASE(.FALSE.)
            ! do nothing
      END SELECT
   END DO
   ! Now, counter i has the correct label
   pot = a*((r-z1)**3.D0)+b*((r-z1)**2.D0)+c*(r-z1)+d
   deriv = 3.D0*(a*((r-z1)**2.D0))+2.D0*b*(r-z1)+c
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_CSPLINES
END MODULE CUBICSPLINES_MOD

!##################################################################################
! MODULE: FOURIER2D_MOD
!> @brief
!! Provides tools to perform 2D interpolations with fourier series
!##################################################################################
module FOURIER2D_MOD
use SURFACE_MOD, only: Surface
implicit none
!/////////////////////////////////////////////////////////////////
! TYPE: TermcalCulator2d
!> @brief
!! Abstract class to calculate terms of the series avoiding the use of
!! unnecessary switches
!----------------------------------------------------------------
type,abstract:: TermCalculator2d
   contains
      procedure(getvalue_termcalculator_example),public,deferred:: getValue
      procedure(getvalue_termcalculator_example),public,deferred:: getDeriv1
      procedure(getvalue_termcalculator_example),public,deferred:: getDeriv2
end type TermCalculator2d
abstract interface
   !###########################################################
   !# FUNCTION: getvalue_termcalculator_example
   !###########################################################
   !> @brief
   !! Just an example that child objects should override
   !-----------------------------------------------------------
   function getvalue_termcalculator_example(this,k,parity,irrep,x) result(answer)
      import TermCalculator2d
      class(TermCalculator2d),intent(in):: this
      integer(kind=4),dimension(2),intent(in):: k
      real(kind=8),dimension(2),intent(in):: x
      character(len=1),intent(in):: parity
      character(len=2),intent(in):: irrep
      real(kind=8):: answer
   end function getvalue_termcalculator_example
   !-------------------------------------------------------------
end interface
!////////////////////////////////////////////////////////////////
! TYPE: Fourier2d
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
!---------------------------------------------------------------
type,abstract :: Fourier2d
   ! public atributes
   integer(kind=4),public :: n
   integer(kind=4),public :: nfunc
   real(kind=8),dimension(:,:),allocatable,public :: xy
   real(kind=8),dimension(:,:),allocatable,public :: f
   integer(kind=4),dimension(:,:),allocatable,public:: kList
   character(len=1),dimension(:),allocatable,public:: parityList
   character(len=2),dimension(:),allocatable,public:: irrepList
   class(termCalculator2d),allocatable:: term
contains
   ! initialize block
   procedure,public,non_overridable :: read => read_FOURIER2D
   procedure(initializeTerms_FOURIER2D),public,deferred:: initializeTerms
   ! tools block
   procedure(interpol_FOURIER2D),public,deferred :: interpol  ! deferred !!!! take a look to interface
   procedure(get_f_and_derivs_FOURIER2D),public,deferred :: get_f_and_derivs ! deferred !!!! take a look to interface
end type Fourier2d

abstract interface
   !###########################################################
   !# SUBROUTINE: INTERPOL_FOURIER2D
   !###########################################################
   !-----------------------------------------------------------
   subroutine interpol_FOURIER2D(this,filename)
      import fourier2d
      import surface
      class(fourier2d),intent(inout) :: this
      character(len=*),intent(in),optional :: filename
   end subroutine interpol_FOURIER2D
   !###########################################################
   !# SUBROUTINE: GET_F_AND_DERIVS
   !###########################################################
   !-----------------------------------------------------------
   subroutine get_f_and_derivs_FOURIER2D(this,r,v,dvdu)
      import fourier2d
      import surface
      class(fourier2d),intent(in):: this
      real(kind=8),dimension(2),intent(in) :: r
      real(kind=8),dimension(:),intent(out) :: v
      real(kind=8),dimension(:,:),intent(out) :: dvdu
   end subroutine get_f_and_derivs_FOURIER2D
   !###########################################################
   !# SUBROUTINE: INITIALIZETERMS_FOURIER2D
   !###########################################################
   !> @brief
   !! Sets terms for this fourier series. Should be overriden by
   !! child non-abstract classes
   !-----------------------------------------------------------
   subroutine initializeTerms_FOURIER2D(this)
      import Fourier2d
      class(fourier2d),intent(inout):: this
   end subroutine initializeTerms_FOURIER2D
end interface
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
subroutine read_FOURIER2D(this,xy,f,kList,irrepList,parityList)
   ! initial declarations
   implicit none
   ! i/o variables
   class(fourier2d),intent(inout):: this
   real(kind=8),dimension(:,:),intent(in):: xy
   real(kind=8),dimension(:,:),intent(in):: f
   integer(kind=4),dimension(:,:),intent(in):: kList
   character(len=1),dimension(:):: parityList
   character(len=2),dimension(:):: irrepList
   ! local variables
   integer(kind=4) :: ndata,nfunc
   ! run section
   ndata=size(f(1,:))
   nfunc=size(f(:,1))
   select case(size(xy(:,1)) == ndata .and. size(xy(1,:))==2)
      case(.false.)
         write(0,*) "READ_FOURIER2D ERR: dimensions mismatch in arrays xy or f"
         call exit(1)
      case(.true.)
         ! do nothing
   end select
   this%n=ndata
   allocate( this%xy(ndata,2),       source=xy(:,:)       )
   allocate( this%f(nfunc,ndata),    source=f(:,:)        )
   allocate( this%kList(ndata,2),    source=kList(:,:)    )
   allocate( this%parityList(ndata), source=parityList(:) )
   allocate( this%irrepList(ndata),  source=irrepList(:)  )
   this%nfunc = nfunc
   return
end subroutine read_fourier2d

end module FOURIER2D_MOD

!#########################################################
! MODULE: FOURIER_P4MM_MOD
!> @brief
!! Provides tools to genererate a symmetry adapted fourier
!! interpolation
!##########################################################
MODULE FOURIER_P4MM_MOD
use SYSTEM_MOD, only: system_surface
use FOURIER2D_MOD, only: Fourier2d,TermCalculator2d
use UNITS_MOD, only: pi
use MATHS_MOD, only: INV_MTRX
IMPLICIT NONE
!////////////////////////////////////////////////////////////////
! TYPE: TermCalculator2d_p4mm
!> @brief
!! Type extension of TermCalculator for p4mm symmetry
!----------------------------------------------------------------º
type,extends(TermCalculator2d):: TermCalculator2d_p4mm
   contains
      procedure,public:: getValue  => termFoup4mm
      procedure,public:: getDeriv1 => termFoup4mm_dx
      procedure,public:: getDeriv2 => termFoup4mm_dy
end type TermCalculator2d_p4mm
!/////////////////////////////////////////////////////////
! TYPE: Fourierp4mm
!> @brief
!! Extends Fourier2d interpolation for p4mm symmetry
!---------------------------------------------------------
TYPE,EXTENDS(Fourier2d):: Fourierp4mm
   PRIVATE
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: coeff
   CONTAINS
      PROCEDURE,PUBLIC:: INTERPOL => INTERPOL_FOURIERP4MM
      PROCEDURE,PUBLIC:: GET_F_AND_DERIVS => GET_F_AND_DERIVS_FOURIERP4MM
      procedure,public:: initializeTerms => initializeTerms_FOURIERP4MM
END TYPE Fourierp4mm
! variables and types, body
CONTAINS
!###########################################################
!# FUNC: initializeTerms_FOURIERP4MM
!###########################################################
!
!-----------------------------------------------------------
subroutine initializeTerms_FOURIERP4MM(this)
   implicit none
   class(Fourierp4mm),intent(inout):: this
   allocate( TermCalculator2d_p4mm::this%term )
end subroutine initializeTerms_FOURIERP4MM
!###########################################################
!# FUNCTION: termfoup4mm
!###########################################################
!
!-----------------------------------------------------------
function termFoup4mm(this,k,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(TermCalculator2d_p4mm),intent(in):: this
   integer(kind=4),dimension(2),intent(in):: k
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   real(kind=8),dimension(2),intent(in):: x
   ! Dummy out variable
   real(kind=8):: answer
   ! Local variables
   real(kind=8):: g
   ! Paramters
   character(len=*),parameter:: routinename='termfoup4mm: '
   ! Run section -------------------------------------------
   g=2.D0*PI/system_surface%norm_s1
   select case( irrep )
   case('A1')
      select case( parity )
      case('+')
         answer=dcos( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) )+&
                dcos( g*dfloat(k(2))*x(1) )*dcos( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('A2')
      select case( parity )
      case('+')
         answer=-dsin( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) )+&
                 dsin( g*dfloat(k(2))*x(1) )*dsin( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('B1')
      select case( parity )
      case('+')
         answer=dcos( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) )-&
                dcos( g*dfloat(k(2))*x(1) )*dcos( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('B2')
      select case( parity )
      case('+')
         answer=dsin( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) )+&
                dsin( g*dfloat(k(2))*x(1) )*dsin( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('EE','E')
      select case( parity )
      case('-')
         answer=dsin( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) )+&
                dcos( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//' irrep not implemented: "'//irrep//'"'
      call exit(1)
   end select
   return
end function termFoup4mm
!###########################################################
!# FUNCTION: termfoup4mm_dx
!###########################################################
!> @brief
!! ! type brief explanation
!-----------------------------------------------------------
function termfoup4mm_dx(this,k,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(TermCalculator2d_p4mm),intent(in):: this
   integer(kind=4),dimension(2),intent(in):: k
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   real(kind=8),dimension(2),intent(in):: x
   ! Dummy out variable
   real(kind=8):: answer
   ! Local variables
   real(kind=8):: g
   ! Paramters
   character(len=*),parameter:: routinename='termfoup4mm_dx: '
   ! Run section -------------------------------------------
   g=2.D0*PI/system_surface%norm_s1
   select case( irrep )
   case('A1')
      select case( parity )
      case('+')
         answer=-g*( dfloat(k(1))*dsin( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) )+&
                     dfloat(k(2))*dsin( g*dfloat(k(2))*x(1) )*dcos( g*dfloat(k(1))*x(2) ) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('A2')
      select case( parity )
      case('+')
         answer=-g*dfloat(k(1))*dcos( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) )+&
                 g*dfloat(k(2))*dcos( g*dfloat(k(2))*x(1) )*dsin( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('B1')
      select case( parity )
      case('+')
         answer=-g*dfloat(k(1))*dsin( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) )-&
                 g*dfloat(k(2))*dsin( g*dfloat(k(2))*x(1) )*dcos( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('B2')
      select case( parity )
      case('+')
         answer=g*dfloat(k(1))*dcos( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) )+&
                g*dfloat(k(2))*dcos( g*dfloat(k(2))*x(1) )*dsin( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('EE','E')
      select case( parity )
      case('-')
         answer=g*dfloat(k(1))*( dcos( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) )-&
                                 dsin( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) ) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select
   case default
      write(0,*) 'ERR '//routinename//' irrep not implemented: "'//irrep//'"'
      call exit(1)
   end select
   return
end function termfoup4mm_dx
!###########################################################
!# FUNCTION: termfoup4mm_dy
!###########################################################
!> @brief
!! ! type brief explanation
!-----------------------------------------------------------
function termfoup4mm_dy(this,k,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(TermCalculator2d_p4mm),intent(in):: this
   integer(kind=4),dimension(2),intent(in):: k
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   real(kind=8),dimension(2),intent(in):: x
   ! Dummy out variable
   real(kind=8):: answer
   ! Local variables
   real(kind=8):: g
   ! Paramters
   character(len=*),parameter:: routinename='termfoup4mm_dy: '
   ! Run section -------------------------------------------
   g=2.D0*PI/system_surface%norm_s1
   select case( irrep )
   case('A1')
      select case( parity )
      case('+')
         answer=-g*( dfloat(k(2))*dcos( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) )+&
                     dfloat(k(1))*dcos( g*dfloat(k(2))*x(1) )*dsin( g*dfloat(k(1))*x(2) ) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('A2')
      select case( parity )
      case('+')
         answer=-g*dfloat(k(2))*dsin( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) )+&
                 g*dfloat(k(1))*dsin( g*dfloat(k(2))*x(1) )*dcos( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('B1')
      select case( parity )
      case('+')
         answer=-g*dfloat(k(2))*dcos( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) )+&
                 g*dfloat(k(1))*dcos( g*dfloat(k(2))*x(1) )*dsin( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('B2')
      select case( parity )
      case('+')
         answer=g*dfloat(k(2))*dsin( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) )+&
                g*dfloat(k(1))*dsin( g*dfloat(k(2))*x(1) )*dcos( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('EE','E')
      select case( parity )
      case('-')
         answer=g*dfloat(k(2))*( -dsin( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) )+&
                                  dcos( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) ) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select
   case default
      write(0,*) 'ERR '//routinename//' irrep not implemented: "'//irrep//'"'
      call exit(1)
   end select
   return
end function termfoup4mm_dy
!###########################################################
!# SUBROUTINE: INTERPOL_FOURIERP4MM
!###########################################################
!> @brief
!! Interpols with a symmetry adapted fourier series for p4mm
!! wallpaper group
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 28/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INTERPOL_FOURIERP4MM(this,filename)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourierp4mm),INTENT(INOUT) :: this
   CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
   ! Local variables
   INTEGER(KIND=4) :: i,j !counters
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: tmtrx, inv_tmtrx
   ! Run section
   ALLOCATE(tmtrx(this%n,this%n))
   ALLOCATE(inv_tmtrx(this%n,this%n))
   ALLOCATE(this%coeff(this%n,this%nfunc))
   DO i = 1, this%n ! loop over points
      DO j = 1, this%n ! loop over terms (one for each k point)
         tmtrx(i,j)=this%term%getValue( k=this%kList(j,:),parity=this%parityList(j),irrep=this%irrepList(j),x=this%xy(i,:) )
      END DO
   END DO
   CALL INV_MTRX(this%n,tmtrx,inv_tmtrx)
   DO i = 1,this%nfunc ! looop over functions
      this%coeff(:,i)=matmul(inv_tmtrx,this%f(i,:))
   END DO
   SELECT CASE(present(filename)) ! Check if we want to print coefficients
      CASE(.TRUE.)
         OPEN (10,FILE=filename,STATUS="replace",ACTION="write")
         DO i = 1,this%nfunc
            WRITE(10,*) this%coeff(:,i)
         END DO
         CLOSE(10)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE INTERPOL_FOURIERP4MM
!###########################################################
!# SUBROUTINE: GET_F_AND_DERIVS_FOURIERP4MM
!###########################################################
!> @brief
!! Gets all values for functions at a given point r
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_F_AND_DERIVS_FOURIERP4MM(this,r,v,dvdu)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   class(Fourierp4mm),INTENT(IN):: this
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: r
   REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: dvdu
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: terms,terms_dx,terms_dy
   ! Run section
   ALLOCATE(terms(this%n))
   ALLOCATE(terms_dx(this%n))
   ALLOCATE(terms_dy(this%n))
   SELECT CASE( size(v) == this%nfunc .AND. size(dvdu(:,1)) == this%nfunc .AND.&
                size(dvdu(1,:)) == 2)
      CASE(.FALSE.)
         WRITE(0,*) "GET_F_AND_DERIVS_FOURIERP4MM ERR: size mismatch in v and stored values of f"
         WRITE(0,*) "GET_F_AND_DERIVS_FOURIERP4MM ERR: nfunc: ", this%nfunc
         WRITE(0,*) "GET_F_AND_DERIVS_FOURIERP4MM ERR: size(v): ",size(v)
         WRITE(0,*) "GET_F_AND_DERIVS_FOURIERP4MM ERR: size(vdvdu): ",size(dvdu(:,1))
         CALL EXIT(1)
      CASE(.TRUE.)
         ! do nothing
   END SELECT
   DO i = 1, this%nfunc
      DO j = 1, this%n
         terms(j)   =this%term%getValue ( k=this%kList(j,:),parity=this%parityList(j),irrep=this%irrepList(j),x=r )
         terms_dx(j)=this%term%getDeriv1( k=this%kList(j,:),parity=this%parityList(j),irrep=this%irrepList(j),x=r )
         terms_dy(j)=this%term%getDeriv2( k=this%kList(j,:),parity=this%parityList(j),irrep=this%irrepList(j),x=r )
      END DO
      v(i)=dot_product(terms,this%coeff(:,i))
      dvdu(i,1)=dot_product(terms_dx,this%coeff(:,i))
      dvdu(i,2)=dot_product(terms_dy,this%coeff(:,i))
   END DO
   DEALLOCATE(terms)
   DEALLOCATE(terms_dx)
   DEALLOCATE(terms_dy)
   RETURN
END SUBROUTINE GET_F_AND_DERIVS_FOURIERP4MM
END MODULE FOURIER_P4MM_MOD

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
   CHARACTER(LEN=10):: id
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: param
CONTAINS
   ! Initialization block
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: READ => READ_FUNCTION1D
   ! Set block
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_ID => SET_ID_FUNCTION1D
   ! Get block
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: getparam => getparam_function1d
   PROCEDURE(getvalue_function1d),PUBLIC,DEFERRED:: getvalue
   PROCEDURE(getvalue_function1d),PUBLIC,DEFERRED:: getderiv
   ! Plot block
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: PLOT => PLOT_FUNCTION1D
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
   !###########################################################
   !# SUBROUTINE: getderiv_function1d
   !###########################################################
   !-----------------------------------------------------------
   REAL(KIND=8) FUNCTION getderiv_function1d(this,x)
      IMPORT Function1d
      CLASS(Function1d),INTENT(IN)::this
      REAL(KIND=8),INTENT(IN) :: x
   END FUNCTION getderiv_function1d
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
   SELECT CASE(allocated(this%param))
   CASE(.TRUE.)
      SELECT CASE(size(param)/=size(this%param))
      CASE(.TRUE.)
         WRITE(0,*) "READ_FUNCTION1D ERR: wrong dimension of parameters"
         CALL EXIT(1)
      CASE(.FALSE.)
               ! do nothing
      END SELECT
   CASE(.FALSE.)
      ALLOCATE(this%param(size(param)))
   END SELECT
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

!#########################################################
! MODULE: LOGISTIC_FUNCTION
!> @brief
!! Compendium of routines and types to manage a logistic function
!##########################################################
MODULE LOGISTIC_FUNCTION_MOD
USE FUNCTION1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: Logistic_func
!> @brief
!! Function @f$ f(x)=\over{1}{1+e^{ax-b}} @f$, where @b a and @b are
!! parameters
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Function1d) :: Logistic_func
CONTAINS
   PROCEDURE,PUBLIC :: getvalue => getvalue_logistic
   PROCEDURE,PUBLIC :: getderiv => getderiv_logistic
END TYPE Logistic_func
!/////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# FUNCTION: getvalue_logistic
!###########################################################
!> @brief
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getvalue_logistic(this,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Logistic_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: a
   ! Run section
   a=this%getparam()
   getvalue_logistic=1.D0/(1.D0+dexp(a(1)*(x-a(2))))
   RETURN
END FUNCTION getvalue_logistic
!###########################################################
!# FUNCTION: getderiv_logistic
!###########################################################
!> @brief
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_logistic(this,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Logistic_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: a
   ! Run section
   a=this%getparam()
   getderiv_logistic=-a(1)*dexp(a(1)*(x-a(2)))/((1.D0+dexp(a(1)*(x-a(2))))**2.D0)
   RETURN
END FUNCTION getderiv_logistic

END MODULE LOGISTIC_FUNCTION_MOD! contains, body

!#########################################################
! MODULE: ONE_FUNCTION
!> @brief
!! Compendium of routines and types to manage a ONE function
!##########################################################
MODULE ONE_FUNCTION_MOD
USE FUNCTION1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: ONE_func
!> @brief
!! Function @f$ f(x)=\over{1}{1+e^{ax-b}} @f$, where @b a and @b are
!! parameters
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Function1d) :: ONE_func
CONTAINS
   PROCEDURE,PUBLIC :: getvalue => getvalue_ONE
   PROCEDURE,PUBLIC :: getderiv => getderiv_ONE
END TYPE ONE_func
!/////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# FUNCTION: getvalue_ONE
!###########################################################
!> @brief
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getvalue_ONE(this,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(ONE_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   getvalue_ONE=1.D0
   RETURN
END FUNCTION getvalue_ONE
!###########################################################
!# FUNCTION: getderiv_ONE
!###########################################################
!> @brief
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_ONE(this,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(ONE_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   getderiv_ONE=0.D0
   RETURN
END FUNCTION getderiv_ONE

END MODULE ONE_FUNCTION_MOD! contains, body

!#################################################################################
! MODULE: HLiF001_WS_MOD
!> @brief
!! CRP3D specific implementation for H/LiF001 without switch function smoothing
!#################################################################################
module PES_HLIF001_NS_MOD
! Initial declarations
use LiF001SURF_MOD, only: LiF001Surf,pi
use PES_MOD, only: PES
use CUBICSPLINES_MOD, only: Csplines
use FOURIER_P4MM_MOD, only: Fourierp4mm
use ONE_FUNCTION_MOD, only: One_func
implicit none
! Local module variable, used to simulate SYSTEM_MOD
type(LiF001Surf),private:: sysLiF001Surf
!///////////////////////////////////////////////////////////////////////////////
! TYPE & SUBTYPES: Symmetric point
!------------------------------------------------------------------------------
type :: Symmpoint
private
   character(len=:),allocatable:: filename
   character(len=:),allocatable:: alias
   integer(kind=4):: n
   real(kind=8):: x
   real(kind=8):: y
   real(kind=8),dimension(:),allocatable:: z
   real(kind=8),dimension(:),allocatable:: v
   real(kind=8):: dz1
   real(kind=8):: dz2
   type(Csplines),public:: interz
   contains
      procedure,public:: PLOT_DATA => PLOT_DATA_SYMMPOINT
      procedure,public:: PLOT => PLOT_INTERPOL_SYMMPOINT
end type Symmpoint

type,extends(Symmpoint) :: Pair_pot
   private
   integer(kind=4):: id
   real(kind=8):: vasint
	real(kind=8):: rumpling
   contains
      procedure,public:: get_v_and_derivs => get_v_and_derivs_PAIRPOT
end type Pair_pot

type,extends(Symmpoint) :: Sitio
   private
	real(kind=8),dimension(:),allocatable:: dvdx,dvdy,dvdz
end type Sitio
!////////////////////////////////////////////////////////////////////////////////
! TYPE: PES_HLIF001_NS
!------------------------------------------------------------------------------
type,extends(PES) :: PES_HLIF001_NS
   integer(kind=4):: max_order=2
   type(Pair_pot),dimension(:),allocatable:: all_pairpots
   type(Sitio),dimension(:),allocatable:: all_sites
   integer(kind=4),dimension(:,:),allocatable:: klist
   type(One_func) dampFunc
   contains
      ! Initialization block
      procedure,public:: initialize                  => initialize_PES_HLIF001_NS
      ! Get block
      procedure,public:: get_v_and_derivs            => GET_V_AND_DERIVS_PES_HLIF001_NS
      procedure,public:: get_v_and_derivs_correction => GET_V_AND_DERIVS_CORRECTION_PES_HLIF001_NS
      procedure,public:: get_repul_corrections       => GET_REPUL_CORRECTIONS_PES_HLIF001_NS
      procedure,public:: getPot                      => getpot_PES_HLIF001_NS
      ! Enquire block
      procedure,public:: is_allowed                  => is_allowed_PES_HLIF001_NS
      ! Tools block
      procedure,public:: extract_vasint              => EXTRACT_VASINT_PES_HLIF001_NS
      procedure,public:: smooth                      => SMOOTH_PES_HLIF001_NS
      procedure,public:: interpol                    => INTERPOL_Z_PES_HLIF001_NS
      ! Plot tools
      procedure,public:: plot_xymap                  => PLOT_XYMAP_PES_HLIF001_NS
      procedure,public:: plot_direction1d            => PLOT_DIRECTION1D_PES_HLIF001_NS
      procedure,public:: plot_sitios                 => PLOT_SITIOS_PES_HLIF001_NS
      procedure,public:: plot_pairpots               => PLOT_PAIRPOTS_PES_HLIF001_NS
      procedure,public:: plot_z                      => PLOT_Z_PES_HLIF001_NS
end type PES_HLIF001_NS

private initialize_PES_HLiF001_NS,get_v_and_derivs_PES_HLiF001_NS,get_v_and_derivs_correction_PES_HLiF001_NS,&
        get_repul_corrections_PES_HLiF001_NS,is_allowed_PES_HLiF001_NS,extract_vasint_PES_HLiF001_NS,&
        smooth_PES_HLiF001_NS,interpol_Z_PES_HLiF001_NS,plot_XYmap_PES_HLiF001_NS,plot_direction1d_PES_HLiF001_NS,&
        plot_sitios_PES_HLiF001_NS,plot_pairpots_PES_HLiF001_NS,plot_Z_PES_HLiF001_NS,interaction_aenv,&
        interaction_AP,get_v_and_derivs_PAIRPOT,plot_data_SYMMPOINT,plot_interpol_SYMMPOINT

!///////////////////////////////////////////////////////////////////////////
contains
!###########################################################
!# SUBROUTINE: GET_REPUL_CORRECTIONS_PES_HLIF001_NS
!###########################################################
SUBROUTINE GET_REPUL_CORRECTIONS_PES_HLIF001_NS(this,P,v,dvdz,dvdx,dvdy)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_NS),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(3),INTENT(IN) :: P
   REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT):: dvdx,dvdy,dvdz ! corrections to the derivatives
   ! Local variables
   INTEGER(KIND=4) :: npairpots
   INTEGER(KIND=4) :: l,k ! counters
   REAL(KIND=8) :: aux1,aux2,aux3,aux4
   ! Run section
   npairpots=size(this%all_pairpots)
   FORALL(l=1:npairpots)
      v(l)=0.D0
      dvdz(l)=0.D0
      dvdx(l)=0.D0
      dvdy(l)=0.D0
   END FORALL
   DO l = 1, npairpots
      DO k = 0, this%max_order
         CALL INTERACTION_AENV(k,P,this%all_pairpots(l),this%dampfunc,aux1,aux2,aux3,aux4)
         v(l)=v(l)+aux1
         dvdz(l)=dvdz(l)+aux2
         dvdx(l)=dvdx(l)+aux3
         dvdy(l)=dvdy(l)+aux4
      END DO
   END DO
   RETURN
END SUBROUTINE GET_REPUL_CORRECTIONS_PES_HLIF001_NS
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_PAIRPOT
!###########################################################
SUBROUTINE GET_V_AND_DERIVS_PAIRPOT(this,x,v,dvdu)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Pair_pot),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   REAL(KIND=8),INTENT(OUT) :: v,dvdu
   ! Run section
   SELECT CASE(X>this%z(this%n)-this%rumpling)
      CASE(.TRUE.)
         v=0.D0
         dvdu=0.D0
      CASE(.FALSE.)
         CALL this%interz%GET_V_AND_DERIVS(x,v,dvdu,this%rumpling)
   END SELECT
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_PAIRPOT
!###########################################################
!# SUBROUTINE: INITIALIZE_PES_HLIF001_NS
!###########################################################
subroutine initialize_PES_HLIF001_NS(this,filename,tablename)
   implicit none
   ! I/O variables
   class(PES_HLIF001_NS),intent(out):: this
   character(len=*),optional,intent(in):: filename
   character(len=*),optional,intent(in):: tablename
   ! Local variables
   real(kind=8),dimension(129):: commonGrid
   integer(kind=4):: i ! counter
   integer(kind=4),save:: invoked
   data invoked/0/
   ! HEY HO!, LET'S GO!! ------------------
   call this%set_pesType('CRP3D')
   call this%set_alias('H_on_LiF001')
   call this%set_dimensions(3)
   allocate( this%all_pairpots(2) )
   ! Create surface if this was the first call.
   select case( invoked )
   case(0)
      call sysLiF001Surf%initialize('dummyString')
      invoked=1
   case default
      ! do nothing
   end select
   ! Pair potential for Li
   this%all_pairpots(1)%alias='Pairpot_Li'
   this%all_pairpots(1)%vasint=-7.1343585831139720d0
   this%all_pairpots(1)%dz1=-7.74124547149313003d-005
   this%all_pairpots(1)%dz2=0.d0
   this%all_pairpots(1)%id=1
   this%all_pairpots(1)%rumpling=-0.12359375826911972d0
   this%all_pairpots(1)%n=89
   allocate(this%all_pairpots(1)%z(89))
   allocate(this%all_pairpots(1)%v(89))
   this%all_pairpots(1)%z(:)=[ -0.12359375826911972451d0,-0.11338356797313858815d0,-0.10469649694026329778d0,-0.09448630664428216142d0,&
                               -0.08579923561140685717d0,-0.07558904531542573468d0,-0.06690197428255041656d0,-0.05669178398656929407d0,&
                               -0.04800471295369398983d0,-0.03779452265771286734d0,-0.02910745162483756310d0,-0.01889726132885643367d0,&
                               -0.01021019029598112943d0, 0.00000000000000000000d0, 0.01021019029598112943d0, 0.01889726132885643367d0,&
                                0.02910745162483756310d0, 0.03779452265771286734d0, 0.04800471295369398983d0, 0.05669178398656929407d0,&
                                0.06690197428255041656d0, 0.07558904531542573468d0, 0.08579923561140685717d0, 0.09448630664428216142d0,&
                                0.10469649694026329778d0, 0.11338356797313858815d0, 0.12359375826911972451d0, 0.13228082930199502876d0,&
                                0.14249101959797616512d0, 0.15117809063085146937d0, 0.16138828092683260573d0, 0.17007535195970788222d0,&
                                0.18028554225568901859d0, 0.18897261328856432283d0, 0.19918280358454545920d0, 0.20786987461742076344d0,&
                                0.21808006491340189981d0, 0.22676713594627717629d0, 0.23697732624225828491d0, 0.24566439727513364466d0,&
                                0.25587458757111475327d0, 0.26456165860399005751d0, 0.27477184889997119388d0, 0.28345891993284649812d0,&
                                0.29366911022882763449d0, 0.30235618126170293873d0, 0.31256637155768401959d0, 0.37794522657712864566d0,&
                                0.50153898484624837018d0, 0.56691783986569299625d0, 0.69051159813481266525d0, 0.75589045315425729132d0,&
                                0.87948421142337707135d0, 0.94486306644282158640d0, 1.06845682471194125540d0, 1.13383567973138599250d0,&
                                1.25742943800050555048d0, 1.32280829301995028757d0, 1.44640205128906984555d0, 1.51178090630851458265d0,&
                                1.63537466457763436267d0, 1.70075351959707887772d0, 1.82434727786619865775d0, 1.88972613288564317280d0,&
                                2.01331989115476295282d0, 2.07869874617420791196d0, 2.20153094481177458164d0, 2.36215766610705424355d0,&
                                2.83458919932846509226d0, 3.30702073254987549689d0, 3.40150703919415775545d0, 3.77945226577128634560d0,&
                                4.15739749234841582393d0, 4.53534271892554396999d0, 5.29123317207980115029d0, 5.66917839865693018453d0,&
                                6.04712362523405833059d0, 6.42506885181118647665d0, 6.80301407838831551089d0, 7.18095930496544365695d0,&
                                7.55890453154257269119d0, 7.93684975811970172543d0, 8.31479498469683164785d0, 8.69274021127395890574d0,&
                                9.07068543785108793998d0, 9.44863066442821697422d0, 9.82657589100534600846d0,10.01554850429390874922d0,&
                               10.48798003751531915384d0 ]

   this%all_pairpots(1)%v(:)=[  9.13637967606230283479d0, 9.06832750843688373266d0, 8.90326498259120313605d0, 8.58331075848343161283d0,&
                                8.20392411180231384549d0, 7.63206806126462566908d0, 7.03835748948212014398d0, 6.21459984256694042415d0,&
                                5.40656554141708056704d0, 4.33090652817683974263d0, 3.33211761650003612800d0, 2.22360275064672885392d0,&
                                1.46479057654297051272d0, 0.75474225051959265009d0, 0.15214536709775466905d0,-0.30267699281657389765d0,&
                               -0.77687704211649644126d0,-1.13507419544155974123d0,-1.51043841948065837855d0,-1.79707707081678691452d0,&
                               -2.10214263141106361132d0,-2.33735240742591665608d0,-2.58787101157876797686d0,-2.78112055534153856939d0,&
                               -2.98755173052881151108d0,-3.14774025503347854027d0,-3.32054352873101787935d0,-3.45657024697155890181d0,&
                               -3.60620495923014505735d0,-3.72661132147541795945d0,-3.86067137452116160290d0,-3.96863141541037478532d0,&
                               -4.08865106449428239443d0,-4.18516561218743277095d0,-4.29232116241032368720d0,-4.37839104506740905975d0,&
                               -4.47385880153010173643d0,-4.55048484731111813062d0,-4.63544111511443102103d0,-4.70362415217937623879d0,&
                               -4.77924523642412868440d0,-4.83998609293299963952d0,-4.90744829872001098181d0,-4.96174780283280369986d0,&
                               -5.02222743526289239213d0,-5.07108641513960467506d0,-5.12575977931359005879d0,-5.42664945665589293355d0,&
                               -5.85159330426782098300d0,-6.01950817493016820947d0,-6.26631361840675893404d0,-6.36922072891826651642d0,&
                               -6.52670514611862007115d0,-6.59596737715163872195d0,-6.70590532582269460704d0,-6.75459931987306472223d0,&
                               -6.83205122595604752433d0,-6.86637019022921268885d0,-6.92105443619809967970d0,-6.94536416651140164902d0,&
                               -6.98400120097157106613d0,-7.00114267339517137856d0,-7.02835930911559181311d0,-7.04042393089925422345d0,&
                               -7.05965099682119578972d0,-7.06817858395823339634d0,-7.08167042202009078267d0,-7.09526035231560658900d0,&
                               -7.11804121391004507302d0,-7.12743723234768022934d0,-7.12851292897997534936d0,-7.13136969960303535032d0,&
                               -7.13284730060497107473d0,-7.13364775057311284456d0,-7.13434548200933971174d0,-7.13448321515789185554d0,&
                               -7.13456043013040375200d0,-7.13459651023333218944d0,-7.13454700870452107608d0,-7.13457967492788203145d0,&
                               -7.13456476333787215083d0,-7.13453759184147884298d0,-7.13450463614975838311d0,-7.13447005570740167002d0,&
                               -7.13443648744154756969d0,-7.13438097583484864828d0,-7.13438475550669082281d0,-7.13436622608774762000d0,&
                               -7.13435858311397197440d0 ]
   call this%all_pairpots(1)%interZ%read( this%all_pairpots(1)%z(:),this%all_pairpots(1)%v(:) )
   ! Pair potential for F
   this%all_pairpots(2)%alias='Pairpot_F'
   this%all_pairpots(2)%vasint=-7.1343585831139720d0
   this%all_pairpots(2)%dz1=-6.94879417153515533d-5
   this%all_pairpots(2)%dz2=0.d0
   this%all_pairpots(2)%id=2
   this%all_pairpots(2)%rumpling=0.d0
   this%all_pairpots(2)%n=76
   allocate(this%all_pairpots(2)%z(76))
   allocate(this%all_pairpots(2)%v(76))
   this%all_pairpots(2)%z(:)=[ 0.00000000000000000000d0, 0.01021019029598112943d0, 0.01889726132885643367d0, 0.02910745162483756310d0,&
                               0.03779452265771286734d0, 0.04800471295369398983d0, 0.05669178398656929407d0, 0.06690197428255041656d0,&
                               0.07558904531542573468d0, 0.08579923561140685717d0, 0.09448630664428216142d0, 0.10469649694026329778d0,&
                               0.11338356797313858815d0, 0.12359375826911972451d0, 0.13228082930199502876d0, 0.14249101959797616512d0,&
                               0.15117809063085146937d0, 0.16138828092683260573d0, 0.17007535195970788222d0, 0.18028554225568901859d0,&
                               0.18897261328856432283d0, 0.19918280358454545920d0, 0.20786987461742076344d0, 0.21808006491340189981d0,&
                               0.22676713594627717629d0, 0.23697732624225828491d0, 0.24566439727513364466d0, 0.25587458757111475327d0,&
                               0.26456165860399005751d0, 0.27477184889997119388d0, 0.28345891993284649812d0, 0.29366911022882763449d0,&
                               0.30235618126170293873d0, 0.31256637155768401959d0, 0.37794522657712864566d0, 0.50153898484624837018d0,&
                               0.56691783986569299625d0, 0.69051159813481266525d0, 0.75589045315425729132d0, 0.87948421142337707135d0,&
                               0.94486306644282158640d0, 1.06845682471194125540d0, 1.13383567973138599250d0, 1.25742943800050555048d0,&
                               1.32280829301995028757d0, 1.44640205128906984555d0, 1.51178090630851458265d0, 1.63537466457763436267d0,&
                               1.70075351959707887772d0, 1.82434727786619865775d0, 1.88972613288564317280d0, 2.01331989115476295282d0,&
                               2.07869874617420791196d0, 2.20153094481177458164d0, 2.36215766610705424355d0, 2.83458919932846509226d0,&
                               3.30702073254987549689d0, 3.40150703919415775545d0, 3.77945226577128634560d0, 4.15739749234841582393d0,&
                               4.53534271892554396999d0, 5.29123317207980115029d0, 5.66917839865693018453d0, 6.04712362523405833059d0,&
                               6.42506885181118647665d0, 6.80301407838831551089d0, 7.18095930496544365695d0, 7.55890453154257269119d0,&
                               7.93684975811970172543d0, 8.31479498469683164785d0, 8.69274021127395890574d0, 9.07068543785108793998d0,&
                               9.44863066442821697422d0, 9.82657589100534600846d0,10.01554850429390874922d0,10.48798003751531915384d0 ]

   this%all_pairpots(2)%v(:)=[35.94183402487214351595d0,35.89857817326345923448d0,35.79366014961507147518d0,35.59028904596869580246d0,&
                              35.34914120263270831401d0,34.98565488951776103477d0,34.60827726271296000959d0,34.08467578269853959227d0,&
                              33.57106840864373964450d0,32.88735180429894455756d0,32.23751471921294609047d0,31.39368303310690322405d0,&
                              30.60761627320849243006d0,29.60366954791030735805d0,28.68137314941830595671d0,27.51731142749706293671d0,&
                              26.45878542663025001502d0,25.13460875065507593717d0,23.93985318363227321470d0,22.45556159617226654746d0,&
                              21.12457649921226021661d0,19.49007498664534310251d0,18.07575364028142672623d0,16.43792379626155053529d0,&
                              15.10737562624412788637d0,13.66979033954163469389d0,12.58416142588480468589d0,11.45376974357187371822d0,&
                              10.59321451689409343544d0, 9.67263148687921692215d0, 8.94402288586876714760d0, 8.12914701145938423110d0,&
                               7.46794715791469609911d0, 6.72871379453375784152d0, 2.96265525418775155231d0,-1.22114850507602268337d0,&
                              -2.63138763162365441062d0,-4.41765828174076613521d0,-5.04604892046661923644d0,-5.84557849100102178141d0,&
                              -6.13938317512004960719d0,-6.52646072102204932719d0,-6.66451051723715881536d0,-6.84111663060843788742d0,&
                              -6.90141040276800232789d0,-6.97758158903419101193d0,-7.00438098661579378046d0,-7.03814683143213049021d0,&
                              -7.05034999134869000414d0,-7.06674205602035154783d0,-7.07314495695511791240d0,-7.08312610202180614749d0,&
                              -7.08757267028558413102d0,-7.09488257990142123788d0,-7.10291428450212425361d0,-7.11880084201481722062d0,&
                              -7.12665054842541234592d0,-7.12763481242552998651d0,-7.13044848486165694368d0,-7.13212299841474983708d0,&
                              -7.13315234284039600965d0,-7.13416050067259011058d0,-7.13437958950106043687d0,-7.13449736358108133061d0,&
                              -7.13456153972916773398d0,-7.13455976722084450614d0,-7.13458681154711271688d0,-7.13451967913154927459d0,&
                              -7.13453963638738120068d0,-7.13450543479165677496d0,-7.13447613339670905219d0,-7.13442351954518816370d0,&
                              -7.13439290160601746749d0,-7.13438831769229597768d0,-7.13437595607961760891d0,-7.13435858311397197440d0 ]
   call this%all_pairpots(2)%interZ%read( this%all_pairpots(2)%z(:),this%all_pairpots(2)%v(:) )
   ! Set sitios
   ! initialize common numbers
   allocate( this%all_sites(6) )
   this%all_sites(:)%n=129
   this%all_sites(:)%dz2=0.d0
   commonGrid(:)=[-2.20153094481177458164d0,-2.07869874617420791196d0,-2.01331989115476295282d0,-1.88972613288564317280d0,&
                  -1.82434727786619865775d0,-1.70075351959707887772d0,-1.63537466457763436267d0,-1.51178090630851458265d0,&
                  -1.44640205128906984555d0,-1.32280829301995028757d0,-1.25742943800050555048d0,-1.13383567973138599250d0,&
                  -1.06845682471194125540d0,-0.94486306644282158640d0,-0.87948421142337707135d0,-0.75589045315425729132d0,&
                  -0.69051159813481266525d0,-0.56691783986569299625d0,-0.50153898484624837018d0,-0.37794522657712864566d0,&
                  -0.31256637155768401959d0,-0.30235618126170293873d0,-0.29366911022882763449d0,-0.28345891993284649812d0,&
                  -0.27477184889997119388d0,-0.26456165860399005751d0,-0.25587458757111475327d0,-0.24566439727513364466d0,&
                  -0.23697732624225828491d0,-0.22676713594627717629d0,-0.21808006491340189981d0,-0.20786987461742076344d0,&
                  -0.19918280358454545920d0,-0.18897261328856432283d0,-0.18028554225568901859d0,-0.17007535195970788222d0,&
                  -0.16138828092683260573d0,-0.15117809063085146937d0,-0.14249101959797616512d0,-0.13228082930199502876d0,&
                  -0.12359375826911972451d0,-0.11338356797313858815d0,-0.10469649694026329778d0,-0.09448630664428216142d0,&
                  -0.08579923561140685717d0,-0.07558904531542573468d0,-0.06690197428255041656d0,-0.05669178398656929407d0,&
                  -0.04800471295369398983d0,-0.03779452265771286734d0,-0.02910745162483756310d0,-0.01889726132885643367d0,&
                  -0.01021019029598112943d0, 0.00000000000000000000d0, 0.01021019029598112943d0, 0.01889726132885643367d0,&
                   0.02910745162483756310d0, 0.03779452265771286734d0, 0.04800471295369398983d0, 0.05669178398656929407d0,&
                   0.06690197428255041656d0, 0.07558904531542573468d0, 0.08579923561140685717d0, 0.09448630664428216142d0,&
                   0.10469649694026329778d0, 0.11338356797313858815d0, 0.12359375826911972451d0, 0.13228082930199502876d0,&
                   0.14249101959797616512d0, 0.15117809063085146937d0, 0.16138828092683260573d0, 0.17007535195970788222d0,&
                   0.18028554225568901859d0, 0.18897261328856432283d0, 0.19918280358454545920d0, 0.20786987461742076344d0,&
                   0.21808006491340189981d0, 0.22676713594627717629d0, 0.23697732624225828491d0, 0.24566439727513364466d0,&
                   0.25587458757111475327d0, 0.26456165860399005751d0, 0.27477184889997119388d0, 0.28345891993284649812d0,&
                   0.29366911022882763449d0, 0.30235618126170293873d0, 0.31256637155768401959d0, 0.37794522657712864566d0,&
                   0.50153898484624837018d0, 0.56691783986569299625d0, 0.69051159813481266525d0, 0.75589045315425729132d0,&
                   0.87948421142337707135d0, 0.94486306644282158640d0, 1.06845682471194125540d0, 1.13383567973138599250d0,&
                   1.25742943800050555048d0, 1.32280829301995028757d0, 1.44640205128906984555d0, 1.51178090630851458265d0,&
                   1.63537466457763436267d0, 1.70075351959707887772d0, 1.82434727786619865775d0, 1.88972613288564317280d0,&
                   2.01331989115476295282d0, 2.07869874617420791196d0, 2.20153094481177458164d0, 2.36215766610705424355d0,&
                   2.83458919932846509226d0, 3.30702073254987549689d0, 3.40150703919415775545d0, 3.77945226577128634560d0,&
                   4.15739749234841582393d0, 4.53534271892554396999d0, 5.29123317207980115029d0, 5.66917839865693018453d0,&
                   6.04712362523405833059d0, 6.42506885181118647665d0, 6.80301407838831551089d0, 7.18095930496544365695d0,&
                   7.55890453154257269119d0, 7.93684975811970172543d0, 8.31479498469683164785d0, 8.69274021127395890574d0,&
                   9.07068543785108793998d0, 9.44863066442821697422d0, 9.82657589100534600846d0,10.01554850429390874922d0,&
                  10.48798003751531915384d0 ]
   allocate( this%all_sites(1)%z(129),source=commonGrid )
   allocate( this%all_sites(2)%z(129),source=commonGrid )
   allocate( this%all_sites(3)%z(129),source=commonGrid )
   allocate( this%all_sites(4)%z(129),source=commonGrid )
   allocate( this%all_sites(5)%z(129),source=commonGrid )
   allocate( this%all_sites(6)%z(129),source=commonGrid )
   allocate( this%all_sites(1)%v(129) )
   allocate( this%all_sites(2)%v(129) )
   allocate( this%all_sites(3)%v(129) )
   allocate( this%all_sites(4)%v(129) )
   allocate( this%all_sites(5)%v(129) )
   allocate( this%all_sites(6)%v(129) )
   ! initiallize one by one
   ! ** SITIO 1 ******************************************************************************************************************
   this%all_sites(1)%alias='topLi'
   this%all_sites(1)%x=0.d0
   this%all_sites(1)%y=0.d0
   this%all_sites(1)%dz1=8.15821691383332159d-002
   this%all_sites(1)%v(:)=[-6.89842276915349739141d0,-6.89161937084925568797d0,-6.88872153608881632181d0,-6.87622416636372957299d0,&
                           -6.86407350011311390148d0,-6.83245899564393699421d0,-6.81072705934956523777d0,-6.75863145794980368919d0,&
                           -6.72446446780923867692d0,-6.64509003831201372492d0,-6.59362504356348644308d0,-6.47422553351983776082d0,&
                           -6.39731081628031539310d0,-6.22024157177208092406d0,-6.10671890208879553086d0,-5.82123720039279124450d0,&
                           -5.61981863977230933216d0,-5.17486061313926448690d0,-4.78564416523587965457d0,-3.43262358816163848374d0,&
                           -2.05846885231100440450d0,-1.74876595967715209312d0,-1.45738648295269945798d0,-1.07517894177084372132d0,&
                           -0.71010027663856134517d0,-0.22659542356497541782d0, 0.23715634365227189484d0, 0.85189653991610558847d0,&
                            1.46479838533047912463d0, 2.37621681672069540170d0, 3.33211527427179188265d0, 4.50000714669734058759d0,&
                            5.40656816103944315444d0, 6.34612651986772924317d0, 7.03836243470419287149d0, 7.72601781684040567200d0,&
                            8.20392882790527089298d0, 8.63968146340184262044d0, 8.90326776642915085347d0, 9.08711788533851283489d0,&
                            9.13637967606230283479d0, 9.06832750843688373266d0, 8.90326498259120313605d0, 8.58331075848343161283d0,&
                            8.20392411180231384549d0, 7.63206806126462566908d0, 7.03835748948212014398d0, 6.21459984256694042415d0,&
                            5.40656554141708056704d0, 4.33090652817683974263d0, 3.33211761650003612800d0, 2.22360275064672885392d0,&
                            1.46479057654297051272d0, 0.75474225051959265009d0, 0.15214536709775466905d0,-0.30267699281657389765d0,&
                           -0.77687704211649644126d0,-1.13507419544155974123d0,-1.51043841948065837855d0,-1.79707707081678691452d0,&
                           -2.10214263141106361132d0,-2.33735240742591665608d0,-2.58787101157876797686d0,-2.78112055534153856939d0,&
                           -2.98755173052881151108d0,-3.14774025503347854027d0,-3.32054352873101787935d0,-3.45657024697155890181d0,&
                           -3.60620495923014505735d0,-3.72661132147541795945d0,-3.86067137452116160290d0,-3.96863141541037478532d0,&
                           -4.08865106449428239443d0,-4.18516561218743277095d0,-4.29232116241032368720d0,-4.37839104506740905975d0,&
                           -4.47385880153010173643d0,-4.55048484731111813062d0,-4.63544111511443102103d0,-4.70362415217937623879d0,&
                           -4.77924523642412868440d0,-4.83998609293299963952d0,-4.90744829872001098181d0,-4.96174780283280369986d0,&
                           -5.02222743526289239213d0,-5.07108641513960467506d0,-5.12575977931359005879d0,-5.42664945665589293355d0,&
                           -5.85159330426782098300d0,-6.01950817493016820947d0,-6.26631361840675893404d0,-6.36922072891826651642d0,&
                           -6.52670514611862007115d0,-6.59596737715163872195d0,-6.70590532582269460704d0,-6.75459931987306472223d0,&
                           -6.83205122595604752433d0,-6.86637019022921268885d0,-6.92105443619809967970d0,-6.94536416651140164902d0,&
                           -6.98400120097157106613d0,-7.00114267339517137856d0,-7.02835930911559181311d0,-7.04042393089925422345d0,&
                           -7.05965099682119578972d0,-7.06817858395823339634d0,-7.08167042202009078267d0,-7.09526035231560658900d0,&
                           -7.11804121391004507302d0,-7.12743723234768022934d0,-7.12851292897997534936d0,-7.13136969960303535032d0,&
                           -7.13284730060497107473d0,-7.13364775057311284456d0,-7.13434548200933971174d0,-7.13448321515789185554d0,&
                           -7.13456043013040375200d0,-7.13459651023333218944d0,-7.13454700870452107608d0,-7.13457967492788203145d0,&
                           -7.13456476333787215083d0,-7.13453759184147884298d0,-7.13450463614975838311d0,-7.13447005570740167002d0,&
                           -7.13443648744154756969d0,-7.13438097583484864828d0,-7.13438475550669082281d0,-7.13436622608774762000d0,&
                           -7.13435858311397197440d0 ]
   ! ** SITIO 2 ******************************************************************************************************************
   this%all_sites(2)%alias='topF'
   this%all_sites(2)%x=2.7216780628885480d0
   this%all_sites(2)%y=2.7216780628885480d0
   this%all_sites(2)%dz1=-0.15738884744871576d0
   this%all_sites(2)%v=[-6.87907498431797836957d0,-6.90118218389560222903d0,-6.91357281358908082325d0,-6.93083748959737988571d0,&
                        -6.93516974660099094763d0,-6.93620262739340098790d0,-6.93254567100649055078d0,-6.91271774527259896814d0,&
                        -6.89236865850650648468d0,-6.82601138191416811907d0,-6.76901355295993312211d0,-6.60134779107096303363d0,&
                        -6.46838559924402023427d0,-6.08818074569892431214d0,-5.79753398418044074702d0,-4.98020683906518435435d0,&
                        -4.35181126919280814036d0,-2.62888018365377318375d0,-1.22625825492317463983d0, 2.96302449612861984463d0,&
                         6.72861752598492746102d0, 7.46788420798443475235d0, 8.12911587850032901770d0, 8.94402288586876537124d0,&
                         9.67264601750633801203d0,10.59323358356773425726d0,11.45378545245022117172d0,12.58416898952393481181d0,&
                        13.66979085473309929455d0,15.10737121412130079534d0,16.43791805939430261674d0,18.07574875543115311416d0,&
                        19.49007210096592856985d0,21.12457649921226732204d0,22.45556384681248829338d0,23.93985767454301338830d0,&
                        25.13461482069859087574d0,26.45879299026938724637d0,27.51731997119913586403d0,28.68138252517930908425d0,&
                        29.60367937710202212998d0,30.60762635806067777366d0,31.39369311719515565073d0,32.23752456770139929176d0,&
                        32.88736127026645306159d0,33.57107723288939382655d0,34.08468391510380968157d0,34.60828443241255314433d0,&
                        34.98566113049513859323d0,35.34914624505879032768d0,35.59029299522835287917d0,35.79366274961602556459d0,&
                        35.89857958809135141109d0,35.94183402487214351595d0,35.89857817326345923448d0,35.79366014961507147518d0,&
                        35.59028904596869580246d0,35.34914120263270831401d0,34.98565488951776103477d0,34.60827726271296000959d0,&
                        34.08467578269853959227d0,33.57106840864373964450d0,32.88735180429894455756d0,32.23751471921294609047d0,&
                        31.39368303310690322405d0,30.60761627320849243006d0,29.60366954791030735805d0,28.68137314941830595671d0,&
                        27.51731142749706293671d0,26.45878542663025001502d0,25.13460875065507593717d0,23.93985318363227321470d0,&
                        22.45556159617226654746d0,21.12457649921226021661d0,19.49007498664534310251d0,18.07575364028142672623d0,&
                        16.43792379626155053529d0,15.10737562624412788637d0,13.66979033954163469389d0,12.58416142588480468589d0,&
                        11.45376974357187371822d0,10.59321451689409343544d0, 9.67263148687921692215d0, 8.94402288586876714760d0,&
                         8.12914701145938423110d0, 7.46794715791469609911d0, 6.72871379453375784152d0, 2.96265525418775155231d0,&
                        -1.22114850507602268337d0,-2.63138763162365441062d0,-4.41765828174076613521d0,-5.04604892046661923644d0,&
                        -5.84557849100102178141d0,-6.13938317512004960719d0,-6.52646072102204932719d0,-6.66451051723715881536d0,&
                        -6.84111663060843788742d0,-6.90141040276800232789d0,-6.97758158903419101193d0,-7.00438098661579378046d0,&
                        -7.03814683143213049021d0,-7.05034999134869000414d0,-7.06674205602035154783d0,-7.07314495695511791240d0,&
                        -7.08312610202180614749d0,-7.08757267028558413102d0,-7.09488257990142123788d0,-7.10291428450212425361d0,&
                        -7.11880084201481722062d0,-7.12665054842541234592d0,-7.12763481242552998651d0,-7.13044848486165694368d0,&
                        -7.13212299841474983708d0,-7.13315234284039600965d0,-7.13416050067259011058d0,-7.13437958950106043687d0,&
                        -7.13449736358108133061d0,-7.13456153972916773398d0,-7.13455976722084450614d0,-7.13458681154711271688d0,&
                        -7.13451967913154927459d0,-7.13453963638738120068d0,-7.13450543479165677496d0,-7.13447613339670905219d0,&
                        -7.13442351954518816370d0,-7.13439290160601746749d0,-7.13438831769229597768d0,-7.13437595607961760891d0,&
                        -7.13435858311397197440d0 ]
   ! ** SITIO 3 ******************************************************************************************************************
   this%all_sites(3)%alias='Hollow'
   this%all_sites(3)%x=2.7216780628885480d0
   this%all_sites(3)%y=0.d0
   this%all_sites(3)%dz1=-4.69246625508239895d-003
   this%all_sites(3)%v=[-7.06087793948749187223d0,-7.06169928158901782922d0,-7.06219152067239797077d0,-7.06254771731697594817d0,&
                        -7.06230569380749084729d0,-7.06130296037807880793d0,-7.06048918499986566388d0,-7.05844921961311921876d0,&
                        -7.05712441466509954324d0,-7.05420750104602145569d0,-7.05247785656989378822d0,-7.04893341285345265135d0,&
                        -7.04695256663319824497d0,-7.04311920223866660962d0,-7.04109375051585306693d0,-7.03739760572271233485d0,&
                        -7.03557201745042970487d0,-7.03250261019674116625d0,-7.03113752572527861417d0,-7.02917545358314743709d0,&
                        -7.02851303679548156111d0,-7.02843474557034308958d0,-7.02837362338697158037d0,-7.02830829317233440889d0,&
                        -7.02825829172660654365d0,-7.02820614161327394953d0,-7.02816744731282039282d0,-7.02812869639159831792d0,&
                        -7.02810149564404884615d0,-7.02807636300574323229d0,-7.02806084221872939821d0,-7.02804954695414352273d0,&
                        -7.02804589253529687909d0,-7.02804865373523579564d0,-7.02805705209218789520d0,-7.02807408884745665745d0,&
                        -7.02809472638783727660d0,-7.02812624659455220666d0,-7.02815924176018924641d0,-7.02820526380243304487d0,&
                        -7.02825056717226370750d0,-7.02831101981917161226d0,-7.02836857078330545789d0,-7.02844338279815072923d0,&
                        -7.02851312074669909435d0,-7.02860222089275499258d0,-7.02868408521582743731d0,-7.02878740225636722272d0,&
                        -7.02888133234407597172d0,-7.02899879504237201644d0,-7.02910473028482662983d0,-7.02923626740415308234d0,&
                        -7.02935414719146312024d0,-7.02949968749509324084d0,-7.02965276407674188874d0,-7.02978892346857708873d0,&
                        -7.02995589350859884092d0,-7.03010384347798833460d0,-7.03028463573965645139d0,-7.03044431567670979888d0,&
                        -7.03063885892329754057d0,-7.03081020821812696653d0,-7.03101843121290759342d0,-7.03120138925562176979d0,&
                        -7.03142321131777325860d0,-7.03161766706620916523d0,-7.03185287691142857369d0,-7.03205861042142554851d0,&
                        -7.03230694270577050276d0,-7.03252372821643678691d0,-7.03278491759596757760d0,-7.03301252934641052406d0,&
                        -7.03328631047718388913d0,-7.03352452270651262722d0,-7.03381063024458796917d0,-7.03405921719190896368d0,&
                        -7.03435738579334479681d0,-7.03461612169776806525d0,-7.03492608601862290385d0,-7.03519474511925491100d0,&
                        -7.03551623981558726939d0,-7.03579459635153625641d0,-7.03612735607940464888d0,-7.03641518428977974509d0,&
                        -7.03675894370524268595d0,-7.03705601782915124431d0,-7.03741051158826635969d0,-7.03981189318787325959d0,&
                        -7.04490278290859262000d0,-7.04783667739777364147d0,-7.05372072931313365274d0,-7.05696307030613301237d0,&
                        -7.06322640526540901362d0,-7.06656469565906686370d0,-7.07283420489514291774d0,-7.07609244097356171466d0,&
                        -7.08207098930207035181d0,-7.08511606054587517178d0,-7.09060273906750815343d0,-7.09334866629683524764d0,&
                        -7.09822272454889890980d0,-7.10062582835543842918d0,-7.10483581515609596835d0,-7.10688803045187977858d0,&
                        -7.11044455233813277317d0,-7.11216182168590460577d0,-7.11509771302703875051d0,-7.11840651508443755802d0,&
                        -7.12530839366201007579d0,-7.12932443281130900914d0,-7.12989992693555585390d0,-7.13167484418369568289d0,&
                        -7.13282403546793641880d0,-7.13355536697153524983d0,-7.13428594963942863671d0,-7.13444636726879277688d0,&
                        -7.13452684829341698958d0,-7.13455128719486975086d0,-7.13455026786885415646d0,-7.13453102855387655268d0,&
                        -7.13450641151416320440d0,-7.13447261401719412532d0,-7.13450772409165523413d0,-7.13446951359362380174d0,&
                        -7.13443884441551734454d0,-7.13440949647940936273d0,-7.13438305676662221089d0,-7.13436315578209434562d0,&
                        -7.13435858311397197440d0 ]
   ! ** SITIO 4 ******************************************************************************************************************
   this%all_sites(4)%alias='Bridge'
   this%all_sites(4)%x=1.3608390314442740d0
   this%all_sites(4)%y=1.3608390314442740d0
   this%all_sites(4)%dz1=-9.34171716476418432d-003
   this%all_sites(4)%v=[-7.01734542360757362900d0,-7.01910913898698396451d0,-7.02018644062349839174d0,-7.02077900195917692372d0,&
                        -7.02000975784999248219d0,-7.01718135589118396922d0,-7.01496951631789844583d0,-7.00950107625823015667d0,&
                        -7.00596464216334169350d0,-6.99817010132702410630d0,-6.99353725531451431152d0,-6.98401109877033920981d0,&
                        -6.97866848849064069782d0,-6.96829902813987800414d0,-6.96280583302785593958d0,-6.95275170578670298482d0,&
                        -6.94776609507806952593d0,-6.93934010648358601259d0,-6.93557474808632790086d0,-6.93014506921935424799d0,&
                        -6.92829533924338214490d0,-6.92807363699524980660d0,-6.92789952433594979908d0,-6.92771203948999847455d0,&
                        -6.92756719420705380230d0,-6.92741428971156558703d0,-6.92729902060994806590d0,-6.92718105941320683172d0,&
                        -6.92709567529788738938d0,-6.92701302034817611997d0,-6.92695783002412657225d0,-6.92691084426973002763d0,&
                        -6.92688615654192041404d0,-6.92687520293112157788d0,-6.92688132660452460243d0,-6.92690676808560557021d0,&
                        -6.92694401196519216057d0,-6.92700619249363924723d0,-6.92707475007426243252d0,-6.92717369208129962743d0,&
                        -6.92727347206327959128d0,-6.92740904594029061059d0,-6.92753993804109313714d0,-6.92771201416951321050d0,&
                        -6.92787390810660408391d0,-6.92808235686787199370d0,-6.92827514235871610992d0,-6.92851983413426886216d0,&
                        -6.92874340089633289352d0,-6.92902420606760660604d0,-6.92927844381835633669d0,-6.92959523276678979187d0,&
                        -6.92988003122369011777d0,-6.93023267433071943344d0,-6.93060463550379335373d0,-6.93093629085830009728d0,&
                        -6.93134390653849763453d0,-6.93170584244843368538d0,-6.93214898301719273377d0,-6.93254108920002298788d0,&
                        -6.93301962503878144162d0,-6.93344179121197168314d0,-6.93395559270216566006d0,-6.93440770858318256131d0,&
                        -6.93495660989368989391d0,-6.93543837182184219614d0,-6.93602170633388492860d0,-6.93653239307327051932d0,&
                        -6.93714928688110177291d0,-6.93768815489207302249d0,-6.93833773408994503029d0,-6.93890403983285342093d0,&
                        -6.93958543051501930421d0,-6.94017843045021631809d0,-6.94089075871092830994d0,-6.94150970929876631743d0,&
                        -6.94225210123227576275d0,-6.94289625893310624605d0,-6.94366784063366804247d0,-6.94433646190784248375d0,&
                        -6.94513635946970797619d0,-6.94582870077757696947d0,-6.94665604029500016736d0,-6.94737135809691608301d0,&
                        -6.94822526566414921945d0,-6.94896281642046353966d0,-6.94984241813175884772d0,-6.95578404806645522740d0,&
                        -6.96825362467727149607d0,-6.97534315145011341031d0,-6.98930870757213540401d0,-6.99685482902732314159d0,&
                        -7.01111953834615064807d0,-7.01854500732761099613d0,-7.03214963133864934264d0,-7.03903653854122612898d0,&
                        -7.05132794945392404884d0,-7.05741256654557247430d0,-7.06805970783239256860d0,-7.07323010665734486224d0,&
                        -7.08214003401021496842d0,-7.08640261928132897395d0,-7.09365139856436055510d0,-7.09707958319132803382d0,&
                        -7.10284741693550714103d0,-7.10554885788336765273d0,-7.11003367320138313090d0,-7.11487183210872853323d0,&
                        -7.12413886648707173066d0,-7.12891881720788500587d0,-7.12956415865976378399d0,-7.13150387492979742632d0,&
                        -7.13273050712688938546d0,-7.13350451343692704143d0,-7.13427180968791940785d0,-7.13444012401582039473d0,&
                        -7.13453257670234997079d0,-7.13455678795222070221d0,-7.13456465454598465215d0,-7.13458469408834439207d0,&
                        -7.13456375501156703223d0,-7.13453593656660700617d0,-7.13450362781379521238d0,-7.13447059314077058900d0,&
                        -7.13443747606051292820d0,-7.13438360160353823858d0,-7.13438656093252809853d0,-7.13436967307073999223d0,&
                        -7.13435858311397197440d0 ]
   ! ** SITIO 5 ******************************************************************************************************************
   this%all_sites(5)%alias='halfHollowF'
   this%all_sites(5)%x=2.7216780628885480d0
   this%all_sites(5)%y=1.3608390314442740d0
   this%all_sites(5)%dz1=-3.73455010716034277E-002
   this%all_sites(5)%v=[-7.00416123562780068568d0,-7.00955367691972863753d0,-7.01260489108868245722d0,-7.01641734487907431372d0,&
                        -7.01700440152958559992d0,-7.01651325928337765703d0,-7.01548043652173447526d0,-7.01222306308404519370d0,&
                        -7.00987367523607929343d0,-7.00435546791209873163d0,-7.00085351035068281789d0,-6.99306728938978938714d0,&
                        -6.98825097231746283910d0,-6.97745570154741390212d0,-6.97069193538925802045d0,-6.95594980396307782655d0,&
                        -6.94725006023632829510d0,-6.92955211827853467810d0,-6.91998249318655922480d0,-6.90273203430610848841d0,&
                        -6.89463036209347368555d0,-6.89346159578449402261d0,-6.89249072252290151397d0,-6.89137851896055053658d0,&
                        -6.89045780299273502578d0,-6.88940693960230454707d0,-6.88854044567453360770d0,-6.88755569988131544079d0,&
                        -6.88674749273985220555d0,-6.88583364196914082811d0,-6.88508778636025020603d0,-6.88424960803733654302d0,&
                        -6.88357016870728521951d0,-6.88281244025746197224d0,-6.88220348195251307999d0,-6.88153098080107294976d0,&
                        -6.88099656826749228600d0,-6.88041391210070507611d0,-6.87995714026848759204d0,-6.87946624259136463309d0,&
                        -6.87908781112477196729d0,-6.87868930671053036008d0,-6.87838975561844190310d0,-6.87808427915665632213d0,&
                        -6.87786414844795324086d0,-6.87765233462819836063d0,-6.87751216431176004562d0,-6.87739464782361142881d0,&
                        -6.87733497790831727059d0,-6.87731239344135136804d0,-6.87733376393608164534d0,-6.87740674617987313155d0,&
                        -6.87750969709350812309d0,-6.87767888073763167256d0,-6.87790030610033209513d0,-6.87812997181308283245d0,&
                        -6.87844857618846905467d0,-6.87876119410468245263d0,-6.87917761218232914899d0,-6.87957372231088370995d0,&
                        -6.88008858878036821949d0,-6.88056873113014422216d0,-6.88118268068103944302d0,-6.88174739526091805430d0,&
                        -6.88246083864992641566d0,-6.88310946964814718285d0,-6.88391972083849257302d0,-6.88464903022272256550d0,&
                        -6.88555212134704142102d0,-6.88635873316202129502d0,-6.88735069635295005241d0,-6.88823123464341957600d0,&
                        -6.88930810203359555999d0,-6.89025919084429627759d0,-6.89141699456635503651d0,-6.89243525794202582802d0,&
                        -6.89367003012860557476d0,-6.89475209211398798459d0,-6.89605986489772515569d0,-6.89720234953755895191d0,&
                        -6.89857915505108998389d0,-6.89977868639011493457d0,-6.90122055676607804031d0,-6.90247375884903391352d0,&
                        -6.90397672622006552956d0,-6.90528022309169386972d0,-6.90684031959042954441d0,-6.91747259100995393055d0,&
                        -6.93968723618298444222d0,-6.95191528833970462387d0,-6.97452606792964413529d0,-6.98587351099715814229d0,&
                        -7.00567168343398094521d0,-7.01511470328137054508d0,-7.03109323158619670124d0,-7.03863076138616872157d0,&
                        -7.05134496708400249076d0,-7.05735395370780960178d0,-7.06757698356927654970d0,-7.07245862109246115068d0,&
                        -7.08083719822429724644d0,-7.08487585467741887157d0,-7.09183152137041084018d0,-7.09517106726483071100d0,&
                        -7.10088438295883417339d0,-7.10360576328006843028d0,-7.10819313166986521679d0,-7.11324455869295046995d0,&
                        -7.12317186210545827407d0,-7.12835846524138982971d0,-7.12905884224952579586d0,-7.13116780839535913117d0,&
                        -7.13251175715024032797d0,-7.13336800341086707533d0,-7.13422453218720509227d0,-7.13441279188392218913d0,&
                        -7.13451349791033795356d0,-7.13454958681170925416d0,-7.13455374752996007715d0,-7.13454102811457069322d0,&
                        -7.13452485239696976294d0,-7.13453792143908316348d0,-7.13450400840355491994d0,-7.13447078765003528389d0,&
                        -7.13443782297732287390d0,-7.13438402962089579518d0,-7.13438660653671785639d0,-7.13437005762999731928d0,&
                        -7.13435858311397197440d0 ]
   ! ** SITIO 6 ******************************************************************************************************************
   this%all_sites(6)%alias='halfHollowLi'
   this%all_sites(6)%x=1.3608390314442740d0
   this%all_sites(6)%y=0.d0
   this%all_sites(6)%dz1=1.09526892073127137d-002
   this%all_sites(6)%v=[-7.00415315077621336570d0,-7.00373007401976011010d0,-7.00371222750970101600d0,-7.00158156764360217750d0,&
                        -6.99884809597901380585d0,-6.99144284982398733774d0,-6.98630178762169151696d0,-6.97424980496590851686d0,&
                        -6.96664818123409190775d0,-6.95001169752298153526d0,-6.94008330528329153708d0,-6.91937461892249938700d0,&
                        -6.90753245408903726599d0,-6.88398200796059622775d0,-6.87114413139800817021d0,-6.84691136678328504672d0,&
                        -6.83449730208988448510d0,-6.81271747551268447296d0,-6.80253598931567982078d0,-6.78688887092610571017d0,&
                        -6.78091823401417048700d0,-6.78014447738610659400d0,-6.77952109958926651956d0,-6.77882999250059015139d0,&
                        -6.77827774306756669631d0,-6.77767118528644463282d0,-6.77719168058485887940d0,-6.77667157187945701224d0,&
                        -6.77626642827693181914d0,-6.77583466841541781633d0,-6.77550550227957337768d0,-6.77516399103011313088d0,&
                        -6.77491241872857319350d0,-6.77466305585933348254d0,-6.77449069375971824059d0,-6.77433537903886584530d0,&
                        -6.77424384350879815742d0,-6.77418436661899470153d0,-6.77417460566897844387d0,-6.77421089268341081180d0,&
                        -6.77428220360623178919d0,-6.77441329934921032674d0,-6.77456486940979107914d0,-6.77478981864798601720d0,&
                        -6.77502083511124908455d0,-6.77533868261132887767d0,-6.77564833274219768811d0,-6.77605812327083167901d0,&
                        -6.77644559433422966066d0,-6.77694637265808630389d0,-6.77741085191893599671d0,-6.77800166280468552316d0,&
                        -6.77854233752791035528d0,-6.77922222574222210767d0,-6.77994984495306063366d0,-6.78060629350228616374d0,&
                        -6.78142151869578757584d0,-6.78215209811647223859d0,-6.78305397405930499133d0,-6.78385787161637043852d0,&
                        -6.78484544307520476281d0,-6.78572184603357442256d0,-6.78679415777507966112d0,-6.78774225339967607340d0,&
                        -6.78889825158255533211d0,-6.78991670056336005246d0,-6.79115396768542023409d0,-6.79224029364168302436d0,&
                        -6.79355584775164356159d0,-6.79470751356879265614d0,-6.79609837271537475800d0,-6.79731284127883750301d0,&
                        -6.79877602351076326670d0,-6.80005075770596789653d0,-6.80158328107195586654d0,-6.80291574378433150372d0,&
                        -6.80451462633310377726d0,-6.80590228044807687979d0,-6.80756454022835288953d0,-6.80900484863135346814d0,&
                        -6.81072750369185442310d0,-6.81221792926830982395d0,-6.81399799765775604499d0,-6.81553600329309450245d0,&
                        -6.81737050306020631041d0,-6.81895355163985605884d0,-6.82083950083335377457d0,-6.83350817178206337843d0,&
                        -6.85963235225001088935d0,-6.87418257984575831188d0,-6.90212676042115180763d0,-6.91682497194263046225d0,&
                        -6.94382795159492527404d0,-6.95746312331058813783d0,-6.98173161042099543749d0,-6.99367754865614532633d0,&
                        -7.01444370813055684266d0,-7.02446678463968154915d0,-7.04160795308638132184d0,-7.04974944331464303815d0,&
                        -7.06349800580577191766d0,-7.06995074865808703635d0,-7.08073054425141457102d0,-7.08574141145314584378d0,&
                        -7.09404216565429646124d0,-7.09786991400715461253d0,-7.10413636116513114160d0,-7.11076011700106569435d0,&
                        -7.12292967161225032413d0,-7.12876940568110395446d0,-7.12951179286446734551d0,-7.13163480918870185121d0,&
                        -7.13287472732416283350d0,-7.13361380054054006195d0,-7.13431574784245281506d0,-7.13446842325632335502d0,&
                        -7.13453236080935759844d0,-7.13459074106095947343d0,-7.13456562046095488228d0,-7.13454334766708697657d0,&
                        -7.13457008178702434265d0,-7.13453431001571392045d0,-7.13450272014070474569d0,-7.13446947796920216689d0,&
                        -7.13443632521088844101d0,-7.13438203038740947903d0,-7.13438533987800926894d0,-7.13436779344102589562d0,&
                        -7.13435858311397197440d0 ]
   ! read interZ
   do i=1,6
      call this%all_sites(i)%interZ%read( this%all_sites(i)%z(:),this%all_sites(i)%v(:) )
   enddo
   ! Set Fourier Coefficients
   allocate( this%klist(6,2) )
   this%klist(:,1)=[0,1,1,2,2,2]
   this%klist(:,2)=[0,0,1,0,1,2]
   call this%interpol()
   return
end subroutine initialize_PES_HLIF001_NS
!#######################################################################
!# SUBROUTINE: EXTRACT_VASINT_PES_HLIF001_NS ######################################
!#######################################################################
SUBROUTINE EXTRACT_VASINT_PES_HLIF001_NS(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_NS),INTENT(INOUT):: this
   ! Local variables
   INTEGER(KIND=4):: npairpots, nsites
   INTEGER(KIND=4):: i,j ! counters
   REAL(KIND=8):: control_vasint
   character(len=*),parameter:: routinename="EXTRACT_VASINT_PES_HLIF001_NS: "
   ! Run section ------------------------
   npairpots=size(this%all_pairpots)
   control_vasint=this%all_pairpots(1)%vasint
   DO i = 1, npairpots
      IF (this%all_pairpots(1)%vasint/=control_vasint) THEN
         WRITE(0,*) "EXTRACT_VASINT_PES_HLIF001_NS ERR: Incoherences in vasint values found"
         WRITE(0,*) "EXTRACT_VASINT_PES_HLIF001_NS ERR: vasint's value at pairpot",1,control_vasint
         WRITE(0,*) "EXTRACT_VASINT_PES_HLIF001_NS ERR: vasint's value at pairpot",i,control_vasint
         CALL EXIT(1)
      END IF
      DO j = 1, this%all_pairpots(i)%n
         this%all_pairpots(i)%interz%f(j)=this%all_pairpots(i)%interz%f(j)-this%all_pairpots(i)%vasint
      END DO
   END DO
   nsites=size(this%all_sites)
   DO i = 1, nsites
      DO j = 1, this%all_sites(i)%n
         this%all_sites(i)%interz%f(j)=this%all_sites(i)%interz%f(j)-this%all_pairpots(1)%vasint
      END DO
   END DO
   RETURN
END SUBROUTINE EXTRACT_VASINT_PES_HLIF001_NS
!############################################################
!# SUBROUTINE: SMOOTH_PES_HLIF001_NS ############################
!############################################################
SUBROUTINE SMOOTH_PES_HLIF001_NS(this)
   ! Initial declaraitons
   IMPLICIT NONE
   CLASS(PES_HLIF001_NS),INTENT(INOUT):: this
   ! Local variables
   REAL(KIND=8),DIMENSION(3):: A
   INTEGER(KIND=4):: i,j ! counters
   INTEGER(KIND=4):: npairpots,nsites
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: v,dvdzr,dummy
   CHARACTER(LEN=*),PARAMETER:: routinename="SMOOTH_PES_HLIF001_NS: "
   ! Run section ----------
   nsites = size(this%all_sites)
   npairpots = size(this%all_pairpots)
   ALLOCATE(v(npairpots))
   ALLOCATE(dvdzr(npairpots))
   ALLOCATE(dummy(npairpots))
   DO i = 1, nsites ! loop over sites
      A(1) = this%all_sites(i)%x
      A(2) = this%all_sites(i)%y
      DO j = 1, this%all_sites(i)%n ! loop over pairs v,z
         A(3)=this%all_sites(i)%z(j)
         CALL this%GET_REPUL_CORRECTIONS(A,v,dvdzr,dummy,dummy)
         this%all_sites(i)%interz%f(j)=this%all_sites(i)%interz%f(j)-sum(v)
         IF (j.EQ.1) THEN
            this%all_sites(i)%dz1=this%all_sites(i)%dz1-sum(dvdzr) ! correct first derivative
         ELSE IF (j.EQ.this%all_sites(i)%n) THEN
            this%all_sites(i)%dz2=this%all_sites(i)%dz2-sum(dvdzr) ! correct first derivative
         END IF
      END DO
   END DO
   RETURN
END SUBROUTINE SMOOTH_PES_HLIF001_NS
!############################################################
!# SUBROUTINE: INTERPOL_Z_PES_HLIF001_NS ########################
!############################################################
SUBROUTINE INTERPOL_Z_PES_HLIF001_NS(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_NS),INTENT(INOUT) :: this
   ! Local variables
   INTEGER(KIND=4) :: nsites,npairpots
   REAL(KIND=8) :: dz1,dz2
   INTEGER(KIND=4) :: i ! counters
   ! Run secton------------------------
   nsites=size(this%all_sites)
   npairpots=size(this%all_pairpots)
   CALL this%EXTRACT_VASINT()
   DO i = 1, npairpots ! loop pairpots
      dz1=this%all_pairpots(i)%dz1
      dz2=this%all_pairpots(i)%dz2
      CALL this%all_pairpots(i)%interz%INTERPOL(dz1,1,dz2,1)
   END DO
   CALL this%SMOOTH()
   DO i = 1, nsites
      dz1=this%all_sites(i)%dz1
      dz2=this%all_sites(i)%dz2
      CALL this%all_sites(i)%interz%INTERPOL(dz1,1,dz2,1)
   END DO
   RETURN
END SUBROUTINE INTERPOL_Z_PES_HLIF001_NS
!##################################################################
!# SUBROUTINE: INTERACTION_AP #####################################
!##################################################################
SUBROUTINE INTERACTION_AP(A,P,pairpot,dampfunc,interac,dvdz_corr,dvdx_corr,dvdy_corr)
   IMPLICIT NONE
   ! I/O VAriables --------------------------------------------
   REAL(KIND=8),DIMENSION(3),INTENT(IN):: A, P
   TYPE(Pair_pot),INTENT(IN):: pairpot
   type(One_func),INTENT(IN):: dampfunc
   REAL(KIND=8),INTENT(OUT):: interac,dvdz_corr,dvdx_corr,dvdy_corr
   ! Local variables ------------------------------------------
   REAL(KIND=8):: r ! distance
   REAL(KIND=8):: v,pre
   REAL(KIND=8):: aux ! dv/dr
   CHARACTER(LEN=*),PARAMETER:: routinename = "INTERACTION_AP: "
   ! GABBA, GABBA HEY! ----------------------------------------
   ! Find the distance between A and P, in a.u.
   r=dsqrt((A(1)-P(1))**2.D0+(A(2)-P(2))**2.D0+(A(3)-P(3))**2.D0)
   CALL pairpot%GET_V_AND_DERIVS(r,v,aux)
   interac=v*dampfunc%getvalue(r)
   pre=aux*dampfunc%getvalue(r)+v*dampfunc%getderiv(r)
   dvdz_corr=pre*(A(3)-P(3))/r
   dvdx_corr=pre*(A(1)-P(1))/r
   dvdy_corr=pre*(A(2)-P(2))/r
   RETURN
END SUBROUTINE INTERACTION_AP
!##################################################################
!# SUBROUTINE: INTERACTION_AENV ###################################
!##################################################################
SUBROUTINE INTERACTION_AENV(n,A,pairpot,dampfunc,interac,dvdz_term,dvdx_term,dvdy_term)
   IMPLICIT NONE
   ! I/O VAriables ------------------------------------------
   INTEGER,INTENT(IN) :: n
   REAL(KIND=8),DIMENSION(3),INTENT(IN) :: A
   TYPE(Pair_pot),INTENT(IN) :: pairpot
   type(One_func),INTENT(IN) :: dampfunc
   REAL(KIND=8),INTENT(OUT) :: interac, dvdz_term, dvdx_term, dvdy_term
   ! Local variables ----------------------------------------
   REAL(KIND=8),DIMENSION(3) :: P
   REAL(KIND=8),DIMENSION(3) :: ghost_A ! A in cartesians, but inside unitcell
   REAL(KIND=8),DIMENSION(2) :: aux
   REAL(KIND=8) :: dummy1, dummy2, dummy3, dummy4
   REAL(KIND=8) :: atomx, atomy
   INTEGER :: pairid
   INTEGER :: i, k ! Counters
   CHARACTER(LEN=18), PARAMETER :: routinename = "INTERACTION_AENV: "
   ! SUSY IS A HEADBANGER !!!! -------------------
   ! Defining some aliases to make the program simpler:
   pairid = pairpot%id
   interac=0.D0
   dvdz_term=0.D0
   dvdx_term=0.D0
   dvdy_term=0.D0
   ! ghost A definition
   ghost_A(1:2)=sysLiF001Surf%project_unitcell(A(1:2))
   ghost_A(3)=A(3)

   SELECT CASE(n)
      CASE(0)
         DO i=1, sysLiF001Surf%atomtype(pairid)%n
            P(:)=sysLiF001Surf%atomtype(pairid)%atom(i,:)
            CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac=interac+dummy1
            dvdz_term=dvdz_term+dummy2
            dvdx_term=dvdx_term+dummy3
            dvdy_term=dvdy_term+dummy4
         END DO
         RETURN

      CASE(1 :)
         DO i=1, sysLiF001Surf%atomtype(pairid)%n
            atomx=sysLiF001Surf%atomtype(pairid)%atom(i,1)
            atomy=sysLiF001Surf%atomtype(pairid)%atom(i,2)
            P(3)=sysLiF001Surf%atomtype(pairid)%atom(i,3)
            DO k= -n,n
               aux(1)=dfloat(n)
               aux(2)=dfloat(k)
               aux=sysLiF001Surf%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac=interac+dummy1
               dvdz_term=dvdz_term+dummy2
               dvdx_term=dvdx_term+dummy3
               dvdy_term=dvdy_term+dummy4
               !
               aux(1)=dfloat(-n)
               aux(2)=dfloat(k)
               aux=sysLiF001Surf%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
            END DO
            DO k= -n+1, n-1
               aux(1)=dfloat(k)
               aux(2)=dfloat(n)
               aux=sysLiF001Surf%surf2cart(aux)
               P(1) =atomx+aux(1)
               P(2) =atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
               !
               aux(1)=dfloat(k)
               aux(2)=dfloat(-n)
               aux=sysLiF001Surf%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
            END DO
         END DO
         RETURN

      CASE DEFAULT
         WRITE(0,*) "INTERACTION_AENV ERR: Wrong environment order."
         CALL EXIT(1)
   END SELECT
END SUBROUTINE INTERACTION_AENV
!############################################################
!# SUBROUTINE: GET_V_AND_DERIVS_PES_HLIF001_NS ##################
!############################################################
!> @brief
!! Subroutine that calculates the 3D potential for a point A and
!! its derivatives in cartesian coordinates.
!
!> @param[in] this - PES_HLIF001_NS PES
!> @param[in] X - Point in space to calculate the potential and it's derivatives. Cartesian's
!> @param[out] v - Value of the potential at X
!> @param[out] dvdu - derivatives, cartesian coordinates
!
!> @warning
!! - All sitios should have been interpolated in Z direction previously. It
!!   is optinal to extract "vasint".
!! - Corrugation, MUST have been extracted from sitios previously (smoothed)
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 06/Feb/2014
!> @version 1.0
!------------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_PES_HLIF001_NS(this,X,v,dvdu,errCode)
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_NS),TARGET,INTENT(IN):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN):: X
   REAL(KIND=8),INTENT(OUT):: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT):: dvdu
   integer(kind=1),optional,intent(out):: errCode
   ! Local variables
   INTEGER(KIND=4):: nsites,npairpots
   type(Fourierp4mm):: interpolxy
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: pot,dvdz,dvdx,dvdy
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: potarr
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: f,derivarr ! arguments to the xy interpolation
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: xy ! arguments to the xy interpolation
   INTEGER :: i ! counters
   ! Pointers
	REAL(KIND=8),POINTER:: zmax
   CHARACTER(LEN=*),PARAMETER:: routinename="GET_V_AND_DERIVS_PES_HLIF001_NS: "
   zmax => this%all_sites(1)%z(this%all_sites(1)%n)
   npairpots = size(this%all_pairpots)
   nsites = size(this%all_sites)
   ! GABBA, GABBA HEY! ----------------------
   SELECT CASE(X(3)>zmax)
      CASE(.TRUE.)
         dvdu = (/0.D0,0.D0,0.D0/)
         v = 0.D0
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   ! Initialization section
   ALLOCATE(pot(npairpots))
   ALLOCATE(dvdz(npairpots))
   ALLOCATE(dvdx(npairpots))
   ALLOCATE(dvdy(npairpots))
   CALL this%GET_REPUL_CORRECTIONS(X,pot,dvdz,dvdx,dvdy)
   v=sum(pot)
   dvdu(1)=sum(dvdx)
   dvdu(2)=sum(dvdy)
   dvdu(3)=sum(dvdz)
   ! Now, we have all the repulsive interaction and corrections to the derivarives
   ! stored in v(:) and dvdu(:) respectively.
   ! Let's get v and derivatives from xy interpolation of the corrugationless function
   ! f(1,i) ==> v values interpolated for Z at site "i"
   ! f(2,i) ==> dvdz values interpolated for Z at site "i"
   ALLOCATE(f(2,nsites))
   ALLOCATE(xy(nsites,2))
   ! potarr(1) ==> v interpolated at X,Y for a given Z
   ! potarr(2) ==> dvdz interpolated at X,Y for a given Z
   ALLOCATE(potarr(2))
   ! derivarr(1,1:2) ==> dvdx, dvdy interpolated at X,Y
   ! derivarr(2,1:2) ==> d(dvdz)/dx, d(dvdz)/dy interpolated at X,Y. This is not needed
   ALLOCATE(derivarr(2,2))
   DO i=1,nsites
      xy(i,1)=this%all_sites(i)%x
      xy(i,2)=this%all_sites(i)%y
      CALL this%all_sites(i)%interz%GET_V_AND_DERIVS(X(3),f(1,i),f(2,i))
   END DO
   CALL interpolxy%READ(xy,f,this%klist)
   CALL interpolxy%INTERPOL(sysLiF001Surf)
   CALL interpolxy%GET_F_AND_DERIVS(sysLiF001Surf,X,potarr,derivarr)
   ! Corrections from the smoothing procedure
   v = v + potarr(1)
   dvdu(1)=dvdu(1)+derivarr(1,1)
   dvdu(2)=dvdu(2)+derivarr(1,2)
   dvdu(3)=dvdu(3)+potarr(2)
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_PES_HLIF001_NS
!############################################################
!# SUBROUTINE: GET_V_AND_DERIVS_CORRECTION_PES_HLIF001_NS ############
!############################################################
!> @brief
!! Subroutine that calculates the correction to the 3D PES for a point A and
!! its derivatives in cartesian coordinates.
!
!> @param[in] this - PES_HLIF001_NS PES
!> @param[in] X - Point in space to calculate the potential and it's derivatives. Cartesian's
!> @param[out] v - Value of the potential at X
!> @param[out] dvdu - derivatives, cartesian coordinates
!
!> @warning
!! - All sitios should have been interpolated in Z direction previously. It
!!   is optinal to extract "vasint".
!! - Corrugation, MUST have been extracted from sitios previously (smoothed)
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 06/Feb/2014
!> @version 1.0
!------------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_CORRECTION_PES_HLIF001_NS(this,X,v,dvdu)
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_NS),TARGET,INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(3), INTENT(IN) :: X
   REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(3),INTENT(OUT) :: dvdu
   ! Local variables
   INTEGER(KIND=4) :: npairpots
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: pot,dvdz,dvdx,dvdy
   INTEGER :: i ! counters
   ! Pointers
	REAL(KIND=8), POINTER :: zmax
   CHARACTER(LEN=24),PARAMETER :: routinename="GET_V_AND_DERIVS_PES_HLIF001_NS: "
   zmax => this%all_sites(1)%interz%x(this%all_sites(1)%n)
   npairpots = size(this%all_pairpots)
   SELECT CASE(size(v)==npairpots+1)
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         WRITE(*,*) "GET_V_AND_DERIVS_CORRECTION_PES_HLIF001_NS ERR: wrong number of dimensions array v"
         CALL EXIT(1)
   END SELECT
   ! GABBA, GABBA HEY! ----------------------
   SELECT CASE(X(3)>zmax)
      CASE(.TRUE.)
         dvdu = (/0.D0,0.D0,0.D0/)
         FORALL(i=1:npairpots+1) v(i) = 0.D0
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   ! Initializing variables
   ALLOCATE(pot(npairpots))
   ALLOCATE(dvdz(npairpots))
   ALLOCATE(dvdx(npairpots))
   ALLOCATE(dvdy(npairpots))
   FORALL(i=1:npairpots+1) v(i) = 0.D0
   ! Compute
   CALL this%GET_REPUL_CORRECTIONS(X,pot,dvdz,dvdx,dvdy)
   v(1)=sum(pot)
   FORALL(i=2:npairpots+1) v(i)=pot(i-1)
   dvdu(1)=sum(dvdx)
   dvdu(2)=sum(dvdy)
   dvdu(3)=sum(dvdz)
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_CORRECTION_PES_HLIF001_NS
!############################################################
!# FUNCTION: getpot_crp3d ###################################
!############################################################
!> @brief
!! Subroutine that calculates the 3D potential for a point A
!
!> @param[in] this - PES_HLIF001_NS PES
!> @param[in] X - Point in space to calculate the potential. Cartesian's
!
!> @warning
!! - All sitios should have been interpolated in Z direction previously. It
!!   is optinal to extract "vasint".
!! - Corrugation, MUST have been extracted from sitios previously (smoothed)
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 06/Feb/2014
!> @version 1.0
!------------------------------------------------------------
function getpot_PES_HLIF001_NS(this,X) result(finalPot)
   implicit none
   ! I/O variables
   class(PES_HLIF001_NS),target,intent(in) :: this
	real(kind=8),dimension(3), intent(in) :: X
	! Function dummy variable
	real(kind=8):: finalPot
   ! Local variables
   class(Fourierp4mm),allocatable:: interpolxy
   real(kind=8):: v
   integer(kind=4) :: nsites,npairpots
   real(kind=8),dimension(:),allocatable :: pot,dummy
   real(kind=8),dimension(:),allocatable :: potarr
   real(kind=8),dimension(:,:),allocatable :: f,derivarr ! arguments to the xy interpolation
   real(kind=8),dimension(:,:),allocatable :: xy ! arguments to the xy interpolation
   INTEGER :: i ! counters
   ! Pointers
   real(kind=8), pointer :: zmax
   character(len=24),parameter :: routinename="GET_V_AND_DERIVS_PES_HLIF001_NS: "
   zmax => this%all_sites(1)%interz%x(this%all_sites(1)%n)
   npairpots = size(this%all_pairpots)
   nsites = size(this%all_sites)
   ! GABBA, GABBA HEY! ----------------------
   allocate(pot(npairpots))
   allocate(dummy(npairpots))
   select case(x(3)>zmax)
      case(.true.)
         finalPot=0.D0
         return
      case(.false.)
         ! do nothing
   end select
   !
   call this%GET_REPUL_CORRECTIONS(X,pot,dummy,dummy,dummy)
   v=sum(pot)
   ! now, we have all the repulsive interaction and corrections to the derivarives
   ! stored in v(:) and dvdu(:) respectively.
   ! let's get v and derivatives from xy interpolation of the corrugationless function
   allocate(f(2,nsites))
   allocate(xy(nsites,2))
   allocate(potarr(2))
   allocate(derivarr(2,2))
   do i=1,nsites
      xy(i,1)=this%all_sites(i)%x
      xy(i,2)=this%all_sites(i)%y
      call this%all_sites(i)%interz%GET_V_AND_DERIVS(X(3),f(1,i),f(2,i))
   end do
   call interpolxy%READ(xy,f,this%klist)
   call interpolxy%INTERPOL(sysLiF001Surf)
   call interpolxy%GET_F_AND_DERIVS(sysLiF001Surf,X,potarr,derivarr)
   ! corrections from the smoothing procedure
   finalPot=v+potarr(1)
   return
End function getpot_PES_HLIF001_NS
!######################################################################
! SUBROUTINE: PLOT_DATA_SYMMPOINT #####################################
!######################################################################
SUBROUTINE PLOT_DATA_SYMMPOINT(this,filename)
   IMPLICIT NONE
   ! I/O Variables ----------------------
   CLASS(Symmpoint), INTENT(IN) :: this
   CHARACTER(LEN=*), INTENT(IN) :: filename
   ! Local variables --------------------
   INTEGER :: i ! Counter
   ! GABBA GABBA HEY!! ------------------
   OPEN(11, FILE=filename, STATUS="replace")
   DO i=1,this%n
      WRITE(11,*) this%z(i),this%v(i)
   END DO
   CLOSE (11)
   WRITE(*,*) "PLOT_DATA_SYMMPOINT: ",this%alias,filename," file created"
END SUBROUTINE PLOT_DATA_SYMMPOINT
!#######################################################################
! SUBROUTINE: PLOT_XYMAP_PES_HLIF001_NS
!#######################################################################
SUBROUTINE PLOT_XYMAP_PES_HLIF001_NS(this,filename,init_xyz,nxpoints,nypoints,Lx,Ly)
   IMPLICIT NONE
   CLASS(PES_HLIF001_NS),INTENT(IN) :: this
   REAL*8,DIMENSION(3),INTENT(IN) :: init_xyz ! Initial position to start the scan (in a.u.)
   INTEGER,INTENT(IN) :: nxpoints, nypoints ! number of points in XY plane
   CHARACTER(LEN=*),INTENT(IN) :: filename ! filename
   REAL*8,INTENT(IN) :: Lx ! Length of X axis
   REAL*8,INTENT(IN) :: Ly ! Length of X axis
   ! Local variables
   REAL*8 :: xmin, ymin, xmax, ymax, z
   REAL*8, DIMENSION(3) :: r, dvdu
   REAL*8 :: xdelta, ydelta
   INTEGER :: xinpoints, nxdelta
   INTEGER :: yinpoints, nydelta
   INTEGER :: i, j ! counters
   REAL*8 :: v ! potential
	! GABBA, GABBA HEY! ---------
   xmin = init_xyz(1)
   ymin = init_xyz(2)
   z = init_xyz(3)
   xmax = init_xyz(1)+Lx
   ymax = init_xyz(2)+Ly
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
   r(3) = z
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v
   END DO
   r(2) = ymax
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v
   ! inpoints in XY
   DO i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      r(3) = z
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v
      DO j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         CALL this%GET_V_AND_DERIVS(r,v,dvdu)
         WRITE(11,*) r(1), r(2), v
      END DO
      r(2) = ymax
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v
   END DO
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   r(3) = z
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v
   END DO
   r(2) = ymax
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_XYMAP_PES_HLIF001_NS
!#######################################################################
! SUBROUTINE: PLOT_DIRECTION1D_PES_HLIF001_NS ###################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. To define
!! the direction, the angle alpha is given.
!
!> @param[in] this - PES_HLIF001_NS PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] angle - Angle between the surface vector S1 and the direction of the
!!                    cut. It should be given in degrees.
!> @param[in] z - Distance to the surface. All points have the same z.
!> @param[in] L - Length of the graphic
!
!> @warning
!! - The graph starts always at 0,0
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 09/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT_DIRECTION1D_PES_HLIF001_NS(this,filename,npoints,angle,z,L)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_HLIF001_NS),INTENT(IN) :: this
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*), INTENT(IN) :: filename
   REAL*8, INTENT(IN) :: z, angle
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL*8 :: delta,L,v,s,alpha
   REAL*8 :: xmax, xmin, ymax, ymin
   REAL*8, DIMENSION(3) :: r, dvdu
   INTEGER :: i ! Counter
   CHARACTER(LEN=24), PARAMETER :: routinename = "PLOT_DIRECTION1D_PES_HLIF001_NS: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT_DIRECTION1D_PES_HLIF001_NS ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   ! Change alpha to radians
   alpha = angle * PI / 180.D0
   !
   xmin = 0.D0
   ymin = 0.D0
   xmax = L*DCOS(alpha)
   ymax = L*DSIN(alpha)
   r(3) = z
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=L/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(1) = xmin
   r(2) = ymin
   s = DSQRT(r(1)**2.D0+r(2)**2.D0)
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) s, v , dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   ! cycle for inpoints
   DO i=1, inpoints
      r(1)=xmin+(DFLOAT(i)*delta)*DCOS(alpha)
      r(2)=ymin+(DFLOAT(i)*delta)*DSIN(alpha)
      s = DSQRT(r(1)**2.D0+r(2)**2.D0)
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) s, v , dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   END DO
   ! Final value
   r(1) = xmax
   r(2) = ymax
   s = DSQRT(r(1)**2.D0+r(2)**2.D0)
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) s, v, dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT_DIRECTION1D_PES_HLIF001_NS
!#######################################################################
!# SUBROUTINE: PLOT_Z_PES_HLIF001_NS #######################################
!#######################################################################
SUBROUTINE PLOT_Z_PES_HLIF001_NS(this,npoints,xyz,L,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_HLIF001_NS),INTENT(IN):: this
   INTEGER,INTENT(IN):: npoints
   CHARACTER(LEN=*),INTENT(IN):: filename
   REAL(KIND=8),DIMENSION(3),INTENT(IN):: xyz
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8):: delta,L
   REAL(KIND=8):: zmax, zmin
   REAL(KIND=8),DIMENSION(3):: r, dvdu
   REAL(KIND=8):: v
   INTEGER:: i ! Counter
   CHARACTER(LEN=*),PARAMETER:: routinename = "PLOT_DIRECTION1D_PES_HLIF001_NS: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT_Z_PES_HLIF001_NS ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   r(1)=xyz(1)
   r(2)=xyz(2)
   zmin=xyz(3)
   zmax=xyz(3)+L
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=L/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(3) = zmin
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)
   ! cycle for inpoints
   DO i=1, inpoints
      r(3)=zmin+(DFLOAT(i)*delta)
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(3),v,dvdu(:)
   END DO
   ! Final value
   r(3) = zmax
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT_Z_PES_HLIF001_NS
!###########################################################
!# SUBROUTINE: PLOT_INTERPOL_SYMMPOINT
!###########################################################
SUBROUTINE PLOT_INTERPOL_SYMMPOINT(this,npoints,filename)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Symmpoint),INTENT(IN):: this
   INTEGER(KIND=4),INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Run section
   CALL this%interz%PLOT(npoints,filename)
   RETURN
END SUBROUTINE PLOT_INTERPOL_SYMMPOINT
!###########################################################
!# SUBROUTINE: PLOT_PAIRPOTS_PES_HLIF001_NS
!###########################################################
SUBROUTINE PLOT_PAIRPOTS_PES_HLIF001_NS(this,npoints)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_NS),INTENT(IN)::this
   INTEGER(KIND=4),INTENT(IN) :: npoints
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   INTEGER(KIND=4) :: npairpots
   CHARACTER(LEN=12) :: stringbase
   CHARACTER(LEN=13) :: filename
   ! Run section
   npairpots=size(this%all_pairpots)
   WRITE(stringbase,'(A12)') "-pairpot.dat"
   DO i = 1, npairpots
      WRITE(filename,'(I1,A12)') i,stringbase
      CALL this%all_pairpots(i)%PLOT(npoints,filename)
   END DO
   RETURN
END SUBROUTINE PLOT_PAIRPOTS_PES_HLIF001_NS
!###########################################################
!# SUBROUTINE: PLOT_SITIOS_PES_HLIF001_NS
!###########################################################
SUBROUTINE PLOT_SITIOS_PES_HLIF001_NS(this,npoints)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_NS),INTENT(IN)::this
   INTEGER(KIND=4),INTENT(IN) :: npoints
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   INTEGER(KIND=4) :: nsitios
   CHARACTER(LEN=10) :: stringbase
   CHARACTER(LEN=11) :: filename
   ! Run section
   nsitios=size(this%all_sites)
   WRITE(stringbase,'(A10)') "-sitio.dat"
   DO i = 1, nsitios
      WRITE(filename,'(I1,A10)') i,stringbase
      CALL this%all_sites(i)%PLOT(npoints,filename)
   END DO
   RETURN
END SUBROUTINE PLOT_SITIOS_PES_HLIF001_NS
!###########################################################
!# FUNCTION: is_allowed_PES_HLIF001_NS
!###########################################################
LOGICAL FUNCTION is_allowed_PES_HLIF001_NS(this,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_NS),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8) :: xmin,xmax
   ! Run section
   xmin=this%all_sites(1)%z(1)
   xmax=this%all_sites(1)%z(this%all_sites(1)%n)
   SELECT CASE(size(x)/=3)
      CASE(.TRUE.)
         WRITE(0,*) "is_allowed_PES_HLIF001_NS ERR: array doesn't have 3 dimensions: 3"
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE( x(3)<xmin )
      CASE(.TRUE.)
         is_allowed_PES_HLIF001_NS=.FALSE.
      CASE(.FALSE.)
         is_allowed_PES_HLIF001_NS=.TRUE.
   END SELECT
   RETURN
END FUNCTION is_allowed_PES_HLIF001_NS
END MODULE PES_HLIF001_NS_MOD

!#########################################################
! MODULE: XEXPONENTIAL_FUNCTION
!> @brief
!! Compendium of routines and types to manage a xexponential function
!##########################################################
MODULE XEXPONENTIAL_FUNCTION_MOD
USE FUNCTION1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: xexponential_func
!> @brief
!! Function @f$ f(x)=axe^{bx} @f$, where @b a and @b b are
!! parameters
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Function1d) :: xexponential_func
CONTAINS
   PROCEDURE,PUBLIC :: getvalue => getvalue_xexponential
   PROCEDURE,PUBLIC :: getderiv => getderiv_xexponential
END TYPE xexponential_func
!/////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# FUNCTION: getvalue_xexponential
!###########################################################
!> @brief
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getvalue_xexponential(this,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Xexponential_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: a
   ! Run section
   a=this%getparam()
   getvalue_xexponential=a(1)*x*dexp(a(2)*x)
   RETURN
END FUNCTION getvalue_xexponential
!###########################################################
!# FUNCTION: getderiv_xexponential
!###########################################################
!> @brief
!! Computes value of the function
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_xexponential(this,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(xexponential_func),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: a
   ! Run section
   a=this%getparam()
   getderiv_xexponential=a(1)*(1.D0+a(2)*x)*dexp(a(2)*x)
   RETURN
END FUNCTION getderiv_xexponential

END MODULE XEXPONENTIAL_FUNCTION_MOD! contains, body

!#########################################################
! MODULE: EXTRAPOL_TO_VACUUM_MOD
!> @brief
!! Provides tools to extrapolate the potential in the vacuum,
!! which is only dependent on the distance between both atoms
!##########################################################
MODULE EXTRAPOL_TO_VACUUM_MOD
! Initial declarations
use CUBICSPLINES_MOD, only: Csplines
use MATHS_MOD, only: ORDER
use UNITS_MOD, only: Energy,Length
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: Vacuumpot
!> @brief
!! Class that stores all information that can be extracted from the
!! monodimensional vacuum potential
!----------------------------------------------------------------
TYPE :: Vacuumpot
   PRIVATE
   INTEGER(KIND=4),PUBLIC:: n
   TYPE(Csplines),PUBLIC:: rpot
   REAL(KIND=8):: surfen
   REAL(KIND=8):: potmin
   REAL(KIND=8):: req
   REAL(KIND=8):: forceConstant
   REAL(KIND=8),DIMENSION(2),PUBLIC:: root

   CONTAINS
      ! Initialization block
      PROCEDURE,PUBLIC:: INITIALIZE => INITIALIZE_VACUUMPOT
      PROCEDURE,PUBLIC:: INITIALIZE_DIRECT => INITIALIZE_DIRECT_VACUUMPLOT
      PROCEDURE,PUBLIC:: READ => READ_VACUUMPOT
      ! Set block
      PROCEDURE,PUBLIC:: SET_ROOTS => SET_ROOTS_VACUUMPOT
      ! Get block
      PROCEDURE,PUBLIC:: getpot => getpot_vacuumpot
      PROCEDURE,PUBLIC:: getderiv => getderiv_vacuumpot
      PROCEDURE,PUBLIC:: getscalefactor => getscalefactor_vacuumpot
      PROCEDURE,PUBLIC:: getreq => getreq_vacuumpot
      PROCEDURE,PUBLIC:: getpotmin => getpotmin_vacuumpot
      PROCEDURE,PUBLIC:: getForceConstant => getForceConstant_vacuumpot
      ! Tools block
      PROCEDURE,PUBLIC:: SHIFTPOT => SHIFTPOT_VACUUMPOT
      PROCEDURE,PUBLIC:: SHIFTPOT_UNDO => SHIFTPOT_UNDO_VACUUMPOT
      ! Enquire block
      PROCEDURE,PUBLIC:: is_allowed => is_allowed_VACUUMPOT
      ! Plot tools block
      PROCEDURE,PUBLIC:: PLOT => PLOT_VACUUMPOT
      PROCEDURE,PUBLIC:: PLOTDATA => PLOTDATA_VACUUMPOT
END TYPE Vacuumpot
!/////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_ROOTS_VACUUMPOT
!###########################################################
!> @brief
!! Calculates roots of the function if any
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_ROOTS_VACUUMPOT(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(INOUT):: this
   ! Local variables
   CHARACTER(LEN=21) :: routinename="SET_ROOTS_VACUUMPOT: "
   ! Run section
   CALL this%rpot%SET_XROOT()
   SELECT CASE(allocated(this%rpot%xroot))
      CASE(.TRUE.)
         SELECT CASE(size(this%rpot%xroot))
            CASE(1)
               WRITE(0,*) "SET_ROOTS_VACUUMPOT ERR: we are in req. Bad. Pure classic not implemented"
               CALL EXIT(1)
            CASE(2)
               this%root=this%rpot%xroot
            CASE DEFAULT
               WRITE(0,*) "SET_ROOTS_VACUUMPOT ERR: too many roots. Check the vacuumpot"
               CALl EXIT(1)
         END SELECT
         ! body
      CASE(.FALSE.)
         WRITE(0,*) "SET_ROOTS_VACUUMPOT ERR: there aren't roots"
         CALL EXIT(1)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE SET_ROOTS_VACUUMPOT
!###########################################################
!# FUNCTION: getpotmin_VACUUMPOT
!###########################################################
!> @brief
!! Standard get function. Gets potential at the minimum prior
!! shift.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jan/2015
!> @version 1.0
!-----------------------------------------------------------
PURE FUNCTION getpotmin_vacuumpot(this) result(potmin)
    ! Initial declarations
    IMPLICIT NONE
    ! I/O variables
    CLASS(Vacuumpot),INTENT(IN):: this
    ! Dummy function variable
    REAL(KIND=8):: potmin
    ! Run section
    potmin=this%potmin
    RETURN
END FUNCTION
!###########################################################
!# FUNCTION: getForceConstant_VACUUMPOT
!###########################################################
!> @brief
!! Standard get function. Gets forceConstant atribute.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jan/2015
!> @version 1.0
!-----------------------------------------------------------
PURE FUNCTION getForceConstant_vacuumpot(this) result(forceConstant)
    ! Initial declarations
    IMPLICIT NONE
    ! I/O variables
    CLASS(Vacuumpot),INTENT(IN):: this
    ! Dummy function variable
    REAL(KIND=8):: forceConstant
    ! Run section
    forceConstant=this%forceConstant
    RETURN
END FUNCTION
!###########################################################
!# FUNCTION: is_allowed_VACUUMPOT
!###########################################################
!> @brief
!! Checks if the potential can be calculated in @b X
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
LOGICAL FUNCTION is_allowed_VACUUMPOT(this,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Local variables
   ! Run section
   SELECT CASE(x< this%rpot%x(1) .OR. x > this%rpot%x(this%n))
      CASE(.TRUE.)
         is_allowed_VACUUMPOT=.FALSE.
      CASE(.FALSE.)
         is_allowed_VACUUMPOT=.TRUE.
   END SELECT
   RETURN
END FUNCTION is_allowed_VACUUMPOT
!###########################################################
!# SUBROUTINE: READ_DIRECT_VACUUMPLOT
!###########################################################
!> @brief
!! Reads information directly from array
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INITIALIZE_DIRECT_VACUUMPLOT(this,x,f,surfaceEnergy)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(OUT):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: f
   REAL(KIND=8),INTENT(IN):: surfaceEnergy
   ! Local variables
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: aux1,aux2
   CHARACTER(LEN=29),PARAMETER :: routinename="INITIALIZE_DIRECT_VACUUMPOT: "
   ! Run section
   this%n=size(x)
   ALLOCATE(aux1(size(x)))
   ALLOCATE(aux2(size(f)))
   aux1=x
   aux2=f
   this%surfEn=surfaceEnergy
   CALL ORDER(aux1,aux2)
   CALL this%rpot%READ(aux1,aux2)
   CALL this%rpot%INTERPOL(0.D0,0,0.D0,0)
   CALL this%rpot%SET_MINIMUM()
   SELECT CASE(ALLOCATED(this%rpot%xmin))
      CASE(.TRUE.)
         SELECT CASE(size(this%rpot%xmin))
            CASE(1)
               this%req=this%rpot%xmin(1)
               this%potmin = this%rpot%getvalue(this%rpot%xmin(1))
               this%forceConstant = this%rpot%getderiv2(this%rpot%xmin(1))
            CASE DEFAULT
               WRITE(0,*) "INITIALIZE_DIRECT_VACUUMPOT ERR: More than one minimum. Something is wrong"
               CALL EXIT(1)
         END SELECT
      CASE(.FALSE.)
        WRITE(0,*) "INITIALIZE_DIRECT_VACUUMPOT ERR: There's not any minimum. Something is wrong"
        CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE INITIALIZE_DIRECT_VACUUMPLOT
!###########################################################
!# FUNCTION: getreq_vacuumpot
!###########################################################
!> @brief
!! Common get function. Gets "req" atribute
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getreq_vacuumpot(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(IN):: this
   ! Run section
   getreq_vacuumpot=this%req
   RETURN
END FUNCTION getreq_vacuumpot
!###########################################################
!# SUBROUTINE: SHIFTPOT_VACUUMPOT
!###########################################################
!> @brief
!! Shift the entire potential so that the equilibrium geomtry has
!! energy 0
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 27/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SHIFTPOT_VACUUMPOT(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(INOUT):: this
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   ! Run section
   DO i = 1, this%n
      this%rpot%f(i)=this%rpot%f(i)-this%potmin
   END DO
   CALL this%rpot%REINTERPOL(0.D0,0,0.D0,0)
   RETURN
END SUBROUTINE SHIFTPOT_VACUUMPOT
!###########################################################
!# SUBROUTINE: SHIFTPOT_UNDO_VACUUMPOT
!###########################################################
!> @brief
!! Shift the entire potential so that the equilibrium geomtry has
!! energy 0
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 27/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SHIFTPOT_UNDO_VACUUMPOT(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(INOUT):: this
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   ! Run section
   DO i = 1, this%n
      this%rpot%f(i)=this%rpot%f(i)+this%potmin
   END DO
   CALL this%rpot%REINTERPOL(0.D0,0,0.D0,0)
   RETURN
END SUBROUTINE SHIFTPOT_UNDO_VACUUMPOT
!###########################################################
!# FUNCTION: getscalefactor_vacuumpot
!###########################################################
!> @brief
!! Common get function. Gets the sum of atributes surfen and potmin
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getscalefactor_vacuumpot(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(IN):: this
   ! Run section
   getscalefactor_vacuumpot=this%surfen+this%potmin
   RETURN
END FUNCTION getscalefactor_vacuumpot
!###########################################################
!# FUNCTION: getderiv_vacuumpot
!###########################################################
!> @brief
!! Common get function. Gets value of the derivative of the potential at R
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_vacuumpot(this,r)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(IN)::this
   REAL(KIND=8),INTENT(IN) :: r
   ! Run section
   getderiv_vacuumpot=this%rpot%getderiv(r)
   RETURN
END FUNCTION getderiv_vacuumpot
!###########################################################
!# FUNCTION: getpot_vacuumpot
!###########################################################
!> @brief
!! Common get function. Gets value of the potential at R
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getpot_vacuumpot(this,r)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(IN)::this
   REAL(KIND=8),INTENT(IN) :: r
   ! Run section
   getpot_vacuumpot=this%rpot%getvalue(r)
   RETURN
END FUNCTION getpot_vacuumpot
!###########################################################
!# SUBROUTINE: PLOTDATA_VACUUMPOT
!###########################################################
!> @brief
!! Plots a graph with values of the potential at grid points
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 26/Mar/2013
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE PLOTDATA_VACUUMPOT(this,filename)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(IN)::this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Run section
   CALL this%rpot%PLOTDATA(filename)
   RETURN
END SUBROUTINE PLOTDATA_VACUUMPOT
!###########################################################
!# SUBROUTINE: PLOT_VACUUMPOT
!###########################################################
!> @brief
!! Plots a graph once interpolation was done
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE PLOT_VACUUMPOT(this,npoints,filename)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(IN)::this
   INTEGER(KIND=4),INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Run section
   CALL this%rpot%PLOT(npoints,filename)
   RETURN
END SUBROUTINE PLOT_VACUUMPOT
!###########################################################
!# SUBROUTINE: INITIALIZE_VACUUMPOT
!###########################################################
!> @brief
!! Reads data from file, interpolates data, finds minimum
!! energy point
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 26/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INITIALIZE_VACUUMPOT(this,filename)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(INOUT) :: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   CHARACTER(LEN=22),PARAMETER :: routinename="INITIALIZE_VACUUMPOT: "
   ! Run section
   CALL this%READ(filename)
   CALL this%rpot%INTERPOL(0.D0,0,0.D0,0)
   CALL this%rpot%SET_MINIMUM()
   SELECT CASE(ALLOCATED(this%rpot%xmin))
      CASE(.TRUE.)
         SELECT CASE(size(this%rpot%xmin))
            CASE(1)
               this%req=this%rpot%xmin(1)
               this%potmin = this%rpot%getvalue(this%rpot%xmin(1))
               this%forceConstant = this%rpot%getderiv2(this%rpot%xmin(1))
            CASE DEFAULT
               WRITE(0,*) "INITIALIZE_VACUUMPOT ERR: More than one minimum. Something is wrong"
               CALL EXIT(1)
         END SELECT
      CASE(.FALSE.)
        WRITE(0,*) "INITIALIZE_VACUUMPOT ERR: There's not any minimum. Something is wrong"
        CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE INITIALIZE_VACUUMPOT
!###########################################################
!# SUBROUTINE: READ_VACUUMPOT
!###########################################################
!> @brief
!! Reads data from file
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 26/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE READ_VACUUMPOT(this,filename)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(INOUT) :: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables
   CHARACTER(LEN=10) :: lenunits,enunits
   REAL(KIND=8) :: aux1,aux2
   TYPE(Energy) :: en
   TYPE(Length) :: len
   INTEGER(KIND=4) :: i ! counter
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: r,f
   ! Run section
   OPEN (111,FILE=filename,STATUS="old",ACTION="read")
   READ(111,*) ! dummy line
   READ(111,*) lenunits,enunits
   READ(111,*) aux1
   CALL en%READ(aux1,enunits)
   CALL en%TO_STD()
   this%surfen=en%getvalue()
   READ(111,*) this%n
   ALLOCATE(r(this%n))
   ALLOCATE(f(this%n))
   DO i = 1, this%n
      READ(111,*) aux1,aux2
      CALL len%READ(aux1,lenunits)
      CALL len%TO_STD()
      r(i)=len%getvalue()
      CALL en%READ(aux2,enunits)
      CALL en%TO_STD()
      f(i)=en%getvalue()
   END DO
   CALL ORDER(r,f) ! order R from low values to high values
   CALL this%rpot%READ(r,f)
   CLOSE(111)
   RETURN
END SUBROUTINE READ_VACUUMPOT
END MODULE EXTRAPOL_TO_VACUUM_MOD

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
use INTERPOL1D_MOD
use UNITS_MOD, only: pi
use MATHS_MOD, only: INV_MTRX
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: Termcalculator
!> @brief
!! Abstract class to calculate terms of the series avoiding the use of
!! unnecessary switches
!----------------------------------------------------------------
type,abstract:: Termcalculator
private
   real(kind=8):: shift=0.D0
   contains
      procedure,public,non_overridable:: getShift => getShift_TERMCALCULATOR
      procedure(getvalue_termcalculator_example),public,deferred:: getvalue
      procedure(getvalue_termcalculator_example),public,deferred:: getderiv
end type Termcalculator
!
ABSTRACT INTERFACE
   !###########################################################
   !# FUNCTION: getvalue_termcalculator_example
   !###########################################################
   !> @brief
   !! Just an example that child objects should override
   !-----------------------------------------------------------
   function getvalue_termcalculator_example(this,kpoint,parity,irrep,x) result(answer)
      import Termcalculator
      class(Termcalculator),intent(in):: this
      integer(kind=4),intent(in):: kpoint
      real(kind=8),intent(in):: x
      character(len=1),intent(in):: parity
      character(len=2),intent(in):: irrep
      real(kind=8):: answer
   end function getvalue_termcalculator_example
   !-------------------------------------------------------------
END INTERFACE
!/////////////////////////////////////////////////////////////////////////////
! TYPE: FOURIER1D
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!----------------------------------------------------------------------------
type,abstract,extends(Interpol1d):: Fourier1d
   ! public atributes
   integer(kind=4),dimension(:),allocatable,public:: kList
   character(len=1),dimension(:),allocatable,public:: parityList
   character(len=2),dimension(:),allocatable,public:: irrepList
   class(TermCalculator),allocatable,public:: term
   ! private atributes
   real(kind=8),dimension(:),allocatable,private:: coeff
   real(kind=8),dimension(:,:),allocatable,private:: extracoeff
   real(kind=8),dimension(:,:),allocatable,private:: extrafuncs
   contains
      ! initialize block
      procedure(initializeTerms_FOURIER1D),public,deferred:: initializeTerms
      ! get block
      procedure,public,non_overridable:: getValue => getvalue_FOURIER1D
      procedure,public,non_overridable:: getDeriv => getderiv_FOURIER1D
      procedure,public,non_overridable:: getKlist => getklist_FOURIER1D
      procedure,public,non_overridable:: getParityList => getParityList_FOURIER1D
      ! set block
      procedure,public,non_overridable:: setShift => set_shift_FOURIER1D
      procedure,public,non_overridable:: setKlist => setKlist_FOURIER1D
      procedure,public,non_overridable:: setParityList => setParityList_FOURIER1D
      procedure,public,non_overridable:: setIrrepList => setIrrepList_FOURIER1D
      ! tools
      procedure,public,non_overridable:: interpol => interpol_FOURIER1D
      procedure,public,non_overridable:: add_morefuncs => add_more_funcs_FOURIER1D
      procedure,public,non_overridable:: get_allfuncs_and_derivs => get_allfunc_and_derivs_FOURIER1D
      ! plotting tools
      procedure,public,non_overridable:: plotCyclic => plotcyclic_interpol_FOURIER1D
      procedure,public,non_overridable:: plotCyclic_all => plotcyclic_all_interpol_FOURIER1D
END TYPE FOURIER1D
ABSTRACT INTERFACE
   !###########################################################
   !# SUBROUTINE: SET_IRREP_FOURIER1D
   !###########################################################
   !> @brief
   !! Sets irrep for this fourier series. Should be overriden by
   !! child non-abstract classes
   !-----------------------------------------------------------
   subroutine initializeTerms_FOURIER1D(this)
      import Fourier1d
      class(fourier1d),intent(inout):: this
   end subroutine initializeTerms_FOURIER1D
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
!# FUNCTION: getParitylist_FOURIER1D
!###########################################################
!> @brief
!! Common get function. Gets parityList atribute
!-----------------------------------------------------------
function getParityList_FOURIER1D(this) result(charArray)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier1d),intent(in):: this
   ! Dummy output variable
   character(len=1),dimension(:),allocatable:: charArray
   ! Run section
   allocate( charArray(size(this%ParityList)),source=this%parityList(:) )
   return
end function getParityList_FOURIER1D
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
!###################################################################
!# SUBROUTINE: setKlist_FOURIER1D
!###################################################################
!> @brief
!! Common set subroutine. Sets Klist atribute of a FOURIER1D object
!-------------------------------------------------------------------
subroutine setKlist_FOURIER1D(this,kList)
   implicit none
   ! I/O variables
   class(Fourier1d),intent(inout):: this
   integer(kind=4),dimension(:):: kList
   ! Run section
   allocate( this%kList(size(kList)),source=kList(:) )
   return
end subroutine setKlist_FOURIER1D
!###################################################################
!# SUBROUTINE: setParityList_FOURIER1D
!###################################################################
!> @brief
!! Common set subroutine. Sets parityList atribute of a FOURIER1D object
!-------------------------------------------------------------------
subroutine setParityList_FOURIER1D(this,parityList)
   implicit none
   ! I/O variables
   class(Fourier1d),intent(inout):: this
   character(len=1),dimension(:):: parityList
   ! Run section
   allocate( this%parityList(size(parityList)),source=parityList(:) )
   return
end subroutine setParityList_FOURIER1D
!###################################################################
!# SUBROUTINE: setIrrepList_FOURIER1D
!###################################################################
!> @brief
!! Common set subroutine. Sets IrrepList atribute of a FOURIER1D object
!-------------------------------------------------------------------
subroutine setIrrepList_FOURIER1D(this,irrepList)
   implicit none
   ! I/O variables
   class(Fourier1d),intent(inout):: this
   character(len=2),dimension(:):: irrepList
   ! Run section
   allocate( this%irrepList(size(irrepList)),source=irrepList(:) )
   return
end subroutine setIrrepList_FOURIER1D
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
      terms(i)=this%term%getvalue( kpoint=this%kList(i),parity=this%parityList(i),irrep=this%irrepList(i),x=x )
      terms_dx(i)=this%term%getderiv( kpoint=this%kList(i),parity=this%parityList(i),irrep=this%irrepList(i),x=x )
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
         terms(i,j)=this%term%getvalue( kpoint=this%kList(j),parity=this%parityList(j),irrep=this%irrepList(j),x=this%x(i))
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
      terms(i)=this%term%getvalue( kpoint=this%kList(i),parity=this%parityList(i),irrep=this%irrepList(i),x=r)
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
      terms(i)=this%term%getvalue( kpoint=this%klist(i),parity=this%parityList(i),irrep=this%irrepList(i),x=r)
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
   CLOSE(11)
END SUBROUTINE PLOTCYCLIC_ALL_INTERPOL_FOURIER1D

END MODULE FOURIER1D_MOD

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
!> @param[in] f(:,:) - stores F falues for each point in the grid
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
   ! Run section
   nx=size(this%x)
   ny=size(this%y)
   OPEN (10,FILE=filename,STATUS="replace",ACTION="write")
   DO i = 1, nx
      DO j = 1, ny
         WRITE(10,*) this%x(i),this%y(j),this%fgrid(i,j)
      END DO
   END DO
   RETURN
END SUBROUTINE PLOTDATA_INTERPOLGRID2D
END MODULE INTERPOLGRID2D_MOD

!##########################################################
! MODULE: BICUBICSPLINES
!
!> @brief
!! Provides tools to perform bicubic splines interpolations on
!! 2D functions
!##########################################################
MODULE BICSPLINES_MOD
use INTERPOLGRID2D_MOD
use CUBICSPLINES_MOD
use UNITS_MOD, only: pi
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
TYPE,EXTENDS(Interpolgrid2d) :: Bicsplines
   PRIVATE
   TYPE(Csplines),DIMENSION(:),ALLOCATABLE :: xcsplines
   TYPE(Csplines),DIMENSION(:),ALLOCATABLE :: ycsplines
   REAL(KIND=8),DIMENSION(:,:,:,:),ALLOCATABLE :: coeff
   CONTAINS
      ! Get block
      PROCEDURE,PUBLIC :: getvalue => getvalue_bicsplines
      PROCEDURE,PUBLIC :: getderivx => getderivx_bicsplines
      PROCEDURE,PUBLIC :: getderivy => getderivy_bicsplines
      PROCEDURE,PUBLIC :: getderivxy => getderivxy_bicsplines
      ! Tools block
      PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_BICSPLINES
      PROCEDURE,PUBLIC :: INTERPOL_NEWGRID => INTERPOL_NEWGRID_BICSPLINES
      PROCEDURE,PUBLIC :: REBOOT => REBOOT_BICSPLINES
      ! Pot tools
      PROCEDURE,PUBLIC :: PLOT_XYMAP => PLOT_XYMAP_BICSPLINES
      PROCEDURE,PUBLIC :: PLOT_SPLINES => PLOT_SPLINES_BICSPLINES
      PROCEDURE,PUBLIC :: PLOT_1D => PLOT_1D_BICSPLINES
      PROCEDURE,PUBLIC :: PLOT_DUALDERIVS_AT_GRID
END TYPE
!//////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: REBOOT_BICSPLINES
!###########################################################
!> @brief
!! Deallocates all information inside a Bicsplines variable
!-----------------------------------------------------------
SUBROUTINE REBOOT_BICSPLINES(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Bicsplines),INTENT(INOUT) :: this
   ! Run section
   DEALLOCATE(this%xcsplines)
   DEALLOCATE(this%ycsplines)
   DEALLOCATE(this%x)
   DEALLOCATE(this%y)
   DEALLOCATE(this%fgrid)
   DEALLOCATE(this%coeff)
   RETURN
END SUBROUTINE REBOOT_BICSPLINES
!###########################################################
!# SUBROUTINE: INTERPOL_NEWGRID_BICSPLINES
!###########################################################
!> @brief
!! Use monodimensional cubic splines interpolation to generate new
!! grid.
!
!> @warnig
!! - Use after we have a correct interpolation of the initial grid
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 21/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INTERPOL_NEWGRID_BICSPLINES(this,nxpoints,nypoints)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Bicsplines),INTENT(INOUT):: this
   INTEGER,INTENT(IN) :: nxpoints,nypoints ! number of points in XY plane
   ! Local variables
   INTEGER(KIND=4) :: oldnx,oldny
   REAL*8 :: xmin, ymin, xmax, ymax
   REAL*8 :: xdelta, ydelta
   INTEGER :: xinpoints, nxdelta
   INTEGER :: yinpoints, nydelta
   INTEGER :: i, j ! counters
   REAL(KIND=8),DIMENSION(nxpoints) :: newxgrid
   REAL(KIND=8),DIMENSION(nypoints) :: newygrid
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: v
   TYPE(Csplines),DIMENSION(:),ALLOCATABLE :: newyspline
   CHARACTER(LEN=29),PARAMETER :: routinename="INTERPOL_NEWGRID_BICSPLINES: "
   ! GABBA, GABBA HEY! ---------
   oldnx=size(this%x)
   oldny=size(this%y)
   !
   newxgrid(1) = this%x(1)
   newygrid(1) = this%y(1)
   newxgrid(nxpoints) = this%x(oldnx)
   newygrid(nypoints) = this%y(oldny)
   ! For X, grid parameters
   xinpoints=nxpoints-2
   nxdelta=nxpoints-1
   xdelta=(newxgrid(nxpoints)-newxgrid(1))/DFLOAT(nxdelta)
   ! For Y, grid parameters
   yinpoints=nypoints-2
   nydelta=nypoints-1
   ydelta=(newygrid(nypoints)-newygrid(1))/DFLOAT(nydelta)
   ! generate new x and y grid
   DO i = 1,xinpoints
      newxgrid(i+1)=newxgrid(1)+i*xdelta
   END DO
   DO j = 1,yinpoints
      newygrid(j+1)=newygrid(1)+j*ydelta
   END DO
   ! generate new cubic splines in Y
   ALLOCATE(v(nxpoints,oldny))
   DO j = 1, oldny
      DO i = 1, nxpoints
         v(i,j)=this%xcsplines(j)%getvalue(newxgrid(i))
      END DO
   END DO
   ALLOCATE(newyspline(nxpoints))
   DO i = 1, nxpoints
      CALL newyspline(i)%READ(this%y,v(i,:))
      CALL newyspline(i)%INTERPOL(0.D0,0,0.D0,0)
   END DO
   DEALLOCATE(v)
   ! store all data in v
   ALLOCATE(v(nxpoints,nypoints))
   DO i = 1, nxpoints
      DO j = 1, nypoints
         v(i,j)=newyspline(i)%getvalue(newygrid(j))
      END DO
   END DO
   ! Now generate new interpolation
   CALL this%REBOOT()
   CALL this%READ(newxgrid,newygrid,v)
   CALL this%INTERPOL()
   RETURN
END SUBROUTINE INTERPOL_NEWGRID_BICSPLINES
!##########################################################
! SUBROUTINE: INTERPOL_BICSPLINES
!> @brief
!! Sets coefficients
!
!> @param[in,out] this - bicubic splines object to set coefficients
!----------------------------------------------------------
SUBROUTINE INTERPOL_BICSPLINES(this,filename)
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
   ALLOCATE(this%xcsplines(ny))
   ALLOCATE(this%ycsplines(nx))
   ALLOCATE(this%coeff(nx-1,ny-1,4,4))
   DO j = 1, ny
      CALL this%xcsplines(j)%READ(this%x,this%fgrid(:,j))
      CALL this%xcsplines(j)%INTERPOL(0.D0,0,0.D0,0) ! last and initial 2 splines are equal
   END DO
   DO i = 1, nx
      CALL this%ycsplines(i)%READ(this%y,this%fgrid(i,:))
      CALL this%ycsplines(i)%INTERPOL(0.D0,0,0.D0,0) ! last and initial 2 splines are equal
   END DO
   !
   IF(present(filename)) OPEN (521,FILE=filename,STATUS="replace",ACTION="write")
   DO i = 1, nx-1
      x0=this%x(i)
      x1=this%x(i+1)

      DO j = 1, ny-1

         y0=this%y(j)
         y1=this%y(j+1)

         ! Create matrices for reduced coordinates (0 to x1-x0)
         xmtrx(1,:)=(/1.D0,0.D0,0.D0,0.D0/)
         xmtrx(2,:)=(/1.D0,x1-x0,(x1-x0)**2.D0,(x1-x0)**3.D0/)
         xmtrx(3,:)=(/0.D0,1.D0,0.D0,0.D0/)
         xmtrx(4,:)=(/0.D0,1.D0,2.D0*(x1-x0),3.D0*((x1-x0)**2.D0)/)
         ymtrx(:,1)=(/1.D0,0.D0,0.D0,0.D0/)
         ymtrx(:,2)=(/1.D0,y1-y0,(y1-y0)**2.D0,(y1-y0)**3.D0/)
         ymtrx(:,3)=(/0.D0,1.D0,0.D0,0.D0/)
         ymtrx(:,4)=(/0.D0,1.D0,2.D0*(y1-y0),3.D0*((y1-y0)**2.D0)/)

         smtrx(1,1)=this%fgrid(i,j)
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
         smtrx(4,2)=this%xcsplines(j+1)%getderiv(x1)

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
END SUBROUTINE INTERPOL_BICSPLINES
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
   ! Check sizes of arrays
   SELECT CASE(nx == size(f(:,1)))
      CASE(.FALSE.)
         WRITE(0,*) "d2fdxdy: size mismatch of arrays x and f"
         CALL EXIT(1)
      CASE (.TRUE.)
         ! do nothing
   END SELECT
   SELECT CASE(ny == size(f(1,:)))
      CASE(.FALSE.)
         WRITE(0,*) "d2fdxdy: size mismatch of arrays y and f"
         CALL EXIT(1)
      CASE (.TRUE.)
         ! do nothing
   END SELECT
   ! initialize variables
   hx0=0.D0
   hx1=0.D0
   hy0=0.D0
   hy1=0.D0
   ! Check if we are in a corner, edge or bulk point of the grid
   IF ( i/=1 .AND. i/=nx .AND. j/=1 .AND. j/=ny) THEN ! we are in the 2D bulk
      hx0=x(i+1)-x(i)
      hx1=x(i)-x(i-1)
      hy0=y(j+1)-y(j)
      hy1=y(j)-y(j-1)
      d2fdxdy_finitdiff=((1.D0/hx1)-(1.D0/hx0))*((1.D0/hy1)-(1.D0/hy0))*f(i,j)
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
      hx0=x(i+1)-x(i)
      hx1=x(i)-x(i-1)
      hy1=y(j)-y(j-1)
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
   REAL(KIND=8),POINTER:: x1,x2,y1,y2
   REAL(KIND=8),DIMENSION(:,:),POINTER :: coeff
   REAL(KIND=8),DIMENSION(4) :: vecy,vecx
   INTEGER(KIND=4) :: nx,ny
   INTEGER(KIND=4) :: i,j
   ! Run section
   nx=size(this%x)
   ny=size(this%y)
   SELECT CASE(x(1)<this%x(1) .OR. x(1)>this%x(nx) .OR. x(2)<this%y(1) .OR. x(2)>this%y(ny))
      CASE(.TRUE.)
         WRITE(0,*) "getvalue_bicsplines ERR: requested X,Y outside interpolation limits:"
         WRITE(0,*) "Your request: ",x(1), x(2)
         WRITE(0,*) "X limits: ", this%x(1),this%x(nx)
         WRITE(0,*) "Y limits: ", this%y(1),this%y(ny)
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   DO i = 1, nx-1
      x1 => this%x(i)
      x2 => this%x(i+1)
      IF((x(1).LE.x2).AND.(x(1).GE.x1)) EXIT
   END DO
   DO j = 1, ny-1
      y1 => this%y(j)
      y2 => this%y(j+1)
      IF((x(2).LE.y2).AND.(x(2).GE.y1)) EXIT
   END DO
   ! Now, i and j have the correct value
   coeff => this%coeff(i,j,:,:)
   vecx=(/1.D0,x(1)-x1,(x(1)-x1)**2.D0,(x(1)-x1)**3.D0/)
   vecy=(/1.D0,x(2)-y1,(x(2)-y1)**2.D0,(x(2)-y1)**3.D0/)
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
   REAL(KIND=8),POINTER :: x1,x2,y1,y2
   REAL(KIND=8),DIMENSION(:,:),POINTER :: coeff
   REAL(KIND=8),DIMENSION(4) :: vecy,vecx
   INTEGER(KIND=4) :: nx,ny
   INTEGER(KIND=4) :: i,j
   ! Run section
   nx=size(this%x)
   ny=size(this%y)
   SELECT CASE(x(1)<this%x(1) .OR. x(1)>this%x(nx) .OR. x(2)<this%y(1) .OR. x(2)>this%y(ny))
      CASE(.TRUE.)
         WRITE(0,*) "getderivx_bicsplines ERR: requested X,Y outside interpolation limits:"
         WRITE(0,*) "Your request: ",x(1), x(2)
         WRITE(0,*) "X limits: ", this%x(1),this%x(nx)
         WRITE(0,*) "Y limits: ", this%y(1),this%y(ny)
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   DO i = 1, nx-1
      x1 => this%x(i)
      x2 => this%x(i+1)
      IF((x(1).LE.x2).AND.(x(1).GE.x1)) EXIT
   END DO
   DO j = 1, ny-1
      y1 => this%y(j)
      y2 => this%y(j+1)
      IF((x(2).LE.y2).AND.(x(2).GE.y1)) EXIT
   END DO
   ! Now, i and j have the correct value
   coeff => this%coeff(i,j,:,:)
   vecx=(/0.D0,1.D0,2.D0*(x(1)-x1),3.D0*(x(1)-x1)**2.D0/)
   vecy=(/1.D0,x(2)-y1,(x(2)-y1)**2.D0,(x(2)-y1)**3.D0/)
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
   REAL(KIND=8),POINTER :: x1,x2,y1,y2
   REAL(KIND=8),DIMENSION(:,:),POINTER :: coeff
   REAL(KIND=8),DIMENSION(4) :: vecy,vecx
   INTEGER(KIND=4) :: nx,ny
   INTEGER(KIND=4) :: i,j
   ! Run section
   nx=size(this%x)
   ny=size(this%y)
   SELECT CASE(x(1)<this%x(1) .OR. x(1)>this%x(nx) .OR. x(2)<this%y(1) .OR. x(2)>this%y(ny))
      CASE(.TRUE.)
         WRITE(0,*) "getderivy_bicsplines ERR: requested X,Y outside interpolation limits:"
         WRITE(0,*) "Your request: ",x(1), x(2)
         WRITE(0,*) "X limits: ", this%x(1),this%x(nx)
         WRITE(0,*) "Y limits: ", this%y(1),this%y(ny)
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   DO i = 1, nx-1
      x1 => this%x(i)
      x2 => this%x(i+1)
      IF((x(1).LE.x2).AND.(x(1).GE.x1)) EXIT
   END DO
   DO j = 1, ny-1
      y1 => this%y(j)
      y2 => this%y(j+1)
      IF((x(2).LE.y2).AND.(x(2).GE.y1)) EXIT
   END DO
   ! Now, i and j have the correct value
   coeff => this%coeff(i,j,:,:)
   vecx=(/1.D0,(x(1)-x1),(x(1)-x1)**2.D0,(x(1)-x1)**3.D0/)
   vecy=(/0.D0,1.D0,2.D0*(x(2)-y1),3.D0*(x(2)-y1)**2.D0/)
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
   REAL(KIND=8),POINTER :: x1,x2,y1,y2
   REAL(KIND=8),DIMENSION(:,:),POINTER :: coeff
   REAL(KIND=8),DIMENSION(4) :: vecy,vecx
   INTEGER(KIND=4) :: nx,ny
   INTEGER(KIND=4) :: i,j
   ! Run section
   nx=size(this%x)
   ny=size(this%y)
   SELECT CASE(x(1)<this%x(1) .OR. x(1)>this%x(nx) .OR. x(2)<this%y(1) .OR. x(2)>this%y(ny))
      CASE(.TRUE.)
         WRITE(0,*) "getderivxy_bicsplines ERR: requested X,Y outside interpolation limits:"
         WRITE(0,*) "Your request: ",x(1), x(2)
         WRITE(0,*) "X limits: ", this%x(1),this%x(nx)
         WRITE(0,*) "Y limits: ", this%y(1),this%y(ny)
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   DO i = 1, nx-1
      x1 => this%x(i)
      x2 => this%x(i+1)
      IF((x(1).LE.x2).AND.(x(1).GE.x1)) EXIT
   END DO
   DO j = 1, ny-1
      y1 => this%y(j)
      y2 => this%y(j+1)
      IF((x(2).LE.y2).AND.(x(2).GE.y1)) EXIT
   END DO
   ! Now, i and j have the correct value
   coeff => this%coeff(i,j,:,:)
   vecx=(/0.D0,1.D0,2.D0*(x(1)-x1),3.D0*(x(1)-x1)**2.D0/)
   vecy=(/0.D0,1.D0,2.D0*(x(2)-y1),3.D0*(x(2)-y1)**2.D0/)
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
   ydelta=Ly/DFLOAT(nydelta)
   ! Let's go!
   ! 1st XY point
   OPEN(11,file=filename,status="replace")
   r(1) = xmin
   r(2) = ymin
   WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r),this%getderivxy(r)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r),this%getderivxy(r)
   END DO
   r(2) = ymax
   WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r),this%getderivxy(r)
   ! inpoints in XY
   DO i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r),this%getderivxy(r)
      DO j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r),this%getderivxy(r)
      END DO
      r(2) = ymax
      WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r),this%getderivxy(r)
   END DO
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r),this%getderivxy(r)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r),this%getderivxy(r)
   END DO
   r(2) = ymax
   WRITE(11,*) r(1),r(2),this%getvalue(r),this%getderivx(r),this%getderivy(r),this%getderivxy(r)
   CLOSE(11)
   WRITE(*,*) "PLOT_XYMAP_BICSPLINES: Graph created: ",filename
   RETURN
END SUBROUTINE PLOT_XYMAP_BICSPLINES
!###############################################################
! SUBROUTINE: PLOT_1D_BICSPLINES
!##############################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut (X,Y)
!
!> @param[in] this - Interpolation 2D object
!> @param[in] filename - Name of the file to print the output
!> @param[in] init_xy - Initial position to start the scan
!> @param[in] npoints - Number of points
!> @param[in] angle - Defines the direction of the 1D cut respect to X axis (rad)
!> @param[in] L - Length of the curve in a.u.
!
!> @author A.S. Muzas
!> @date 20/Feb/2014
!> @version 1.0
!---------------------------------------------------------------
SUBROUTINE PLOT_1D_BICSPLINES(this,filename,init_xy,npoints,angle,L)
   IMPLICIT NONE
   CLASS(Bicsplines),INTENT(IN) :: this
   REAL*8,DIMENSION(2),INTENT(IN) :: init_xy
   INTEGER,INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),INTENT(IN) :: angle
   REAL*8,INTENT(IN) :: L
   ! Local variables
   REAL*8 :: xmin, ymin, xmax, ymax
   REAL*8, DIMENSION(2) :: r
   REAL(KIND=8) :: alpha
   REAL*8 :: delta,s
   INTEGER :: inpoints, ndelta
   INTEGER :: i, j ! counters
   ! GABBA, GABBA HEY! ---------
   alpha = angle*PI/180.D0
   xmin = init_xy(1)
   ymin = init_xy(2)
   xmax = init_xy(1)+L*dcos(alpha)
   ymax = init_xy(2)+L*dsin(alpha)
   ! For X, grid parameters
   inpoints=npoints-2
   ndelta=npoints-1
   delta=L/DFLOAT(ndelta)
   ! Let's go!
   ! 1st XY point
   OPEN(11,file=filename,status="replace")
   r(1) = xmin
   r(2) = ymin
   s = dsqrt(r(1)**2.D0+r(2)**2.D0)
   WRITE(11,*) s,this%getvalue(r),dcos(alpha)*this%getderivx(r)+dsin(alpha)*this%getderivy(r),&
      this%getderivx(r),this%getderivy(r)
   DO i =1, inpoints
      r(1)=xmin+(DFLOAT(i)*delta)*DCOS(alpha)
      r(2)=ymin+(DFLOAT(i)*delta)*DSIN(alpha)
      s = DSQRT(r(1)**2.D0+r(2)**2.D0)
      WRITE(11,*) s,this%getvalue(r),dcos(alpha)*this%getderivx(r)+dsin(alpha)*this%getderivy(r),&
         this%getderivx(r),this%getderivy(r)
   END DO
   r(1) = xmax
   r(2) = ymax
   s = DSQRT(r(1)**2.D0+r(2)**2.D0)
   WRITE(11,*) s,this%getvalue(r),dcos(alpha)*this%getderivx(r)+dsin(alpha)*this%getderivy(r),&
      this%getderivx(r),this%getderivy(r)
   CLOSE(11)
   WRITE(*,*) "PLOT_1D_BICSPLINES: Graph created: ",filename
   RETURN
END SUBROUTINE PLOT_1D_BICSPLINES
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
      CALL this%xcsplines(i)%PLOT(npoints,filename)
   END DO
   DO i = 1, nx
      WRITE(filename,'(I4,A12)') i,"-yspline.dat"
      filename=adjustl(filename)
      CALL this%ycsplines(i)%PLOT(npoints,filename)
   END DO
   RETURN
END SUBROUTINE PLOT_SPLINES_BICSPLINES
!###########################################################
!# SUBROUTINE: PLOT_FINITEDIFF_AT_GRID
!###########################################################
!> @brief
!! Creates a file with finite differences results at gridpoints
!
!> @details
!! - Useful to check if finite differences give good results in the
!!   proposed grid
!! - FORMAT:
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 21/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE PLOT_DUALDERIVS_AT_GRID(this,filename)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Bicsplines),TARGET,INTENT(IN):: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   INTEGER(KIND=4) :: nx,ny
   REAL(KIND=8),DIMENSION(:),POINTER :: x,y
   ! Run section
   nx=size(this%x)
   ny=size(this%y)
   OPEN (10,FILE=filename,STATUS="replace",ACTION="write")
   DO i = 1, nx
      DO j = 1, ny
         WRITE(10,*) this%x(i),this%y(j),d2fdxdy_finitdiff(i,j,this%x,this%y,this%fgrid)
      END DO
   END DO
   CLOSE(10)
   RETURN
END SUBROUTINE PLOT_DUALDERIVS_AT_GRID
END MODULE BICSPLINES_MOD

!########################################################
! MODULE : FOURIER1D_2_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes FOURIER1D_mod in its scope
!########################################################
MODULE FOURIER1D_2_MOD
USE FOURIER1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: term_A1
!> @brief
!! Child class of abstract termcalculator. Strcture to avoid unnecessary switches
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Termcalculator) :: term_2
PRIVATE
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_2
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_2
END TYPE term_2
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_2
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
TYPE,EXTENDS(Fourier1d):: Fourier1d_2
   CONTAINS
      ! Set block
      PROCEDURE,PUBLIC :: initializeTerms => initializeTerms_FOURIER1D_2
END TYPE FOURIER1D_2
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: initializeTerm_FOURIER1D_2
!###########################################################
!> @brief
!! Sets Term for this fourier series
!-----------------------------------------------------------
subroutine initializeTerms_FOURIER1D_2(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier1d_2),intent(inout):: this
   ! Run section
   allocate( Term_2::this%term )
   return
end subroutine InitializeTerms_FOURIER1D_2
!###########################################################
!# FUNCTION: termfou1d_2
!###########################################################
!-----------------------------------------------------------
function termfou1d_2(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term_2),intent(in):: this
   integer(kind=4),intent(in):: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Local
   character(len=1):: auxChar
   ! Parameters
   character(len=*),parameter:: routinename='termfou_2: '
   ! Run section
   auxChar=trim(irrep)
   select case( auxChar )
   case('A')
      ! check parity
      select case( parity )
      case('+')
         select case( mod(kpoint,2)==0 )
         case(.true.) ! is par
            answer=dcos(dfloat(kpoint)*x)
         case(.false.) ! is odd
            write(0,*) 'ERR '//routinename//'bad combination of parity and Kpoint'
            call exit(1)
         end select
      case('-')
         select case( mod(kpoint,2)==0 )
         case(.true.) ! is par
            answer=dsin(dfloat(kpoint)*x)
         case(.false.) ! is odd
            write(0,*) 'ERR '//routinename//'bad combination of parity and Kpoint'
            call exit(1)
         end select
      case('o')
         select case( mod(kpoint,2)==0 )
         case(.true.) ! is par
            answer=dsin(dfloat(kpoint)*x)+dcos(dfloat(kpoint)*x)
         case(.false.) ! is odd
            write(0,*) 'ERR '//routinename//'bad combination of parity and Kpoint'
            call exit(1)
         end select
      case default
         write(0,*) 'ERR '//routinename//'bad parity symbol: "'//parity//'"'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep not implemented: "'//auxChar//'"'
      write(0,*) 'Implemented ones: A'
      call exit(1)
   end select

   return
end function termfou1d_2
!###########################################################
!# FUNCTION: termfou1d_dx_2
!###########################################################
!-----------------------------------------------------------
function termfou1d_dx_2(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term_2),intent(in):: this
   integer(kind=4),intent(in):: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Local
   character(len=1):: auxChar
   ! Parameters
   character(len=*),parameter:: routinename='termfou_dx_2: '
   ! Run section
   auxChar=trim(irrep)
   select case( auxChar )
   case('A')
      ! check parity
      select case( parity )
      case('+')
         select case( mod(kpoint,2)==0 )
         case(.true.) ! kpoint is par
            answer=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
         case(.false.) ! kpoint is odd
            write(0,*) 'ERR '//routinename//'bad combination of parity and Kpoint'
            call exit(1)
         end select
      case('-')
         select case( mod(kpoint,2)==0 )
         case(.true.) ! is par
            answer=dfloat(kpoint)*dcos(dfloat(kpoint)*x)
         case(.false.) ! is odd
            write(0,*) 'ERR '//routinename//'bad combination of parity and Kpoint'
            call exit(1)
         end select
      case('o')
         select case( mod(kpoint,2)==0 )
         case(.true.) ! is par
            answer=dfloat(kpoint)*( dcos(dfloat(kpoint)*x)-dsin(dfloat(kpoint)*x) )
         case(.false.) ! is odd
            write(0,*) 'ERR '//routinename//'bad combination of parity and Kpoint'
            call exit(1)
         end select
      case default
         write(0,*) 'ERR '//routinename//'bad parity symbol: "'//parity//'"'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep not implemented: "'//auxChar//'"'
      write(0,*) 'Implemented ones: A'
      call exit(1)
   end select

   return
end function termfou1d_dx_2
END MODULE FOURIER1D_2_MOD

!########################################################
! MODULE : FOURIER1D_4MM_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes FOURIER1D_mod in its scope
!########################################################
MODULE FOURIER1D_4MM_MOD
USE FOURIER1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: term_A1
!> @brief
!! Child class of abstract termcalculator. Strcture to avoid unnecessary switches
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Termcalculator) :: term_4mm
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_4mm
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_4mm
END TYPE term_4mm
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_4MM
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
type,extends(Fourier1d):: Fourier1d_4mm
   contains
      ! Set block
      procedure,public:: initializeTerms => initializeTerms_FOURIER1D_4MM
end type Fourier1d_4mm
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_IRREP_FOURIER1D_4MM
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!-----------------------------------------------------------
subroutine initializeTerms_FOURIER1D_4MM(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier1d_4mm),intent(inout):: this
   ! Run section
   allocate(Term_4mm::this%term)
   return
end subroutine initializeTerms_FOURIER1D_4MM
!###########################################################
!# FUNCTION: termfou1d_4mm
!###########################################################
!-----------------------------------------------------------
function termfou1d_4mm(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term_4mm),intent(in):: this
   integer(kind=4),intent(in):: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou1d_4mm: '
   ! Run section
   select case( irrep )
   case('A1')
      ! kpoint check
      select case( mod(kpoint,4)/=0 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('+')
         answer=dcos(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case('A2')
      ! kpoint check
      select case( mod(kpoint,4)/=0 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('-')
         answer=dsin(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case('B1')
      ! kpoint check
      select case( mod(kpoint,4)/=2 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('+')
         answer=dcos(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case('B2')
      ! kpoint check
      select case( mod(kpoint,4)/=2 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('-')
         answer=dsin(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case('E','EE')
      ! kpoint check
      select case( mod(kpoint,2)==0 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('+')
         answer=dcos(dfloat(kpoint)*x)
      case('-')
         answer=dsin(dfloat(kpoint)*x)
      case('o')
         answer=dsin(dfloat(kpoint)*x)+dcos(dfloat(kpoint)*x)
      case('n')
         answer=dsin(dfloat(kpoint)*x)-dcos(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep "'//irrep//'" not implemented yet'
      write(0,*) 'Implemented ones: A1'
      call exit(1)
   end select

   return
end function termfou1d_4mm
!###########################################################
!# FUNCTION: termfou1d_dx_4mm_A1
!###########################################################
!-----------------------------------------------------------
function termfou1d_dx_4mm(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term_4mm),intent(in) :: this
   integer(kind=4),intent(in) :: kpoint
   real(kind=8),intent(in) :: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routineName='termfou1d_dx_4mm: '
    select case( irrep )
   case('A1')
      select case( parity )
      case('+')
         answer=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case('A2')
      ! kpoint check
      select case( mod(kpoint,4)/=0 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('-')
         answer=dfloat(kpoint)*dcos(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case('B1')
      ! kpoint check
      select case( mod(kpoint,4)/=2 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('+')
         answer=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case('B2')
      ! kpoint check
      select case( mod(kpoint,4)/=2 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('-')
         answer=dfloat(kpoint)*dcos(dfloat(kpoint)*x)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case('E','EE')
      ! kpoint check
      select case( mod(kpoint,2)==0 )
      case(.true.)
         write(0,*) 'ERR '//routinename//' bad kpoint number'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! parity check
      select case( parity )
      case('+')
         answer=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
      case('-')
         answer=dfloat(kpoint)*dcos(dfloat(kpoint)*x)
      case('o')
         answer=dfloat(kpoint)*( dcos(dfloat(kpoint)*x)-dsin(dfloat(kpoint)*x) )
      case('n')
         answer=dfloat(kpoint)*( dcos(dfloat(kpoint)*x)+dsin(dfloat(kpoint)*x) )
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep "'//irrep//'" not implemented yet'
      write(0,*) 'Implemented ones: A1'
      call exit(1)
   end select

   return
end function termfou1d_dx_4mm
END MODULE FOURIER1D_4MM_MOD

!########################################################
! MODULE : FOURIER1D_M_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes FOURIER1D_mod in its scope
!########################################################
MODULE FOURIER1D_M_MOD
USE FOURIER1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: term_Ap
!> @brief
!! Child class of abstract termcalculator. Strcture to avoid unnecessary switches
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Termcalculator) :: term_M
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_M
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_M
END TYPE term_M
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_M
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
TYPE,EXTENDS(FOURIER1D):: Fourier1d_M
   CONTAINS
      ! Set block
      PROCEDURE,PUBLIC :: initializeTerms => initializeTerms_FOURIER1D_M
END TYPE FOURIER1D_M
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_IRREP_FOURIER1D_M
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine initializeTerms_FOURIER1D_M(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier1d_M),intent(inout):: this
   ! Run section
   allocate(term_M::this%term)
   return
end subroutine initializeTerms_FOURIER1D_M
!###########################################################
!# FUNCTION: termfou1d_M
!###########################################################
!-----------------------------------------------------------
function termfou1d_M(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(term_M),intent(in):: this
   integer(kind=4),intent(in):: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou_M: '
   ! Run section
   select case( irrep )
   case('Ap')
      select case( parity )
      case('+')
         answer=dcos(dfloat(kpoint)*(x+this%getShift()))
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep "'//irrep//'" not implemented yet'
      write(0,*) 'Implemented ones: A1'
      call exit(1)
   end select
   return
end function termfou1d_M
!###########################################################
!# FUNCTION: termfou1d_dx_M
!###########################################################
!-----------------------------------------------------------
function termfou1d_dx_M(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(term_M),intent(in):: this
   integer(kind=4),intent(in):: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou_M: '
   ! Run section
   select case( irrep )
   case('Ap')
      select case( parity )
      case('+')
         answer=-dsin(dfloat(kpoint)*(x+this%getshift()))*dfloat(kpoint)
      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep "'//irrep//'" not implemented yet'
      write(0,*) 'Implemented ones: A1'
      call exit(1)
   end select
   return
end function termfou1d_dx_M

END MODULE FOURIER1D_M_MOD

!########################################################
! MODULE : FOURIER1D_M45_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes FOURIER1D_mod in its scope
!########################################################
MODULE FOURIER1D_M45_MOD
USE FOURIER1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: term_Ap
!> @brief
!! Child class of abstract termcalculator. Strcture to avoid unnecessary switches
!----------------------------------------------------------------
TYPE,EXTENDS(Termcalculator):: term_M45
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_M45
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_M45
END TYPE term_M45
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_M45
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
TYPE,EXTENDS(FOURIER1D):: Fourier1d_M45
   CONTAINS
      ! Set block
      PROCEDURE,PUBLIC :: initializeTerms => initializeTerms_FOURIER1D_M45
END TYPE FOURIER1D_M45
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_IRREP_FOURIER1D_M45
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine initializeTerms_FOURIER1D_M45(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier1d_M45),intent(inout):: this
   ! Run section
   allocate(term_M45::this%term)
   return
end subroutine initializeTerms_FOURIER1D_M45
!###########################################################
!# FUNCTION: termfou1d_M45_Ap
!###########################################################
!-----------------------------------------------------------
function termfou1d_M45(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(term_M45),intent(in) :: this
   integer(kind=4),intent(in) :: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou_M45: '
   ! Run section
   select case( irrep )
   case('Ap')
      select case( parity )
      case('+')
         select case( mod(kpoint,4) )
         case(0)
            answer=dcos(dfloat(kpoint)*x)
         case default
            write(0,*) "termfou1d_m45 ERR: incorrect Kpoint-parity combination"
            call exit(1)
         end select

      case('-')
         select case( mod(kpoint,4) )
         case(2)
            answer=dsin(dfloat(kpoint)*x)
         case default
            write(0,*) "termfou1d_m45 ERR: incorrect point-parity combination"
            call exit(1)
         end select

      case('o')
         select case( mod(kpoint,4) )
         case(1)
            answer=dcos(dfloat(kpoint)*x)+dsin(dfloat(kpoint)*x)
         case(3)
            answer=dcos(dfloat(kpoint)*x)-dsin(dfloat(kpoint)*x)
         case default
            write(0,*) "termfou1d_m45 ERR: incorrect Kpoint-parity combination"
            call exit(1)
         end select

      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep "'//irrep//'" not implemented yet'
      write(0,*) 'Implemented ones: Ap'
      call exit(1)
   end select
   return
end function termfou1d_M45
!###########################################################
!# FUNCTION: termfou1d_dx_M45_Ap
!###########################################################
!-----------------------------------------------------------
function termfou1d_dx_M45(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(term_M45),intent(in) :: this
   integer(kind=4),intent(in) :: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou_dx_M45: '
   ! Run section
   select case( irrep )
   case('Ap')
      select case( parity )
      case('+')
         select case( mod(kpoint,4) )
         case(0)
            answer=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
         case default
            write(0,*) "termfou1d_m45 ERR: incorrect Kpoint-parity combination"
            call exit(1)
         end select

      case('-')
         select case( mod(kpoint,4) )
         case(2)
            answer=dfloat(kpoint)*dcos(dfloat(kpoint)*x)
         case default
            write(0,*) "termfou1d_m45 ERR: incorrect Kpoint-parity combination"
            call exit(1)
         end select

      case('o')
         select case( mod(kpoint,4) )
         case(1)
            answer=dfloat(kpoint)*( dcos(dfloat(kpoint)*x)-dsin(dfloat(kpoint)*x) )
         case(3)
            answer=-dfloat(kpoint)*( dcos(dfloat(kpoint)*x)+dsin(dfloat(kpoint)*x) )
         case default
            write(0,*) "termfou1d_m45 ERR: incorrect Kpoint-parity combination"
            call exit(1)
         end select

      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep "'//irrep//'" not implemented yet'
      write(0,*) 'Implemented ones: Ap'
      call exit(1)
   end select
   return
end function termfou1d_dx_M45
END MODULE FOURIER1D_M45_MOD

!########################################################
! MODULE : FOURIER1D_MM2_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes FOURIER1D_mod in its scope
!########################################################
MODULE FOURIER1D_MM2_MOD
USE FOURIER1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: term_A1
!> @brief
!! Child class of abstract termcalculator. Strcture to avoid unnecessary switches
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Termcalculator) :: term_mm2
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_MM2
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_MM2
END TYPE term_mm2
!///////////////////////////////////////////////////////////////////////////////
! TYPE: FOURIER1D_MM2
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!-------------------------------------------------------------------------------
TYPE,EXTENDS(FOURIER1D):: Fourier1d_MM2
   CONTAINS
      PROCEDURE,PUBLIC :: initializeTerms => initializeTerms_FOURIER1D_MM2
END TYPE FOURIER1D_MM2
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: initializeTerms_FOURIER1D_MM2
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE initializeTerms_FOURIER1D_MM2(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d_MM2),INTENT(INOUT):: this
   ! Run section
   ALLOCATE(Term_mm2::this%term)
   RETURN
END SUBROUTINE initializeTerms_FOURIER1D_MM2
!###########################################################
!# FUNCTION: termfou1d_MM2_A1
!###########################################################
!-----------------------------------------------------------
function termfou1d_MM2(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term_mm2),intent(in):: this
   integer(kind=4),intent(in):: kpoint
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   real(kind=8),intent(in):: x
   ! Dummy out variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou1d_MM2: '
   ! Run section
   select case( irrep )
   case('A1')
      ! check KPOINT
      select case( mod(kpoint,2)==1 )
      case(.true.)
         write(0,*) 'ERR '//routinename//'bad irrep. Given: "'//irrep//'"'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! check parity
      select case( parity )
      case('+')
         answer=dcos(dfloat(kpoint)*x)
      case default
         write(0,*) 'ERR '//routinename//'bad parity. Given: "'//parity//'"'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep not implemented. Given: "'//irrep//'"'
      call exit(1)
   end select
   return
end function termfou1d_MM2
!###########################################################
!# FUNCTION: termfou1d_dx_MM2_A1
!###########################################################
!-----------------------------------------------------------
function termfou1d_dx_MM2(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term_mm2),intent(in):: this
   integer(kind=4),intent(in):: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy out variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou1d_MM2: '
   ! Run section
   select case( irrep )
   case('A1')
      ! check KPOINT
      select case( mod(kpoint,2)==1 )
      case(.true.)
         write(0,*) 'ERR '//routinename//'bad irrep'
         call exit(1)
      case(.false.)
         ! do nothing
      end select
      ! check parity
      select case( parity )
      case('+')
         answer=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
      case default
         write(0,*) 'ERR '//routinename//'bad parity'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep not implemented'
      call exit(1)
   end select
   return
end function termfou1d_dx_MM2

END MODULE FOURIER1D_MM2_MOD

!########################################################
! MODULE : FOURIER1D_M45M1352_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes FOURIER1D_mod in its scope
!########################################################
MODULE FOURIER1D_m45m135_2_MOD
USE FOURIER1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: term_A1
!> @brief
!! Child class of abstract termcalculator. Strcture to avoid unnecessary switches
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Termcalculator) :: term_m45m135_2
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_m45m135_2
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_m45m135_2
END TYPE term_m45m135_2
!///////////////////////////////////////////////////////////////////////////////
! TYPE: FOURIER1D_m45m135_2
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!-------------------------------------------------------------------------------
TYPE,EXTENDS(FOURIER1D):: Fourier1d_m45m135_2
   CONTAINS
      PROCEDURE,PUBLIC :: initializeTerms => initializeTerms_FOURIER1D_m45m135_2
END TYPE FOURIER1D_m45m135_2
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_IRREP_FOURIER1D_m45m135_2
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine initializeTerms_FOURIER1D_m45m135_2(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier1d_m45m135_2),intent(inout):: this
   ! Run section
   allocate(Term_m45m135_2::this%term)
   return
end subroutine initializeTerms_FOURIER1D_m45m135_2
!###########################################################
!# FUNCTION: termfou1d_m45m135_2_A1
!###########################################################
!-----------------------------------------------------------
function termfou1d_m45m135_2(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term_m45m135_2),intent(in):: this
   integer(kind=4),intent(in):: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou_M45M135_2: '
   ! Run section
   select case( irrep )
   case('A1')
      select case( parity )
      case('+')
         select case( mod(kpoint,4) )
         case(0)
            answer=dcos(dfloat(kpoint)*x)
         case default
            write(0,*) "termfou1d_m45 ERR: incorrect Kpoint-parity combination"
            call exit(1)
         end select

      case('-')
         select case( mod(kpoint,4) )
         case(2)
            answer=dsin(dfloat(kpoint)*x)
         case default
            write(0,*) "termfou1d_m45 ERR: incorrect Kpoint-parity combination"
            call exit(1)
         end select

      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep "'//irrep//'" not implemented yet'
      write(0,*) 'Implemented ones: Ap'
      call exit(1)
   end select
   return
end function termfou1d_m45m135_2
!###########################################################
!# FUNCTION: termfou1d_dx_m45m135_2
!###########################################################
!-----------------------------------------------------------
function termfou1d_dx_m45m135_2(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term_m45m135_2),intent(in):: this
   integer(kind=4),intent(in):: kpoint
   real(kind=8),intent(in):: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou_dx_M45M135_2: '
   ! Run section
   select case( irrep )
   case('A1')
      select case( parity )
      case('+')
         select case( mod(kpoint,4) )
         case(0)
            answer=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
         case default
            write(0,*) "termfou1d_m45 ERR: incorrect Kpoint-parity combination"
            call exit(1)
         end select

      case('-')
         select case( mod(kpoint,4) )
         case(2)
            answer=dfloat(kpoint)*dcos(dfloat(kpoint)*x)
         case default
            write(0,*) "termfou1d_m45 ERR: incorrect Kpoint-parity combination"
            call exit(1)
         end select

      case default
         write(0,*) routinename//'ERR: parity "'//parity//'" not implemented and/or does not exist'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep "'//irrep//'" not implemented yet'
      write(0,*) 'Implemented ones: Ap'
      call exit(1)
   end select
   return
end function termfou1d_dx_m45m135_2
END MODULE FOURIER1D_m45m135_2_MOD

!########################################################
! MODULE : FOURIER1D_E_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes FOURIER1D_mod in its scope
!########################################################
MODULE FOURIER1D_E_MOD
USE FOURIER1D_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: term_A
!> @brief
!! Child class of abstract termcalculator. Strcture to avoid unnecessary switches
!----------------------------------------------------------------
TYPE,EXTENDS(Termcalculator) :: term_E
PRIVATE
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_E
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_E
END TYPE term_E
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_E
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
TYPE,EXTENDS(FOURIER1D):: Fourier1d_E
   CONTAINS
      ! Set block
      PROCEDURE,PUBLIC :: initializeTerms => initializeTerms_FOURIER1D_E
END TYPE FOURIER1D_E
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: initializeTerms_FOURIER1D_E
!###########################################################
!> @brief
!! Sets Terms for this fourier series
!-----------------------------------------------------------
subroutine initializeTerms_FOURIER1D_E(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier1d_E),intent(inout):: this
   ! Run section
   allocate(Term_E::this%term)
   return
end subroutine initializeTerms_FOURIER1D_E
!###########################################################
!# FUNCTION: termfou1d_E
!###########################################################
!-----------------------------------------------------------
function termfou1d_E(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term_E),intent(in) :: this
   integer(kind=4),intent(in) :: kpoint
   real(kind=8),intent(in) :: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou_E: '
   ! Run section
   select case( irrep )
   case('A')
      ! check parity
      select case( parity )
      case('+')
            answer=dcos(dfloat(kpoint)*x)
      case('-')
            answer=dsin(dfloat(kpoint)*x)
      case('o')
            answer=dsin(dfloat(kpoint)*x)+dcos(dfloat(kpoint)*x)
      case default
         write(0,*) 'ERR '//routinename//'bad parity symbol'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep not implemented'
      write(0,*) 'Implemented ones: A'
      call exit(1)
   end select

   return
end function termfou1d_E
!###########################################################
!# FUNCTION: termfou1d_dx_E
!###########################################################
!-----------------------------------------------------------
function termfou1d_dx_E(this,kpoint,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term_E),intent(in) :: this
   integer(kind=4),intent(in) :: kpoint
   real(kind=8),intent(in) :: x
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou_dx_E: '
   ! Run section
   select case( irrep )
   case('A')
      ! check parity
      select case( parity )
      case('+')
            answer=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
      case('-')
            answer=dfloat(kpoint)*dcos(dfloat(kpoint)*x)
      case('o')
            answer=dfloat(kpoint)*( dcos(dfloat(kpoint)*x)-dsin(dfloat(kpoint)*x) )
      case default
         write(0,*) 'ERR '//routinename//'bad parity symbol'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//'irrep not implemented'
      write(0,*) 'Implemented ones: A'
      call exit(1)
   end select

   return
end function termfou1d_dx_E
END MODULE FOURIER1D_E_MOD

!###################################################################################################
! MODULE: INTERPOL3D
!
!> @brief
!! Module that manages different interpolation schemes for one variable
!
!> @details
!! All types and subprograms intended to create interpolations in 1D should be placed inside this module
!##################################################################################################
module INTERPOL3D_MOD
! Initial declarations
#ifdef DEBUG
   use DEBUG_MOD
#endif
implicit none
!//////////////////////////////////////////////////////////////////////
! TYPE: interpol3d
!> @brief
!! Generic type of three dimensional interpolations
!-----------------------------------------------------------------------
type,abstract:: Interpol3d
   ! public atributes
   integer(kind=4),public:: n
   real(kind=8),dimension(:,:),allocatable,public:: x
   real(kind=8),dimension(:),allocatable,public:: f
   contains
      ! Initialization block
      procedure,non_overridable,public:: read => read_INTERPOL3D
      ! Get functions block
      procedure(getValue_INTERPOL3D),public,deferred:: getValue  ! child types, override this
      procedure(getValue_INTERPOL3D),public,deferred:: getDeriv1 ! child types, override this
      procedure(getValue_INTERPOL3D),public,deferred:: getDeriv2 ! child types, override this
      procedure(getValue_INTERPOL3D),public,deferred:: getDeriv3 ! child types, override this
end type Interpol3d
!//////////////////////////////////////////////////////////////////////
abstract interface
   !###########################################################
   !# FUNCTION: getvalue_interpol1d
   !###########################################################
   !> @brief
   !! Dummy function. Override it!!
   !-----------------------------------------------------------
   function getValue_INTERPOL3D(this,x) result(answer)
      import Interpol3d
      class(Interpol3d),target,intent(in):: this
      real(kind=8),dimension(3),intent(in):: x
      real(kind=8):: answer
   end function getValue_INTERPOL3D
end interface
contains
!######################################################################
! SUBROUTINE: READ_INTERPOL3D #########################################
!######################################################################
!> @brief
!! Read main parameters for a 3D interpolation from arguments
!----------------------------------------------------------------------
subroutine read_INTERPOL3D(this,x,f)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Interpol3d),intent(inout):: this
   real(kind=8),dimension(:,:),intent(in)::x
   real(kind=8), dimension(:),intent(in):: f
   ! Run section -------
   if ( size(x(:,1))/=size(f) ) then
      write(0,*) "READ_INTERPOL3D: dimensions mismatch between x and f"
      call exit(1)
   end if
   this%n=size(x(:,1))
   allocate( this%x(this%n,3), source=x(:,:) )
   allocate( this%f(this%n),   source=f(:)   )
   return
end subroutine READ_INTERPOL3D
end module INTERPOL3D_MOD

!########################################################
! MODULE : FOURIER3D_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes Interpol1d_mod in its scope
!########################################################
module FOURIER3D_MOD
use INTERPOL3D_MOD, only: Interpol3d
use MATHS_MOD, only: inv_mtrx
use FOURIER1D_MOD, only: TermCalculator
use FOURIER2D_MOD, only: TermCalculator2d
#ifdef DEBUG
use DEBUG_MOD, only: verbose_write, debug_write
#endif
implicit none
!/////////////////////////////////////////////////////////////////
! TYPE: Termcalculator
!> @brief
!! Abstract class to calculate terms of the series avoiding the use of
!! unnecessary switches
!----------------------------------------------------------------
type,abstract:: TermCalculator3d
   ! public atributes
   class(TermCalculator),allocatable,public:: angleFourier
   class(TermCalculator2d),allocatable,public:: xyFourier
   contains
      procedure(getvalue_termcalculator_example),public,deferred:: getValue  ! deferred
      procedure(getvalue_termcalculator_example),public,deferred:: getDeriv1 ! deferred
      procedure(getvalue_termcalculator_example),public,deferred:: getDeriv2 ! deferred
      procedure(getvalue_termcalculator_example),public,deferred:: getDeriv3 ! deferred
end type TermCalculator3d
!
abstract interface
   !###########################################################
   !# FUNCTION: getvalue_termcalculator_example
   !###########################################################
   !> @brief
   !! Just an example that child objects should override
   !-----------------------------------------------------------
   function getvalue_termcalculator_example(this,k,parityXY,irrepXY,l,parityAngle,irrepAngle,x) result(answer)
      import TermCalculator3d
      class(TermCalculator3d),intent(in):: this
      integer(kind=4),dimension(2),intent(in):: k
      integer(kind=4),intent(in):: l
      character(len=1),intent(in):: parityXY
      character(len=1),intent(in):: parityAngle
      character(len=2),intent(in):: irrepXY
      character(len=2),intent(in):: irrepAngle
      real(kind=8),dimension(3),intent(in):: x
      real(kind=8):: answer
   end function getvalue_termcalculator_example
   !-------------------------------------------------------------
end interface
!
!/////////////////////////////////////////////////////////////////////////////
! TYPE: FOURIER3D
!> @brief
!! Class to store all information needed for a 3D REAL combined fourier interpolation
!----------------------------------------------------------------------------
type,abstract,extends(Interpol3d):: Fourier3d
   ! public atributes
   class(TermCalculator3d),allocatable,public:: term
   integer(kind=4),dimension(:,:),allocatable,public:: kListXY
   character(len=1),dimension(:),allocatable,public:: parityListXY
   character(len=2),dimension(:),allocatable,public:: irrepListXY
   integer(kind=4),dimension(:),allocatable,public:: kListAngle
   character(len=1),dimension(:),allocatable,public:: parityListAngle
   character(len=2),dimension(:),allocatable,public:: irrepListAngle
   ! private atributes
   real(kind=8),dimension(:),allocatable,private:: coeff
   real(kind=8),dimension(:,:),allocatable,private:: extracoeff
   real(kind=8),dimension(:,:),allocatable,private:: extrafuncs
   contains
      ! initialize block
      procedure(initializeTerms_FOURIER3D),public,deferred:: initializeTerms
      ! get block
      procedure,public,non_overridable:: getValue  => getvalue_FOURIER3D
      procedure,public,non_overridable:: getDeriv1 => getDeriv1_FOURIER3D
      procedure,public,non_overridable:: getDeriv2 => getDeriv2_FOURIER3D
      procedure,public,non_overridable:: getDeriv3 => getDeriv3_FOURIER3D
      ! set block
      procedure,public,non_overridable:: setKlist => setKlist_FOURIER3D
      procedure,public,non_overridable:: setParityList => setParityList_FOURIER3D
      procedure,public,non_overridable:: setIrrepList => setIrrepList_FOURIER3D
      ! tools
      procedure,public,non_overridable:: interpol => interpol_FOURIER3D
      procedure,public,non_overridable:: add_morefuncs => add_more_funcs_FOURIER3D
      procedure,public,non_overridable:: get_allfuncs_and_derivs => get_allfunc_and_derivs_FOURIER3D
end type Fourier3d
abstract interface
   !###########################################################
   !# SUBROUTINE: SET_IRREP_FOURIER3D
   !###########################################################
   !> @brief
   !! Sets irrep for this fourier series. Should be overriden by
   !! child non-abstract classes
   !-----------------------------------------------------------
   subroutine initializeTerms_FOURIER3D(this)
      import Fourier3d
      class(Fourier3d),intent(inout):: this
   end subroutine initializeTerms_FOURIER3D
end interface
!/////////////////////////////////////////////////////////////////////////////
contains
!###################################################################
!# SUBROUTINE: setKlist_FOURIER3D
!###################################################################
!> @brief
!! Common set subroutine. Sets Klist atribute of a FOURIER3D object
!-------------------------------------------------------------------
subroutine setKlist_FOURIER3D(this,kListXY,kListAngle)
   implicit none
   ! I/O variables
   class(Fourier3d),intent(inout):: this
   integer(kind=4),dimension(:,:),intent(in):: kListXY
   integer(kind=4),dimension(:),intent(in):: kListAngle
   ! Local variables
   integer(kind=4):: n
   ! Run section
   n=size(kListXY(:,1))
   allocate( this%kListXY(n,2),  source=kListXY(:,1:2)  )
   allocate( this%kListAngle(n), source=kListAngle(:) )
   return
end subroutine setKlist_FOURIER3D
!###################################################################
!# SUBROUTINE: setParityList_FOURIER3D
!###################################################################
!> @brief
!! Common set subroutine. Sets parityList atribute of a FOURIER3D object
!-------------------------------------------------------------------
subroutine setParityList_FOURIER3D(this,parityListXY,parityListAngle)
   implicit none
   ! I/O variables
   class(Fourier3d),intent(inout):: this
   character(len=1),dimension(:),intent(in):: parityListXY
   character(len=1),dimension(:),intent(in):: parityListAngle
   ! Local variables
   integer(kind=4):: n
   ! Run section
   n=size(parityListXY)
   allocate( this%parityListXY(n),    source=parityListXY(:)    )
   allocate( this%parityListAngle(n), source=parityListAngle(:) )
   return
end subroutine setParityList_FOURIER3D
!###################################################################
!# SUBROUTINE: setIrrepList_FOURIER3D
!###################################################################
!> @brief
!! Common set subroutine. Sets IrrepList atribute of a FOURIER3D object
!-------------------------------------------------------------------
subroutine setIrrepList_FOURIER3D(this,irrepListXY,irrepListAngle)
   implicit none
   ! I/O variables
   class(Fourier3d),intent(inout):: this
   character(len=2),dimension(:),intent(in):: irrepListXY
   character(len=2),dimension(:),intent(in):: irrepListAngle
   ! Local variables
   integer(kind=4):: n
   ! Run section
   n=size( irrepListXY )
   allocate( this%irrepListXY(n),    source=irrepListXY(:)    )
   allocate( this%irrepListAngle(n), source=irrepListAngle(:) )
   return
end subroutine setIrrepList_FOURIER3D
!###########################################################
!# SUBROUTINE: GET_ALLFUNC_AND_DERIVS_FOURIER3D
!###########################################################
!> @brief
!! Get value of the potential and derivs for an specific point x
!! for the main function and extra ones
!-----------------------------------------------------------
subroutine GET_ALLFUNC_AND_DERIVS_FOURIER3D(this,x,f,dfdx)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier3d),intent(in) :: this
   real(kind=8),dimension(3),intent(in) :: x
   real(kind=8),dimension(:),intent(out) :: f
   real(kind=8),dimension(:,:),intent(out):: dfdx
   ! Local variables
   integer(kind=4) :: nfuncs
   integer(kind=4) :: i ! counters
   real(kind=8),dimension(:),allocatable :: terms
   real(kind=8),dimension(:),allocatable :: terms_dx
   real(kind=8),dimension(:),allocatable :: terms_dy
   real(kind=8),dimension(:),allocatable :: terms_dz
   ! Check section
   select case(allocated(this%extracoeff))
      case(.false.)
         write(0,*) "GET_ALLFUNCS_AND_DERIVS ERR: extra coefficients are not allocated"
         write(0,*) "GET_ALLFUNCS_AND_DERIVS ERR: did you use ADD_MOREFUNCS and INTERPOL before this?"
         call EXIT(1)
      case(.true.)
         ! do nothing
   end select
   nfuncs=size(this%extrafuncs(:,1))+1
   select case( size(f)/=nfuncs .or. size(dfdx(:,1))/=nfuncs .or. size(dfdx(1,:))/=3 )
      case(.true.)
         write(0,*) "GET_ALLFUNCS_AND DERIVS ERR: size mismatch of output arguments"
         write(0,*) "nfuncs: ",nfuncs
         write(0,*) "size f: ", size(f)
         write(0,*) "size dfdx: ",size(dfdx(:,1)),size(dfdx(1,:))
         call EXIT(1)
      case(.false.)
         ! do nothing
   end select
   ! Run section
   allocate(terms(this%n))
   allocate(terms_dx(this%n))
   allocate(terms_dy(this%n))
   allocate(terms_dz(this%n))
   do i = 1, this%n
      terms(i)=this%term%getValue( k=this%kListXY(i,1:2), parityXY=this%parityListXY(i),       irrepXY=this%irrepListXY(i),&
                                   l=this%kListAngle(i),  parityAngle=this%parityListAngle(i), irrepAngle=this%irrepListAngle(i),&
                                   x=x )
      terms_dx(i)=this%term%getDeriv1( k=this%kListXY(i,1:2), parityXY=this%parityListXY(i),       irrepXY=this%irrepListXY(i),&
                                       l=this%kListAngle(i),  parityAngle=this%parityListAngle(i), irrepAngle=this%irrepListAngle(i),&
                                       x=x )
      terms_dy(i)=this%term%getDeriv2( k=this%kListXY(i,1:2), parityXY=this%parityListXY(i),       irrepXY=this%irrepListXY(i),&
                                       l=this%kListAngle(i),  parityAngle=this%parityListAngle(i), irrepAngle=this%irrepListAngle(i),&
                                       x=x )
      terms_dz(i)=this%term%getDeriv3( k=this%kListXY(i,1:2), parityXY=this%parityListXY(i),       irrepXY=this%irrepListXY(i),&
                                       l=this%kListAngle(i),  parityAngle=this%parityListAngle(i), irrepAngle=this%irrepListAngle(i),&
                                       x=x )
   end do
   f(1)=dot_product(terms,this%coeff)
   dfdx(1,1)=dot_product(terms_dx,this%coeff)
   dfdx(1,2)=dot_product(terms_dy,this%coeff)
   dfdx(1,3)=dot_product(terms_dz,this%coeff)
   do i = 2, nfuncs
      f(i)=dot_product(terms,this%extracoeff(:,i-1))
      dfdx(i,1)=dot_product(terms_dx,this%extracoeff(:,i-1))
      dfdx(i,1)=dot_product(terms_dy,this%extracoeff(:,i-1))
      dfdx(i,1)=dot_product(terms_dz,this%extracoeff(:,i-1))
   end do
   return
end subroutine GET_ALLFUNC_AND_DERIVS_FOURIER3D
!###########################################################
!# SUBROUTINE: ADD_MORE_FUNCS_FOURIER3D
!###########################################################
!> @brief
!! Adds a new set of functions to interpolate at the same time
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine ADD_MORE_FUNCS_FOURIER3D(this,f)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(FOURIER3D),intent(inout) :: this
   real(kind=8),dimension(:,:),intent(in) :: f
   ! Local variables
   integer(kind=4) :: nfuncs, ndata
   ! Run section
   nfuncs=size(f(:,1)) ! number of rows
   ndata=size(f(1,:)) ! number of columns
   select case(ndata == this%n)
      case(.false.)
         write(0,*) "ADD_MORE_FUNCS_FOURIER3D ERR: size mismatch between extra functions and the original one"
         call EXIT(1)
      case(.true.)
         ! donothing
   end select
   allocate(this%extrafuncs(nfuncs,ndata))
   this%extrafuncs=f
   return
end subroutine ADD_MORE_FUNCS_FOURIER3D
!###########################################################
!# SUBROUTINE: INTERPOL_FOURIER3D
!###########################################################
!> @brief
!! Performs a generic FOURIER3D interpolation
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Mar/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine INTERPOL_FOURIER3D(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(FOURIER3D),intent(inout) :: this
   ! Local variables
   real(kind=8),dimension(:,:),allocatable :: terms,inv_terms
   integer(kind=4) :: i,j ! counters
   ! Run section
   allocate(this%coeff(this%n))
   allocate(terms(this%n,this%n))
   allocate(inv_terms(this%n,this%n))
   do i = 1, this%n ! loop over eq for different points
      do j = 1, this%n ! loop over coefficients
         terms(i,j)=this%term%getValue( k=this%kListXY(j,1:2), parityXY=this%parityListXY(j),       irrepXY=this%irrepListXY(j),&
                                        l=this%kListAngle(j),  parityAngle=this%parityListAngle(j), irrepAngle=this%irrepListAngle(j),&
                                        x=this%x(i,:) )
      end do
   end do
   call INV_MTRX(this%n,terms,inv_terms)
   this%coeff=matmul(inv_terms,this%f)
   deallocate(terms)
   ! Check if there are extra functions to be interpolated
   select case(allocated(this%extrafuncs))
      case(.true.)
         allocate(this%extracoeff(this%n,size(this%extrafuncs(:,1))))
         do i = 1, size(this%extrafuncs(:,1))
            this%extracoeff(:,i)=matmul(inv_terms,this%extrafuncs(i,:))
         end do
      case(.false.)
         ! do nothing
   end select
   return
end subroutine INTERPOL_FOURIER3D
!###########################################################
!# FUNCTION: getvalue_FOURIER3D
!###########################################################
!> @brief
!! Get's fourier 3D interpolation value for a given point
!! inside the correct range
!-----------------------------------------------------------
function getValue_FOURIER3D(this,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier3d),target,intent(in) :: this
   real(kind=8),dimension(3),intent(in) :: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Local variables
   integer(kind=4):: i ! counters
   real(kind=8),dimension(:),allocatable:: terms
   ! Run section
   allocate(terms(this%n))
   do i = 1, this%n
      terms(i)=this%term%getValue( k=this%kListXY(i,1:2), parityXY=this%parityListXY(i),       irrepXY=this%irrepListXY(i),&
                                   l=this%kListAngle(i),  parityAngle=this%parityListAngle(i), irrepAngle=this%irrepListAngle(i),&
                                   x=x )
   end do
   answer=dot_product(terms,this%coeff)
   deallocate(terms)
   return
end function getValue_FOURIER3D
!###########################################################
!# FUNCTION: getDeriv1_FOURIER3D
!###########################################################
!> @brief
!! Get's derivative value at a given point X using the interpolation
!-----------------------------------------------------------
function getDeriv1_FOURIER3D(this,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier3d),target,intent(in) :: this
   real(kind=8),dimension(3),intent(in) :: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Local variables
   integer(kind=4):: i ! counters
   real(kind=8),dimension(:),allocatable:: terms
   ! Run section
   allocate(terms(this%n))
   do i = 1, this%n
      terms(i)=this%term%getDeriv1( k=this%kListXY(i,1:2), parityXY=this%parityListXY(i),       irrepXY=this%irrepListXY(i),&
                                    l=this%kListAngle(i),  parityAngle=this%parityListAngle(i), irrepAngle=this%irrepListAngle(i),&
                                    x=x )
   end do
   answer=dot_product(terms,this%coeff)
   deallocate(terms)
   return
end function getDeriv1_FOURIER3D
!###########################################################
!# FUNCTION: getDeriv2_FOURIER3D
!###########################################################
function getDeriv2_FOURIER3D(this,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier3d),target,intent(in) :: this
   real(kind=8),dimension(3),intent(in) :: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Local variables
   integer(kind=4):: i ! counters
   real(kind=8),dimension(:),allocatable:: terms
   ! Run section
   allocate(terms(this%n))
   do i = 1, this%n
      terms(i)=this%term%getDeriv1( k=this%kListXY(i,1:2), parityXY=this%parityListXY(i),  irrepXY=this%irrepListXY(i),&
                                    l=this%kListAngle(i),  parityAngle=this%parityListAngle(i), irrepAngle=this%irrepListAngle(i),&
                                    x=x )
   end do
   answer=dot_product(terms,this%coeff)
   deallocate(terms)
   return
end function getDeriv2_FOURIER3D
!###########################################################
!# FUNCTION: getDeriv3_FOURIER3D
!###########################################################
function getDeriv3_FOURIER3D(this,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier3d),target,intent(in) :: this
   real(kind=8),dimension(3),intent(in) :: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Local variables
   integer(kind=4):: i ! counters
   real(kind=8),dimension(:),allocatable:: terms
   ! Run section
   allocate(terms(this%n))
   do i = 1, this%n
      terms(i)=this%term%getDeriv1( k=this%kListXY(i,1:2), parityXY=this%parityListXY(i),       irrepXY=this%irrepListXY(i),&
                                    l=this%kListAngle(i),  parityAngle=this%parityListAngle(i), irrepAngle=this%irrepListAngle(i),&
                                    x=x )
   end do
   answer=dot_product(terms,this%coeff)
   deallocate(terms)
   return
end function getDeriv3_FOURIER3D

end module FOURIER3D_MOD

!########################################################
! MODULE : FOURIER3D_P4MM_MOD
!
!> @brief
!! Provides tools to perform 1D periodical interpolations with
!! symmetry discrimination.
!
!> @warning
!! - Includes FOURIER1D_mod in its scope
!########################################################
module FOURIER3D_P4MM_MOD
use FOURIER3D_MOD, only: TermCalculator3d,Fourier3d
use FOURIER1D_4MM_MOD, only: Term_4mm
use FOURIER_P4MM_MOD, only: TermCalculator2d_p4mm
implicit none
!/////////////////////////////////////////////////////////////////
! TYPE: term_A1
!> @brief
!! Child class of abstract termcalculator. Structure to avoid unnecessary switches
!----------------------------------------------------------------
type,extends(TermCalculator3d) :: Term3d_p4mm
private
   ! some atributes
contains
   procedure,public:: getValue  => termFou3d_p4mm
   procedure,public:: getDeriv1 => termFou3d_dx_p4mm
   procedure,public:: getDeriv2 => termFou3d_dy_p4mm
   procedure,public:: getDeriv3 => termFou3d_dz_p4mm
end type term3d_p4mm
!/////////////////////////////////////////////////
! TYPE: FOURIER3D_P4MM
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
type,extends(Fourier3d):: Fourier3d_p4mm
   contains
      ! Set block
      procedure,public:: initializeTerms => initializeTerms_FOURIER3D_P4MM
end type Fourier3d_p4mm
!//////////////////////////////////////////////////
contains
!###########################################################
!# SUBROUTINE: SET_IRREP_FOURIER3D_P4MM
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!-----------------------------------------------------------
subroutine initializeTerms_FOURIER3D_P4MM(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Fourier3d_p4mm),intent(inout):: this
   ! Run section
   allocate( Term3d_p4mm::   this%term )
   allocate( TermCalculator2d_p4mm::   this%term%xyFourier )
   allocate( Term_4mm:: this%term%angleFourier )
   return
end subroutine initializeTerms_FOURIER3D_P4MM
!###########################################################
!# FUNCTION: termfou1d_4mm
!###########################################################
!-----------------------------------------------------------
function termFou3d_p4mm(this,k,parityXY,irrepXY,l,parityAngle,irrepAngle,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term3d_p4mm),intent(in):: this
   integer(kind=4),dimension(2),intent(in):: k
   integer(kind=4),intent(in):: l
   character(len=1),intent(in):: parityXY
   character(len=1),intent(in):: parityAngle
   character(len=2),intent(in):: irrepXY
   character(len=2),intent(in):: irrepAngle
   real(kind=8),dimension(3),intent(in):: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou1d_4mm: '
   ! Run section
   answer=this%xyFourier%getValue( x=x(1:2),k=k,irrep=irrepXY,parity=parityXY )*&
          this%angleFourier%getValue( x=x(3),kpoint=l,irrep=irrepAngle,parity=parityAngle )
   return
end function termFou3d_p4mm
!###########################################################
!# FUNCTION: termfou1d_dx_4mm_A1
!###########################################################
!-----------------------------------------------------------
function termFou3d_dx_p4mm(this,k,parityXY,irrepXY,l,parityAngle,irrepAngle,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term3d_p4mm),intent(in):: this
   integer(kind=4),dimension(2),intent(in):: k
   integer(kind=4),intent(in):: l
   character(len=1),intent(in):: parityXY
   character(len=1),intent(in):: parityAngle
   character(len=2),intent(in):: irrepXY
   character(len=2),intent(in):: irrepAngle
   real(kind=8),dimension(3),intent(in):: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou1d_4mm: '
   ! Run section
   answer=this%xyFourier%getDeriv1( x=x(1:2),k=k,irrep=irrepXY,parity=parityXY )*&
          this%angleFourier%getValue( x=x(3),kpoint=l,irrep=irrepAngle,parity=parityAngle )
   return
end function termFou3d_dx_p4mm
!###########################################################
!# FUNCTION: termfou1d_dy_4mm_A1
!###########################################################
!-----------------------------------------------------------
function termFou3d_dy_p4mm(this,k,parityXY,irrepXY,l,parityAngle,irrepAngle,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term3d_p4mm),intent(in):: this
   integer(kind=4),dimension(2),intent(in):: k
   integer(kind=4),intent(in):: l
   character(len=1),intent(in):: parityXY
   character(len=1),intent(in):: parityAngle
   character(len=2),intent(in):: irrepXY
   character(len=2),intent(in):: irrepAngle
   real(kind=8),dimension(3),intent(in):: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou1d_4mm: '
   ! Run section
   answer=this%xyFourier%getDeriv2( x=x(1:2),k=k,irrep=irrepXY,parity=parityXY )*&
          this%angleFourier%getValue( x=x(3),kpoint=l,irrep=irrepAngle,parity=parityAngle )
   return
end function termFou3d_dy_p4mm
!###########################################################
!# FUNCTION: termfou1d_dz_4mm_A1
!###########################################################
!-----------------------------------------------------------
function termFou3d_dz_p4mm(this,k,parityXY,irrepXY,l,parityAngle,irrepAngle,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Term3d_p4mm),intent(in):: this
   integer(kind=4),dimension(2),intent(in):: k
   integer(kind=4),intent(in):: l
   character(len=1),intent(in):: parityXY
   character(len=1),intent(in):: parityAngle
   character(len=2),intent(in):: irrepXY
   character(len=2),intent(in):: irrepAngle
   real(kind=8),dimension(3),intent(in):: x
   ! Dummy output variable
   real(kind=8):: answer
   ! Parameters
   character(len=*),parameter:: routinename='termfou1d_4mm: '
   ! Run section
   answer=this%xyFourier%getValue( x=x(1:2),k=k,irrep=irrepXY,parity=parityXY )*&
          this%angleFourier%getDeriv( x=x(3),kpoint=l,irrep=irrepAngle,parity=parityAngle )
   return
end function termFou3d_dz_p4mm

end module FOURIER3D_P4MM_MOD

!#########################################################
! MODULE: INTERPOL_WYCKOFF_GENERIC
!> @brief
!! Provides generic routines to interpol a CRP6D PES
!##########################################################
MODULE WYCKOFF_GENERIC_MOD
   use BICSPLINES_MOD
   use FOURIER1D_MOD
   use MATHS_MOD
   use UNITS_MOD
#ifdef DEBUG
   use DEBUG_MOD
#endif
! Initial declarations
IMPLICIT NONE
!//////////////////////////////////////////////////
! TYPE: Fouklist
!> @brief
!! Auxiliar type to get rid of some technical problems
!--------------------------------------------------
TYPE Fouklist
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE:: k
END TYPE
!//////////////////////////////////////////////////////////////
! TTPE: TermsInfo
!> @brief
!! Used to store information of terms in a phiCut or a thetaCut
!--------------------------------------------------------------
type TermsInfo
   character(len=2),dimension(:),allocatable,public:: irrepList
   character(len=1),dimension(:),allocatable,public:: parityList
   integer(kind=4),dimension(:),allocatable,public:: kpointList
   contains
      procedure,public:: addTerm => addTerm_TERMSINFO
      procedure,public:: reboot  => reboot_TERMSINFO
end type
!--------------------------------------------------
!//////////////////////////////////////////////////
! TYPE: Cut2d
!> @brief
!! Contains data for a Z,R interpolation letting X,Y,Theta,Phi
!! fixed
!-------------------------------------------------
TYPE Cut2d
   CHARACTER(LEN=:),ALLOCATABLE:: alias
   CHARACTER(LEN=:),ALLOCATABLE:: filename
   REAL(KIND=8):: x,y,phi,theta
   TYPE(Bicsplines):: interrz
   CONTAINS
      ! Initialize block
      PROCEDURE,PUBLIC:: READ => READ_CUT2D
      ! Get block
      PROCEDURE,PUBLIC:: gettheta => gettheta_cut2d
      PROCEDURE,PUBLIC:: getphi => getphi_cut2d
      PROCEDURE,PUBLIC:: getgridsizer => getgridsizer_cut2d
      PROCEDURE,PUBLIC:: getgridsizez => getgridsizez_cut2d
      PROCEDURE,PUBLIC:: getfirstr => getfirstr_cut2d
      PROCEDURE,PUBLIC:: getfirstz => getfirstz_cut2d
      PROCEDURE,PUBLIC:: getlastz => getlastz_cut2d
      PROCEDURE,PUBLIC:: getlastr => getlastr_cut2d
      PROCEDURE,PUBLIC:: getgridvaluer => getgridvaluer_cut2d
      PROCEDURE,PUBLIC:: getgridvaluez => getgridvaluez_cut2d
      PROCEDURE,PUBLIC:: getalias => getalias_cut2d
      PROCEDURE,PUBLIC:: getpotatgridpoint => getpotatgridpoint_cut2d
      ! Tools block
      PROCEDURE,PUBLIC:: PRINT_INPUT => PRINT_INPUT_CUT2D
      PROCEDURE,PUBLIC:: INTERPOL => INTERPOL_CUT2D
      PROCEDURE,PUBLIC:: CHANGEPOT_AT_GRIDPOINT => CHANGEPOT_AT_GRIDPOINT_CUT2D
END TYPE Cut2d
!////////////////////////////////////////////////////////////////////////////////////////////////////
! TYPE: WYCKOFFSITIO
!
!> @brief
!! Stores all data related to a single X,Y position
!! of the molecule on the surface.
!
!> @param id - Letter that specifies the symmetry subgroup
!!             related to this site. Its meaning changes for
!!             each wallpaper group
!> @param mynumber - Number that identifies this wyckoff site
!> @param x,y - Position in XY plane of this wyckoffsite
!> @param is_homonuclear - True if the molecule is homonuclear
!> @param  n2dcuts - Number of Z,R cuts for this site
!> @param nphicuts - Number of Phi interpolations that will exist
!!                   for different Theta values.
!
!> @brief nphipoints - Array of integer numbers storing the number of points
!!                     in which each phi interpolation is based
!>
!---------------------------------------------------------------------------------------------------
type,abstract:: WyckoffSitio
   ! public atributes
   character(len=1),public:: id
   integer(kind=4),public:: mynumber
   logical,public:: is_homonucl=.FALSE.
   real(kind=8),public:: x,y
   integer(kind=4),public:: n2dcuts
   integer(kind=4),public:: nphicuts
   integer(kind=4),dimension(:),allocatable,public:: nphipoints
   type(Cut2d),dimension(:),allocatable,public:: zrcut
   type(TermsInfo),dimension(:),allocatable:: phiTerms
   type(TermsInfo):: thetaTerms
   character(len=2),dimension(:,:),allocatable,public:: irrepList
   character(len=1),dimension(:,:),allocatable,public:: parityList
   integer(kind=4),dimension(:,:),allocatable,public:: kpointList
   real(kind=8),dimension(:),allocatable,public:: phiList
   contains
      ! Initialization block
      procedure,public,non_overridable:: INITIALIZE => INITIALIZE_WYCKOFFSITIO
      procedure,public,non_overridable:: SET_ID => SET_ID_WYCKOFFSITIO
      procedure,public,non_overridable:: SET_HOMONUCL => SET_HOMONUCL_WYCKOFFSITIO
      procedure,public,non_overridable:: SET_MYNUMBER => SET_MYNUMBER_WYCKOFFSITIO
      ! interface procedures
      procedure(get_v_and_derivs_WYCKOFFSITIO),public,deferred:: GET_V_AND_DERIVS
end type WyckoffSitio

ABSTRACT INTERFACE
   !###########################################################
   !# SUBROUTINE: GET_V_AND_DERIVS_WYCKOFFSITIO
   !###########################################################
   !-----------------------------------------------------------
   SUBROUTINE GET_V_AND_DERIVS_WYCKOFFSITIO(this,x,v,dvdu)
      IMPORT Wyckoffsitio
      CLASS(Wyckoffsitio),INTENT(IN) :: this
      REAL(KIND=8),DIMENSION(4),INTENT(IN) ::x
      REAL(KIND=8),INTENT(OUT) :: v
      REAL(KIND=8),DIMENSION(4),INTENT(OUT) :: dvdu
   END SUBROUTINE GET_V_AND_DERIVS_WYCKOFFSITIO
END INTERFACE
!//////////////////////////////////////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: PRINT_INPUT_CUT2D
!###########################################################
!> @brief
!! Prints an input file from data stored in RAM
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE PRINT_INPUT_CUT2D(this,filename)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables
   INTEGER(KIND=4) :: nr,nz
   INTEGER(KIND=4) :: i,j ! counters
   ! Run section
   nr=this%getgridsizeR()
   nz=this%getgridsizeZ()
   OPEN (10,FILE=filename,STATUS="replace",ACTION="write")
   WRITE(10,*) "# FILE generated by PRINT_CUT2D_INPUT CRP6D"
   WRITE(10,*) this%alias
   WRITE(10,*) "au au rad"
   WRITE(10,*) this%x,this%y
   WRITE(10,*) this%theta,this%phi
   WRITE(10,*) nr,nz
   DO j = 1, nz
      DO i = 1, nr
         WRITE(10,*) this%getgridvalueR(i),this%getgridvalueZ(j),this%getpotatgridpoint(i,j)
      END DO
   END DO
   CLOSE(10)
   RETURN
END SUBROUTINE PRINT_INPUT_CUT2D
!###########################################################
!# FUNCTION: gettheta_cut2d
!###########################################################
!> @brief
!! Common get functions. Gets theta atribute
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION gettheta_cut2d(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Run section
   gettheta_cut2d=this%theta
   RETURN
END FUNCTION gettheta_cut2d
!###########################################################
!# FUNCTION: getphi_cut2d
!###########################################################
!> @brief
!! Common get functions. Gets phi atribute
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getphi_cut2d(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Run section
   getphi_cut2d=this%phi
   RETURN
END FUNCTION getphi_cut2d
!###########################################################
!# SUBROUTINE: SET_MYNUMBER_WYCKOFFSITIO
!###########################################################
!> @brief
!! Common set function. Sets_mynumber atribute
!-----------------------------------------------------------
SUBROUTINE SET_MYNUMBER_WYCKOFFSITIO(this,mynumber)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Wyckoffsitio),INTENT(INOUT):: this
   INTEGER(KIND=4),INTENT(IN) :: mynumber
   ! Run section
   this%mynumber=mynumber
   RETURN
END SUBROUTINE SET_MYNUMBER_WYCKOFFSITIO
!###########################################################
!# SUBROUTINE: SET_HOMONUCL_WYCKOFFSITIO
!###########################################################
!> @brief
!! Common set function. Sets homonucl atribute
!-----------------------------------------------------------
SUBROUTINE SET_HOMONUCL_WYCKOFFSITIO(this,is_homonucl)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Wyckoffsitio),INTENT(INOUT):: this
   LOGICAL,INTENT(IN) :: is_homonucl
   ! Run section
   this%is_homonucl=is_homonucl
   RETURN
END SUBROUTINE SET_HOMONUCL_WYCKOFFSITIO
!###########################################################
!# FUNCTION: getpotatgridpoint_cut2d
!###########################################################
!> @brief
!! Common get function. Gets the potential at a given grid point
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getpotatgridpoint_cut2d(this,i,j)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: i,j
   ! Run section
   getpotatgridpoint_cut2d=this%interrz%fgrid(i,j)
   RETURN
END FUNCTION getpotatgridpoint_cut2d
!###########################################################
!# SUBROUTINE: CHANGEPOT_AT_GRIDPOINT_CUT2D
!###########################################################
!> @brief
!! Changes the value of the potential at a given point of the
!! grid.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 26/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE CHANGEPOT_AT_GRIDPOINT_CUT2D(this,i,j,newpot)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(INOUT):: this
   INTEGER(KIND=4),INTENT(IN) :: i,j
   REAL(KIND=8),INTENT(IN) :: newpot
   ! Run section
   this%interrz%fgrid(i,j)=newpot
   RETURN
END SUBROUTINE CHANGEPOT_AT_GRIDPOINT_CUT2D
!###########################################################
!# FUNCTION: getgridvaluer_cut2d
!###########################################################
!> @brief
!! Common get function. Gets ith component of the grid in R
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getgridvaluer_cut2d(this,i)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   INTEGER(KIND=4),INTENT(IN) :: i
   ! Run section
   getgridvaluer_cut2d=this%interrz%x(i)
   RETURN
END FUNCTION getgridvaluer_cut2d
!###########################################################
!# FUNCTION: getgridvaluez_cut2d
!###########################################################
!> @brief
!! Common get function. Gets ith component of the grid in Z
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getgridvaluez_cut2d(this,i)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   INTEGER(KIND=4),INTENT(IN) :: i
   ! Run section
   getgridvaluez_cut2d=this%interrz%y(i)
   RETURN
END FUNCTION getgridvaluez_cut2d
!###########################################################
!# FUNCTION: getalias_cut2d
!###########################################################
!> @brief
!! Common get function. Gets alias of a Cut2d object
!-----------------------------------------------------------
CHARACTER(LEN=30) FUNCTION getalias_cut2d(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Run section
   getalias_cut2d=this%alias
   RETURN
END FUNCTION getalias_cut2d
!###########################################################
!# FUNCTION: getgridsizer_cut2d
!###########################################################
!> @brief
!! Common get function. Gets R grid size
!-----------------------------------------------------------
INTEGER(KIND=4) FUNCTION getgridsizer_cut2d(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Run section
   getgridsizer_cut2d=size(this%interrz%x)
   RETURN
END FUNCTION getgridsizer_cut2d
!###########################################################
!# FUNCTION: getgridsizez_cut2d
!###########################################################
!> @brief
!! Common get function. Gets Z grid size
!-----------------------------------------------------------
INTEGER(KIND=4) FUNCTION getgridsizez_cut2d(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Run section
   getgridsizez_cut2d=size(this%interrz%y)
   RETURN
END FUNCTION getgridsizez_cut2d
!###########################################################
!# FUNCTION: getfirstr_cut2d
!###########################################################
!> @brief
!! Just common get function to get first R value
!
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getfirstr_cut2d(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Local variables

   ! Run section
   getfirstr_cut2d=this%interrz%x(1)
   RETURN
END FUNCTION getfirstr_cut2d
!###########################################################
!# FUNCTION: getlastr_cut2d
!###########################################################
!> @brief
!! Just common get function to get last R value
!
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getlastr_cut2d(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Run section
   getlastr_cut2d=this%interrz%x(size(this%interrz%x))
   RETURN
END FUNCTION getlastr_cut2d
!###########################################################
!# FUNCTION: getfirstz_cut2d
!###########################################################
!> @brief
!! Just common get function to get first Z value
!
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getfirstz_cut2d(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Run section
   getfirstz_cut2d=this%interrz%y(1)
   RETURN
END FUNCTION getfirstz_cut2d
!###########################################################
!# FUNCTION: getlastz_cut2d
!###########################################################
!> @brief
!! Just common get function to get last Z value
!
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getlastz_cut2d(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Run section
   getlastz_cut2d=this%interrz%y(size(this%interrz%y))
   RETURN
END FUNCTION getlastz_cut2d
!###########################################################
!# SUBROUTINE: SET_ID_WYCKOFFSITIO
!###########################################################
!> @brief
!! Set the correct id (letter) for this wyckoff sitio
!-----------------------------------------------------------
SUBROUTINE SET_ID_WYCKOFFSITIO(this,id)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Wyckoffsitio),INTENT(INOUT):: this
   CHARACTER,INTENT(IN):: id
   ! Run section
   this%id=id
   RETURN
END SUBROUTINE SET_ID_WYCKOFFSITIO
!#####################################################
! SUBROUTINE: READ_CUT2D
!> @brief
!! Set up a Cut2d object from file
!-----------------------------------------------------
SUBROUTINE READ_CUT2D(this,filename)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(OUT):: this
   CHARACTER(LEN=*),INTENT(IN):: filename
   ! Local variables
   INTEGER(KIND=4):: i,j ! counters
   INTEGER(KIND=4):: nx,ny
   REAL(KIND=8):: auxr1,auxr2
   CHARACTER(LEN=10):: units1,units2,units3
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: z,r
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: f
   TYPE(Length):: len1,len2
   TYPE(Energy):: en
   TYPE(Angle):: angl
   ! Run section ----------------------
   this%filename = filename
   OPEN (UNIT=10,FILE=filename,STATUS="old",ACTION="read")
   READ(10,*)
   READ(10,*) this%alias
   READ(10,*) units1,units2,units3
   READ(10,*) auxr1,auxr2
   CALL len1%READ(auxr1,units1)
   CALL len2%READ(auxr2,units1)
   CALL len1%TO_STD()
   CALL len2%TO_STD()
   this%x=len1%getvalue()
   this%y=len2%getvalue()
   READ(10,*) auxr1,auxr2
   CALL angl%READ(auxr1,units3)
   CALL angl%TO_STD()
   this%theta=angl%getvalue()
   CALL angl%READ(auxr2,units3)
   CALL angl%TO_STD()
   this%phi=angl%getvalue()
   READ(10,*) nx,ny
   ALLOCATE(r(nx))
   ALLOCATE(z(ny))
   ALLOCATE(f(nx,ny))
   DO j = 1, ny
      DO i = 1, nx
         READ(10,*) r(i),z(j),f(i,j)
         CALL len1%READ(r(i),units1)
         CALL len2%READ(z(j),units1)
         CALL en%READ(f(i,j),units2)
         CALL len1%TO_STD()
         CALL len2%TO_STD()
         CALL en%TO_STD()
         r(i)=len1%getvalue()
         z(j)=len2%getvalue()
         f(i,j)=en%getvalue()
      END DO
   END DO
   ! order array in a good way:
   DO j = 1, ny
      CALL ORDER(r,f(:,j))
   END DO
   DO i = 1, nx
      CALL ORDER(z,f(i,:))
   END DO
   !
   CALL this%interrz%READ(r,z,f)
   CLOSE(10)
   RETURN
END SUBROUTINE READ_CUT2D
!###########################################################
!# SUBROUTINE: INTERPOl_CUT2D
!###########################################################
!> @brief
!! Interpolates a RZ-2dcur of the potential
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 25/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INTERPOl_CUT2D(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(INOUT) :: this
   ! Local variables
   ! Run section
   CALL this%interrz%INTERPOL()
   RETURN
END SUBROUTINE INTERPOl_CUT2D
!###########################################################
!# SUBROUTINE: INITIALIZE_WYCKOFFSITIO
!###########################################################
!> @brief
!! Initialize a Wyckoffsitio from file.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 20/03/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine initialize_WYCKOFFSITIO(this,nphiPoints,fileNames,letter,myNumber,phiTerms,thetaTerms,phiList)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Wyckoffsitio),intent(inout):: this
   integer(kind=4),dimension(:),intent(in):: nphipoints
   character(len=*),dimension(:),intent(in):: filenames
   character,intent(in):: letter
   integer(kind=4),intent(in):: mynumber
   type(TermsInfo),dimension(:),intent(in):: phiTerms
   type(TermsInfo),intent(in):: thetaTerms
   real(kind=8),dimension(:),intent(in):: phiList
   ! Local variables
   integer(kind=4):: i ! counters
   ! Parameters
   character(len=*),parameter:: routinename="INITIALIZE_WYCKOFFSITIO: "
   ! Run section
   this%mynumber=mynumber
   this%id=letter
   this%n2dcuts=sum(nphipoints(:))
   this%nphicuts=size(nphipoints(:))
   allocate( this%zrcut(this%n2dcuts) )
   allocate( this%nphipoints(this%nphicuts),source=nphipoints(:) )
   allocate( this%phiTerms(size(phiTerms(:))),source=phiTerms(:) )
   allocate( this%phiList(size(phiList)),source=phiList(:) )
   this%thetaTerms=thetaTerms
   SELECT CASE(this%n2dcuts==size(filenames(:)))
      CASE(.true.)
         ! do nothing
      CASE(.false.)
         WRITE(0,*) "INITIALIZE_WYCKOFFSITIO ERR: mismatch between number of n2dcuts and number of files to open"
         CALL EXIT(1)
   END SELECT
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"n2dcuts: ",this%n2dcuts)
   CALL VERBOSE_WRITE(routinename,"nphicuts: ",this%nphicuts)
   CALL VERBOSE_WRITE(routinename,'theta structure: ',this%nphipoints(:))
   DO i = 1, this%n2dcuts
      CALL VERBOSE_WRITE(routinename,trim(filenames(i)))
   END DO
#endif
   DO i = 1, this%n2dcuts
      CALL this%zrcut(i)%READ(trim(filenames(i)))
   END DO
   ! All zrcuts should belong to the same XY position (center of mass)
   this%x=this%zrcut(1)%x
   this%y=this%zrcut(1)%y
   ! DEBUGGING PART
   RETURN
END SUBROUTINE INITIALIZE_WYCKOFFSITIO
!############################################################
! SUBROUTINE: addTerm_TERMSINFO
!############################################################
!> @brief
!! Adds a new term to the already existing ones
!------------------------------------------------------------
subroutine addTerm_TERMSINFO(this,irrep,parity,kpoint)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(TermsInfo),intent(inout):: this
   character(len=2),intent(in):: irrep
   character(len=1),intent(in):: parity
   integer(kind=4),intent(in):: kpoint
   ! Local variables
   integer(kind=4):: n
   character(len=2),dimension(:),allocatable:: auxIrrepList
   character(len=1),dimension(:),allocatable:: auxParityList
   integer(kind=4),dimension(:),allocatable::  auxKpointList

   ! Run section
   select case( allocated(this%kpointList) )
   case(.true.)
      n=size( this%kpointList(:) )
      call move_alloc( from=this%kpointList,to=auxKpointList )
      call move_alloc( from=this%irrepList,to=auxIrrepList )
      call move_alloc( from=this%parityList,to=auxParityList )
      allocate( this%kpointList(n+1) )
      allocate( this%irrepList(n+1) )
      allocate( this%ParityList(n+1) )
      this%kpointList(1:n)=auxKpointList(:)
      this%irrepList(1:n)=auxIrrepList(:)
      this%parityList(1:n)=auxParityList(:)
      this%kpointList(n+1)=kpoint
      this%irrepList(n+1)=trim(irrep)
      this%parityList(n+1)=trim(parity)

   case(.false.)
      allocate( this%kpointList(1) )
      allocate( this%parityList(1) )
      allocate( this%irrepList(1) )
      this%kpointList(1)=kpoint
      this%parityList(1)=trim(parity)
      this%irrepList(1)=trim(irrep)
   end select
   return
end subroutine addTerm_TERMSINFO
!########################################################
! SUBROUTINE: reboot_TERMSINFO
!########################################################
!> @brief
!! Deallocates all atributes of TERMSINFO
!--------------------------------------------------------
subroutine reboot_TERMSINFO(this)
   ! initial declaration
   implicit none
   ! I/O variables
   class(TermsInfo),intent(inout):: this
   ! GET TO THE CHOPPAHH !!!
   deallocate(this%kpointList)
   deallocate(this%parityList)
   deallocate(this%irrepList)
   return
end subroutine reboot_TERMSINFO

END MODULE WYCKOFF_GENERIC_MOD

!#########################################################
! MODULE: WYCKOFF_P4MM
!> @brief
!! Provides tools to interpolate through wyckoff sites belonging to
!! p4mm wallpaper symmetry
!##########################################################
module WYCKOFF_P4MM_MOD
! Initial declarations
use WYCKOFF_GENERIC_MOD
use FOURIER1D_2_MOD
use FOURIER1D_4MM_MOD
use FOURIER1D_M_MOD
use FOURIER1D_M45_MOD
use FOURIER1D_MM2_MOD
use FOURIER1D_M45M135_2_MOD
use FOURIER1D_E_MOD
#ifdef DEBUG
use DEBUG_MOD
use UNITS_MOD
#endif
implicit none
!/////////////////////////////////////////////////////////////////
! TYPE: Wyckoffp4mm
!> @brief
!! Subclass of Wyckoffsitio for generic p4mm symmetry.
!! in p4mm symmetry does not matter
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 20/03/2014
!> @version 1.0
!----------------------------------------------------------------
type,extends(Wyckoffsitio) :: Wyckoffp4mm
   contains
      procedure,public :: GET_V_AND_DERIVS => GET_V_AND_DERIVS_WYCKOFFP4MM
end type Wyckoffp4mm
!/////////////////////////////////////////////////////////////////
contains
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_WYCKOFFP4MM
!###########################################################
!> @brief
!! Symmetry adapted interpolation for p4mm wallpaper group.
!
!> @param[in] x - Array which stands for Z,r,theta,phi
!> @param[out] v - potential at X
!> @param[out] dvdu - Array with derivatives: dvdz, dvdr, dvdtheta dvdphi
!
!> @warning
!! - Only stands for a,b,c, and f p4mm Wyckoff sites.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.2
!-----------------------------------------------------------
subroutine GET_V_AND_DERIVS_WYCKOFFP4MM(this,x,v,dvdu)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Wyckoffp4mm),intent(in) :: this
   real(kind=8),dimension(4),intent(in) ::x
   real(kind=8),intent(out) :: v
   real(kind=8),dimension(4),intent(out) :: dvdu
   ! Local variables
   integer(kind=4) :: i,j,h ! counters
   real(kind=8) :: z,r,theta,phi
   class(Fourier1d),dimension(:),allocatable :: phicut
   class(Fourier1d),allocatable :: thetacut
   real(kind=8) :: aux_theta
   real(kind=8),dimension(:,:),allocatable :: aux
   real(kind=8),dimension(:),allocatable :: f
   real(kind=8),dimension(:),allocatable :: dfdz
   real(kind=8),dimension(:),allocatable :: dfdr
   real(kind=8),dimension(:),allocatable :: dfdphi
   real(kind=8),dimension(:),allocatable :: philist
   real(kind=8),dimension(:),allocatable :: thetalist
   character(len=2) :: theta_irrep, phi_irrep
   character(len=30),parameter :: routinename="GET_V_AND_DERIVS_WYCKOFFP4MM: "
#ifdef DEBUG
   type(Angle),dimension(:),allocatable :: beta
   real(kind=8),dimension(:),allocatable :: philistdeg
   character(len=18) :: filename
#endif
   ! Run section
   z=x(1)
   r=x(2)
   theta=x(3)
   phi=x(4)
   ! CONDITIONS: edit this part to include/edit symmetries
   select case( this%id )
      case("a" : "b")
         allocate( Fourier1d_4mm::phiCut(this%nphicuts) )
         allocate( Fourier1d_mm2::thetaCut )

      case("c")
         allocate( Fourier1d_mm2::phiCut(this%nphicuts) )
         allocate( Fourier1d_mm2::thetaCut )

      case("f")
         allocate( Fourier1d_m45::phiCut(this%nphicuts) )
         allocate( Fourier1d_2::thetaCut )

      case default
         write(0,*) "GET_V_AND_DERIVS_WYCKOFFP4MM: Unexpected error with Wyckoff id"
         call EXIT(1)
   end select
   ! PHI INTERPOLATION ----------------------------------------------------
   h=0 ! initialize h
   do i = 1, this%nphicuts ! loop over specific phi cuts (each one for a different theta value)
      allocate(f(this%nphipoints(i)))
      allocate(dfdr(this%nphipoints(i)))
      allocate(dfdz(this%nphipoints(i)))
      allocate(philist(this%nphipoints(i)))
      allocate(aux(2,this%nphipoints(i)))
      do j = 1, this%nphipoints(i) ! loop over number of zrcuts inside
         h=h+1 ! numbering of zrcuts
         f(j)=this%zrcut(h)%interrz%getvalue((/r,z/)) ! storing potential at this site
         dfdr(j)=this%zrcut(h)%interrz%getderivx((/r,z/)) ! storing d/dr at this site
         dfdz(j)=this%zrcut(h)%interrz%getderivy((/r,z/)) ! storing d/dz at this site
         philist(j)=this%zrcut(h)%phi
         aux_theta=this%zrcut(h)%theta
      end do
      aux(1,:)=dfdz(:)
      aux(2,:)=dfdr(:)
      call phiCut(i)%read(philist,f)
      call phiCut(i)%add_morefuncs(aux)
      call phiCut(i)%setKlist( this%phiTerms(i)%kpointList(:) )
      call phiCut(i)%setParityList( this%phiTerms(i)%parityList(:) )
      call phiCut(i)%setIrrepList( this%phiTerms(i)%irrepList(:) )
      call phiCut(i)%initializeTerms()
      call phiCut(i)%interpol()
#ifdef DEBUG
      call DEBUG_WRITE(routinename,"NEW PHICUT")
      call DEBUG_WRITE(routinename,"For theta: ",this%zrcut(h)%theta)
      call DEBUG_WRITE(routinename,"Is homonuclear: ",this%is_homonucl)
      allocate(beta(size(philist)))
      allocate(philistdeg(size(philist)))
      do j = 1, size(philist)
         call beta(j)%READ(philist(j),"rad")
         call beta(j)%TO_DEG()
         philistdeg(j)=beta(j)%getvalue()
      end do
      call DEBUG_WRITE(routinename,"At Phi: (deg)",philistdeg)
      call DEBUG_WRITE(routinename,"At Phi: (rad)",philist)
      call DEBUG_WRITE(routinename,"Klist: ",phicut(i)%getklist())
      call DEBUG_WRITE(routinename,"f: ",f)
      call DEBUG_WRITE(routinename,"dfdz: ",aux(1,:))
      call DEBUG_WRITE(routinename,"dfdr: ",aux(2,:))
      select case(get_debugmode())
         case(.true.)
            write(filename,'(I1,A1,I1,A1,A14)') this%mynumber,"-",i,"-","wyckoffphi.raw"
            call phicut(i)%PLOTDATA(filename)
            write(filename,'(I1,A1,I1,A1,A14)') this%mynumber,"-",i,"-","wyckoffphi.cyc"
            call phicut(i)%PLOTCYCLIC_ALL(300,filename)
         case(.false.)
            ! do nothing
      end select
      deallocate(beta)
      deallocate(philistdeg)
#endif
      deallocate(philist)
      deallocate(f)
      deallocate(dfdr)
      deallocate(dfdz)
      deallocate(aux)
   end do
   ! THETA INTERPOLATION --------------------------
   allocate(thetalist(this%nphicuts))
   allocate(f(this%nphicuts))
   allocate(dfdz(this%nphicuts))
   allocate(dfdr(this%nphicuts))
   allocate(dfdphi(this%nphicuts))
   h=0 ! reboot h
   do i = 1, this%nphicuts
      do j = 1, this%nphipoints(i)
         h=h+1
         thetalist(i)=this%zrcut(h)%theta
      end do
      allocate(aux(2,3))
      call phicut(i)%GET_ALLFUNCS_AND_DERIVS(phi,aux(1,:),aux(2,:))
      f(i)=aux(1,1)
      dfdz(i)=aux(1,2)
      dfdr(i)=aux(1,3)
      dfdphi(i)=aux(2,1)
      deallocate(aux)
   end do
   allocate(aux(3,this%nphicuts))
   aux(1,:)=dfdz(:)
   aux(2,:)=dfdr(:)
   aux(3,:)=dfdphi(:)
   call thetaCut%read(thetalist,f)
   call thetaCut%add_morefuncs(aux)
   call thetaCut%setKlist     ( this%thetaTerms%kpointList(:)    )
   call thetaCut%setParityList( this%thetaTerms%parityList(:) )
   call thetaCut%setIrrepList ( this%thetaTerms%irrepList(:)  )
   call thetaCut%initializeTerms()
   call thetaCut%interpol()
#ifdef DEBUG
   call DEBUG_WRITE(routinename,"NEW THETACUT")
   allocate(beta(size(thetalist)))
   allocate(philistdeg(size(thetalist)))
   do j = 1, size(thetalist)
      call beta(j)%READ(thetalist(j),"rad")
      call beta(j)%TO_DEG()
      philistdeg(j)=beta(j)%getvalue()
   end do
   call DEBUG_WRITE(routinename,"At Theta: (deg) ",philistdeg)
   call DEBUG_WRITE(routinename,"At Theta: (rad) ",thetalist)
   call DEBUG_WRITE(routinename,"Klist:          ",thetacut%getklist())
   call DEBUG_WRITE(routinename,"f:              ",f)
   call DEBUG_WRITE(routinename,"dfdz:           ",dfdz)
   call DEBUG_WRITE(routinename,"dfdr:           ",dfdr)
   call DEBUG_WRITE(routinename,"dfdphi:         ",dfdphi)
   select case(get_debugmode())
      case(.true.)
         write(filename,'(I1,A1,A16)') this%mynumber,"-","wyckofftheta.raw"
         call thetacut%PLOTDATA(filename)
         write(filename,'(I1,A1,A16)') this%mynumber,"-","wyckofftheta.cyc"
         call thetacut%PLOTCYCLIC(300,filename)
      case(.false.)
         ! do nothing
   end select
   deallocate(beta)
   deallocate(philistdeg)
#endif
   deallocate(f)
   deallocate(dfdr)
   deallocate(dfdz)
   deallocate(dfdphi)
   deallocate(aux)
   deallocate(thetalist)
   allocate(aux(2,4))
   call thetacut%GET_ALLFUNCS_AND_DERIVS(theta,aux(1,:),aux(2,:))
   v=aux(1,1) ! value of the potential
   dvdu(1)=aux(1,2) ! dvdz
   dvdu(2)=aux(1,3) ! dvdr
   dvdu(3)=aux(2,1) ! dvdtheta
   dvdu(4)=aux(1,4) ! dvdphi
   deallocate(aux)
   return
end subroutine GET_V_AND_DERIVS_WYCKOFFP4MM
end module WYCKOFF_P4MM_MOD

MODULE PES_H2LiF001_MOD
! Massive name spaces
use SYSTEM_MOD
use LINK_FUNCTION1D_MOD
! Selective name spaces
use UNITS_MOD, only: Length, pi
use PES_MOD, only: PES
use CRP3D_MOD, only: CRP3D
use EXTRAPOL_TO_VACUUM_MOD, only: Vacuumpot
use FOURIER_P4MM_MOD, only: Fourierp4mm
use FOURIER3D_P4MM_MOD, only: Fourier3d_p4mm
use WYCKOFF_P4MM_MOD, only: WyckoffSitio,Wyckoffp4mm,TermsInfo
use FOURIER3D_P4MM_MOD, only: Fourier3d_p4mm
use AOTUS_MODULE, only: flu_State, OPEN_CONFIG_FILE, CLOSE_CONFIG, AOT_GET_VAL
use AOT_TABLE_MODULE, only: AOT_TABLE_OPEN, AOT_TABLE_CLOSE, AOT_TABLE_LENGTH, AOT_TABLE_GET_VAL
#if DEBUG
use DEBUG_MOD, only: verbose_write, debug_write
#endif
IMPLICIT NONE
!/////////////////////////////////////////////////
! TYPE: CRP6D
!
!-------------------------------------------------
type,extends(PES):: PES_H2LiF001
   integer(kind=4):: nsites=4
   integer(kind=4):: natomic=1
   logical:: is_interpolated=.false.
   logical:: is_homonucl=.false.
   logical:: is_smooth=.false.
   logical:: is_shifted=.false.
   logical:: is_resized=.false.
   real(kind=8):: zvacuum
   integer(kind=4),dimension(2):: grid=[25,50]
   type(Wyckoffp4mm),dimension(:),allocatable:: wyckoffSite
   type(Pes_HLiF001_NS),dimension(:),allocatable:: atomicCrp
   type(VacuumPot):: farpot
   type(Logistic_func):: dumpfunc
   character(len=30):: extrapol2vac_flag
   integer(kind=4),dimension(:,:),allocatable:: kListXY
   character(len=1),dimension(:),allocatable:: parityListXY
   character(len=2),dimension(:),allocatable:: irrepListXY
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE:: kListAngle
   character(len=1),dimension(:),allocatable:: parityListAngle
   character(len=2),dimension(:),allocatable:: irrepListAngle
   integer(kind=4):: totalTerms=0

   CONTAINS
      ! Initialization block
      procedure,public:: initialize => INITIALIZE_CRP6D
      ! Set block
      procedure,public:: SET_SMOOTH => SET_SMOOTH_CRP6D
      ! Get block
      procedure,public:: GET_V_AND_DERIVS => GET_V_AND_DERIVS_CRP6D
      procedure,public:: GET_V_AND_DERIVS_PURE => GET_V_AND_DERIVS_PURE_CRP6D
      procedure,public:: GET_V_AND_DERIVS_SMOOTH => GET_V_AND_DERIVS_SMOOTH_CRP6D
      ! Tools block
      procedure,public:: SMOOTH => SMOOTH_CRP6D
      procedure,public:: INTERPOL => INTERPOL_CRP6D
      procedure,public:: EXTRACT_VACUUMSURF => EXTRACT_VACUUMSURF_CRP6D
      procedure,public:: ADD_VACUUMSURF => ADD_VACUUMSURF_CRP6D
      procedure,public:: INTERPOL_NEW_RZGRID => INTERPOL_NEW_RZGRID_CRP6D
      ! Plot toolk
      procedure,public:: PLOT1D_THETA => PLOT1D_THETA_CRP6D
      procedure,public:: PLOT1D_ATOMIC_INTERAC_THETA => PLOT1D_ATOMIC_INTERAC_THETA_CRP6D
      procedure,public:: PLOT1D_PHI => PLOT1D_PHI_CRP6D
      procedure,public:: PLOT1D_ATOMIC_INTERAC_PHI => PLOT1D_ATOMIC_INTERAC_PHI_CRP6D
      procedure,public:: PLOT1D_R => PLOT1D_R_CRP6D
      procedure,public:: PLOT1D_Z => PLOT1D_Z_CRP6D
      procedure,public:: PLOT_XYMAP => PLOT_XYMAP_CRP6D
      procedure,public:: PLOT_RZMAP => PLOT_RZMAP_CRP6D
      procedure,public:: PLOT_ATOMIC_INTERAC_RZ => PLOT_ATOMIC_INTERAC_RZ_CRP6D
      ! Enquire block
      procedure,public:: is_allowed => is_allowed_CRP6D
end type PES_H2LiF001
CONTAINS
!###########################################################
!# SUBROUTINE: INITIALIZE_CRP6D
!###########################################################
!> @brief
!! Specific implementation of initialize PES from file
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INITIALIZE_CRP6D(this,filename,tablename)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(OUT)::this
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: filename,tablename
   ! Local variables
   CHARACTER(LEN=:),ALLOCATABLE:: auxstring
   ! Run section
   SELECT CASE(allocated(system_inputfile) .or. .not.present(filename))
      CASE(.TRUE.)
         auxstring=system_inputfile
      CASE(.FALSE.)
         auxstring=filename
   END SELECT
   SELECT CASE(present(tablename))
      CASE(.TRUE.)
         CALL this%READ(filename=auxstring,tablename=tablename)
      CASE(.FALSE.)
         CALL this%READ(filename=auxstring,tablename='pes')
   END SELECT
   CALL this%INTERPOL()
   SELECT CASE(this%is_resized)
      CASE(.TRUE.)
         CALL this%INTERPOL_NEW_RZGRID(this%grid(1),this%grid(2))
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE INITIALIZE_CRP6D
!###########################################################
!# SUBROUTINE: GET_ATOMICPOT_AND_DERIVS_CRP6D
!###########################################################
!> @brief
!! Get atomic potentials and derivatives
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_ATOMICPOT_AND_DERIVS_CRP6D(this,molecx,atomicx,v,dvdu)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(IN):: this
   REAL(KIND=8),DIMENSION(6),INTENT(IN):: molecx
   REAL(KIND=8),DIMENSION(2),INTENT(OUT):: v
   REAL(KIND=8),DIMENSION(6),INTENT(OUT):: dvdu
   REAL(KIND=8),DIMENSION(6),INTENT(OUT):: atomicx
   ! Local variables
   REAL(KIND=8):: vcorra,vcorrb
   ! Run section
   atomicx(:)=from_molecular_to_atomic(molecx)
   SELECT CASE(this%natomic)
      CASE(1)
         CALL this%atomiccrp(1)%GET_V_AND_DERIVS(atomicx(1:3),v(1),dvdu(1:3))
         CALL this%atomiccrp(1)%GET_V_AND_DERIVS(atomicx(4:6),v(2),dvdu(4:6))
      CASE(2)
         CALL this%atomiccrp(1)%GET_V_AND_DERIVS(atomicx(1:3),v(1),dvdu(1:3))
         CALL this%atomiccrp(2)%GET_V_AND_DERIVS(atomicx(4:6),v(2),dvdu(4:6))
      CASE DEFAULT
         WRITE(0,*) "GET_ATOMICPOT_AND_DERIVS_CRP6D ERR: wrong number of atomic potentials"
         CALL EXIT(1)
   END SELECT
   vcorra=this%dumpfunc%getvalue(atomicx(3))
   vcorrb=this%dumpfunc%getvalue(atomicx(6))
   dvdu(1)=dvdu(1)*vcorra
   dvdu(2)=dvdu(2)*vcorra
   dvdu(3)=dvdu(3)*vcorra+v(1)*this%dumpfunc%getderiv(atomicx(3))
   dvdu(4)=dvdu(4)*vcorrb
   dvdu(5)=dvdu(5)*vcorrb
   dvdu(6)=dvdu(6)*vcorrb+v(2)*this%dumpfunc%getderiv(atomicx(6))
   v(1)=v(1)*vcorra
   v(2)=v(2)*vcorrb
   RETURN
END SUBROUTINE GET_ATOMICPOT_AND_DERIVS_CRP6D
!###########################################################
!# SUBROUTINE: SET_SMOOTH_CRP6D
!###########################################################
!> @brief
!! Common set function. Sets is_smooth atribute
!-----------------------------------------------------------
SUBROUTINE SET_SMOOTH_CRP6D(this,is_smooth)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(INOUT):: this
   LOGICAL,INTENT(IN) :: is_smooth
   ! Run section
   this%is_smooth=is_smooth
   RETURN
END SUBROUTINE SET_SMOOTH_CRP6D
!###########################################################
!# SUBROUTINE: CHEAT_CARTWHEEL_ONTOP_CRP6D
!###########################################################
!> @brief
!! This subroutine changes data from a raw cut2d input by the
!! value interpolated from atomic potentials. Only data too close
!! to the plane @f$ \pi: z=\over{r}{2}+rumpling@f$
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE CHEAT_CARTWHEEL_ONTOP_CRP6D(this,wyckoff,cut2d,toptype,dmax)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O
   CLASS(CRP6D),INTENT(INOUT) :: this
   INTEGER(KIND=4),INTENT(IN) :: wyckoff,cut2d,toptype
   REAL(KIND=8),INTENT(IN) :: dmax
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   REAL(KIND=8) :: r,z,rump,interac
   REAL(KIND=8) :: x,y ! x,y position
   REAL(KIND=8) :: m1, m2 ! masses
   REAL(KIND=8) :: dist ! distance to the plane
   INTEGER(KIND=4) :: nr,nz
   REAL(KIND=8),DIMENSION(3) :: pos1,pos2
   ! Run section
   nr=this%wyckoffsite(wyckoff)%zrcut(cut2d)%getgridsizeR()
   nz=this%wyckoffsite(wyckoff)%zrcut(cut2d)%getgridsizeZ()
   x=this%wyckoffsite(wyckoff)%x
   y=this%wyckoffsite(wyckoff)%y
   m1=system_mass(1)
   m2=system_mass(2)
   rump=this%atomiccrp(1)%getrumpling(toptype)
   DO i = 1, nr
      DO j = 1, nz
        r=this%wyckoffsite(wyckoff)%zrcut(cut2d)%getgridvalueR(i)
        z=this%wyckoffsite(wyckoff)%zrcut(cut2d)%getgridvalueZ(j)
        dist=(2.D0/sqrt(5.D0))*(z+rump-r/2.D0)
        pos1=(/x,y,z+(m2/(m1+m2))*r/)
        pos2=(/x,y,z-(m1/(m1+m2))*r/)
        SELECT CASE(dist <= dmax)
           CASE(.TRUE.)
              SELECT CASE(this%is_homonucl)
                 CASE(.TRUE.)
                    interac=this%atomiccrp(1)%getpot(pos1)+this%atomiccrp(1)%getpot(pos2)!+this%farpot%getpot(r)
                 CASE(.FALSE.)
                    interac=this%atomiccrp(1)%getpot(pos1)+this%atomiccrp(2)%getpot(pos2)!+this%farpot%getpot(r)
              END SELECT
              CALL this%wyckoffsite(wyckoff)%zrcut(cut2d)%CHANGEPOT_AT_GRIDPOINT(i,j,interac)
           CASE(.FALSE.)
              ! do nothing
        END SELECT
      END DO
   END DO
   RETURN
END SUBROUTINE CHEAT_CARTWHEEL_ONTOP_CRP6D
!###########################################################
!# SUBROUTINE: INTERPOL_NEW_RZGRID_CRP6D
!###########################################################
!> @brief
!! Creates a different R-Z grid using the ones stored in files
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INTERPOL_NEW_RZGRID_CRP6D(this,nRpoints,nZpoints)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(INOUT) :: this
   INTEGER(KIND=4),INTENT(IN) :: nrpoints,nzpoints
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   ! Run section
   DO i = 1, this%nsites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         CALL this%wyckoffsite(i)%zrcut(j)%interrz%INTERPOL_NEWGRID(nrpoints,nzpoints)
      END DO
   END DO
   RETURN
END SUBROUTINE INTERPOL_NEW_RZGRID_CRP6D
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_PURE_CRP6D
!###########################################################
!> @brief
!! Gets the potential and derivatives respect to all DOFs if we
!! are in the pure CRP6D region. Otherwise, there'd be errors and
!! unexpected behaviors.
!
!> @warning
!! - Assumed a Fourier interpolation in XY (symmetry adapted)
!! - Assumed a Fourier interpolation in theta (symmetry adapted)
!! - Assumed a Fourier interpolation in phi (symmetry adapted)
!! - Assumed a bicubic splines interpolation in R-Z (symmetry independent)
!! - Assumed that Z-R grids have the same size
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_PURE_CRP6D(this,x,v,dvdu)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),TARGET,INTENT(IN):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN):: x
   REAL(KIND=8),INTENT(OUT):: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT):: dvdu
   ! Local variables
   INTEGER(KIND=4):: i,j,h ! counters
   REAL(KIND=8):: ma,mb
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: f ! smooth function and derivs
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: xyList
   real(kind=8),dimension(:),allocatable:: phiList
   real(kind=8),dimension(:,:),allocatable:: completeGeomList
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: aux1
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: aux2
   REAL(KIND=8),DIMENSION(6):: atomicx
   REAL(KIND=8),DIMENSION(2):: atomic_v
   REAL(KIND=8),DIMENSION(6):: atomic_dvdu
   REAL(KIND=8),DIMENSION(3):: dvdu_atomicA,dvdu_atomicB
   TYPE(Fourier3d_p4mm):: fouInterpol
   CHARACTER(LEN=*),PARAMETER:: routinename="GET_V_AND_DERIVS_CRP6D: "
   ! Run section
   SELECT CASE(this%is_smooth)
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         WRITE(0,*) "GET_V_AND_DERIVS_CRP6D ERR: smooth the PES first (CALL thispes%SMOOTH())"
         CALl EXIT(1)
   END SELECT
   ! f(1,:) smooth potential values
   ! f(2,:) smooth dvdz
   ! f(3,:) smooth dvdr
   ! f(4,:) smooth dvdtheta
   ! f(5,:) smooth dvdphi
   allocate( f(5,this%totalTerms)      )
   allocate( xyList(this%totalTerms,2) )
   allocate( phiList(this%totalTerms)  )
   allocate( completeGeomList(this%totalTerms,3)  )
   h=0
   do i = 1, this%nsites
      do j=1,size(this%wyckoffsite(i)%phiList)
         h=h+1
         xyList(h,1)=this%wyckoffSite(i)%x
         xyList(h,2)=this%wyckoffSite(i)%y
         phiList(h)=this%wyckoffSite(i)%phiList(j)
         call this%wyckoffsite(i)%get_v_and_derivs( [x(3),x(4),x(5),phiList(h)],f(1,h),f(2:5,h) )
      enddo
   end do
   completeGeomList(:,1:2)=xyList(:,:)
   completeGeomList(:,3)=phiList(:)
   call fouInterpol%read( x=completeGeomList(:,:),f=f(1,:) )
   call fouInterpol%add_moreFuncs( f=f(2:4,1:this%totalTerms) )
   call fouInterpol%initializeTerms()
   call fouInterpol%setKlist( kListXY=this%kListXY,kListAngle=this%kListAngle )
   call fouInterpol%setParityList( parityListXY=this%parityListXY,parityListAngle=this%parityListAngle )
   call fouInterpol%setIrrepList( irrepListXY=this%irrepListXY,irrepListAngle=this%irrepListAngle )
   CALL fouInterpol%interpol()
   allocate(aux1(4))
   allocate(aux2(4,3))
   call fouInterpol%get_allFuncs_and_derivs( x=[x(1),x(2),x(6)],f=aux1,dfdx=aux2 )
#ifdef DEBUG
   !-------------------------------------
   ! Results for the smooth potential
   !-------------------------------------
   ! v=aux1(1)         ! v
   ! dvdu(1)=aux2(1,1) ! dvdx
   ! dvdu(2)=aux2(1,2) ! dvdy
   ! dvdu(3)=aux1(2)   ! dvdz
   ! dvdu(4)=aux1(3)   ! dvdr
   ! dvdu(5)=aux1(4)   ! dvdtheta
   ! dvdu(6)=aux2(1,3) ! dvdphi
   CALL DEBUG_WRITE(routinename,"Smooth V at each thetaCut:")
   CALL DEBUG_WRITE(routinename,f(1,:))
   CALL DEBUG_WRITE(routinename,"Smooth dv/dz at each thetaCut:")
   CALL DEBUG_WRITE(routinename,f(2,:))
   CALL DEBUG_WRITE(routinename,"Smooth dv/dr at each thetaCut:")
   CALL DEBUG_WRITE(routinename,f(3,:))
   CALL DEBUG_WRITE(routinename,"Smooth dv/dtheta at each thetaCut:")
   CALL DEBUG_WRITE(routinename,f(4,:))
   CALL DEBUG_WRITE(routinename,"Smooth interpolated values:")
   CALL DEBUG_WRITE(routinename,"v: ",aux1(1))   ! v
   CALL DEBUG_WRITE(routinename,"dvdx: ",aux2(1,1)) ! dvdx
   CALL DEBUG_WRITE(routinename,"dvdy: ",aux2(1,2)) ! dvdy
   CALL DEBUG_WRITE(routinename,"dvdz: ",aux1(2))   ! dvdz
   CALL DEBUG_WRITE(routinename,"dvdr: ",aux1(3))   ! dvdr
   CALL DEBUG_WRITE(routinename,"dvdtheta: ",aux1(4))   ! dvdtheta
   CALL DEBUG_WRITE(routinename,"dvdphi: ",aux2(1,3))   ! dvdphi
#endif
   !--------------------------------------
   ! Results for the real potential
   !-------------------------------------
   CALL this%GET_ATOMICPOT_AND_DERIVS(x,atomicx,atomic_v,atomic_dvdu)
   dvdu_atomicA=atomic_dvdu(1:3)
   dvdu_atomicB=atomic_dvdu(4:6)
#ifdef DEBUG
   CALL DEBUG_WRITE(routinename,"Contributions of the atomic potential: ")
   CALL DEBUG_WRITE(routinename, "Position Atom A: ",atomicx(1:3))
   CALL DEBUG_WRITE(routinename,"Va: ",atomic_v(1))
   CALL DEBUG_WRITE(routinename,"dVa/dxa; dVa/dya: ",dvdu_atomicA)
   CALL DEBUG_WRITE(routinename, "Position Atom B: ",atomicx(4:6))
   CALL DEBUG_WRITE(routinename,"Vb: ",atomic_v(2))
   CALL DEBUG_WRITE(routinename,"dVa/dxa; dVb/dyb: ",dvdu_atomicB)
#endif
   v=aux1(1)+sum(atomic_v)

   ma=system_mass(1)
   mb=system_mass(2)
   dvdu(1)=aux2(1,1)+dvdu_atomicA(1)+dvdu_atomicB(1)
   dvdu(2)=aux2(1,2)+dvdu_atomicA(2)+dvdu_atomicB(2)
   dvdu(3)=aux1(2)+dvdu_atomicA(3)+dvdu_atomicB(3)
   dvdu(4)=aux1(3)&
      +(mb/(ma+mb))*dvdu_atomicA(1)*dcos(x(6))*dsin(x(5))&
      +(mb/(ma+mb))*dvdu_atomicA(2)*dsin(x(6))*dsin(x(5))&
      +(mb/(ma+mb))*dvdu_atomicA(3)*dcos(x(5))&
      -(ma/(ma+mb))*dvdu_atomicB(1)*dcos(x(6))*dsin(x(5))&
      -(ma/(ma+mb))*dvdu_atomicB(2)*dsin(x(6))*dsin(x(5))&
      -(ma/(ma+mb))*dvdu_atomicB(3)*dcos(x(5))
   dvdu(5)=aux1(4)&
      +(mb/(ma+mb))*x(4)*dvdu_atomicA(1)*dcos(x(6))*dcos(x(5))&
      +(mb/(ma+mb))*x(4)*dvdu_atomicA(2)*dsin(x(6))*dcos(x(5))&
      -(mb/(ma+mb))*x(4)*dvdu_atomicA(3)*dsin(x(5))&
      -(ma/(ma+mb))*x(4)*dvdu_atomicB(1)*dcos(x(6))*dcos(x(5))&
      -(ma/(ma+mb))*x(4)*dvdu_atomicB(2)*dsin(x(6))*dcos(x(5))&
      +(ma/(ma+mb))*x(4)*dvdu_atomicB(3)*dsin(x(5))
   dvdu(6)=aux2(1,3)&
      -(mb/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicA(1)*dsin(x(6))&
      +(mb/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicA(2)*dcos(x(6))&
      +(ma/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicB(1)*dsin(x(6))&
      -(ma/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicB(2)*dcos(x(6))
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_PURE_CRP6D
!###########################################################
!# SUBROUTINE: GET_V AND DERIVS_CRP6D
!###########################################################
!> @brief
!! Gets potential for any configuration. Discriminates between 3 regions:
!! - Pure CRP6D region: zmin to zcrp
!! - Extrapolation region zcrp to zvacuum
!! - Vacuum region: further than zvacuum
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_CRP6D(this,X,v,dvdu,errCode)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),TARGET,INTENT(IN):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   REAL(KIND=8),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: dvdu
   integer(kind=1),optional,intent(out):: errCode
   ! Local variables
   REAL(KIND=8) :: zcrp, zvac ! last CRP6D z value and Z infinity
   REAL(KIND=8) :: vzcrp, vzvac ! potentials at zcrp and zvac
   REAL(KIND=8),DIMENSION(6) :: dvducrp ! derivatives at zcrp
   REAL(KIND=8),DIMENSION(6) :: dvduvac ! derivatives at vacuum
   REAL(KIND=8) :: alpha,beta,gama ! parameters
   CLASS(Function1d),ALLOCATABLE:: extrapolfunc
   INTEGER(KIND=4) :: i !counter
   ! Local Parameter
   REAL(KIND=8),PARAMETER :: zero=0.D0 ! what we will consider zero (a.u.)
   REAL(KIND=8),PARAMETER :: dz=0.5D0 ! 0.25 Angstroems approx
   character(len=*),parameter:: routinename='GET_V_AND_DERIVS_CRP6D: '
   ! Run section
   zcrp=this%wyckoffsite(1)%zrcut(1)%getlastZ()
   zvac=this%zvacuum
   ! Check if we are in the pure CRP6D region
   SELECT CASE(x(3)<= zcrp) !easy
      CASE(.TRUE.)
         call this%get_v_and_derivs_pure(x,v,dvdu)
         RETURN
      CASE(.FALSE.)
         ! do nothing, next switch
   END SELECT
   ! Check if we are in the extrapolation region
   SELECT CASE(x(3)>zcrp .AND. x(3)<zvac)
      CASE(.TRUE.) ! uff
         ! Set potential and derivs
         vzvac=this%farpot%getpot(x(4))
         CALL this%GET_V_AND_DERIVS_PURE([x(1),x(2),zcrp,x(4),x(5),x(6)],vzcrp,dvducrp)
         dvduvac(1:3)=zero
         dvduvac(4)=this%farpot%getderiv(x(4))
         dvduvac(5:6)=zero
         ! Check kind of extrapolation
         SELECT CASE(this%extrapol2vac_flag)
            CASE("Xexponential")
               ALLOCATE(Xexponential_func::extrapolfunc)
               ! Extrapol potential
               beta=-1.D0/zvac
               alpha=(vzcrp-vzvac)/(zcrp*dexp(beta*zcrp)-zvac*dexp(beta*zvac))
               gama=vzvac-alpha*zvac*dexp(beta*zvac)
               CALL extrapolfunc%READ([alpha,beta])
               v=extrapolfunc%getvalue(x(3))+gama
               dvdu(3)=extrapolfunc%getderiv(x(3))
               ! Extrapol derivatives
               DO i = 1, 6
                  SELECT CASE(i)
                     CASE(3)
                        ! Skip dvdz
                     CASE DEFAULT
                        beta=-1.D0/zvac
                        alpha=(dvducrp(i)-dvduvac(i))/(zcrp*dexp(beta*zcrp)-zvac*dexp(beta*zvac))
                        gama=dvduvac(i)-alpha*zvac*dexp(beta*zvac)
                        CALL extrapolfunc%READ([alpha,beta])
                        dvdu(i)=extrapolfunc%getvalue(x(3))+gama
                  END SELECT
               END DO
               DEALLOCATE(extrapolfunc)
               RETURN

            CASE("Linear")
               ALLOCATE(Linear_func::extrapolfunc)
               ! Extrapol potential
               beta=(vzcrp*zvac-vzvac*zcrp)/(zvac-zcrp)
               alpha=(vzvac-beta)/zvac
               CALL extrapolfunc%READ([alpha,beta])
               v=extrapolfunc%getvalue(x(3))
               dvdu(3)=extrapolfunc%getderiv(x(3))
               ! Extrapol derivatives
               DO i = 1, 6
                  SELECT CASE(i)
                     CASE(3)
                        ! do nothing
                     CASE DEFAULT
                        beta=(dvducrp(i)*zvac-dvduvac(i)*zcrp)/(zvac-zcrp)
                        alpha=(dvduvac(i)-beta)/zvac
                        CALL extrapolfunc%READ([alpha,beta])
                        dvdu(i)=extrapolfunc%getvalue(x(3))
                  END SELECT
               END DO
               DEALLOCATE(extrapolfunc)
               RETURN

            CASE DEFAULT
               WRITE(0,*) "GET_V_AND_DERIVS_CRP6D ERR: type of extrapolation function isn't implemented yet"
               WRITE(0,*) "Implemented ones: Linear, Xexponential, None"
               CALL EXIT(1)
         END SELECT
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   ! Check if we are in the Vacuum region
   SELECT CASE(x(3)>=zvac) !easy
      CASE(.TRUE.)
         v=this%farpot%getpot(x(4))
         dvdu(1:3)=0.D0
         dvdu(4)=this%farpot%getderiv(x(4))
         dvdu(5:6)=0.D0
      CASE(.FALSE.) ! this's the last switch!
#ifdef DEBUG
         call debug_write(routinename,'Unclassificable point')
         call debug_write(routinename,'Asking potential at Z: ',x(3))
#endif
          errCode=1_1
   END SELECT
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_CRP6D
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_SMOOTH_CRP6D
!###########################################################
!> @brief
!! Gets the SMOOTH potential and the derivatives respect to all DOFs
!
!> @details
!! This routine doesn't check PES smoothness, so it can be used to get
!! interpolated values for the raw PES
!> @warning
!! - Assumed a 2D Fourier interpolation in XY (symmetry adapted)
!! - Assumed a 1D Fourier interpolation in theta (symmetry adapted)
!! - Assumed a 1D Fourier interpolation in phi (symmetry adapted)
!! - Assumed a bicubic splines interpolation in R-Z (symmetry independent)
!! - Assumed that Z-R grids have the same size. This can be ensured using
!!   using subroutine INTERPOL_NEWGRID
!
!> @see
!! interpol_newgrid_crp6d
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_SMOOTH_CRP6D(this,x,v,dvdu)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),TARGET,INTENT(IN):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN):: x
   REAL(KIND=8),INTENT(OUT):: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT):: dvdu
   ! Local variables
   INTEGER(KIND=4):: i,j,h ! counters
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: f ! smooth function and derivs
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: xyList
   real(kind=8),dimension(:),allocatable:: phiList
   real(kind=8),dimension(:,:),allocatable:: completeGeomList
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: aux1
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: aux2
   TYPE(Fourier3d_p4mm):: fouInterpol
   CHARACTER(LEN=*),PARAMETER:: routinename="GET_V_AND_DERIVS_CRP6D: "
   ! Run section
   ! f(1,:) smooth potential values
   ! f(2,:) smooth dvdz
   ! f(3,:) smooth dvdr
   ! f(4,:) smooth dvdtheta
   ! f(5,:) smooth dvdphi
   allocate( f(5,this%totalTerms)      )
   allocate( xyList(this%totalTerms,2) )
   allocate( phiList(this%totalTerms)  )
   allocate( completeGeomList(this%totalTerms,3)  )
   h=0
   do i = 1, this%nsites
      do j=1,size(this%wyckoffsite(i)%phiList)
         h=h+1
         xyList(h,1)=this%wyckoffSite(i)%x
         xyList(h,2)=this%wyckoffSite(i)%y
         phiList(h)=this%wyckoffSite(i)%phiList(j)
         call this%wyckoffsite(i)%get_v_and_derivs( [x(3),x(4),x(5),phiList(h)],f(1,h),f(2:5,h) )
      enddo
   end do
   completeGeomList(:,1:2)=xyList(:,:)
   completeGeomList(:,3)=phiList(:)
   call fouInterpol%read( x=completeGeomList(:,:),f=f(1,:) )
   call fouInterpol%add_moreFuncs( f=f(2:4,1:this%totalTerms) )
   call fouInterpol%initializeTerms()
   call fouInterpol%setKlist( kListXY=this%kListXY,kListAngle=this%kListAngle )
   call fouInterpol%setParityList( parityListXY=this%parityListXY,parityListAngle=this%parityListAngle )
   call fouInterpol%setIrrepList( irrepListXY=this%irrepListXY,irrepListAngle=this%irrepListAngle )
   CALL fouInterpol%interpol()
   allocate(aux1(4))
   allocate(aux2(4,3))
   call fouInterpol%get_allFuncs_and_derivs( x=[x(1),x(2),x(6)],f=aux1,dfdx=aux2 )
#ifdef DEBUG
   !-------------------------------------
   ! Results for the smooth potential
   !-------------------------------------
   ! v=aux1(1)         ! v
   ! dvdu(1)=aux2(1,1) ! dvdx
   ! dvdu(2)=aux2(1,2) ! dvdy
   ! dvdu(3)=aux1(2)   ! dvdz
   ! dvdu(4)=aux1(3)   ! dvdr
   ! dvdu(5)=aux1(4)   ! dvdtheta
   ! dvdu(6)=aux2(1,3) ! dvdphi
   CALL DEBUG_WRITE(routinename,"Smooth V at each thetaCut:")
   CALL DEBUG_WRITE(routinename,f(1,:))
   CALL DEBUG_WRITE(routinename,"Smooth dv/dz at each thetaCut:")
   CALL DEBUG_WRITE(routinename,f(2,:))
   CALL DEBUG_WRITE(routinename,"Smooth dv/dr at each thetaCut:")
   CALL DEBUG_WRITE(routinename,f(3,:))
   CALL DEBUG_WRITE(routinename,"Smooth dv/dtheta at each thetaCut:")
   CALL DEBUG_WRITE(routinename,f(4,:))
   CALL DEBUG_WRITE(routinename,"Smooth interpolated values:")
   CALL DEBUG_WRITE(routinename,"v: ",aux1(1))   ! v
   CALL DEBUG_WRITE(routinename,"dvdx: ",aux2(1,1)) ! dvdx
   CALL DEBUG_WRITE(routinename,"dvdy: ",aux2(1,2)) ! dvdy
   CALL DEBUG_WRITE(routinename,"dvdz: ",aux1(2))   ! dvdz
   CALL DEBUG_WRITE(routinename,"dvdr: ",aux1(3))   ! dvdr
   CALL DEBUG_WRITE(routinename,"dvdtheta: ",aux1(4))   ! dvdtheta
   CALL DEBUG_WRITE(routinename,"dvdphi: ",aux2(1,3))   ! dvdphi
#endif
   v=aux1(1)
   dvdu(1)=aux2(1,1)
   dvdu(2)=aux2(1,2)
   dvdu(3)=aux1(2)
   dvdu(4)=aux1(3)
   dvdu(5)=aux1(4)
   dvdu(6)=aux2(1,3)
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_SMOOTH_CRP6D
!###########################################################
!# SUBROUTINE: INTERPOL_CRP6D
!###########################################################
!> @brief
!! Interpolates CRP6D smooth 2dcuts in R and Z.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INTERPOL_CRP6D(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(INOUT)::this
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   ! Run section
   CALL this%EXTRACT_VACUUMSURF()
   CALL this%SMOOTH()
   !CALL this%SMOOTH_EXTRA()
   DO i = 1, this%nsites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         CALL this%wyckoffsite(i)%zrcut(j)%INTERPOL()
      END DO
   END DO
   this%is_interpolated=.TRUE.
   RETURN
END SUBROUTINE INTERPOL_CRP6D
!###########################################################
!# SUBROUTINE: RAWINTERPOL_CRP6D
!###########################################################
!> @brief
!! Interpolates CRP6D 2dcuts in R and Z. There is not smoothing
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE RAWINTERPOL_CRP6D(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(INOUT)::this
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   ! Run section
   CALL this%EXTRACT_VACUUMSURF()
   DO i = 1, this%nsites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         CALL this%wyckoffsite(i)%zrcut(j)%INTERPOL()
      END DO
   END DO
   this%is_interpolated=.TRUE.
   RETURN
END SUBROUTINE RAWINTERPOL_CRP6D
!###########################################################
!# SUBROUTINE: SMOOTH_CRP6D
!###########################################################
!> @brief
!! Smooths a RZ-2dcut of the potential using atomic potentials
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 25/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SMOOTH_CRP6D(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
    CLASS(CRP6D),INTENT(INOUT):: this
   ! Local variables
   REAL(KIND=8),DIMENSION(6):: molcoord,atomcoord,dummy
   INTEGER(KIND=4):: nr,nz
   INTEGER(KIND=4):: i,j,k,l ! counters
   CHARACTER(LEN=*),PARAMETER:: routinename="SMOOTH_CRP6D: "
   REAL(KIND=8):: newpot
   REAL(KIND=8),DIMENSION(2):: atomic_v
   ! Run section
   DO i = 1, this%nsites ! cycle wyckoff sites
      molcoord(1)=this%wyckoffsite(i)%x
      molcoord(2)=this%wyckoffsite(i)%y
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         molcoord(5)=this%wyckoffsite(i)%zrcut(j)%theta
         molcoord(6)=this%wyckoffsite(i)%zrcut(j)%phi
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizeR()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizeZ()
         DO k = 1, nr
            DO l = 1, nz
               molcoord(3)=this%wyckoffsite(i)%zrcut(j)%getgridvalueZ(l)
               molcoord(4)=this%wyckoffsite(i)%zrcut(j)%getgridvalueR(k)
               CALL this%GET_ATOMICPOT_AND_DERIVS(molcoord,atomcoord,atomic_v,dummy)
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)-sum(atomic_v)
               CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
            END DO
         END DO
      END DO
   END DO
   this%is_smooth=.TRUE.
   RETURN
END SUBROUTINE SMOOTH_CRP6D
!###########################################################
!# SUBROUTINE: ROUGH_CRP6D
!###########################################################
!> @brief
!! ROUGHs a RZ-2dcut of the potential using atomic potentials
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE ROUGH_CRP6D(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
    CLASS(CRP6D),INTENT(INOUT):: this
   ! Local variables
   REAL(KIND=8),DIMENSION(6):: molcoord,atomcoord
   INTEGER(KIND=4):: nr,nz
   INTEGER(KIND=4):: i,j,k,l ! counters
   REAL(KIND=8):: newpot
   ! Parameters
   CHARACTER(LEN=*),PARAMETER:: routinename="ROUGH_CRP6D: "
   ! Run section
   DO i = 1, this%nsites ! cycle wyckoff sites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         molcoord(1)=this%wyckoffsite(i)%zrcut(j)%x
         molcoord(2)=this%wyckoffsite(i)%zrcut(j)%y
         molcoord(5)=this%wyckoffsite(i)%zrcut(j)%theta
         molcoord(6)=this%wyckoffsite(i)%zrcut(j)%phi
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizer()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizez()
         DO k = 1, nr
            DO l = 1, nz
               molcoord(3)=this%wyckoffsite(i)%zrcut(j)%getgridvalueZ(l)
               molcoord(4)=this%wyckoffsite(i)%zrcut(j)%getgridvalueR(k)
               atomcoord(:)=from_molecular_to_atomic(molcoord)
#ifdef DEBUG
               CALL DEBUG_WRITE(routinename,"Molecular coords:")
               CALL DEBUG_WRITE(routinename,molcoord)
               CALL DEBUG_WRITE(routinename,"Atomic coords:")
               CALL DEBUG_WRITE(routinename,atomcoord)
#endif
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)
               SELECT CASE(this%natomic)
                  CASE(1)
                     newpot=newpot+this%atomiccrp(1)%getpot(atomcoord(1:3))+&
                        this%atomiccrp(1)%getpot(atomcoord(4:6))
                     CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  CASE(2)
                     newpot=newpot+this%atomiccrp(1)%getpot(atomcoord(1:3))+&
                        this%atomiccrp(2)%getpot(atomcoord(4:6))
                     CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  CASE DEFAULT
                     WRITE(0,*) "SMOOTH_CRP6D ERR: Something is wrong with number of atomic crp potentials"
                     CALL EXIT(1)
               END SELECT
            END DO
         END DO
      END DO
   END DO
   this%is_smooth=.FALSE.
   RETURN
END SUBROUTINE ROUGH_CRP6D
!###########################################################
!# SUBROUTINE: SMOOTH_EXTRA_CRP6D
!###########################################################
!> @brief
!! Smooths a RZ-2dcut of the potential using atomic potentials as well
!! as the vacuum potential
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SMOOTH_EXTRA_CRP6D(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
    CLASS(CRP6D),INTENT(INOUT) :: this
   ! Local variables
   REAL(KIND=8),DIMENSION(6) :: molcoord,atomcoord
   INTEGER(KIND=4) :: nr,nz
   INTEGER(KIND=4) :: i,j,k,l ! counters
   CHARACTER(LEN=20),PARAMETER :: routinename="SMOOTH_EXTRA_CRP6D: "
   REAL(KIND=8) :: newpot
   ! Run section
   DO i = 1, this%nsites ! cycle wyckoff sites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         molcoord(1)=this%wyckoffsite(i)%zrcut(j)%x
         molcoord(2)=this%wyckoffsite(i)%zrcut(j)%y
         molcoord(5)=this%wyckoffsite(i)%zrcut(j)%theta
         molcoord(6)=this%wyckoffsite(i)%zrcut(j)%phi
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizer()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizez()
         DO k = 1, nr
            DO l = 1, nz
               molcoord(3)=this%wyckoffsite(i)%zrcut(j)%getgridvalueZ(l)
               molcoord(4)=this%wyckoffsite(i)%zrcut(j)%getgridvalueR(k)
               atomcoord(:)=from_molecular_to_atomic(molcoord)
#ifdef DEBUG
               CALL DEBUG_WRITE(routinename,"Molecular coords:")
               CALL DEBUG_WRITE(routinename,molcoord)
               CALL DEBUG_WRITE(routinename,"Atomic coords:")
               CALL DEBUG_WRITE(routinename,atomcoord)
#endif
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)
               SELECT CASE(this%natomic)
                  CASE(1)
                     newpot=newpot-this%atomiccrp(1)%getpot(atomcoord(1:3))&
                        -this%atomiccrp(1)%getpot(atomcoord(4:6))-this%farpot%getpot(molcoord(4))&
                        -0.8*dexp(-molcoord(4))
                     CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  CASE(2)
                     newpot=newpot-this%atomiccrp(1)%getpot(atomcoord(1:3))-&
                        this%atomiccrp(2)%getpot(atomcoord(4:6))-this%farpot%getpot(molcoord(4))
                     CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  CASE DEFAULT
                     WRITE(0,*) "SMOOTH_EXTRA_CRP6D ERR: Something is wrong with number of atomic crp potentials"
                     CALL EXIT(1)
               END SELECT
            END DO
         END DO
      END DO
   END DO
   this%is_smooth=.TRUE.
   RETURN
END SUBROUTINE SMOOTH_EXTRA_CRP6D
!###########################################################
!# SUBROUTINE: EXTRACT_VACUUMSURF_CRP6D
!###########################################################
!> @brief
!! Extracts energy at equilibrium distance in the vacuum and surface energy
!! to an Rz-2dcut of the potential. Shifts the vacuum potential as well.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 25/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE EXTRACT_VACUUMSURF_CRP6D(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
    CLASS(CRP6D),INTENT(INOUT) :: this
   ! Local variables
   INTEGER(KIND=4) :: nr,nz
   INTEGER(KIND=4) :: i,j,k,l ! counters
   REAL(KIND=8) :: newpot
   CHARACTER(LEN=26),PARAMETER :: routinename="EXTRACT_VACUUMSURF_CRP6D: "
   ! Run section
   DO i = 1, this%nsites ! cycle wyckoff sites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizeR()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizeZ()
         DO k = 1, nr
            DO l = 1, nz
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)
               newpot=newpot-this%farpot%getscalefactor()
               CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
            END DO
         END DO
      END DO
   END DO
   CALL this%farpot%SHIFTPOT()
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Potential shifted: ",this%farpot%getscalefactor())
#endif
   this%is_shifted=.TRUE.
   RETURN
END SUBROUTINE EXTRACT_VACUUMSURF_CRP6D
!###########################################################
!# SUBROUTINE: ADD_VACUUMSURF_CRP6D
!###########################################################
!> @brief
!! Adds energy at equilibrium distance in the vacuum and surface energy
!! to an Rz-2dcut of the potential. Shifts the vacuum potential as well.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE ADD_VACUUMSURF_CRP6D(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
    CLASS(CRP6D),INTENT(INOUT) :: this
   ! Local variables
   INTEGER(KIND=4) :: nr,nz
   INTEGER(KIND=4) :: i,j,k,l ! counters
   REAL(KIND=8) :: newpot
   CHARACTER(LEN=22),PARAMETER :: routinename="ADD_VACUUMSURF_CRP6D: "
   ! Run section
   DO i = 1, this%nsites ! cycle wyckoff sites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizeR()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizeZ()
         DO k = 1, nr
            DO l = 1, nz
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)
               newpot=newpot+this%farpot%getscalefactor()
               CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
            END DO
         END DO
      END DO
   END DO
   CALL this%farpot%SHIFTPOT_UNDO()
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Potential re-shifted: ",this%farpot%getscalefactor())
#endif
   this%is_shifted=.FALSE.
   RETURN
END SUBROUTINE ADD_VACUUMSURF_CRP6D
!###########################################################
!# SUBROUTINE: READ_CRP6D
!###########################################################
!> @brief
!! Sets up a CRP6D Object
!-----------------------------------------------------------
SUBROUTINE READ_CRP6D(this,filename,tablename)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(OUT):: this
   CHARACTER(LEN=*),INTENT(IN):: filename,tablename
   ! Local variables
   INTEGER(KIND=4):: i,j,k,n ! counters
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: param_damp
   ! Wyckoff variables
   CHARACTER(LEN=1),DIMENSION(:),ALLOCATABLE:: wyckoff_letters
   CHARACTER(LEN=1024),DIMENSION(:),ALLOCATABLE:: cuts2d_files
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE:: thetablocks_data
   type(TermsInfo),dimension(:),allocatable:: phiTerms
   type(TermsInfo):: thetaTerms
   INTEGER(KIND=4):: n2dcuts
   INTEGER(KIND=4):: nthetablocks
   real(kind=8),dimension(:),allocatable:: phiList

   ! Lua specifications
   TYPE(flu_State):: conf
   INTEGER(KIND=4):: ierr
   INTEGER(KIND=4):: pes_table,crp3d_table,vacfunc_table,dampfunc_table,param_table,extrapol_table
   INTEGER(KIND=4):: resize_table,wyckoff_table,inwyckoff_table,cut2d_table,files_table,kpoints_table
   integer(kind=4):: phi_table,term_table,theta_table,fourier_table,magnitude_table,phiList_table
   INTEGER(KIND=4):: inkpoints_table
   ! Auxiliar, dummy variables
   INTEGER(KIND=4):: auxInt,auxInt2
   REAL(KIND=8):: auxReal
   CHARACTER(LEN=1024):: auxString,auxString2
   TYPE(Length):: len
   type(Angle):: angl
   ! Parameters
   CHARACTER(LEN=*),PARAMETER:: routinename="READ_CRP6D: "
   ! Run section -----------------------
   ! Open Lua config file
   CALL OPEN_CONFIG_FILE(L=conf,ErrCode=ierr,filename=filename)
   ! Open PES table
   CALL AOT_TABLE_OPEN(L=conf,thandle=pes_table,key=tablename)
   ! get pes.kind
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='kind',val=auxstring)
   CALL this%SET_PESTYPE(trim(auxstring))
   SELECT CASE(trim(auxstring))
      CASE('CRP6D')
         ! do nothing
      CASE DEFAULT
         WRITE(0,*) 'READ CRP6D ERR: wrong kind of PES. Expected: CRP6D. Encountered: '//trim(auxstring)
   END SELECT
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'PES type: '//trim(auxstring))
#endif
   ! get pes.name
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='name',val=auxstring)
   CALL this%SET_ALIAS(trim(auxstring))
#ifdef DEBUG
   CALl VERBOSE_WRITE(routinename,'PES name: '//trim(auxstring))
#endif
   ! get pes.dimensions
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='dimensions',val=auxint)
   CALL this%SET_DIMENSIONS(auxint)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'PES dimensions: ',auxint)
#endif
   ! get crp3d subpes
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=crp3d_table,key='crp3dPes')
   this%natomic=aot_table_length(L=conf,thandle=crp3d_table)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Atomic potentials found: ",this%natomic)
#endif
   SELECT CASE(this%natomic)
      CASE(1)
         this%is_homonucl=.TRUE.
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,'This PES is for homonuclear projectiles')
#endif
      CASE(2)
         this%is_homonucl=.FALSE.
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,'This PES is for heteronuclear projectiles')
#endif
      CASE DEFAULT
         WRITE(0,*) "READ_CRP6D ERR: Wrong number of atomic potentials. Allowed number: 1 or 2."
         CALL EXIT(1)
   END SELECT
   ALLOCATE(this%atomiccrp(this%natomic))
   DO i = 1, this%natomic
      CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=crp3d_table,pos=i,val=auxstring)
#ifdef DEBUG
      CALl VERBOSE_WRITE(routinename,'Atomic potential keyword: '//trim(auxstring))
#endif
      CALL this%atomiccrp(i)%INITIALIZE(filename=filename,tablename=trim(auxstring))
   END DO
   CALL AOT_TABLE_CLOSE(L=conf,thandle=crp3d_table)
   ! get pes.vacuumFunction
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=vacfunc_table,key='vacuumFunction')
   CALL AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=vacfunc_table,key='kind',val=auxstring)
   SELECT CASE(trim(auxstring))
      CASE('Numerical')
         CALL AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=vacfunc_table,key='source',val=auxstring)
         CALL this%farpot%INITIALIZE(trim(auxstring))
      CASE DEFAULT
         WRITE(0,*) 'READ_CRP6D ERR: wrong kind of vacuuum function: '//trim(auxstring)
         WRITE(0,*) 'Implemented ones: Numerical'
         WRITE(0,*) 'Case sensitive'
         CALL EXIT(1)
   END SELECT
   CALL AOT_TABLE_CLOSE(L=conf,thandle=vacfunc_table)
   ! get pes.dampFunction
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=dampfunc_table,key='dampFunction')
   CALL AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=dampfunc_table,key='kind',val=auxstring)
   SELECT CASE(trim(auxstring))
      CASE("Logistic")
         ALLOCATE(Logistic_func::this%dumpfunc)
         ALLOCATE(param_damp(2))
         CALL AOT_TABLE_OPEN(L=conf,parent=dampfunc_table,thandle=param_table,key='param')
         auxint=aot_table_length(L=conf,thandle=param_table)
         SELECT CASE(auxint)
            CASE(2)
               ! do nothing
            CASE DEFAULT
               WRITE(0,*) 'READ_CRP6D ERR: incorrect number of parameters in damp function: ',auxint
               CALL EXIT(1)
         END SELECT
         CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=param_table,pos=1,val=param_damp(1))
         CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=param_table,pos=2,val=param_damp(2))
         CALL AOT_TABLE_CLOSE(L=conf,thandle=param_table)
         CALL this%dumpfunc%READ(param_damp)
      CASE("None")
         ALLOCATE(One_func::this%dumpfunc)
      CASE DEFAULT
         WRITE(0,*) "READ_CRP6D ERR: Keyword for dumping function needed"
         WRITE(*,*) "Currently implemented: Logistic, None"
         CALL EXIT(1)
   END SELECT
   CALL AOT_TABLE_CLOSE(L=conf,thandle=dampfunc_table)
   ! get pes.extrapolFunction
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=extrapol_table,key='extrapolFunction')
   CALL AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=extrapol_table,key='kind',val=auxstring)
   SELECT CASE(trim(auxstring))
      CASE("Xexponential","Linear")
         this%extrapol2vac_flag=trim(auxstring)
         CALL AOT_TABLE_OPEN(L=conf,parent=extrapol_table,thandle=param_table,key='upToZ')
         auxint=aot_table_length(L=conf,thandle=param_table)
         CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=param_table,pos=1,val=auxreal)
         CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=param_table,pos=2,val=auxstring)
         CALL len%READ(auxreal,trim(auxstring))
         CALL len%TO_STD()
         this%zvacuum=len%getvalue()
         CALl AOT_TABLE_CLOSE(L=conf,thandle=param_table)
      CASE("None")
         this%zvacuum=0.D0
      CASE DEFAULT
         WRITE(0,*) "READ_CRP6D ERR: Keyword for extrapolation function needed"
         WRITE(0,*) "Currently implemented: None, Xexponential, Linear"
         WRITE(0,*) "Case sensitive."
         CALL EXIT(1)
   END SELECT
   CALL AOT_TABLE_CLOSE(L=conf,thandle=extrapol_table)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Z vacuum: ",this%zvacuum)
#endif
   ! get pes.resize
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=resize_table,key='resize')
   auxint=aot_table_length(L=conf,thandle=resize_table)
   SELECT CASE(auxint)
      CASE(0)
         ! do nothing, resize is not required
      CASE DEFAULT
         CALL AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=resize_table,key='r',val=auxint)
         this%grid(1)=auxint
         CALL AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=resize_table,key='z',val=auxint)
         this%grid(2)=auxint
         if( this%grid(1)/=0 .and. this%grid(2)/=0 ) this%is_resized=.true.
   END SELECT
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,'New grid (R,Z):',this%grid(:))
         call verbose_write(routinename,'Is this PES going to be resized?: ',this%is_resized)
#endif
   CALL AOT_TABLE_CLOSE(L=conf,thandle=resize_table)
   ! get wyckoff sites
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=wyckoff_table,key='wyckoffSite')
   this%nsites=aot_table_length(L=conf,thandle=wyckoff_table)
   ALLOCATE(wyckoff_letters(this%nsites))
   SELECT CASE(this%nsites)
      CASE(0)
         WRITE(0,*) "CRP3D_READ ERR: there aren't Wyckoff sites"
         CALL EXIT(1)
      CASE DEFAULT
         ! do nothing
   END SELECT
   ! Allocate with the correct type (symmetry) all wyckoff sites
   SELECT CASE(system_surface%getsymmlabel())
      CASE("p4mm")
         ALLOCATE(Wyckoffp4mm::this%wyckoffsite(this%nsites))
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Allocated p4mm Wyckoff sites")
#endif
      CASE DEFAULT
         WRITE(0,*) "READ_CRP6D ERR: surface symmetry is not implemented yet"
         WRITE(0,*) "Good luck!"
         CALL EXIT(1)
   END SELECT
   DO i = 1, this%nsites
      CALL AOT_TABLE_OPEN(L=conf,parent=wyckoff_table,thandle=inwyckoff_table,pos=i)
      call aot_table_open( L=conf,parent=inWyckoff_table,thandle=phiList_table,key='phiList' )
      allocate( phiList(aot_table_length(L=conf,thandle=phiList_table)) )
      do k=1,size(phiList)
         call aot_table_open( L=conf,parent=phiList_table,thandle=magnitude_table,pos=k )
         call aot_get_val( L=conf,errCode=iErr,thandle=magnitude_table,pos=1,val=auxReal )
         call aot_get_val( L=conf,errCode=iErr,thandle=magnitude_table,pos=2,val=auxString )
         call angl%read( auxReal,trim(auxString) )
         call angl%to_std()
         phiList(k)=angl%getValue()
         call aot_table_close( L=conf,thandle=magnitude_table )
      enddo
      call aot_table_close( L=conf,thandle=phiList_table )
      CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=inwyckoff_table,key='kind',val=wyckoff_letters(i))
      CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=inwyckoff_table,key='n2dcuts',val=n2dcuts)
      ALLOCATE(cuts2d_files(n2dcuts))
      nthetablocks=aot_table_length(L=conf,thandle=inwyckoff_table)-7
      ALLOCATE(thetablocks_data(nthetablocks))
      call aot_table_open(L=conf,parent=inwyckoff_table,thandle=theta_table,key='thetaFourierTerms')
      do k=1,nThetaBlocks
         call aot_table_open( L=conf,parent=theta_table,thandle=term_table,pos=k )
         call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=1,val=auxString  )
         call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=2,val=auxString2 )
         call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=3,val=auxInt     )
         call thetaTerms%addTerm( irrep=trim(auxString),parity=trim(auxString2),kpoint=auxInt )
         call aot_table_close( L=conf,thandle=term_table )
      enddo
      call aot_table_close(L=conf,thandle=theta_table)
#ifdef DEBUG
      CALL VERBOSE_WRITE( routinename,'Wyckoff Site number: ',i)
      CALL VERBOSE_WRITE( routinename,'Wyckoff Letter: '//trim(wyckoff_letters(i)))
      CALL VERBOSE_WRITE( routinename,'Wyckoff number of cut2ds: ',n2dcuts)
      CALl VERBOSE_WRITE( routinename,'Wyckoff number of theta blocks: ',nthetablocks)
      call verbose_write( routinename,'Wyckoff theta terms irreps: ',thetaTerms%irrepList(:) )
      call verbose_write( routinename,'Wyckoff theta terms parities: ',thetaTerms%parityList(:) )
      call verbose_write( routinename,'Wyckoff theta terms kpoints: ',thetaTerms%kpointList(:) )
#endif
      allocate( phiTerms(nThetaBlocks) )
      DO j = 1, nThetaBlocks
         CALL AOT_TABLE_OPEN(L=conf,parent=inwyckoff_table,thandle=cut2d_table,pos=j)
         CALL AOT_TABLE_OPEN(L=conf,parent=cut2d_table,thandle=files_table,key='files')
         thetablocks_data(j)=aot_table_length(L=conf,thandle=files_table)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,'Wyckoff thetablock: ',j)
         CALL VERBOSE_WRITE(routinename,'Wyckoff cut2d inside: ',thetablocks_data(j))
#endif
         SELECT CASE(j)
            CASE(1)
               auxint=1
               auxint2=thetablocks_data(1)
            CASE DEFAULT
               auxint=sum(thetablocks_data(1:j-1))+1
               auxint2=sum(thetablocks_data(1:j-1))+thetablocks_data(j)
         END SELECT
         n=0
         DO k = auxint, auxint2
            n=n+1
            CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=files_table,pos=n,val=cuts2d_files(k))
#ifdef DEBUG
            CALL VERBOSE_WRITE(routinename,'Wyckoff pos array: ',k)
            CALL VERBOSE_WRITE(routinename,'Wyckoff cut2d filename: '//trim(cuts2d_files(k)))
#endif
         END DO
         CALL AOT_TABLE_CLOSE(L=conf,thandle=files_table)
         call aot_table_open( L=conf,parent=cut2d_table,thandle=phi_table,key='phiFourierTerms' )
         do k=1,thetaBlocks_data(j)
            call aot_table_open( L=conf,parent=phi_table,thandle=term_table,pos=k )
            call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=1,val=auxString )
            call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=2,val=auxString2 )
            call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=3,val=auxInt )
            call phiTerms(j)%addTerm( irrep=trim(auxString),parity=trim(auxString2),kpoint=auxInt )
            call aot_table_close( L=conf,thandle=term_table )
         enddo
#ifdef DEBUG
         call verbose_write( routinename,'Wyckoff phi terms irreps: ',  phiTerms(j)%irrepList(:)  )
         call verbose_write( routinename,'Wyckoff phi terms parities: ',phiTerms(j)%parityList(:) )
         call verbose_write( routinename,'Wyckoff phi terms kpoints: ', phiTerms(j)%kpointList(:) )
#endif
         call aot_table_close(L=conf,thandle=phi_table)
         CALL AOT_TABLE_CLOSE(L=conf,thandle=cut2d_table)
      END DO
      CALL this%wyckoffsite(i)%INITIALIZE(mynumber=i,letter=wyckoff_letters(i),nphipoints=thetablocks_data(:),&
                                          filenames=cuts2d_files(:),phiTerms=phiTerms,thetaTerms=thetaTerms,  &
                                          phiList=phiList(:) )
      CALL AOT_TABLE_CLOSE(L=conf,thandle=inwyckoff_table)
      deallocate(cuts2d_files)
      deallocate(thetablocks_data)
      deallocate(phiTerms)
      deallocate(phiList)
      call thetaTerms%reboot()
   END DO
   CALL AOT_TABLE_CLOSE(L=conf,thandle=wyckoff_table)
   ! Get kpoints for Fourier interpolation
   call aot_table_open( L=conf,parent=pes_table,thandle=fourier_table,key='fourierTerms' )
   this%totalTerms=aot_table_length( L=conf,thandle=fourier_table )
   allocate( this%kListXY(this%totalTerms,2)    )
   allocate( this%irrepListXY(this%totalTerms)  )
   allocate( this%parityListXY(this%totalTerms) )
   allocate( this%kListAngle(this%totalTerms)    )
   allocate( this%irrepListAngle(this%totalTerms)  )
   allocate( this%parityListAngle(this%totalTerms) )
   do i=1,this%totalTerms
      call aot_table_open( L=conf,parent=fourier_table,thandle=term_table,pos=i )
      call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=1,val=auxString )
      call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=2,val=auxString2 )
      call aot_table_open( L=conf,parent=term_table,thandle=kpoints_table,pos=3 )
      call aot_get_val( L=conf,errCode=iErr,thandle=kpoints_table,pos=1,val=auxInt )
      call aot_get_val( L=conf,errCode=iErr,thandle=kpoints_table,pos=2,val=auxInt2 )
      call aot_table_close( L=conf,thandle=kpoints_table )
      this%kListXY(i,1:2)=[auxInt,auxInt2]
      this%irrepListXY(i)=trim(auxString)
      this%parityListXY(i)=trim(auxString2)
      call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=4,val=auxString )
      call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=5,val=auxString2 )
      call aot_get_val( L=conf,ErrCode=iErr,thandle=term_table,pos=6,val=auxInt )
      this%kListAngle(i)=auxInt
      this%irrepListAngle(i)=trim(auxString)
      this%parityListAngle(i)=trim(auxString2)
      call aot_table_close( L=conf,thandle=term_table )
   enddo
   call aot_table_close( L=conf,thandle=fourier_table )
   ! ENDDDDDD!!!!!!
   CALL CLOSE_CONFIG(conf)
   RETURN
END SUBROUTINE READ_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_PHI_CRP6D #######################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES.
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0. Initial PHI value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 09/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_PHI_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER,INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8),DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=18),PARAMETER :: routinename = "PLOT1D_PHI_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_PHI_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:5)=x(1:5)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(6)=xmin
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(6),v,dvdu(:)
   ! cycle for inpoints
   DO i=1, inpoints
      r(6)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(6),v,dvdu(:)
   END DO
   ! Final value
   r(6) = xmax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(6),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_PHI_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_ATOMIC_INTERAC_PHI_CRP6D #######################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the atomic interaction
!! PES
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0. Initial PHI value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_ATOMIC_INTERAC_PHI_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER,INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta
   REAL(KIND=8),DIMENSION(2) :: v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8),DIMENSION(6) :: r, dvdu, ratom
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=18),PARAMETER :: routinename = "PLOT1D_PHI_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_PHI_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:5)=x(1:5)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(6)=xmin
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
   WRITE(11,*) r(6),sum(v),v(1),v(2),ratom(:),dvdu(:)
   ! cycle for inpoints
   DO i=1, inpoints
      r(6)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
      WRITE(11,*) r(6),sum(v),v(1),v(2),ratom(:),dvdu(:)
   END DO
   ! Final value
   r(6) = xmax
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
   WRITE(11,*) r(6),sum(v),v(1),v(2),ratom(:),dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_ATOMIC_INTERAC_PHI_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_PHI_SMOOTH_CRP6D #######################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the smooth PES.
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0. Initial PHI value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 09/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_PHI_SMOOTH_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER,INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8),DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=25),PARAMETER :: routinename = "PLOT1D_PHI_SMOOTH_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_PHI_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:5)=x(1:5)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(6)=xmin
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(6),v,dvdu(:)
   ! cycle for inpoints
   DO i=1, inpoints
      r(6)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(6),v,dvdu(:)
   END DO
   ! Final value
   r(6) = xmax
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(6),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_PHI_SMOOTH_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_THETA_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES.
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial THETA value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_THETA_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=20),PARAMETER :: routinename = "PLOT1D_THETA_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_THETA_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:4)=x(1:4)
   r(6)=x(6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(5)=xmin
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(5),v,dvdu(:)
   ! cycle for inpoints
   DO i=1, inpoints
      r(5)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(5),v,dvdu(:)
   END DO
   ! Final value
   r(5) = xmax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(5),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_THETA_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_ATOMIC_INTERAC_THETA_CRP6D ########################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the atomic PES
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial THETA value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_ATOMIC_INTERAC_THETA_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8), DIMENSION(6) :: r,dvdu,ratom
   REAL(KIND=8),DIMENSION(2) :: v
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=20),PARAMETER :: routinename = "PLOT1D_THETA_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_THETA_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:4)=x(1:4)
   r(6)=x(6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(5)=xmin
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
   WRITE(11,*) r(5),sum(v),v(:),ratom(:),dvdu(:)

   ! cycle for inpoints
   DO i=1, inpoints
      r(5)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
      WRITE(11,*) r(5),sum(v),v(:),ratom(:),dvdu(:)
   END DO
   ! Final value
   r(5) = xmax
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
   WRITE(11,*) r(5),sum(v),v(:),ratom(:),dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_ATOMIC_INTERAC_THETA_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_THETA_SMOOTH_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the smooth PES.
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial THETA value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_THETA_SMOOTH_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=27),PARAMETER :: routinename = "PLOT1D_THETA_SMOOTH_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_THETA_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:4)=x(1:4)
   r(6)=x(6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(5)=xmin
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(5),v,dvdu(:)
   ! cycle for inpoints
   DO i=1, inpoints
      r(5)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(5),v,dvdu(:)
   END DO
   ! Final value
   r(5) = xmax
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(5),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_THETA_SMOOTH_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_R_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES.
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial R value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_R_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=16),PARAMETER :: routinename = "PLOT1D_R_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_R_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = thispes%wyckoffsite(1)%zrcut(1)%getfirstr()
   xmax = thispes%wyckoffsite(1)%zrcut(1)%getlastr()
   r(1:3)=x(1:3)
   r(5:6)=x(5:6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=(xmax-xmin)/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(4)=xmin
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),v,dvdu(:)
   ! cycle for inpoints
   DO i=1, inpoints
      r(4)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(4),v,dvdu(:)
   END DO
   ! Final value
   r(4) = xmax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_R_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_R_SMOOTH_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the smooth PES.
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial R value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_R_SMOOTH_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=23),PARAMETER :: routinename = "PLOT1D_R_SMOOTH_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_R_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = thispes%wyckoffsite(1)%zrcut(1)%getfirstr()
   xmax = thispes%wyckoffsite(1)%zrcut(1)%getlastr()
   r(1:3)=x(1:3)
   r(5:6)=x(5:6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=(xmax-xmin)/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(4)=xmin
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(4),v,dvdu(:)
   ! cycle for inpoints
   DO i=1, inpoints
      r(4)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(4),v,dvdu(:)
   END DO
   ! Final value
   r(4) = xmax
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(4),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_R_SMOOTH_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_Z_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES.
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!> @param[in] L - Length to plot from minimim Z
!
!> @warning
!! - The graph starts always at 0,0, Initial Z value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_Z_CRP6D(thispes,npoints,X,L,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   REAL(KIND=8),INTENT(IN) :: L
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=16),PARAMETER :: routinename = "PLOT1D_Z_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_Z_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = thispes%wyckoffsite(1)%zrcut(1)%getfirstz()
   xmax = xmin+L
   r(1:2)=x(1:2)
   r(4:6)=x(4:6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=(xmax-xmin)/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(3)=xmin
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)
   ! cycle for inpoints
   DO i=1, inpoints
      r(3)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(3),v,dvdu(:)
   END DO
   ! Final value
   r(3) = xmax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_Z_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_Z_SMOOTH_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the smooth PES.
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial Z value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_Z_SMOOTH_CRP6D(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=23),PARAMETER :: routinename = "PLOT1D_Z_SMOOTH_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   SELECT CASE(npoints)
      CASE(: 1)
         WRITE(0,*) "PLOT1D_Z_CRP6D ERR: Less than 2 points"
         CALL EXIT(1)
      CASE DEFAULT
         ! do nothing
   END SELECT
   !
   xmin = thispes%wyckoffsite(1)%zrcut(1)%getfirstz()
   xmax = thispes%wyckoffsite(1)%zrcut(1)%getlastz()
   r(1:2)=x(1:2)
   r(4:6)=x(4:6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=(xmax-xmin)/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(3)=xmin
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)
   ! cycle for inpoints
   DO i=1, inpoints
      r(3)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(3),v,dvdu(:)
   END DO
   ! Final value
   r(3) = xmax
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_Z_SMOOTH_CRP6D
!#######################################################################
! SUBROUTINE: PLOT_XYMAP_CRP6D
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (X,Y) of the PES.
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the file to print the output
!> @param[in] init_point - Initial position to start the scan (a.u. and radians). 6 DIMENSIONAL
!> @param[in] nxpoints - Number of points in X axis (auxiliar cartesian coordinates)
!> @param[in] nypoints - Number of points in Y axis (auxiliar cartesian coordinates)
!> @param[in] Lx - Length of X axis (a.u.)
!> @param[in] Ly - Length of Y axis (a.u.)
!
!> @warning
!! - Z,R,THETA,PHI parameters are taken from @b init_point
!
!> @author A.S. Muzas
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT_XYMAP_CRP6D(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(CRP6D),INTENT(IN) :: thispes
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: init_point ! Initial position to start the scan (in a.u. and radians)
   INTEGER,INTENT(IN) :: nxpoints, nypoints ! number of points in XY plane
   CHARACTER(LEN=*),INTENT(IN) :: filename ! filename
   REAL(KIND=8),INTENT(IN) :: Lx ! Length of X axis
   REAL(KIND=8),INTENT(IN) :: Ly ! Length of X axis
   ! Local variables
   REAL(KIND=8) :: xmin, ymin, xmax, ymax
   REAL(KIND=8),DIMENSION(6) :: r,dvdu
   REAL(KIND=8) :: xdelta, ydelta
   INTEGER :: xinpoints, nxdelta
   INTEGER :: yinpoints, nydelta
   INTEGER :: i, j ! counters
   REAL(KIND=8) :: v ! potential
   ! GABBA, GABBA HEY! ---------
   xmin = init_point(1)
   ymin = init_point(2)
   xmax = init_point(1)+Lx
   ymax = init_point(2)+Ly
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
   r(3:6)=init_point(3:6)
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
   END DO
   r(2) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   ! inpoints in XY
   DO i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      r(3:6) = init_point(3:6)
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
      DO j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
         WRITE(11,*) r(1:2),v,dvdu(:)
      END DO
      r(2) = ymax
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
   END DO
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   r(3:6) = init_point(3:6)
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
   END DO
   r(2) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_XYMAP_CRP6D
!######################################################################################
! SUBROUTINE: PLOT_XYMAP_SMOOTH_CRP6D
!######################################################################################
!> @brief
!! Same as plot_xymap_crp6d but calling get_v_and_derivs_smooth, i.e. potential
!! and derivatives are not corrected. Thus, we get smooth corrugationless potential.
!-------------------------------------------------------------------------------------
SUBROUTINE PLOT_XYMAP_SMOOTH_CRP6D(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(CRP6D),INTENT(IN) :: thispes
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: init_point ! Initial position to start the scan (in a.u. and radians)
   INTEGER,INTENT(IN) :: nxpoints, nypoints ! number of points in XY plane
   CHARACTER(LEN=*),INTENT(IN) :: filename ! filename
   REAL(KIND=8),INTENT(IN) :: Lx ! Length of X axis
   REAL(KIND=8),INTENT(IN) :: Ly ! Length of X axis
   ! Local variables
   REAL(KIND=8) :: xmin, ymin, xmax, ymax
   REAL(KIND=8),DIMENSION(6) :: r,dvdu
   REAL(KIND=8) :: xdelta, ydelta
   INTEGER :: xinpoints, nxdelta
   INTEGER :: yinpoints, nydelta
   INTEGER :: i, j ! counters
   REAL(KIND=8) :: v ! potential
   ! GABBA, GABBA HEY! ---------
   xmin = init_point(1)
   ymin = init_point(2)
   xmax = init_point(1)+Lx
   ymax = init_point(2)+Ly
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
   r(3:6)=init_point(3:6)
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
   END DO
   r(2) = ymax
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   ! inpoints in XY
   DO i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      r(3:6) = init_point(3:6)
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
      DO j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
         WRITE(11,*) r(1:2),v,dvdu(:)
      END DO
      r(2) = ymax
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
   END DO
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   r(3:6) = init_point(3:6)
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
   END DO
   r(2) = ymax
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_XYMAP_SMOOTH_CRP6D
!#######################################################################
! SUBROUTINE: PLOT_RZMAP_CRP6D
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (R,Z) of the PES.
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the file to print the output
!> @param[in] init_point - Initial position to start the scan (a.u. and radians). 6 DIMENSIONAL
!> @param[in] nxpoints - Number of points in R axis
!> @param[in] nypoints - Number of points in Z axis
!> @param[in] Lx - Length of R axis (a.u.)
!> @param[in] Ly - Length of Z axis (a.u.)
!
!> @warning
!! - X,Y,THETA,PHI parameters are taken from @b init_point
!
!> @author A.S. Muzas
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT_RZMAP_CRP6D(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(CRP6D),INTENT(IN) :: thispes
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: init_point
   INTEGER,INTENT(IN) :: nxpoints, nypoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),INTENT(IN) :: Lx
   REAL(KIND=8),INTENT(IN) :: Ly
   ! Local variables
   REAL(KIND=8) :: xmin, ymin, xmax, ymax
   REAL(KIND=8),DIMENSION(6) :: r,dvdu
   REAL(KIND=8) :: xdelta, ydelta
   INTEGER :: xinpoints, nxdelta
   INTEGER :: yinpoints, nydelta
   INTEGER :: i, j ! counters
   REAL(KIND=8) :: v ! potential
   ! GABBA, GABBA HEY! ---------
   xmin = init_point(4)
   ymin = init_point(3)
   xmax = init_point(4)+Lx
   ymax = init_point(3)+Ly
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
   r(4) = xmin
   r(3) = ymin
   r(1:2)=init_point(1:2)
   r(5:6)=init_point(5:6)
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),r(3),v,dvdu(:)
   DO i =1, yinpoints
      r(3) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(4),r(3),v,dvdu(:)
   END DO
   r(3) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),r(3),v,dvdu(:)
   ! inpoints in XY
   DO i = 1, xinpoints
      r(4) = xmin+DFLOAT(i)*xdelta
      r(3) = ymin
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(4),r(3),v,dvdu(:)
      DO j = 1, yinpoints
         r(3) = ymin + DFLOAT(j)*ydelta
         CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
         WRITE(11,*) r(4),r(3),v,dvdu(:)
      END DO
      r(3) = ymax
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(4),r(3),v,dvdu(:)
   END DO
   ! Last point in XY plane
   r(4) = xmax
   r(3) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),r(3),v,dvdu(:)
   DO i =1, yinpoints
      r(3) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(4),r(3),v,dvdu(:)
   END DO
   r(3) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),r(3),v,dvdu(:)
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_RZMAP_CRP6D
!#######################################################################
! SUBROUTINE: PLOT_ATOMIC_INTERAC_RZ_CRP6D
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (R,Z) of the sub of atomic potentials.
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the file to print the output
!> @param[in] init_point - Initial position to start the scan (a.u. and radians). 6 DIMENSIONAL
!> @param[in] nxpoints - Number of points in R axis
!> @param[in] nypoints - Number of points in Z axis
!> @param[in] Lx - Length of R axis (a.u.)
!> @param[in] Ly - Length of Z axis (a.u.)
!
!> @warning
!! - X,Y,THETA,PHI parameters are taken from @b init_point
!
!> @author A.S. Muzas
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT_ATOMIC_INTERAC_RZ_CRP6D(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(CRP6D),INTENT(IN) :: thispes
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: init_point
   INTEGER,INTENT(IN) :: nxpoints, nypoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),INTENT(IN) :: Lx
   REAL(KIND=8),INTENT(IN) :: Ly
   ! Local variables
   REAL(KIND=8) :: xmin, ymin, xmax, ymax
   REAL(KIND=8),DIMENSION(6) :: r
   REAL(KIND=8) :: xdelta, ydelta
   INTEGER :: xinpoints, nxdelta
   INTEGER :: yinpoints, nydelta
   REAL(KIND=8),DIMENSION(2) :: v
   REAL(KIND=8),DIMENSION(6) :: atomicx
   REAL(KIND=8),DIMENSION(6) :: dvdu
   INTEGER :: i, j ! counters
   ! GABBA, GABBA HEY! ---------
   xmin = init_point(4)
   ymin = init_point(3)
   xmax = init_point(4)+Lx
   ymax = init_point(3)+Ly
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
   r(4) = xmin
   r(3) = ymin
   r(1:2)=init_point(1:2)
   r(5:6)=init_point(5:6)
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
   WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   DO i =1, yinpoints
      r(3) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
      WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   END DO
   r(3) = ymax
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
   WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   ! inpoints in XY
   DO i = 1, xinpoints
      r(4) = xmin+DFLOAT(i)*xdelta
      r(3) = ymin
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
      WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
      DO j = 1, yinpoints
         r(3) = ymin + DFLOAT(j)*ydelta
         CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
         WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
      END DO
      r(3) = ymax
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
      WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   END DO
   ! Last point in XY plane
   r(4) = xmax
   r(3) = ymax
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
   WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   DO i =1, yinpoints
      r(3) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
      WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   END DO
   r(3) = ymax
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
   WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_ATOMIC_INTERAC_RZ_CRP6D
!###########################################################
!# FUNCTION: is_allowed_CRP6D
!###########################################################
!> @brief
!! Enquires if the CRP6D potential can be calculated at @b X
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
LOGICAL FUNCTION is_allowed_CRP6D(this,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8) :: zmin, rmin, rmax
   ! Run section
   zmin=this%wyckoffsite(1)%zrcut(1)%getfirstZ()
   rmin=this%wyckoffsite(1)%zrcut(1)%getfirstR()
   rmax=this%wyckoffsite(1)%zrcut(1)%getlastR()
   SELECT CASE(size(x)/=6)
      CASE(.TRUE.)
         WRITE(0,*) "is_allowed_CRP6D ERR: wrong number of dimensions"
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(x(3)<zmin .OR. x(4)<rmin .OR. x(4)>rmax)
      CASE(.TRUE.)
         is_allowed_CRP6D=.FALSE.
      CASE(.FALSE.)
         is_allowed_CRP6D=.TRUE.
   END SELECT
   RETURN
END FUNCTION is_allowed_CRP6D

END MODULE PES_H2LiF001_MOD
