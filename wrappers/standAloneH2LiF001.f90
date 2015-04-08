!#########################################################
! MODULE MATHS
!> @brief
!! Contains some mathematical operations
!##########################################################
MODULE MATHS_MOD
! Initial declarations
IMPLICIT NONE
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
REAL(KIND=8),PARAMETER:: dmass2au = 3671.482934845D0
REAL(KIND=8),PARAMETER:: pmass2au = 1836.15267247D0
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
      CASE("angst")
         ! do nothing
      CASE("au")
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
      case('Kelvin')
         ! do nothing
      case('Celsius')
         this%mag=this%mag+kelvinParam
      case('Fahrenheit')
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
SUBROUTINE MASS_AU(this)
   ! Initial declarations
        IMPLICIT NONE
        ! I/O variables
        CLASS(mass), INTENT(INOUT) :: this
        ! Run section
        IF (this%units.EQ."hmass") THEN
                this%mag = this%mag*hmass2au
        ELSE IF (this%units.EQ."dmass") THEN
                this%mag = this%mag*dmass2au
        ELSE IF (this%units.EQ."pmass") THEN
                this%mag = this%mag*pmass2au
        ELSE IF (this%units.EQ."au") THEN
                RETURN
        ELSE
                WRITE(0,*) "MASS_AU ERR: incorrect units"
                WRITE(0,*) "Supported ones: hmass, dmass, pmass, au"
                CALL EXIT(1)
        END IF
        this%units = "au"
        RETURN
END SUBROUTINE MASS_AU
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
   logical:: initialized=.false.
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
   real(kind=8),public:: norm_s1, norm_s2
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
      READ(10,*) surf%symmlabel
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
      ! Set primitive vectors norms
      surf%norm_s1=DSQRT(DOT_PRODUCT(surf%surf2cart_mtrx(1:2,1),surf%surf2cart_mtrx(1:2,1)))
      surf%norm_s2=DSQRT(DOT_PRODUCT(surf%surf2cart_mtrx(1:2,2),surf%surf2cart_mtrx(1:2,2)))
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
!! Projects 2D point into the C4v unit cell
!
!> @warning
!! - Input/output in cartesian coordinates (r)
!----------------------------------------------------------------
FUNCTION project_unitcell_SURFACE(surf,r)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface),INTENT(IN) :: surf
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: project_unitcell_SURFACE
   REAL(KIND=8), DIMENSION(2) :: aux
   INTEGER :: i ! counters
   ! HEY, HO! LET'S GO !!! ----------------------
   project_unitcell_SURFACE = surf%cart2surf(r)
   FORALL (i=1:2)
      aux(i)=dfloat(int(project_unitcell_SURFACE(i),8))
      project_unitcell_SURFACE(i)=project_unitcell_SURFACE(i)-aux(i)
   END FORALL
   project_unitcell_SURFACE = surf%surf2cart(project_unitcell_SURFACE)
   RETURN
END FUNCTION project_unitcell_SURFACE
!################################################################
! SUBROUTINE: project_iwscell ###################################
!################################################################
!> @brief
!! Projects R into Irreducible WS cell. Needs information from
!! surface main vectors.
!
!> @param[in] surf - Surface specifications
!> @param[in] x - 2D point
!
!> @warning
!! - r is in cartesian coordinates (Input and output)
!! - Only C4v symmetry
!----------------------------------------------------------------
FUNCTION project_iwscell_SURFACE(surf,x)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface) , INTENT(IN) :: surf
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: x
   ! Local variables
   REAL*8,DIMENSION(2) :: r
   REAL*8,DIMENSION(2,2) :: Proj_x, Proj_y, Rot
   REAL*8,DIMENSION(2) :: aux
   REAL*8 :: angle, radius, alpha
   INTEGER :: i ! counters
   REAL(KIND=8),DIMENSION(2) :: project_iwscell_SURFACE
   ! HEY, HO! LET'S GO! ------------------
   ! Go to surface coordinates
   r = x
   r = surf%cart2surf(r)
   FORALL (i=1:2)
      aux(i)=DFLOAT(INT(r(i),8))
      r(i)=r(i)-aux(i)
   END FORALL
   ! Now, r vector is inside the unit cell. Let's define this vector
   ! taking as the origin the center of the cell (in surface units is 0.5,0.5):
   FORALL (i=1:2) r(i)=r(i)-0.5D0
   ! Calculate angle
   IF (r(1).EQ.0.D0) THEN
      angle = PI/2.D0
   ELSE
      angle = DATAN(DABS(r(2))/DABS(r(1))) ! Only angles from 0 to pi rad.
   END IF
   radius = DSQRT(r(1)**2.D0+r(2)**2.D0) ! Radius
   r(1)=radius*DCOS(angle)
   r(2)=radius*DSIN(angle)
   alpha = PI/2.D0 ! 90 deg.
   ! 2D Projectors ========================
   !Rotation
   Rot(1,1)=DCOS(alpha)
   Rot(1,2)=-DSIN(alpha)
   Rot(2,1)=DSIN(alpha)
   Rot(2,2)=DCOS(alpha)
   ! Projector X
   Proj_x(1,1)=1.D0
   Proj_x(1,2)=0.D0
   Proj_x(2,1)=0.D0
   Proj_x(2,2)=-1.D0
   ! Projector Y
   Proj_y(1,1)=-1.D0
   Proj_y(1,2)=0.D0
   Proj_y(2,1)=0.D0
   Proj_y(2,2)=1.D0
   ! Project onto IWS cell
   IF ((angle.LE.(PI/4.D0)).AND.(angle.GE.0.D0)) THEN
      r = MATMUL(Proj_y,r)
      r = MATMUL(Rot, r)
   ELSE IF ((angle.LE.(PI/2.D0)).AND.(angle.GT.(PI/4.D0))) THEN
      r = MATMUL(Proj_y, r)
      r = MATMUL(Proj_x, r)
   ELSE
      WRITE(*,*) "ERR GET_V_XYZ: Incorrect angle value."
      WRITE(*,*) "angle : ", angle
      STOP
   END IF
   FORALL (i=1:2) r(i)=r(i)+0.5D0
   ! Go to cartesian coordinates
   r = surf%surf2cart(r)
   project_iwscell_SURFACE = r
   RETURN
END FUNCTION project_iwscell_SURFACE
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
   CHARACTER(LEN=2) :: irrep="NO"
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
type,abstract :: Fourier2d
   integer(kind=4) :: n
   integer(kind=4) :: nfunc
   real(kind=8),dimension(:,:),allocatable :: xy
   real(kind=8),dimension(:,:),allocatable :: f
   integer(kind=4),dimension(:,:),allocatable :: klist
contains
   procedure,public,non_overridable :: read => read_FOURIER2D
   procedure(interpol_FOURIER2D),public,deferred :: interpol  ! deferred !!!! take a look to interface
   procedure(get_f_and_derivs_FOURIER2D),public,deferred :: get_f_and_derivs ! deferred !!!! take a look to interface
end type Fourier2d

abstract interface
   !###########################################################
   !# SUBROUTINE: INTERPOL_FOURIER2D
   !###########################################################
   !-----------------------------------------------------------
   subroutine interpol_FOURIER2D(this,surf,filename)
      import fourier2d
      import surface
      class(fourier2d),intent(inout) :: this
      class(surface),intent(in) :: surf
      character(len=*),intent(in),optional :: filename
   end subroutine interpol_FOURIER2D
   !###########################################################
   !# SUBROUTINE: GET_F_AND_DERIVS
   !###########################################################
   !-----------------------------------------------------------
   subroutine get_f_and_derivs_FOURIER2D(this,surf,r,v,dvdu)
      import fourier2d
      import surface
      class(fourier2d),intent(in):: this
      class(surface),intent(in):: surf
      real(kind=8),dimension(2),intent(in) :: r
      real(kind=8),dimension(:),intent(out) :: v
      real(kind=8),dimension(:,:),intent(out) :: dvdu
   end subroutine get_f_and_derivs_FOURIER2D
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
subroutine read_FOURIER2D(this,xy,f,klist)
   ! initial declarations
   implicit none
   ! i/o variables
   class(fourier2d),intent(inout) :: this
   real(kind=8),dimension(:,:),intent(in) :: xy
   real(kind=8),dimension(:,:),intent(in) :: f
   integer(kind=4),dimension(:,:),intent(in) :: klist
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
   allocate(this%xy(ndata,2))
   allocate(this%f(nfunc,ndata))
   allocate(this%klist(ndata,2))
   this%f = f
   this%nfunc = nfunc
   this%xy = xy
   this%klist=klist
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
use FOURIER2D_MOD
use UNITS_MOD, only: pi
use MATHS_MOD, only: INV_MTRX
IMPLICIT NONE

TYPE,EXTENDS(Fourier2d) :: Fourierp4mm
   PRIVATE
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: coeff
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE :: termmap
   CONTAINS
      PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_FOURIERP4MM
      PROCEDURE,PUBLIC :: GET_F_AND_DERIVS => GET_F_AND_DERIVS_FOURIERP4MM
      PROCEDURE,PUBLIC :: SET_TERMMAP => SET_TERMMAP_FOURIERP4MM
END TYPE Fourierp4mm
! variables and types, body
CONTAINS
!###########################################################
!# FUNCTION: termfoup4mm
!###########################################################
!
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfoup4mm(id,surf,k,r)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: id
   class(Surface),INTENT(IN) :: surf
   INTEGER(KIND=4),DIMENSION(2),INTENT(IN) :: k
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8) :: g
    INTEGER(KIND=4) ::i ! counters
   ! Run section
   g=2.D0*PI/surf%norm_s1
   SELECT CASE(id)
      CASE(0)
         termfoup4mm=1.D0
      CASE(1)
         termfoup4mm=dcos(g*dfloat(k(1))*r(1))+dcos(g*dfloat(k(1))*r(2))
      CASE(2)
         termfoup4mm=dcos(g*dfloat(k(1))*r(1))*dcos(g*dfloat(k(1))*r(2))
      CASE(3)
         termfoup4mm=dcos(g*dfloat(k(1))*r(1))*dcos(g*dfloat(k(2))*r(2))+&
            dcos(g*dfloat(k(2))*r(1))*dcos(g*dfloat(k(1))*r(2))
      CASE(4)
         termfoup4mm=0.D0
         DO i = 0, -k(1) ! k(1) is negative
            termfoup4mm=termfoup4mm+dcos(g*dfloat(-k(1))*r(1))*dcos(g*dfloat(i)*r(2))&
               +dcos(g*dfloat(i)*r(1))*dcos(g*dfloat(-k(1))*r(2))
         END DO
      CASE DEFAULT
         WRITE(0,*) "termfoup4mm ERR: Incorrect fourier term id: ", id
         CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfoup4mm
!###########################################################
!# FUNCTION: termfoup4mm_dx
!###########################################################
!> @brief
!! ! type brief explanation
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfoup4mm_dx(id,surf,k,r)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: id
   class(Surface),INTENT(IN) :: surf
   INTEGER(KIND=4),DIMENSION(2),INTENT(IN) :: k
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8):: g
   INTEGER(KIND=4) :: i ! counter
   ! Run section
   g=2.D0*PI/surf%norm_s1
   SELECT CASE(id)
      CASE(0)
         termfoup4mm_dx=0.D0
      CASE(1)
         termfoup4mm_dx=-g*dfloat(k(1))*dsin(g*dfloat(k(1))*r(1))
      CASE(2)
         termfoup4mm_dx=-g*dfloat(k(1))*dsin(g*dfloat(k(1))*r(1))*dcos(g*dfloat(k(1))*r(2))
      CASE(3)
         termfoup4mm_dx=-g*dfloat(k(1))*dsin(g*dfloat(k(1))*r(1))*dcos(g*dfloat(k(2))*r(2))-&
            g*dfloat(k(2))*dsin(g*dfloat(k(2))*r(1))*dcos(g*dfloat(k(1))*r(2))
      CASE(4)
         termfoup4mm_dx=0.D0
          DO i = 0, -k(1)
            termfoup4mm_dx=termfoup4mm_dx&
               -g*dfloat(-k(1))*dsin(g*dfloat(-k(1))*r(1))*dcos(g*dfloat(i)*r(2))&
               -g*dfloat(i)*dsin(g*dfloat(i)*r(1))*dcos(g*dfloat(-k(1))*r(2))
         END DO
      CASE DEFAULT
         WRITE(0,*) "termfoup4mm_dx ERR: Incorrect fourier term id: ", id
         CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfoup4mm_dx
!###########################################################
!# FUNCTION: termfoup4mm_dy
!###########################################################
!> @brief
!! ! type brief explanation
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfoup4mm_dy(id,surf,k,r)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: id
   class(Surface),INTENT(IN) :: surf
   INTEGER(KIND=4),DIMENSION(2),INTENT(IN) :: k
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8):: g
   INTEGER(KIND=4) :: i ! counters
   ! Run section
   g=2.D0*PI/surf%norm_s1
   SELECT CASE(id)
      CASE(0)
         termfoup4mm_dy=0.D0
      CASE(1)
         termfoup4mm_dy=-g*dfloat(k(1))*dsin(g*dfloat(k(1))*r(2))
      CASE(2)
         termfoup4mm_dy=-g*dfloat(k(1))*dcos(g*dfloat(k(1))*r(1))*dsin(g*dfloat(k(1))*r(2))
      CASE(3)
         termfoup4mm_dy=-g*dfloat(k(2))*dcos(g*dfloat(k(1))*r(1))*dsin(g*dfloat(k(2))*r(2))-&
            g*dfloat(k(1))*dcos(g*dfloat(k(2))*r(1))*dsin(g*dfloat(k(1))*r(2))
      CASE(4)
         termfoup4mm_dy=0.D0
            DO i = 0, -k(1)
               termfoup4mm_dy=termfoup4mm_dy&
               -g*dfloat(i)*dcos(g*dfloat(-k(1))*r(1))*dsin(g*dfloat(i)*r(2))&
               -g*dfloat(-k(1))*dcos(g*dfloat(i)*r(1))*dsin(g*dfloat(-k(1))*r(2))
            END DO
      CASE DEFAULT
         WRITE(0,*) "termfoup4mm_dy ERR: Incorrect fourier term id: ", id
         CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfoup4mm_dy
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
SUBROUTINE INTERPOL_FOURIERP4MM(this,surf,filename)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourierp4mm),INTENT(INOUT) :: this
   class(Surface),INTENT(IN) :: surf
   CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
   ! Local variables
   INTEGER(KIND=4) :: i,j !counters
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: tmtrx, inv_tmtrx
   ! Run section
   ALLOCATE(tmtrx(this%n,this%n))
   ALLOCATE(inv_tmtrx(this%n,this%n))
   ALLOCATE(this%coeff(this%n,this%nfunc))
   CALL this%SET_TERMMAP()
   DO i = 1, this%n ! loop over points
      DO j = 1, this%n ! loop over terms (one for each k point)
         tmtrx(i,j)=termfoup4mm(this%termmap(j),surf,this%klist(j,:),this%xy(i,:))
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
!# SUBROUTINE: SET_TERMMAP_FOURIERP4MM
!###########################################################
!> @brief
!! Sets the map that identifies the kind of terms that should
!! be used during the fourier interpolation. These terms depend
!! on the kpoints chosen.
!
!> @warning
!! - If the first coordinate of the Kpoint is negative, an average will be done instead for
!!   that order
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 28/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_TERMMAP_FOURIERP4MM(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourierp4mm),INTENT(INOUT)::this
   ! Local variables
   INTEGER(KIND=4) :: i !counters
   CHARACTER(LEN=25),PARAMETER :: routinename="SET_TERMMAP_FOURIERP4MM: "
   ! Run section
   ALLOCATE(this%termmap(this%n))
   DO i = 1, this%n
      SELECT CASE(this%klist(i,1)==0 .AND. this%klist(i,2)==0)
         CASE(.TRUE.)
            this%termmap(i)=0
            CYCLE
         CASE DEFAULT
            ! do nothing
      END SELECT
      !
      SELECT CASE(this%klist(i,1)>0 .AND. this%klist(i,2)==0)
         CASE(.TRUE.)
            this%termmap(i)=1
            CYCLE
         CASE DEFAULT
            ! do nothing
      END SELECT
      !
      SELECT CASE(this%klist(i,1)>0 .AND. this%klist(i,2)==this%klist(i,1))
         CASE(.TRUE.)
            this%termmap(i)=2
            CYCLE
         CASE DEFAULT
            ! do nothing
      END SELECT
      !
      SELECT CASE(this%klist(i,2)<this%klist(i,1) .AND. this%klist(i,2)>0)
         CASE(.TRUE.)
            this%termmap(i)=3
            CYCLE
         CASE DEFAULT
            ! do nothing
      END SELECT
      SELECT CASE(this%klist(i,1)<0)
         CASE(.TRUE.)
            this%termmap(i)=4
            CYCLE
         CASE DEFAULT
            WRITE(0,*) "INTERPOL_FOURIER4MM ERR: Kpoint list has a bad item at position ",i
            CALL EXIT(1)
      END SELECT
   END DO
   RETURN
END SUBROUTINE SET_TERMMAP_FOURIERP4MM
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
SUBROUTINE GET_F_AND_DERIVS_FOURIERP4MM(this,surf,r,v,dvdu)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   class(Fourierp4mm),INTENT(IN):: this
   class(Surface),INTENT(IN):: surf
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
   SELECT CASE(size(v) == this%nfunc .AND. size(dvdu(:,1)) == this%nfunc .AND.&
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
         terms(j)=termfoup4mm(this%termmap(j),surf,this%klist(j,:),r)
         terms_dx(j)=termfoup4mm_dx(this%termmap(j),surf,this%klist(j,:),r)
         terms_dy(j)=termfoup4mm_dy(this%termmap(j),surf,this%klist(j,:),r)
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
TYPE,EXTENDS(Termcalculator) :: term_A
PRIVATE
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_2_A
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_2_A
END TYPE term_A
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_2
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
TYPE,EXTENDS(Fourier1d):: Fourier1d_2
   CONTAINS
      ! Set block
      PROCEDURE,PUBLIC :: SET_IRREP => SET_IRREP_FOURIER1D_2
END TYPE FOURIER1D_2
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_IRREP_FOURIER1D_2
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_IRREP_FOURIER1D_2(this,irrep)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d_2),INTENT(INOUT):: this
   CHARACTER(LEN=2),INTENT(IN) :: irrep
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   INTEGER(KIND=4) :: npar
   ! Run section
   ALLOCATE(this%klist(this%n))
   SELECT CASE(irrep)
      CASE("A")
         ALLOCATE(Term_A::this%term)
         this%irrep=irrep
         SELECT CASE(mod(this%n,2) == 0)
            CASE(.TRUE.) ! case is even
               CALL this%SET_AVERAGE_LASTKPOINT(.TRUE.)
               this%klist(1)=0
               npar=(this%n-2)/2
                DO i = 1, npar
                   this%klist(i+1)=2*i
                   this%klist(i+1+npar)=-2*i
                END DO
                this%klist(this%n)=npar+1
                CALL this%SET_LASTKPOINT(npar+1)
            CASE(.FALSE.) ! case is odd
               this%klist(1)=0
               npar=(this%n-1)/2
                DO i = 1, npar
                   this%klist(i+1)=2*i
                   this%klist(i+1+npar)=-2*i
                END DO
                CALL this%SET_LASTKPOINT(npar)
         END SELECT
      CASE DEFAULT
         WRITE(0,*) "SET_IRREP_FOURIER1D_2 ERR: irrep used is not implemented or does not exist"
         WRITE(0,*) "List of irreps implemented: A"
         CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE SET_IRREP_FOURIER1D_2
!###########################################################
!# FUNCTION: termfou1d_2_A1
!###########################################################
REAL(KIND=8) FUNCTION termfou1d_2_A(this,kpoint,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_A),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   SELECT CASE(this%getaveragelast() .AND. kpoint == this%getlastkpoint())
      CASE(.TRUE.)
         termfou1d_2_A=dcos(dfloat(kpoint)*x)+dsin(dfloat(kpoint)*x)
      CASE(.FALSE.)
         SELECT CASE(kpoint)
            CASE(0)
               termfou1d_2_A=1.D0
            CASE(: -1)
               termfou1d_2_A=dsin(dfloat(-kpoint)*x)
            CASE(1 :)
               termfou1d_2_A=dcos(dfloat(kpoint)*x)
            CASE DEFAULT
               WRITE(0,*) "Termfou1d_2_A ERR: Something went really wrong with kpoints of this interpolation"
               CALL EXIT(1)
         END SELECT
   END SELECT
   RETURN
END FUNCTION termfou1d_2_A
!###########################################################
!# FUNCTION: termfou1d_dx_2_A1
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_dx_2_A(this,kpoint,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_A),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   SELECT CASE(this%getaveragelast() .AND. kpoint == this%getlastkpoint())
      CASE(.TRUE.)
         termfou1d_dx_2_A=(-dsin(dfloat(kpoint)*x)+dcos(dfloat(kpoint)*x))*dfloat(kpoint)
      CASE(.FALSE.)
         SELECT CASE(kpoint)
            CASE(0)
               termfou1d_dx_2_A=0.D0
            CASE(: -1)
               termfou1d_dx_2_A=dfloat(-kpoint)*dcos(dfloat(-kpoint)*x)
            CASE(1 :)
               termfou1d_dx_2_A=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
            CASE DEFAULT
               WRITE(0,*) "Termfou1d_2_A ERR: Something went really wrong with kpoints of this interpolation"
               CALL EXIT(1)
         END SELECT
   END SELECT
   RETURN
END FUNCTION termfou1d_dx_2_A
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
TYPE,EXTENDS(Fourier1d):: Fourier1d_4mm
   CONTAINS
      ! Set block
      PROCEDURE,PUBLIC :: SET_IRREP => SET_IRREP_FOURIER1D_4MM
END TYPE FOURIER1D_4MM
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_IRREP_FOURIER1D_4MM
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_IRREP_FOURIER1D_4MM(this,irrep)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d_4mm),INTENT(INOUT):: this
   CHARACTER(LEN=2),INTENT(IN) :: irrep
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   ! Run section
   SELECT CASE(irrep)
      CASE("A1")
         ALLOCATE(Term_A1::this%term)
         this%irrep=irrep
         ALLOCATE(this%klist(this%n))
         FORALL (i=1:this%n) this%klist(i)=(i-1)*4
      CASE DEFAULT
         WRITE(0,*) "SET_IRREP_FOURIER1D_4MM ERR: irrep used is not implemented or does not exist"
         WRITE(0,*) "List of irreps implemented: A1"
         CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE SET_IRREP_FOURIER1D_4MM
!###########################################################
!# FUNCTION: termfou1d_4mm_A1
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_4mm_A1(this,kpoint,x)
   ! Initial declarations
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
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_A1),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   termfou1d_dx_4mm_A1=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
   RETURN
END FUNCTION termfou1d_dx_4mm_A1
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
TYPE,EXTENDS(Termcalculator) :: term_Ap
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_M_Ap
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_M_Ap
END TYPE term_Ap
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_M
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
TYPE,EXTENDS(FOURIER1D):: Fourier1d_M
   CONTAINS
      ! Set block
      PROCEDURE,PUBLIC :: SET_IRREP => SET_IRREP_FOURIER1D_M
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
SUBROUTINE SET_IRREP_FOURIER1D_M(this,irrep)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d_M),INTENT(INOUT):: this
   CHARACTER(LEN=2),INTENT(IN) :: irrep
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   ! Run section
   SELECT CASE(irrep)
      CASE("Ap")
         ALLOCATE(term_Ap::this%term)
         this%irrep=irrep
         ALLOCATE(this%klist(this%n))
         FORALL(i=1:this%n) this%klist(i)=i-1

      CASE DEFAULT
         WRITE(0,*) "SET_IRREP_FOURIER1D_M ERR: irrep used is not implemented or does not exist"
         WRITE(0,*) "List of irreps implemented: Ap"
         CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE SET_IRREP_FOURIER1D_M
!###########################################################
!# FUNCTION: termfou1d_M_Ap
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_M_Ap(this,kpoint,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(term_Ap),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   termfou1d_M_Ap=dcos(dfloat(kpoint)*(x+this%getshift()))
   RETURN
END FUNCTION termfou1d_M_Ap
!###########################################################
!# FUNCTION: termfou1d_dx_M_Ap
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_dx_M_Ap(this,kpoint,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(term_Ap),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   termfou1d_dx_M_Ap=-dsin(dfloat(kpoint)*(x+this%getshift()))*dfloat(kpoint)
   RETURN
END FUNCTION termfou1d_dx_M_Ap
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
TYPE,EXTENDS(Termcalculator) :: term_Ap
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_M45_Ap
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_M45_Ap
END TYPE term_Ap
!/////////////////////////////////////////////////////////////////
! TYPE: term_expanded_m45m1352_A1
!> @brief
!! Child class of abstract termcalculator. Strcture to avoid unnecessary switches
!----------------------------------------------------------------
TYPE,EXTENDS(Termcalculator) :: term_expanded_m45m1352_A1
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_expanded_m45m1352_A1
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_expanded_m45m1352_A1
END TYPE term_expanded_m45m1352_A1
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_M45
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
TYPE,EXTENDS(FOURIER1D):: Fourier1d_M45
   CONTAINS
      ! Set block
      PROCEDURE,PUBLIC :: SET_IRREP => SET_IRREP_FOURIER1D_M45
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
SUBROUTINE SET_IRREP_FOURIER1D_M45(this,irrep)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d_M45),INTENT(INOUT):: this
   CHARACTER(LEN=2),INTENT(IN) :: irrep
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   ! Run section
   SELECT CASE(irrep)
      CASE("Ap")
         ALLOCATE(term_Ap::this%term)
         this%irrep=irrep
         ALLOCATE(this%klist(this%n))
         FORALL(i=1:this%n) this%klist(i)=i-1
      CASE("A1") ! expanded symmetry
         ALLOCATE(term_expanded_m45m1352_A1::this%term)
         this%irrep=irrep
         ALLOCATE(this%klist(this%n))
         FORALL(i=1:this%n) this%klist(i)=(i-1)*2
      CASE DEFAULT
         WRITE(0,*) "SET_IRREP_FOURIER1D_M45 ERR: irrep used is not implemented or does not exist"
         WRITE(0,*) "List of irreps implemented: Ap, A1"
         WRITE(0,*) "Be careful with A1, it is an expanded symmetry flag. In this case we're using symmetry m45m1352_A1"
         CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE SET_IRREP_FOURIER1D_M45
!###########################################################
!# FUNCTION: termfou1d_expanded_m45m1352_A1
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_expanded_m45m1352_A1(this,kpoint,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_expanded_m45m1352_A1),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   SELECT CASE(mod(kpoint,4))
      CASE(0)
         termfou1d_expanded_m45m1352_A1=dcos(dfloat(kpoint)*x)
      CASE(2)
         termfou1d_expanded_m45m1352_A1=dsin(dfloat(kpoint)*x)
      CASE DEFAULT
          WRITE(0,*) "termfou1d_expanded_m45m1352_A1 ERR: Unclassificable kpoint"
          CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfou1d_expanded_m45m1352_A1
!###########################################################
!# FUNCTION: termfou1d_dx_expanded_m45m1352_A1
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_dx_expanded_m45m1352_A1(this,kpoint,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_expanded_m45m1352_A1),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   SELECT CASE(mod(kpoint,4))
      CASE(0)
         termfou1d_dx_expanded_m45m1352_A1=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
      CASE(2)
         termfou1d_dx_expanded_m45m1352_A1=dfloat(kpoint)*dcos(dfloat(kpoint)*x)
      CASE DEFAULT
          WRITE(0,*) "termfou1d_dx_expanded_m45m1352_A1 ERR: Unclassificable kpoint"
          CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfou1d_dx_expanded_m45m1352_A1
!###########################################################
!# FUNCTION: termfou1d_M45_Ap
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_M45_Ap(this,kpoint,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(term_Ap),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   SELECT CASE(mod(kpoint,4))
      CASE(0)
         termfou1d_M45_Ap=dcos(dfloat(kpoint)*x)
      CASE(1)
         termfou1d_M45_Ap=dcos(dfloat(kpoint)*x)+dsin(dfloat(kpoint)*x)
      CASE(2)
         termfou1d_M45_Ap=dsin(dfloat(kpoint)*x)
      CASE(3)
         termfou1d_M45_Ap=dcos(dfloat(kpoint)*x)-dsin(dfloat(kpoint)*x)
      CASE DEFAULT
          WRITE(0,*) "termfou1d_m45_Ap ERR: Unclassificable kpoint"
          CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfou1d_M45_Ap
!###########################################################
!# FUNCTION: termfou1d_dx_M45_Ap
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_dx_M45_Ap(this,kpoint,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(term_Ap),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   SELECT CASE(mod(kpoint,4))
      CASE(0)
         termfou1d_dx_M45_Ap=-dsin(dfloat(kpoint)*x)*dfloat(kpoint)
      CASE(1)
         termfou1d_dx_M45_Ap=(-dsin(dfloat(kpoint)*x)+dcos(dfloat(kpoint)*x))*dfloat(kpoint)
      CASE(2)
         termfou1d_dx_M45_Ap=-dcos(dfloat(kpoint)*x)*dfloat(kpoint)
      CASE(3)
         termfou1d_dx_M45_Ap=(-dsin(dfloat(kpoint)*x)-cos(dfloat(kpoint)*x))*dfloat(kpoint)
      CASE DEFAULT
          WRITE(0,*) "termfou1d_m45_Ap ERR: Unclassificable kpoint"
          CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfou1d_dx_M45_Ap
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
TYPE,EXTENDS(Termcalculator) :: term_A1
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_MM2_A1
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_MM2_A1
END TYPE term_A1
!///////////////////////////////////////////////////////////////////////////////
! TYPE: FOURIER1D_MM2
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!-------------------------------------------------------------------------------
TYPE,EXTENDS(FOURIER1D):: Fourier1d_MM2
   CONTAINS
      PROCEDURE,PUBLIC :: SET_IRREP => SET_IRREP_FOURIER1D_MM2
END TYPE FOURIER1D_MM2
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_IRREP_FOURIER1D_MM2
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_IRREP_FOURIER1D_MM2(this,irrep)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d_MM2),INTENT(INOUT):: this
   CHARACTER(LEN=2),INTENT(IN) :: irrep
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   ! Run section
   SELECT CASE(irrep)
      CASE("A1")
         ALLOCATE(Term_A1::this%term)
         this%irrep=irrep
         ALLOCATE(this%klist(this%n))
         FORALL(i=1:this%n) this%klist(i)=(i-1)*2

      CASE DEFAULT
         WRITE(0,*) "SET_IRREP_FOURIER1D_MM2 ERR: irrep used is not implemented or does not exist"
         WRITE(0,*) "List of irreps implemented: A1"
         CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE SET_IRREP_FOURIER1D_MM2
!###########################################################
!# FUNCTION: termfou1d_MM2_A1
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_MM2_A1(this,kpoint,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_A1),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   termfou1d_MM2_A1=dcos(dfloat(kpoint)*x)
   RETURN
END FUNCTION termfou1d_MM2_A1
!###########################################################
!# FUNCTION: termfou1d_dx_MM2_A1
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_dx_MM2_A1(this,kpoint,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_A1),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   termfou1d_dx_MM2_A1=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
   RETURN
END FUNCTION termfou1d_dx_MM2_A1
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
MODULE FOURIER1D_m45m1352_MOD
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
TYPE,EXTENDS(Termcalculator) :: term_A1
PRIVATE
   ! some atributes
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_m45m1352_A1
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_m45m1352_A1
END TYPE term_A1
!///////////////////////////////////////////////////////////////////////////////
! TYPE: FOURIER1D_m45m1352
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!-------------------------------------------------------------------------------
TYPE,EXTENDS(FOURIER1D):: Fourier1d_m45m1352
   CONTAINS
      PROCEDURE,PUBLIC :: SET_IRREP => SET_IRREP_FOURIER1D_m45m1352
END TYPE FOURIER1D_m45m1352
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_IRREP_FOURIER1D_m45m1352
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_IRREP_FOURIER1D_m45m1352(this,irrep)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d_m45m1352),INTENT(INOUT):: this
   CHARACTER(LEN=2),INTENT(IN) :: irrep
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   ! Run section
   SELECT CASE(irrep)
      CASE("A1")
         ALLOCATE(Term_A1::this%term)
         this%irrep=irrep
         ALLOCATE(this%klist(this%n))
         FORALL(i=1:this%n) this%klist(i)=(i-1)*2

      CASE DEFAULT
         WRITE(0,*) "SET_IRREP_FOURIER1D_m45m1352 ERR: irrep used is not implemented or does not exist"
         WRITE(0,*) "List of irreps implemented: A1"
         CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE SET_IRREP_FOURIER1D_m45m1352
!###########################################################
!# FUNCTION: termfou1d_m45m1352_A1
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_m45m1352_A1(this,kpoint,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_A1),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   SELECT CASE(mod(kpoint,4))
      CASE(0)
         termfou1d_m45m1352_A1=dcos(dfloat(kpoint)*x)
      CASE(2)
         termfou1d_m45m1352_A1=dsin(dfloat(kpoint)*x)
      CASE DEFAULT
         WRITE(0,*) "termfou1d_m45m1352_A1 ERR: wring kpoint"
         CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfou1d_m45m1352_A1
!###########################################################
!# FUNCTION: termfou1d_dx_m45m1352_A1
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_dx_m45m1352_A1(this,kpoint,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_A1),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   SELECT CASE(mod(kpoint,4))
      CASE(0)
         termfou1d_dx_m45m1352_A1=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
      CASE(2)
         termfou1d_dx_m45m1352_A1=dfloat(kpoint)*dcos(dfloat(kpoint)*x)
      CASE DEFAULT
         WRITE(0,*) "termfou1d_m45m1352_A1 ERR: wring kpoint"
         CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfou1d_dx_m45m1352_A1
END MODULE FOURIER1D_m45m1352_MOD
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
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Termcalculator) :: term_A
PRIVATE
CONTAINS
   PROCEDURE,PUBLIC:: getvalue => termfou1d_E_A
   PROCEDURE,PUBLIC:: getderiv => termfou1d_dx_E_A
END TYPE term_A
!/////////////////////////////////////////////////
! TYPE: FOURIER1D_E
!> @brief
!! Class to store all information needed for a 1D REAL fourier interpolation
!------------------------------------------------
TYPE,EXTENDS(FOURIER1D):: Fourier1d_E
   CONTAINS
      ! Set block
      PROCEDURE,PUBLIC :: SET_IRREP => SET_IRREP_FOURIER1D_E
END TYPE FOURIER1D_E
!//////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_IRREP_FOURIER1D_E
!###########################################################
!> @brief
!! Sets irrep for this fourier series
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_IRREP_FOURIER1D_E(this,irrep)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier1d_E),INTENT(INOUT):: this
   CHARACTER(LEN=2),INTENT(IN) :: irrep
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   INTEGER(KIND=4) :: npar
   ! Run section
   ALLOCATE(this%klist(this%n))
   SELECT CASE(irrep)
      CASE("A")
         ALLOCATE(Term_A::this%term)
         this%irrep=irrep
         SELECT CASE(mod(this%n,2) == 0)
            CASE(.TRUE.) ! case is even
               CALL this%SET_AVERAGE_LASTKPOINT(.TRUE.)
               this%klist(1)=0
               npar=(this%n-2)/2
                DO i = 1, npar
                   this%klist(i+1)=i
                   this%klist(i+1+npar)=-i
                END DO
                this%klist(this%n)=npar+1
                CALl this%SET_LASTKPOINT(npar+1)
            CASE(.FALSE.) ! case is odd
               this%klist(1)=0
               npar=(this%n-1)/2
                DO i = 1, npar
                   this%klist(i+1)=i
                   this%klist(i+1+npar)=-i
                END DO
                CALL this%SET_LASTKPOINT(npar)
         END SELECT
      CASE DEFAULT
         WRITE(0,*) "SET_IRREP_FOURIER1D_E ERR: irrep used is not implemented or does not exist"
         WRITE(0,*) "List of irreps implemented: A"
         CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE SET_IRREP_FOURIER1D_E
!###########################################################
!# FUNCTION: termfou1d_E_A
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_E_A(this,kpoint,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_A),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   SELECT CASE(this%getaveragelast() .AND. kpoint == this%getlastkpoint())
      CASE(.TRUE.)
         termfou1d_E_A=dcos(dfloat(kpoint)*x)+dsin(dfloat(kpoint)*x)
      CASE(.FALSE.)
         SELECT CASE(kpoint)
            CASE(0)
               termfou1d_E_A=1.D0
            CASE(: -1)
               termfou1d_E_A=dsin(dfloat(-kpoint)*x)
            CASE(1 :)
               termfou1d_E_A=dcos(dfloat(kpoint)*x)
            CASE DEFAULT
               WRITE(0,*) "Termfou1d_E_A ERR: Something went really wrong with kpoints of this interpolation"
               CALL EXIT(1)
         END SELECT
   END SELECT
   RETURN
END FUNCTION termfou1d_E_A
!###########################################################
!# FUNCTION: termfou1d_dx_E_A
!###########################################################
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfou1d_dx_E_A(this,kpoint,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Term_A),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: kpoint
   REAL(KIND=8),INTENT(IN) :: x
   ! Run section
   SELECT CASE(this%getaveragelast() .AND. kpoint == this%getlastkpoint())
      CASE(.TRUE.)
         termfou1d_dx_E_A=(-dsin(dfloat(kpoint)*x)+dcos(dfloat(kpoint)*x))*dfloat(kpoint)
      CASE(.FALSE.)
         SELECT CASE(kpoint)
            CASE(0)
               termfou1d_dx_E_A=0.D0
            CASE(: -1)
               termfou1d_dx_E_A=dfloat(-kpoint)*dcos(dfloat(-kpoint)*x)
            CASE(1 :)
               termfou1d_dx_E_A=-dfloat(kpoint)*dsin(dfloat(kpoint)*x)
            CASE DEFAULT
               WRITE(0,*) "Termfou1d_E_A ERR: Something went really wrong with kpoints of this interpolation"
               CALL EXIT(1)
         END SELECT
   END SELECT
   RETURN
END FUNCTION termfou1d_dx_E_A
END MODULE FOURIER1D_E_MOD
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
TYPE,ABSTRACT:: Wyckoffsitio
   CHARACTER:: id
   INTEGER(KIND=4):: mynumber
   LOGICAL:: is_homonucl=.FALSE.
   REAL(KIND=8):: x,y
   INTEGER(KIND=4):: n2dcuts
   INTEGER(KIND=4):: nphicuts
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE:: nphipoints
   TYPE(Cut2d),DIMENSION(:),ALLOCATABLE:: zrcut
   CONTAINS
      ! Initialization block
      PROCEDURE,PUBLIC,NON_OVERRIDABLE:: INITIALIZE => INITIALIZE_WYCKOFFSITIO
      PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_ID => SET_ID_WYCKOFFSITIO
      PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_HOMONUCL => SET_HOMONUCL_WYCKOFFSITIO
      PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_MYNUMBER => SET_MYNUMBER_WYCKOFFSITIO
      ! Interface procedures
      PROCEDURE(GET_V_AND_DERIVS_WYCKOFFSITIO),PUBLIC,DEFERRED:: GET_V_AND_DERIVS
END TYPE Wyckoffsitio

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
SUBROUTINE INITIALIZE_WYCKOFFSITIO(this,nphipoints,filenames,letter,is_homonucl,mynumber)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Wyckoffsitio),INTENT(INOUT):: this
   INTEGER(KIND=4),DIMENSION(:),INTENT(IN):: nphipoints
   CHARACTER(LEN=*),DIMENSION(:),INTENT(IN):: filenames
   CHARACTER,INTENT(IN):: letter
   LOGICAL,INTENT(IN):: is_homonucl
   INTEGER(KIND=4),INTENT(IN):: mynumber
   ! Local variables
   INTEGER(KIND=4):: i ! counters
   ! Parameters
   CHARACTER(LEN=*),PARAMETER:: routinename="INITIALIZE_WYCKOFFSITIO: "
   ! Run section
   this%mynumber=mynumber
   this%is_homonucl=is_homonucl
   this%id=letter
   this%n2dcuts=sum(nphipoints(:))
   this%nphicuts=size(nphipoints(:))
   ALLOCATE(this%zrcut(this%n2dcuts))
   ALLOCATE(this%nphipoints(this%nphicuts))
   this%nphipoints(:)=nphipoints(:)
   SELECT CASE(this%n2dcuts==size(filenames(:)))
      CASE(.true.)
         ! do nothing
      CASE(.false.)
         WRITE(0,*) "INITIALIZE_WYCKOFFSITIO ERR: mismatch between number of n2dcuts and number of files to open"
         CALL EXIT(1)
   END SELECT
   DO i = 1, this%n2dcuts
      CALL this%zrcut(i)%READ(trim(filenames(i)))
   END DO
   ! All zrcuts should belong to the same XY position (center of mass)
   this%x=this%zrcut(1)%x
   this%y=this%zrcut(1)%y
   RETURN
END SUBROUTINE INITIALIZE_WYCKOFFSITIO

END MODULE WYCKOFF_GENERIC_MOD
!#########################################################
! MODULE: WYCKOFF_P4MM
!> @brief
!! Provides tools to interpolate through wyckoff sites belonging to
!! p4mm wallpaper symmetry
!##########################################################
MODULE WYCKOFF_P4MM_MOD
! Initial declarations
use WYCKOFF_GENERIC_MOD
use FOURIER1D_2_MOD
use FOURIER1D_4MM_MOD
use FOURIER1D_M_MOD
use FOURIER1D_M45_MOD
use FOURIER1D_MM2_MOD
use FOURIER1D_M45M1352_MOD
use FOURIER1D_E_MOD
IMPLICIT NONE
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
TYPE,EXTENDS(Wyckoffsitio) :: Wyckoffp4mm
   CONTAINS
      PROCEDURE,PUBLIC :: GET_V_AND_DERIVS => GET_V_AND_DERIVS_WYCKOFFP4MM
END TYPE Wyckoffp4mm
!/////////////////////////////////////////////////////////////////
CONTAINS
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
SUBROUTINE GET_V_AND_DERIVS_WYCKOFFP4MM(this,x,v,dvdu)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Wyckoffp4mm),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(4),INTENT(IN) ::x
   REAL(KIND=8),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(4),INTENT(OUT) :: dvdu
   ! Local variables
   INTEGER(KIND=4) :: i,j,h ! counters
   REAL(KIND=8) :: z,r,theta,phi
   CLASS(Fourier1d),DIMENSION(:),ALLOCATABLE :: phicut
   CLASS(Fourier1d),ALLOCATABLE :: thetacut
   REAL(KIND=8) :: aux_theta
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: aux
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: f
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: dfdz
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: dfdr
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: dfdphi
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: philist
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: thetalist
   CHARACTER(LEN=2) :: theta_irrep, phi_irrep
   REAL(KIND=8) :: phi_shift
   CHARACTER(LEN=30),PARAMETER :: routinename="GET_V_AND_DERIVS_WYCKOFFP4MM: "
   ! Run section
   z=x(1)
   r=x(2)
   theta=x(3)
   phi=x(4)
   ! CONDITIONS: edit this part to include/edit symmetries
   SELECT CASE(this%id)
      CASE("a" : "b")
         ALLOCATE(Fourier1d_4mm::phicut(this%nphicuts))
         phi_irrep="A1"
         phi_shift=0.D0
         SELECT CASE(this%is_homonucl)
            CASE(.TRUE.)
               ALLOCATE(Fourier1d_mm2::thetacut)
               theta_irrep="A1"
            CASE(.FALSE.)
               ALLOCATE(Fourier1d_m::thetacut)
               theta_irrep="Ap"
         END SELECT

      CASE("c")
         ALLOCATE(Fourier1d_mm2::phicut(this%nphicuts))
         phi_irrep="A1"
         phi_shift=0.D0
         SELECT CASE(this%is_homonucl)
            CASE(.TRUE.)
               ALLOCATE(Fourier1d_mm2::thetacut)
               theta_irrep="A1"
            CASE(.FALSE.)
               ALLOCATE(Fourier1d_m::thetacut)
               theta_irrep="Ap"
         END SELECT

      CASE("f")
         ALLOCATE(Fourier1d_m45::phicut(this%nphicuts))
         !phi_irrep="Ap" ! should be decided later. There are special symmetries depending upon theta
         phi_shift=0.D0
         SELECT CASE(this%is_homonucl)
            CASE(.TRUE.)
               ALLOCATE(Fourier1d_2::thetacut)
               theta_irrep="A"
            CASE(.FALSE.)
               ALLOCATE(Fourier1d_e::thetacut)
               theta_irrep="A"
         END SELECT

      CASE DEFAULT
         WRITE(0,*) "GET_V_AND_DERIVS_WYCKOFFP4MM: Unexpected error with Wyckoff id"
         CALL EXIT(1)
   END SELECT
   ! PHI INTERPOLATION ----------------------------------------------------
   h=0 ! initialize h
   DO i = 1, this%nphicuts ! loop over specific phi cuts (each one for a different theta value)
      ALLOCATE(f(this%nphipoints(i)))
      ALLOCATE(dfdr(this%nphipoints(i)))
      ALLOCATE(dfdz(this%nphipoints(i)))
      ALLOCATE(philist(this%nphipoints(i)))
      ALLOCATE(aux(2,this%nphipoints(i)))
      DO j = 1, this%nphipoints(i) ! loop over number of zrcuts inside
         h=h+1 ! numbering of zrcuts
         f(j)=this%zrcut(h)%interrz%getvalue((/r,z/)) ! storing potential at this site
         dfdr(j)=this%zrcut(h)%interrz%getderivx((/r,z/)) ! storing d/dr at this site
         dfdz(j)=this%zrcut(h)%interrz%getderivy((/r,z/)) ! storing d/dz at this site
         philist(j)=this%zrcut(h)%phi
         aux_theta=this%zrcut(h)%theta
      END DO
      SELECT CASE(this%id=="f" .AND. this%is_homonucl) ! "f" has special simmetry if theta=90 and we've a homonuclear molecule
         CASE(.TRUE.)
            SELECT CASE(aux_theta>=dacos(0.D0)-1.D-6 .AND. aux_theta<=dacos(0.D0)+1.D-6) ! check if theta is pi/2
               CASE(.TRUE.)
                  phi_irrep="A1" ! expanded symmetry m45 -> m45m1352
               CASE(.FALSE.)
                  phi_irrep="Ap"
            END SELECT
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      aux(1,:)=dfdz(:)
      aux(2,:)=dfdr(:)
      CALL phicut(i)%READ(philist,f)
      CALL phicut(i)%ADD_MOREFUNCS(aux)
      CALL phicut(i)%SET_IRREP(phi_irrep)
      CALL phicut(i)%SET_SHIFT(phi_shift)
      CALL phicut(i)%INTERPOL()
      DEALLOCATE(philist)
      DEALLOCATE(f)
      DEALLOCATE(dfdr)
      DEALLOCATE(dfdz)
      DEALLOCATE(aux)
   END DO
   ! THETA INTERPOLATION --------------------------
   ALLOCATE(thetalist(this%nphicuts))
   ALLOCATE(f(this%nphicuts))
   ALLOCATE(dfdz(this%nphicuts))
   ALLOCATE(dfdr(this%nphicuts))
   ALLOCATE(dfdphi(this%nphicuts))
   h=0 ! reboot h
   DO i = 1, this%nphicuts
      DO j = 1, this%nphipoints(i)
         h=h+1
         thetalist(i)=this%zrcut(h)%theta
      END DO
      ALLOCATE(aux(2,3))
      CALL phicut(i)%GET_ALLFUNCS_AND_DERIVS(phi,aux(1,:),aux(2,:))
      f(i)=aux(1,1)
      dfdz(i)=aux(1,2)
      dfdr(i)=aux(1,3)
      dfdphi(i)=aux(2,1)
      DEALLOCATE(aux)
   END DO
   ALLOCATE(aux(3,this%nphicuts))
   aux(1,:)=dfdz(:)
   aux(2,:)=dfdr(:)
   aux(3,:)=dfdphi(:)
   CALL thetacut%READ(thetalist,f)
   CALL thetacut%ADD_MOREFUNCS(aux)
   CALL thetacut%SET_IRREP(theta_irrep)
   CALL thetacut%INTERPOL()
   DEALLOCATE(f)
   DEALLOCATE(dfdr)
   DEALLOCATE(dfdz)
   DEALLOCATE(dfdphi)
   DEALLOCATE(aux)
   DEALLOCATE(thetalist)
   ALLOCATE(aux(2,4))
   CALL thetacut%GET_ALLFUNCS_AND_DERIVS(theta,aux(1,:),aux(2,:))
   v=aux(1,1) ! value of the potential
   dvdu(1)=aux(1,2) ! dvdz
   dvdu(2)=aux(1,3) ! dvdr
   dvdu(3)=aux(2,1) ! dvdtheta
   dvdu(4)=aux(1,4) ! dvdphi
   DEALLOCATE(aux)
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_WYCKOFFP4MM
END MODULE WYCKOFF_P4MM_MOD
!#######################################################
! MODULE: PES_H2LiF001_MOD
!#######################################################
!> @brief
!! CRP6D PES explicit implementation
!#######################################################
module PES_H2LiF001_MOD
! Initial declarations
use PES_MOD, only: PES
use LiF001SURF_MOD, only: LiF001Surf,pi
use LOGISTIC_FUNCTION_MOD, only: Logistic_func
use XEXPONENTIAL_FUNCTION_MOD, only: Xexponential_func
use PES_HLIF001_NS_MOD, only: Pes_HLiF001_NS
use EXTRAPOL_TO_VACUUM_MOD, only: Vacuumpot
use FOURIER_P4MM_MOD, only: Fourierp4mm
use WYCKOFF_P4MM_MOD, only: WyckoffSitio, Wyckoffp4mm
implicit none
! Local module variable, used to simulate SYSTEM_MOD
type(LiF001Surf),private:: sysLiF001Surf
!/////////////////////////////////////////////////
! TYPE: PES_H2LiF001
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
   character(len=30):: extrapol2vac_flag='Xexponential'
   integer(kind=4),dimension(:,:),allocatable:: xyklist
   contains
      ! Initialization block
      procedure,public:: initialize                  => INITIALIZE_PES_H2LiF001
      ! Set block
      procedure,public:: set_smooth                  => SET_SMOOTH_PES_H2LiF001
      ! Get block
      procedure,public:: get_v_and_derivs            => GET_V_AND_DERIVS_PES_H2LiF001
      procedure,public:: get_v_and_derivs_pure       => GET_V_AND_DERIVS_PURE_PES_H2LiF001
      procedure,public:: get_atomicpot_and_derivs    => GET_ATOMICPOT_AND_DERIVS_PES_H2LiF001
      ! Tools block
      procedure,public:: smooth                      => SMOOTH_PES_H2LiF001
      procedure,public:: interpol                    => INTERPOL_PES_H2LiF001
      procedure,public:: extract_vacuumsurf          => EXTRACT_VACUUMSURF_PES_H2LiF001
      procedure,public:: add_vacuumsurf              => ADD_VACUUMSURF_PES_H2LiF001
      procedure,public:: interpol_new_rzgrid         => INTERPOL_NEW_RZGRID_PES_H2LiF001
      ! Plot toolk
      procedure,public:: plot1d_theta                => PLOT1D_THETA_PES_H2LiF001
      procedure,public:: plot1d_atomic_interac_theta => PLOT1D_ATOMIC_INTERAC_THETA_PES_H2LiF001
      procedure,public:: plot1d_phi                  => PLOT1D_PHI_PES_H2LiF001
      procedure,public:: plot1d_atomic_interac_phi   => PLOT1D_ATOMIC_INTERAC_PHI_PES_H2LiF001
      procedure,public:: plot1d_r                    => PLOT1D_R_PES_H2LiF001
      procedure,public:: plot1d_z                    => PLOT1D_Z_PES_H2LiF001
      procedure,public:: plot_xymap                  => PLOT_XYMAP_PES_H2LiF001
      procedure,public:: plot_rzmap                  => PLOT_RZMAP_PES_H2LiF001
      procedure,public:: plot_atomic_interac_rz      => PLOT_ATOMIC_INTERAC_RZ_PES_H2LiF001
      ! Enquire block
      procedure,public:: is_allowed                  => is_allowed_PES_H2LiF001
end type PES_H2LiF001

private initialize_PES_H2LiF001,set_smooth_PES_H2LiF001,get_v_and_derivs_PES_H2LiF001,get_v_and_derivs_pure_PES_H2LiF001,&
        get_atomicPot_and_derivs_PES_H2LiF001,smooth_PES_H2LiF001,interpol_PES_H2LiF001,extract_vacuumSurf_PES_H2LiF001,&
        add_vacuumSurf_PES_H2LiF001,interpol_new_RZgrid_PES_H2LiF001,plot1d_theta_PES_H2LiF001,&
        plot1d_atomic_interac_theta_PES_H2LiF001,plot1d_phi_PES_H2LiF001,plot1d_atomic_interac_phi_PES_H2LiF001,&
        plot1d_R_PES_H2LiF001,plot1d_Z_PES_H2LiF001,plot_XYmap_PES_H2LiF001,plot_RZmap_PES_H2LiF001,&
        plot_atomic_interac_RZ_PES_H2LiF001,is_allowed_PES_H2LiF001,from_molecular_to_atomic

contains
!###########################################################
!# SUBROUTINE: GET_ATOMICPOT_AND_DERIVS_PES_H2LiF001
!###########################################################
SUBROUTINE GET_ATOMICPOT_AND_DERIVS_PES_H2LiF001(this,molecx,atomicx,v,dvdu)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),INTENT(IN):: this
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
         WRITE(0,*) "GET_ATOMICPOT_AND_DERIVS_PES_H2LiF001 ERR: wrong number of atomic potentials"
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
END SUBROUTINE GET_ATOMICPOT_AND_DERIVS_PES_H2LiF001
!###########################################################
!# SUBROUTINE: SET_SMOOTH_PES_H2LiF001
!###########################################################
SUBROUTINE SET_SMOOTH_PES_H2LiF001(this,is_smooth)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),INTENT(INOUT):: this
   LOGICAL,INTENT(IN) :: is_smooth
   ! Run section
   this%is_smooth=is_smooth
   RETURN
END SUBROUTINE SET_SMOOTH_PES_H2LiF001
!###########################################################
!# SUBROUTINE: INTERPOL_NEW_RZGRID_PES_H2LiF001
!###########################################################
SUBROUTINE INTERPOL_NEW_RZGRID_PES_H2LiF001(this,nRpoints,nZpoints)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),INTENT(INOUT) :: this
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
END SUBROUTINE INTERPOL_NEW_RZGRID_PES_H2LiF001
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_PURE_PES_H2LiF001
!###########################################################
SUBROUTINE GET_V_AND_DERIVS_PURE_PES_H2LiF001(this,x,v,dvdu)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),TARGET,INTENT(IN):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN):: x
   REAL(KIND=8),INTENT(OUT):: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT):: dvdu
   ! Local variables
   INTEGER(KIND=4):: i ! counters
   REAL(KIND=8):: ma,mb
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: f ! smooth function and derivs
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: xy
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: aux1
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: aux2
   REAL(KIND=8),DIMENSION(6):: atomicx
   REAL(KIND=8),DIMENSION(2):: atomic_v
   REAL(KIND=8),DIMENSION(6):: atomic_dvdu
   REAL(KIND=8),DIMENSION(3):: dvdu_atomicA,dvdu_atomicB
   TYPE(Fourierp4mm):: xyinterpol
   CHARACTER(LEN=*),PARAMETER:: routinename="GET_V_AND_DERIVS_PES_H2LiF001: "
   ! Run section
   SELECT CASE(this%is_smooth)
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         WRITE(0,*) "GET_V_AND_DERIVS_PES_H2LiF001 ERR: smooth the PES first (CALL thispes%SMOOTH())"
         CALl EXIT(1)
   END SELECT
   ALLOCATE(f(5,this%nsites))
   ALLOCATE(xy(this%nsites,2))
   DO i = 1, this%nsites
      CALL this%wyckoffsite(i)%GET_V_AND_DERIVS(x(3:6),f(1,i),f(2:5,i))
      xy(i,1)=this%wyckoffsite(i)%x
      xy(i,2)=this%wyckoffsite(i)%y
   END DO
   ! f(1,:) smooth potential values
   ! f(2,:) smooth dvdz
   ! f(3,:) smooth dvdr
   ! f(4,:) smooth dvdtheta
   ! f(5,:) smooth dvdphi
   CALL xyinterpol%READ(xy,f,this%xyklist)
   CALL xyinterpol%INTERPOL(sysLiF001Surf)
   ALLOCATE(aux1(5))
   ALLOCATE(aux2(5,2))
   CALL xyinterpol%GET_F_AND_DERIVS(sysLiF001Surf,x(1:2),aux1,aux2)
   !--------------------------------------
   ! Results for the real potential
   !-------------------------------------
   CALL this%GET_ATOMICPOT_AND_DERIVS(x,atomicx,atomic_v,atomic_dvdu)
   dvdu_atomicA=atomic_dvdu(1:3)
   dvdu_atomicB=atomic_dvdu(4:6)
   v=aux1(1)+sum(atomic_v)
   ! It does not matter ma and mb values inside this subroutine as long ma=mb for D2 and H2. I let
   ! equations to depend on explicit values of Ma and Mb just in case someone want to
   ! implement an explicit HD PES in the future
   ma=1.d0
   mb=1.d0
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
   dvdu(6)=aux1(5)&
      -(mb/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicA(1)*dsin(x(6))&
      +(mb/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicA(2)*dcos(x(6))&
      +(ma/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicB(1)*dsin(x(6))&
      -(ma/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicB(2)*dcos(x(6))
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_PURE_PES_H2LiF001
!###########################################################
!# SUBROUTINE: GET_V AND DERIVS_PES_H2LiF001
!###########################################################
subroutine get_v_and_derivs_PES_H2LiF001(this,x,v,dvdu,errCode)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   class(PES_H2LiF001),target,intent(in):: this
   real(kind=8),dimension(:),intent(in) :: x
   real(kind=8),intent(out):: v
   real(kind=8),dimension(:),intent(out):: dvdu
   integer(kind=1),optional,intent(out):: errCode
   ! Local variables
   real(kind=8),dimension(6):: geom
   real(kind=8):: zCrpMin,zCrpMax,zvac ! first and last PES_H2LiF001 z value and Z infinity
   real(kind=8):: rCrpMin,rCrpMax ! first and last PES_H2LiF001 r value
   real(kind=8):: vzCrpMax, vzvac ! potentials at zCrpMax and zvac
   real(kind=8),dimension(6):: dvducrp ! derivatives at zCrpMax
   real(kind=8),dimension(6):: dvduvac ! derivatives at vacuum
   real(kind=8):: alpha,beta,gama ! parameters
   type(Xexponential_func):: extrapolfunc
   integer(kind=4):: i !counter
   ! local parameter
   real(kind=8),parameter:: zero=0.D0 ! what we will consider zero (a.u.)
   real(kind=8),parameter:: dz=0.5D0 ! 0.25 Angstroems approx
   character(len=*),parameter:: routinename='GET_V_AND_DERIVS_PES_H2LiF001: '
   ! Run section
   zCrpMin=this%wyckoffSite(1)%zrcut(1)%getFirstZ()
   zCrpMax=this%wyckoffSite(1)%zrcut(1)%getLastZ()
   zvac=this%zvacuum
   rCrpMin=this%wyckoffSite(1)%zrcut(1)%getFirstR()
   rCrpMax=this%wyckoffSite(1)%zrcut(1)%getLastR()
   geom(:)=x(1:6)
   if( geom(4)>rCrpMax ) geom(4)=rCrpMax
   if( geom(4)<rCrpMin ) geom(4)=rCrpMin
   if( geom(3)<zCrpMin ) geom(3)=zCrpMin
   ! *************************************************************************
   ! SWITCH 1: Check if we are inside Z's range
   select case( geom(3) >= zCrpMin .and. geom(3)<= zCrpMax ) !easy
   case(.true.)
      call this%get_v_and_derivs_pure(geom,v,dvdu)
      return
   case(.false.)
      ! do nothing, next switch
   end select
   ! *************************************************************************
   ! SWITCH 2: Check if we are inside the extrapolation region
   select case(geom(3)>zCrpMax .AND. geom(3)<zvac)
   case(.true.) ! uff
      ! Set potential and derivs
      vzvac=this%farpot%getpot(geom(4))
      CALL this%GET_V_AND_DERIVS_PURE([geom(1),geom(2),zCrpMax,geom(4),geom(5),geom(6)],vzCrpMax,dvducrp)
      dvduvac(1:3)=zero
      dvduvac(4)=this%farpot%getderiv(geom(4))
      dvduvac(5:6)=zero
      ! Extrapol potential
      beta=-1.D0/zvac
      alpha=(vzCrpMax-vzvac)/(zCrpMax*dexp(beta*zCrpMax)-zvac*dexp(beta*zvac))
      gama=vzvac-alpha*zvac*dexp(beta*zvac)
      CALL extrapolfunc%READ([alpha,beta])
      v=extrapolfunc%getvalue(geom(3))+gama
      dvdu(3)=extrapolfunc%getderiv(geom(3))
      ! extrapol derivatives
      do i = 1, 6
         select case(i)
         case(3)
            ! skip dvdz
         case default
            beta=-1.D0/zvac
            alpha=(dvducrp(i)-dvduvac(i))/(zCrpMax*dexp(beta*zCrpMax)-zvac*dexp(beta*zvac))
            gama=dvduvac(i)-alpha*zvac*dexp(beta*zvac)
            call extrapolfunc%READ([alpha,beta])
            dvdu(i)=extrapolfunc%getvalue(geom(3))+gama
         end select
      end do
      return

   case(.false.)
      ! do nothing
   end select
   ! *************************************************************************
   ! SWITCH 3: Check if we are inside the Vacuum region
   select case(geom(3)>=zvac) !easy
      case(.true.)
         v=this%farpot%getpot(geom(4))
         dvdu(1:3)=0.D0
         dvdu(4)=this%farpot%getderiv(geom(4))
         dvdu(5:6)=0.D0
      case(.false.) ! this's the last switch!
          errCode=1_1
   end select
   return
end subroutine get_v_and_derivs_PES_H2LiF001
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_SMOOTH_PES_H2LiF001
!###########################################################
SUBROUTINE GET_V_AND_DERIVS_SMOOTH_PES_H2LiF001(this,x,v,dvdu)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),INTENT(IN):: this
   REAL(KIND=8),DIMENSION(6),INTENT(IN):: x
   REAL(KIND=8),INTENT(OUT):: v
   REAL(KIND=8),DIMENSION(6),INTENT(OUT):: dvdu
   ! Local variables
   INTEGER(KIND=4):: i ! counters
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: f ! smooth function and derivs
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: xy
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: aux1
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: aux2
   TYPE(Fourierp4mm):: xyinterpol
   CHARACTER(LEN=*),PARAMETER:: routinename="GET_V_AND_DERIVS_PES_H2LiF001: "
   ! Run section
   SELECT CASE(x(3) >= this%wyckoffsite(1)%zrcut(1)%getlastz())
      CASE(.TRUE.)
         v=this%farpot%getpot(x(4))
         dvdu(1)=0.D0
         dvdu(2)=0.D0
         dvdu(3)=0.D0
         dvdu(4)=this%farpot%getderiv(x(4))
         dvdu(5)=0.D0
         dvdu(6)=0.D0
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   ALLOCATE(f(5,this%nsites))
   ALLOCATE(xy(this%nsites,2))
   DO i = 1, this%nsites
      CALL this%wyckoffsite(i)%GET_V_AND_DERIVS(x(3:6),f(1,i),f(2:5,i))
      xy(i,1)=this%wyckoffsite(i)%x
      xy(i,2)=this%wyckoffsite(i)%y
   END DO
   ! f(1,:) smooth potential values
   ! f(2,:) smooth dvdz
   ! f(3,:) smooth dvdr
   ! f(4,:) smooth dvdtheta
   ! f(5,:) smooth dvdphi
   CALL xyinterpol%READ(xy,f,this%xyklist)
   CALL xyinterpol%INTERPOL(sysLiF001Surf)
   ALLOCATE(aux1(5))
   ALLOCATE(aux2(5,2))
   CALL xyinterpol%GET_F_AND_DERIVS(sysLiF001Surf,x(1:2),aux1,aux2)
   v=aux1(1)
   dvdu(1)=aux2(1,1)
   dvdu(2)=aux2(1,2)
   dvdu(3)=aux1(2)
   dvdu(4)=aux1(3)
   dvdu(5)=aux1(4)
   dvdu(6)=aux1(5)
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_SMOOTH_PES_H2LiF001
!###########################################################
!# SUBROUTINE: INTERPOL_PES_H2LiF001
!###########################################################
SUBROUTINE INTERPOL_PES_H2LiF001(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),INTENT(INOUT)::this
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
END SUBROUTINE INTERPOL_PES_H2LiF001
!###########################################################
!# SUBROUTINE: SMOOTH_PES_H2LiF001
!###########################################################
SUBROUTINE SMOOTH_PES_H2LiF001(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
    CLASS(PES_H2LiF001),INTENT(INOUT):: this
   ! Local variables
   REAL(KIND=8),DIMENSION(6):: molcoord,atomcoord,dummy
   INTEGER(KIND=4):: nr,nz
   INTEGER(KIND=4):: i,j,k,l ! counters
   CHARACTER(LEN=*),PARAMETER:: routinename="SMOOTH_PES_H2LiF001: "
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
END SUBROUTINE SMOOTH_PES_H2LiF001
!###########################################################
!# SUBROUTINE: EXTRACT_VACUUMSURF_PES_H2LiF001
!###########################################################
SUBROUTINE EXTRACT_VACUUMSURF_PES_H2LiF001(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
    CLASS(PES_H2LiF001),INTENT(INOUT) :: this
   ! Local variables
   INTEGER(KIND=4) :: nr,nz
   INTEGER(KIND=4) :: i,j,k,l ! counters
   REAL(KIND=8) :: newpot
   CHARACTER(LEN=*),PARAMETER :: routinename="EXTRACT_VACUUMSURF_PES_H2LiF001: "
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
   this%is_shifted=.TRUE.
   RETURN
END SUBROUTINE EXTRACT_VACUUMSURF_PES_H2LiF001
!###########################################################
!# SUBROUTINE: ADD_VACUUMSURF_PES_H2LiF001
!###########################################################
SUBROUTINE ADD_VACUUMSURF_PES_H2LiF001(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
    CLASS(PES_H2LiF001),INTENT(INOUT) :: this
   ! Local variables
   INTEGER(KIND=4) :: nr,nz
   INTEGER(KIND=4) :: i,j,k,l ! counters
   REAL(KIND=8) :: newpot
   CHARACTER(LEN=*),PARAMETER :: routinename="ADD_VACUUMSURF_PES_H2LiF001: "
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
   this%is_shifted=.FALSE.
   RETURN
END SUBROUTINE ADD_VACUUMSURF_PES_H2LiF001
!###########################################################
!# SUBROUTINE: READ_PES_H2LiF001
!###########################################################
subroutine initialize_PES_H2LiF001(this,filename,tablename)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Pes_H2LiF001),intent(out):: this
   character(len=*),optional,intent(in):: filename,tablename
   ! local variables
   integer(kind=4):: i
   real(kind=8),dimension(14,14):: gridPot1414
   real(kind=8),dimension(14,16):: gridPot1416
   integer(kind=4),save:: invoked ! 0 if this was the first call. 1 if not.
   ! Parameters
   data invoked/0/
   character(len=*),parameter:: routinename="READ_PES_H2LiF001: "
   real(kind=8),dimension(14),parameter:: gridR14=[ 0.7558904532d0,0.9448630664d0,1.1338356797d0,1.2283219864d0,&
                                                    1.3228082930d0,1.4172945997d0,1.5117809063d0,1.6062672130d0,&
                                                    1.8897261329d0,2.3621576661d0,2.8345891993d0,3.3070207325d0,&
                                                    3.7794522658d0,4.3463701056d0 ]
   real(kind=8),dimension(14),parameter:: gridZ14=[ 0.4724315332d0,0.9448630664d0,1.8897261329d0,2.8345891993d0,&
                                                    3.4015070392d0,3.9684248791d0,4.7243153322d0,5.4802057854d0,&
                                                    6.0471236252d0,6.4250688518d0,6.8030140784d0,7.1809593050d0,&
                                                    7.3699319183d0,7.5589045315d0 ]
   real(kind=8),dimension(16),parameter:: gridZ16=[ 0.4724315332d0,0.9448630664d0,1.3983973383d0,1.8897261329d0,&
                                                    2.3054658821d0,2.8345891993d0,3.4015070392d0,3.9684248791d0,&
                                                    4.7243153322d0,5.4802057854d0,6.0471236252d0,6.4250688518d0,&
                                                    6.8030140784d0,7.1809593050d0,7.3699319183d0,7.5589045315d0 ]
   ! Run section -----------------------
   call this%set_pesType('PES_H2LiF001')
   call this%set_alias('H2_on_LiF001')
   call this%set_dimensions(6)
   this%is_homonucl=.true.
   call this%farPot%initialize_direct( x=[ 0.7558904532d0,0.9448630664d0,1.1338356797d0,1.2283219864d0,1.3228082930d0,&
                                           1.4172945997d0,1.5117809063d0,1.6062672130d0,1.7007535196d0,1.8897261329d0,&
                                           2.3621576661d0,2.8345891993d0,3.3070207325d0,3.7794522658d0,4.3463701056d0 ],&
                                       f=[ -0.0412144309d0,-0.1744799612d0,-0.2298169263d0,-0.2422003318d0,-0.2483344814d0,&
                                           -0.2500191955d0,-0.2485072327d0,-0.2446916237d0,-0.2392079371d0,-0.2250493114d0,&
                                           -0.1828462465d0,-0.1423982815d0,-0.1082070069d0,-0.0810442795d0,-0.0564546003d0 ],&
                                       surfaceEnergy=-7.09306998104d0 )
   call this%dumpFunc%read([3.d0,4.d0])
   this%zvacuum=13.2280829302d0
   allocate( this%wyckoffSite(4) )
   allocate( this%atomicCrp(1) )
   call this%atomicCrp(1)%initialize()
   allocate( this%xyKlist(4,2) )
   this%xyklist(:,1)=[0,1,1,2]
   this%xyklist(:,2)=[0,0,1,0]
   select case( invoked ) ! check if this routine has been initialized before
   case(0)
      call sysLiF001Surf%initialize('dummyString')
      invoked=1
   case default
      ! do nothing
   end select
   ! Create wyckoff sites//////////////////////////////////
   ! Wickoff Top Li -----------> 'a'
   this%wyckoffSite(1)%id='a'
   this%wyckoffSite(1)%myNumber=1
   this%wyckoffSite(1)%is_homonucl=.true.
   this%wyckoffSite(1)%x=0.d0
   this%wyckoffSite(1)%y=0.d0
   this%wyckoffSite(1)%n2dcuts=5
   this%wyckoffSite(1)%nphicuts=3
   allocate( this%wyckoffSite(1)%nPhiPoints(3) )
   this%wyckoffSite(1)%nPhiPoints(:)=[1,2,2]
   allocate( this%wyckoffSite(1)%zrCut(5) )
   this%wyckoffSite(1)%zrCut(:)%x=0.d0
   this%wyckoffSite(1)%zrCut(:)%y=0.d0
   ! Reading zrcuts
   this%wyckoffSite(1)%zrcut(1)%alias='Top_Li_0_0'
   this%wyckoffSite(1)%zrCut(1)%theta=0.d0
   this%wyckoffSite(1)%zrCut(1)%phi=0.d0
   gridPot1416(:,:)=reshape( [ -2.3486329404d0, 1.0844017627d0, 8.8266781296d0, 9.1272988309d0, 6.5171540893d0, 1.6638643184d0,&
                               -0.9649307703d0,-2.4497027752d0,-4.5344580512d0,-5.7406206931d0,-6.3350386134d0,-6.6717466581d0,&
                               -6.8772623685d0,-7.0181287846d0,-5.8769513827d0,-5.7776624228d0,-5.5109630136d0,-5.2525253008d0,&
                               -5.0118020103d0,-4.6926259863d0,-4.2606813963d0,-3.6820328918d0, 0.6399452089d0, 1.3764483020d0,&
                               -4.7204565840d0,-5.8608453541d0,-6.4126558964d0,-6.7673050127d0,-6.6586128413d0,-6.7073472242d0,&
                               -6.6589171257d0,-6.6112389184d0,-6.5509718623d0,-6.4789630291d0,-6.3952271187d0,-6.2990585441d0,&
                               -5.9175156555d0,-4.6496075711d0, 2.0317609625d0, 0.0824348596d0,-4.9097362666d0,-6.0776842388d0,&
                               -6.9680495110d0,-7.0663885008d0,-7.0788310874d0,-7.0666236965d0,-7.0459661657d0,-7.0185276494d0,&
                               -6.9853878427d0,-6.9472269758d0,-6.8055553864d0,-6.4702699762d0,-5.9208833637d0,-4.8334825967d0,&
                                0.5489994402d0,-1.2821154148d0,-7.0648608314d0,-7.1814932679d0,-7.2164459187d0,-7.2171338661d0,&
                               -7.2105138426d0,-7.1983207840d0,-7.1817347109d0,-7.1615556564d0,-7.0839759930d0,-6.9074087142d0,&
                               -6.6625938957d0,-6.3011745703d0,-5.6766355873d0,-3.6644919811d0,-7.1108893614d0,-7.2375201869d0,&
                               -7.2847423350d0,-7.2924879904d0,-7.2935801803d0,-7.2897972047d0,-7.2823650212d0,-7.2721420938d0,&
                               -7.2300049498d0,-7.1361504802d0,-7.0179293704d0,-6.8652991322d0,-6.6558110727d0,-6.2631022691d0,&
                               -7.1266893664d0,-7.2574589009d0,-7.3097704631d0,-7.3204310749d0,-7.3246987741d0,-7.3243661927d0,&
                               -7.3206776629d0,-7.3145136985d0,-7.2867102614d0,-7.2261021715d0,-7.1583944794d0,-7.0828874780d0,&
                               -6.9915426198d0,-6.8415160713d0,-7.1319312902d0,-7.2641987272d0,-7.3183562080d0,-7.3300766703d0,&
                               -7.3355008708d0,-7.3364269538d0,-7.3341054989d0,-7.3294243698d0,-7.3068566091d0,-7.2583361073d0,&
                               -7.2088160238d0,-7.1609661972d0,-7.1120951068d0,-7.0441860284d0,-7.1342082784d0,-7.2671530054d0,&
                               -7.3221252188d0,-7.3343068852d0,-7.3402286715d0,-7.3416912946d0,-7.3399471717d0,-7.3358863712d0,&
                               -7.3154735909d0,-7.2718069400d0,-7.2295404384d0,-7.1928925411d0,-7.1615273594d0,-7.1278002985d0,&
                               -7.1347778929d0,-7.2679750878d0,-7.3232067515d0,-7.3355295352d0,-7.3416001563d0,-7.3432193316d0,&
                               -7.3416376406d0,-7.3377473570d0,-7.3178946365d0,-7.2753620697d0,-7.2346438172d0,-7.2003467743d0,&
                               -7.1728178546d0,-7.1471344861d0,-7.1348653563d0,-7.2681114278d0,-7.3234217350d0,-7.3357867805d0,&
                               -7.3419011333d0,-7.3435644078d0,-7.3420293884d0,-7.3381850415d0,-7.3184723359d0,-7.2761819472d0,&
                               -7.2357139575d0,-7.2017344288d0,-7.1747170598d0,-7.1503221226d0,-7.1348403668d0,-7.2681033430d0,&
                               -7.3234437846d0,-7.3358253673d0,-7.3419547873d0,-7.3436338640d0,-7.3421146469d0,-7.3382868371d0,&
                               -7.3186222731d0,-7.2764006057d0,-7.2359888425d0,-7.2020449606d0,-7.1750764682d0,-7.1508513129d0,&
                               -7.1347775254d0,-7.2680559364d0,-7.3234151201d0,-7.3358058902d0,-7.3419456000d0,-7.3436345990d0,&
                               -7.3421253042d0,-7.3383085192d0,-7.3186737222d0,-7.2764921115d0,-7.2360924756d0,-7.2021552086d0,&
                               -7.1751793663d0,-7.1509336313d0,-7.1346952069d0,-7.2679905226d0,-7.3233588937d0,-7.3357551761d0,&
                               -7.3419007658d0,-7.3435956447d0,-7.3420925973d0,-7.3382820597d0,-7.3186663723d0,-7.2765101187d0,&
                               -7.2361288574d0,-7.2021702758d0,-7.1751576842d0,-7.1508505779d0,-7.1346639700d0,-7.2679545082d0,&
                               -7.3233258193d0,-7.3357235717d0,-7.3418706313d0,-7.3435677152d0,-7.3420668727d0,-7.3382578052d0,&
                               -7.3186487326d0,-7.2765012988d0,-7.2361237125d0,-7.2021577810d0,-7.1751290197d0,-7.1507881040d0,&
                               -7.1346239132d0,-7.2679181264d0,-7.3232905399d0,-7.3356901298d0,-7.3418379244d0,-7.3435364783d0,&
                               -7.3420374733d0,-7.3382295082d0,-7.3186248456d0,-7.2764843941d0,-7.2361093803d0,-7.2021383039d0,&
                               -7.1750937403d0,-7.1507212202d0 ],shape( gridPot1416 ) )
   call this%wyckoffSite(1)%zrCut(1)%interRZ%read( gridR14,gridZ16,gridPot1416 )
   ! Second cut2d
   this%wyckoffSite(1)%zrcut(2)%alias='Top_Li_45_0'
   this%wyckoffSite(1)%zrCut(2)%theta=0.785398163397d0
   this%wyckoffSite(1)%zrCut(2)%phi=0.d0
   gridPot1414(:,:)=reshape( [  -4.3930691146d0,-4.5492390478d0,-4.7432240368d0,-4.8587091594d0,&
                                -4.9825143292d0,-5.1107735173d0,-5.2401098281d0,-5.3680834739d0,&
                                -5.7370930050d0,-6.2247672104d0,-6.5372845760d0,-6.7355163171d0,&
                                -6.8625587350d0,-6.9560923830d0,-6.0979945325d0,-6.1737356271d0,&
                                -6.1832555399d0,-6.1793263020d0,-6.1744048323d0,-6.1707174050d0,&
                                -6.1697615551d0,-6.1724810052d0,-6.2061441222d0,-6.3326963041d0,&
                                -6.4939343368d0,-6.6527818258d0,-6.7889858506d0,-6.9054363779d0,&
                                -6.9945082903d0,-7.1105093734d0,-7.1472392217d0,-7.1501071390d0,&
                                -7.1467232611d0,-7.1389735634d0,-7.1281839615d0,-7.1153139802d0,&
                                -7.0705907862d0,-6.9974831482d0,-6.9453715023d0,-6.9210739508d0,&
                                -6.9218636939d0,-6.9449492526d0,-7.1160849811d0,-7.2459383549d0,&
                                -7.2974362870d0,-7.3077808546d0,-7.3118203405d0,-7.3113694263d0,&
                                -7.3076919213d0,-7.3016881840d0,-7.2755064946d0,-7.2230618998d0,&
                                -7.1741801521d0,-7.1342861869d0,-7.1044027380d0,-7.0809349863d0,&
                                -7.1288891810d0,-7.2609574367d0,-7.3149429307d0,-7.3266031241d0,&
                                -7.3319935152d0,-7.3329199657d0,-7.3306393025d0,-7.3260441669d0,&
                                -7.3040971022d0,-7.2583669767d0,-7.2151321304d0,-7.1791530708d0,&
                                -7.1508605002d0,-7.1258511143d0,-7.1329521864d0,-7.2657940154d0,&
                                -7.3206644331d0,-7.3328001629d0,-7.3386826274d0,-7.3401143811d0,&
                                -7.3383496785d0,-7.3342800582d0,-7.3139334267d0,-7.2708294079d0,&
                                -7.2299994374d0,-7.1960985523d0,-7.1694585988d0,-7.1456619406d0,&
                                -7.1346040686d0,-7.2677667192d0,-7.3229965453d0,-7.3353237390d0,&
                                -7.3414024449d0,-7.3430344825d0,-7.3414726362d0,-7.3376066071d0,&
                                -7.3178681770d0,-7.2756876688d0,-7.2355974622d0,-7.2022717039d0,&
                                -7.1761469760d0,-7.1529993109d0,-7.1349616395d0,-7.2682202058d0,&
                                -7.3235529301d0,-7.3359323079d0,-7.3420624628d0,-7.3437452144d0,&
                                -7.3422318772d0,-7.3384128873d0,-7.3187949949d0,-7.2767273072d0,&
                                -7.2365823441d0,-7.2030360899d0,-7.1766008302d0,-7.1531158063d0,&
                                -7.1349311376d0,-7.2682345381d0,-7.3235962943d0,-7.3359900043d0,&
                                -7.3421341240d0,-7.3438297379d0,-7.3423281604d0,-7.3385201953d0,&
                                -7.3189250876d0,-7.2768500499d0,-7.2366176235d0,-7.2028868876d0,&
                                -7.1761833578d0,-7.1523242258d0,-7.1348921833d0,-7.2681808841d0,&
                                -7.3235540326d0,-7.3359528875d0,-7.3421017846d0,-7.3438018084d0,&
                                -7.3423057434d0,-7.3385014532d0,-7.3189125928d0,-7.2768283678d0,&
                                -7.2365544146d0,-7.2027380529d0,-7.1759014905d0,-7.1518387672d0,&
                                -7.1348006775d0,-7.2681037105d0,-7.3234834739d0,-7.3358852687d0,&
                                -7.3420367383d0,-7.3437411720d0,-7.3422473119d0,-7.3384448592d0,&
                                -7.3188607762d0,-7.2767717739d0,-7.2364698912d0,-7.2025943630d0,&
                                -7.1756651924d0,-7.1514495919d0,-7.1347231364d0,-7.2680191870d0,&
                                -7.3234011554d0,-7.3358036852d0,-7.3419569923d0,-7.3436614260d0,&
                                -7.3421694034d0,-7.3383684206d0,-7.3187876451d0,-7.2766964377d0,&
                                -7.2363776504d0,-7.2024650054d0,-7.1754740959d0,-7.1511522898d0,&
                                -7.1346746273d0,-7.2679761903d0,-7.3233585262d0,-7.3357614235d0,&
                                -7.3419147305d0,-7.3436195317d0,-7.3421267741d0,-7.3383261589d0,&
                                -7.3187464858d0,-7.2766545435d0,-7.2363306112d0,-7.2024051040d0,&
                                -7.1753925124d0,-7.1510321195d0,-7.1346327331d0,-7.2679342961d0,&
                                -7.3233155295d0,-7.3357187943d0,-7.3418724688d0,-7.3435769025d0,&
                                -7.3420841449d0,-7.3382838972d0,-7.3187042241d0,-7.2766111793d0,&
                                -7.2362828371d0,-7.2023477750d0,-7.1753186462d0,-7.1509273840d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(1)%zrCut(2)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! third cut2d
   this%wyckoffSite(1)%zrcut(3)%alias='Top_Li_45_45'
   this%wyckoffSite(1)%zrCut(3)%theta=0.785398163397d0
   this%wyckoffSite(1)%zrCut(3)%phi=0.785398163397d0
   gridPot1414(:,:)=reshape( [ -4.3930823443d0,-4.5492658748d0,-4.7432806308d0,-4.8587903754d0,&
                               -4.9826212697d0,-5.1109150022d0,-5.2403027620d0,-5.3683429241d0,&
                               -5.7369761422d0,-6.2233861707d0,-6.5348999123d0,-6.7319990392d0,&
                               -6.8575891237d0,-6.9488557059d0,-6.0980114372d0,-6.1737826663d0,&
                               -6.1833514556d0,-6.1794549246d0,-6.1745764517d0,-6.1709390035d0,&
                               -6.1700426874d0,-6.1728290213d0,-6.2067321114d0,-6.3337142604d0,&
                               -6.4951518419d0,-6.6510373354d0,-6.7852300695d0,-6.8995230440d0,&
                               -6.9945185801d0,-7.1105380379d0,-7.1472983881d0,-7.1501868851d0,&
                               -7.1468305692d0,-7.1391135783d0,-7.1283636657d0,-7.1155414585d0,&
                               -7.0710111985d0,-6.9983882841d0,-6.9468738147d0,-6.9230687042d0,&
                               -6.9239664903d0,-6.9462244542d0,-7.1160890235d0,-7.2459512171d0,&
                               -7.2974627465d0,-7.3078176040d0,-7.3118703196d0,-7.3114355750d0,&
                               -7.3077775472d0,-7.3017980645d0,-7.2757185382d0,-7.2235749204d0,&
                               -7.1751823062d0,-7.1359446840d0,-7.1067690271d0,-7.0839734205d0,&
                               -7.1288910185d0,-7.2609651540d0,-7.3149572629d0,-7.3266240712d0,&
                               -7.3320225471d0,-7.3329581849d0,-7.3306889141d0,-7.3261084782d0,&
                               -7.3042231524d0,-7.2586852259d0,-7.2157932507d0,-7.1803220669d0,&
                               -7.1527163411d0,-7.1286506779d0,-7.1329525539d0,-7.2657984253d0,&
                               -7.3206728855d0,-7.3328119226d0,-7.3386991646d0,-7.3401364307d0,&
                               -7.3383783430d0,-7.3343168075d0,-7.3140047204d0,-7.2710116846d0,&
                               -7.2303930227d0,-7.1968364787d0,-7.1706915387d0,-7.1476894009d0,&
                               -7.1346073760d0,-7.2677689241d0,-7.3230009552d0,-7.3353296189d0,&
                               -7.3414105298d0,-7.3430451398d0,-7.3414862334d0,-7.3376242468d0,&
                               -7.3179016188d0,-7.2757707222d0,-7.2357767989d0,-7.2026167801d0,&
                               -7.1767503999d0,-7.1540587940d0,-7.1349664170d0,-7.2682216758d0,&
                               -7.3235551351d0,-7.3359352478d0,-7.3420668727d0,-7.3437510943d0,&
                               -7.3422395946d0,-7.3384224421d0,-7.3188130021d0,-7.2767673639d0,&
                               -7.2366624576d0,-7.2031856597d0,-7.1768617504d0,-7.1535861977d0,&
                               -7.1349458373d0,-7.2682352731d0,-7.3235981318d0,-7.3359925767d0,&
                               -7.3421374314d0,-7.3438341478d0,-7.3423340403d0,-7.3385271777d0,&
                               -7.3189375823d0,-7.2768754070d0,-7.2366639276d0,-7.2029670012d0,&
                               -7.1763160229d0,-7.1525579515d0,-7.1348844660d0,-7.2681816190d0,&
                               -7.3235555026d0,-7.3359547249d0,-7.3421047245d0,-7.3438058508d0,&
                               -7.3423108883d0,-7.3385062306d0,-7.3189225151d0,-7.2768482125d0,&
                               -7.2365889590d0,-7.2027924419d0,-7.1759863815d0,-7.1519791497d0,&
                               -7.1347981050d0,-7.2681044455d0,-7.3234849438d0,-7.3358867387d0,&
                               -7.3420396782d0,-7.3437437445d0,-7.3422513544d0,-7.3384503716d0,&
                               -7.3188699636d0,-7.2767883111d0,-7.2364963507d0,-7.2026329498d0,&
                               -7.1757184789d0,-7.1515297054d0,-7.1347227689d0,-7.2680199220d0,&
                               -7.3234022579d0,-7.3358055227d0,-7.3419595647d0,-7.3436647334d0,&
                               -7.3421734458d0,-7.3383735655d0,-7.3187957299d0,-7.2767107700d0,&
                               -7.2363993325d0,-7.2024933024d0,-7.1755090077d0,-7.1511956540d0,&
                               -7.1346757298d0,-7.2679769253d0,-7.3233596286d0,-7.3357628935d0,&
                               -7.3419173030d0,-7.3436224717d0,-7.3421315516d0,-7.3383320388d0,&
                               -7.3187545707d0,-7.2766681408d0,-7.2363504559d0,-7.2024300935d0,&
                               -7.1754211768d0,-7.1510629890d0,-7.1346378780d0,-7.2679350311d0,&
                               -7.3233166319d0,-7.3357206317d0,-7.3418746738d0,-7.3435798425d0,&
                               -7.3420889223d0,-7.3382894096d0,-7.3187119415d0,-7.2766240416d0,&
                               -7.2363015793d0,-7.2023705596d0,-7.1753421658d0,-7.1509483311d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(1)%zrCut(3)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Fourth cut2d
   this%wyckoffSite(1)%zrcut(4)%alias='Top_Li_90_0'
   this%wyckoffSite(1)%zrCut(4)%theta=1.57079632679d0
   this%wyckoffSite(1)%zrCut(4)%phi=0.d0
   gridPot1414(:,:)=reshape( [ -5.0184217276d0,-5.2816556754d0,-5.4847890415d0,-5.5759494181d0,&
                               -5.6628281308d0,-5.7462241099d0,-5.8264578096d0,-5.9036012584d0,&
                               -6.1161872860d0,-6.4275860022d0,-6.6500858953d0,-6.7907248286d0,&
                               -6.8786692718d0,-6.9416892199d0,-6.2526276089d0,-6.4004646326d0,&
                               -6.4832391806d0,-6.5138748882d0,-6.5408786600d0,-6.5656657125d0,&
                               -6.5890941425d0,-6.6116655781d0,-6.6763572531d0,-6.7753040767d0,&
                               -6.8581234589d0,-6.9217192690d0,-6.9668467056d0,-7.0012264345d0,&
                               -7.0178705715d0,-7.1470955318d0,-7.1996522120d0,-7.2112785961d0,&
                               -7.2171221063d0,-7.2189841947d0,-7.2181095607d0,-7.2153695310d0,&
                               -7.2013492959d0,-7.1744274751d0,-7.1527652177d0,-7.1371669666d0,&
                               -7.1257772481d0,-7.1153966662d0,-7.1208634958d0,-7.2533448139d0,&
                               -7.3079951032d0,-7.3201080483d0,-7.3260346120d0,-7.3275802887d0,&
                               -7.3260015377d0,-7.3221906326d0,-7.3030306368d0,-7.2630286286d0,&
                               -7.2260657897d0,-7.1959989616d0,-7.1726399879d0,-7.1519269656d0,&
                               -7.1309383234d0,-7.2640822318d0,-7.3193319026d0,-7.3316918032d0,&
                               -7.3378226932d0,-7.3395252894d0,-7.3380542139d0,-7.3343010053d0,&
                               -7.3150167968d0,-7.2739395033d0,-7.2351939546d0,-7.2031580977d0,&
                               -7.1780414037d0,-7.1557624927d0,-7.1339095064d0,-7.2672309140d0,&
                               -7.3226455893d0,-7.3350745786d0,-7.3412620625d0,-7.3430076555d0,&
                               -7.3415637745d0,-7.3378215907d0,-7.3184671910d0,-7.2769338384d0,&
                               -7.2373967092d0,-7.2044329318d0,-7.1784511587d0,-7.1553832397d0,&
                               -7.1349785442d0,-7.2683197965d0,-7.3237454966d0,-7.3361741184d0,&
                               -7.3423564574d0,-7.3440913931d0,-7.3426309749d0,-7.3388656390d0,&
                               -7.3194017263d0,-7.2775390998d0,-7.2374955648d0,-7.2039129288d0,&
                               -7.1773075197d0,-7.1536185371d0,-7.1351152517d0,-7.2684395993d0,&
                               -7.3238465572d0,-7.3362630518d0,-7.3424314261d0,-7.3441505595d0,&
                               -7.3426721341d0,-7.3388869536d0,-7.3193528497d0,-7.2773344060d0,&
                               -7.2370821349d0,-7.2032411511d0,-7.1763524047d0,-7.1523444380d0,&
                               -7.1350259509d0,-7.2683414786d0,-7.3237399842d0,-7.3361513338d0,&
                               -7.3423134607d0,-7.3440256118d0,-7.3425405715d0,-7.3387473062d0,&
                               -7.3191871102d0,-7.2771164825d0,-7.2368010026d0,-7.2028839477d0,&
                               -7.1759099429d0,-7.1517972405d0,-7.1349362825d0,-7.2682466653d0,&
                               -7.3236414960d0,-7.3360502732d0,-7.3422101951d0,-7.3439201412d0,&
                               -7.3424317935d0,-7.3386363232d0,-7.3190665725d0,-7.2769786726d0,&
                               -7.2366437155d0,-7.2027038760d0,-7.1757034117d0,-7.1515587374d0,&
                               -7.1348275045d0,-7.2681433997d0,-7.3235349229d0,-7.3359418627d0,&
                               -7.3420999471d0,-7.3438080558d0,-7.3423178706d0,-7.3385205628d0,&
                               -7.3189456672d0,-7.2768478450d0,-7.2365037006d0,-7.2025528363d0,&
                               -7.1755406121d0,-7.1513805031d0,-7.1347371012d0,-7.2680408691d0,&
                               -7.3234283499d0,-7.3358338197d0,-7.3419904341d0,-7.3436967053d0,&
                               -7.3422050502d0,-7.3384059049d0,-7.3188273343d0,-7.2767239997d0,&
                               -7.2363747104d0,-7.2024201712d0,-7.1754028022d0,-7.1512368133d0,&
                               -7.1346856521d0,-7.2679916250d0,-7.3233776358d0,-7.3357820031d0,&
                               -7.3419375151d0,-7.3436430513d0,-7.3421506612d0,-7.3383511485d0,&
                               -7.3187703729d0,-7.2766652008d0,-7.2363144415d0,-7.2023580649d0,&
                               -7.1753406958d0,-7.1511736044d0,-7.1346474328d0,-7.2679438509d0,&
                               -7.3233276567d0,-7.3357312890d0,-7.3418856986d0,-7.3435904998d0,&
                               -7.3420977422d0,-7.3382971269d0,-7.3187148814d0,-7.2766078719d0,&
                               -7.2362560101d0,-7.2022988984d0,-7.1752804269d0,-7.1511129681d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(1)%zrCut(4)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Fifth cut2d
   this%wyckoffSite(1)%zrcut(5)%alias='Top_Li_90_45'
   this%wyckoffSite(1)%zrCut(5)%theta=1.57079632679d0
   this%wyckoffSite(1)%zrCut(5)%phi=0.785398163397d0
   gridpot1414(:,:)=reshape( [ -5.0184786890d0,-5.2817479162d0,-5.4849353038d0,-5.5761136876d0,&
                               -5.6629935027d0,-5.7463681673d0,-5.8265401280d0,-5.9035718590d0,&
                               -6.1151924818d0,-6.4139604549d0,-6.6235473700d0,-6.7500183359d0,&
                               -6.8211885494d0,-6.8491095844d0,-6.2527040475d0,-6.4006380894d0,&
                               -6.4835713945d0,-6.5143088978d0,-6.5414254900d0,-6.5663378577d0,&
                               -6.5898974827d0,-6.6125960710d0,-6.6775887230d0,-6.7761860605d0,&
                               -6.8564260076d0,-6.9137641425d0,-6.9476279109d0,-6.9587846386d0,&
                               -7.0179179781d0,-7.1472057798d0,-7.1998734430d0,-7.2115773681d0,&
                               -7.2175142216d0,-7.2194894979d0,-7.2187482640d0,-7.2161600090d0,&
                               -7.2027163708d0,-7.1770506419d0,-7.1566264693d0,-7.1416062851d0,&
                               -7.1294477708d0,-7.1155524833d0,-7.1208844430d0,-7.2533918530d0,&
                               -7.3080946939d0,-7.3202432859d0,-7.3262128463d0,-7.3278125444d0,&
                               -7.3262984722d0,-7.3225621683d0,-7.3037049869d0,-7.2644835344d0,&
                               -7.2285959807d0,-7.1997238732d0,-7.1774140928d0,-7.1572916321d0,&
                               -7.1309511857d0,-7.2641105288d0,-7.3193892315d0,-7.3317697118d0,&
                               -7.3379248563d0,-7.3396586894d0,-7.3382250983d0,-7.3345141514d0,&
                               -7.3154041347d0,-7.2747895152d0,-7.2367223590d0,-7.2055085845d0,&
                               -7.1812091956d0,-7.1595711928d0,-7.1339168562d0,-7.2672470837d0,&
                               -7.3226782962d0,-7.3351190453d0,-7.3413201265d0,-7.3430833591d0,&
                               -7.3416607927d0,-7.3379417610d0,-7.3186843795d0,-7.2774068022d0,&
                               -7.2382533359d0,-7.2057687697d0,-7.1802794377d0,-7.1575665171d0,&
                               -7.1349825867d0,-7.2683286163d0,-7.3237627687d0,-7.3361976380d0,&
                               -7.3423862244d0,-7.3441299799d0,-7.3426805865d0,-7.3389262754d0,&
                               -7.3195082993d0,-7.2777621682d0,-7.2378854752d0,-7.2045034905d0,&
                               -7.1780866054d0,-7.1544086476d0,-7.1351178242d0,-7.2684443767d0,&
                               -7.3238572145d0,-7.3362773840d0,-7.3424494332d0,-7.3441737116d0,&
                               -7.3427019011d0,-7.3389229680d0,-7.3194142211d0,-7.2774560463d0,&
                               -7.2372805813d0,-7.2035167711d0,-7.1766680814d0,-7.1524800430d0,&
                               -7.1350156611d0,-7.2683458885d0,-7.3237488040d0,-7.3361627261d0,&
                               -7.3423277930d0,-7.3440443539d0,-7.3425637236d0,-7.3387752357d0,&
                               -7.3192334144d0,-7.2772046809d0,-7.2369377101d0,-7.2030588745d0,&
                               -7.1760749473d0,-7.1517145545d0,-7.1349355475d0,-7.2682503403d0,&
                               -7.3236492133d0,-7.3360605630d0,-7.3422230574d0,-7.3439363109d0,&
                               -7.3424527407d0,-7.3386609453d0,-7.3191073642d0,-7.2770540087d0,&
                               -7.2367572709d0,-7.2028416860d0,-7.1758143946d0,-7.1513963054d0,&
                               -7.1348341194d0,-7.2681467072d0,-7.3235422728d0,-7.3359514175d0,&
                               -7.3421117069d0,-7.3438231230d0,-7.3423373478d0,-7.3385429799d0,&
                               -7.3189827840d0,-7.2769158312d0,-7.2366029237d0,-7.2026671267d0,&
                               -7.1756166832d0,-7.1511680920d0,-7.1347374687d0,-7.2680438091d0,&
                               -7.3234353323d0,-7.3358430070d0,-7.3420014589d0,-7.3437106701d0,&
                               -7.3422234249d0,-7.3384275870d0,-7.3188618787d0,-7.2767868411d0,&
                               -7.2364643788d0,-7.2025193944d0,-7.1754571912d0,-7.1509931653d0,&
                               -7.1346878571d0,-7.2679949325d0,-7.3233842507d0,-7.3357908230d0,&
                               -7.3419485399d0,-7.3436570161d0,-7.3421686684d0,-7.3383720956d0,&
                               -7.3188045498d0,-7.2767262047d0,-7.2364008024d0,-7.2024525106d0,&
                               -7.1753870000d0,-7.1509185641d0,-7.1346459628d0,-7.2679471583d0,&
                               -7.3233342716d0,-7.3357401089d0,-7.3418967234d0,-7.3436044645d0,&
                               -7.3421150144d0,-7.3383180741d0,-7.3187483233d0,-7.2766670383d0,&
                               -7.2363394311d0,-7.2023893018d0,-7.1753212187d0,-7.1508494754d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(1)%zrCut(5)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Wickoff Top F -----------> 'b'
   this%wyckoffSite(2)%id='b'
   this%wyckoffSite(2)%myNumber=2
   this%wyckoffSite(2)%is_homonucl=.true.
   this%wyckoffSite(2)%x=2.7216780628885480d0
   this%wyckoffSite(2)%y=2.7216780628885480d0
   this%wyckoffSite(2)%n2dcuts=5
   this%wyckoffSite(2)%nphicuts=3
   allocate( this%wyckoffSite(2)%nPhiPoints(3) )
   this%wyckoffSite(2)%nPhiPoints(:)=[1,2,2]
   allocate( this%wyckoffSite(2)%zrCut(5) )
   this%wyckoffSite(2)%zrCut(:)%x=2.7216780628885480d0
   this%wyckoffSite(2)%zrCut(:)%y=2.7216780628885480d0
   ! reading zrcuts
   this%wyckoffSite(2)%zrcut(1)%alias='Top_F_0_0'
   this%wyckoffSite(2)%zrCut(1)%theta=0.d0
   this%wyckoffSite(2)%zrCut(1)%phi=0.d0
   gridPot1416(:,:)=reshape( [ 33.4715168470d0,36.7280777391d0,32.7111077495d0,27.9642199020d0,&
                               21.3856928713d0,13.9533148245d0, 9.0618268879d0, 5.5913042852d0,&
                               -0.5049412292d0,-4.6567383910d0,-6.2356988133d0,-6.8413948663d0,&
                               -7.0578498187d0,-7.1242317001d0,-1.9033033846d0, 0.0034275311d0,&
                                3.2801723154d0, 5.7512859620d0, 9.1276728621d0,13.9376318364d0,&
                               20.9998533972d0,27.4752624360d0,35.7943159572d0,13.6062627973d0,&
                               -0.6447695306d0,-4.7270460195d0,-6.2813546984d0,-6.9360302369d0,&
                               -5.9031036726d0,-5.7674850646d0,-5.4007425994d0,-5.1231609799d0,&
                               -4.7715419645d0,-4.3329174513d0,-3.7898079910d0,-3.1189144748d0,&
                                0.0390016364d0,16.3720911156d0,35.6009409801d0,11.3358053823d0,&
                               -1.1717378205d0,-5.3378967244d0,-6.7182110598d0,-6.8361484648d0,&
                               -6.8604724758d0,-6.8476532086d0,-6.8195539394d0,-6.7756550328d0,&
                               -6.7145070953d0,-6.6338603332d0,-6.2405863250d0,-4.7196967512d0,&
                               -0.6132108564d0,13.6361947599d0,35.7370124293d0, 8.7378711640d0,&
                               -6.9382990946d0,-7.0706065884d0,-7.1257691633d0,-7.1377939101d0,&
                               -7.1429116211d0,-7.1424430672d0,-7.1370585561d0,-7.1269477142d0,&
                               -7.0648972132d0,-6.7898784917d0,-6.0360156927d0,-4.0851696691d0,&
                                1.2479442124d0,28.4739811097d0,-7.0611300399d0,-7.1933041336d0,&
                               -7.2486306104d0,-7.2614649448d0,-7.2683367012d0,-7.2710050697d0,&
                               -7.2706780007d0,-7.2681731667d0,-7.2523099529d0,-7.2062387936d0,&
                               -7.1133137144d0,-6.8802009836d0,-6.2704458868d0,-4.1779815604d0,&
                               -7.1104149276d0,-7.2431556959d0,-7.2984219039d0,-7.3109839258d0,&
                               -7.3174554820d0,-7.3196405969d0,-7.3187920550d0,-7.3157947800d0,&
                               -7.2994578674d0,-7.2642814131d0,-7.2290340327d0,-7.1898500645d0,&
                               -7.1240533373d0,-6.9118171633d0,-7.1272883804d0,-7.2604888828d0,&
                               -7.3159524346d0,-7.3284905695d0,-7.3348522452d0,-7.3368418537d0,&
                               -7.3357165893d0,-7.3323687258d0,-7.3146584909d0,-7.2771436770d0,&
                               -7.2422424752d0,-7.2129385631d0,-7.1874631958d0,-7.1537780291d0,&
                               -7.1336412363d0,-7.2670148280d0,-7.3225349738d0,-7.3350418717d0,&
                               -7.3413252714d0,-7.3431884622d0,-7.3418853311d0,-7.3383096217d0,&
                               -7.3196236922d0,-7.2798127805d0,-7.2426720748d0,-7.2122925100d0,&
                               -7.1884785796d0,-7.1664249420d0,-7.1348840985d0,-7.2682492378d0,&
                               -7.3237212420d0,-7.3361833057d0,-7.3424038641d0,-7.3441884113d0,&
                               -7.3427845871d0,-7.3390865025d0,-7.3199011496d0,-7.2788267961d0,&
                               -7.2400467030d0,-7.2081508610d0,-7.1834244449d0,-7.1616776641d0,&
                               -7.1349855266d0,-7.2683315563d0,-7.3237682811d0,-7.3362064578d0,&
                               -7.3423979842d0,-7.3441465171d0,-7.3427019011d0,-7.3389549399d0,&
                               -7.3195733457d0,-7.2779903815d0,-7.2384620721d0,-7.2056768964d0,&
                               -7.1801070833d0,-7.1576940372d0,-7.1349381200d0,-7.2682712874d0,&
                               -7.3236947825d0,-7.3361211994d0,-7.3423005985d0,-7.3440333291d0,&
                               -7.3425718085d0,-7.3388053701d0,-7.3193488073d0,-7.2775648243d0,&
                               -7.2377241456d0,-7.2045200277d0,-7.1784779857d0,-7.1555412618d0,&
                               -7.1348473491d0,-7.2681797816d0,-7.3235911494d0,-7.3360109514d0,&
                               -7.3421826331d0,-7.3439058090d0,-7.3424328960d0,-7.3386543304d0,&
                               -7.3191496259d0,-7.2772366528d0,-7.2371945879d0,-7.2036975778d0,&
                               -7.1772895125d0,-7.1538857047d0,-7.1347558433d0,-7.2680798234d0,&
                               -7.3234820039d0,-7.3358970285d0,-7.3420624628d0,-7.3437790238d0,&
                               -7.3422991285d0,-7.3385124780d0,-7.3189794766d0,-7.2769889624d0,&
                               -7.2368168048d0,-7.2031268607d0,-7.1764615502d0,-7.1526825317d0,&
                               -7.1347099067d0,-7.2680298443d0,-7.3234283499d0,-7.3358411695d0,&
                               -7.3420047664d0,-7.3437187549d0,-7.3422362871d0,-7.3384463292d0,&
                               -7.3189026705d0,-7.2768838593d0,-7.2366698075d0,-7.2029122447d0,&
                               -7.1761499160d0,-7.1522239002d0,-7.1346647050d0,-7.2679820702d0,&
                               -7.3233761658d0,-7.3357878830d0,-7.3419485399d0,-7.3436606910d0,&
                               -7.3421760182d0,-7.3383834879d0,-7.3188310093d0,-7.2767927210d0,&
                               -7.2365433898d0,-7.2027325405d0,-7.1758937732d0,-7.1518435446d0 ],shape( gridPot1416 ) )
   call this%wyckoffSite(2)%zrCut(1)%interRZ%read( gridR14,gridZ16,gridPot1416 )
   ! Second cut2d
   this%wyckoffSite(2)%zrcut(2)%alias='Top_F_45_0'
   this%wyckoffSite(2)%zrCut(2)%theta=0.785398163397d0
   this%wyckoffSite(2)%zrCut(2)%phi=0.d0
   gridPot1414(:,:)=reshape( [ 7.3744382510d0,   5.3167479737d0,   2.8013896828d0,   1.5772557258d0,&
                               0.4404425759d0,  -0.5879814523d0,  -1.5016393216d0,  -2.3030362705d0,&
                              -4.1194461324d0,  -5.7720787303d0,  -6.5066907626d0,  -6.8351525581d0,&
                              -6.9866299699d0,  -7.0681296839d0,  -3.9607838600d0,  -3.8982769324d0,&
                              -3.8325867783d0,  -3.8276432590d0,  -3.8507247753d0,  -3.9053151632d0,&
                              -3.9921137623d0,  -4.1093819622d0,  -4.6009152142d0,  -5.5587819707d0,&
                              -6.2990761838d0,  -6.7339651281d0,  -6.9516990012d0,  -7.0609808376d0,&
                              -6.7478251362d0,  -6.8910640842d0,  -6.9587118749d0,  -6.9776399825d0,&
                              -6.9906194767d0,  -6.9990527119d0,  -7.0041998224d0,  -7.0068615761d0,&
                              -7.0055080984d0,  -6.9933573015d0,  -6.9929074897d0,  -7.0134132458d0,&
                              -7.0455703755d0,  -7.0787094472d0,  -7.0617937327d0,  -7.1940872617d0,&
                              -7.2495136967d0,  -7.2623947028d0,  -7.2693355479d0,  -7.2721314365d0,&
                              -7.2720208211d0,  -7.2698724555d0,  -7.2566155039d0,  -7.2277540537d0,&
                              -7.2007793138d0,  -7.1777385893d0,  -7.1580736578d0,  -7.1380107311d0,&
                              -7.1098879423d0,  -7.2422928218d0,  -7.2970941508d0,  -7.3093669555d0,&
                              -7.3155092378d0,  -7.3173279619d0,  -7.3160773823d0,  -7.3126449953d0,&
                              -7.2948961736d0,  -7.2577554679d0,  -7.2234165308d0,  -7.1949799028d0,&
                              -7.1719230085d0,  -7.1496892992d0,  -7.1267279532d0,  -7.2596098389d0,&
                              -7.3146621658d0,  -7.3269548152d0,  -7.3330412384d0,  -7.3347261950d0,&
                              -7.3332643068d0,  -7.3295471126d0,  -7.3105451389d0,  -7.2704394976d0,&
                              -7.2329629030d0,  -7.2020243810d0,  -7.1774750966d0,  -7.1548195050d0,&
                              -7.1333101249d0,  -7.2665032774d0,  -7.3217944749d0,  -7.3341665028d0,&
                              -7.3403003327d0,  -7.3419981515d0,  -7.3405131113d0,  -7.3367356481d0,&
                              -7.3173297993d0,  -7.2759213945d0,  -7.2367550659d0,  -7.2042899769d0,&
                              -7.1787583831d0,  -7.1558448112d0,  -7.1347264439d0,  -7.2680044873d0,&
                              -7.3233695510d0,  -7.3357684059d0,  -7.3419206104d0,  -7.3436276166d0,&
                              -7.3421414739d0,  -7.3383522509d0,  -7.3188376242d0,  -7.2769896974d0,&
                              -7.2371093294d0,  -7.2038273029d0,  -7.1775842421d0,  -7.1541616921d0,&
                              -7.1348980632d0,  -7.2681955838d0,  -7.3235768172d0,  -7.3359797145d0,&
                              -7.3421352265d0,  -7.3438433351d0,  -7.3423549875d0,  -7.3385609871d0,&
                              -7.3190118160d0,  -7.2770249767d0,  -7.2368958159d0,  -7.2032698156d0,&
                              -7.1766493393d0,  -7.1528243841d0,  -7.1348789536d0,  -7.2681801491d0,&
                              -7.3235657924d0,  -7.3359697922d0,  -7.3421256717d0,  -7.3438334128d0,&
                              -7.3423454326d0,  -7.3385488598d0,  -7.3189871939d0,  -7.2769525805d0,&
                              -7.2367322814d0,  -7.2029684712d0,  -7.1761756405d0,  -7.1521342318d0,&
                              -7.1348186847d0,  -7.2681184102d0,  -7.3235047885d0,  -7.3359095233d0,&
                              -7.3420657703d0,  -7.3437738789d0,  -7.3422844287d0,  -7.3384871210d0,&
                              -7.3189184727d0,  -7.2768566648d0,  -7.2365819766d0,  -7.2027318055d0,&
                              -7.1758228470d0,  -7.1516179038d0,  -7.1347345287d0,  -7.2680371942d0,&
                              -7.3234232050d0,  -7.3358283073d0,  -7.3419838193d0,  -7.3436904579d0,&
                              -7.3422002728d0,  -7.3384029650d0,  -7.3188302743d0,  -7.2767533992d0,&
                              -7.2364478416d0,  -7.2025462214d0,  -7.1755633967d0,  -7.1512434282d0,&
                              -7.1346941045d0,  -7.2679941975d0,  -7.3233794733d0,  -7.3357842081d0,&
                              -7.3419397201d0,  -7.3436459913d0,  -7.3421550711d0,  -7.3383566609d0,&
                              -7.3187836027d0,  -7.2767015826d0,  -7.2363861027d0,  -7.2024675778d0,&
                              -7.1754597636d0,  -7.1510997383d0,  -7.1346485353d0,  -7.2679519358d0,&
                              -7.3233353741d0,  -7.3357401089d0,  -7.3418948859d0,  -7.3436007896d0,&
                              -7.3421095020d0,  -7.3383103567d0,  -7.3187358285d0,  -7.2766497661d0,&
                              -7.2363273038d0,  -7.2023962842d0,  -7.1753686253d0,  -7.1509755256d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(2)%zrCut(2)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Third cut2d
   this%wyckoffSite(2)%zrcut(3)%alias='Top_F_45_45'
   this%wyckoffSite(2)%zrCut(3)%theta=0.785398163397d0
   this%wyckoffSite(2)%zrCut(3)%phi=0.785398163397d0
   gridPot1414(:,:)=reshape( [ 7.3744812477d0, 5.3168402145d0, 2.8015803383d0, 1.5775181527d0,&
                               0.4407941934d0,-0.5875210935d0,-1.5010470327d0,-2.3022891935d0,&
                              -4.1180518630d0,-5.7690799853d0,-6.5014565562d0,-6.8269313665d0,&
                              -6.9744773355d0,-7.0498873187d0,-3.9607515206d0,-3.8982133561d0,&
                              -3.8324555832d0,-3.8274635548d0,-3.8504888446d0,-3.9050020589d0,&
                              -3.9917109897d0,-4.1088751890d0,-4.5999582617d0,-5.5565281346d0,&
                              -6.2944865605d0,-6.7255278505d0,-6.9374490828d0,-7.0368005165d0,&
                              -6.7478137439d0,-6.8910317448d0,-6.9586633658d0,-6.9776921665d0,&
                              -6.9905220910d0,-6.9989365840d0,-7.0040340830d0,-7.0066451225d0,&
                              -7.0050560817d0,-6.9920865098d0,-6.9898841227d0,-7.0070067359d0,&
                              -7.0331186016d0,-7.0534803003d0,-7.0617874853d0,-7.1940758694d0,&
                              -7.2494905447d0,-7.2623634658d0,-7.2692940211d0,-7.2720770475d0,&
                              -7.2719502624d0,-7.2697824197d0,-7.2564380046d0,-7.2272704325d0,&
                              -7.1995985580d0,-7.1750970478d0,-7.1525983758d0,-7.1259330653d0,&
                              -7.1098838999d0,-7.2422865744d0,-7.2970812885d0,-7.3093489483d0,&
                              -7.3154853507d0,-7.3172970924d0,-7.3160369581d0,-7.3125939138d0,&
                              -7.2947973180d0,-7.2575015301d0,-7.2228380964d0,-7.1937462279d0,&
                              -7.1694189095d0,-7.1441272888d0,-7.1267253807d0,-7.2596065315d0,&
                              -7.3146551834d0,-7.3269452603d0,-7.3330291111d0,-7.3347100253d0,&
                              -7.3332440947d0,-7.3295213880d0,-7.3104933223d0,-7.2703072000d0,&
                              -7.2326722159d0,-7.2014393317d0,-7.1763524047d0,-7.1524418237d0,&
                              -7.1333090224d0,-7.2665021749d0,-7.3217922699d0,-7.3341635628d0,&
                              -7.3402962902d0,-7.3419930066d0,-7.3405061289d0,-7.3367264608d0,&
                              -7.3173106897d0,-7.2758703129d0,-7.2366400405d0,-7.2040621311d0,&
                              -7.1783401758d0,-7.1550286087d0,-7.1347257089d0,-7.2680041198d0,&
                              -7.3233691835d0,-7.3357676709d0,-7.3419202429d0,-7.3436268816d0,&
                              -7.3421407389d0,-7.3383507810d0,-7.3188335817d0,-7.2769746301d0,&
                              -7.2370689052d0,-7.2037402070d0,-7.1774177677d0,-7.1538397680d0,&
                              -7.1348833635d0,-7.2681959513d0,-7.3235771847d0,-7.3359800820d0,&
                              -7.3421359615d0,-7.3438440701d0,-7.3423560899d0,-7.3385620896d0,&
                              -7.3190129184d0,-7.2770224043d0,-7.2368818511d0,-7.2032319638d0,&
                              -7.1765673883d0,-7.1526542348d0,-7.1348833635d0,-7.2681805166d0,&
                              -7.3235661599d0,-7.3359701597d0,-7.3421267741d0,-7.3438348828d0,&
                              -7.3423469026d0,-7.3385510648d0,-7.3189901339d0,-7.2769547855d0,&
                              -7.2367286064d0,-7.2029489940d0,-7.1761256614d0,-7.1520210439d0,&
                              -7.1348139073d0,-7.2681187777d0,-7.3235055235d0,-7.3359106258d0,&
                              -7.3420672402d0,-7.3437760839d0,-7.3422870012d0,-7.3384900609d0,&
                              -7.3189228826d0,-7.2768618097d0,-7.2365845491d0,-7.2027244556d0,&
                              -7.1757938150d0,-7.1515414652d0,-7.1347334262d0,-7.2680379292d0,&
                              -7.3234239400d0,-7.3358294097d0,-7.3419856567d0,-7.3436926629d0,&
                              -7.3422028452d0,-7.3384051700d0,-7.3188354192d0,-7.2767603816d0,&
                              -7.2364544565d0,-7.2025462214d0,-7.1755472270d0,-7.1511905091d0,&
                              -7.1346918995d0,-7.2679949325d0,-7.3233802083d0,-7.3357853106d0,&
                              -7.3419415575d0,-7.3436481962d0,-7.3421580111d0,-7.3383599683d0,&
                              -7.3187887476d0,-7.2767093000d0,-7.2363941876d0,-7.2024701503d0,&
                              -7.1754480038d0,-7.1510541691d0,-7.1346503728d0,-7.2679526707d0,&
                              -7.3233364766d0,-7.3357412114d0,-7.3418967234d0,-7.3436029945d0,&
                              -7.3421124419d0,-7.3383140316d0,-7.3187413409d0,-7.2766582184d0,&
                              -7.2363364911d0,-7.2024006941d0,-7.1753605404d0,-7.1509373063d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(2)%zrCut(3)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Fourth cut2d
   this%wyckoffSite(2)%zrcut(4)%alias='Top_F_90_0'
   this%wyckoffSite(2)%zrCut(4)%theta=1.57079632679d0
   this%wyckoffSite(2)%zrCut(4)%phi=0.d0
   gridPot1414(:,:)=reshape( [ 0.9741794146d0,-0.9041065359d0,-2.4291775949d0,-3.0667038997d0,&
                              -3.6277203507d0,-4.1181084570d0,-4.5440853224d0,-4.9120283881d0,&
                              -5.7275583925d0,-6.4391914392d0,-6.7431337173d0,-6.8884126204d0,&
                              -6.9689340673d0,-7.0192523461d0,-5.1308742958d0,-5.5135512685d0,&
                              -5.8165858383d0,-5.9475990208d0,-6.0670850424d0,-6.1758244588d0,&
                              -6.2743773297d0,-6.3632111065d0,-6.5765523327d0,-6.7934101019d0,&
                              -6.9072381973d0,-6.9740712555d0,-7.0171866665d0,-7.0473416930d0,&
                              -6.7654780421d0,-6.9156820898d0,-6.9907973434d0,-7.0138200608d0,&
                              -7.0309430415d0,-7.0438607969d0,-7.0537309307d0,-7.0613593557d0,&
                              -7.0757500240d0,-7.0863525719d0,-7.0913372504d0,-7.0948104291d0,&
                              -7.0976735690d0,-7.1001504736d0,-7.0621351339d0,-7.1942206617d0,&
                              -7.2490638850d0,-7.2614899344d0,-7.2678648399d0,-7.2699841735d0,&
                              -7.2690941048d0,-7.2660762502d0,-7.2498238610d0,-7.2158435973d0,&
                              -7.1855169515d0,-7.1621311508d0,-7.1452110264d0,-7.1316078961d0,&
                              -7.1093502997d0,-7.2414038556d0,-7.2957072313d0,-7.3076665642d0,&
                              -7.3134498056d0,-7.3148595097d0,-7.3131506661d0,-7.3092151808d0,&
                              -7.2897101088d0,-7.2492542465d0,-7.2121241981d0,-7.1824083261d0,&
                              -7.1600621638d0,-7.1413644745d0,-7.1261903106d0,-7.2587840816d0,&
                              -7.3134806750d0,-7.3255671606d0,-7.3314275755d0,-7.3328633717d0,&
                              -7.3311299060d0,-7.3271194521d0,-7.3071039321d0,-7.2649770779d0,&
                              -7.2254700831d0,-7.1930773902d0,-7.1680977713d0,-7.1466108082d0,&
                              -7.1329995931d0,-7.2660317835d0,-7.3211340895d0,-7.3333991769d0,&
                              -7.3394198188d0,-7.3409956299d0,-7.3393801296d0,-7.3354644889d0,&
                              -7.3155897188d0,-7.2732162766d0,-7.2329114540d0,-7.1993148532d0,&
                              -7.1729170778d0,-7.1497058364d0,-7.1345768741d0,-7.2677766415d0,&
                              -7.3230516693d0,-7.3354009126d0,-7.3415009331d0,-7.3431520803d0,&
                              -7.3416067712d0,-7.3377554419d0,-7.3180350189d0,-7.2757659448d0,&
                              -7.2353582241d0,-7.2014639538d0,-7.1746167341d0,-7.1507576021d0,&
                              -7.1348091298d0,-7.2680673286d0,-7.3233982154d0,-7.3357735508d0,&
                              -7.3419003983d0,-7.3435776375d0,-7.3420569504d0,-7.3382295082d0,&
                              -7.3185711915d0,-7.2763664288d0,-7.2359642205d0,-7.2020126212d0,&
                              -7.1750529486d0,-7.1510089675d0,-7.1348161122d0,-7.2680926857d0,&
                              -7.3234426821d0,-7.3358286748d0,-7.3419650771d0,-7.3436518711d0,&
                              -7.3421407389d0,-7.3383224840d0,-7.3186887894d0,-7.2765130586d0,&
                              -7.2361204051d0,-7.2021570460d0,-7.1751646665d0,-7.1510622540d0,&
                              -7.1347620907d0,-7.2680577738d0,-7.3234195300d0,-7.3358114026d0,&
                              -7.3419544198d0,-7.3436470937d0,-7.3421414739d0,-7.3383290989d0,&
                              -7.3187119415d0,-7.2765564228d0,-7.2361729566d0,-7.2022062901d0,&
                              -7.1751981084d0,-7.1510644589d0,-7.1347054967d0,-7.2679938300d0,&
                              -7.3233622011d0,-7.3357573811d0,-7.3419037057d0,-7.3435996871d0,&
                              -7.3420981097d0,-7.3382894096d0,-7.3186825420d0,-7.2765409881d0,&
                              -7.2361656068d0,-7.2021989403d0,-7.1751841437d0,-7.1510343245d0,&
                              -7.1346720549d0,-7.2679574482d0,-7.3233272892d0,-7.3357235717d0,&
                              -7.3418709988d0,-7.3435677152d0,-7.3420672402d0,-7.3382596426d0,&
                              -7.3186560825d0,-7.2765200410d0,-7.2361475996d0,-7.2021820356d0,&
                              -7.1751650340d0,-7.1510108049d0,-7.1346327331d0,-7.2679192289d0,&
                              -7.3232898049d0,-7.3356868224d0,-7.3418346170d0,-7.3435324358d0,&
                              -7.3420326959d0,-7.3382258333d0,-7.3186248456d0,-7.2764924790d0,&
                              -7.2361229775d0,-7.2021581485d0,-7.1751404120d0,-7.1509828754d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(2)%zrCut(4)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Fifth cut2d
   this%wyckoffSite(2)%zrcut(5)%alias='Top_F_90_45'
   this%wyckoffSite(2)%zrCut(5)%theta=1.57079632679d0
   this%wyckoffSite(2)%zrCut(5)%phi=0.785398163397d0
   gridPot1414(:,:)=reshape( [ 0.9743341660d0,-0.9037241960d0,-2.4283785544d0,-3.0655893662d0,&
                              -3.6261854049d0,-4.1160920215d0,-4.5414735478d0,-4.9086048210d0,&
                              -5.7208854500d0,-6.4226847447d0,-6.7098175139d0,-6.8282418475d0,&
                              -6.8621427326d0,-6.8124943940d0,-5.1307636803d0,-5.5132782211d0,&
                              -5.8160221036d0,-5.9468214051d0,-6.0660284993d0,-6.1744349668d0,&
                              -6.2725773477d0,-6.3609102312d0,-6.5721475585d0,-6.7828579006d0,&
                              -6.8859570305d0,-6.9347939440d0,-6.9476834024d0,-6.9162590542d0,&
                              -6.7654368829d0,-6.9156012413d0,-6.9906360139d0,-7.0136025048d0,&
                              -7.0306549268d0,-7.0434877912d0,-7.0532557619d0,-7.0607603416d0,&
                              -7.0746335795d0,-7.0836871433d0,-7.0857983920d0,-7.0843280515d0,&
                              -7.0792331250d0,-7.0668537473d0,-7.0621171267d0,-7.1941787675d0,&
                              -7.2489830365d0,-7.2613818914d0,-7.2677233550d0,-7.2698044693d0,&
                              -7.2688695664d0,-7.2657998953d0,-7.2493420774d0,-7.2148443831d0,&
                              -7.1837232169d0,-7.1591739326d0,-7.1405813464d0,-7.1241062563d0,&
                              -7.1093403774d0,-7.2413807035d0,-7.2956623971d0,-7.3076062953d0,&
                              -7.3133700595d0,-7.3147577141d0,-7.3130235134d0,-7.3090564237d0,&
                              -7.2894308140d0,-7.2486750771d0,-7.2111055068d0,-7.1807983381d0,&
                              -7.1576752951d0,-7.1376660224d0,-7.1261859006d0,-7.2587734242d0,&
                              -7.3134589929d0,-7.3255377612d0,-7.3313875188d0,-7.3328122901d0,&
                              -7.3310663297d0,-7.3270397061d0,-7.3069613447d0,-7.2646738959d0,&
                              -7.2249269281d0,-7.1922079012d0,-7.1667964777d0,-7.1445263865d0,&
                              -7.1329981231d0,-7.2660284761d0,-7.3211282096d0,-7.3333910920d0,&
                              -7.3394080591d0,-7.3409805627d0,-7.3393617549d0,-7.3354409694d0,&
                              -7.3155463546d0,-7.2731166860d0,-7.2327218275d0,-7.1989852118d0,&
                              -7.1723713503d0,-7.1486489258d0,-7.1345765066d0,-7.2677773765d0,&
                              -7.3230538743d0,-7.3354020151d0,-7.3415016681d0,-7.3431531828d0,&
                              -7.3416086087d0,-7.3377565444d0,-7.3180357539d0,-7.2757571250d0,&
                              -7.2353266197d0,-7.2013786953d0,-7.1744142453d0,-7.1501699804d0,&
                              -7.1348157447d0,-7.2680691661d0,-7.3234007879d0,-7.3357779607d0,&
                              -7.3419055432d0,-7.3435842524d0,-7.3420657703d0,-7.3382394305d0,&
                              -7.3185858913d0,-7.2763888459d0,-7.2359862700d0,-7.2020104162d0,&
                              -7.1749673227d0,-7.1505808378d0,-7.1348238296d0,-7.2680945231d0,&
                              -7.3234470920d0,-7.3358345546d0,-7.3419720595d0,-7.3436606910d0,&
                              -7.3421521312d0,-7.3383353462d0,-7.3187093690d0,-7.2765468680d0,&
                              -7.2361619318d0,-7.2021842405d0,-7.1751205674d0,-7.1506896158d0,&
                              -7.1347800979d0,-7.2680599788d0,-7.3234246750d0,-7.3358180175d0,&
                              -7.3419621372d0,-7.3436570161d0,-7.3421547036d0,-7.3383441661d0,&
                              -7.3187358285d0,-7.2765975821d0,-7.2362266106d0,-7.2022511243d0,&
                              -7.1751782638d0,-7.1507241602d0,-7.1347143166d0,-7.2679964024d0,&
                              -7.3233673460d0,-7.3357647309d0,-7.3419121581d0,-7.3436107119d0,&
                              -7.3421120744d0,-7.3383055793d0,-7.3187086340d0,-7.2765865573d0,&
                              -7.2362262431d0,-7.2022547993d0,-7.1751782638d0,-7.1507116654d0,&
                              -7.1346713199d0,-7.2679600206d0,-7.3233328016d0,-7.3357309216d0,&
                              -7.3418794512d0,-7.3435791075d0,-7.3420819400d0,-7.3382761798d0,&
                              -7.3186829095d0,-7.2765670801d0,-7.2362108084d0,-7.2022412020d0,&
                              -7.1751642990d0,-7.1506940257d0,-7.1346331006d0,-7.2679221688d0,&
                              -7.3232956848d0,-7.3356945397d0,-7.3418438043d0,-7.3435441956d0,&
                              -7.3420477631d0,-7.3382431054d0,-7.3186524076d0,-7.2765406206d0,&
                              -7.2361880238d0,-7.2022202549d0,-7.1751429844d0,-7.1506697712d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(2)%zrCut(5)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Wickoff Hollow -----------> 'c'
   this%wyckoffSite(3)%id='c'
   this%wyckoffSite(3)%myNumber=3
   this%wyckoffSite(3)%is_homonucl=.true.
   this%wyckoffSite(3)%x=2.7216780628885480d0
   this%wyckoffSite(3)%y=0.d0
   this%wyckoffSite(3)%n2dcuts=4
   this%wyckoffSite(3)%nphicuts=2
   allocate( this%wyckoffSite(3)%nPhiPoints(2) )
   this%wyckoffSite(3)%nPhiPoints(:)=[1,3]
   allocate( this%wyckoffSite(3)%zrCut(4) )
   this%wyckoffSite(3)%zrCut(:)%x=2.7216780628885480d0
   this%wyckoffSite(3)%zrCut(:)%y=0.d0
   ! Reading zrcuts
   this%wyckoffSite(3)%zrcut(1)%alias='Hollow_0_0'
   this%wyckoffSite(3)%zrCut(1)%theta=0.d0
   this%wyckoffSite(3)%zrCut(1)%phi=0.d0
   gridPot1414(:,:)=reshape( [ -6.9155832341d0,-7.0431757894d0,-7.0949357443d0,-7.1065338314d0,&
                               -7.1126007775d0,-7.1149446495d0,-7.1148035321d0,-7.1130384620d0,&
                               -7.1032602015d0,-7.0866983830d0,-7.0775919002d0,-7.0752053990d0,&
                               -7.0765959935d0,-7.0795432894d0,-6.9649287583d0,-7.0927587143d0,&
                               -7.1442276144d0,-7.1554471835d0,-7.1609845718d0,-7.1626478463d0,&
                               -7.1616850140d0,-7.1589677689d0,-7.1456549582d0,-7.1216381716d0,&
                               -7.1043222570d0,-7.0941989203d0,-7.0893777763d0,-7.0877354490d0,&
                               -7.0637296871d0,-7.1931130371d0,-7.2449909573d0,-7.2559113868d0,&
                               -7.2608001495d0,-7.2614708247d0,-7.2591846492d0,-7.2548324266d0,&
                               -7.2349852184d0,-7.1962117402d0,-7.1619216796d0,-7.1347734830d0,&
                               -7.1143239534d0,-7.0973924367d0,-7.1139211808d0,-7.2455660842d0,&
                               -7.2993351247d0,-7.3109945831d0,-7.3164617803d0,-7.3175466204d0,&
                               -7.3155074003d0,-7.3112382312d0,-7.2907258602d0,-7.2484402489d0,&
                               -7.2086991610d0,-7.1749118312d0,-7.1468114595d0,-7.1193957278d0,&
                               -7.1260984372d0,-7.2585944550d0,-7.3131826380d0,-7.3252162046d0,&
                               -7.3310288454d0,-7.3324242173d0,-7.3306598821d0,-7.3266292161d0,&
                               -7.3066232509d0,-7.2646639736d0,-7.2250246813d0,-7.1914078684d0,&
                               -7.1634706637d0,-7.1357473401d0,-7.1316064261d0,-7.2645588705d0,&
                               -7.3195983352d0,-7.3318457829d0,-7.3378583400d0,-7.3394396635d0,&
                               -7.3378440078d0,-7.3339640140d0,-7.3143060648d0,-7.2726654042d0,&
                               -7.2333421561d0,-7.2004088806d0,-7.1736645591d0,-7.1478499955d0,&
                               -7.1342582575d0,-7.2674587598d0,-7.3227451799d0,-7.3351083880d0,&
                               -7.3412289881d0,-7.3429087998d0,-7.3414009750d0,-7.3375966848d0,&
                               -7.3180864679d0,-7.2764487473d0,-7.2370380357d0,-7.2043381185d0,&
                               -7.1784996678d0,-7.1548683816d0,-7.1349043106d0,-7.2681768416d0,&
                               -7.3235360254d0,-7.3359301029d0,-7.3420778975d0,-7.3437812288d0,&
                               -7.3422895737d0,-7.3384952058d0,-7.3189710242d0,-7.2771392671d0,&
                               -7.2373522425d0,-7.2042657223d0,-7.1783008540d0,-7.1551998605d0,&
                               -7.1349465723d0,-7.2682411530d0,-7.3236165064d0,-7.3360186688d0,&
                               -7.3421719758d0,-7.3438771445d0,-7.3423869594d0,-7.3385903866d0,&
                               -7.3190386430d0,-7.2770812032d0,-7.2370453856d0,-7.2036148918d0,&
                               -7.1772990674d0,-7.1539459736d0,-7.1348936533d0,-7.2681988912d0,&
                               -7.3235804921d0,-7.3359852269d0,-7.3421392689d0,-7.3438448051d0,&
                               -7.3423538850d0,-7.3385562097d0,-7.3189916038d0,-7.2769742626d0,&
                               -7.2368226847d0,-7.2032099142d0,-7.1766677140d0,-7.1530525975d0,&
                               -7.1348230946d0,-7.2681268625d0,-7.3235117709d0,-7.3359172406d0,&
                               -7.3420727526d0,-7.3437782888d0,-7.3422866337d0,-7.3384882234d0,&
                               -7.3189151652d0,-7.2768588698d0,-7.2366253408d0,-7.2028777003d0,&
                               -7.1761510184d0,-7.1522830666d0,-7.1347323238d0,-7.2680423391d0,&
                               -7.3234283499d0,-7.3358341872d0,-7.3419896992d0,-7.3436952354d0,&
                               -7.3422032127d0,-7.3384040675d0,-7.3188258644d0,-7.2767453143d0,&
                               -7.2364584989d0,-7.2026171476d0,-7.1757548607d0,-7.1516767027d0,&
                               -7.1346893270d0,-7.2679989749d0,-7.3233846182d0,-7.3357900880d0,&
                               -7.3419456000d0,-7.3436507687d0,-7.3421587461d0,-7.3383592333d0,&
                               -7.3187795602d0,-7.2766901904d0,-7.2363850002d0,-7.2025102071d0,&
                               -7.1755961036d0,-7.1514326872d0,-7.1346470653d0,-7.2679563457d0,&
                               -7.3233405190d0,-7.3357463563d0,-7.3419007658d0,-7.3436059345d0,&
                               -7.3421139119d0,-7.3383132966d0,-7.3187321536d0,-7.2766369038d0,&
                               -7.2363170140d0,-7.2024153938d0,-7.1754593961d0,-7.1512235835d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(3)%zrCut(1)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! 2nd cut2d
   this%wyckoffSite(3)%zrcut(2)%alias='Hollow_90_0'
   this%wyckoffSite(3)%zrCut(2)%theta=1.57079632679d0
   this%wyckoffSite(3)%zrCut(2)%phi=0.d0
   gridPot1414(:,:)=reshape( [ -6.8925656617d0,-7.0062897566d0,-7.0400972985d0,-7.0410130916d0,&
                               -7.0351728889d0,-7.0243105233d0,-7.0095927859d0,-6.9918072149d0,&
                               -6.9244964155d0,-6.7748156782d0,-6.5731195782d0,-6.3010488876d0,&
                               -5.9367693621d0,-5.3591787448d0,-6.9498236831d0,-7.0685791281d0,&
                               -7.1083084563d0,-7.1125651307d0,-7.1103590687d0,-7.1034671002d0,&
                               -7.0931085678d0,-7.0801298086d0,-7.0305788557d0,-6.9265415155d0,&
                               -6.8019491751d0,-6.6550058950d0,-6.4870559656d0,-6.2760931557d0,&
                               -7.0595226244d0,-7.1864867662d0,-7.2353225772d0,-7.2444584595d0,&
                               -7.2473785609d0,-7.2458938882d0,-7.2412638407d0,-7.2343828970d0,&
                               -7.2058959224d0,-7.1498528337d0,-7.0959875100d0,-7.0482593236d0,&
                               -7.0073962787d0,-6.9684544886d0,-7.1126202547d0,-7.2435981578d0,&
                               -7.2965730454d0,-7.3077830596d0,-7.3127673706d0,-7.3133333102d0,&
                               -7.3107413803d0,-7.3058882644d0,-7.2834737483d0,-7.2379409667d0,&
                               -7.1959019434d0,-7.1618176791d0,-7.1358506057d0,-7.1139513152d0,&
                               -7.1253333163d0,-7.2574423637d0,-7.3115737525d0,-7.3233526463d0,&
                               -7.3288922396d0,-7.3299961893d0,-7.3279216899d0,-7.3235632199d0,&
                               -7.3024830719d0,-7.2586312044d0,-7.2176347594d0,-7.1841175372d0,&
                               -7.1584154265d0,-7.1365396556d0,-7.1311393422d0,-7.2638514460d0,&
                               -7.3186057359d0,-7.3306929565d0,-7.3365338943d0,-7.3379300012d0,&
                               -7.3361366341d0,-7.3320464342d0,-7.3116770181d0,-7.2686351057d0,&
                               -7.2279385353d0,-7.1942588810d0,-7.1680599195d0,-7.1453458964d0,&
                               -7.1340234293d0,-7.2671008214d0,-7.3222380392d0,-7.3345167238d0,&
                               -7.3405465532d0,-7.3421275091d0,-7.3405127438d0,-7.3365937957d0,&
                               -7.3166833787d0,-7.2741684516d0,-7.2336368857d0,-7.1997595201d0,&
                               -7.1730611352d0,-7.1495048176d0,-7.1347863453d0,-7.2680044873d0,&
                               -7.3232901724d0,-7.3356438256d0,-7.3417478886d0,-7.3434019757d0,&
                               -7.3418585041d0,-7.3380082772d0,-7.3182860168d0,-7.2759959956d0,&
                               -7.2355493206d0,-7.2015991913d0,-7.1746839854d0,-7.1507355525d0,&
                               -7.1348697662d0,-7.2681397248d0,-7.3234750215d0,-7.3358521943d0,&
                               -7.3419805118d0,-7.3436588535d0,-7.3421389014d0,-7.3383118267d0,&
                               -7.3186524076d0,-7.2764406624d0,-7.2360285318d0,-7.2020626003d0,&
                               -7.1750830830d0,-7.1510100699d0,-7.1348532290d0,-7.2681268625d0,&
                               -7.3234794314d0,-7.3358668941d0,-7.3420040314d0,-7.3436911929d0,&
                               -7.3421804282d0,-7.3383618058d0,-7.3187262737d0,-7.2765453980d0,&
                               -7.2361490696d0,-7.2021809331d0,-7.1751815712d0,-7.1510666639d0,&
                               -7.1347782604d0,-7.2680743110d0,-7.3234382722d0,-7.3358312472d0,&
                               -7.3419742644d0,-7.3436669384d0,-7.3421616860d0,-7.3383489435d0,&
                               -7.3187292136d0,-7.2765692851d0,-7.2361843489d0,-7.2022165800d0,&
                               -7.1752072958d0,-7.1510696039d0,-7.1347161541d0,-7.2680019148d0,&
                               -7.3233717559d0,-7.3357676709d0,-7.3419143630d0,-7.3436103444d0,&
                               -7.3421087670d0,-7.3382993319d0,-7.3186898919d0,-7.2765439280d0,&
                               -7.2361674442d0,-7.2022022477d0,-7.1751881861d0,-7.1510383669d0,&
                               -7.1346702174d0,-7.2679633281d0,-7.3233346391d0,-7.3357316565d0,&
                               -7.3418790837d0,-7.3435761675d0,-7.3420753251d0,-7.3382669925d0,&
                               -7.3186608599d0,-7.2765200410d0,-7.2361464971d0,-7.2021820356d0,&
                               -7.1751676065d0,-7.1510141124d0,-7.1346360405d0,-7.2679236388d0,&
                               -7.3232956848d0,-7.3356930697d0,-7.3418412319d0,-7.3435390507d0,&
                               -7.3420393107d0,-7.3382317131d0,-7.3186277855d0,-7.2764906415d0,&
                               -7.2361196701d0,-7.2021563111d0,-7.1751411470d0,-7.1509854479d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(3)%zrCut(2)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! third cut2d
   this%wyckoffSite(3)%zrcut(3)%alias='Hollow_90_45'
   this%wyckoffSite(3)%zrCut(3)%theta=1.57079632679d0
   this%wyckoffSite(3)%zrCut(3)%phi=0.785398163397d0
   gridPot1414(:,:)=reshape( [-6.9006196438d0,-7.0199344136d0,-7.0617540434d0,-7.0677937950d0,&
                              -7.0679308700d0,-7.0639924448d0,-7.0572294665d0,-7.0485217138d0,&
                              -7.0162065620d0,-6.9582091442d0,-6.9074149616d0,-6.8693941095d0,&
                              -6.8470615445d0,-6.8419978549d0,-6.9566623651d0,-7.0798291991d0,&
                              -7.1255971764d0,-7.1335912572d0,-7.1356613467d0,-7.1336261691d0,&
                              -7.1287407137d0,-7.1218876995d0,-7.0950257801d0,-7.0458228434d0,&
                              -7.0029690875d0,-6.9706870101d0,-6.9502775373d0,-6.9413033520d0,&
                              -7.0635521879d0,-7.1927352540d0,-7.2443246920d0,-7.2550529225d0,&
                              -7.2597149420d0,-7.2601269019d0,-7.2575489367d0,-7.2528725851d0,&
                              -7.2318997451d0,-7.1910605872d0,-7.1550694004d0,-7.1271204360d0,&
                              -7.1068954448d0,-7.0912597093d0,-7.1144422862d0,-7.2463558272d0,&
                              -7.3004295196d0,-7.3122484701d0,-7.3178788343d0,-7.3191297813d0,&
                              -7.3172588731d0,-7.3131587509d0,-7.2931858600d0,-7.2521210614d0,&
                              -7.2145650883d0,-7.1844163092d0,-7.1614744404d0,-7.1417330702d0,&
                              -7.1263681773d0,-7.2589986976d0,-7.3137335104d0,-7.3258413106d0,&
                              -7.3317259801d0,-7.3331922782d0,-7.3314955618d0,-7.3275270021d0,&
                              -7.3076871438d0,-7.2660402359d0,-7.2272009763d0,-7.1954870435d0,&
                              -7.1709803884d0,-7.1495140049d0,-7.1317052818d0,-7.2646988855d0,&
                              -7.3197762019d0,-7.3320379818d0,-7.3380600938d0,-7.3396458272d0,&
                              -7.3380483340d0,-7.3341576829d0,-7.3144141079d0,-7.2724522581d0,&
                              -7.2327949586d0,-7.1999818535d0,-7.1743253120d0,-7.1516046740d0,&
                              -7.1342685473d0,-7.2674672122d0,-7.3227426075d0,-7.3350955257d0,&
                              -7.3412014261d0,-7.3428617607d0,-7.3413289463d0,-7.3374919492d0,&
                              -7.3178351026d0,-7.2757380153d0,-7.2355901123d0,-7.2020188686d0,&
                              -7.1754821807d0,-7.1517399115d0,-7.1348962257d0,-7.2681614069d0,&
                              -7.3235069934d0,-7.3358926186d0,-7.3420290209d0,-7.3437172849d0,&
                              -7.3422094601d0,-7.3383941452d0,-7.3187791927d0,-7.2766615259d0,&
                              -7.2363614807d0,-7.2025102071d0,-7.1756056584d0,-7.1513841781d0,&
                              -7.1349359150d0,-7.2682235133d0,-7.3235907819d0,-7.3359863294d0,&
                              -7.3421322865d0,-7.3438297379d0,-7.3423296304d0,-7.3385216653d0,&
                              -7.3189225151d0,-7.2768074207d0,-7.2364739336d0,-7.2025506313d0,&
                              -7.1755409796d0,-7.1511677246d0,-7.1348829960d0,-7.2681827215d0,&
                              -7.3235577075d0,-7.3359569299d0,-7.3421061945d0,-7.3438069533d0,&
                              -7.3423101533d0,-7.3385047606d0,-7.3189122253d0,-7.2768004383d0,&
                              -7.2364599689d0,-7.2025164544d0,-7.1754748308d0,-7.1510523317d0,&
                              -7.1348164797d0,-7.2681117953d0,-7.3234911912d0,-7.3358926186d0,&
                              -7.3420440882d0,-7.3437466844d0,-7.3422517218d0,-7.3384485342d0,&
                              -7.3188607762d0,-7.2767530317d0,-7.2364103573d0,-7.2024572880d0,&
                              -7.1753991272d0,-7.1509501685d0,-7.1347304863d0,-7.2680272719d0,&
                              -7.3234085053d0,-7.3358106676d0,-7.3419636071d0,-7.3436669384d0,&
                              -7.3421730783d0,-7.3383709931d0,-7.3187861751d0,-7.2766817380d0,&
                              -7.2363394311d0,-7.2023830544d0,-7.1753171762d0,-7.1508542528d0,&
                              -7.1346882246d0,-7.2679842752d0,-7.3233655085d0,-7.3357680384d0,&
                              -7.3419206104d0,-7.3436243092d0,-7.3421308166d0,-7.3383287314d0,&
                              -7.3187446484d0,-7.2766416812d0,-7.2362990068d0,-7.2023422626d0,&
                              -7.1752741795d0,-7.1508075812d0,-7.1346444929d0,-7.2679412785d0,&
                              -7.3233217768d0,-7.3357243067d0,-7.3418772462d0,-7.3435805775d0,&
                              -7.3420870849d0,-7.3382853672d0,-7.3187016517d0,-7.2765990520d0,&
                              -7.2362567451d0,-7.2022996334d0,-7.1752300803d0,-7.1507605420d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(3)%zrCut(3)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! third cut2d
   this%wyckoffSite(3)%zrcut(4)%alias='Hollow_90_90'
   this%wyckoffSite(3)%zrCut(4)%theta=1.57079632679d0
   this%wyckoffSite(3)%zrCut(4)%phi=1.57079632679d0
   gridPot1414(:,:)=reshape( [-6.9064025177d0,-7.0281412729d0,-7.0723304993d0,-7.0794484761d0,&
                              -7.0805373586d0,-7.0773614819d0,-7.0711074817d0,-7.0625684085d0,&
                              -7.0274220886d0,-6.9381322527d0,-6.7570407645d0,-6.3143201714d0,&
                              -5.2611778488d0,-2.3857877032d0,-6.9621012652d0,-7.0876817950d0,&
                              -7.1359285143d0,-7.1451173156d0,-7.1483089945d0,-7.1472829534d0,&
                              -7.1432526548d0,-7.1370596585d0,-7.1104990836d0,-7.0518978743d0,&
                              -6.9716663796d0,-6.8353256474d0,-6.5833230284d0,-6.0721726190d0,&
                              -7.0672947392d0,-7.1982932220d0,-7.2518888056d0,-7.2636651270d0,&
                              -7.2693903044d0,-7.2708613799d0,-7.2693241556d0,-7.2656547354d0,&
                              -7.2473230694d0,-7.2088042640d0,-7.1716738481d0,-7.1390496345d0,&
                              -7.1113230035d0,-7.0848671641d0,-7.1162290384d0,-7.2490282381d0,&
                              -7.3041066571d0,-7.3164647202d0,-7.3226551441d0,-7.3244797481d0,&
                              -7.3231931542d0,-7.3196843286d0,-7.3014621756d0,-7.2628959636d0,&
                              -7.2267963662d0,-7.1967898071d0,-7.1728149147d0,-7.1509769956d0,&
                              -7.1273927485d0,-7.2605318795d0,-7.3158432891d0,-7.3282616212d0,&
                              -7.3344693172d0,-7.3362681967d0,-7.3349114116d0,-7.3312908681d0,&
                              -7.3125064504d0,-7.2725187744d0,-7.2349481016d0,-7.2039033740d0,&
                              -7.1794658076d0,-7.1576418532d0,-7.1322664440d0,-7.2655375051d0,&
                              -7.3209279258d0,-7.3333572826d0,-7.3395539539d0,-7.3413186565d0,&
                              -7.3399041749d0,-7.3362013129d0,-7.3170288224d0,-7.2759996705d0,&
                              -7.2371486512d0,-7.2049139804d0,-7.1795793630d0,-7.1571009031d0,&
                              -7.1345107253d0,-7.2678277230d0,-7.3232346809d0,-7.3356577904d0,&
                              -7.3418364545d0,-7.3435702876d0,-7.3421120744d0,-7.3383515159d0,&
                              -7.3189225151d0,-7.2771925536d0,-7.2373739246d0,-7.2040823432d0,&
                              -7.1777859959d0,-7.1544273897d0,-7.1349884666d0,-7.2683124466d0,&
                              -7.3237127897d0,-7.3361252418d0,-7.3422906761d0,-7.3440079721d0,&
                              -7.3425284443d0,-7.3387425288d0,-7.3192117323d0,-7.2772204831d0,&
                              -7.2370262760d0,-7.2032712856d0,-7.1764880098d0,-7.1526086656d0,&
                              -7.1349899365d0,-7.2683014218d0,-7.3236958850d0,-7.3361039272d0,&
                              -7.3422642166d0,-7.3439748977d0,-7.3424880200d0,-7.3386940197d0,&
                              -7.3191316188d0,-7.2770665034d0,-7.2367697657d0,-7.2028839477d0,&
                              -7.1759503671d0,-7.1518916863d0,-7.1349252577d0,-7.2682323331d0,&
                              -7.3236231213d0,-7.3360304286d0,-7.3421885130d0,-7.3438966217d0,&
                              -7.3424068040d0,-7.3386094962d0,-7.3190364380d0,-7.2769478031d0,&
                              -7.2366205634d0,-7.2026961587d0,-7.1757170089d0,-7.1516013666d0,&
                              -7.1348396318d0,-7.2681430323d0,-7.3235316155d0,-7.3359374528d0,&
                              -7.3420940672d0,-7.3438003384d0,-7.3423086833d0,-7.3385095380d0,&
                              -7.3189302325d0,-7.2768287353d0,-7.2364864284d0,-7.2025429140d0,&
                              -7.1755413471d0,-7.1513963054d0,-7.1347455535d0,-7.2680467490d0,&
                              -7.3234327598d0,-7.3358374946d0,-7.3419926391d0,-7.3436974403d0,&
                              -7.3422046827d0,-7.3384040675d0,-7.3188207195d0,-7.2767122399d0,&
                              -7.2363625831d0,-7.2024102489d0,-7.1753987597d0,-7.1512408557d0,&
                              -7.1346959419d0,-7.2679993424d0,-7.3233838832d0,-7.3357878830d0,&
                              -7.3419422925d0,-7.3436467262d0,-7.3421532337d0,-7.3383518834d0,&
                              -7.3187666980d0,-7.2766556460d0,-7.2363037842d0,-7.2023492450d0,&
                              -7.1753355509d0,-7.1511743394d0,-7.1346533127d0,-7.2679526707d0,&
                              -7.3233353741d0,-7.3357386389d0,-7.3418923134d0,-7.3435960122d0,&
                              -7.3421021521d0,-7.3383004344d0,-7.3187134114d0,-7.2766005220d0,&
                              -7.2362460878d0,-7.2022908136d0,-7.1752752820d0,-7.1511126006d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(3)%zrCut(4)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Wickoff Bridge -----------> 'f'
   this%wyckoffSite(4)%id='f'
   this%wyckoffSite(4)%myNumber=4
   this%wyckoffSite(4)%is_homonucl=.true.
   this%wyckoffSite(4)%x=1.36083903144d0
   this%wyckoffSite(4)%y=1.36083903144d0
   this%wyckoffSite(4)%n2dcuts=7
   this%wyckoffSite(4)%nphicuts=3
   allocate( this%wyckoffSite(4)%nPhiPoints(3) )
   this%wyckoffSite(4)%nPhiPoints(:)=[1,3,3]
   allocate( this%wyckoffSite(4)%zrCut(7) )
   this%wyckoffSite(4)%zrCut(:)%x=1.36083903144d0
   this%wyckoffSite(4)%zrCut(:)%y=1.36083903144d0
   ! First zrcuts
   this%wyckoffSite(4)%zrcut(1)%alias='Bridge_0_0'
   this%wyckoffSite(4)%zrCut(1)%theta=0.d0
   this%wyckoffSite(4)%zrCut(1)%phi=0.d0
   gridPot1414(:,:)=reshape( [-6.6746098227d0,-6.8113661897d0,-6.8747940554d0,-6.8930176784d0,&
                              -6.9061375550d0,-6.9158812712d0,-6.9234034906d0,-6.9294829315d0,&
                              -6.9436067997d0,-6.9648464398d0,-6.9878140332d0,-7.0109010619d0,&
                              -7.0315031012d0,-7.0497109220d0,-6.8048167250d0,-6.9354557993d0,&
                              -6.9912331904d0,-7.0050615941d0,-7.0134922568d0,-7.0183001711d0,&
                              -7.0206954921d0,-7.0215127971d0,-7.0195106939d0,-7.0150743153d0,&
                              -7.0165020266d0,-7.0234755786d0,-7.0339050371d0,-7.0478271516d0,&
                              -7.0202012137d0,-7.1482479906d0,-7.1989161230d0,-7.2092731854d0,&
                              -7.2136151182d0,-7.2137437409d0,-7.2109063254d0,-7.2059793434d0,&
                              -7.1841307670d0,-7.1406596224d0,-7.1002507992d0,-7.0671889011d0,&
                              -7.0430262197d0,-7.0263350436d0,-7.1061104792d0,-7.2372419945d0,&
                              -7.2905292513d0,-7.3019446943d0,-7.3071700808d0,-7.3080101704d0,&
                              -7.3057206874d0,-7.3011883931d0,-7.2797308295d0,-7.2348826878d0,&
                              -7.1904781104d0,-7.1494077994d0,-7.1114538311d0,-7.0698249302d0,&
                              -7.1236255751d0,-7.2559290265d0,-7.3103279504d0,-7.3222714811d0,&
                              -7.3279973935d0,-7.3293089769d0,-7.3274630583d0,-7.3233511763d0,&
                              -7.3030813509d0,-7.2603488678d0,-7.2190055093d0,-7.1820992642d0,&
                              -7.1488407573d0,-7.1112994839d0,-7.1308677647d0,-7.2637540603d0,&
                              -7.3187314186d0,-7.3309502018d0,-7.3369370344d0,-7.3384952058d0,&
                              -7.3368797055d0,-7.3329828070d0,-7.3132870061d0,-7.2715566771d0,&
                              -7.2318762255d0,-7.1979650505d0,-7.1692509651d0,-7.1391918544d0,&
                              -7.1341200800d0,-7.2673055151d0,-7.3225798080d0,-7.3349378711d0,&
                              -7.3410544288d0,-7.3427313005d0,-7.3412220058d0,-7.3374173481d0,&
                              -7.3179177885d0,-7.2763373968d0,-7.2369954065d0,-7.2042741747d0,&
                              -7.1781832561d0,-7.1537067354d0,-7.1348774836d0,-7.2681456047d0,&
                              -7.3235033185d0,-7.3358966610d0,-7.3420440882d0,-7.3437466844d0,&
                              -7.3422568668d0,-7.3384636014d0,-7.3189467697d0,-7.2771462495d0,&
                              -7.2374176563d0,-7.2044072072d0,-7.1785037102d0,-7.1553869146d0,&
                              -7.1349395899d0,-7.2682290257d0,-7.3236040117d0,-7.3360061740d0,&
                              -7.3421602160d0,-7.3438657522d0,-7.3423763021d0,-7.3385804642d0,&
                              -7.3190327631d0,-7.2770900230d0,-7.2370836049d0,-7.2037008852d0,&
                              -7.1774464322d0,-7.1541605896d0,-7.1348899784d0,-7.2681904389d0,&
                              -7.3235724072d0,-7.3359771420d0,-7.3421319190d0,-7.3438378227d0,&
                              -7.3423472701d0,-7.3385503298d0,-7.3189886639d0,-7.2769801425d0,&
                              -7.2368462043d0,-7.2032624658d0,-7.1767639972d0,-7.1532095171d0,&
                              -7.1348190522d0,-7.2681198802d0,-7.3235055235d0,-7.3359106258d0,&
                              -7.3420665052d0,-7.3437724089d0,-7.3422811213d0,-7.3384834460d0,&
                              -7.3189129603d0,-7.2768621772d0,-7.2366297507d0,-7.2029085698d0,&
                              -7.1762094499d0,-7.1523848622d0,-7.1347271788d0,-7.2680368267d0,&
                              -7.3234224700d0,-7.3358283073d0,-7.3419838193d0,-7.3436897230d0,&
                              -7.3421980678d0,-7.3383992901d0,-7.3188232919d0,-7.2767471518d0,&
                              -7.2364665837d0,-7.2026355223d0,-7.1757894051d0,-7.1517384416d0,&
                              -7.1346896945d0,-7.2679934625d0,-7.3233791058d0,-7.3357845756d0,&
                              -7.3419400876d0,-7.3436456238d0,-7.3421539686d0,-7.3383544559d0,&
                              -7.3187762528d0,-7.2766909253d0,-7.2363912476d0,-7.2025241718d0,&
                              -7.1756218282d0,-7.1514800938d0,-7.1346481678d0,-7.2679515683d0,&
                              -7.3233353741d0,-7.3357412114d0,-7.3418956209d0,-7.3436011571d0,&
                              -7.3421091345d0,-7.3383088867d0,-7.3187288462d0,-7.2766365363d0,&
                              -7.2363217914d0,-7.2024264186d0,-7.1754792408d0,-7.1512599654d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(4)%zrCut(1)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Second cuts2d
   this%wyckoffSite(4)%zrcut(2)%alias='Bridge_90_0'
   this%wyckoffSite(4)%zrCut(2)%theta=1.57079632679d0
   this%wyckoffSite(4)%zrCut(2)%phi=0.d0
   gridPot1414(:,:)=reshape( [-6.6458520056d0,-6.7681765451d0,-6.8148717079d0,-6.8238907273d0,&
                              -6.8272709303d0,-6.8267626871d0,-6.8235563085d0,-6.8184977638d0,&
                              -6.7975223514d0,-6.7619155625d0,-6.7443820919d0,-6.7518649895d0,&
                              -6.7767203957d0,-6.8098524850d0,-6.7978108336d0,-6.9250487579d0,&
                              -6.9767300692d0,-6.9882950819d0,-6.9942833845d0,-6.9964604145d0,&
                              -6.9960293449d0,-6.9938210779d0,-6.9816030297d0,-6.9576579043d0,&
                              -6.9401854375d0,-6.9320693490d0,-6.9309404097d0,-6.9331009026d0,&
                              -7.0248044342d0,-7.1553163559d0,-7.2088935649d0,-7.2208580428d0,&
                              -7.2269036743d0,-7.2288289714d0,-7.2278775314d0,-7.2249239881d0,&
                              -7.2094822891d0,-7.1778910990d0,-7.1495481818d0,-7.1266504121d0,&
                              -7.1087391584d0,-7.0929325386d0,-7.1077520715d0,-7.2397670406d0,&
                              -7.2940575540d0,-7.3060282793d0,-7.3118405526d0,-7.3132976634d0,&
                              -7.3116571735d0,-7.3078109891d0,-7.2886936225d0,-7.2491697230d0,&
                              -7.2129532629d0,-7.1837290968d0,-7.1613101709d0,-7.1419517287d0,&
                              -7.1242609710d0,-7.2568903888d0,-7.3116538660d0,-7.3237888608d0,&
                              -7.3297095446d0,-7.3312188394d0,-7.3295721021d0,-7.3256619739d0,&
                              -7.3060264418d0,-7.2647716491d0,-7.2263006178d0,-7.1948744322d0,&
                              -7.1706051777d0,-7.1495489168d0,-7.1310765009d0,-7.2640598147d0,&
                              -7.3191327212d0,-7.3313948686d0,-7.3374206555d0,-7.3390104313d0,&
                              -7.3374188181d0,-7.3335369868d0,-7.3138250162d0,-7.2719285803d0,&
                              -7.2323374295d0,-7.1995959856d0,-7.1740515295d0,-7.1516961799d0,&
                              -7.1341498469d0,-7.2673441019d0,-7.3226147198d0,-7.3349650656d0,&
                              -7.3410694961d0,-7.3427268906d0,-7.3411911363d0,-7.3373526693d0,&
                              -7.3176866353d0,-7.2755730109d0,-7.2354082032d0,-7.2018362244d0,&
                              -7.1753594380d0,-7.1519361529d0,-7.1348532290d0,-7.2681342124d0,&
                              -7.3234768590d0,-7.3358610142d0,-7.3419959465d0,-7.3436820056d0,&
                              -7.3421712408d0,-7.3383537209d0,-7.3187281112d0,-7.2765869247d0,&
                              -7.2362626250d0,-7.2024036340d0,-7.1755475945d0,-7.1516311335d0,&
                              -7.1349134979d0,-7.2682106510d0,-7.3235753472d0,-7.3359697922d0,&
                              -7.3421146469d0,-7.3438098932d0,-7.3423072133d0,-7.3384970433d0,&
                              -7.3188879707d0,-7.2767504592d0,-7.2363938201d0,-7.2024631679d0,&
                              -7.1755016579d0,-7.1514312172d0,-7.1348587414d0,-7.2681724317d0,&
                              -7.3235448452d0,-7.3359429652d0,-7.3420914948d0,-7.3437896811d0,&
                              -7.3422899411d0,-7.3384830785d0,-7.3188809884d0,-7.2767482543d0,&
                              -7.2363850002d0,-7.2024352384d0,-7.1754421239d0,-7.1513228067d0,&
                              -7.1348072924d0,-7.2681029755d0,-7.3234797989d0,-7.3358793888d0,&
                              -7.3420301234d0,-7.3437301472d0,-7.3422322447d0,-7.3384272196d0,&
                              -7.3188299068d0,-7.2767012152d0,-7.2363379611d0,-7.2023797470d0,&
                              -7.1753715652d0,-7.1512276260d0,-7.1347161541d0,-7.2680191870d0,&
                              -7.3233974805d0,-7.3357981728d0,-7.3419496424d0,-7.3436504012d0,&
                              -7.3421539686d0,-7.3383496785d0,-7.3187556732d0,-7.2766310239d0,&
                              -7.2362685049d0,-7.2023080858d0,-7.1752929217d0,-7.1511364876d0,&
                              -7.1346812422d0,-7.2679765578d0,-7.3233544837d0,-7.3357551761d0,&
                              -7.3419066457d0,-7.3436077720d0,-7.3421113394d0,-7.3383074168d0,&
                              -7.3187137789d0,-7.2765909672d0,-7.2362284481d0,-7.2022680290d0,&
                              -7.1752517624d0,-7.1510923884d0,-7.1346319981d0,-7.2679335611d0,&
                              -7.3233107520d0,-7.3357114444d0,-7.3418629140d0,-7.3435640403d0,&
                              -7.3420676077d0,-7.3382636851d0,-7.3186707822d0,-7.2765487055d0,&
                              -7.2361869214d0,-7.2022265023d0,-7.1752091332d0,-7.1510475543d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(4)%zrCut(2)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! third cut2d
   this%wyckoffSite(4)%zrcut(3)%alias='Bridge_90_54'
   this%wyckoffSite(4)%zrCut(3)%theta=1.57079632679d0
   this%wyckoffSite(4)%zrCut(3)%phi=0.942477796077d0
   gridPot1414(:,:)=reshape( [-6.6059602453d0,-6.6999852318d0,-6.7044785719d0,-6.6853839899d0,&
                              -6.6547934839d0,-6.6134240333d0,-6.5613171649d0,-6.4980344590d0,&
                              -6.2319002544d0,-5.4521493932d0,-4.1165565330d0,-2.3944015246d0,&
                              -1.3930856731d0,-2.5257612532d0,-6.7789617371d0,-6.8945967969d0,&
                              -6.9304090145d0,-6.9319638784d0,-6.9262486234d0,-6.9147740140d0,&
                              -6.8984782606d0,-6.8779552324d0,-6.7935721664d0,-6.5820665690d0,&
                              -6.3018904471d0,-6.0311177426d0,-5.8992990149d0,-6.0217481346d0,&
                              -7.0222896778d0,-7.1517156570d0,-7.2041786264d0,-7.2155856170d0,&
                              -7.2210711888d0,-7.2224312814d0,-7.2209069194d0,-7.2173642844d0,&
                              -7.2000171328d0,-7.1646561969d0,-7.1321749382d0,-7.1059245276d0,&
                              -7.0866524463d0,-7.0723352767d0,-7.1074573419d0,-7.2393683105d0,&
                              -7.2935853252d0,-7.3055398808d0,-7.3113495816d0,-7.3128217596d0,&
                              -7.3112143441d0,-7.3074199763d0,-7.2885734522d0,-7.2497860092d0,&
                              -7.2143835466d0,-7.1857572921d0,-7.1635621695d0,-7.1438546088d0,&
                              -7.1241768150d0,-7.2567805083d0,-7.3115259784d0,-7.3236598706d0,&
                              -7.3295849644d0,-7.3311049165d0,-7.3294761864d0,-7.3255917827d0,&
                              -7.3060859757d0,-7.2652097011d0,-7.2272142061d0,-7.1961889556d0,&
                              -7.1720979353d0,-7.1507866340d0,-7.1310562887d0,-7.2640300478d0,&
                              -7.3190989119d0,-7.3313625292d0,-7.3373912561d0,-7.3389861768d0,&
                              -7.3374022809d0,-7.3335307394d0,-7.3138694829d0,-7.2721244542d0,&
                              -7.2327438771d0,-7.2001990420d0,-7.1747291870d0,-7.1521529740d0,&
                              -7.1341483770d0,-7.2673407945d0,-7.3226121474d0,-7.3349646981d0,&
                              -7.3410709660d0,-7.3427313005d0,-7.3411992212d0,-7.3373647966d0,&
                              -7.3177178722d0,-7.2756542269d0,-7.2355552005d0,-7.2020409182d0,&
                              -7.1755622942d0,-7.1518813965d0,-7.1348561690d0,-7.2681356824d0,&
                              -7.3234801664d0,-7.3358661591d0,-7.3420032964d0,-7.3436919279d0,&
                              -7.3421837356d0,-7.3383695231d0,-7.3187564081d0,-7.2766435187d0,&
                              -7.2363497209d0,-7.2025076346d0,-7.1756163158d0,-7.1514231324d0,&
                              -7.1349134979d0,-7.2682128560d0,-7.3235797571d0,-7.3359760396d0,&
                              -7.3421230992d0,-7.3438205506d0,-7.3423204431d0,-7.3385135805d0,&
                              -7.3189155327d0,-7.2768022758d0,-7.2364695237d0,-7.2025476914d0,&
                              -7.1755428171d0,-7.1511890392d0,-7.1348804235d0,-7.2681750042d0,&
                              -7.3235492552d0,-7.3359492125d0,-7.3420999471d0,-7.3438003384d0,&
                              -7.3423035384d0,-7.3384992482d0,-7.3189078154d0,-7.2767982334d0,&
                              -7.2364573964d0,-7.2025142495d0,-7.1754755658d0,-7.1510696039d0,&
                              -7.1348072924d0,-7.2681055479d0,-7.3234842089d0,-7.3358860037d0,&
                              -7.3420385758d0,-7.3437408045d0,-7.3422458420d0,-7.3384433893d0,&
                              -7.3188563663d0,-7.2767500918d0,-7.2364074173d0,-7.2024547156d0,&
                              -7.1753994947d0,-7.1509663382d0,-7.1347220339d0,-7.2680221270d0,&
                              -7.3234022579d0,-7.3358047877d0,-7.3419580947d0,-7.3436614260d0,&
                              -7.3421675659d0,-7.3383662157d0,-7.3187821327d0,-7.2766791656d0,&
                              -7.2363364911d0,-7.2023797470d0,-7.1753164413d0,-7.1508696875d0,&
                              -7.1346794047d0,-7.2679794978d0,-7.3233592612d0,-7.3357621585d0,&
                              -7.3419154655d0,-7.3436187968d0,-7.3421249367d0,-7.3383235865d0,&
                              -7.3187406059d0,-7.2766387413d0,-7.2362956994d0,-7.2023389552d0,&
                              -7.1752734445d0,-7.1508226484d0,-7.1346367755d0,-7.2679361336d0,&
                              -7.3233158970d0,-7.3357184268d0,-7.3418717338d0,-7.3435750651d0,&
                              -7.3420812050d0,-7.3382802223d0,-7.3186972417d0,-7.2765961121d0,&
                              -7.2362538051d0,-7.2022963260d0,-7.1752293454d0,-7.1507756092d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(4)%zrCut(3)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Fourh cut2d
   this%wyckoffSite(4)%zrcut(4)%alias='Bridge_90_144'
   this%wyckoffSite(4)%zrCut(4)%theta=1.57079632679d0
   this%wyckoffSite(4)%zrCut(4)%phi=2.51327412287d0
   gridPot1414(:,:)=reshape( [-6.6769162103d0,-6.8143105457d0,-6.8780765051d0,-6.8963177678d0,&
                              -6.9093439337d0,-6.9188741362d0,-6.9260608343d0,-6.9316820111d0,&
                              -6.9436865458d0,-6.9589536855d0,-6.9728078137d0,-6.9831847207d0,&
                              -6.9870246577d0,-6.9796001915d0,-6.8134903007d0,-6.9478138625d0,&
                              -7.0070769271d0,-7.0225675027d0,-7.0325515595d0,-7.0387871850d0,&
                              -7.0424702024d0,-7.0444271040d0,-7.0449041102d0,-7.0413250934d0,&
                              -7.0393380574d0,-7.0382565247d0,-7.0354992229d0,-7.0268263821d0,&
                              -7.0271527161d0,-7.1585418442d0,-7.2128746193d0,-7.2251544064d0,&
                              -7.2314664705d0,-7.2336100587d0,-7.2328320755d0,-7.2300126672d0,&
                              -7.2148267435d0,-7.1838044329d0,-7.1570226271d0,-7.1367344270d0,&
                              -7.1216907232d0,-7.1080640733d0,-7.1080434937d0,-7.2401624634d0,&
                              -7.2945224330d0,-7.3065137379d0,-7.3123322586d0,-7.3137820195d0,&
                              -7.3121191125d0,-7.3082350763d0,-7.2889218358d0,-7.2489157852d0,&
                              -7.2123439590d0,-7.1831547049d0,-7.1611253218d0,-7.1422740203d0,&
                              -7.1243451269d0,-7.2570050467d0,-7.3117887361d0,-7.3239299782d0,&
                              -7.3298524995d0,-7.3313581193d0,-7.3297021947d0,-7.3257769993d0,&
                              -7.3060609862d0,-7.2645687929d0,-7.2258250816d0,-7.1942228667d0,&
                              -7.1699179654d0,-7.1488264250d0,-7.1311011229d0,-7.2640950941d0,&
                              -7.3191742480d0,-7.3314400703d0,-7.3374673272d0,-7.3390574705d0,&
                              -7.3374643872d0,-7.3335792485d0,-7.3138452283d0,-7.2718694138d0,&
                              -7.2321569903d0,-7.1992861887d0,-7.1736274423d0,-7.1510288121d0,&
                              -7.1341549919d0,-7.2673532892d0,-7.3226264796d0,-7.3349793978d0,&
                              -7.3410852983d0,-7.3427448978d0,-7.3412109810d0,-7.3373739839d0,&
                              -7.3177123598d0,-7.2755961629d0,-7.2354129806d0,-7.2017983726d0,&
                              -7.1752352252d0,-7.1515135357d0,-7.1348734412d0,-7.2681386223d0,&
                              -7.3234845764d0,-7.3358698340d0,-7.3420069713d0,-7.3436952354d0,&
                              -7.3421870430d0,-7.3383720956d0,-7.3187567756d0,-7.2766321264d0,&
                              -7.2363192189d0,-7.2024521431d0,-7.1755369372d0,-7.1513268491d0,&
                              -7.1349274627d0,-7.2682139585d0,-7.3235815946d0,-7.3359778770d0,&
                              -7.3421249367d0,-7.3438223880d0,-7.3423226480d0,-7.3385157854d0,&
                              -7.3189170027d0,-7.2767997033d0,-7.2364599689d0,-7.2025289492d0,&
                              -7.1755148876d0,-7.1511537598d0,-7.1348800560d0,-7.2681757392d0,&
                              -7.3235507251d0,-7.3359506825d0,-7.3421014171d0,-7.3438021759d0,&
                              -7.3423053759d0,-7.3385010857d0,-7.3189100203d0,-7.2767982334d0,&
                              -7.2364533540d0,-7.2025046947d0,-7.1754604986d0,-7.1510508617d0,&
                              -7.1348083949d0,-7.2681059154d0,-7.3234849438d0,-7.3358871062d0,&
                              -7.3420400457d0,-7.3437422745d0,-7.3422476794d0,-7.3384452267d0,&
                              -7.3188589388d0,-7.2767519292d0,-7.2364066823d0,-7.2024499382d0,&
                              -7.1753910424d0,-7.1509556809d0,-7.1347249739d0,-7.2680224945d0,&
                              -7.3234029929d0,-7.3358055227d0,-7.3419591972d0,-7.3436628959d0,&
                              -7.3421690359d0,-7.3383676856d0,-7.3187847051d0,-7.2766817380d0,&
                              -7.2363372261d0,-7.2023779095d0,-7.1753116638d0,-7.1508630726d0,&
                              -7.1346816097d0,-7.2679794978d0,-7.3233599961d0,-7.3357628935d0,&
                              -7.3419165680d0,-7.3436198992d0,-7.3421264067d0,-7.3383254239d0,&
                              -7.3187428109d0,-7.2766413138d0,-7.2362971693d0,-7.2023378527d0,&
                              -7.1752697696d0,-7.1508171360d0,-7.1346356730d0,-7.2679365010d0,&
                              -7.3233162644d0,-7.3357191618d0,-7.3418724688d0,-7.3435761675d0,&
                              -7.3420826750d0,-7.3382816922d0,-7.3186994467d0,-7.2765986845d0,&
                              -7.2362552751d0,-7.2022955910d0,-7.1752264054d0,-7.1507711993d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(4)%zrCut(4)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Fifth cut2d
   this%wyckoffSite(4)%zrcut(5)%alias='Bridge_135_54'
   this%wyckoffSite(4)%zrCut(5)%theta=2.35619449019d0
   this%wyckoffSite(4)%zrCut(5)%phi=0.942477796077d0
   gridPot1414(:,:)=reshape( [-6.6524169050d0,-6.7781863263d0,-6.8274892213d0,-6.8369484976d0,&
                              -6.8398568392d0,-6.8377022263d0,-6.8314265440d0,-6.8216280714d0,&
                              -6.7745113938d0,-6.6526102065d0,-6.5226774542d0,-6.4649644761d0,&
                              -6.5300482664d0,-6.7171592941d0,-6.8002509888d0,-6.9286968634d0,&
                              -6.9801139471d0,-6.9902630083d0,-6.9934620370d0,-6.9910229843d0,&
                              -6.9836290200d0,-6.9715513542d0,-6.9059141191d0,-6.6704281169d0,&
                              -6.2204734192d0,-5.5837439500d0,-5.0443605039d0,-5.1022491438d0,&
                              -7.0260785333d0,-7.1577601860d0,-7.2129686976d0,-7.2258658733d0,&
                              -7.2328570650d0,-7.2356783107d0,-7.2354901542d0,-7.2330624937d0,&
                              -7.2163794025d0,-7.1587663826d0,-7.0204808761d0,-6.6657627901d0,&
                              -5.7925348748d0,-3.1559945441d0,-7.1090416054d0,-7.2418514624d0,&
                              -7.2972157911d0,-7.3098443292d0,-7.3163971015d0,-7.3186825420d0,&
                              -7.3179538029d0,-7.3150983803d0,-7.2993468844d0,-7.2658902986d0,&
                              -7.2340874324d0,-7.2030827616d0,-7.1616262151d0,-7.0569255495d0,&
                              -7.1252370331d0,-7.2584173233d0,-7.3138867550d0,-7.3264414271d0,&
                              -7.3328310323d0,-7.3348618000d0,-7.3337901897d0,-7.3305099449d0,&
                              -7.3130830473d0,-7.2763318844d0,-7.2426015161d0,-7.2150950135d0,&
                              -7.1926900523d0,-7.1682153691d0,-7.1317181441d0,-7.2650527815d0,&
                              -7.3205622700d0,-7.3330790902d0,-7.3393856420d0,-7.3412852146d0,&
                              -7.3400331651d0,-7.3365217670d0,-7.3181232173d0,-7.2790572144d0,&
                              -7.2427992275d0,-7.2132965016d0,-7.1902866465d0,-7.1692652973d0,&
                              -7.1344710360d0,-7.2678383803d0,-7.3233210419d0,-7.3357941304d0,&
                              -7.3420312259d0,-7.3438345153d0,-7.3424571506d0,-7.3387899354d0,&
                              -7.3197291628d0,-7.2789653411d0,-7.2405703809d0,-7.2090174101d0,&
                              -7.1844909103d0,-7.1627716915d0,-7.1350016963d0,-7.2683624257d0,&
                              -7.3238017230d0,-7.3362417372d0,-7.3424369384d0,-7.3441876763d0,&
                              -7.3427474702d0,-7.3390056539d0,-7.3196435369d0,-7.2780888697d0,&
                              -7.2385649702d0,-7.2057195257d0,-7.1799865455d0,-7.1572320982d0,&
                              -7.1350141911d0,-7.2683385387d0,-7.3237572563d0,-7.3361811008d0,&
                              -7.3423590299d0,-7.3440895556d0,-7.3426247275d0,-7.3388557167d0,&
                              -7.3193859241d0,-7.2775519620d0,-7.2376069153d0,-7.2042058209d0,&
                              -7.1778561871d0,-7.1544046051d0,-7.1349351800d0,-7.2682591601d0,&
                              -7.3236683230d0,-7.3360870225d0,-7.3422576017d0,-7.3439800426d0,&
                              -7.3425063947d0,-7.3387263591d0,-7.3192172447d0,-7.2772807520d0,&
                              -7.2371773157d0,-7.2035535204d0,-7.1769301041d0,-7.1531242586d0,&
                              -7.1348436742d0,-7.2681628769d0,-7.3235650574d0,-7.3359793470d0,&
                              -7.3421451488d0,-7.3438624448d0,-7.3423829169d0,-7.3385962665d0,&
                              -7.3190625300d0,-7.2770639310d0,-7.2368631090d0,-7.2030989313d0,&
                              -7.1762906659d0,-7.1522272076d0,-7.1347481260d0,-7.2680621837d0,&
                              -7.3234588518d0,-7.3358705690d0,-7.3420330634d0,-7.3437463169d0,&
                              -7.3422627466d0,-7.3384716862d0,-7.3189225151d0,-7.2768864318d0,&
                              -7.2366279133d0,-7.2027784771d0,-7.1758537164d0,-7.1516120239d0,&
                              -7.1347029243d0,-7.2680133071d0,-7.3234074028d0,-7.3358176500d0,&
                              -7.3419786744d0,-7.3436904579d0,-7.3422050502d0,-7.3384125198d0,&
                              -7.3188578363d0,-7.2768085232d0,-7.2365308951d0,-7.2026527944d0,&
                              -7.1756883444d0,-7.1513841781d0,-7.1346599276d0,-7.2679662680d0,&
                              -7.3233570562d0,-7.3357662009d0,-7.3419261228d0,-7.3436364364d0,&
                              -7.3421495587d0,-7.3383555584d0,-7.3187960974d0,-7.2767368620d0,&
                              -7.2364441666d0,-7.2025447514d0,-7.1755494320d0,-7.1511963890d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(4)%zrCut(5)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Sixth cut2d
   this%wyckoffSite(4)%zrcut(6)%alias='Bridge_135_144'
   this%wyckoffSite(4)%zrCut(6)%theta=2.35619449019d0
   this%wyckoffSite(4)%zrCut(6)%phi=2.51327412287d0
   gridPot1414(:,:)=reshape( [-6.6740618902d0,-6.8098815170d0,-6.8718555793d0,-6.8891685540d0,&
                              -6.9012737818d0,-6.9099135482d0,-6.9162675066d0,-6.9211371597d0,&
                              -6.9316298271d0,-6.9481552637d0,-6.9695275689d0,-6.9941073552d0,&
                              -7.0183938819d0,-7.0433815857d0,-6.8076059988d0,-6.9389466178d0,&
                              -6.9949077555d0,-7.0086891200d0,-7.0169621281d0,-7.0215116946d0,&
                              -7.0235615720d0,-7.0239614047d0,-7.0204643389d0,-7.0137564845d0,&
                              -7.0144201773d0,-7.0223764063d0,-7.0347561515d0,-7.0511419407d0,&
                              -7.0229416109d0,-7.1522345575d0,-7.2042286055d0,-7.2152754527d0,&
                              -7.2203273825d0,-7.2211954016d0,-7.2191344994d0,-7.2150329072d0,&
                              -7.1960820151d0,-7.1593966335d0,-7.1283408811d0,-7.1059399623d0,&
                              -7.0915511314d0,-7.0825082249d0,-7.1067510199d0,-7.2382243040d0,&
                              -7.2918672942d0,-7.3034811836d0,-7.3089149388d0,-7.3099766268d0,&
                              -7.3079241770d0,-7.3036513329d0,-7.2831933510d0,-7.2413299895d0,&
                              -7.2027843570d0,-7.1713611114d0,-7.1469713191d0,-7.1254119598d0,&
                              -7.1237880072d0,-7.2561756144d0,-7.3106642067d0,-7.3226533066d0,&
                              -7.3284229507d0,-7.3297793683d0,-7.3279786514d0,-7.3239130735d0,&
                              -7.3038251572d0,-7.2618813147d0,-7.2227462231d0,-7.1905483016d0,&
                              -7.1652666032d0,-7.1424765091d0,-7.1308699697d0,-7.2637555303d0,&
                              -7.3187233338d0,-7.3309336646d0,-7.3369080024d0,-7.3384496366d0,&
                              -7.3368128217d0,-7.3328894637d0,-7.3130900297d0,-7.2712071910d0,&
                              -7.2318049318d0,-7.1992924361d0,-7.1738240512d0,-7.1510413069d0,&
                              -7.1340689984d0,-7.2672525961d0,-7.3224985920d0,-7.3348393829d0,&
                              -7.3409349935d0,-7.3425876107d0,-7.3410500189d0,-7.3372126544d0,&
                              -7.3175818997d0,-7.2756571668d0,-7.2358852094d0,-7.2028762303d0,&
                              -7.1770105852d0,-7.1540742287d0,-7.1348491866d0,-7.2681092229d0,&
                              -7.3234496645d0,-7.3358323497d0,-7.3419687520d0,-7.3436577510d0,&
                              -7.3421517637d0,-7.3383408586d0,-7.3187512632d0,-7.2767453143d0,&
                              -7.2366738499d0,-7.2031996244d0,-7.1768150787d0,-7.1533521045d0,&
                              -7.1349138654d0,-7.2682051386d0,-7.3235709373d0,-7.3359668522d0,&
                              -7.3421146469d0,-7.3438132007d0,-7.3423160332d0,-7.3385110080d0,&
                              -7.3189280275d0,-7.2768750395d0,-7.2366698075d0,-7.2029622238d0,&
                              -7.1762715562d0,-7.1524080143d0,-7.1348782186d0,-7.2681724317d0,&
                              -7.3235481527d0,-7.3359481101d0,-7.3420992121d0,-7.3438014409d0,&
                              -7.3423064783d0,-7.3385043931d0,-7.3189239851d0,-7.2768526224d0,&
                              -7.2365918989d0,-7.2027850920d0,-7.1759507346d0,-7.1518773540d0,&
                              -7.1348069249d0,-7.2681062829d0,-7.3234871488d0,-7.3358896786d0,&
                              -7.3420426182d0,-7.3437470519d0,-7.3422539268d0,-7.3384533116d0,&
                              -7.3188743735d0,-7.2767923535d0,-7.2364970857d0,-7.2026248650d0,&
                              -7.1756923869d0,-7.1514620866d0,-7.1347187265d0,-7.2680254344d0,&
                              -7.3234081378d0,-7.3358121376d0,-7.3419665471d0,-7.3436713483d0,&
                              -7.3421793257d0,-7.3383794454d0,-7.3188019773d0,-7.2767144449d0,&
                              -7.2363982300d0,-7.2024852175d0,-7.1754884281d0,-7.1511504524d0,&
                              -7.1346768323d0,-7.2679835402d0,-7.3233662435d0,-7.3357702433d0,&
                              -7.3419250203d0,-7.3436301890d0,-7.3421381664d0,-7.3383382862d0,&
                              -7.3187611856d0,-7.2766714482d0,-7.2363489859d0,-7.2024223762d0,&
                              -7.1754031697d0,-7.1510255047d0,-7.1346312631d0,-7.2679420134d0,&
                              -7.3233239818d0,-7.3357279816d0,-7.3418831261d0,-7.3435879273d0,&
                              -7.3420955372d0,-7.3382960245d0,-7.3187185563d0,-7.2766273490d0,&
                              -7.2363001093d0,-7.2023632098d0,-7.1753259961d0,-7.1509170942d0 ],shape( gridPot1414 ) )
   call this%wyckoffSite(4)%zrCut(6)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   ! Sixth cut2d, THE LAST ONE (thank goodness (>_<)Uu )
   this%wyckoffSite(4)%zrcut(7)%alias='Bridge_135_234'
   this%wyckoffSite(4)%zrCut(7)%theta=2.35619449019d0
   this%wyckoffSite(4)%zrCut(7)%phi=4.08407044967d0
   gridPot1414(:,:)=reshape( [-6.6335189320d0,-6.7466719423d0,-6.7809414232d0,-6.7825547186d0,&
                              -6.7777637091d0,-6.7683389771d0,-6.7554906779d0,-6.7400603712d0,&
                              -6.6836281071d0,-6.5747909375d0,-6.4733613294d0,-6.4107511363d0,&
                              -6.4123273148d0,-6.4961907454d0,-6.7815885789d0,-6.8965412037d0,&
                              -6.9309933288d0,-6.9318922173d0,-6.9257462601d0,-6.9142606259d0,&
                              -6.8985885086d0,-6.8795170787d0,-6.8064597873d0,-6.6435941270d0,&
                              -6.4368394824d0,-6.1999474510d0,-5.9866477515d0,-5.9159890885d0,&
                              -7.0153969744d0,-7.1397228821d0,-7.1849440295d0,-7.1918088035d0,&
                              -7.1920825860d0,-7.1875241997d0,-7.1793342450d0,-7.1683454617d0,&
                              -7.1235502390d0,-7.0221548078d0,-6.8890998328d0,-6.7113896500d0,&
                              -6.4628594747d0,-5.9963936726d0,-7.1044207452d0,-7.2345192370d0,&
                              -7.2863780475d0,-7.2969210615d0,-7.3011509088d0,-7.3008638966d0,&
                              -7.2973072969d0,-7.2913623585d0,-7.2647613593d0,-7.2080847122d0,&
                              -7.1471863026d0,-7.0828007496d0,-7.0102976380d0,-6.9030274596d0,&
                              -7.1225451450d0,-7.2542451724d0,-7.3078775054d0,-7.3193675494d0,&
                              -7.3245870561d0,-7.3253378448d0,-7.3228727001d0,-7.3180798531d0,&
                              -7.2953919220d0,-7.2474564695d0,-7.1997073360d0,-7.1553549427d0,&
                              -7.1129807655d0,-7.0605369057d0,-7.1301992945d0,-7.2627265492d0,&
                              -7.3172614456d0,-7.3292255560d0,-7.3349330937d0,-7.3361855107d0,&
                              -7.3342366940d0,-7.3299774472d0,-7.3090130595d0,-7.2645845951d0,&
                              -7.2217330442d0,-7.1846504024d0,-7.1531088239d0,-7.1201814284d0,&
                              -7.1337981559d0,-7.2668093992d0,-7.3218749559d0,-7.3341146862d0,&
                              -7.3401026213d0,-7.3416402131d0,-7.3399798786d0,-7.3360127889d0,&
                              -7.3159465547d0,-7.2731299157d0,-7.2322268141d0,-7.1978283430d0,&
                              -7.1702781087d0,-7.1447968615d0,-7.1347216664d0,-7.2679206988d0,&
                              -7.3231843344d0,-7.3355269628d0,-7.3416185310d0,-7.3432604908d0,&
                              -7.3417048919d0,-7.3378421703d0,-7.3180813230d0,-7.2757380153d0,&
                              -7.2352597359d0,-7.2013033592d0,-7.1743491990d0,-7.1501343335d0,&
                              -7.1348539640d0,-7.2681062829d0,-7.3234323923d0,-7.3358073601d0,&
                              -7.3419316352d0,-7.3436063020d0,-7.3420830424d0,-7.3382522928d0,&
                              -7.3185829513d0,-7.2763645913d0,-7.2359660579d0,-7.2020346708d0,&
                              -7.1750838180d0,-7.1508696875d0,-7.1348388968d0,-7.2681081204d0,&
                              -7.3234581168d0,-7.3358444770d0,-7.3419808793d0,-7.3436676734d0,&
                              -7.3421572761d0,-7.3383382862d0,-7.3187020191d0,-7.2765270233d0,&
                              -7.2361483346d0,-7.2022066576d0,-7.1752179531d0,-7.1509325289d0,&
                              -7.1347789954d0,-7.2680643887d0,-7.3234287174d0,-7.3358224274d0,&
                              -7.3419658121d0,-7.3436599560d0,-7.3421565411d0,-7.3383449011d0,&
                              -7.3187310511d0,-7.2765839848d0,-7.2362163208d0,-7.2022632516d0,&
                              -7.1752392677d0,-7.1508836523d0,-7.1347099067d0,-7.2679986074d0,&
                              -7.3233702860d0,-7.3357684059d0,-7.3419165680d0,-7.3436147543d0,&
                              -7.3421161168d0,-7.3383088867d0,-7.3187086340d0,-7.2765799424d0,&
                              -7.2362199958d0,-7.2022588417d0,-7.1752083982d0,-7.1507969239d0,&
                              -7.1346632350d0,-7.2679614906d0,-7.3233357416d0,-7.3357349640d0,&
                              -7.3418845961d0,-7.3435842524d0,-7.3420878199d0,-7.3382813247d0,&
                              -7.3186858494d0,-7.2765634052d0,-7.2362063985d0,-7.2022426720d0,&
                              -7.1751826737d0,-7.1507498847d0,-7.1346400829d0,-7.2679243738d0,&
                              -7.3232993598d0,-7.3356993171d0,-7.3418500517d0,-7.3435508105d0,&
                              -7.3420540105d0,-7.3382500878d0,-7.3186579200d0,-7.2765402531d0,&
                              -7.2361854514d0,-7.2022202549d0,-7.1751529068d0,-7.1507024781d0 ],shape( gridPot1414  ) )
   call this%wyckoffSite(4)%zrCut(7)%interRZ%read( gridR14,gridZ14,gridPot1414 )
   call this%interpol()
   return
end subroutine initialize_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT1D_PHI_PES_H2LiF001 #######################################
!#######################################################################
SUBROUTINE PLOT1D_PHI_PES_H2LiF001(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
   INTEGER,INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8),DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=*),PARAMETER :: routinename = "PLOT1D_PHI_PES_H2LiF001: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_PHI_PES_H2LiF001 ERR: Less than 2 points"
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
END SUBROUTINE PLOT1D_PHI_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT1D_ATOMIC_INTERAC_PHI_PES_H2LiF001 #######################################
!#######################################################################
SUBROUTINE PLOT1D_ATOMIC_INTERAC_PHI_PES_H2LiF001(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
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
   CHARACTER(LEN=*),PARAMETER :: routinename = "PLOT1D_PHI_PES_H2LiF001: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_PHI_PES_H2LiF001 ERR: Less than 2 points"
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
END SUBROUTINE PLOT1D_ATOMIC_INTERAC_PHI_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT1D_THETA_PES_H2LiF001 #########################################
!#######################################################################
SUBROUTINE PLOT1D_THETA_PES_H2LiF001(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=*),PARAMETER :: routinename = "PLOT1D_THETA_PES_H2LiF001: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_THETA_PES_H2LiF001 ERR: Less than 2 points"
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
END SUBROUTINE PLOT1D_THETA_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT1D_ATOMIC_INTERAC_THETA_PES_H2LiF001 ########################
!#######################################################################
SUBROUTINE PLOT1D_ATOMIC_INTERAC_THETA_PES_H2LiF001(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
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
   CHARACTER(LEN=*),PARAMETER :: routinename = "PLOT1D_THETA_PES_H2LiF001: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_THETA_PES_H2LiF001 ERR: Less than 2 points"
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
END SUBROUTINE PLOT1D_ATOMIC_INTERAC_THETA_PES_H2LiF001
!#######################################################################
!# SUBROUTINE: PLOT1D_R_PES_H2LiF001 ###################################
!#######################################################################
SUBROUTINE PLOT1D_R_PES_H2LiF001(thispes,npoints,X,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=*),PARAMETER :: routinename = "PLOT1D_R_PES_H2LiF001: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_R_PES_H2LiF001 ERR: Less than 2 points"
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
END SUBROUTINE PLOT1D_R_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT1D_Z_PES_H2LiF001 #########################################
!#######################################################################
SUBROUTINE PLOT1D_Z_PES_H2LiF001(thispes,npoints,X,L,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
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
   CHARACTER(LEN=*),PARAMETER :: routinename = "PLOT1D_Z_PES_H2LiF001: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_Z_PES_H2LiF001 ERR: Less than 2 points"
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
END SUBROUTINE PLOT1D_Z_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT_XYMAP_PES_H2LiF001
!#######################################################################
SUBROUTINE PLOT_XYMAP_PES_H2LiF001(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
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
END SUBROUTINE PLOT_XYMAP_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT_RZMAP_PES_H2LiF001
!#######################################################################
SUBROUTINE PLOT_RZMAP_PES_H2LiF001(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
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
END SUBROUTINE PLOT_RZMAP_PES_H2LiF001
!#######################################################################
! SUBROUTINE: PLOT_ATOMIC_INTERAC_RZ_PES_H2LiF001
!#######################################################################
SUBROUTINE PLOT_ATOMIC_INTERAC_RZ_PES_H2LiF001(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(PES_H2LiF001),INTENT(IN) :: thispes
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
END SUBROUTINE PLOT_ATOMIC_INTERAC_RZ_PES_H2LiF001
!###########################################################
!# FUNCTION: is_allowed_PES_H2LiF001
!###########################################################
LOGICAL FUNCTION is_allowed_PES_H2LiF001(this,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_H2LiF001),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8) :: zmin, rmin, rmax
   ! Run section
   zmin=this%wyckoffsite(1)%zrcut(1)%getfirstZ()
   rmin=this%wyckoffsite(1)%zrcut(1)%getfirstR()
   rmax=this%wyckoffsite(1)%zrcut(1)%getlastR()
   SELECT CASE(size(x)/=6)
      CASE(.TRUE.)
         WRITE(0,*) "is_allowed_PES_H2LiF001 ERR: wrong number of dimensions"
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(x(3)<zmin .OR. x(4)<rmin .OR. x(4)>rmax)
      CASE(.TRUE.)
         is_allowed_PES_H2LiF001=.FALSE.
      CASE(.FALSE.)
         is_allowed_PES_H2LiF001=.TRUE.
   END SELECT
   RETURN
END FUNCTION is_allowed_PES_H2LiF001
!###########################################################
! FUNCTION: from_molecular_to_atomic
!###########################################################
pure function from_molecular_to_atomic(molcoord) result(atomcoord)
   ! Initial declarations
   implicit none
   ! i/o variables
   real(kind=8),dimension(6),intent(in):: molcoord
   ! dymmy function variable
   real(kind=8),dimension(6):: atomcoord
   real(kind=8),dimension(2),parameter:: masa=[ 1.d0,1.d0 ]
   real(kind=8):: mTot
   ! run section
   mTot=sum(masa(:))
   atomcoord(1)=molcoord(1)+(masa(2)/(mTot))*molcoord(4)*dcos(molcoord(6))*dsin(molcoord(5))
   atomcoord(2)=molcoord(2)+(masa(2)/(mTot))*molcoord(4)*dsin(molcoord(6))*dsin(molcoord(5))
   atomcoord(3)=molcoord(3)+(masa(2)/(mTot))*molcoord(4)*dcos(molcoord(5))
   atomcoord(4)=molcoord(1)-(masa(1)/(mTot))*molcoord(4)*dcos(molcoord(6))*dsin(molcoord(5))
   atomcoord(5)=molcoord(2)-(masa(1)/(mTot))*molcoord(4)*dsin(molcoord(6))*dsin(molcoord(5))
   atomcoord(6)=molcoord(3)-(masa(1)/(mTot))*molcoord(4)*dcos(molcoord(5))
   return
end function from_molecular_to_atomic

end module PES_H2LiF001_MOD
