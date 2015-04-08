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
! MODULE: HLIF001_WS_WS_MOD
!> @brief
!! CRP3D specific implementation for H/LiF001
!##########################################################
module PES_HLIF001_WS_MOD
! Initial declarations
use LiF001SURF_MOD, only: LiF001Surf,pi
use PES_MOD, only: PES
use CUBICSPLINES_MOD, only: Csplines
use FOURIER_P4MM_MOD, only: Fourierp4mm
use LOGISTIC_FUNCTION_MOD, only: Logistic_func
implicit none
! Local module variable, used to simulate SYSTEM_MOD
type(LiF001Surf):: sysLiF001Surf
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
! TYPE: PES_HLIF001_WS
!------------------------------------------------------------------------------
type,extends(PES) :: PES_HLIF001_WS
   integer(kind=4):: max_order=2
   type(Pair_pot),dimension(:),allocatable:: all_pairpots
   type(Sitio),dimension(:),allocatable:: all_sites
   integer(kind=4),dimension(:,:),allocatable:: klist
   type(Logistic_func) dampFunc
   contains
      ! Initialization block
      procedure,public:: initialize => initialize_PES_HLIF001_WS
      ! Get block
      procedure,public:: get_v_and_derivs => GET_V_AND_DERIVS_PES_HLIF001_WS
      procedure,public:: get_v_and_derivs_correction => GET_V_AND_DERIVS_CORRECTION_PES_HLIF001_WS
      procedure,public:: get_repul_corrections => GET_REPUL_CORRECTIONS_PES_HLIF001_WS
      procedure,public:: getpot => getpot_crp3d
      ! Enquire block
      procedure,public:: is_allowed => is_allowed_PES_HLIF001_WS
      ! Tools block
      procedure,public:: extract_vasint => EXTRACT_VASINT_PES_HLIF001_WS
      procedure,public:: smooth => SMOOTH_PES_HLIF001_WS
      procedure,public:: interpol => INTERPOL_Z_PES_HLIF001_WS
      ! Plot tools
      procedure,public:: plot_xymap => PLOT_XYMAP_PES_HLIF001_WS
      procedure,public:: plot_direction1d => PLOT_DIRECTION1D_PES_HLIF001_WS
      procedure,public:: plot_sitios => PLOT_SITIOS_PES_HLIF001_WS
      procedure,public:: plot_pairpots => PLOT_PAIRPOTS_PES_HLIF001_WS
      procedure,public:: plot_z => PLOT_Z_PES_HLIF001_WS
end type PES_HLIF001_WS
!///////////////////////////////////////////////////////////////////////////
contains
!###########################################################
!# SUBROUTINE: GET_REPUL_CORRECTIONS_PES_HLIF001_WS
!###########################################################
SUBROUTINE GET_REPUL_CORRECTIONS_PES_HLIF001_WS(this,P,v,dvdz,dvdx,dvdy)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_WS),INTENT(IN) :: this
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
END SUBROUTINE GET_REPUL_CORRECTIONS_PES_HLIF001_WS
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
!# SUBROUTINE: INITIALIZE_PES_HLIF001_WS
!###########################################################
subroutine initialize_PES_HLIF001_WS(this,filename,tablename)
   implicit none
   ! I/O variables
   class(PES_HLIF001_WS),intent(out):: this
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
   ! Create surface
   select case( invoke )
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
   ! Set switch function
   call this%dampFunc%read([5.d0,5.d0])
   ! Set Fourier Coefficients
   allocate( this%klist(6,2) )
   this%klist(:,1)=[0,1,1,2,2,2]
   this%klist(:,2)=[0,0,1,0,1,2]
   call this%interpol()
   return
end subroutine initialize_PES_HLIF001_WS
!#######################################################################
!# SUBROUTINE: EXTRACT_VASINT_PES_HLIF001_WS ######################################
!#######################################################################
SUBROUTINE EXTRACT_VASINT_PES_HLIF001_WS(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_WS),INTENT(INOUT):: this
   ! Local variables
   INTEGER(KIND=4):: npairpots, nsites
   INTEGER(KIND=4):: i,j ! counters
   REAL(KIND=8):: control_vasint
   character(len=*),parameter:: routinename="EXTRACT_VASINT_PES_HLIF001_WS: "
   ! Run section ------------------------
   npairpots=size(this%all_pairpots)
   control_vasint=this%all_pairpots(1)%vasint
   DO i = 1, npairpots
      IF (this%all_pairpots(1)%vasint/=control_vasint) THEN
         WRITE(0,*) "EXTRACT_VASINT_PES_HLIF001_WS ERR: Incoherences in vasint values found"
         WRITE(0,*) "EXTRACT_VASINT_PES_HLIF001_WS ERR: vasint's value at pairpot",1,control_vasint
         WRITE(0,*) "EXTRACT_VASINT_PES_HLIF001_WS ERR: vasint's value at pairpot",i,control_vasint
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
END SUBROUTINE EXTRACT_VASINT_PES_HLIF001_WS
!############################################################
!# SUBROUTINE: SMOOTH_PES_HLIF001_WS ############################
!############################################################
SUBROUTINE SMOOTH_PES_HLIF001_WS(this)
   ! Initial declaraitons
   IMPLICIT NONE
   CLASS(PES_HLIF001_WS),INTENT(INOUT):: this
   ! Local variables
   REAL(KIND=8),DIMENSION(3):: A
   INTEGER(KIND=4):: i,j ! counters
   INTEGER(KIND=4):: npairpots,nsites
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: v,dvdzr,dummy
   CHARACTER(LEN=*),PARAMETER:: routinename="SMOOTH_PES_HLIF001_WS: "
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
END SUBROUTINE SMOOTH_PES_HLIF001_WS
!############################################################
!# SUBROUTINE: INTERPOL_Z_PES_HLIF001_WS ########################
!############################################################
SUBROUTINE INTERPOL_Z_PES_HLIF001_WS(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_WS),INTENT(INOUT) :: this
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
END SUBROUTINE INTERPOL_Z_PES_HLIF001_WS
!##################################################################
!# SUBROUTINE: INTERACTION_AP #####################################
!##################################################################
SUBROUTINE INTERACTION_AP(A,P,pairpot,dampfunc,interac,dvdz_corr,dvdx_corr,dvdy_corr)
   IMPLICIT NONE
   ! I/O VAriables --------------------------------------------
   REAL(KIND=8),DIMENSION(3),INTENT(IN):: A, P
   TYPE(Pair_pot),INTENT(IN):: pairpot
   type(Logistic_func),INTENT(IN):: dampfunc
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
   type(Logistic_func),INTENT(IN) :: dampfunc
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
!# SUBROUTINE: GET_V_AND_DERIVS_PES_HLIF001_WS ##################
!############################################################
!> @brief
!! Subroutine that calculates the 3D potential for a point A and
!! its derivatives in cartesian coordinates.
!
!> @param[in] this - PES_HLIF001_WS PES
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
SUBROUTINE GET_V_AND_DERIVS_PES_HLIF001_WS(this,X,v,dvdu,errCode)
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_WS),TARGET,INTENT(IN):: this
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
   CHARACTER(LEN=*),PARAMETER:: routinename="GET_V_AND_DERIVS_PES_HLIF001_WS: "
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
END SUBROUTINE GET_V_AND_DERIVS_PES_HLIF001_WS
!############################################################
!# SUBROUTINE: GET_V_AND_DERIVS_CORRECTION_PES_HLIF001_WS ############
!############################################################
!> @brief
!! Subroutine that calculates the correction to the 3D PES for a point A and
!! its derivatives in cartesian coordinates.
!
!> @param[in] this - PES_HLIF001_WS PES
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
SUBROUTINE GET_V_AND_DERIVS_CORRECTION_PES_HLIF001_WS(this,X,v,dvdu)
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_WS),TARGET,INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(3), INTENT(IN) :: X
   REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(3),INTENT(OUT) :: dvdu
   ! Local variables
   INTEGER(KIND=4) :: npairpots
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: pot,dvdz,dvdx,dvdy
   INTEGER :: i ! counters
   ! Pointers
   REAL(KIND=8), POINTER :: zmax
   CHARACTER(LEN=24),PARAMETER :: routinename="GET_V_AND_DERIVS_PES_HLIF001_WS: "
   zmax => this%all_sites(1)%interz%x(this%all_sites(1)%n)
   npairpots = size(this%all_pairpots)
   SELECT CASE(size(v)==npairpots+1)
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         WRITE(*,*) "GET_V_AND_DERIVS_CORRECTION_PES_HLIF001_WS ERR: wrong number of dimensions array v"
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
END SUBROUTINE GET_V_AND_DERIVS_CORRECTION_PES_HLIF001_WS
!############################################################
!# FUNCTION: getpot_crp3d ###################################
!############################################################
!> @brief
!! Subroutine that calculates the 3D potential for a point A
!
!> @param[in] this - PES_HLIF001_WS PES
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
REAL(KIND=8) FUNCTION getpot_crp3d(this,X)
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_WS),TARGET,INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(3), INTENT(IN) :: X
   ! Local variables
   CLASS(Fourierp4mm),ALLOCATABLE:: interpolxy
   REAL(KIND=8):: v
   INTEGER(KIND=4) :: nsites,npairpots
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: pot,dummy
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: potarr
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: f,derivarr ! arguments to the xy interpolation
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: xy ! arguments to the xy interpolation
   INTEGER :: i ! counters
   ! Pointers
   REAL(KIND=8), POINTER :: zmax
   CHARACTER(LEN=24),PARAMETER :: routinename="GET_V_AND_DERIVS_PES_HLIF001_WS: "
   zmax => this%all_sites(1)%interz%x(this%all_sites(1)%n)
   npairpots = size(this%all_pairpots)
   nsites = size(this%all_sites)
   ! GABBA, GABBA HEY! ----------------------
   ALLOCATE(pot(npairpots))
   ALLOCATE(dummy(npairpots))
   SELECT CASE(X(3)>zmax)
      CASE(.TRUE.)
         getpot_crp3d=0.D0
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   !
   CALL this%GET_REPUL_CORRECTIONS(X,pot,dummy,dummy,dummy)
   v=sum(pot)
   ! Now, we have all the repulsive interaction and corrections to the derivarives
   ! stored in v(:) and dvdu(:) respectively.
   ! Let's get v and derivatives from xy interpolation of the corrugationless function
   ALLOCATE(f(2,nsites))
   ALLOCATE(xy(nsites,2))
   ALLOCATE(potarr(2))
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
   getpot_crp3d=v+potarr(1)
   RETURN
END FUNCTION getpot_crp3d
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
! SUBROUTINE: PLOT_XYMAP_PES_HLIF001_WS
!#######################################################################
SUBROUTINE PLOT_XYMAP_PES_HLIF001_WS(this,filename,init_xyz,nxpoints,nypoints,Lx,Ly)
   IMPLICIT NONE
   CLASS(PES_HLIF001_WS),INTENT(IN) :: this
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
END SUBROUTINE PLOT_XYMAP_PES_HLIF001_WS
!#######################################################################
! SUBROUTINE: PLOT_DIRECTION1D_PES_HLIF001_WS ###################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. To define
!! the direction, the angle alpha is given.
!
!> @param[in] this - PES_HLIF001_WS PES used
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
SUBROUTINE PLOT_DIRECTION1D_PES_HLIF001_WS(this,filename,npoints,angle,z,L)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_HLIF001_WS),INTENT(IN) :: this
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*), INTENT(IN) :: filename
   REAL*8, INTENT(IN) :: z, angle
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL*8 :: delta,L,v,s,alpha
   REAL*8 :: xmax, xmin, ymax, ymin
   REAL*8, DIMENSION(3) :: r, dvdu
   INTEGER :: i ! Counter
   CHARACTER(LEN=24), PARAMETER :: routinename = "PLOT_DIRECTION1D_PES_HLIF001_WS: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT_DIRECTION1D_PES_HLIF001_WS ERR: Less than 2 points"
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
END SUBROUTINE PLOT_DIRECTION1D_PES_HLIF001_WS
!#######################################################################
!# SUBROUTINE: PLOT_Z_PES_HLIF001_WS #######################################
!#######################################################################
SUBROUTINE PLOT_Z_PES_HLIF001_WS(this,npoints,xyz,L,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_HLIF001_WS),INTENT(IN):: this
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
   CHARACTER(LEN=*),PARAMETER:: routinename = "PLOT_DIRECTION1D_PES_HLIF001_WS: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT_Z_PES_HLIF001_WS ERR: Less than 2 points"
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
END SUBROUTINE PLOT_Z_PES_HLIF001_WS
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
!# SUBROUTINE: PLOT_PAIRPOTS_PES_HLIF001_WS
!###########################################################
SUBROUTINE PLOT_PAIRPOTS_PES_HLIF001_WS(this,npoints)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_WS),INTENT(IN)::this
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
END SUBROUTINE PLOT_PAIRPOTS_PES_HLIF001_WS
!###########################################################
!# SUBROUTINE: PLOT_SITIOS_PES_HLIF001_WS
!###########################################################
SUBROUTINE PLOT_SITIOS_PES_HLIF001_WS(this,npoints)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_WS),INTENT(IN)::this
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
END SUBROUTINE PLOT_SITIOS_PES_HLIF001_WS
!###########################################################
!# FUNCTION: is_allowed_PES_HLIF001_WS
!###########################################################
LOGICAL FUNCTION is_allowed_PES_HLIF001_WS(this,x)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_HLIF001_WS),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8) :: xmin,xmax
   ! Run section
   xmin=this%all_sites(1)%z(1)
   xmax=this%all_sites(1)%z(this%all_sites(1)%n)
   SELECT CASE(size(x)/=3)
      CASE(.TRUE.)
         WRITE(0,*) "is_allowed_PES_HLIF001_WS ERR: array doesn't have 3 dimensions: 3"
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE( x(3)<xmin )
      CASE(.TRUE.)
         is_allowed_PES_HLIF001_WS=.FALSE.
      CASE(.FALSE.)
         is_allowed_PES_HLIF001_WS=.TRUE.
   END SELECT
   RETURN
END FUNCTION is_allowed_PES_HLIF001_WS
END MODULE PES_HLIF001_WS_MOD
