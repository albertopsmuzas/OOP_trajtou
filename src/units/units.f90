!#########################################################
! RAISON D'ÊTRE:
! - This module manages all input/output units
! - (Other) routines inside the program will use always a.u. for
!   calculations.
! UPDATES:
! - Last big update 06/Nov/2013 --> Alberto Muzas
! FUNCTIONALITY:
! - Numbers with a magnitude associated will be declared as
!   a "Quantity" subtype. Only I/O data recomended.
! - Several Routines to change units to a.u.
! IDEAS FOR THE FUTURE:
! - Routines to go to non-au from au units.
! - Units control of the output
!##########################################################
MODULE UNITS_MOD
IMPLICIT NONE
!=========================================================
! Quantity derived data
!----------------------
TYPE :: Quantity ! Some physical quantity
PRIVATE
   CHARACTER(LEN=10) :: units ! units for this physical quantity
   REAL(KIND=8),PUBLIC :: mag ! Magnitude
CONTAINS
   PROCEDURE, PUBLIC :: READ
   PROCEDURE, PUBLIC :: GET
END TYPE Quantity
!==========================================================
! Length derived data
!--------------------
TYPE, EXTENDS(Quantity) :: length
CONTAINS
   PROCEDURE, PUBLIC :: TO_STD => LENGTH_AU
END TYPE length
!============================================================
! Energy derived data
!--------------------
TYPE, EXTENDS(Quantity) :: energy
CONTAINS
   PROCEDURE,PUBLIC :: TO_STD => ENERGY_AU
END TYPE
!==========================================================
! Mass derived data
!------------------
TYPE, EXTENDS(Quantity) :: mass
CONTAINS
   PROCEDURE,PUBLIC :: TO_STD => MASS_AU
END TYPE
!==========================================================
! Time derived data
!------------------
TYPE, EXTENDS(Quantity) :: time
CONTAINS
   PROCEDURE :: TO_STD => TIME_AU
END TYPE
!==========================================================
! Angle derived data
!------------------
TYPE, EXTENDS(Quantity) :: angle
CONTAINS
   PROCEDURE :: TO_STD => TO_RAD
   PROCEDURE :: TO_DEG
END TYPE
!==========================================================
! Conversion factors:

REAL*8, PARAMETER, PRIVATE :: au2ev = 27.21138386D0
REAL*8, PARAMETER, PRIVATE :: au2kcalmol = 627.503D0
REAL*8, PARAMETER, PRIVATE :: au2kjmol = 2.6255D3
REAL*8, PARAMETER, PRIVATE :: au2angst = 0.52917720859D0
REAL*8, PARAMETER, PRIVATE :: au2fs = 0.02418884326505D0
REAL*8, PARAMETER, PRIVATE :: au2ps = 2.418884326505D-5
REAL*8, PARAMETER, PRIVATE :: hmass2au = 1837.15264409D0
REAL*8, PARAMETER, PRIVATE :: dmass2au = 3671.482934845D0
REAL*8, PARAMETER, PRIVATE :: pmass2au = 1836.15267247D0

!==========================================================

! MODULE CONTAINS:
CONTAINS
!###########################################################
!# SUBROUTINE: READ     ####################################
!###########################################################
! - Reads a quantity from unit "u"
! - Stores Kind and units
! - Kind is only needed if "this" is a pure Quantity type variable
!-----------------------------------------------------------
SUBROUTINE READ(this,u)
	! Initial declarations
	IMPLICIT NONE
	! I/O variables
	CLASS(Quantity), INTENT(OUT) :: this
	INTEGER, INTENT(IN) :: u ! unit to read
	! Run section
	! Check variable type
	SELECT TYPE (this)
	TYPE IS (Quantity)
		WRITE(0,*) "READ_QUANT ERR: Quantities not allowed. Allocate it with a specific sub-type before (length, time, etc.)"
		CALL EXIT(1)
	END SELECT
	READ(u,*) this%mag, this%units
# ifdef VERBOSE
   WRITE(*,*) "Quantity read from unit:", u
#endif
	RETURN
END SUBROUTINE READ
!###########################################################
!# SUBROUTINE: GET #########################################
!###########################################################
! - Reads a quantity from arguments
! - Stores units
!-----------------------------------------------------------
SUBROUTINE GET(this,mag,units)
	! Initial declarations
	IMPLICIT NONE
	! I/O variables
	CLASS(Quantity), INTENT(OUT) :: this
	REAL(KIND=8), INTENT(IN) :: mag ! magnitude
	CHARACTER(LEN=*) :: units
	! Run section
	this%mag=mag
	this%units=units
	RETURN
END SUBROUTINE GET
!###########################################################
!# SUBROUTINE: GO TO RAD ####################################
!############################################################
! - Changes degrees to radians
!------------------------------------------------------------
SUBROUTINE TO_RAD(this)
   ! Initial declarations
   USE CONSTANTS_MOD
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
      WRITE(0,*) "Supported ones: deg"
      CALL EXIT(1)
   END IF
   this%units = "rad"
   RETURN
END SUBROUTINE TO_RAD
!############################################################
!# SUBROUTINE: TO DEG ####################################
!############################################################
! - Changes radians to degrees
!------------------------------------------------------------
SUBROUTINE TO_DEG(this)
	! Initial declarations
	USE CONSTANTS_MOD
        IMPLICIT NONE
        ! I/O variables
        CLASS(angle), INTENT(INOUT) :: this
        ! GO, GO, GO !!!------
        IF (this%units.EQ."deg") THEN
                RETURN
        ELSE IF (this%units.EQ."rad") THEN
                this%mag = this%mag*180.D0/pi
        ELSE
                WRITE(0,*) "TO_DEG ERR: incorrect kind"
                WRITE(0,*) "Supported ones: rad"
                CALL EXIT(1)
        END IF
        this%units = "deg"
        RETURN
END SUBROUTINE TO_DEG
!############################################################
!# SUBROUTINE: LENGTH_AU ####################################
!############################################################
! - Manage length unit changes to au
!------------------------------------------------------------
SUBROUTINE LENGTH_AU(this)
        IMPLICIT NONE
        ! I/O variables
        CLASS(length), INTENT(INOUT) :: this
        ! Run section
        IF (this%units.EQ."angst") THEN
                this%mag = this%mag / au2angst
        ELSE
                WRITE(0,*) "LENGTH_AU ERR: incorrect units"
                WRITE(0,*) "Supported ones: angst"
                CALL EXIT(1)
        END IF
        this%units = "au"
        RETURN
END SUBROUTINE LENGTH_AU
!############################################################
!# SUBROUTINE: MASS_AU ######################################
!############################################################
! - Manage mass unit changes to au
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
        ELSE
                WRITE(0,*) "MASS_AU ERR: incorrect units"
                WRITE(0,*) "Supported ones: hmass, dmass, pmass"
                CALL EXIT(1)
        END IF
        this%units = "au"
        RETURN
END SUBROUTINE MASS_AU
!############################################################
!# SUBROUTINE: ENERGY_AU ####################################
!############################################################
! - Manage energy unit changes to au
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
        ELSE
                WRITE(0,*) "ENERGY_AU ERR: incorrect units"
                WRITE(0,*) "Supported ones: ev, kcalmol, kjmol"
                CALL EXIT(1)
        END IF
        this%units = "au"
        RETURN
END SUBROUTINE ENERGY_AU
!############################################################
!# SUBROUTINE: TIME_AU ######################################
!############################################################
! - Manage energy unit changes to au
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
        ELSE
                WRITE(0,*) "TIME_AU ERR: incorrect units"
                WRITE(0,*) "Supported ones: fs, ps"
                CALL EXIT(1)
        END IF
        this%units = "au"
        RETURN
END SUBROUTINE TIME_AU
END MODULE UNITS_MOD
