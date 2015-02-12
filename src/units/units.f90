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
   PROCEDURE,PUBLIC:: TO_CELSIUS => TO_CELSIUS_UNITS
   PROCEDURE,PUBLIC:: TO_FARENHEIT => TO_FARENHEIT_UNITS
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
!! simple units change function. From au to angstroem units
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
      case('Farenheit')
         
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
