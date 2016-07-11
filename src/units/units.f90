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
module UNITS_MOD
implicit none
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
type :: Quantity
   private
   character(len=10) :: units
   real(kind=8) :: mag
contains
   procedure,public :: READ => READ_QUANTITY_FROM_ARGUMENTS
   procedure,public :: getvalue => get_magnitude_quantity
   procedure,public :: getunits => get_units_quantity
end type Quantity
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
type, extends(Quantity) :: Length
contains
   procedure,public :: TO_STD => LENGTH_AU
   procedure,public :: TO_ANGST => TO_ANGST_LENGTH
end type Length
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
type,extends(Quantity) :: Temperature
contains
   procedure,public:: TO_STD => TO_KELVIN_UNITS
end type Temperature
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
type, extends(Quantity) :: Energy
contains
   procedure,public :: TO_STD => ENERGY_AU
end type Energy
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
type, extends(Quantity) :: Mass
contains
   procedure,public :: TO_STD => MASS_AU
end type Mass
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
type, extends(Quantity) :: Time
contains
   procedure,public :: TO_STD => TIME_AU
end type Time
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
type, extends(Quantity) :: angle
contains
   procedure,public :: TO_STD => TO_RAD
   procedure,public :: TO_DEG
end type
!//////////////////////////////////////////////////////////////////////
! Conversion factors:
real(kind=8),parameter:: au2ev = 27.21138386D0
real(kind=8),parameter:: au2kcalmol = 627.503D0
real(kind=8),parameter:: au2kjmol = 2.6255D3
real(kind=8),parameter:: au2angst = 0.52917720859D0
real(kind=8),parameter:: au2fs = 0.02418884326505D0
real(kind=8),parameter:: au2ps = 2.418884326505D-5
real(kind=8),parameter:: hmass2au = 1837.15264409D0
real(kind=8),parameter:: dmass2au = 3671.482934845d0
real(kind=8),parameter:: pmass2au = 1836.15267247d0
real(kind=8),parameter:: dalton2au=1822.88848325d0
real(kind=8),parameter:: kelvinParam=273.15d0
real(kind=8),parameter:: fahrenheitParam=459.67d0
real(kind=8),parameter:: boltzmann=8.617332478d-5/au2ev
real(kind=8),parameter:: pi=dacos(-1.D0)
!//////////////////////////////////////////////////////////////////////

! MODULE CONTAINS:
contains
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
subroutine READ_QUANTITY_FROM_ARGUMENTS(quant,mag,units)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(Quantity),intent(out) :: quant
   real(kind=8),intent(in) :: mag
   character(len=*),intent(in) :: units
   ! Run section
   ! Check variable type
   select type (quant)
   type is (Quantity)
      write(0,*) "READ_QUANT_FROM_ARGUMENTS ERR: Pure Quantities're not allowed. &
      & Declare it with a specific sub-type instead (length, time, etc.)"
      call EXIT(1)
   end select
   quant%mag=mag
   quant%units=units
   return
end subroutine READ_QUANTITY_FROM_ARGUMENTS
!###########################################################
!# FUNCTION: get_magnitude_quantity ########################
!###########################################################
! - Typical get function
!-----------------------------------------------------------
real(kind=8) function get_magnitude_quantity(quant) 
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(Quantity),intent(in) :: quant
   ! Run section
   get_magnitude_quantity=quant%mag
   return
end function get_magnitude_quantity
!###########################################################
!# FUNCTION: get_units_quantity ############################
!###########################################################
! - Typical get function
!-----------------------------------------------------------
character(len=10) function get_units_quantity(quant) 
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(Quantity),intent(in) :: quant
   ! Run section
   get_units_quantity=quant%units
   return
end function get_units_quantity
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
subroutine TO_RAD(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(angle), intent(inout) :: this
   ! GO, GO, GO !!!------
   if (this%units.eq."rad") then
      return
   else if (this%units.eq."deg") then
      this%mag = this%mag*(pi/180.D0)
   else
      write(0,*) "TO_RAD ERR: incorrect kind"
      write(0,*) "Supported ones: deg, rad"
      call EXIT(1)
   end if
   this%units = "rad"
   return
end subroutine TO_RAD
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
subroutine TO_DEG(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(angle), intent(inout) :: this
   ! GO, GO, GO !!!------
   select case(this%units)
      case("deg")
         return
      case("rad")
         this%mag = this%mag*180.D0/pi
         this%units = "deg"
      case default
         write(0,*) "TO_DEG ERR: incorrect kind"
         write(0,*) "Supported ones: rad, deg"
         call EXIT(1)
   end select
   return
end subroutine TO_DEG
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
subroutine LENGTH_AU(this)
   implicit none
   ! I/O variables
   class(length), intent(inout) :: this
   ! Run section
   select case(this%units)
      case("angst")
         this%mag = this%mag / au2angst
         this%units = "au"
      case("au")
         ! do nothing
      case default
         write(0,*) "TO_ANGST ERR: incorrect units"
         write(0,*) "Supported ones: angst, au"
         write(0,*) 'Encountered one: '//this%units
         call EXIT(1)
   end select
   return
end subroutine LENGTH_AU
!###########################################################
!# SUBROUTINE: TO_ANGST_LENGTH 
!###########################################################
!> @brief
!! simple units change function. From au to angstroem units
!-----------------------------------------------------------
subroutine TO_ANGST_LENGTH(this)
   ! Initial declarations     
   implicit none
   ! I/O variables
   class(Length),intent(inout):: this
   ! Run section
   select case(this%units)
      case('angst','Angst')
         ! do nothing
      case('au','bohr','bohrs','Bohr','Bohrs')
         this%mag = this%mag*au2angst
         this%units = "angst"
      case default
         write(0,*) "TO_ANGST ERR: incorrect units"
         write(0,*) "Supported ones: angst, au"
         call EXIT(1)
   end select
   return
end subroutine TO_ANGST_LENGTH
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
subroutine ENERGY_AU(this)
   implicit none
   ! I/O variables
   class(energy), intent(inout) :: this
   ! GO, GO, GO !!!
   select case ( this%units )
   case('ev','eV','electronVolts','electronvolts')
      this%mag = this%mag/au2ev
   case('kcalmol')
      this%mag = this%mag/au2kcalmol
   case('kjmol')
      this%mag = this%mag/au2kjmol
   case('au')
      return
   case default
      write(0,*) "ENERGY_AU ERR: incorrect units"
      write(0,*) "Supported ones: ev, kcalmol, kjmol, au"
      call EXIT(1)
   end select
   this%units = "au"
   return
end subroutine ENERGY_AU
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
subroutine TIME_AU(this)
   ! Initial declarations
        implicit none
        ! I/O variables
        class(time), intent(inout) :: this
        ! GO, GO, GO !!!
        if (this%units.eq."ps") then
                this%mag = this%mag/au2ps
        else if (this%units.eq."fs") then
                this%mag = this%mag/au2fs
        else if (this%units.eq."au") then
                return
        else
                write(0,*) "TIME_AU ERR: incorrect units"
                write(0,*) "Supported ones: fs, ps, au"
                call EXIT(1)
        end if
        this%units = "au"
        return
end subroutine TIME_AU
end module UNITS_MOD
