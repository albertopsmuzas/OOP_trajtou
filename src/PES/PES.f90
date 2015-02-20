MODULE PES_MOD
   USE SURFACE_MOD
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
!> @param surf - Surface specifications
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
