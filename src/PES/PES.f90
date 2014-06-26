MODULE PES_MOD
   USE SURFACE_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: Atomdata
!> @brief
!! Type to store atoms characteristics
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 25/Mar/2013
!> @version 1.0
!----------------------------------------------------------------
TYPE :: Atomdata
PRIVATE
   REAL(KIND=8) :: mass
   CHARACTER(LEN=2) :: myname
   CONTAINS
      PROCEDURE,PUBLIC :: getname => getname_atomdata
      PROCEDURE,PUBLIC :: getmass => getmass_atomdata
END TYPE Atomdata
!////////////////////////////////////////////////////////////////////////////////////////
! TYPE: PES
!> @brief
!! Type for a generic Potential energy surface
!
!> param alias - An alias for the PES
!> param dimensions - Number of coordinates which this PES depends on
!> param r - last point qhere the potential was calculated
!> param dvdu - last calculated derivatives
!> param v - last value of the potential (at r, obviously)
!> param initialized - Controls status of the PES
!> param atomdat - Contains array of atoms (with some info) related to this PES
!----------------------------------------------------------------------------------------
TYPE,ABSTRACT :: PES
PRIVATE
   CHARACTER(LEN=30) :: alias
   INTEGER :: dimensions
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: r
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: dvdu 
	REAL(KIND=8) :: v
	LOGICAL :: initialized = .FALSE.
   TYPE(Atomdata),DIMENSION(:),ALLOCATABLE,PUBLIC:: atomdat
   TYPE(Surface),PUBLIC:: surf
   CONTAINS
      ! Initialization block
      PROCEDURE(INITIALIZE_PES),DEFERRED,PUBLIC :: INITIALIZE
      ! Set block
      PROCEDURE,NON_OVERRIDABLE,PUBLIC :: SET_ALIAS => SET_ALIAS_PES
      PROCEDURE,NON_OVERRIDABLE,PUBLIC :: SET_DIMENSIONS => SET_DIMENSIONS_PES
      PROCEDURE,NON_OVERRIDABLE,PUBLIC :: SET_ATOMS => SET_ATOMS_PES
      ! Get block
      PROCEDURE,NON_OVERRIDABLE,PUBLIC :: getalias => getalias_PES
      PROCEDURE,NON_OVERRIDABLE,PUBLIC :: getdimensions => getdimensions_PES
      PROCEDURE,NON_OVERRIDABLE,PUBLIC :: getlastvalue => getlastvalue_PES
      PROCEDURE(GET_V_AND_DERIVS_PES),DEFERRED :: GET_V_AND_DERIVS
      ! Enquire block
      PROCEDURE,NON_OVERRIDABLE,PUBLIC :: is_initialized => is_initialized_PES
      PROCEDURE(is_allowed_PES),DEFERRED,PUBLIC :: is_allowed
END TYPE PES
ABSTRACT INTERFACE
   !###########################################################
   !# SUBROUTINE: INITIALIZE_PES 
   !###########################################################
   !-----------------------------------------------------------
   SUBROUTINE INITIALIZE_PES(this,filename)
      IMPORT PES
      CLASS(PES),INTENT(OUT) :: this
      CHARACTER(LEN=*),INTENT(IN) :: filename
   END SUBROUTINE INITIALIZE_PES
   !###########################################################
   !# SUBROUTINE: GET_V_AND_DERIVS_PES 
   !###########################################################
   !-----------------------------------------------------------
   SUBROUTINE GET_V_AND_DERIVS_PES(this,x,v,dvdu)
      IMPORT PES
      CLASS(PES),TARGET,INTENT(IN):: this
      REAL(KIND=8),DIMENSION(:),INTENT(IN):: x
      REAL(KIND=8),INTENT(OUT):: v
      REAL(KIND=8),DIMENSION(:),INTENT(OUT):: dvdu
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
!###########################################################
!# FUNCTION: getname_atomdata 
!###########################################################
!> @brief
!! Common get function. Gets name of Atomdata type variable
!-----------------------------------------------------------
CHARACTER(LEN=2) FUNCTION getname_atomdata(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Atomdata),INTENT(IN) :: this
   ! Run section
   getname_atomdata=this%myname
   RETURN
END FUNCTION getname_atomdata  
!###########################################################
!# FUNCTION: getmass_atomdata 
!###########################################################
!> @brief 
!! Common get function. Gets mass of Atomdata
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getmass_atomdata(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Atomdata),INTENT(IN):: this
   ! Run section
   getmass_atomdata=this%mass
   RETURN
END FUNCTION getmass_atomdata
!###########################################################
!# SUBROUTINE: SET_ATOMS_PES 
!###########################################################
!> @brief
!! Sets atoms of the PES
!-----------------------------------------------------------
SUBROUTINE SET_ATOMS_PES(this,n,names,masses)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES),INTENT(INOUT) :: this
   INTEGER(KIND=4),INTENT(IN) :: n
   CHARACTER(LEN=2),DIMENSION(:) :: names
   REAL(KIND=8),DIMENSION(:) :: masses
   ! Local variables
   INTEGER(KIND=4) :: i ! counter
   ! Run section
   ! Check data
   SELECT CASE(n==size(names) .AND. n==size(masses))
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         WRITE(0,*) "SET_ATOMS_PES ERR: Size mismatch between declared atoms, names and masses"
         CALL EXIT(1)
   END SELECT
   ! Set information
   ALLOCATE(this%atomdat(n))
   DO i = 1, n
      this%atomdat(i)%mass=masses(i)
      this%atomdat(i)%myname=names(i)
   END DO
   RETURN
END SUBROUTINE SET_ATOMS_PES
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
	CLASS(PES), INTENT(INOUT) :: this
	CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: alias
	! Run section
	IF (PRESENT(alias)) THEN
		this%alias=alias
	ELSE
		this%alias="Mysterious PES"
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
CHARACTER(LEN=30) FUNCTION getalias_PES(this)
	! Initial declarations
	IMPLICIT NONE
	! I/O VAriables
	CLASS(PES), INTENT(IN) :: this
	! Run section
	getalias_PES=this%alias
	RETURN
END FUNCTION getalias_PES
!##########################################################
! FUNCTION: get_dimensions ################################
!##########################################################
!> @brief
!! Common get function. Gets dimensions from PES derived type
!----------------------------------------------------------
INTEGER FUNCTION getdimensions_PES(this)
	! Initial declarations
	IMPLICIT NONE
	! I/O variables
	CLASS(PES), INTENT(IN) :: this
	! Run section
	getdimensions_PES=this%dimensions
	RETURN
END FUNCTION getdimensions_PES
!##########################################################
! FUNCTION: get_last_value ################################
!##########################################################
!> @brief
!! Common get function. Gets last value calculated of the PES
!----------------------------------------------------------
REAL(KIND=8) FUNCTION getlastvalue_PES(this)
	! Initial declarations
	IMPLICIT NONE
	! I/O variables
	CLASS(PES), INTENT(IN) :: this
	! Run part
	getlastvalue_PES=this%v
END FUNCTION getlastvalue_PES
END MODULE PES_MOD
