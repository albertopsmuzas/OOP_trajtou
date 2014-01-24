MODULE PES_MOD
IMPLICIT NONE
!====================================================================================================
! Type for a generic Potential energy surface
!--------------------------------------------
TYPE PES
PRIVATE
	CHARACTER(LEN=30) :: alias ! Just an alias for the PES
	INTEGER :: dimensions ! number of coordinates which this PES depends on
	REAL(KIND=8), DIMENSION(:),ALLOCATABLE :: r ! last point where the potential was calculated
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: dvdu ! last calculated derivatives
	REAL(KIND=8) :: v ! last value of the potential (obviously at point r)
	LOGICAL :: initialized = .FALSE.
CONTAINS
	! Set block
   PROCEDURE,PUBLIC :: INITIALIZE => INIT_PES! very simple one, SHOULD be overridden by child types
	PROCEDURE,PUBLIC :: SET_ALIAS
	PROCEDURE,PUBLIC :: SET_DIMENSIONS
	! Get block
	PROCEDURE,PUBLIC :: get_alias
	PROCEDURE,PUBLIC :: get_dimensions
	PROCEDURE,PUBLIC :: get_last_value
	! Enquire block
	PROCEDURE,PUBLIC :: is_initialized
END TYPE PES
! MODULE CONTAINS:
CONTAINS
!###############################################################
! SUBROUTINE: INIT_PES #########################################
!###############################################################
! - Gives initial values to some variables.
! - Some important arrays are allocated as well
! - If no argument is given, default initialization will be loaded.
!---------------------------------------------------------------
SUBROUTINE INIT_PES(this,alias,dimensions)
	! Initial declarations
	IMPLICIT NONE
	! I/O variables
	CLASS(PES), INTENT(INOUT) :: this
	CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: alias
	INTEGER, INTENT(IN), OPTIONAL :: dimensions
	! Local variables
	! Run section
	IF(.NOT.this%is_initialized()) THEN
		CALL this%SET_ALIAS(alias)
		CALL this%SET_DIMENSIONS(dimensions)
		this%initialized=.TRUE.
	END IF
   RETURN
END SUBROUTINE INIT_PES
!###############################################################
! SUBROUTINE: SET_ALIAS ########################################
!###############################################################
! - Sets an alias for PES.
! - If no argument is given a default string will be loaded
!---------------------------------------------------------------
SUBROUTINE SET_ALIAS(this,alias)
	! Initial declarations
	IMPLICIT NONE
	! I/O variables
	CLASS(PES), INTENT(INOUT) :: this
	CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: alias
	! Run section
	IF (PRESENT(alias)) THEN
		this%alias=alias
	ELSE
		this%alias="Testing PES (1D)"
	END IF
	RETURN
END SUBROUTINE SET_ALIAS
!###############################################################
! SUBROUTINE: SET_DIMENSIONS ###################################
!###############################################################
! - Sets number of dimensions and allocates some arrays
! - If no argument is given, default value will be 1
!---------------------------------------------------------------
SUBROUTINE SET_DIMENSIONS(this,dimensions)
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
END SUBROUTINE SET_DIMENSIONS
!#########################################################
! FUNCTION: is_initiallized ##############################
!#########################################################
! - Enquires whether the PES is initiallized or not
!---------------------------------------------------------
LOGICAL FUNCTION is_initialized(this)
	! Initial declarations
	IMPLICIT NONE
	! I/O variables
	CLASS(PES), INTENT(IN) :: this
	! Run section
	IF(this%initialized) THEN
		is_initialized=.TRUE.
	ELSE
		is_initialized=.FALSE.
	END IF
	RETURN
END FUNCTION is_initialized
!##########################################################
! FUNCTION: get_alias #####################################
!##########################################################
CHARACTER(LEN=30) FUNCTION get_alias(this)
	! Initial declarations
	IMPLICIT NONE
	! I/O VAriables
	CLASS(PES), INTENT(IN) :: this
	! Run section
	IF (this%is_initialized()) THEN
		get_alias=this%alias
	ELSE
		WRITE(0,*) "PES_get_alias: argument is not initialized."
		CALL EXIT(1)
	END IF
	RETURN
END FUNCTION get_alias
!##########################################################
! FUNCTION: get_dimensions ################################
!##########################################################
INTEGER FUNCTION get_dimensions(this)
	! Initial declarations
	IMPLICIT NONE
	! I/O variables
	CLASS(PES), INTENT(IN) :: this
	! Run section
	IF(this%is_initialized()) THEN
		get_dimensions=this%dimensions
	ELSE
		WRITE(0,*) "PES_get_dimensions: argument is not initialized."
		CALL EXIT(1)
	END IF
	RETURN
END FUNCTION get_dimensions
!##########################################################
! FUNCTION: get_last_value ################################
!##########################################################
REAL(KIND=8) FUNCTION get_last_value(this)
	! Initial declarations
	IMPLICIT NONE
	! I/O variables
	CLASS(PES), INTENT(IN) :: this
	! Run part
	IF(this%is_initialized()) THEN
		get_last_value=this%v
	ELSE
		WRITE(0,*) "PES_get_value: argument is not initialized."
		CALL EXIT(1)
	END IF
END FUNCTION get_last_value
END MODULE PES_MOD
