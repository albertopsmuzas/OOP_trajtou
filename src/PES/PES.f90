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
	PROCEDURE,PUBLIC :: SET_ALIAS
	PROCEDURE,PUBLIC :: SET_DIMENSIONS
	! Get block
	PROCEDURE,PUBLIC :: getalias
	PROCEDURE,PUBLIC :: getdimensions
	PROCEDURE,PUBLIC :: getlastvalue
	! Enquire block
	PROCEDURE,PUBLIC :: is_initialized
END TYPE PES
! MODULE CONTAINS:
CONTAINS
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
		this%alias="Mysterious PES"
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
CHARACTER(LEN=30) FUNCTION getalias(this)
	! Initial declarations
	IMPLICIT NONE
	! I/O VAriables
	CLASS(PES), INTENT(IN) :: this
	! Run section
	getalias=this%alias
	RETURN
END FUNCTION getalias
!##########################################################
! FUNCTION: get_dimensions ################################
!##########################################################
INTEGER FUNCTION getdimensions(this)
	! Initial declarations
	IMPLICIT NONE
	! I/O variables
	CLASS(PES), INTENT(IN) :: this
	! Run section
	getdimensions=this%dimensions
	RETURN
END FUNCTION getdimensions
!##########################################################
! FUNCTION: get_last_value ################################
!##########################################################
REAL(KIND=8) FUNCTION getlastvalue(this)
	! Initial declarations
	IMPLICIT NONE
	! I/O variables
	CLASS(PES), INTENT(IN) :: this
	! Run part
	getlastvalue=this%v
END FUNCTION getlastvalue
END MODULE PES_MOD
