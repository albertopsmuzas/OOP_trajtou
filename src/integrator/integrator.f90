!#########################################################
! MODULE: INTEGRATOR_MOD
!> @brief
!! Contains generic type for integration schemes
!##########################################################
MODULE INTEGRATOR_MOD
! Initial declarations
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: INTEGRATOR
!> @brief
!! Generic type to implement different integration schemes
!
!> @param nv - @b integer(kind=4). Number of dimensions of the function to integrate. It is size(errscal) by definition.
!> @param dt - @b real(kind=8). Time step to be used in the next integration. Mandatory to set it before integration
!> @param dt_old - @b real(kind=8). Time step used in the previous integration
!> @param err - @b real(kind=8). Error threshold during integration. Mandatory to set it before integration.
!> @param errscal(:) - @b real(kind=8) @b array. Array against which the integration error is scaled. This 
!!                     parameter can be used to scale in a different way errors in each component of the function
!!                     to integrate. Mandatory to set it before integration.
!> @param realparam(:) - @b real(kind=8) @b array. Real internal parameters of the integrator. Its meaning can change
!!                       from each specific implementation.
!> @param intparam(:) - @b integer(kind=8) @b array. Integer internal parameters of the integrator. Its meaning can change
!!                      from each specific implementation.
!> @param stringparam(:) - @b character(len=*) @b array. String internal parameters of the integrator.
!!                         Its meaning can change from each specific implementation. Max number of characters: 20.
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!----------------------------------------------------------------
TYPE,ABSTRACT :: INTEGRATOR
   PRIVATE
   INTEGER(KIND=4):: nv
   REAL(KIND=8):: dt
   REAL(KIND=8) :: err
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: errscal
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: realparam
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE:: intparam
   CHARACTER(LEN=20),DIMENSION(:),ALLOCATABLE:: stringparam
CONTAINS
   ! Set block
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_TIMESTEP => SET_TIMESTEP_INTEGRATOR ! Mandatory
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_ERROR => SET_ERR_INTEGRATOR ! Mandatory
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_ERRSCALING => SET_ERRSCAL_INTEGRATOR ! Mandatory
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_REALPARAM => SET_REALPARAM_INTEGRATOR
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_INTPARAM => SET_INTPARAM_INTEGRATOR
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_STRINGPARAM => SET_STRINGPARAM_INTEGRATOR
   ! Tools block
   PROCEDURE(INTEGRATE_INTEGRATOR),PUBLIC,DEFERRED:: INTEGRATE ! Deferred subroutine!!!!!
END TYPE INTEGRATOR

ABSTRACT INTERFACE
   !###########################################################
   !# SUBROUTINE: INTEGRATE_INTEGRATOR
   !###########################################################
   !-----------------------------------------------------------
   SUBROUTINE INTEGRATE_INTEGRATOR(this)
      IMPORT Integrator
      CLASS(Integrator),INTENT(INOUT):: this
   END SUBROUTINE INTEGRATE_INTEGRATOR
END INTERFACE
CONTAINS
!###########################################################
!# SUBROUTINE: SET_ERRSCAL_INTEGRATOR
!###########################################################
!> @brief
!! Standard set subroutine. Sets errscal allocatable atribute
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_ERRSCAL_INTEGRATOR(this,errscal)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Integrator),INTENT(INOUT):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN):: errscal
   ! Run section
   SELECT CASE(allocated(this%errscal))
      CASE(.TRUE.)
         DEALLOCATE(this%errscal)
      CASE(.FALSE.)
         ! nothing special 
   END SELECT
   this%nv=size(errscal)
   ALLOCATE(this%errscal(this%nv))
   this%errscal=errscal
   RETURN
END SUBROUTINE SET_ERRSCAL_INTEGRATOR
!###########################################################
!# SUBROUTINE: SET_ERR_INTEGRATOR
!###########################################################
!> @brief
!! Standard set subroutine. Sets eps atribute
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_ERR_INTEGRATOR(this,err)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Integrator),INTENT(INOUT):: this
   REAL(KIND=8),INTENT(IN):: err
   ! Run section
   this%err=err
   RETURN
END SUBROUTINE SET_ERR_INTEGRATOR
!###########################################################
!# SUBROUTINE: SET_TIMESTEP_INTEGRATOR
!###########################################################
!> @brief
!! Standard set subroutine. Sets timestep atribute
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_TIMESTEP_INTEGRATOR(this,timestep)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Integrator),INTENT(INOUT):: this
   REAL(KIND=8),INTENT(IN):: timestep
   ! Run section
   this%dt=timestep
   RETURN
END SUBROUTINE SET_TIMESTEP_INTEGRATOR
!###########################################################
!# SUBROUTINE: SET_REALPARAM_INTEGRATOR
!###########################################################
!> @brief
!! Standard set subroutine. Sets realparam allocatable atribute.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_REALPARAM_INTEGRATOR(this,realparam)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Integrator),INTENT(INOUT):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN):: realparam
   ! Run section
   SELECT CASE(allocated(this%realparam))
      CASE(.TRUE.)
         DEALLOCATE(this%realparam)
      CASE(.FALSE.)
         ! nothing special 
   END SELECT
   ALLOCATE(this%realparam(size(realparam)))
   this%realparam=realparam
   RETURN
END SUBROUTINE SET_REALPARAM_INTEGRATOR
!###########################################################
!# SUBROUTINE: SET_INTPARAM_INTEGRATOR
!###########################################################
!> @brief
!! Standard set subroutine. Sets intparam allocatable atribute.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_INTPARAM_INTEGRATOR(this,intparam)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Integrator),INTENT(INOUT):: this
   INTEGER(KIND=4),DIMENSION(:),INTENT(IN):: intparam
   ! Run section
   SELECT CASE(allocated(this%intparam))
      CASE(.TRUE.)
         DEALLOCATE(this%intparam)
      CASE(.FALSE.)
         ! nothing special 
   END SELECT
   ALLOCATE(this%intparam(size(intparam)))
   this%intparam=intparam
   RETURN
END SUBROUTINE SET_INTPARAM_INTEGRATOR
!###########################################################
!# SUBROUTINE: SET_REALPARAM_INTEGRATOR
!###########################################################
!> @brief
!! Standard set subroutine. Sets stringparam allocatable atribute.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_STRINGPARAM_INTEGRATOR(this,stringparam)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Integrator),INTENT(INOUT):: this
   CHARACTER(LEN=*),DIMENSION(:),INTENT(IN):: stringparam
   ! Run section
   SELECT CASE(allocated(this%stringparam))
      CASE(.TRUE.)
         DEALLOCATE(this%stringparam)
      CASE(.FALSE.)
         ! nothing special 
   END SELECT
   ALLOCATE(this%stringparam(size(stringparam)))
   this%stringparam=stringparam
   RETURN
END SUBROUTINE SET_STRINGPARAM_INTEGRATOR
END MODULE INTEGRATOR_MOD
