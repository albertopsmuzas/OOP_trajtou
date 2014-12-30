!#########################################################
! MODULE: INTEGRATOR_O1_MOD
!> @brief
!! Contains an abstract type for integration of first order
!! differential equations. 
!! @details
!! - Any system of ordinary differential
!!   equations can be reduced to a set of first order differential
!!   equations.
!##########################################################
MODULE INTEGRATOR_O1_MOD
! Initial declarations
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: INTEGRATOR_O1
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
!!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!----------------------------------------------------------------
TYPE,ABSTRACT :: Integrator_o1
   PRIVATE
   INTEGER(KIND=4):: nv
   REAL(KIND=8):: err
   REAL(KIND=8):: dt
   REAL(KIND=8):: dt_old
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE,PUBLIC:: errscal
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE,PUBLIC:: realparam
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE,PUBLIC:: intparam
   CHARACTER(LEN=20),DIMENSION(:),ALLOCATABLE,PUBLIC:: stringparam
CONTAINS
   ! Set block
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_TIMESTEP_NEXT => SET_TIMESTEP_NEXT_INTEGRATOR_O1 ! Mandatory
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_TIMESTEP_USED => SET_TIMESTEP_USED_INTEGRATOR_O1 ! Mandatory
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_ERROR => SET_ERR_INTEGRATOR_O1 ! Mandatory
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_ERRSCALING => SET_ERRSCAL_INTEGRATOR_O1 ! Mandatory
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_REALPARAM => SET_REALPARAM_INTEGRATOR_O1
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_INTPARAM => SET_INTPARAM_INTEGRATOR_O1
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_STRINGPARAM => SET_STRINGPARAM_INTEGRATOR_O1
   ! Get block
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: getnv => getnv_INTEGRATOR_O1
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: geterr => geterr_INTEGRATOR_O1
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: gettimestep_next => gettimestep_next_INTEGRATOR_O1
   PROCEDURE,PUBLIC,NON_OVERRIDABLE:: gettimestep_used => gettimestep_used_INTEGRATOR_O1
   ! Tools block
   PROCEDURE(INTEGRATE_INTEGRATOR_O1),PUBLIC,DEFERRED:: INTEGRATE ! Deferred subroutine!!!!!
END TYPE Integrator_o1

ABSTRACT INTERFACE
   !###########################################################
   !# SUBROUTINE: INTEGRATE_INTEGRATOR_O1
   !###########################################################
   !-----------------------------------------------------------
   SUBROUTINE INTEGRATE_INTEGRATOR_O1(this,f,dfdx,x,DERIV,switch)
      IMPORT integrator_o1
      CLASS(integrator_o1),INTENT(INOUT):: this
      REAL(KIND=8),DIMENSION(:),INTENT(INOUT):: f
      REAL(KIND=8),DIMENSION(:),INTENT(IN):: dfdx
      REAL(KIND=8),INTENT(INOUT):: x
      EXTERNAL DERIV
      LOGICAL,INTENT(INOUT) :: switch
   END SUBROUTINE INTEGRATE_INTEGRATOR_O1
END INTERFACE
CONTAINS
!###########################################################
!# FUNCTION: geterr_INTEGRATOR_O1
!###########################################################
!> @brief
!! Standard get function. Gets @b err atribute
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION geterr_INTEGRATOR_O1(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Integrator_o1),INTENT(IN):: this
   ! Run section
   geterr_INTEGRATOR_O1=this%err
   RETURN
END FUNCTION geterr_INTEGRATOR_O1
!###########################################################
!# FUNCTION: getnv_INTEGRATOR_O1
!###########################################################
!> @brief
!! Standard get function. Gets @b nv atribute
!-----------------------------------------------------------
INTEGER(KIND=4) FUNCTION getnv_INTEGRATOR_O1(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Integrator_o1),INTENT(IN):: this
   ! Run section
   getnv_INTEGRATOR_O1=this%nv
   RETURN
END FUNCTION getnv_INTEGRATOR_O1
!###########################################################
!# FUNCTION: gettimestep_next_INTEGRATOR_O1
!###########################################################
!> @brief
!! Standard get function. Gets dt atribute.
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION gettimestep_next_INTEGRATOR_O1(this)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Integrator_o1),INTENT(IN):: this
   ! Run section
   gettimestep_next_INTEGRATOR_O1=this%dt
   RETURN
END FUNCTION gettimestep_next_INTEGRATOR_O1
!###########################################################
!# FUNCTION: gettimestep_used_INTEGRATOR_O1
!###########################################################
!> @brief
!! Standard get function. Gets dt_old atribute.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION gettimestep_used_INTEGRATOR_O1(this)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Integrator_o1),INTENT(IN):: this
   ! Run section
   gettimestep_used_INTEGRATOR_O1=this%dt_old
   RETURN
END FUNCTION gettimestep_used_INTEGRATOR_O1
!###########################################################
!# SUBROUTINE: SET_ERRSCAL_INTEGRATOR_O1
!###########################################################
!> @brief
!! Standard set subroutine. Sets errscal allocatable atribute
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_ERRSCAL_INTEGRATOR_O1(this,errscal)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(integrator_o1),INTENT(INOUT):: this
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
END SUBROUTINE SET_ERRSCAL_INTEGRATOR_O1
!###########################################################
!# SUBROUTINE: SET_ERR_INTEGRATOR_O1
!###########################################################
!> @brief
!! Standard set subroutine. Sets eps atribute
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_ERR_INTEGRATOR_O1(this,err)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(integrator_o1),INTENT(INOUT):: this
   REAL(KIND=8),INTENT(IN):: err
   ! Run section
   this%err=err
   RETURN
END SUBROUTINE SET_ERR_INTEGRATOR_O1
!###########################################################
!# SUBROUTINE: SET_TIMESTEP_NEXT_INTEGRATOR_O1
!###########################################################
!> @brief
!! Standard set subroutine. Sets dt atribute
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_TIMESTEP_NEXT_INTEGRATOR_O1(this,timestep)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(integrator_o1),INTENT(INOUT):: this
   REAL(KIND=8),INTENT(IN):: timestep
   ! Run section
   this%dt=timestep
   RETURN
END SUBROUTINE SET_TIMESTEP_NEXT_INTEGRATOR_O1
!###########################################################
!# SUBROUTINE: SET_TIMESTEP_USED_INTEGRATOR_O1
!###########################################################
!> @brief
!! Standard set subroutine. Sets dt_old atribute
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_TIMESTEP_USED_INTEGRATOR_O1(this,timestep)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(integrator_o1),INTENT(INOUT):: this
   REAL(KIND=8),INTENT(IN):: timestep
   ! Run section
   this%dt_old=timestep
   RETURN
END SUBROUTINE SET_TIMESTEP_USED_INTEGRATOR_O1
!###########################################################
!# SUBROUTINE: SET_REALPARAM_INTEGRATOR_O1
!###########################################################
!> @brief
!! Standard set subroutine. Sets realparam allocatable atribute.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_REALPARAM_INTEGRATOR_O1(this,realparam)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(integrator_o1),INTENT(INOUT):: this
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
END SUBROUTINE SET_REALPARAM_INTEGRATOR_O1
!###########################################################
!# SUBROUTINE: SET_INTPARAM_INTEGRATOR_O1
!###########################################################
!> @brief
!! Standard set subroutine. Sets intparam allocatable atribute.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_INTPARAM_INTEGRATOR_O1(this,intparam)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(integrator_o1),INTENT(INOUT):: this
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
END SUBROUTINE SET_INTPARAM_INTEGRATOR_O1
!###########################################################
!# SUBROUTINE: SET_REALPARAM_INTEGRATOR_O1
!###########################################################
!> @brief
!! Standard set subroutine. Sets stringparam allocatable atribute.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_STRINGPARAM_INTEGRATOR_O1(this,stringparam)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(integrator_o1),INTENT(INOUT):: this
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
END SUBROUTINE SET_STRINGPARAM_INTEGRATOR_O1
END MODULE INTEGRATOR_O1_MOD
