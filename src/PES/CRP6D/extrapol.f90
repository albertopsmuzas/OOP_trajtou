!#########################################################
! MODULE: EXTRAPOL_TO_VACUUM_MOD
!> @brief
!! Provides tools to extrapolate the potential in the vacuum,
!! which is only dependent on the distance between both atoms
!##########################################################
MODULE EXTRAPOL_TO_VACUUM_MOD
! Initial declarations
USE CUBICSPLINES_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: Vacuumpot
!> @brief
!! Class that stores all information that can be extracted from the
!! monodimensional vacuum potential
!----------------------------------------------------------------
TYPE :: Vacuumpot
   PRIVATE
   INTEGER(KIND=4) :: n
   TYPE(Csplines) :: rpot
   REAL(KIND=8) :: surfen
   REAL(KIND=8) :: potmin
   CONTAINS
      ! Initialization block
      PROCEDURE,PUBLIC :: INITIALIZE => INITIALIZE_VACUUMPOT
      PROCEDURE,PUBLIC :: READ => READ_VACUUMPOT
      ! Get block
      PROCEDURE,PUBLIC :: getpot => getpot_vacuumpot
      PROCEDURE,PUBLIC :: getderiv => getderiv_vacuumpot
      PROCEDURE,PUBLIC :: getscalefactor => getscalefactor_vacuumpot
      ! Tools block
      PROCEDURE,PUBLIC :: SHIFTPOT => SHIFTPOT_VACUUMPOT
      ! Plot tools block
      PROCEDURE,PUBLIC :: PLOT => PLOT_VACUUMPOT
      PROCEDURE,PUBLIC :: PLOTDATA => PLOTDATA_VACUUMPOT
END TYPE Vacuumpot
!/////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SHIFTPOT_VACUUMPOT 
!###########################################################
!> @brief
!! Shift the entire potential so that the equilibrium geomtry has
!! energy 0
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 27/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SHIFTPOT_VACUUMPOT(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(INOUT):: this
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   ! Run section
   DO i = 1, this%n
      this%rpot%f(i)=this%rpot%f(i)-this%potmin
   END DO
   CALL this%rpot%REINTERPOL(0.D0,0,0.D0,0)
   RETURN
END SUBROUTINE SHIFTPOT_VACUUMPOT
!###########################################################
!# FUNCTION: getscalefactor_vacuumpot 
!###########################################################
!> @brief 
!! Common get function. Gets the sum of atributes surfen and potmin
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getscalefactor_vacuumpot(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(IN):: this
   ! Run section
   getscalefactor_vacuumpot=this%surfen+this%potmin
   RETURN
END FUNCTION getscalefactor_vacuumpot
!###########################################################
!# FUNCTION: getderiv_vacuumpot 
!###########################################################
!> @brief 
!! Common get function. Gets value of the derivative of the potential at R
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getderiv_vacuumpot(this,r) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(IN)::this
   REAL(KIND=8),INTENT(IN) :: r
   ! Run section
   getderiv_vacuumpot=this%rpot%getderiv(r)
   RETURN
END FUNCTION getderiv_vacuumpot
!###########################################################
!# FUNCTION: getpot_vacuumpot 
!###########################################################
!> @brief 
!! Common get function. Gets value of the potential at R
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getpot_vacuumpot(this,r) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(IN)::this
   REAL(KIND=8),INTENT(IN) :: r
   ! Run section
   getpot_vacuumpot=this%rpot%getvalue(r)
   RETURN
END FUNCTION getpot_vacuumpot
!###########################################################
!# SUBROUTINE: PLOTDATA_VACUUMPOT 
!###########################################################
!> @brief
!! Plots a graph with values of the potential at grid points
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 26/Mar/2013
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE PLOTDATA_VACUUMPOT(this,filename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(IN)::this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Run section
   CALL this%rpot%PLOTDATA(filename)
   RETURN
END SUBROUTINE PLOTDATA_VACUUMPOT
!###########################################################
!# SUBROUTINE: PLOT_VACUUMPOT 
!###########################################################
!> @brief
!! Plots a graph once interpolation was done
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE PLOT_VACUUMPOT(this,npoints,filename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(IN)::this
   INTEGER(KIND=4),INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Run section
   CALL this%rpot%PLOT(npoints,filename)
   RETURN
END SUBROUTINE PLOT_VACUUMPOT
!###########################################################
!# SUBROUTINE: INITIALIZE_VACUUMPOT 
!###########################################################
!> @brief
!! Reads data from file, interpolates data, finds minimum
!! energy point
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 26/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INITIALIZE_VACUUMPOT(this,filename)
   ! Initial declarations   
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(INOUT) :: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   CHARACTER(LEN=22),PARAMETER :: routinename="INITIALIZE_VACUUMPOT: "
   ! Run section
   CALL this%READ(filename)
   CALL this%rpot%INTERPOL(0.D0,0,0.D0,0)
   CALL this%rpot%SET_MINIMUM()
   SELECT CASE(ALLOCATED(this%rpot%xmin))
      CASE(.TRUE.)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Minimums found: ",size(this%rpot%xmin))
         CALL VERBOSE_WRITE(routinename,this%rpot%xmin)
#endif
         SELECT CASE(size(this%rpot%xmin))
            CASE(1)
               !do nothing
               this%potmin = this%rpot%getvalue(this%rpot%xmin(1))
            CASE DEFAULT
               WRITE(0,*) "INITIALIZE_VACUUMPOT ERR: More than one minimum. Something is wrong"
               CALL EXIT(1)
         END SELECT
      CASE(.FALSE.)
        WRITE(0,*) "INITIALIZE_VACUUMPOT ERR: There's not any minimum. Something is wrong"
        CALL EXIT(1) 
   END SELECT
   RETURN
END SUBROUTINE INITIALIZE_VACUUMPOT
!###########################################################
!# SUBROUTINE: READ_VACUUMPOT 
!###########################################################
!> @brief
!! Reads data from file
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 26/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE READ_VACUUMPOT(this,filename)
   ! Initial declarations   
   USE UNITS_MOD
   USE MATHS_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(INOUT) :: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables
   CHARACTER(LEN=10) :: lenunits,enunits
   REAL(KIND=8) :: aux1,aux2
   TYPE(Energy) :: en
   TYPE(Length) :: len
   INTEGER(KIND=4) :: i ! counter
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: r,f
   ! Run section
   OPEN (111,FILE=filename,STATUS="old",ACTION="read")
   READ(111,*) ! dummy line
   READ(111,*) lenunits,enunits
   READ(111,*) aux1
   CALL en%READ(aux1,enunits)
   CALL en%TO_STD()
   this%surfen=en%getvalue()
   READ(111,*) this%n
   ALLOCATE(r(this%n))
   ALLOCATE(f(this%n))
   DO i = 1, this%n
      READ(111,*) aux1,aux2
      CALL len%READ(aux1,lenunits)
      CALL len%TO_STD()
      r(i)=len%getvalue()
      CALL en%READ(aux2,enunits)
      CALL en%TO_STD()
      f(i)=en%getvalue()
   END DO
   CALL ORDER(r,f) ! order R from low values to high values
   CALL this%rpot%READ(r,f)
   CLOSE(111)
   RETURN
END SUBROUTINE READ_VACUUMPOT
END MODULE EXTRAPOL_TO_VACUUM_MOD
