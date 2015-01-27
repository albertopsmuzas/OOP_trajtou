!#########################################################
! MODULE: EXTRAPOL_TO_VACUUM_MOD
!> @brief
!! Provides tools to extrapolate the potential in the vacuum,
!! which is only dependent on the distance between both atoms
!##########################################################
MODULE EXTRAPOL_TO_VACUUM_MOD
! Initial declarations
use CUBICSPLINES_MOD, only: Csplines
use MATHS_MOD
use UNITS_MOD, only: Energy,Length
#ifdef DEBUG
use DEBUG_MOD, only: VERBOSE_WRITE, DEBUG_WRITE
#endif
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: Vacuumpot
!> @brief
!! Class that stores all information that can be extracted from the
!! monodimensional vacuum potential
!----------------------------------------------------------------
TYPE :: Vacuumpot
   PRIVATE
   INTEGER(KIND=4),PUBLIC :: n
   TYPE(Csplines),PUBLIC :: rpot
   REAL(KIND=8) :: surfen
   REAL(KIND=8) :: potmin
   REAL(KIND=8) :: req
   REAL(KIND=8),DIMENSION(2),PUBLIC :: root

   CONTAINS
      ! Initialization block
      PROCEDURE,PUBLIC :: INITIALIZE => INITIALIZE_VACUUMPOT
      PROCEDURE,PUBLIC :: INITIALIZE_DIRECT => INITIALIZE_DIRECT_VACUUMPLOT
      PROCEDURE,PUBLIC :: READ => READ_VACUUMPOT
      ! Set block
      PROCEDURE,PUBLIC :: SET_ROOTS => SET_ROOTS_VACUUMPOT
      ! Get block
      PROCEDURE,PUBLIC :: getpot => getpot_vacuumpot
      PROCEDURE,PUBLIC :: getderiv => getderiv_vacuumpot
      PROCEDURE,PUBLIC :: getscalefactor => getscalefactor_vacuumpot
      PROCEDURE,PUBLIC :: getreq => getreq_vacuumpot 
      ! Tools block
      PROCEDURE,PUBLIC :: SHIFTPOT => SHIFTPOT_VACUUMPOT
      PROCEDURE,PUBLIC :: SHIFTPOT_UNDO => SHIFTPOT_UNDO_VACUUMPOT
      ! Enquire block
      PROCEDURE,PUBLIC :: is_allowed => is_allowed_VACUUMPOT
      ! Plot tools block
      PROCEDURE,PUBLIC :: PLOT => PLOT_VACUUMPOT
      PROCEDURE,PUBLIC :: PLOTDATA => PLOTDATA_VACUUMPOT
END TYPE Vacuumpot
!/////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_ROOTS_VACUUMPOT
!###########################################################
!> @brief
!! Calculates roots of the function if any
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_ROOTS_VACUUMPOT(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(INOUT):: this
   ! Local variables
   CHARACTER(LEN=21) :: routinename="SET_ROOTS_VACUUMPOT: "
   ! Run section
   CALL this%rpot%SET_XROOT()
   SELECT CASE(allocated(this%rpot%xroot))
      CASE(.TRUE.)
         SELECT CASE(size(this%rpot%xroot))
            CASE(1)
               WRITE(0,*) "SET_ROOTS_VACUUMPOT ERR: we are in req. Bad. Pure classic not implemented"
               CALL EXIT(1)
            CASE(2)
               this%root=this%rpot%xroot
#ifdef DEBUG
               CALL VERBOSE_WRITE(routinename,"Found roots at x: ",this%root)
#endif
            CASE DEFAULT
               WRITE(0,*) "SET_ROOTS_VACUUMPOT ERR: too many roots. Check the vacuumpot"
               CALl EXIT(1)
         END SELECT
         ! body
      CASE(.FALSE.)
         WRITE(0,*) "SET_ROOTS_VACUUMPOT ERR: there aren't roots"
         CALL EXIT(1)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE SET_ROOTS_VACUUMPOT
!###########################################################
!# FUNCTION: is_allowed_VACUUMPOT 
!###########################################################
!> @brief 
!! Checks if the potential can be calculated in @b X
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
LOGICAL FUNCTION is_allowed_VACUUMPOT(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   ! Local variables
   ! Run section
   SELECT CASE(x< this%rpot%x(1) .OR. x > this%rpot%x(this%n))
      CASE(.TRUE.)
         is_allowed_VACUUMPOT=.FALSE.
      CASE(.FALSE.)
         is_allowed_VACUUMPOT=.TRUE.
   END SELECT
   RETURN
END FUNCTION is_allowed_VACUUMPOT
!###########################################################
!# SUBROUTINE: READ_DIRECT_VACUUMPLOT 
!###########################################################
!> @brief
!! Reads information directly from array
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INITIALIZE_DIRECT_VACUUMPLOT(this,x,f,shift)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(OUT):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: f
   REAL(KIND=8),INTENT(IN):: shift
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: aux1,aux2
   CHARACTER(LEN=29),PARAMETER :: routinename="INITIALIZE_DIRECT_VACUUMPOT: "
   ! Run section
   this%n=size(x)
   ALLOCATE(aux1(size(x)))
   ALLOCATE(aux2(size(f)))
   aux1=x
   aux2=f
   CALL ORDER(aux1,aux2)
   aux2=aux2-shift
   CALL this%rpot%READ(aux1,aux2)
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
               this%req=this%rpot%xmin(1)
               this%potmin = this%rpot%getvalue(this%rpot%xmin(1))
#ifdef DEBUG
                CALL VERBOSE_WRITE(routinename,'Potmin: ',this%potmin)
#endif
            CASE DEFAULT
               WRITE(0,*) "INITIALIZE_DIRECT_VACUUMPOT ERR: More than one minimum. Something is wrong"
               CALL EXIT(1)
         END SELECT
      CASE(.FALSE.)
        WRITE(0,*) "INITIALIZE_DIRECT_VACUUMPOT ERR: There's not any minimum. Something is wrong"
        CALL EXIT(1) 
   END SELECT
   RETURN
END SUBROUTINE INITIALIZE_DIRECT_VACUUMPLOT
!###########################################################
!# FUNCTION: getreq_vacuumpot 
!###########################################################
!> @brief 
!! Common get function. Gets "req" atribute
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getreq_vacuumpot(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(IN):: this
   ! Run section
   getreq_vacuumpot=this%req
   RETURN
END FUNCTION getreq_vacuumpot  
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
!# SUBROUTINE: SHIFTPOT_UNDO_VACUUMPOT 
!###########################################################
!> @brief
!! Shift the entire potential so that the equilibrium geomtry has
!! energy 0
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 27/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SHIFTPOT_UNDO_VACUUMPOT(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Vacuumpot),INTENT(INOUT):: this
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   ! Run section
   DO i = 1, this%n
      this%rpot%f(i)=this%rpot%f(i)+this%potmin
   END DO
   CALL this%rpot%REINTERPOL(0.D0,0,0.D0,0)
   RETURN
END SUBROUTINE SHIFTPOT_UNDO_VACUUMPOT
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
               this%req=this%rpot%xmin(1)
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
