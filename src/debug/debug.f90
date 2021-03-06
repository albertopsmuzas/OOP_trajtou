!#########################################################
! MODULE: DEBUG_MOD
!
!> @brief
!! Manages debug messages 
!
!> @details
!! - Provides 3 global variables: debugmode, verbosemode and silentmode.
!!   Silentmode=.true. redirects all output to some specific file.
!! - There are 2 different types of subroutines: those that 
!!   print when verbosemode=.true. and those that print when
!!   debugmode=.true.
!! - It is hardly recommended to use printing subroutines 
!!   inside #ifdef DEBUG #endif blocks. Otherwise, the 
!!   standard program will check verbosemode and debugmode too
!!   often.
!! - The philosophy under this module is to set 3 different
!!   levels of verbosity:
!!   . @b Not @b using @b this @b module: the program will only output
!!     in some specific files, avoiding unnecessary standard 
!!     printing or checking debugging variables.
!!   . @b Using @b this @b module, @b verbosemode=.true. : subroutines
!!     print some std output. Just to know which routines
!!     were executed, maybe some @b very useful information.
!!   . @b Using @b this @b module, @b debugmode=.true. : subroutines
!!     print whatever the developer thinks it is useful to
!!     follow possible bugs. This option forces verbosemode=.true.
!!     Not recomended for real calculations (e.g. run long dynamics)
!
!> @todo
!! - Use unrestricted polymorphism instead of interfaces. NOT
!!   all compilers have this feature (gfortran doesn't). This
!!   is why the INTERFACE structures are kept for the moment.
!! - Implement silentmode properly.
!##########################################################
MODULE DEBUG_MOD
! Initial declarations
IMPLICIT NONE
! Global variables
LOGICAL,PRIVATE :: debugmode=.FALSE.
LOGICAL,PRIVATE :: verbosemode=.FALSE.
LOGICAL,PRIVATE :: silentmode=.FALSE.
CHARACTER(LEN=20),PRIVATE :: outputdebug_msg
INTERFACE DEBUG_WRITE
   MODULE PROCEDURE DEBUG_WRITE_STRING , DEBUG_WRITE_STRING_REAL, &
      DEBUG_WRITE_STRING_INT, DEBUG_WRITE_STRING_REAL_INT, &
      DEBUG_WRITE_STRING_INT_REAL, DEBUG_WRITE_STRING_INT_INT, &
      DEBUG_WRITE_STRING_REAL_REAL, DEBUG_WRITE_REAL, &
      DEBUG_WRITE_INT, DEBUG_WRITE_STRING_STRING, DEBUG_WRITE_REAL_REAL,&
      DEBUG_WRITE_VECT, DEBUG_WRITE_VECTINT, DEBUG_WRITE_STRING_VECT, &
      DEBUG_WRITE_STRING_INTVECT, DEBUG_WRITE_MSG_MULTISTRING, DEBUG_WRITE_STRING_LOG1
END INTERFACE DEBUG_WRITE
INTERFACE VERBOSE_WRITE
   MODULE PROCEDURE VERBOSE_WRITE_STRING1, VERBOSE_WRITE_STRING2, &
      VERBOSE_WRITE_STRING_STRING1, VERBOSE_WRITE_STRING_REAL1, &
      VERBOSE_WRITE_STRING_INT1, VERBOSE_WRITE_REAL_REAL_REAL1, &
      VERBOSE_WRITE_STRING_LOG1, VERBOSE_WRITE_INT1, VERBOSE_WRITE_VECT, &
      VERBOSE_WRITE_VECTCHAR, VERBOSE_WRITE_VECTINT, VERBOSE_WRITE_STRING_VECT, &
      VERBOSE_WRITE_STRING_INTVECT, VERBOSE_WRITE_MSG_MULTISTRING
END INTERFACE VERBOSE_WRITE
CONTAINS
!###########################################################
!# SUBROUTINE: SET_SILENT_MODE
!###########################################################
!> @brief
!! Common set subroutine. Sets silent mode.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2004
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_SILENT_MODE(is_silent,filename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   LOGICAL,INTENT(IN) :: is_silent
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Run section
   silentmode=is_silent
   outputdebug_msg=filename
   RETURN
END SUBROUTINE SET_SILENT_MODE
!###########################################################
!# SUBROUTINE: DEBUG_SEPARATOR1
!###########################################################
!> @brief
!! Prints a line "----------------------------------" when debugmode is
!! activated
!-----------------------------------------------------------
SUBROUTINE DEBUG_SEPARATOR1()
   ! Initial declarations   
   IMPLICIT NONE
   ! Run section
   SELECT CASE(debugmode)
      CASE(.TRUE.)
         WRITE(*,*) "---------------------------------------------------------"
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE DEBUG_SEPARATOR1
!###########################################################
!# SUBROUTINE: VERBOSE_SEPARATOR1
!###########################################################
!> @brief
!! Prints a line "----------------------------------" when verbosemode is
!! activated
!-----------------------------------------------------------
SUBROUTINE VERBOSE_SEPARATOR1()
   ! Initial declarations   
   IMPLICIT NONE
   ! Run section
   SELECT CASE(verbosemode)
      CASE(.TRUE.)
         WRITE(*,*) "---------------------------------------------------------"
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE VERBOSE_SEPARATOR1
!###########################################################
!# SUBROUTINE: SET_VERBOSE_MODE 
!###########################################################
! - Sets verbosemode .true. or .false.
! - Cannot set verbosemode=.false. while debugmode=.true.
!-----------------------------------------------------------
SUBROUTINE SET_VERBOSE_MODE(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   LOGICAL, INTENT(IN) :: this
   ! Run section
   IF ((debugmode.EQV..true.).AND.(this.EQV..false.)) THEN
      WRITE(0,*) "SET_VERBOSE_MODE ERR: debugmode=.true. Cannot set verbosemode=.false."
      CALL EXIT(1)
   END IF
   verbosemode=this
   RETURN
END SUBROUTINE SET_VERBOSE_MODE
!###########################################################
!# SUBROUTINE: SET_DEBUG_MODE 
!###########################################################
! - Sets debugmode .true. or .false.
! - Interacts also with verbosemode
!-----------------------------------------------------------
SUBROUTINE SET_DEBUG_MODE(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   LOGICAL, INTENT(IN) :: this
   ! Run section
   IF (this.EQV..true.) THEN
      debugmode=.true.
      verbosemode=.true.
   ELSE
      debugmode=.false.
   END IF
   RETURN
END SUBROUTINE SET_DEBUG_MODE
!##############################################
! SUBROUTINE: VERBOSE_WRITE ################### <--- for interface
!##############################################
! - Prints only if VERBOSE = .TRUE.
!----------------------------------------------
SUBROUTINE VERBOSE_WRITE_MSG_MULTISTRING(routinename,msg,x)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*),INTENT(IN) :: msg, routinename
   CHARACTER(LEN=*),DIMENSION(:) :: x
	! RUN ---
	IF(verbosemode.EQV..true.) WRITE(*,*) routinename, msg, x
END SUBROUTINE VERBOSE_WRITE_MSG_MULTISTRING

SUBROUTINE VERBOSE_WRITE_STRING_VECT(routinename,msg,x)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*), INTENT(IN) :: msg, routinename
   REAL(KIND=8),DIMENSION(:) :: x
	! RUN ---
	IF(verbosemode.EQV..true.) WRITE(*,*) routinename, msg, x
END SUBROUTINE VERBOSE_WRITE_STRING_VECT

SUBROUTINE VERBOSE_WRITE_STRING_INTVECT(routinename,msg,x)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*), INTENT(IN) :: msg, routinename
   INTEGER(KIND=4),DIMENSION(:) :: x
	! RUN ---
	IF(verbosemode.EQV..true.) WRITE(*,*) routinename, msg, x
END SUBROUTINE VERBOSE_WRITE_STRING_INTVECT

SUBROUTINE VERBOSE_WRITE_STRING1(routinename,msg)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*), INTENT(IN) :: msg, routinename
	! RUN ---
	IF(verbosemode.EQV..true.) WRITE(*,*) routinename, msg
END SUBROUTINE VERBOSE_WRITE_STRING1

SUBROUTINE VERBOSE_WRITE_STRING2(string)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*), INTENT(IN) :: string
	! RUN ---
	IF(verbosemode.EQV..true.) WRITE(*,*) string
END SUBROUTINE VERBOSE_WRITE_STRING2

SUBROUTINE VERBOSE_WRITE_STRING_STRING1(routinename,msg1,msg2)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*), INTENT(IN) :: msg1,msg2,routinename
	! RUN ---
	IF(verbosemode.EQV..true.) WRITE(*,*) routinename, msg1, msg2
END SUBROUTINE VERBOSE_WRITE_STRING_STRING1

SUBROUTINE VERBOSE_WRITE_STRING_REAL1(routinename,msg,a)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*), INTENT(IN) :: msg, routinename
	REAL*8, INTENT(IN) :: a
	! RUN ---
	IF(verbosemode.EQV..true.) WRITE(*,*) routinename, msg, a
END SUBROUTINE VERBOSE_WRITE_STRING_REAL1

SUBROUTINE VERBOSE_WRITE_STRING_LOG1(routinename,msg,a)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*), INTENT(IN) :: msg, routinename
	LOGICAL, INTENT(IN) :: a
	! RUN ---
	IF(verbosemode.EQV..true.) WRITE(*,*) routinename, msg, a
END SUBROUTINE VERBOSE_WRITE_STRING_LOG1

SUBROUTINE VERBOSE_WRITE_STRING_INT1(routinename,msg,a)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*), INTENT(IN) :: msg, routinename
	INTEGER, INTENT(IN) :: a
	! RUN ---
	IF(verbosemode.EQV..true.) WRITE(*,*) routinename, msg, a
END SUBROUTINE VERBOSE_WRITE_STRING_INT1

SUBROUTINE VERBOSE_WRITE_INT1(routinename,a)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*), INTENT(IN) :: routinename
	INTEGER, INTENT(IN) :: a
	! RUN ---
	IF(verbosemode.EQV..true.) WRITE(*,*) routinename, a
END SUBROUTINE VERBOSE_WRITE_INT1

SUBROUTINE VERBOSE_WRITE_REAL_REAL_REAL1(routinename,a,b,c)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*), INTENT(IN) :: routinename
	REAL*8, INTENT(IN) :: a,b,c
	! RUN ---
	IF(verbosemode.EQV..true.) WRITE(*,*) routinename, a,b,c
END SUBROUTINE VERBOSE_WRITE_REAL_REAL_REAL1

SUBROUTINE VERBOSE_WRITE_VECT(routinename,a)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*),INTENT(IN) :: routinename
	REAL*8,DIMENSION(:),INTENT(IN) :: a
	! RUN ---
	IF(verbosemode.EQV..true.) WRITE(*,*) routinename,a(:)
END SUBROUTINE VERBOSE_WRITE_VECT

SUBROUTINE VERBOSE_WRITE_VECTCHAR(routinename,a)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*),INTENT(IN) :: routinename
	CHARACTER,DIMENSION(:),INTENT(IN) :: a
	! RUN ---
	IF(verbosemode.EQV..true.) WRITE(*,*) routinename,a(:)
END SUBROUTINE VERBOSE_WRITE_VECTCHAR

SUBROUTINE VERBOSE_WRITE_VECTINT(routinename,a)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*),INTENT(IN) :: routinename
	INTEGER(KIND=4),DIMENSION(:),INTENT(IN) :: a
	! RUN ---
	IF(verbosemode.EQV..true.) WRITE(*,*) routinename,a(:)
END SUBROUTINE VERBOSE_WRITE_VECTINT
!##############################################
! SUBROUTINE: DEBUG_WRITE ##################### <--- for interface
!##############################################
! - Prints only if DEBUG = .TRUE.
!----------------------------------------------
SUBROUTINE DEBUG_WRITE_STRING_VECT(routinename,msg,x)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*), INTENT(IN) :: msg, routinename
   REAL(KIND=8),DIMENSION(:) :: x
	! RUN ---
	IF(debugmode.EQV..true.) WRITE(*,*) routinename, msg, x
END SUBROUTINE DEBUG_WRITE_STRING_VECT

SUBROUTINE DEBUG_WRITE_STRING_INTVECT(routinename,msg,x)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*), INTENT(IN) :: msg, routinename
   INTEGER(KIND=4),DIMENSION(:) :: x
	! RUN ---
	IF(debugmode.EQV..true.) WRITE(*,*) routinename, msg, x
END SUBROUTINE DEBUG_WRITE_STRING_INTVECT

SUBROUTINE DEBUG_WRITE_VECTINT(routinename,a)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*),INTENT(IN) :: routinename
	INTEGER(KIND=4),DIMENSION(:),INTENT(IN) :: a
	! RUN ---
	IF(debugmode.EQV..true.) WRITE(*,*) routinename,a(:)
END SUBROUTINE DEBUG_WRITE_VECTINT

SUBROUTINE DEBUG_WRITE_REAL(subroutinename,a)
	IMPLICIT NONE
	! I/O variables
   CHARACTER(LEN=*),INTENT(IN) :: subroutinename
	REAL*8, INTENT(IN) :: a
	! RUN ---
	IF(debugmode.EQV..true.) WRITE(*,*) subroutinename, a
	RETURN
END SUBROUTINE DEBUG_WRITE_REAL

SUBROUTINE DEBUG_WRITE_INT(subroutinename,a)
	IMPLICIT NONE
	! I/O variables
   CHARACTER(LEN=*),INTENT(IN) :: subroutinename
	INTEGER, INTENT(IN) :: a
	! RUN ---
	IF(debugmode.EQV..true.) WRITE(*,*) subroutinename, a
	RETURN
END SUBROUTINE DEBUG_WRITE_INT

SUBROUTINE DEBUG_WRITE_STRING(subroutinename,string)
	IMPLICIT NONE
	! I/O variables
   CHARACTER(LEN=*),INTENT(IN) :: subroutinename
	CHARACTER(LEN=*), INTENT(IN) :: string
	! RUN ---
	IF(debugmode.EQV..true.) WRITE(*,*) subroutinename, string
	RETURN
END SUBROUTINE DEBUG_WRITE_STRING

SUBROUTINE DEBUG_WRITE_STRING_STRING(subroutinename,string1,string2)
	IMPLICIT NONE
	! I/O variables
   CHARACTER(LEN=*),INTENT(IN) :: subroutinename
	CHARACTER(LEN=*), INTENT(IN) :: string1, string2
	! RUN ---
	IF(debugmode.EQV..true.) WRITE(*,*) subroutinename, string1, string2
	RETURN
END SUBROUTINE DEBUG_WRITE_STRING_STRING

SUBROUTINE DEBUG_WRITE_STRING_REAL(subroutinename,string,a)
	IMPLICIT NONE
	! I/O variables
   CHARACTER(LEN=*),INTENT(IN) :: subroutinename
	CHARACTER(LEN=*), INTENT(IN) :: string
	REAL*8, INTENT(IN) :: a
	! RUN ---
	IF(debugmode.EQV..true.) WRITE(*,*) subroutinename, string, a
	RETURN
END SUBROUTINE DEBUG_WRITE_STRING_REAL

SUBROUTINE DEBUG_WRITE_STRING_REAL_REAL(subroutinename,string,a,b)
	IMPLICIT NONE
	! I/O variables
   CHARACTER(LEN=*),INTENT(IN) :: subroutinename
	CHARACTER(LEN=*), INTENT(IN) :: string
	REAL*8, INTENT(IN) :: a,b
	! RUN ---
	IF(debugmode.EQV..true.) WRITE(*,*) subroutinename, string, a, b
	RETURN
END SUBROUTINE DEBUG_WRITE_STRING_REAL_REAL

SUBROUTINE DEBUG_WRITE_REAL_REAL(subroutinename,a,b)
	IMPLICIT NONE
	! I/O variables
   CHARACTER(LEN=*),INTENT(IN) :: subroutinename
	REAL*8, INTENT(IN) :: a,b
	! RUN ---
	IF(debugmode.EQV..true.) WRITE(*,*) subroutinename, a, b
	RETURN
END SUBROUTINE DEBUG_WRITE_REAL_REAL

SUBROUTINE DEBUG_WRITE_STRING_INT(subroutinename,string,a)
	IMPLICIT NONE
	! I/O variables
   CHARACTER(LEN=*),INTENT(IN) :: subroutinename
	CHARACTER(LEN=*), INTENT(IN) :: string
	INTEGER, INTENT(IN) :: a
	! RUN ---
	IF(debugmode.EQV..true.) WRITE(*,*) subroutinename, string, a
	RETURN
END SUBROUTINE DEBUG_WRITE_STRING_INT

SUBROUTINE DEBUG_WRITE_STRING_INT_INT(subroutinename,string,a,b)
	IMPLICIT NONE
	! I/O variables
   CHARACTER(LEN=*),INTENT(IN) :: subroutinename
	CHARACTER(LEN=*), INTENT(IN) :: string
	INTEGER, INTENT(IN) :: a,b
	! RUN ---
	IF(debugmode.EQV..true.) WRITE(*,*) subroutinename, string, a, b
	RETURN
END SUBROUTINE DEBUG_WRITE_STRING_INT_INT

SUBROUTINE DEBUG_WRITE_STRING_REAL_INT(subroutinename,string,a,b)
	IMPLICIT NONE
	! I/O variables
   CHARACTER(LEN=*),INTENT(IN) :: subroutinename
	CHARACTER(LEN=*), INTENT(IN) :: string
	REAL*8, INTENT(IN) :: a
	INTEGER, INTENT(IN) :: b
	! RUN ---
	IF(debugmode.EQV..true.) WRITE(*,*) subroutinename, string, a, b
	RETURN
END SUBROUTINE DEBUG_WRITE_STRING_REAL_INT

SUBROUTINE DEBUG_WRITE_STRING_INT_REAL(subroutinename,string,a,b)
	IMPLICIT NONE
	! I/O variables
   CHARACTER(LEN=*),INTENT(IN) :: subroutinename
	CHARACTER(LEN=*), INTENT(IN) :: string
	REAL*8, INTENT(IN) :: b
	INTEGER, INTENT(IN) :: a
	! RUN ---
	IF(debugmode.EQV..true.) WRITE(*,*) subroutinename, string, a, b
	RETURN
END SUBROUTINE DEBUG_WRITE_STRING_INT_REAL

SUBROUTINE DEBUG_WRITE_VECT(subroutinename,a)
	IMPLICIT NONE
	! I/O variables
   CHARACTER(LEN=*),INTENT(IN) :: subroutinename
	REAL(KIND=8),DIMENSION(:),INTENT(IN) :: a
	! RUN ---
	IF(debugmode.EQV..true.) WRITE(*,*) subroutinename,a(:)
	RETURN
END SUBROUTINE DEBUG_WRITE_VECT

SUBROUTINE DEBUG_WRITE_MSG_MULTISTRING(routinename,msg,x)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*),INTENT(IN) :: msg, routinename
   CHARACTER(LEN=*),DIMENSION(:) :: x
	! RUN ---
	IF(debugMode.EQV..true.) WRITE(*,*) routinename, msg, x
END SUBROUTINE DEBUG_WRITE_MSG_MULTISTRING

SUBROUTINE DEBUG_WRITE_STRING_LOG1(routinename,msg,a)
	IMPLICIT NONE
	! I/O variables
	CHARACTER(LEN=*), INTENT(IN) :: msg, routinename
	LOGICAL, INTENT(IN) :: a
	! RUN ---
	IF(debugMode.EQV..true.) WRITE(*,*) routinename, msg, a
END SUBROUTINE DEBUG_WRITE_STRING_LOG1
!
!###################################################
! FUNCTION get_debugmode
!> @brief
!! Usual get function
!---------------------------------------------------
LOGICAL FUNCTION get_debugmode()
   IMPLICIT NONE
   get_debugmode=debugmode
END FUNCTION get_debugmode
!###################################################
! FUNCTION get_verbosemode
!> @brief
!! Usual get function
!---------------------------------------------------
LOGICAL FUNCTION get_verbosemode()
   IMPLICIT NONE
   get_verbosemode=verbosemode
END FUNCTION get_verbosemode
END MODULE DEBUG_MOD
