!##################################################
! PROGRAM: TEST_SYSTEMINPUT
!> @brief
!! Test if system part of Lua conf file is read correctly
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!##################################################
PROGRAM TEST_SYSTEMINPUT
! Initial declarations
USE SYSTEM_MOD
USE DEBUG_MOD
! use some modules?
IMPLICIT NONE
! Variables
REAL(KIND=4),DIMENSION(2):: timearr
REAL(KIND=4):: timer
CHARACTER(LEN=1024):: luafile
! GABBA GABBA HEY! ===============================
!
! STEP 0: HELLO! & system specifications
WRITE(*,*) "******************************************************" 
WRITE(*,*) "************** TEST SYSTEM INPUT *********************"
WRITE(*,*) "******************************************************" 
CALL SET_VERBOSE_MODE(.TRUE.)
! STEP 1: INITIALIZE SYSTEM VIA LUA CONFIG FILE
SELECT CASE(command_argument_count())
   CASE(1)
      CALL GET_COMMAND_ARGUMENT(1,luafile)
   CASE DEFAULT
      WRITE(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      WRITE(0,*) "It is only needed one string, which is a config lua file"
      CALL EXIT(1)
END SELECT
CALL ETIME(timearr,timer)
CALL INITIALIZE_SYSTEM(trim(luafile))
CALL ETIME(timearr,timer)
! STEP 2: PRINT TIME 
CALL VERBOSE_WRITE("****************** RUN TIME ***************************")
CALL VERBOSE_WRITE('',"User time: ",real(timearr(1),kind=8))
CALL VERBOSE_WRITE('',"System time: ",real(timearr(2),kind=8))
CALL VERBOSE_WRITE('',"Total time: ",real(timer,kind=8))
CALL VERBOSE_WRITE("*******************************************************")
END PROGRAM TEST_SYSTEMINPUT
