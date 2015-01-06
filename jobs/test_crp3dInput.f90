!##################################################
! PROGRAM: TEST_SYSTEMINPUT
!> @brief
!! Test if system part of Lua conf file is read correctly
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!##################################################
PROGRAM TEST_CRP3DINPUT
! Initial declarations
USE DEBUG_MOD
USE CRP3D_MOD
USE SYSTEM_MOD
! use some modules?
IMPLICIT NONE
! Variables
REAL(KIND=4),DIMENSION(2):: timearr
TYPE(CRP3D):: thispes
REAL(KIND=4):: timer
CHARACTER(LEN=1024):: luafile
! GABBA GABBA HEY! ===============================
!
! STEP 0: HELLO! & system specifications
WRITE(*,*) "******************************************************" 
WRITE(*,*) "************** TEST CRP3D INPUT **********************"
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
CALL thispes%INITIALIZE()
CALL ETIME(timearr,timer)
! STEP 2: PRINT TIME 
CALL VERBOSE_WRITE("****************** RUN TIME ***************************")
CALL VERBOSE_WRITE('',"User time: ",real(timearr(1),kind=8))
CALL VERBOSE_WRITE('',"System time: ",real(timearr(2),kind=8))
CALL VERBOSE_WRITE('',"Total time: ",real(timer,kind=8))
CALL VERBOSE_WRITE("*******************************************************")
END PROGRAM TEST_CRP3DINPUT
