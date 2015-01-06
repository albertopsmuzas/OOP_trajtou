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
REAL(KIND=8),DIMENSION(3):: x,dvdx
REAL(KIND=8):: v
TYPE(CRP3D):: thispes
REAL(KIND=4):: timer
CHARACTER(LEN=1024):: luafile
CHARACTER(LEN=1024):: auxstring
! GABBA GABBA HEY! ===============================
!
! STEP 1: INITIALIZE SYSTEM VIA LUA CONFIG FILE
SELECT CASE(command_argument_count())
   CASE(4)
      CALL GET_COMMAND_ARGUMENT(1,luafile)
      CALL GET_COMMAND_ARGUMENT(2,auxstring)
      READ(auxstring,*) x(1)
      CALL GET_COMMAND_ARGUMENT(3,auxstring)
      READ(auxstring,*) x(2)
      CALL GET_COMMAND_ARGUMENT(4,auxstring)
      READ(auxstring,*) x(3)
   CASE DEFAULT
      WRITE(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      WRITE(0,*) "Four arguments needed: lua config file, 3 real numbers"
      CALL EXIT(1)
END SELECT
CALL ETIME(timearr,timer)
CALL INITIALIZE_SYSTEM(trim(luafile))
CALL thispes%INITIALIZE()
CALL thispes%GET_V_AND_DERIVS(x,v,dvdx)
CALL ETIME(timearr,timer)

   ! STEP 2: GET VALUES
CALL VERBOSE_WRITE('Potential calculated at:',x(:))
WRITE(*,*) "Potential (a.u.): ", v
WRITE(*,*) "Derivatives in auxiliar cartesian coord. (a.u.): ", dvdx(:)

! STEP 3: PRINT TIME 
CALL VERBOSE_WRITE("****************** RUN TIME ***************************")
CALL VERBOSE_WRITE('',"User time: ",real(timearr(1),kind=8))
CALL VERBOSE_WRITE('',"System time: ",real(timearr(2),kind=8))
CALL VERBOSE_WRITE('',"Total time: ",real(timer,kind=8))
CALL VERBOSE_WRITE("******************************************************")
END PROGRAM TEST_CRP3DINPUT
