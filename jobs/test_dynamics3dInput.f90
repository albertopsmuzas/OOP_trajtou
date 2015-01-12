PROGRAM  TEST_DYNAMICS
use DEBUG_MOD, only: VERBOSE_WRITE, SET_VERBOSE_MODE
use DYNATOM_MOD, only: Dynatom
use SYSTEM_MOD
IMPLICIT NONE
TYPE(Dynatom) :: thisdynamics
REAL(KIND=4),DIMENSION(2) :: timearr
REAL(KIND=4) :: timer
CHARACTER(LEN=1024):: luafile
! RUN SECTION: GABBA, GABBA HEY! =========================================================================
!
! STEP 0: HELLO! & system specifications -----------------------------------------------------------------
WRITE(*,*) "******************************************************" 
WRITE(*,*) "************** TEST 3D DYNAMICS INPUT ****************"
WRITE(*,*) "******************************************************" 
! STEP 1: START UP DYNAMICS
SELECT CASE(command_argument_count())
   CASE(1)
      CALL GET_COMMAND_ARGUMENT(1,luafile)
   CASE DEFAULT
      WRITE(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      WRITE(0,*) "It is only needed one string, which is a config lua file"
      CALL EXIT(1)
END SELECT
CALL SET_VERBOSE_MODE(.true.)
CALL ETIME(timearr,timer)
CALL INITIALIZE_SYSTEM(trim(luafile))
CALL thisdynamics%INITIALIZE()
CALL ETIME(timearr,timer)
! STEP 4: PRINT TIME 
CALL VERBOSE_WRITE("****************** RUN TIME ***************************")
CALL VERBOSE_WRITE('',"User time: ",real(timearr(1),kind=8))
CALL VERBOSE_WRITE('',"System time: ",real(timearr(2),kind=8))
CALL VERBOSE_WRITE('',"Total time: ",real(timer,kind=8))
CALL VERBOSE_WRITE("*******************************************************")
END PROGRAM TEST_DYNAMICS
