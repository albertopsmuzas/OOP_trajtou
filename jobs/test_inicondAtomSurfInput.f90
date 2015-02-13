PROGRAM test_inicondAtomSurfInput
! Initial declarations
use SYSTEM_MOD
use INITATOMSURF_MOD, only: InitAtomSurf
use CRP3D_MOD, only: Crp3d
use DEBUG_MOD, only: VERBOSE_WRITE,SET_VERBOSE_MODE
IMPLICIT NONE
! Some objects
TYPE(InitAtomSurf) :: inicond
TYPE(Crp3d):: thispes
REAL(KIND=4),DIMENSION(2):: timeArr
REAL(KIND=4):: timer
CHARACTER(LEN=1024):: luafile
! STEP 0: HELLO! & system specifications
CALL SET_VERBOSE_MODE(.TRUE.)
call verbose_write("******************************************************************")
call verbose_write("*************** TEST INICOND ATOM SURF INPUT *********************")
call verbose_write("******************************************************************")
! STEP 1: READ PES, INTERPOLATION NOT NEEDED
SELECT CASE(command_argument_count())
   CASE(1)
      CALL GET_COMMAND_ARGUMENT(1,luafile)
   CASE DEFAULT
      WRITE(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      WRITE(0,*) "It is only needed one string, which is a config lua file"
      CALL EXIT(1)
END SELECT
! STEP 2: GENERATE NEW INITIAL CONDITIONS
CALL ETIME(timearr,timer)
CALL INITIALIZE_SYSTEM(trim(luafile))
CALL thispes%INITIALIZE()
CALL inicond%INITIALIZE()
CALL inicond%GENERATE_TRAJS(thispes)
CALL ETIME(timearr,timer)
! STEP 3: PRINT TIME 
CALL VERBOSE_WRITE("************************* RUN TIME *******************************")
CALL VERBOSE_WRITE('',"User time: ",real(timeArr(1),kind=8))
CALL VERBOSE_WRITE('',"System time: ",real(timeArr(2),kind=8))
CALL VERBOSE_WRITE('',"Total time: ",real(timer,kind=8))
CALL VERBOSE_WRITE("******************************************************************")
END PROGRAM test_inicondAtomSurfInput
