PROGRAM TEST_INICOND6D_INPUT
! Initial declarations
use DEBUG_MOD, only: VERBOSE_WRITE,SET_VERBOSE_MODE
use INITDIATOMIC_MOD, only: Initdiatomic
use CRP6D_MOD, only: CRP6D
use SYSTEM_MOD
IMPLICIT NONE
! Some objects
TYPE(Initdiatomic) :: inicondat
REAL(KIND=4),DIMENSION(2):: timearr
TYPE(CRP6D):: thispes
REAL(KIND=4):: timer
CHARACTER(LEN=1024):: luafile
! STEP 0: HELLO! & system specifications
CALL SET_VERBOSE_MODE(.TRUE.)
WRITE(*,*) "**************************************************************"
WRITE(*,*) "******************* TEST INICOND6D INPUT *********************"
WRITE(*,*) "**************************************************************"
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
CALL inicondat%INITIALIZE()
CALL inicondat%GENERATE_TRAJS(thispes)
CALL ETIME(timearr,timer)
! STEP 3: PRINT TIME 
CALL VERBOSE_WRITE("****************** RUN TIME ***************************")
CALL VERBOSE_WRITE('',"User time: ",real(timearr(1),kind=8))
CALL VERBOSE_WRITE('',"System time: ",real(timearr(2),kind=8))
CALL VERBOSE_WRITE('',"Total time: ",real(timer,kind=8))
CALL VERBOSE_WRITE("*******************************************************")
END PROGRAM TEST_INICOND6D_INPUT
