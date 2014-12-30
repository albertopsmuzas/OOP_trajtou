PROGRAM  TEST_DYNAMICS
USE DEBUG_MOD
USE DYNATOM_MOD
IMPLICIT NONE
TYPE(Dynatom) :: this
REAL(KIND=4),DIMENSION(2) :: timearr
REAL(KIND=4) :: timer
CHARACTER(LEN=1024):: luafile
! RUN SECTION: GABBA, GABBA HEY! =========================================================================
!
! STEP 0: HELLO! & system specifications -----------------------------------------------------------------
WRITE(*,*) "******************************************************" 
WRITE(*,*) "****************** 3D DYNAMICS ***********************"
WRITE(*,*) "******************************************************" 
CALL SET_VERBOSE_MODE(.TRUE.)
! STEP 1: START UP DYNAMICS
SELECT CASE(command_argument_count())
   CASE(1)
      CALL GET_COMMAND_ARGUMENT(1,luafile)
   CASE DEFAULT
      WRITE(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      WRITE(0,*) "It is only needed one string, which is a config lua file"
      CALL EXIT(1)
END SELECT
WRITE(*,*) "Input file: ", trim(luafile)
CALL this%INITIALIZE(trim(luafile))
! STEP 3: RUN DYNAMICS
CALL ETIME(timearr,timer)
CALL this%RUN()
CALL ETIME(timearr,timer)
! STEP 4: PRINT TIME 
WRITE(*,*) "****************** RUN TIME ***************************"
WRITE(*,*) "User time: ", timearr(1)
WRITE(*,*) "System time: ", timearr(2)
WRITE(*,*) "Total time: ",timer
WRITE(*,*) "******************************************************" 
END PROGRAM TEST_DYNAMICS
