!##################################################
! PROGRAM: DRAWTRAJ6D
!> @brief
!! Creates a XYZ file
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jul/2014
!> @version 1.0
!##################################################
PROGRAM DRAWTRAJ6D
   USE SYSTEM_MOD, only: INITIALIZE_SYSTEM
   USE DRAWTRAJ_MOD, only: Drawtraj
   USE DEBUG_MOD, only: VERBOSE_WRITE, SET_VERBOSE_MODE
! Initial declarations
IMPLICIT NONE
TYPE(Drawtraj):: this
CHARACTER(LEN=1024):: trajFile,outFile,auxString
INTEGER(KIND=4):: nOrder
REAL(KIND=4),DIMENSION(2):: timearr
REAL(KIND=4):: timer
! GABBA, GABBA HEY! =================================
SELECT CASE(command_argument_count())
   CASE(4)
      CALL GET_COMMAND_ARGUMENT(1,auxString)
      CALL INITIALIZE_SYSTEM(trim(auxString))
      CALL VERBOSE_WRITE("******************************************************")
      CALL VERBOSE_WRITE("**************** DRAW TRAJECTORY 6D ******************")
      CALL VERBOSE_WRITE("******************************************************")
      CALL GET_COMMAND_ARGUMENT(2,auxString)
      READ(auxString,*) nOrder
      CALL VERBOSE_WRITE('','Surface order: ',nOrder)
      CALL GET_COMMAND_ARGUMENT(3,trajFile)
      CALL VERBOSE_WRITE('','Trajectory file: ',trim(trajFile))
      CALL GET_COMMAND_ARGUMENT(4,outFile)
      CALL VERBOSE_WRITE('','Output file: ',trim(trajFile))
      CALL VERBOSE_WRITE('','Output type: XYZ')
      CALL ETIME(timearr,timer)
      CALL this%INITIALIZE("XYZ")
      CALL this%DRAW(nOrder,trim(trajFile),trim(outFile))
      CALL ETIME(timearr,timer)

   CASE DEFAULT
      WRITE(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      WRITE(0,*) "Four arguments needed: lua config file, order of surface (int), trajectory file and output file"
      CALL EXIT(1)
END SELECT
! STEP 3: PRINT TIME
CALL VERBOSE_WRITE("****************** RUN TIME ***************************")
CALL VERBOSE_WRITE('',"User time: ",real(timearr(1),kind=8))
CALL VERBOSE_WRITE('',"System time: ",real(timearr(2),kind=8))
CALL VERBOSE_WRITE('',"Total time: ",real(timer,kind=8))
CALL VERBOSE_WRITE("******************************************************")
END PROGRAM DRAWTRAJ6D
