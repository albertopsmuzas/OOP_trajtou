PROGRAM ANALYSIS
! Initial declarations
use SYSTEM_MOD, only: INITIALIZE_SYSTEM
use DEBUG_MOD, only: VERBOSE_WRITE
use DIFFRACTIONCRP3D_MOD, only: Allowed_peaksCRP3D
IMPLICIT NONE
! Variables
TYPE(Allowed_peaksCRP3D):: anlyDiff
CHARACTER(LEN=1024):: auxString
REAL(KIND=4),DIMENSION(2):: timearr
REAL(KIND=4):: timer
! HEY HO, LET'S GO --------------------------------------
! Get input
SELECT CASE(command_argument_count())
   CASE(1)
      CALL GET_COMMAND_ARGUMENT(1,auxString)
      CALL INITIALIZE_SYSTEM(trim(auxString))
      CALL VERBOSE_WRITE("******************************************************")
      CALL VERBOSE_WRITE("************** DIFFRACTION ANALYSIS 6D ***************")
      CALL VERBOSE_WRITE("******************************************************")
      CALL ETIME(timearr,timer)
      CALL anlyDiff%INITIALIZE()
      CALL anlyDiff%SETUP()
      CALL anlyDiff%ASSIGN_PEAKS()
      CALL anlyDiff%PRINT_LABMOMENTA_AND_ANGLES()
      CALL ETIME(timearr,timer)

   CASE DEFAULT
      WRITE(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      WRITE(0,*) "ONE argument needed: lua config file"
      CALL EXIT(1)
END SELECT
! Time
CALL VERBOSE_WRITE("****************** RUN TIME ***************************")
CALL VERBOSE_WRITE('',"User time: ",real(timearr(1),kind=8))
CALL VERBOSE_WRITE('',"System time: ",real(timearr(2),kind=8))
CALL VERBOSE_WRITE('',"Total time: ",real(timer,kind=8))
CALL VERBOSE_WRITE("******************************************************")
END PROGRAM ANALYSIS
