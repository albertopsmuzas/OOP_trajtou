PROGRAM TEST_CRP6D
! Initial declarations
use SYSTEM_MOD, only: INITIALIZE_SYSTEM
use CRP6D_MOD, only: CRP6D
use DEBUG_MOD, only: VERBOSE_WRITE
IMPLICIT NONE
! Variables
TYPE(CRP6D):: thispes
TYPE(CRP6D):: thisrawpes
CHARACTER(LEN=1024):: luaFile
REAL(KIND=4),DIMENSION(2):: timearr
REAL(KIND=4):: timer
REAL(KIND=8),DIMENSION(6) :: r
real(kind=8),parameter:: ucell=2.d0*5.44335612578d0
! STEP 1: READ CRP6D INPUT FILES
SELECT CASE(command_argument_count())
   CASE(1)
      CALL GET_COMMAND_ARGUMENT(1,luaFile)
   CASE DEFAULT
      WRITE(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      WRITE(0,*) "It is only needed one string, which is a config lua file"
      CALL EXIT(1)
END SELECT
CALL INITIALIZE_SYSTEM(trim(luaFile))
CALL VERBOSE_WRITE('##############################################')
CALL VERBOSE_WRITE('######### TEST SOME PRE-DEFINED CUTS #########')
CALL VERBOSE_WRITE('##############################################')
CALL thispes%INITIALIZE()
CALL thisrawpes%READ(trim(luaFile),'pes')
CALL thisrawpes%RAWINTERPOL()
!CALL thisrawpes%INTERPOL_NEW_RZGRID(25,50)
! corrugation extraction + addition
r=[0.d0,0.d0,13.6060281568d0,1.42d0,0.d0,0.d0]
CALL thispes%PLOT_XYMAP(r,150,150,ucell,ucell,"xycut_Z7.2.dat")
r=[0.d0,0.d0,13.0391103169d0,1.42d0,0.d0,0.d0]
CALL thispes%PLOT_XYMAP(r,150,150,ucell,ucell,"xycut_Z6.9.dat")
r=[0.d0,0.d0,11.3383567973d0,1.42d0,0.d0,0.d0]
CALL thispes%PLOT_XYMAP(r,150,150,ucell,ucell,"xycut_Z6.dat")
r=[0.d0,0.d0,9.44863066443d0,1.42d0,0.d0,0.d0]
CALL thispes%PLOT_XYMAP(r,150,150,ucell,ucell,"xycut_Z5.dat")
r=[0.d0,0.d0,7.55890453154d0,1.42d0,0.d0,0.d0]
CALL thispes%PLOT_XYMAP(r,150,150,ucell,ucell,"xycut_Z4.dat")
r=[0.d0,0.d0,5.66917839866d0,1.42d0,0.d0,0.d0]
CALL thispes%PLOT_XYMAP(r,150,150,ucell,ucell,"xycut_Z3.dat")
r=[0.d0,0.d0,3.77945226577d0,1.42d0,0.d0,0.d0]
CALL thispes%PLOT_XYMAP(r,150,150,ucell,ucell,"xycut_Z2.dat")
CALL thispes%atomiccrp(1)%PLOT_XYMAP(init_xyz=r(1:3),nxpoints=150,nypoints=150,Lx=ucell,Ly=ucell,filename="xycut_Z2.atom.dat")
r=[0.d0,0.d0,1.88972613289d0,1.42d0,0.d0,0.d0]
CALL thispes%PLOT_XYMAP(r,150,150,ucell,ucell,"xycut_Z1.dat")
! raw
r=[0.d0,0.d0,13.6060281568d0,1.42d0,0.d0,0.d0]
CALL thisrawpes%PLOT_XYMAP_SMOOTH(r,150,150,ucell,ucell,"xycut_Z7.2.raw.dat")
r=[0.d0,0.d0,13.0391103169d0,1.42d0,0.d0,0.d0]
CALL thisrawpes%PLOT_XYMAP_SMOOTH(r,150,150,ucell,ucell,"xycut_Z6.9.raw.dat")
r=[0.d0,0.d0,11.3383567973d0,1.42d0,0.d0,0.d0]
CALL thisrawpes%PLOT_XYMAP_SMOOTH(r,150,150,ucell,ucell,"xycut_Z6.raw.dat")
r=[0.d0,0.d0,9.44863066443d0,1.42d0,0.d0,0.d0]
CALL thisrawpes%PLOT_XYMAP_SMOOTH(r,150,150,ucell,ucell,"xycut_Z5.raw.dat")
r=[0.d0,0.d0,7.55890453154d0,1.42d0,0.d0,0.d0]
CALL thisrawpes%PLOT_XYMAP_SMOOTH(r,150,150,ucell,ucell,"xycut_Z4.raw.dat")
r=[0.d0,0.d0,5.66917839866d0,1.42d0,0.d0,0.d0]
CALL thisrawpes%PLOT_XYMAP_SMOOTH(r,150,150,ucell,ucell,"xycut_Z3.raw.dat")
r=[0.d0,0.d0,3.77945226577d0,1.42d0,0.d0,0.d0]
CALL thisrawpes%PLOT_XYMAP_SMOOTH(r,150,150,ucell,ucell,"xycut_Z2.raw.dat")
r=[0.d0,0.d0,1.88972613289d0,1.42d0,0.d0,0.d0]
CALL thisrawpes%PLOT_XYMAP_SMOOTH(r,150,150,ucell,ucell,"xycut_Z1.raw.dat")


CALL VERBOSE_WRITE("****************** RUN TIME ***************************")
CALL VERBOSE_WRITE('',"User time: ",real(timearr(1),kind=8))
CALL VERBOSE_WRITE('',"System time: ", real(timearr(2),kind=8))
CALL VERBOSE_WRITE('',"Total time: ",real(timer,kind=8))
CALL VERBOSE_WRITE("******************************************************")
END PROGRAM TEST_CRP6D
