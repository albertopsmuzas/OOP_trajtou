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
CALL thisrawpes%INTERPOL_NEW_RZGRID(25,50)
! STEP 2: DO GRAPHS

! Randomconf 1
r=(/4.43141D0,2.55113D0,3.77945D0,1.13384D0,1.91986D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,r,"randomconf.1.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,r,"randomconfsmooth.1.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,r,"randomconfraw.1.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,r,"randomconfatomic.1.dat")
r=(/2.5511300000000001D0,1.0119461257770961D0,3.77945D0,1.0D0,1.91986D0,0.D0/)
CALL thispes%PLOT1D_Z(1000,r,14.D0,"zcut1.dat")
! Randomconf 2
r=(/10.1006D0,6.29279D0,5.669180D0,1.51178D0,1.91986D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,r,"randomconf.2.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,r,"randomconfsmooth.2.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,r,"randomconfraw.2.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,r,"randomconfatomic.2.dat")
r=(/0.84943387422290395D0,0.78611225155419195D0,5.669180D0,1.5D0,1.91986D0,0.D0/)
CALL thispes%PLOT1D_Z(1000,r,14.D0,"zcut2.dat")
! Randomconf 3
r=(/10.1006D0,6.29279D0,5.669180D0,1.3606D0,0.D0,0.558505D0/)
CALL thispes%PLOT1D_THETA(1000,r,"randomconf.3.dat")
CALL thispes%PLOT1D_THETA_SMOOTH(1000,r,"randomconfsmooth.3.dat")
CALL thisrawpes%PLOT1D_THETA_SMOOTH(1000,r,"randomconfraw.3.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_THETA(1000,r,"randomconfatomic.3.dat")
r=(/0.84943387422290395D0,0.78611225155419195D0,5.669180D0,1.6D0,0.D0,2.58308765359D0/)
CALL thispes%PLOT1D_Z(1000,r,14.D0,"zcut3.dat")
! Randomconf 4
r=(/10.1006D0,6.29279D0,7.5589D0,1.3606D0,0.D0,0.558505D0/)
CALL thispes%PLOT1D_THETA(1000,r,"randomconf.4.dat")
CALL thispes%PLOT1D_THETA_SMOOTH(1000,r,"randomconfsmooth.4.dat")
CALL thisrawpes%PLOT1D_THETA_SMOOTH(1000,r,"randomconfraw.4.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_THETA(1000,r,"randomconfatomic.4.dat")
r=(/0.84943387422290395D0,0.78611225155419195D0,7.5589D0,1.7D0,0.D0,2.58308765359D0/)
CALL thispes%PLOT1D_Z(1000,r,14.D0,"zcut4.dat")

r=(/0.84943387422290395D0,0.78611225155419195D0,10.D0,1.3606D0,1.23D0,2.58308765359D0/)
CALL thispes%PLOT_XYMAP(r,100,100,8.D0,8.D0,"xycut10.dat")
r=(/0.84943387422290395D0,0.78611225155419195D0,8.D0,1.3606D0,1.23D0,2.58308765359D0/)
CALL thispes%PLOT_XYMAP(r,100,100,8.D0,8.D0,"xycut8.dat")
r=(/0.84943387422290395D0,0.78611225155419195D0,7.5D0,1.3606D0,1.23D0,2.58308765359D0/)
CALL thispes%PLOT_XYMAP(r,100,100,8.D0,8.D0,"xycut7.5.dat")
r=(/0.84943387422290395D0,0.78611225155419195D0,7.D0,1.3606D0,1.23D0,2.58308765359D0/)
CALL thispes%PLOT_XYMAP(r,100,100,8.D0,8.D0,"xycut7.dat")
r=(/0.84943387422290395D0,0.78611225155419195D0,5.D0,1.3606D0,1.23D0,2.58308765359D0/)
CALL thispes%PLOT_XYMAP(r,100,100,8.D0,8.D0,"xycut5.dat")
r=(/0.84943387422290395D0,0.78611225155419195D0,4.D0,1.3606D0,1.23D0,2.58308765359D0/)
CALL thispes%PLOT_XYMAP(r,100,100,8.D0,8.D0,"xycut4.dat")
CALL VERBOSE_WRITE("****************** RUN TIME ***************************")
CALL VERBOSE_WRITE('',"User time: ",real(timearr(1),kind=8))
CALL VERBOSE_WRITE('',"System time: ", real(timearr(2),kind=8))
CALL VERBOSE_WRITE('',"Total time: ",real(timer,kind=8))
CALL VERBOSE_WRITE("******************************************************")
END PROGRAM TEST_CRP6D
