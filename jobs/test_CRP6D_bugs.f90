PROGRAM TEST_CRP6D
use SYSTEM_MOD
use CRP6D_MOD
use DEBUG_MOD
IMPLICIT NONE
! Variables
TYPE(CRP6D) :: thispes
TYPE(CRP6D) :: thisrawpes
REAL(KIND=8),DIMENSION(6) :: pos
CHARACTER(LEN=1024):: luaFile
REAL(KIND=4),DIMENSION(2):: timearr
REAL(KIND=4):: timer
! STEP 1: READ CRP6D INPUT FILES & INITIALIZE EVERYTHING
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
! datagroup 1
pos=(/0.0D0,0.D0,4.724315D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.1.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.1.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.1.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.1.dat")
pos=(/0.0D0,0.D0,2.834589D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.close.1.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.close.1.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.close.1.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.close.1.dat")

! datagroup 2
pos=(/0.D0,0.D0,5.669178D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.2.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.2.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.2.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.2.dat")
pos=(/0.D0,0.D0,3.779452D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.close.2.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.close.2.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.close.2.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.close.2.dat")

! datagroup 3
pos=(/0.D0,0.D0,4.724315D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.3.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.3.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.3.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.3.dat")
pos=(/0.D0,0.D0,2.834589D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.close.3.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.close.3.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.close.3.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.close.3.dat")

! datagroup 4
pos=(/0.D0,0.D0,5.669178D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.4.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.4.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.4.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.4.dat")
pos=(/0.D0,0.D0,3.779452D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.close.4.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.close.4.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.close.4.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.close.4.dat")

! datagroup 5
pos=(/2.721678D0,2.721678D0,2.834589D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.close.5.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.close.5.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.close.5.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.close.5.dat")

! datagroup 6
pos=(/2.721678D0,2.721678D0,3.779452D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.close.6.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.close.6.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.close.6.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.close.6.dat")

! datagroup 7
pos=(/2.721678D0,2.721678D0,5.669178D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.7.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.7.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.7.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.7.dat")
pos=(/2.721678D0,2.721678D0,2.834589D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.close.7.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.close.7.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.close.7.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.close.7.dat")

! datagroup 8
pos=(/2.721678D0,2.721678D0,4.724315D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.8.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.8.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.8.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.8.dat")
pos=(/2.721678D0,2.721678D0,3.779452D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.close.8.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.close.8.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.close.8.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.close.8.dat")

! datagroup 9
pos=(/2.721678D0,2.721678D0,5.669178D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.9.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.9.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.9.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.9.dat")
pos=(/1.814453D0,0.907225D0,2.834589D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.close.9.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.close.9.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.close.9.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.close.9.dat")

! datagroup 10
pos=(/1.814453D0,0.907225D0,4.724315D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.10.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.10.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.10.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.10.dat")
pos=(/1.814453D0,0.907225D0,3.779452D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.close.10.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.close.10.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.close.10.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.close.10.dat")

! datagroup 11
pos=(/1.814453D0,0.907225D0,5.669178D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.11.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.11.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.11.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.11.dat")
pos=(/1.814453D0,0.907225D0,2.834589D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.close.11.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.close.11.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.close.11.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.close.11.dat")

! datagroup 12
pos=(/1.814453D0,0.907225D0,4.724315D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.12.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.12.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.12.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.12.dat")
pos=(/1.814453D0,0.907225D0,3.779452D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.close.12.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.close.12.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.close.12.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.close.12.dat")
! datagroup 13
pos=(/0.D0,0.D0,4.724315D0,1.511781D0,0.D0,0.785398D0/)
CALL thispes%PLOT1D_THETA(1000,pos,"datagroup.13.dat")
CALL thispes%PLOT1D_THETA_SMOOTH(1000,pos,"datagroupsmooth.13.dat")
CALL thisrawpes%PLOT1D_THETA_SMOOTH(1000,pos,"datagroupraw.13.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_THETA(1000,pos,"datagroupatomic.13.dat")
pos=(/0.D0,0.D0,2.834589D0,1.511781D0,0.D0,0.785398D0/)
CALL thispes%PLOT1D_THETA(1000,pos,"datagroup.close.13.dat")
CALL thispes%PLOT1D_THETA_SMOOTH(1000,pos,"datagroupsmooth.close.13.dat")
CALL thisrawpes%PLOT1D_THETA_SMOOTH(1000,pos,"datagroupraw.close.13.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_THETA(1000,pos,"datagroupatomic.close.13.dat")
! datagroup 14
pos=(/2.721678D0,2.721678D0,4.724315D0,1.511781D0,0.D0,0.785398D0/)
CALL thispes%PLOT1D_THETA(1000,pos,"datagroup.14.dat")
CALL thispes%PLOT1D_THETA_SMOOTH(1000,pos,"datagroupsmooth.14.dat")
CALL thisrawpes%PLOT1D_THETA_SMOOTH(1000,pos,"datagroupraw.14.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_THETA(1000,pos,"datagroupatomic.14.dat")
pos=(/2.721678D0,2.721678D0,2.834589D0,1.511781D0,0.D0,0.785398D0/)
CALL thispes%PLOT1D_THETA(1000,pos,"datagroup.close.14.dat")
CALL thispes%PLOT1D_THETA_SMOOTH(1000,pos,"datagroupsmooth.close.14.dat")
CALL thisrawpes%PLOT1D_THETA_SMOOTH(1000,pos,"datagroupraw.close.14.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_THETA(1000,pos,"datagroupatomic.close.14.dat")
! datagroup 15
pos=(/1.360839D0,1.360839D0,5.669178D0,1.511781D0,1.570796D0,0D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.15.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.15.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.15.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.15.dat")
pos=(/1.360839D0,1.360839D0,3.779452D0,1.511781D0,1.570796D0,0D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.close.15.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.close.15.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.close.15.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.close.15.dat")
! datagroup 16
pos=(/1.360839D0,1.360839D0,5.669178D0,1.511781D0,0.523599D0,0D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.16.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.16.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.16.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.16.dat")
pos=(/1.360839D0,1.360839D0,3.779452D0,1.511781D0,0.523599D0,0D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.close.16.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.close.16.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.close.16.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.close.16.dat")
! datagroup 17
pos=(/1.360839D0,1.360839D0,5.669178D0,1.511781D0,0.D0,0.785398D0/)
CALL thispes%PLOT1D_THETA(1000,pos,"datagroup.17.dat")
CALL thispes%PLOT1D_THETA_SMOOTH(1000,pos,"datagroupsmooth.17.dat")
CALL thisrawpes%PLOT1D_THETA_SMOOTH(1000,pos,"datagroupraw.17.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_THETA(1000,pos,"datagroupatomic.17.dat")
pos=(/1.360839D0,1.360839D0,3.779452D0,1.511781D0,0.D0,0.785398D0/)
CALL thispes%PLOT1D_THETA(1000,pos,"datagroup.close.17.dat")
CALL thispes%PLOT1D_THETA_SMOOTH(1000,pos,"datagroupsmooth.close.17.dat")
CALL thisrawpes%PLOT1D_THETA_SMOOTH(1000,pos,"datagroupraw.close.17.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_THETA(1000,pos,"datagroupatomic.close.17.dat")
! datagroup 18
pos=(/2.721678D0,0.D0,5.669178D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.18.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.18.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.18.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.18.dat")
pos=(/2.721678D0,0.D0,3.779452D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.close.18.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.close.18.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.close.18.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.close.18.dat")
! datagroup 19
pos=(/2.721678D0,0.D0,5.669178D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.19.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.19.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.19.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.19.dat")
pos=(/2.721678D0,0.D0,3.779452D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.close.19.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.close.19.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.close.19.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,pos,"datagroupatomic.close.19.dat")
! datagroup 20
pos=(/2.721678D0,0.D0,5.669178D0,1.511781D0,0.D0,0.785398D0/)
CALL thispes%PLOT1D_THETA(1000,pos,"datagroup.20.dat")
CALL thispes%PLOT1D_THETA_SMOOTH(1000,pos,"datagroupsmooth.20.dat")
CALL thisrawpes%PLOT1D_THETA_SMOOTH(1000,pos,"datagroupraw.20.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_THETA(1000,pos,"datagroupatomic.20.dat")
pos=(/2.721678D0,0.D0,3.779452D0,1.511781D0,0.D0,0.785398D0/)
CALL thispes%PLOT1D_THETA(1000,pos,"datagroup.close.20.dat")
CALL thispes%PLOT1D_THETA_SMOOTH(1000,pos,"datagroupsmooth.close.20.dat")
CALL thisrawpes%PLOT1D_THETA_SMOOTH(1000,pos,"datagroupraw.close.20.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_THETA(1000,pos,"datagroupatomic.close.20.dat")

CALL EXIT(0)
END PROGRAM TEST_CRP6D