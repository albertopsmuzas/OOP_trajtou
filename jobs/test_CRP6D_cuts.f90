PROGRAM TEST_CRP6D
! Initial declarations
USE CRP6D_MOD
USE DEBUG_MOD
IMPLICIT NONE
! Variables
TYPE(CRP6D) :: thispes
TYPE(CRP6D) :: thisrawpes
REAL(KIND=8),DIMENSION(6) :: r
CALL SET_DEBUG_MODE(.FALSE.)
CALL SET_VERBOSE_MODE(.FALSE.)
! STEP 1: READ CRP6D INPUT FILES
CALL thispes%READ("INcrp6d.inp")
CALL thisrawpes%READ("INcrp6d.inp")
CALL thispes%INTERPOL()
CALL thispes%INTERPOL_NEW_RZGRID(50,100)
CALL thisrawpes%RAWINTERPOL()
CALL thisrawpes%INTERPOL_NEW_RZGRID(50,100)
! STEP 2: DO GRAPHS

! Randomconf 1
r=(/4.43141D0,2.55113D0,3.77945D0,1.13384D0,1.91986D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,r,"randomconf.1.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,r,"randomconfsmooth.1.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,r,"randomconfraw.1.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,r,"randomconfatomic.1.dat")
! Randomconf 2
r=(/10.1006D0,6.29279D0,5.669180D0,1.51178D0,1.91986D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,r,"randomconf.2.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,r,"randomconfsmooth.2.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,r,"randomconfraw.2.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_PHI(1000,r,"randomconfatomic.2.dat")

! Randomconf 3
r=(/10.1006D0,6.29279D0,5.669180D0,1.3606D0,0.D0,0.558505D0/)
CALL thispes%PLOT1D_THETA(1000,r,"randomconf.3.dat")
CALL thispes%PLOT1D_THETA_SMOOTH(1000,r,"randomconfsmooth.3.dat")
CALL thisrawpes%PLOT1D_THETA_SMOOTH(1000,r,"randomconfraw.3.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_THETA(1000,r,"randomconfatomic.3.dat")

! Randomconf 4
r=(/10.1006D0,6.29279D0,7.5589D0,1.3606D0,0.D0,0.558505D0/)
CALL thispes%PLOT1D_THETA(1000,r,"randomconf.4.dat")
CALL thispes%PLOT1D_THETA_SMOOTH(1000,r,"randomconfsmooth.4.dat")
CALL thisrawpes%PLOT1D_THETA_SMOOTH(1000,r,"randomconfraw.4.dat")
CALL thispes%PLOT1D_ATOMIC_INTERAC_THETA(1000,r,"randomconfatomic.4.dat")

CALL EXIT(0)
END PROGRAM TEST_CRP6D
