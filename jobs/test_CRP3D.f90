PROGRAM CRP_TEST
! Initial declarations
USE DEBUG_MOD
USE CRP3D_MOD
USE UNITS_MOD
IMPLICIT NONE
! Variables
TYPE(CRP3D) :: crp_pes
CHARACTER(LEN=40) :: filename1,filename2,filename3
INTEGER(KIND=4) :: nsites, npairpots
REAL(KIND=8),DIMENSION(3) :: r
TYPE(Length) :: aux
REAL(KIND=8) :: l,l2
INTEGER(KIND=4) :: i ! counters
REAL(KIND=4) :: tiempo
REAL(KIND=4),DIMENSION(2) :: tiempoarr
! Run section===================================================================================0
CALL SET_DEBUG_MODE(.FALSE.)
CALL SET_VERBOSE_MODE(.TRUE.)
CALL ETIME(tiempoarr,tiempo)
! STEP 1: HELLO!
WRITE(*,*) "***************************************" 
WRITE(*,*) "CRP_PES program executed"
WRITE(*,*) "***************************************" 

! STEP 2: INITIALIZE CRP PES:
CALL crp_pes%READ("INcrp3d.inp")

! STEP 3: DO Z INTERPOLATION EXTRACTING VASINT AND SMOOTHING SITES
CALL crp_pes%INTERPOL()

! STEP 4: PLOT SOME GRAPHS
nsites=size(crp_pes%all_sites)
npairpots=size(crp_pes%all_pairpots)
r=(/0.D0,0.D0,3.34D0/)
l=5.44335612578
l2=12.D0
!CALL crp_pes%PLOT_XYMAP("xymap.dat",r,10,10,l,l)
WRITE(*,*) "popodolto"
CALL crp_pes%PLOT_DIRECTION1D("zscan.dat",500,45.D0,r(3),l2)
WRITE(*,*) "papadelta"
!CALL crp_pes%all_sites(3)%interz%PLOT_INTERPOL(50,"Site3graph")
CALL ETIME(tiempoarr,tiempo)
WRITE(*,*) "Total tile: ", tiempo
END PROGRAM CRP_TEST
