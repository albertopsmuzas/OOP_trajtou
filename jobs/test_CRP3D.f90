PROGRAM CRP_TEST
! Initial declarations
USE DEBUG_MOD
USE CRP3D_MOD
USE UNITS_MOD
IMPLICIT NONE
! Variables
TYPE(CRP3D) :: thispes, thisrawpes
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
WRITE(*,*) "thispes program executed"
WRITE(*,*) "***************************************" 

! STEP 2: INITIALIZE CRP PES:
CALL thispes%READ("INcrp3d.inp")
CALL thisrawpes%READ("INcrp3d.inp")

! STEP 3: DO Z INTERPOLATION EXTRACTING VASINT AND SMOOTHING SITES
CALL thispes%INTERPOL()
CALL thisrawpes%RAWINTERPOL()
! STEP 4: PLOT SOME GRAPHS
nsites=size(thispes%all_sites)
npairpots=size(thispes%all_pairpots)
l=5.44335612578D0
l2=12.D0
r=(/0.D0,0.D0,3.D0/)
CALL thispes%PLOT_XYMAP("xymap.3.dat",r,200,200,l,l)
CALL thispes%PLOT_XYMAP_CORRECTION("xymapcorrection.3.dat",r,200,200,l,l)
CALL thisrawpes%PLOT_XYMAP_SMOOTH("xymapsmooth.3.dat",r,200,200,l,l)
CALL thispes%PLOT_DIRECTION1D("zscan.3.dat",500,45.D0,r(3),l2)
CALL thispes%PLOT_DIRECTION1D_CORRECTION("zscancorrection.3.dat",500,45.D0,r(3),l2)
CALL thisrawpes%PLOT_DIRECTION1D_SMOOTH("zscansmooth.3.dat",500,45.D0,r(3),l2)
r=(/0.D0,0.D0,5.D0/)
CALL thispes%PLOT_XYMAP("xymap.5.dat",r,200,200,l,l)
CALL thispes%PLOT_XYMAP_CORRECTION("xymapcorrection.5.dat",r,200,200,l,l)
CALL thisrawpes%PLOT_XYMAP_SMOOTH("xymapsmooth.5.dat",r,200,200,l,l)
CALL thispes%PLOT_DIRECTION1D("zscan.5.dat",500,45.D0,r(3),l2)
CALL thispes%PLOT_DIRECTION1D_CORRECTION("zscancorrection.5.dat",500,45.D0,r(3),l2)
CALL thisrawpes%PLOT_DIRECTION1D_SMOOTH("zscansmooth.5.dat",500,45.D0,r(3),l2)
r=(/0.D0,0.D0,7.D0/)
CALL thispes%PLOT_XYMAP("xymap.7.dat",r,200,200,l,l)
CALL thispes%PLOT_XYMAP_CORRECTION("xymapcorrection.7.dat",r,200,200,l,l)
CALL thisrawpes%PLOT_XYMAP_SMOOTH("xymapsmooth.7.dat",r,200,200,l,l)
CALL thispes%PLOT_DIRECTION1D("zscan.7.dat",500,45.D0,r(3),l2)
CALL thispes%PLOT_DIRECTION1D_CORRECTION("zscancorrection.7.dat",500,45.D0,r(3),l2)
CALL thisrawpes%PLOT_DIRECTION1D_SMOOTH("zscansmooth.7.dat",500,45.D0,r(3),l2)
r=(/0.D0,0.D0,9.D0/)
CALL thispes%PLOT_XYMAP("xymap.9.dat",r,200,200,l,l)
CALL thispes%PLOT_XYMAP_CORRECTION("xymapcorrection.9.dat",r,200,200,l,l)
CALL thisrawpes%PLOT_XYMAP_SMOOTH("xymapsmooth.9.dat",r,200,200,l,l)
CALL thispes%PLOT_DIRECTION1D("zscan.9.dat",500,45.D0,r(3),l2)
CALL thispes%PLOT_DIRECTION1D_CORRECTION("zscancorrection.9.dat",500,45.D0,r(3),l2)
CALL thisrawpes%PLOT_DIRECTION1D_SMOOTH("zscansmooth.9.dat",500,45.D0,r(3),l2)
r=(/0.D0,0.D0,10.D0/)
CALL thispes%PLOT_XYMAP("xymap.10.dat",r,200,200,l,l)
CALL thispes%PLOT_XYMAP_CORRECTION("xymapcorrection.10.dat",r,200,200,l,l)
CALL thisrawpes%PLOT_XYMAP_SMOOTH("xymapsmooth.10.dat",r,200,200,l,l)
CALL thispes%PLOT_DIRECTION1D("zscan.10.dat",500,45.D0,r(3),l2)
CALL thispes%PLOT_DIRECTION1D_CORRECTION("zscancorrection.10.dat",500,45.D0,r(3),l2)
CALL thisrawpes%PLOT_DIRECTION1D_SMOOTH("zscansmooth.10.dat",500,45.D0,r(3),l2)
! CRP3D sites and pairpots
CALL thispes%PLOT_SITIOS(1000)
CALL thispes%PLOT_PAIRPOTS(1000)

CALL ETIME(tiempoarr,tiempo)
WRITE(*,*) "Total time: ", tiempo
END PROGRAM CRP_TEST
