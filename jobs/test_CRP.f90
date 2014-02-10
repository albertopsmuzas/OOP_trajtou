PROGRAM CRP_TEST
! Initial declarations
USE DEBUG_MOD
USE CRP_MOD
USE UNITS_MOD
IMPLICIT NONE
! Variables
TYPE(CRP) :: crp_pes
CHARACTER(LEN=40) :: filename1,filename2,filename3
INTEGER(KIND=4) :: nsites, npairpots
REAL(KIND=8),DIMENSION(3) :: r
TYPE(Length) :: aux
REAL(KIND=8) :: l
INTEGER(KIND=4) :: i ! counters
! Run section===================================================================================0
CALL SET_DEBUG_MODE(.TRUE.)
! STEP 1: HELLO!
WRITE(*,*) "***************************************" 
WRITE(*,*) "CRP_PES program executed"
WRITE(*,*) "***************************************" 

! STEP 2: INITIALIZE CRP PES:
CALL crp_pes%READ("crp.inp")

! STEP 3: DO Z INTERPOLATION EXTRACTING VASINT AND SMOOTHING SITES
CALL crp_pes%INTERPOL_Z()

! STEP 4: PLOT SOME GRAPHS
nsites=size(crp_pes%all_sites)
npairpots=size(crp_pes%all_pairpots)
!DO i = 1, nsites
   !WRITE(filename1,'(A7,I1,A4)') "Siteraw",i,".dat"
   !WRITE(filename2,'(A18,I1,A4)') "Sitenovasint_dense",i,".dat"
   !WRITE(filename3,'(A12,I1,A4)') "Siteinterpol",i,".dat"
   !CALL crp_pes%all_sites(i)%PLOT_DATA(filename1)
   !CALL crp_PES%all_sites(i)%interz%PLOT_INTERPOL(1000,filename2)
   !CALL crp_PES%all_sites(i)%interz%PLOT_DATA(filename3)
!END DO
!
!DO i = 1, npairpots
   !WRITE(filename1,'(A10,I1,A4)') "Pairpotraw",i,".dat"
   !WRITE(filename2,'(A21,I1,A4)') "Pairpotnovasint_dense",i,".dat"
   !WRITE(filename3,'(A15,I1,A4)') "Pairpotinterpol",i,".dat"
   !CALL crp_pes%all_pairpots(i)%PLOT_DATA(filename1)
   !CALL crp_PES%all_pairpots(i)%interz%PLOT_INTERPOL(1000,filename2)
   !CALL crp_PES%all_pairpots(i)%interz%PLOT_DATA(filename3)
!END DO
!
!CALL aux%READ(3.34D0,"angst")
!CALL aux%TO_STD()
r=(/0.D0,0.D0,3.34D0/)
CALL aux%READ(2.8805D0,"angst")
CALL aux%TO_STD()
l=aux%getvalue()
!CALL crp_pes%PLOT_XYMAP("xymap.dat",r,10,10,l,l)
CALL crp_pes%PLOT_DIRECTION1D("zscan.dat",50,45.D0,r(3),l)
CALL crp_pes%all_sites(3)%interz%PLOT_INTERPOL(50,"Site3graph")
END PROGRAM CRP_TEST
