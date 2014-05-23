PROGRAM TEST_CRP6D
! Initial declarations
USE CRP6D_MOD
USE DEBUG_MOD
IMPLICIT NONE
! Variables
TYPE(CRP6D) :: thispes
TYPE(CRP6D) :: thisrawpes
REAL(KIND=8),DIMENSION(6) :: pos
CALL SET_DEBUG_MODE(.FALSE.)
CALL SET_VERBOSE_MODE(.FALSE.)
! STEP 1: READ CRP6D INPUT FILES
CALL thispes%READ("INcrp6d.inp")
CALL thisrawpes%READ("INcrp6d.inp")
CALL thispes%INTERPOL()
CALL thispes%INTERPOL_NEW_RZGRID(200,200)
CALL thisrawpes%RAWINTERPOL()
CALL thisrawpes%INTERPOL_NEW_RZGRID(200,200)
! STEP 2: DO GRAPHS
! datagroup 1
pos=(/0.0D0,0.D0,4.724315D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.1.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.1.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.1.dat")

! datagroup 2
pos=(/0.D0,0.D0,5.669178D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.2.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.2.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.2.dat")

! datagroup 3
pos=(/0.D0,0.D0,4.724315D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.3.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.3.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.3.dat")

! datagroup 4
pos=(/0.D0,0.D0,5.669178D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.4.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.4.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.4.dat")

! datagroup 7
pos=(/2.721678D0,2.721678D0,5.669178D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.7.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.7.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.7.dat")

! datagroup 8
pos=(/2.721678D0,2.721678D0,4.724315D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.8.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.8.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.8.dat")

! datagroup 9
pos=(/2.721678D0,2.721678D0,5.669178D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.9.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.9.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.9.dat")

! datagroup 10
pos=(/1.814453D0,0.907225D0,4.724315D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.10.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.10.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.10.dat")

! datagroup 11
pos=(/1.814453D0,0.907225D0,5.669178D0,1.511781D0,1.570796D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.11.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.11.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.11.dat")

! datagroup 12
pos=(/1.814453D0,0.907225D0,4.724315D0,1.511781D0,0.523599D0,0.D0/)
CALL thispes%PLOT1D_PHI(1000,pos,"datagroup.12.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupsmooth.12.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,pos,"datagroupraw.12.dat")

CALL EXIT(0)
END PROGRAM TEST_CRP6D
