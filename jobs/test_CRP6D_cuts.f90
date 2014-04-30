PROGRAM TEST_CRP6D
! Initial declarations
USE CRP6D_MOD
USE DEBUG_MOD
IMPLICIT NONE
! Variables
TYPE(CRP6D) :: thispes
TYPE(CRP6D) :: thisrawpes
CALL SET_DEBUG_MODE(.FALSE.)
CALL SET_VERBOSE_MODE(.FALSE.)
! STEP 1: READ CRP6D INPUT FILES
CALL thispes%READ("INcrp6d.inp")
CALL thisrawpes%READ("INcrp6d.inp")
CALL thispes%INTERPOL()
CALL thisrawpes%RAWINTERPOL()
! STEP 2: DO GRAPHS
CALL thispes%PLOT1D_PHI(1000,&
   (/4.43140778D0,2.55113028D0,3.77945227D0,1.13383568D0,1.91986218D0,0.D0/),"phicut.close.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,&
   (/4.43140778D0,2.55113028D0,3.77945227D0,1.13383568D0,1.91986218D0,0.D0/),"phicutsmooth.close.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,&
   (/4.43140778D0,2.55113028D0,3.77945227D0,1.13383568D0,1.91986218D0,0.D0/),"phicutraw.close.dat")
CALL thispes%PLOT1D_PHI(1000,&
   (/10.10058618D0,6.29278802D0,5.66917840D0,1.51178091D0,1.91986218D0,0.D0/),"phicut.far.dat")
CALL thispes%PLOT1D_PHI_SMOOTH(1000,&
   (/10.10058618D0,6.29278802D0,5.66917840D0,1.51178091D0,1.91986218D0,0.D0/),"phicutsmooth.far.dat")
CALL thisrawpes%PLOT1D_PHI_SMOOTH(1000,&
   (/10.10058618D0,6.29278802D0,5.66917840D0,1.51178091D0,1.91986218D0,0.D0/),"phicutraw.far.dat")
CALL thispes%PLOT1D_THETA(1000,&
   (/10.10058618D0,6.29278802D0,5.66917840D0,1.36060282D0,0.D0,0.55850536D0/),"thetacut.close.dat")
CALL thispes%PLOT1D_THETA_SMOOTH(1000,&
   (/10.10058618D0,6.29278802D0,5.66917840D0,1.36060282D0,0.D0,0.55850536D0/),"thetacutsmooth.close.dat")
CALL thisrawpes%PLOT1D_THETA_SMOOTH(1000,&
   (/10.10058618D0,6.29278802D0,5.66917840D0,1.36060282D0,0.D0,0.55850536D0/),"thetacutraw.close.dat")
CALL thispes%PLOT1D_THETA(1000,&
   (/10.10058618D0,6.29278802D0,7.55890453D0,1.36060282D0,0.D0,0.55850536D0/),"thetacut.far.dat")
CALL thispes%PLOT1D_THETA_SMOOTH(1000,&
   (/10.10058618D0,6.29278802D0,7.55890453D0,1.36060282D0,0.D0,0.55850536D0/),"thetacutsmooth.far.dat")
CALL thisrawpes%PLOT1D_THETA_SMOOTH(1000,&
   (/10.10058618D0,6.29278802D0,7.55890453D0,1.36060282D0,0.D0,0.55850536D0/),"thetacutraw.far.dat")
CALL thispes%PLOT1D_Z(3000,&
   (/10.10058618D0,6.29278802D0,0.D0,1.36060282D0,1.91986218D0,0.55850536D0/),"zcut.dat")
CALL thispes%PLOT1D_Z_SMOOTH(3000,&
   (/10.10058618D0,6.29278802D0,0.D0,1.36060282D0,1.91986218D0,0.55850536D0/),"zcutsmooth.dat")
CALL thispes%PLOT1D_Z(3000,&
   (/0.D0,0.D0,0.D0,1.36060282D0,1.91986218D0,0.55850536D0/),"zcut00.dat")
CALL thispes%PLOT1D_Z_SMOOTH(3000,&
   (/0.D0,0.D0,0.D0,1.36060282D0,1.91986218D0,0.55850536D0/),"zcutsmooth00.dat")
CALL thisrawpes%PLOT1D_Z_SMOOTH(3000,&
   (/0.D0,0.D0,0.D0,1.36060282D0,1.91986218D0,0.55850536D0/),"zcutraw00.dat")
CALL thispes%PLOT1D_Z(3000,&
   (/0.D0,0.D0,0.D0,1.36060282D0,0.D0,0.55850536D0/),"zcut00-0.dat")
CALL thispes%PLOT1D_Z_SMOOTH(3000,&
   (/0.D0,0.D0,0.D0,1.36060282D0,0.D0,0.55850536D0/),"zcutsmooth00-0.dat")
CALL thisrawpes%PLOT1D_Z_SMOOTH(3000,&
   (/0.D0,0.D0,0.D0,1.36060282D0,0.D0,0.55850536D0/),"zcutraw00-0.dat")
CALL thisrawpes%PLOT1D_Z_SMOOTH(3000,&
   (/10.10058618D0,6.29278802D0,0.D0,1.36060282D0,1.91986218D0,0.55850536D0/),"zcutraw.dat")
CALL thispes%PLOT1D_R(1000,&
   (/10.10058618D0,6.29278802D0,7.55890453D0,0.D0,1.91986218D0,0.55850536D0/),"rcut.far.dat")
CALL thispes%PLOT1D_R_SMOOTH(1000,&
   (/10.10058618D0,6.29278802D0,7.55890453D0,0.D0,1.91986218D0,0.55850536D0/),"rcutsmooth.far.dat")
CALL thisrawpes%PLOT1D_R_SMOOTH(1000,&
   (/10.10058618D0,6.29278802D0,7.55890453D0,0.D0,1.91986218D0,0.55850536D0/),"rcutraw.far.dat")
CALL thispes%PLOT1D_R(1000,&
   (/10.10058618D0,6.29278802D0,5.66917840D0,0.D0,1.91986218D0,0.55850536D0/),"rcut.close.dat")
CALL thispes%PLOT1D_R_SMOOTH(1000,&
   (/10.10058618D0,6.29278802D0,5.66917840D0,0.D0,1.91986218D0,0.55850536D0/),"rcutsmooth.close.dat")
CALL thisrawpes%PLOT1D_R_SMOOTH(1000,&
   (/10.10058618D0,6.29278802D0,5.66917840D0,0.D0,1.91986218D0,0.55850536D0/),"rcutraw.close.dat")
! STEP 3: WITH DIFFERENT GRID
!CALL thispes%INTERPOL_NEW_RZGRID(50,100)
!CALL thispes%PLOT1D_PHI(1000,&
   !(/4.43140778D0,2.55113028D0,3.77945227D0,1.13383568D0,1.91986218D0,0.D0/),"phicut.close.50-100.dat")
!CALL thispes%PLOT1D_PHI(1000,&
   !(/10.10058618D0,6.29278802D0,5.66917840D0,1.51178091D0,1.91986218D0,0.D0/),"phicut.far.50-100.dat")
!CALL thispes%PLOT1D_THETA(1000,&
   !(/10.10058618D0,6.29278802D0,5.66917840D0,1.36060282D0,0.D0,0.55850536D0/),"thetacut.close.50-100.dat")
!CALL thispes%PLOT1D_THETA(1000,&
   !(/10.10058618D0,6.29278802D0,7.55890453D0,1.36060282D0,0.D0,0.55850536D0/),"thetacut.far.50-100.dat")
!CALL thispes%INTERPOL_NEW_RZGRID(100,200)
!CALL thispes%PLOT1D_PHI(1000,&
   !(/4.43140778D0,2.55113028D0,3.77945227D0,1.13383568D0,1.91986218D0,0.D0/),"phicut.close.100-200.dat")
!CALL thispes%PLOT1D_PHI(1000,&
   !(/10.10058618D0,6.29278802D0,5.66917840D0,1.51178091D0,1.91986218D0,0.D0/),"phicut.far.100-200.dat")
!CALL thispes%PLOT1D_THETA(1000,&
   !(/10.10058618D0,6.29278802D0,5.66917840D0,1.36060282D0,0.D0,0.55850536D0/),"thetacut.close.100-200.dat")
!CALL thispes%PLOT1D_THETA(1000,&
   !(/10.10058618D0,6.29278802D0,7.55890453D0,1.36060282D0,0.D0,0.55850536D0/),"thetacut.far.100-200.dat")
CALL EXIT(0)
END PROGRAM TEST_CRP6D
