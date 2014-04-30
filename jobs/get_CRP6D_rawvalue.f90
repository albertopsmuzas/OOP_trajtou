PROGRAM TEST_CRP6D
! Initial declarations
USE CRP6D_MOD
USE DEBUG_MOD
IMPLICIT NONE
! Variables
TYPE(CRP6D) :: thispes
REAL(KIND=8),DIMENSION(6) :: x
REAL(KIND=8),DIMENSION(6) :: dvdu
REAL(KIND=8) :: v
CALL SET_VERBOSE_MODE(.TRUE.)
READ(*,*) x
WRITE(*,*) "*****************************************"
WRITE(*,*) "********* GET CRP6D POTENTIAL************"
WRITE(*,*) "*****************************************"
! STEP 1: READ CRP6D INPUT FILES
CALL thispes%READ("INcrp6d.inp")
CALL thispes%RAWINTERPOL()
!CALL thispes%INTERPOL_NEW_RZGRID(200,400)
CALL thispes%GET_V_AND_DERIVS_SMOOTH(x,v,dvdu)
WRITE(*,*) "Potential (au): ", v
WRITE(*,*) "dvdx (au): ", dvdu(1)
WRITE(*,*) "dvdy (au): ", dvdu(2)
WRITE(*,*) "dvdz (au): ", dvdu(3)
WRITE(*,*) "dvdr (au): ", dvdu(4)
WRITE(*,*) "dvdtheta (au): ", dvdu(5)
WRITE(*,*) "dvdphi (au): ", dvdu(6)
CALL EXIT(0)
END PROGRAM TEST_CRP6D
