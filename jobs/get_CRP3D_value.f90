PROGRAM GET_CRP_VALUE
! Initial declarations
USE DEBUG_MOD
USE CRP3D_MOD
USE UNITS_MOD
IMPLICIT NONE
! Variables
TYPE(CRP3D) :: crp_pes
REAL(KIND=8),DIMENSION(3) :: r,dvdu
REAL(KIND=8) :: v
! Run section===================================================================================0
CALL SET_DEBUG_MODE(.FALSE.)
CALL SET_VERBOSE_MODE(.FALSE.)
! STEP 1: INSERT VALUE FROM STD INPUT
READ(*,*) r
! STEP 2: INITIALIZE CRP PES:
CALL crp_pes%READ("crp.inp")

! STEP 3: DO Z INTERPOLATION EXTRACTING VASINT AND SMOOTHING SITES
CALL crp_pes%INTERPOL_Z()

! STEP 4: GIVE VALUE
CALL crp_pes%GET_V_AND_DERIVS(r,v,dvdu)
WRITE(*,*) "Potential (a.u.): ", v
WRITE(*,*) "Derivatives in auxiliar cartesian coord. (a.u.): ", dvdu(:)
CALL EXIT(0)
END PROGRAM GET_CRP_VALUE
