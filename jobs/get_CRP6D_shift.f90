PROGRAM GET_CRP_DETAILS
! Initial declarations
USE DEBUG_MOD
USE CRP6D_MOD
IMPLICIT NONE
! Variables
TYPE(CRP6D) :: crp_pes
! Run section===================================================================================0
CALL SET_VERBOSE_MODE(.FALSE.)
CALL crp_pes%READ("INcrp6d.inp")
CALL crp_pes%INTERPOL()
WRITE(*,*)  crp_pes%farpot%getscalefactor()
CALL EXIT(0)
END PROGRAM GET_CRP_DETAILS
