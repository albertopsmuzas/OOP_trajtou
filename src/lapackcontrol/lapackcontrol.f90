!#########################################################
! RAISON D'ÊTRE:
! - Compendium of control routines for lapack 
! UPDATES:
! - Created 12/Dec/2012 
! FUNCTIONALITY:
! - Routine to stop program if "info" /= 0
! IDEAS FOR THE FUTURE:
! - None
!##########################################################
MODULE LAPACKCONTROL_MOD
! Initial declarations
IMPLICIT NONE

CONTAINS
!###########################################################
!# SUBROUTINE: LAPACK_CHECK 
!###########################################################
! - Stop the program if last info value was non equal to zero
!-----------------------------------------------------------
SUBROUTINE LAPACK_CHECK(lapackname,info)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CHARACTER(LEN=*),INTENT(IN) :: lapackname  ! Lapack routine name
   INTEGER(KIND=4), INTENT(IN) :: info
   ! Run section
   IF (info/=0) THEN
      WRITE(0,*) "LAPACK ERR: ", lapackname
      WRITE(0,*) "LAPACK ERR: exit status is ", info
      CALL EXIT(info)
   END IF
   RETURN
END SUBROUTINE LAPACK_CHECK
END MODULE LAPACKCONTROL_MOD
