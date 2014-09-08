PROGRAM  TEST_DYNAMICS
USE DEBUG_MOD
USE DYNDIATOMIC_MOD
IMPLICIT NONE
TYPE(Dyndiatomic) :: this
REAL(KIND=4),DIMENSION(2) :: timearr
REAL(KIND=4) :: timer
! RUN SECTION: GABBA, GABBA HEY! =========================================================================
!
! STEP 0: HELLO! -------------------------------------------------------------------------------------
WRITE(*,*) "*******************************************" 
WRITE(*,*) "**** TEST DYNAMICS 6D program executed ****"
WRITE(*,*) "*******************************************" 
! STEP 0: System specifications
CALL SET_VERBOSE_MODE(.TRUE.)
CALL ETIME(timearr,timer)
! STEP 1: START UP DYNAMICS
CALL this%INITIALIZE("INdynamics6d.inp")
!
! STEP 3: RUN DYNAMICS
CALL this%RUN()
!
CALL ETIME(timearr,timer)
WRITE(*,*) "User time: ", timearr(1)
WRITE(*,*) "System time: ", timearr(2)
WRITE(*,*) "Total time: ",timer
END PROGRAM TEST_DYNAMICS
