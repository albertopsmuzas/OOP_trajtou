PROGRAM  TEST_DYNAMICS
USE DEBUG_MOD
USE CRP_MOD
USE INITATOM_MOD
USE DYNATOM_MOD
IMPLICIT NONE
TYPE(Dynatom) :: dinamica
TYPE(Initatom) :: initialcond
TYPE(CRP) :: thispes
REAL(KIND=4),DIMENSION(2) :: timearr
REAL(KIND=4) :: timer
TYPE(Atom_trajs) :: trajs
! RUN SECTION: GABBA, GABBA HEY! =========================================================================
!
! STEP 0: HELLO! -------------------------------------------------------------------------------------
WRITE(*,*) "***************************************" 
WRITE(*,*) "TEST DYNAMICS program executed"
WRITE(*,*) "***************************************" 
!CALL SET_VERBOSE_MODE(.TRUE.)

CALL ETIME(timearr,timer)
! STEP 1: START UP OUR PES
CALL thispes%READ("crp.inp")
CALL thispes%INTERPOL_Z()
!
! STEP 2: INITIAL CONDITIONS
CALL initialcond%READ("inicond.inp")
CALL initialcond%GENERATE_TRAJS(thispes,trajs)
!
! STEP 3: RUN DYNAMICS
CALL dinamica%READ("dynamics.inp")
CALL dinamica%RUN(initialcond,thispes,trajs)
!
CALL ETIME(timearr,timer)
WRITE(*,*) "User time: ", timearr(1)
WRITE(*,*) "System time: ", timearr(2)
WRITE(*,*) "Total time: ",timer
END PROGRAM TEST_DYNAMICS
