PROGRAM TEST_INICOND
! Initial declarations
USE INITATOM_MOD
USE DEBUG_MOD
USE CRP_MOD
IMPLICIT NONE
! Some objects
TYPE(Initatom) :: inicondat
TYPE(Atom_trajs) :: trajs
TYPE(CRP) :: thispes
WRITE(*,*) "##############################"
WRITE(*,*) "#### TEST INICOND PROGRAM ####"
WRITE(*,*) "##############################"
! STEP 1: SET DEBUG/VERBOSE MODE
CALL SET_DEBUG_MODE(.TRUE.)
CALL SET_VERBOSE_MODE(.TRUE.)
! STEP 2: READ SURFACE
CALL thispes%READ("crp.inp")
! STEP 3: GENERATE NEW INITIAL CONDITIONS
CALL inicondat%READ("inicond.inp")
CALL inicondat%GENERATE_TRAJS(thispes,trajs)
END PROGRAM TEST_INICOND
