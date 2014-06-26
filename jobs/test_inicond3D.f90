PROGRAM TEST_INICOND
! Initial declarations
USE DEBUG_MOD
USE INITATOM_MOD
USE CRP3D_MOD
IMPLICIT NONE
! Some objects
TYPE(Initatom) :: inicondat
TYPE(CRP3D) :: thispes
WRITE(*,*) "##############################"
WRITE(*,*) "#### TEST INICOND PROGRAM ####"
WRITE(*,*) "##############################"
! STEP 1: SET DEBUG/VERBOSE MODE
CALL SET_VERBOSE_MODE(.TRUE.)
! STEP 2: READ PES, INTERPOLATION NOT NEEDED
CALL thispes%READ("INcrp3d.inp")
! STEP 3: GENERATE NEW INITIAL CONDITIONS
CALL inicondat%INITIALIZE("INinicond3d.inp")
CALL inicondat%GENERATE_TRAJS(thispes)
END PROGRAM TEST_INICOND
