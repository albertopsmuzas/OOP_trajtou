PROGRAM TEST
USE DEBUG_MOD
USE SURFACE_MOD 
USE PES_MOD
IMPLICIT NONE
TYPE(Surface) :: newsurf
TYPE(PES) :: newpes
CALL SET_DEBUG_MODE(.TRUE.)
CALL newsurf%INITIALIZE("LiF(001)","surface.inp")
WRITE(*,*) "Is the surface initialized?"
WRITE(*,*) newsurf%is_initialized()
CALL newpes%INITIALIZE()
WRITE(*,*) "Is the new PES initialized?"
WRITE(*,*) newpes%is_initialized()
WRITE(*,*) newpes%get_alias()
WRITE(*,*) newpes%get_dimensions()
END PROGRAM TEST
