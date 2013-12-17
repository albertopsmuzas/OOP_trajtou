PROGRAM TEST
USE DEBUG_MOD
USE SURFACE_MOD 
IMPLICIT NONE
TYPE(Surface) :: newsurf
CALL SET_DEBUG_MODE(.TRUE.)
CALL newsurf%INITIALIZE("LiF(001)","surface.inp")
WRITE(*,*) "Is the surface initialized?"
WRITE(*,*) newsurf%is_initialized()
END PROGRAM TEST
