PROGRAM TEST
USE SURFACE_MOD 
IMPLICIT NONE
TYPE(Surface) :: newsurf
CALL newsurf%INITIALIZE("LiF(001)","Surface.in")
WRITE(*,*) "Is the surface initialized?"
WRITE(*,*) newsurf%is_initialized()
END PROGRAM TEST
