!##################################################
! PROGRAM: cart2surf_program
!> @brief
!! Program to go from cartesian to surface coordinates: X,Y
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014 
!> @version 1.0
!##################################################
PROGRAM cart2surf_program
! Initial declarations
USE SURFACE_MOD
USE DEBUG_MOD
IMPLICIT NONE
! variables
TYPE(Surface) :: surf
REAL(KIND=8),DIMENSION(2) :: cart, aux
! Run yeah, tun
CALL SET_DEBUG_MODE(.TRUE.)
CALL surf%INITIALIZE("INsurface.inp")
READ(*,*) cart(:)
WRITE(*,*) "*******************************************"
WRITE(*,*) "** FROM CARTESIAN TO SURFACE COORDINATES **"
WRITE(*,*) "*******************************************"
WRITE(*,*) " Cartesian a.u. :     ",cart
aux=surf%cart2surf(cart)
WRITE(*,*) " Surface coordinates: ",aux
aux=surf%project_unitcell(cart)
WRITE(*,*) " Cartesian unit cell: ",aux
aux=surf%project_iwscell(cart)
WRITE(*,*) " Cartesian IW cell:   ",aux
CALL EXIT(0)
END PROGRAM cart2surf_program
