!##################################################
! PROGRAM: molec2atom
!> @brief
!! Program to go from molecular coordinates: X,Y,Z,R,Theta,Phi
!! to atomic coordinates. Only works with diatomic molecules
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014 
!> @version 1.0
!##################################################
PROGRAM molec2atom
! Initial declarations
USE DEBUG_MOD
USE SYSTEM_MOD
USE CRP6D_MOD
IMPLICIT NONE
! variables
REAL(KIND=8) :: ma,mb
REAL(KIND=8),DIMENSION(6) :: molec_coord, atomic_coord
! Run yeah, tun
CALL INITIALIZE_SYSTEM('blabla.inp')
READ(*,*) molec_coord(:)
WRITE(*,*) "***********************************************"
WRITE(*,*) "** FROM MOLECULAR COORDINATES TO ATOMIC ONES **"
WRITE(*,*) "***********************************************"
ma=system_mass(1)
mb=system_mass(2)
CALL FROM_MOLECULAR_TO_ATOMIC(ma,mb,molec_coord,atomic_coord)
WRITE(*,*) "Atom A: ",atomic_coord(1:3)
WRITE(*,*) "Atom B: ",atomic_coord(4:6)
CALL EXIT(0)
END PROGRAM molec2atom
