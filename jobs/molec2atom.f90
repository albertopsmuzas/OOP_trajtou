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
USE CRP6D_MOD
USE DEBUG_MOD
IMPLICIT NONE
! variables
TYPE(CRP6D) :: thispes
REAL(KIND=8) :: ma,mb
REAL(KIND=8),DIMENSION(6) :: molec_coord, atomic_coord
! Run yeah, tun
CALL thispes%READ("INcrp6d.inp")
READ(*,*) molec_coord(:)
WRITE(*,*) "***********************************************"
WRITE(*,*) "** FROM MOLECULAR COORDINATES TO ATOMIC ONES **"
WRITE(*,*) "***********************************************"
ma=thispes%atomdat(1)%getmass()
mb=thispes%atomdat(2)%getmass()
CALL FROM_MOLECULAR_TO_ATOMIC(ma,mb,molec_coord,atomic_coord)
WRITE(*,*) "Atom A: ",atomic_coord(1:3)
WRITE(*,*) "Atom B: ",atomic_coord(4:6)
CALL EXIT(0)
END PROGRAM molec2atom
