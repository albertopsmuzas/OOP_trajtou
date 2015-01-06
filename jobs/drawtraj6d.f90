!##################################################
! PROGRAM: DRAWTRAJ6D
!> @brief
!! Creates a XYZ file
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jul/2014
!> @version 1.0
!##################################################
PROGRAM DRAWTRAJ6D
   USE DRAWTRAJ_MOD
! Initial declarations
IMPLICIT NONE
TYPE(Drawtraj):: this
CHARACTER(LEN=40):: dynamicsfile,outputfile
INTEGER(KIND=4):: norder
WRITE(*,*) "TYPE order,input file and output file"
READ(*,*) norder,dynamicsfile,outputfile
CALL this%INITIALIZE("XYZ","INsurface.inp")
CALL this%DRAW(norder,dynamicsfile,outputfile)
END PROGRAM DRAWTRAJ6D
