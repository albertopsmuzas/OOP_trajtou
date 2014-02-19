PROGRAM TEST_CRP6D
! Initial declarations
USE CRP6D_MOD
USE DEBUG_MOD
USE UNITS_MOD
IMPLICIT NONE
! Variables
TYPE(CRP6D) :: thispes
INTEGER(KIND=4) :: i
CHARACTER(LEN=30) :: filename
TYPE(Length) :: r,z
REAL(KIND=8) :: r1,r2,z1,z2
CALL SET_DEBUG_MODE(.TRUE.)
! STEP 1: READ CRP6D INPUT FILES
CALL thispes%READ("crp6d.inp")

! STEP 2: PLOT SOME GRAPHS
! Prepare some units
CALL r%READ(0.4D0,"angst")
CALL z%READ(0.25D0,"angst")
CALL r%TO_STD()
CALL z%TO_STD()
z1=z%getvalue()
r1=r%getvalue()
CALL r%READ(2.3D0,"angst")
CALL z%READ(4D0,"angst")
CALL r%TO_STD()
CALL z%TO_STD()
z2=z%getvalue()
r2=r%getvalue()

! Print graphs
DO i = 1, thispes%n
   WRITE(filename,'(I2,A9)') i,"cut2d.dat"
   filename=adjustl(filename)
   CALL thispes%corte2d(i)%interrz%PLOT_XYMAP(filename,(/r1,z1/),100,100,r2-r1,z2-z1)
END DO
   CALL thispes%corte2d(1)%interrz%PLOT_SPLINES(150)
CALL EXIT(0)
END PROGRAM TEST_CRP6D
