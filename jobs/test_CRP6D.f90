PROGRAM TEST_CRP6D
! Initial declarations
USE CRP6D_MOD
USE DEBUG_MOD
USE UNITS_MOD
IMPLICIT NONE
! Variables
TYPE(CRP6D) :: thispes
TYPE(CRP6D) :: thispesraw
INTEGER(KIND=4) :: i,j
CHARACTER(LEN=11) :: filename
CHARACTER(LEN=15) :: filename1
CHARACTER(LEN=13) :: realname
CHARACTER(LEN=17) :: realname1
TYPE(Length) :: r,z
REAL(KIND=8) :: r1,r2,z1,z2
REAL(KIND=8),DIMENSION(6) :: pos
CALL SET_VERBOSE_MODE(.TRUE.)
! STEP 1: READ CRP6D INPUT FILES
CALL thispes%READ("INcrp6d.inp")
CALL thispesraw%READ("INcrp6d.inp")
CALL thispes%INTERPOL()
CALL thispesraw%RAWINTERPOL()

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
DO i=1,thispes%nsites
   DO j = 1, thispes%wyckoffsite(i)%n2dcuts
      WRITE(filename,'(I1,A10)') j,"-cut2d.dat"
      WRITE(realname,'(I1,A1,A11)') i,"-",filename
      CALL thispes%wyckoffsite(i)%zrcut(j)%interrz%PLOT_XYMAP(realname,(/r1,z1/),300,300,r2-r1,z2-z1)
      WRITE(filename,'(I1,A10)') j,"-dualg.dat"
      WRITE(realname,'(I1,A1,A11)') i,"-",filename
      CALL thispes%wyckoffsite(i)%zrcut(j)%interrz%PLOT_DUALDERIVS_AT_GRID(realname)
      pos(1)=thispes%wyckoffsite(i)%zrcut(j)%x
      pos(2)=thispes%wyckoffsite(i)%zrcut(j)%y
      pos(3)=z1
      pos(4)=r1
      pos(5)=thispes%wyckoffsite(i)%zrcut(j)%theta
      pos(6)=thispes%wyckoffsite(i)%zrcut(j)%phi
      WRITE(filename,'(I1,A10)') j,"-inter.dat"
      WRITE(realname,'(I1,A1,A11)') i,"-",filename
      CALL thispes%PLOT_RZMAP(pos,300,300,r2-r1,z2-z1,realname)
   END DO
END DO
DO i=1,thispesraw%nsites
   DO j = 1, thispesraw%wyckoffsite(i)%n2dcuts
      WRITE(filename1,'(I1,A14)') j,"-cut2d.raw.dat"
      WRITE(realname1,'(I1,A1,A15)') i,"-",filename1
      CALL thispesraw%wyckoffsite(i)%zrcut(j)%interrz%PLOT_XYMAP(realname1,(/r1,z1/),300,300,r2-r1,z2-z1)
      WRITE(filename1,'(I1,A14)') j,"-dualg.raw.dat"
      WRITE(realname1,'(I1,A1,A15)') i,"-",filename1
      CALL thispesraw%wyckoffsite(i)%zrcut(j)%interrz%PLOT_DUALDERIVS_AT_GRID(realname1)
   END DO
END DO

CALL thispes%farpot%PLOT(100,"vacuumpot.dat")

CALL EXIT(0)
END PROGRAM TEST_CRP6D
