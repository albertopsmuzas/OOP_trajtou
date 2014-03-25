PROGRAM TEST_CRP6D
! Initial declarations
USE CRP6D_MOD
USE DEBUG_MOD
IMPLICIT NONE
! Variables
TYPE(CRP6D) :: thispes
INTEGER(KIND=4) :: i
INTEGER(KIND=4) :: nwyckoff,nzrcut,nxgrid,nygrid,nxpoints,nypoints
CHARACTER(LEN=19) :: filename1,filename2,filename3,filename4
CHARACTER(LEN=21) :: realname1,realname2,realname3,realname4
REAL(KIND=8) :: r1,r2,z1,z2
REAL(KIND=8) :: r,z
REAL(KIND=8) :: f1,dfdr1,dfdz1,d2fdrdz1
REAL(KIND=8) :: f2,dfdr2,dfdz2,d2fdrdz2
REAL(KIND=8) :: maxf,maxdfdr,maxdfdz,maxd2fdrdz
REAL(KIND=8) :: avf,avdfdr,avdfdz,avd2fdrdz
READ(*,*) nwyckoff,nzrcut,nxgrid,nygrid,nxpoints,nypoints
CALL SET_VERBOSE_MODE(.FALSE.)
CALL SET_DEBUG_MODE(.FALSE.)
CALL thispes%READ("INcrp6d.inp")
WRITE(*,*) "******************************************"
WRITE(*,*) "******* TESTING CRP6D 2DCUT GRID *********"
WRITE(*,*) "******************************************"
WRITE(*,*) "Selected wyckoff site: ", nwyckoff
WRITE(*,*) "Selected ZR-cut: ", nzrcut
WRITE(*,*) "Alias: ", thispes%wyckoffsite(nwyckoff)%zrcut(nzrcut)%getalias()
WRITE(*,*) "Original R grid size: ",&
            thispes%wyckoffsite(nwyckoff)%zrcut(nzrcut)%getgridsizeR()
WRITE(*,*) "Original Z grid size: ",&
            thispes%wyckoffsite(nwyckoff)%zrcut(nzrcut)%getgridsizeZ()
WRITE(*,*) "New grid size (R,Z): ", nxgrid, nygrid

r1=thispes%wyckoffsite(nwyckoff)%zrcut(nzrcut)%getfirstr()
r2=thispes%wyckoffsite(nwyckoff)%zrcut(nzrcut)%getlastr()
z1=thispes%wyckoffsite(nwyckoff)%zrcut(nzrcut)%getfirstz()
z2=thispes%wyckoffsite(nwyckoff)%zrcut(nzrcut)%getlastz()

WRITE(filename1,'(I1,A10)') nzrcut,"-xymap.oldgrid.dat"
WRITE(realname1,'(I1,A1,A11)') nwyckoff,"-",filename1
CALL thispes%wyckoffsite(nwyckoff)%zrcut(nzrcut)%interrz%PLOT_XYMAP(realname1,(/r1,z1/),nxpoints,nypoints,r2-r1,z2-z1)

CALL thispes%wyckoffsite(nwyckoff)%zrcut(nzrcut)%interrz%INTERPOL_NEWGRID(nxgrid,nygrid)
WRITE(filename2,'(I1,A10)') nzrcut,"-xymap.newgrid.dat"
WRITE(realname2,'(I1,A1,A11)') nwyckoff,"-",filename2
CALL thispes%wyckoffsite(nwyckoff)%zrcut(nzrcut)%interrz%PLOT_XYMAP(realname2,(/r1,z1/),nxpoints,nypoints,r2-r1,z2-z1)


WRITE(filename3,'(I1,A10)') nzrcut,"-absolut_error.dat"
WRITE(realname3,'(I1,A1,A11)') nwyckoff,"-",filename3
WRITE(filename4,'(I1,A10)') nzrcut,"-relativ_error.dat"
WRITE(realname4,'(I1,A1,A11)') nwyckoff,"-",filename4
OPEN (11,FILE=realname1,STATUS="old",ACTION="read")
OPEN (12,FILE=realname2,STATUS="old",ACTION="read")
OPEN (13,FILE=realname3,STATUS="replace",ACTION="write")
OPEN (14,FILE=realname4,STATUS="replace",ACTION="write")
maxf=0.D0
maxdfdr=0.D0
maxdfdz=0.D0
maxd2fdrdz=0.D0
avf=0.D0
avdfdr=0.D0
avdfdz=0.D0
avd2fdrdz=0.D0
DO i = 1, nxpoints*nypoints
   READ(11,*) r,z,f1,dfdr1,dfdz1,d2fdrdz1
   READ(12,*) r,z,f2,dfdr2,dfdz2,d2fdrdz2
   WRITE(13,*) r,z,dabs(f1-f2),dabs(dfdr1-dfdr2),dabs(dfdz1-dfdz2),dabs(d2fdrdz1-d2fdrdz2)
   WRITE(14,*) r,z,dabs(f1-f2)/f1,dabs(dfdr1-dfdr2)/dfdr1,dabs(dfdz1-dfdz2)/dfdz1,dabs(d2fdrdz1-d2fdrdz2)/d2fdrdz1
   ! Some switches
   SELECT CASE(dabs(f1-f2)>=maxf)
      CASE(.TRUE.)
         maxf=dabs(f1-f2)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(dabs(dfdr1-dfdr2)>=maxdfdr)
      CASE(.TRUE.)
         maxdfdr=dabs(dfdr1-dfdr2)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(dabs(dfdz1-dfdz2)>=maxdfdz)
      CASE(.TRUE.)
         maxdfdz=dabs(dfdz1-dfdz2)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(dabs(d2fdrdz1-d2fdrdz2)>=maxd2fdrdz)
      CASE(.TRUE.)
         maxd2fdrdz=dabs(d2fdrdz1-d2fdrdz2)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
  ! Calculate averages
   avf=avf+dabs(f1-f2)
   avdfdr=avdfdr+dabs(dfdr1-dfdr2)
   avdfdz=avdfdz+dabs(dfdz1-dfdz2)
   avd2fdrdz=avd2fdrdz+dabs(d2fdrdz1-d2fdrdz2)
END DO
CLOSE(11)
CLOSE(12)
CLOSE(13)
CLOSE(14)
avf=avf/(nxpoints*nypoints)
avdfdr=avdfdr/(nxpoints*nypoints)
avdfdz=avdfdz/(nxpoints*nypoints)
avd2fdrdz=avd2fdrdz/(nxpoints*nypoints)
WRITE(*,*) "Max abs. err. f: ",maxf
WRITE(*,*) "Max abs. err. dfdr: ",maxdfdr
WRITE(*,*) "Max abs. err. dfdz: ",maxdfdz
WRITE(*,*) "Max abs. err. d2fdrdz: ",maxd2fdrdz
WRITE(*,*) "Average abs. err. f: ",avf
WRITE(*,*) "Average abs. err. dfdr: ",avdfdr
WRITE(*,*) "Average abs. err. dfdz: ",avdfdz
WRITE(*,*) "Average abs. err. d2fdrdz: ",avd2fdrdz
CALL EXIT(0)
END PROGRAM TEST_CRP6D
