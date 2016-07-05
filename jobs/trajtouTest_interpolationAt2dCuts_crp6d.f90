PROGRAM TEST_CRP6D
! Initial declarations
use SYSTEM_MOD
use CRP6D_MOD, only: CRP6D
use DEBUG_MOD, only: VERBOSE_WRITE
use UNITS_MOD, only: Length
IMPLICIT NONE
! Variables
CHARACTER(LEN=1024):: luaFile
REAL(KIND=4),DIMENSION(2):: timearr
REAL(KIND=4):: timer
TYPE(CRP6D):: thispes
TYPE(CRP6D):: thispesraw
INTEGER(KIND=4):: i,j
CHARACTER(LEN=11):: filename
CHARACTER(LEN=15):: filename1
CHARACTER(LEN=13):: realname
CHARACTER(LEN=17):: realname1
TYPE(Length):: r,z
REAL(KIND=8):: r1,r2,z1,z2
REAL(KIND=8),DIMENSION(6):: pos
! STEP 1: INITIALIZE SYSTEM VIA LUA CONFIG FILE
SELECT CASE(command_argument_count())
   CASE(1)
      CALL GET_COMMAND_ARGUMENT(1,luaFile)
   CASE DEFAULT
      WRITE(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      WRITE(0,*) "It is only needed one string, which is a config lua file"
      CALL EXIT(1)
END SELECT
CALL INITIALIZE_SYSTEM(trim(luaFile))
CALL VERBOSE_WRITE('##############################################')
CALL VERBOSE_WRITE('########## TEST CRP6D INTERPOLATION ##########')
CALL VERBOSE_WRITE('##############################################')
CALL ETIME(timearr,timer)
CALL thispes%READ(fileName=trim(luaFile), tableName='pes')
CALL thispesraw%READ(fileName=trim(luaFile), tableName='pes')
CALL thispes%INTERPOL()
CALL thispesraw%RAWINTERPOL()
! STEP 2: PLOT SOME GRAPHS
CALL r%READ(0.4D0,"angst")
CALL z%READ(0.25D0,"angst")
CALL r%TO_STD()
CALL z%TO_STD()
z1=z%getvalue()
r1=r%getvalue()
CALL r%READ(2.28D0,"angst")
CALL z%READ(3.8D0,"angst")
CALL r%TO_STD()
CALL z%TO_STD()
z2=z%getvalue()
r2=r%getvalue()
! STEP 3: PRINT GRAPHS
DO i=1,thispes%nsites
   DO j = 1, thispes%wyckoffsite(i)%n2dcuts
      WRITE(filename,'(I1,A10)') j,"-cut2d.dat"
      WRITE(realname,'(I1,A1,A11)') i,"-",filename
      CALL thispes%wyckoffsite(i)%zrcut(j)%interrz%PLOT_XYMAP(realname,(/r1,z1/),100,100,r2-r1,z2-z1)
      pos(1)=thispes%wyckoffsite(i)%zrcut(j)%x
      pos(2)=thispes%wyckoffsite(i)%zrcut(j)%y
      pos(3)=z1
      pos(4)=r1
      pos(5)=thispes%wyckoffsite(i)%zrcut(j)%theta
      pos(6)=thispes%wyckoffsite(i)%zrcut(j)%phi
      WRITE(filename,'(I1,A10)') j,"-inter.dat"
      WRITE(realname,'(I1,A1,A11)') i,"-",filename
      CALL thispes%PLOT_RZMAP(pos,100,100,r2-r1,z2-z1,realname)
      WRITE(filename,'(I1,A10)') j,"-intat.dat"
      WRITE(realname,'(I1,A1,A11)') i,"-",filename
      CALL thispes%PLOT_ATOMIC_INTERAC_RZ(pos,100,100,r2-r1,z2-z1,realname)
      WRITE(filename1,'(I1,A14)') j,"-cut2d.raw.dat"
      WRITE(realname1,'(I1,A1,A15)') i,"-",filename1
      CALL thispesraw%wyckoffsite(i)%zrcut(j)%interrz%PLOTDATA(realname1)
      WRITE(filename1,'(I1,A14)') j,"-cut2d.smt.dat"
      WRITE(realname1,'(I1,A1,A15)') i,"-",filename1
      CALL thispes%wyckoffsite(i)%zrcut(j)%interrz%PLOTDATA(realname1)
   END DO
END DO
CALL ETIME(timearr,timer)
! PRINT TIMES
CALL VERBOSE_WRITE("****************** RUN TIME ***************************")
CALL VERBOSE_WRITE('',"User time: ",real(timearr(1),kind=8))
CALL VERBOSE_WRITE('',"System time: ", real(timearr(2),kind=8))
CALL VERBOSE_WRITE('',"Total time: ",real(timer,kind=8))
CALL VERBOSE_WRITE("******************************************************")
END PROGRAM TEST_CRP6D
