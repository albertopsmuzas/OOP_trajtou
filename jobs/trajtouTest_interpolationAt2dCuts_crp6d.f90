program TEST_CRP6D
! Initial declarations
use SYSTEM_MOD
use CRP6D_MOD, only: CRP6D
use DEBUG_MOD, only: VERBOSE_WRITE
implicit none
! Variables
character(len=1024):: luaFile
character(len=1024):: auxString1,auxString2,fileName
real(kind=4),dimension(2):: timearr
real(kind=4):: timer
type(CRP6D):: thispes
type(CRP6D):: thispesraw
integer(kind=4):: i,j
real(kind=8):: r1,r2,z1,z2
integer(kind=4):: nx,ny
real(kind=8),dimension(6):: pos
! STEP 1: INITIALIZE SYSTEM VIA LUA CONFIG FILE
select case(command_argument_count())
   case(1)
      call GET_COMMAND_ARGUMENT(1,luaFile)
   case default
      write(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      write(0,*) "It is only needed one string, which is a config lua file"
      call EXIT(1)
end select
call INITIALIZE_SYSTEM(trim(luaFile))
call VERBOSE_WRITE('##############################################')
call VERBOSE_WRITE('########## TEST CRP6D INTERPOLATION ##########')
call VERBOSE_WRITE('##############################################')
call ETIME(timearr,timer)
call thispes%READ(fileName=trim(luaFile), tableName='pes')
call thispesraw%READ(fileName=trim(luaFile), tableName='pes')
call thispes%INTERPOL()
call thispesraw%RAWINTERPOL()
! STEP 2: PRINT GRAPHS
do i=1,thispes%nsites
   do j = 1, thispes%wyckoffsite(i)%n2dcuts
      ! Get some information useful to the graphs
      write(auxString1,'(I2)') i
      write(auxString2,'(I2)') j
      nx=size( thisPes%wyckoffSite(i)%zrcut(j)%interrz%x(:) )
      ny=size( thisPes%wyckoffSite(i)%zrcut(j)%interrz%y(:) )
      r1=thisPes%wyckoffSite(i)%zrcut(j)%interrz%x(1)
      z1=thisPes%wyckoffSite(i)%zrcut(j)%interrz%y(1)
      r2=thisPes%wyckoffSite(i)%zrcut(j)%interrz%x(nx)
      z2=thisPes%wyckoffSite(i)%zrcut(j)%interrz%y(ny)
      pos(1)=thispes%wyckoffsite(i)%zrcut(j)%x
      pos(2)=thispes%wyckoffsite(i)%zrcut(j)%y
      pos(3)=z1
      pos(4)=r1
      pos(5)=thispes%wyckoffsite(i)%zrcut(j)%theta
      pos(6)=thispes%wyckoffsite(i)%zrcut(j)%phi
      ! Print interpolations performed at wyckoff sites
      fileName=trim(adjustL(auxString1))//"-"//trim(adjustL(auxString2))//"-RZ-crp6d_on_cut2d.dat"
      call thispes%wyckoffsite(i)%zrcut(j)%interrz%PLOT_XYMAP(trim(fileName),(/r1,z1/),100,100,r2-r1,z2-z1)
      ! Print interpolations from final potential
      fileName=trim(adjustL(auxString1))//"-"//trim(adjustL(auxString2))//"-RZ-crp6d_finalPot.dat"
      call thispes%PLOT_RZMAP(pos,100,100,r2-r1,z2-z1,trim(fileName))
      ! Print atomic potential
      fileName=trim(adjustL(auxString1))//"-"//trim(adjustL(auxString2))//"-RZ-atomicRepulsion.dat"
      call thispes%PLOT_ATOMIC_INTERAC_RZ(pos,100,100,r2-r1,z2-z1,trim(fileName))
      ! Print RZ raw info at grid points
      fileName=trim(adjustL(auxString1))//"-"//trim(adjustL(auxString2))//"-RZ-rawGrid.dat"
      call thispesraw%wyckoffsite(i)%zrcut(j)%interrz%PLOTDATA(trim(fileName))
      fileName=trim(adjustL(auxString1))//"-"//trim(adjustL(auxString2))//"-RZ-smoothGrid.dat"
      call thispes%wyckoffsite(i)%zrcut(j)%interrz%PLOTDATA(trim(fileName))
   end do
end do
call ETIME(timearr,timer)
! PRINT TIMES
call VERBOSE_WRITE("****************** RUN TIME ***************************")
call VERBOSE_WRITE('',"User time: ",real(timearr(1),kind=8))
call VERBOSE_WRITE('',"System time: ", real(timearr(2),kind=8))
call VERBOSE_WRITE('',"Total time: ",real(timer,kind=8))
call VERBOSE_WRITE("******************************************************")
end program TEST_CRP6D
