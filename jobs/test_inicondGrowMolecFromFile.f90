program TEST_INICOND6D_INPUT
! Initial declarations
use DEBUG_MOD, only: verbose_write,set_verbose_mode
use INITSURFACEGROW_MOD, only: InitSurfaceGrowDiatomic
use SYSTEM_MOD, only: initialize_system
implicit none
! Some objects
type(InitSurfaceGrowDiatomic) :: iniconDat
real(kind=4),dimension(2):: timearr
real(kind=4):: timer
character(len=1024):: luaFile
! STEP 0: HELLO! & system specifications
write(*,*) "**************************************************************"
write(*,*) "******************* TEST INITSURFACEGROW *********************"
write(*,*) "**************************************************************"
! STEP 1: READ PES, INTERPOLATION NOT NEEDED
select case(command_argument_count())
case(1)
   CALL GET_COMMAND_ARGUMENT(1,luaFile)
case default
   write(0,*) "ERR: Bad number of arguments: ",command_argument_count()
   write(0,*) "It is only needed one string, which is a config lua file"
   call exit(1)
end select
! STEP 2: GENERATE NEW INITIAL CONDITIONS
call etime(timearr,timer)
call initialize_system(trim(luaFile))
call set_verbose_mode(.true.) ! we want verbose output only for inicond
call inicondat%initialize()
call inicondat%generate_trajs_from_file('init.raw.dat')
call etime(timearr,timer)
! STEP 3: PRINT TIME 
call verbose_write("****************** RUN TIME ***************************")
call verbose_write('',"User time: ",real(timearr(1),kind=8))
call verbose_write('',"System time: ",real(timearr(2),kind=8))
call verbose_write('',"Total time: ",real(timer,kind=8))
call verbose_write("*******************************************************")
end program TEST_INICOND6D_INPUT
