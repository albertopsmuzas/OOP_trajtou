!##################################################
! PROGRAM: getCut_Z
!> @brief
!! Get a Z cut
!
!> @date Dec/2016
!> @version 1.0
!##################################################
program getCut_Z
! Initial declarations
use SYSTEM_MOD
use DEBUG_MOD, only: verbose_write
use CRP3D_MOD, only: CRP3D
! use some modules?
implicit none
! Variables
real(kind=4),dimension(2):: timearr
real(kind=8),dimension(3):: x
type(CRP3D):: thispes
real(kind=4):: timer
real(kind=8):: len
integer(kind=4):: nPoints
character(len=1024):: luafile
character(len=1024):: auxstring
! GABBA GABBA HEY! ===============================
!
! STEP 1: INITIALIZE SYSTEM VIA LUA CONFIG FILE
select case(command_argument_count())
case(7)
      call get_command_argument(1,luaFile)
      call get_command_argument(2,auxString)
      read(auxString,*) x(1)
      call get_command_argument(3,auxString)
      read(auxString,*) x(2)
      call get_command_argument(4,auxString)
      read(auxString,*) x(3)
      call get_command_argument(5,auxString)
      read(auxString,*) len
      call get_command_argument(6,auxString)
      read(auxString,*) nPoints
      call get_command_argument(7,auxString)
      
case default
   write(0,*) "ERR: Bad number of arguments: ",command_argument_count()
   write(0,*) "Four arguments needed: lua config file, 4 real numbers, 1 integer and a fileName"
   call exit(1)
end select
call etime(timearr,timer)
call initialize_system(trim(luaFile))
call thisPes%initialize()
call thisPes%plot_Z(npoints=nPoints,xyz=x(:),L=len,fileName=trim(auxString))
call etime(timearr,timer)
! STEP 2: PRINT TIME 
CALL VERBOSE_WRITE("****************** RUN TIME ***************************")
CALL VERBOSE_WRITE('',"User time: ",real(timearr(1),kind=8))
CALL VERBOSE_WRITE('',"System time: ",real(timearr(2),kind=8))
CALL VERBOSE_WRITE('',"Total time: ",real(timer,kind=8))
CALL VERBOSE_WRITE("******************************************************")
end program getCut_z
