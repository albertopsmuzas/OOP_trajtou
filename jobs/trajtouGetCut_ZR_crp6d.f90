!###############################################################
! PROGRAM: trajtouGetCut_ZR_crp6d
!> @brief
!! Prints a file with a nice RZ cut (Z on Y axis, R on X axis)
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Sep/2015
!###############################################################
program trajtouGetCut_ZR_crp6d
! Initial declarations
use DEBUG_MOD
use CRP6D_MOD
use SYSTEM_MOD
! use some modules?
implicit none
! Variables
real(kind=4),dimension(2):: timeArr
real(kind=8),dimension(6):: x
real(kind=8):: lr,lz
real(kind=4) timer
type(Crp6d):: thisPes
character(len=1024):: luaFile
character(len=1024):: auxString
integer(kind=4):: nr,nz
! GABBA GABBA HEY! ===============================
!
! STEP 1: INITIALIZE SYSTEM VIA LUA CONFIG FILE
select case(command_argument_count())
   case(12)
      call GET_COMMAND_ARGUMENT(1,luaFile)   ! get Lua input file
      call GET_COMMAND_ARGUMENT(2,auxString) ! get initial X value (au)
      read(auxstring,*) x(1)
      call GET_COMMAND_ARGUMENT(3,auxString) ! get initial Y value (au)
      read(auxstring,*) x(2)
      call GET_COMMAND_ARGUMENT(4,auxString) ! get initial Z value (au)
      read(auxstring,*) x(3)
      call GET_COMMAND_ARGUMENT(5,auxString) ! get initial R value (au)
      read(auxstring,*) x(4)
      call GET_COMMAND_ARGUMENT(6,auxString) ! get initial THETA value (rad)
      read(auxstring,*) x(5)
      call GET_COMMAND_ARGUMENT(7,auxString) ! get initial PHI value (rad)
      read(auxstring,*) x(6)
      call get_command_argument(8,auxString) ! get length of R grid in au
      read(auxString,*) lr
      call get_command_argument(9,auxString) ! get length of Z grid in au
      read(auxString,*) lz
      call get_command_argument(10,auxString) ! get number of points in R grid
      read(auxString,*) nr
      call get_command_argument(11,auxString) ! get number of points in Z grid
      read(auxString,*) nz
      call get_command_argument(12,auxString) ! get output filename
   case default
      write(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      write(0,*) "12 arguments needed: lua config file, 6 real numbers (geometry), 2 real numbers (grid length in R and Z),"
      write(0,*) '2 integer numbers (grid points in R and Z) and a string with the name of the output'
      call EXIT(1)
end select
! STEP 2: START COUNTING TIME, EVALUATE PES POINTS AND PRINT OUTPUT
call ETIME(timearr,timer)
call INITIALIZE_SYSTEM(trim(luafile))
call thispes%INITIALIZE()
call thisPes%PLOT_RZMAP(init_point=x(:),nxpoints=nr,nypoints=nz,Lx=Lr,Ly=Lz,fileName=trim(auxString))
call ETIME(timearr,timer)
! STEP 3: PRINT TIME IF VERBOSE MODE IS DECLARED IN INPUT FILE
call VERBOSE_WRITE("****************** RUN TIME ***************************")
call VERBOSE_WRITE('',"User time: ",real(timearr(1),kind=8))
call VERBOSE_WRITE('',"System time: ",real(timearr(2),kind=8))
call VERBOSE_WRITE('',"Total time: ",real(timer,kind=8))
call VERBOSE_WRITE("******************************************************")
end program trajtouGetCut_ZR_crp6d
