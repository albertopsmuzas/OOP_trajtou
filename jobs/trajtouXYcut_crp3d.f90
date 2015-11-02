program trajtouZcut_crp6d
   ! initial declarations
   use SYSTEM_MOD
   use DEBUG_MOD, only: verbose_write
   use CRP3D_MOD, only: Crp3d
   implicit none
   ! Variables
   character(len=1024):: auxString
   type(Crp3d):: myPes
   real(kind=8),dimension(3):: initPoint
   integer(kind=4),dimension(2):: gridPoint
   real(kind=8),dimension(2):: gridLength
   real(kind=4),dimension(2):: timeArr
   real(kind=4):: timer
   ! Run section
   select case(command_argument_count())
      case(9)
         call get_command_argument(1,auxString) ! Read Lua file
         call initialize_system(trim(auxString))
         call myPes%initialize()
         call get_command_argument(2,auxString) ! Read init X
         read(auxString,*) initPoint(1)
         call get_command_argument(3,auxString) ! Read init Y
         read(auxString,*) initPoint(2)
         call get_command_argument(4,auxString) ! Read init Z
         read(auxString,*) initPoint(3)
         call get_command_argument(5,auxString) ! Read points in X direction
         read(auxString,*) gridPoint(1)
         call get_command_argument(6,auxString) ! Read points in Y direction
         read(auxString,*) gridPoint(2)
         call get_command_argument(7,auxString) ! Read length of X grid in a.u.
         read(auxString,*) gridLength(1)
         call get_command_argument(8,auxString) ! Read length of Y grid in a.u.
         read(auxString,*) gridLength(2)
         call get_command_argument(9,auxString) ! Read output filename
         call etime(timeArr,timer)
         call myPes%plot_xymap( filename=trim(auxString),&
                                init_xyz=initPoint(:),&
                                nxpoints=gridPoint(1),&
                                nyPoints=gridPoint(2),&
                                Lx=gridLength(1),&
                                Ly=gridLength(2) )
         call etime(timeArr,timer)
      case default
         write(0,*) 'ERR: bad number of arguments: ',command_argument_count()
         write(0,*) 'Expected number of arguments: 9'
         write(0,*) 'Lua config file; initial X,Y,Z, points in X, points in Y, length of X (bohrs), length of Y (bohrs)'
         call exit(1)
   end select
   call verbose_write("****************** RUN TIME ***************************")
   call verbose_write('',"User time: ",real(timeArr(1),kind=8))
   call verbose_write('',"System time: ",real(timeArr(2),kind=8))
   call verbose_write('',"Total time: ",real(timer,kind=8))
   call verbose_write("******************************************************")

end program trajtouZcut_crp6d
