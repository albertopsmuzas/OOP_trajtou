program trajtouGetGraph_Thetacut_crp6d
   ! initial declarations
   use SYSTEM_MOD
   use DEBUG_MOD, only: verbose_write
   use CRP6D_MOD, only: Crp6d
   implicit none
   ! Variables
   type(Crp6d):: myPes
   real(kind=8),dimension(6):: initPoint
   real(kind=4),dimension(2):: timeArr
   character(len=1024):: auxString
   real(kind=4):: timer
   ! Run section
   select case(command_argument_count())
      case(7)
         call get_command_argument(1,auxString) ! Read Lua file
         call initialize_system(trim(auxString))
         call myPes%initialize()
         call get_command_argument(2,auxString) ! Read init X
         read(auxString,*) initPoint(1)
         call get_command_argument(3,auxString) ! Read init Y
         read(auxString,*) initPoint(2)
         call get_command_argument(4,auxString) ! Read init Z
         read(auxString,*) initPoint(3)
         call get_command_argument(5,auxString) ! Read init R
         read(auxString,*) initPoint(4)
         call get_command_argument(6,auxString) ! Read init PHI
         read(auxString,*) initPoint(6)
         call get_command_argument(7,auxString) ! Read output fileName
         call etime(timeArr,timer)
         call myPes%plot1d_theta(npoints=1000,X=initPoint,filename=trim(auxString))
         call etime(timeArr,timer)

      case(6)
         call get_command_argument(1,auxString) ! Read Lua file
         call initialize_system(trim(auxString))
         call myPes%initialize()
         call get_command_argument(2,auxString) ! Read init X
         read(auxString,*) initPoint(1)
         call get_command_argument(3,auxString) ! Read init Y
         read(auxString,*) initPoint(2)
         call get_command_argument(4,auxString) ! Read init Z
         read(auxString,*) initPoint(3)
         call get_command_argument(5,auxString) ! Read init R
         read(auxString,*) initPoint(4)
         call get_command_argument(6,auxString) ! Read init PHI
         read(auxString,*) initPoint(6)
         call etime(timeArr,timer)
         call myPes%PLOT1D_PHI(npoints=1000,X=initPoint,filename='thetaCut_crp6d.dat')
         call etime(timeArr,timer)

      case default
         write(0,*) 'ERR: bad number of arguments: ',command_argument_count()
         write(0,*) 'Expected number of arguments: 7 or 6'
         write(0,*) 'Case 7: Lua config file; initial X,Y,Z,R,PHI and graph name'
         write(0,*) 'Case 6: Lua config file; initial X,Y,Z,R and PHI'
         call exit(1)
   end select
   call verbose_write("****************** RUN TIME ***************************")
   call verbose_write('',"User time: ",real(timeArr(1),kind=8))
   call verbose_write('',"System time: ",real(timeArr(2),kind=8))
   call verbose_write('',"Total time: ",real(timer,kind=8))
   call verbose_write("******************************************************")

end program trajtouGetGraph_Thetacut_crp6d
