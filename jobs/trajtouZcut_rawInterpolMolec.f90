program trajtouZcut_crp6d
   ! initial declarations
   use SYSTEM_MOD
   use DEBUG_MOD, only: verbose_write
   use CRP6D_MOD, only: Crp6d
   implicit none
   ! Variables
   character(len=1024):: auxString
   type(Crp6d):: myPes
   real(kind=8),dimension(6):: initPoint
   real(kind=8):: lastZ,intervalZ
   real(kind=4),dimension(2):: timeArr
   real(kind=4):: timer
   ! Run section
   select case(command_argument_count())
      case(8)
         call get_command_argument(1,auxString) ! Read Lua file
         call initialize_system(trim(auxString))
         call myPes%read(filename=trim(auxString), tablename='pes')
         call myPes%rawInterpol
         call get_command_argument(2,auxString) ! Read init X
         read(auxString,*) initPoint(1)
         call get_command_argument(3,auxString) ! Read init Y
         read(auxString,*) initPoint(2)
         call get_command_argument(4,auxString) ! Read init Z
         read(auxString,*) initPoint(3)
         call get_command_argument(5,auxString) ! Read init R
         read(auxString,*) initPoint(4)
         call get_command_argument(6,auxString) ! Read init THETA
         read(auxString,*) initPoint(5)
         call get_command_argument(7,auxString) ! Read init PHI
         read(auxString,*) initPoint(6)
         call get_command_argument(8,auxString) ! Read last Z
         read(auxString,*) lastZ
         intervalZ=lastZ-initPoint(6)
         call etime(timeArr,timer)
         call myPes%PLOT1D_Z_SMOOTH(npoints=1000, X=initPoint, filename='OUTPLOTzcutRaw.out')
         call etime(timeArr,timer)

      case default
         write(0,*) 'ERR: bad number of arguments: ',command_argument_count()
         write(0,*) 'Expected number of arguments: 8'
         write(0,*) 'Lua config file; initial X,Y,Z,R,THETA,PHI and last Z point'
         call exit(1)
   end select
   call verbose_write("****************** RUN TIME ***************************")
   call verbose_write('',"User time: ",real(timeArr(1),kind=8))
   call verbose_write('',"System time: ",real(timeArr(2),kind=8))
   call verbose_write('',"Total time: ",real(timer,kind=8))
   call verbose_write("******************************************************")

end program trajtouZcut_crp6d
