! ############################################################################
! PROGRAM: trajtouGetGraph_vacuumPot
!> @brief
!! Program designed to obtain an ascii file which contains the interpolated
!! vacuum potential given as input file. If no option is given as output file, information will be printed
!! inside a file called interpolatedVacuumPot.dat. This file will be replaced if
!! it exits.
!-----------------------------------------------------------------------------
program trajtouGetGraph_vacuumPot
   ! Initial declarations
   use EXTRAPOL_TO_VACUUM_MOD, only: VacuumPot
   implicit none
   ! Variables 
   character(len=1024):: inputFile
   character(len=1024):: auxString
   type(VacuumPot):: thisPot
   integer(kind=4):: nPoints
   ! GET TO THE CHOPPAAAAAH ! ==========================
   select case(command_argument_count())
   case(2)
         call get_command_argument(1,inputFile)
         call get_command_argument(2,auxString)
         read(auxString,*) nPoints
         auxString='interpolatedVacuumPot.dat'
   case(3)
         call get_command_argument(1,inputFile)
         call get_command_argument(2,auxString)
         read(auxString,*) nPoints
         call get_command_argument(3,auxString)
   case default
      write(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      write(0,*) "Expected number of arguments: 2, vacuumpot file, number of points in the graph"
      write(0,*) "Expected number of arguments: 3, vacuumpot file, number of points in the graph, outputfile"
      call exit(1)
   end select
   call thisPot%initialize( trim(inputFile) ) 
   call thisPot%plot( npoints=nPoints,filename=trim(auxString) )
end program trajtouGetGraph_vacuumPot
