program trajtouGetInfo_vacuumPotCoeff
   ! Initial declarations
   use EXTRAPOL_TO_VACUUM_MOD, only: VacuumPot
   implicit none
   ! Variables 
   character(len=1024):: inputFile
   character(len=1024):: outputFile
   type(VacuumPot):: thisPot
   ! GET TO THE CHOPPAAAAAH ! ==========================
   select case(command_argument_count())
   case(1)
         CALL GET_COMMAND_ARGUMENT(1,inputFile)
         call thisPot%initialize( trim(inputFile) )
         ! call thisPot%rPot%printCoeffInfo()
         call thisPot%rPot%printCoeffInfo()
   case(2)
         CALL GET_COMMAND_ARGUMENT(1,inputFile)
         CALL GET_COMMAND_ARGUMENT(2,outputFile)
         call thisPot%initialize( trim(inputFile) )
         call thisPot%rPot%printCoeffInfo( fileName=trim(outputFile) )
   case default
      write(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      write(0,*) "Expected number of arguments: 1, vacuumpot file: output given in standard channel"
      write(0,*) "Expected number of arguments: 2, vacuumpot file,outputfile: output given in output File"
      call exit(1)
   end select
end program trajtouGetInfo_vacuumPotCoeff
