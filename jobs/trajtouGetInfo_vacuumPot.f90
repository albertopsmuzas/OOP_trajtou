program trajtouGetInfo_vacuumPot
   ! Initial declarations
   use EXTRAPOL_TO_VACUUM_MOD, only: VacuumPot
   implicit none
   ! Variables 
   character(len=1024):: inputFile
   type(VacuumPot):: thisPot
   ! GET TO THE CHOPPAAAAAH ! ==========================
   select case(command_argument_count())
   case(1)
         CALL GET_COMMAND_ARGUMENT(1,inputFile)
   case default
      write(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      write(0,*) "Expected number of arguments: 1, Lua config file"
      call exit(1)
   end select
   call thisPot%initialize( trim(inputFile) ) 
   write(*,*) '####################################'
   write(*,*) '# SOME INFO FROM VACUUMPOT'
   write(*,*) '####################################'
   write(*,*) 'Equilibrium radius, Req. (bohr): ',thisPot%getReq()
   write(*,*) 'Potential at Req., V(Req) (hartree): ',thisPot%getPotMin()
   write(*,*) 'Energy shift to put V=0 at vacuum level, (hartree): ',thisPot%getScaleFactor()
   write(*,*) 'First grid point: ',thisPot%rPot%x(1)
   write(*,*) 'Last grid point: ',thisPot%rPot%x(thisPot%rPot%n)
   write(*,*) 'Potential at first grid point: ',thisPot%getPot(thisPot%rPot%x(1))
   write(*,*) 'Potential at last grid point: ',thisPot%getPot(thisPot%rPot%x(thisPot%rPot%n))
   write(*,*) 'First derivative at first grid point: ',thisPot%getDeriv(thisPot%rPot%x(1))
   write(*,*) 'First derivative at last grid point: ',thisPot%getDeriv(thisPot%rPot%x(thisPot%rPot%n))
end program trajtouGetInfo_vacuumPot
