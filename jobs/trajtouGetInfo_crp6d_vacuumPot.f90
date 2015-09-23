program trajtouGetInfo_crp6d_vacuumPot
   ! Initial declarations
   use SYSTEM_MOD
   use CRP6D_MOD, only: Crp6d
   implicit none
   ! Variables 
   character(len=1024):: luaFile 
   type(Crp6d):: thisPes
   ! GET TO THE CHOPPAAAAAH ! ==========================
   select case(command_argument_count())
   case(1)
         CALL GET_COMMAND_ARGUMENT(1,luaFile)
   case default
      write(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      write(0,*) "Expected number of arguments: 1, Lua config file"
      call exit(1)
   end select
   call initialize_system(trim(luaFile))
   call thisPes%initialize() 
   write(*,*) '####################################'
   write(*,*) '# SOME INFO FROM VACUUMPOT'
   write(*,*) '####################################'
   write(*,*) 'Equilibrium radius, Req. (bohr): ',thisPes%farPot%getReq()
   write(*,*) 'Potential at Req., V(Req) (hartree): ',thisPes%farPot%getPotMin()
   write(*,*) 'Energy shift to put V=0 at vacuum level, (hartree): ',thisPes%farPot%getScaleFactor()
end program trajtouGetInfo_crp6d_vacuumPot
