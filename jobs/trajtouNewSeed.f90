program trajtouNewSeed
   ! initial declarations
   use SYSTEM_MOD, only: system_iSeed, generate_seed, system_seedFilename
   implicit none
   ! Local variables
   integer(kind=1):: iStat
   ! Run section
   call generate_seed(iStat)
   select case(iStat)
      case(0)
         write(*,*) 'Seed read from '//system_seedFilename
      case(1)
         write(*,*) 'Seed generated from /dev/urandom'
      case(2)
         write(*,*) 'Seed generated from system clock'
      case default
         write(*,*) 'Unknown iStat specificator'
         call exit(1)
   end select
   write(*,*) system_iSeed(:)
end program trajtouNewSeed
