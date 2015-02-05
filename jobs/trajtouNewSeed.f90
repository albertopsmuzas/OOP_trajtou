program trajtouNewSeed
   ! initial declarations
   implicit none
   ! Variables
   integer(kind=4):: nseed,clock
   integer(kind=4),dimension(:),allocatable:: seed
   integer(kind=4):: i ! counters
   ! Run section
   call random_seed(size=nseed)
   allocate(seed(nseed))
   call random_seed(put=seed(:))
   call system_clock(count=clock)
   seed(:) = clock+ 37*(/ (i - 1,i = 1,nseed) /)
   write(*,*) seed(:)
end program trajtouNewSeed
