program trajtouTest_cheatCartwheelInput_crp6d
! Initial declarations
use CRP6D_MOD, only: Crp6d
use DEBUG_MOD, only: set_verbose_mode, set_debug_mode
use SYSTEM_MOD, only: initialize_system
implicit none
! Variables
type(Crp6d) :: thisPes
real(kind=8) :: maxDist
integer(kind=4):: wyckoffId,cut2dId,rumplingId
character(len=1024):: luaFile, auxString, auxString2
! Hey, Ho, Let's go! -------------------------
! Avoid verbosity
call SET_VERBOSE_MODE(.false.)
call SET_DEBUG_MODE(.false.)
! read arguments
select case ( command_argument_count() )
case(5)
   call get_command_argument(1,luaFile)
   call get_command_argument(2,auxString)
   read(auxString,*) wyckoffId
   call get_command_argument(3,auxString)
   read(auxString,*) cut2dId
   call get_command_argument(4,auxString)
   read(auxString,*) rumplingId
   call get_command_argument(5,auxString)
   read(auxString,*) maxDist
case default
   write(0,*) 'ERR: Bad number of arguments. Expected 5.'
   write(0,*) 'Number of arguments encountered: ',command_argument_count()
   write(0,*) 'Arguments needed: Lua conf. file(string), wyckoffId(int), cut2dId(int)'
   write(0,*) 'rumplingId(int) and maxDist(real)'
   call exit(1)
end select
! STEP 1: READ CRP6D INPUT FILES
call initialize_system( trim(luaFile) )
call thispes%READ(trim(luaFile),'pes')
! STEP 2: EXTRACT VACUUM POTENTIAL & CHEAT
call thispes%EXTRACT_VACUUMSURF()
call thispes%CHEAT_CARTWHEEL_ONTOP(wyckoff=wyckoffId,cut2d=cut2dId,topType=rumplingId,dmax=maxDist)
! STEP 3: ADD VACUUM POTENTIAL & PRINT
call thispes%ADD_VACUUMSURF()
write(auxString,'(I2)')  wyckoffId
write(auxString2,'(I2)') cut2dId
call thispes%wyckoffsite(wyckoffId)%zrcut(cut2dId)%PRINT_INPUT('NewInput_wyckoff-'//trim(adjustL(auxString))//'-'//trim(adjustL(auxString2))//'.inp')
! STEP 4: GOODBYE! 
call EXIT(0)
end program trajtouTest_cheatCartwheelInput_crp6d
