PROGRAM TEST_CRP6D
! Initial declarations
use SYSTEM_MOD, only: evaluateEnergyRovibrState,initialize_system
use DEBUG_MOD, only: verbose_write
IMPLICIT NONE
! Variables
CHARACTER(LEN=1024):: auxString
REAL(KIND=4),DIMENSION(2):: timearr
REAL(KIND=4):: timer
integer(kind=4),dimension(2):: rovibrState
real(kind=8):: eVibr,eRot,eTot
! STEP 1: INITIALIZE SYSTEM VIA LUA CONFIG FILE
call verbose_write('##################################################')
call verbose_write('########## EVALUATE ROVIBRATIONAL ENERGY #########')
call verbose_write('##################################################')
call etime(timearr,timer)
select case(command_argument_count())
case(3)
   call get_command_argument(1,auxString)
   call initialize_system(trim(auxString))
   call get_command_argument(2,auxString)
   read(auxString,*) rovibrState(1)
   call get_command_argument(3,auxString)
   read(auxString,*) rovibrState(2)

case default
   write(0,*) "ERR: Bad number of arguments: ",command_argument_count()
   write(0,*) "Expected number of arguments: 3: lua config file, "
   call exit(1)

end select
eTot=evaluateEnergyRovibrState([rovibrState(:),0],eVibr=eVibr,eRot=eRot)
write(*,'("Vibrational energy(au): ",F10.6)') eVibr
write(*,'("Rotational energy(au): ",F10.6)') eRot
write(*,'("Total Rovibr. energy(au): ",F10.6)') eTot
call etime(timearr,timer)
! PRINT TIMES
call verbose_write("****************** RUN TIME ***************************")
call verbose_write('',"User time: ",real(timearr(1),kind=8))
call verbose_write('',"System time: ", real(timearr(2),kind=8))
call verbose_write('',"Total time: ",real(timer,kind=8))
call verbose_write("******************************************************")
END PROGRAM TEST_CRP6D
