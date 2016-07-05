program analysis
! Initial declarations
use SYSTEM_MOD, only: initialize_system
use DEBUG_MOD, only: verbose_write
use DIFFRACTIONATOMSURF_MOD, only: Allowed_peaksAtomSurf
implicit none
! Variables
type(Allowed_peaksAtomSurf):: anlyDiff
character(len=1024):: auxString
real(kind=4),dimension(2):: timearr
real(kind=4):: timer
! HEY HO, LET'S GO --------------------------------------
! Get input
select case(command_argument_count())
case(1)
   call get_command_argument(1,auxString)
   call initialize_system(trim(auxString))
   call verbose_write("******************************************************")
   call verbose_write("************** DIFFRACTION ANALYSIS 6D ***************")
   call verbose_write("******************************************************")
   call etime(timearr,timer)
   call anlyDiff%initialize()
   call anlyDiff%setup()
   call anlyDiff%assign_peaks()
   call anlyDiff%print_labmomenta_and_angles()
   call etime(timearr,timer)

case default
   write(0,*) "ERR: Bad number of arguments: ",command_argument_count()
   write(0,*) "ONE argument needed: lua config file"
   call exit(1)

end select
! Time
call verbose_write("****************** RUN TIME ***************************")
call verbose_write('',"User time: ",real(timearr(1),kind=8))
call verbose_write('',"System time: ",real(timearr(2),kind=8))
call verbose_write('',"Total time: ",real(timer,kind=8))
call verbose_write("******************************************************")
end program analysis
