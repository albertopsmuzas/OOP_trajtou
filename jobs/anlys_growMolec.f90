program analysiS
! Initial declarations
use SYSTEM_MOD, only: INITIALIZE_SYSTEM
use DEBUG_MOD, only: VERBOSE_WRITE
use DIFFRACTIONGROW_MOD, only: Allowed_peaksGrow
implicit nonE
! Variables
type(Allowed_peaksGrow):: anlyDiff
character(LEN=1024):: auxString
real(kind=8):: morseEd, morseWidth
integer(kind=4):: dJ
real(kind=4),dimension(2):: timearr
real(kind=4):: timer
! HEY HO, LET'S GO --------------------------------------
! Get input
select case(command_argument_count())
   case(1)
      call get_command_argument(1,auxString)
      call initialize_systeM(trim(auxString))
      call verbose_write("**************************************************************************")
      call verbose_write("********* DIFFRACTION FOR GROW DIATOMICS (WITH QUANTUM BINNING) **********")
      call verbose_write("**************************************************************************")
      call etime(timearr,timer)
      call anlyDiff%initialize()
      call anlyDiff%assignTrajsToPeaks()
      call anlyDiff%sortByDiffOrder()
      !CALL anlyDiff%print_labmomenta_and_angles()
      call anlyDiff%printSeenPeaks()
      call etime(timearr,timer)

   case default
      write(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      write(0,*) "FOUR arguments needed: lua config file, dJ, morseEd, morseWidth"
      call exit(1)
end selecT
! Time
call verbose_write("****************** RUN TIME ***************************")
call verbose_write('',"User time: ",real(timearr(1),kind=8))
call verbose_write('',"System time: ",real(timearr(2),kind=8))
call verbose_write('',"Total time: ",real(timer,kind=8))
call verbose_write("******************************************************")
end program analysis
