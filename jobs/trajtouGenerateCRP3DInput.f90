!##################################################
! PROGRAM: trajtouGenerateCRP3DInput
!> @brief
!! Generate correct CRP3D input files from raw data
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2016
!> @version 1.0
!##################################################
program trajtouGenerateCRP3DInput
! Initial declarations
use DEBUG_MOD
use NEWINPUT3D_MOD, only: NewInput3d
implicit none
! Variables
real(kind=4),dimension(2):: timeArr
type(NewInput3d):: thisInput
real(kind=4):: timer
character(len=1024):: luaFile
! GABBA GABBA HEY! ===============================
!
! STEP 1: INITIALIZE SYSTEM VIA LUA CONFIG FILE
select case(command_argument_count())
   case(1)
      call GET_COMMAND_ARGUMENT(1,luaFile)
   case default
      write(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      write(0,*) "One argument needed: lua config file"
      call EXIT(1)
end select
call etime(timearr,timer)
call thisInput%initialize(fileName=trim(luaFile))
call etime(timearr,timer)

! STEP 2: PRINT TIME
call verbose_write("****************** RUN TIME ***************************")
call verbose_write('',"User time: ",real(timearr(1),kind=8))
call verbose_write('',"System time: ",real(timearr(2),kind=8))
call verbose_write('',"Total time: ",real(timer,kind=8))
call verbose_write("******************************************************")
end program trajtouGenerateCRP3DInput
