!#############################################
! PROGRAM: trajtouSymmetrizeCRP3DInput
!
!> @brief
!! This program symmetrizes a raw CRP3D input file.
!> @warning
!! - 6 argumets needed: zero(real8), unitsZero(string), vtop(real8), unitsVtop(string), inputFile(string)
!!   outputFile(string)
!> @param zero - @b real8, length at which this program is going to symmetrize points.
!> @param unitsZero - @b string, units of "zero" argument.
!> @param vtop - @b real8, only points with energy higher than "vtop" will be symmetrized.
!> @param unitsVtop - @b string, units of "vtop" argument.
!> @param inputFile - @b string, file that contains the raw data to be symmetrized.
!> @param outputFile - @b string, output file with the symmetrized data inside.
!---------------------------------------------
program trajtouSymmetrizeCRP3DInput
! Initial declarations
use CRP3D_MOD, only: SymmPoint
use UNITS_MOD, only: Length, Energy
use DEBUG_MOD, only: verbose_write
implicit none
! Variables -----------------------------------
real(kind=4),dimension(2):: timeArr
real(kind=4):: timer
real(kind=8):: auxReal
type(Length):: leng
type(Energy):: en
type(SymmPoint):: thisPoint
character(len=1012):: auxString
! Rock the Casbah! ----------------------------
call etime(timeArr,timer)
select case(command_argument_count())
   case(6)
      call get_command_argument(number=1, value=auxString)
      read(auxString,'(f15.10)') auxReal
      call get_command_argument(number=2, value=auxString)
      call leng%read(mag=auxReal,units=trim(auxString))
      call leng%to_std()
      call get_command_argument(number=3, value=auxString)
      read(auxString,'(f15.10)') auxReal
      call get_command_argument(number=4, value=auxString)
      call en%read(mag=auxReal,units=trim(auxString))
      call en%to_std()
      call get_command_argument(number=5, value=auxString)
      call thisPoint%initializeRaw(fileName=trim(auxString))
      call get_command_argument(number=6, value=auxString)
      call thisPoint%printSymmetrizedRawInput(zero=leng, vtop=en, fileName=trim(auxString))

   case default
      write(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      write(0,*) "Six argument needed: zero(real8), unitsZero(string), vPot(real8), unitsVPot(string), inputFile(string), outFile(string)"
      call exit(1)
end select
call verbose_write("****************** RUN TIME ***************************")
call verbose_write('',"User time: ",real(timearr(1),kind=8))
call verbose_write('',"System time: ",real(timearr(2),kind=8))
call verbose_write('',"Total time: ",real(timer,kind=8))
call verbose_write("******************************************************")

end program trajtouSymmetrizeCRP3DInput
