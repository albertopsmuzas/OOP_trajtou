!############################################
! PROGRAM: trajtouCheck_interpol_crp3d
!
!> @brief
!! Prints sitio graphs as read from the input in conjunction with
!! the crp3d and the smooth potential (for each sitio)
!
!> @details
!! Useful to check that sitio inputs and the interpolated potential
!! have the same values.
!--------------------------------------------
program trajtouCheckCRP3DInterpol
! Initial declarations
use SYSTEM_MOD, only: initialize_system
use CRP3D_MOD, only: CRP3D
use DEBUG_MOD, only: verbose_write
implicit none
! Variables
real(kind=4),dimension(2):: timeArr
real(kind=4):: timer
character(1012):: luaFile, auxString
type(CRP3D):: thisPes
integer(kind=4):: i ! counters
integer(kind=4):: nPoints
real(kind=8),dimension(3):: xyz
real(kind=8):: L
! Hey, ho, let's go !
call etime(timeArr,timer)
select case(command_argument_count())
   case(2)
      call get_command_argument(number=1, value=luaFile)
      call get_command_argument(number=2, value=auxString)
      read(auxString,'(I10)') nPoints

   case default
      write(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      write(0,*) "Two argument needed: Lua input file, nPoints"
      call exit(1)
end select
call initialize_system(trim(luaFile))
call thisPes%initialize()
do i=1,size(thisPes%all_sites)
   xyz(1)=thisPes%all_sites(i)%x
   xyz(2)=thisPes%all_sites(i)%y
   xyz(3)=thisPes%all_sites(i)%z(1)
   L=thisPes%all_sites(i)%z(thisPes%all_sites(i)%n)-thisPes%all_sites(i)%z(1)
   write(auxString,'(I10)') i
   call thisPes%plot_z(npoints=nPoints, xyz=xyz(:), L=L, filename='sitio'//trim(adjustL(auxString))//'.fullInterpolation.dat')
   call thisPes%plot_z_smooth(npoints=nPoints, xyz=xyz(:), L=L, filename='sitio'//trim(adjustL(auxString))//'.smoothInterpolation.dat')
enddo
do i=1,size(thisPes%all_sites)
   write(auxString,'(I10)') i
   call thisPes%all_sites(i)%plot_data(filename='sitio'//trim(adjustL(auxString))//'.dat')
enddo

call etime(timeArr,timer)
! Last step, verbose messages
call verbose_write("****************** RUN TIME ***************************")
call verbose_write('',"User time: ",real(timearr(1),kind=8))
call verbose_write('',"System time: ",real(timearr(2),kind=8))
call verbose_write('',"Total time: ",real(timer,kind=8))
call verbose_write("******************************************************")
end program trajtouCheckCRP3DInterpol
