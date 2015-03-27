!##################################################
! PROGRAM: TEST_SYSTEMINPUT
!> @brief
!! Test if system part of Lua conf file is read correctly
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!##################################################
program getPot
! Initial declarations
use SYSTEM_MOD
use DEBUG_MOD
use LINK_PES_MOD
implicit none
! Variables
real(kind=4),dimension(2):: timeArr
real(kind=8),dimension(:),allocatable:: x,dvdx
real(kind=8):: v
class(PES),allocatable:: thispes
real(kind=4):: timer
character(len=1024):: auxstring
! GABBA GABBA HEY! ===============================
!
! STEP 1: INITIALIZE SYSTEM VIA LUA CONFIG FILE
call etime(timeArr,timer)
call get_command_argument(1,auxString)
select case(trim(auxString))
case('Crp3d')
   select case( command_argument_count() )
   case(5)
      call get_command_argument( 2,auxString )
      call initialize_system( trim(auxString) )
      allocate( Crp3d:: thisPes )
      call thisPes%initialize()
      allocate( x(3) )
      allocate( dvdx(3) )
      call get_command_argument(3,auxstring)
      read(auxstring,*) x(1)
      call get_command_argument(4,auxstring)
      read(auxstring,*) x(2)
      call get_command_argument(5,auxstring)
      read(auxstring,*) x(3)
   case default
      write(0,*) "ERR: bad number of arguments: ",command_argument_count()
      write(0,*) "Expected number of arguments: 5: Pestype, Lua file, x, y, z"
      call exit(1)
   end select

case('Crp6d')
   select case( command_argument_count() )
   case(8)
      call get_command_argument( 2,auxString )
      call initialize_system( trim(auxString) )
      allocate( Crp6d:: thisPes )
      call thisPes%initialize()
      allocate( x(6) )
      allocate( dvdx(6) )
      call get_command_argument(3,auxstring)
      read(auxstring,*) x(1)
      call get_command_argument(4,auxstring)
      read(auxstring,*) x(2)
      call get_command_argument(5,auxstring)
      read(auxstring,*) x(3)
      call get_command_argument(6,auxstring)
      read(auxstring,*) x(4)
      call get_command_argument(7,auxstring)
      read(auxstring,*) x(5)
      call get_command_argument(8,auxstring)
      read(auxstring,*) x(6)
   case default
      write(0,*) "ERR: bad number of arguments: ",command_argument_count()
      write(0,*) "Expected number of arguments: 8: Pestype, Lua file, x, y, z, r, theta, phi"
      call exit(1)
   end select

case('H2LiF001')
   select case( command_argument_count() )
   case(7)
      allocate( PES_H2LiF001:: thisPes )
      call thisPes%initialize()
      allocate( x(6) )
      allocate( dvdx(6) )
      call get_command_argument(2,auxstring)
      read(auxstring,*) x(1)
      call get_command_argument(3,auxstring)
      read(auxstring,*) x(2)
      call get_command_argument(4,auxstring)
      read(auxstring,*) x(3)
      call get_command_argument(5,auxstring)
      read(auxstring,*) x(4)
      call get_command_argument(6,auxstring)
      read(auxstring,*) x(5)
      call get_command_argument(7,auxstring)
      read(auxstring,*) x(6)
   case default
      write(0,*) "ERR: bad number of arguments: ",command_argument_count()
      write(0,*) "Expected number of arguments: 7: Pestype, x, y, z, r, theta, phi"
      call exit(1)
   end select

case('HLiF001_WS')
   select case( command_argument_count() )
   case(4)
      allocate( PES_HLiF001_WS:: thisPes )
      call thisPes%initialize()
      allocate( x(3) )
      allocate( dvdx(3) )
      call get_command_argument(2,auxstring)
      read(auxstring,*) x(1)
      call get_command_argument(3,auxstring)
      read(auxstring,*) x(2)
      call get_command_argument(4,auxstring)
      read(auxstring,*) x(3)
   case default
      write(0,*) "ERR: bad number of arguments: ",command_argument_count()
      write(0,*) "Expected number of arguments: 4: Pestype, x, y, z"
      call exit(1)
   end select

case('HLiF001_NS')
   select case( command_argument_count() )
   case(4)
      allocate( PES_HLiF001_NS:: thisPes )
      call thisPes%initialize()
      allocate( x(3) )
      allocate( dvdx(3) )
      call get_command_argument(2,auxstring)
      read(auxstring,*) x(1)
      call get_command_argument(3,auxstring)
      read(auxstring,*) x(2)
      call get_command_argument(4,auxstring)
      read(auxstring,*) x(3)
   case default
      write(0,*) "ERR: bad number of arguments: ",command_argument_count()
      write(0,*) "Expected number of arguments: 4: Pestype, x, y, z"
      call exit(1)
   end select

case default
   write(0,*) 'ERR: wrong specified PES: '//trim(auxString)
   write(0,*) 'Implemented ones: Crp3d, Crp6d, HLiF001_WS, HLiF001_NS, H2LiF001'
   write(0,*) 'Warning: Case-sensitive'
   call exit(1)
end select
call thisPes%get_v_and_derivs(x,v,dvdx)
call etime(timeArr,timer)
   ! STEP 2: GET VALUES
call verbose_write('Potential calculated at:',x(:))
write(*,*) "Potential (a.u.): ", v
write(*,*) "Derivatives in auxiliar cartesian coord. (a.u.): ", dvdx(:)

! STEP 3: PRINT TIME 
call verbose_write("****************** RUN TIME ***************************")
call verbose_write('',"User time: ",real(timearr(1),kind=8))
call verbose_write('',"System time: ",real(timearr(2),kind=8))
call verbose_write('',"Total time: ",real(timer,kind=8))
call verbose_write("******************************************************")
end program getPot
