program trajtouPrintPolygon
   ! initial declarations
   use MATHS_MOD, only: parametricPolygonEquation
   use UNITS_MOD, only: pi
   implicit none
   ! Variables
   integer(kind=4):: nArguments
   character(len=1024):: auxString
   integer(kind=4):: N
   real(kind=8):: tiltingAngle
   real(kind=8):: radius
   real(kind=8),dimension(2):: center,position
   integer(kind=4):: i
   integer(kind=4):: nPoints
   ! Run section ......................
   nArguments=command_argument_count()
   select case( nArguments )
   case(6)
      call get_command_argument( number=1,value=auxString )
      read(auxString,'(I10)') nPoints
      call get_command_argument( number=2,value=auxString )
      read(auxString,'(I10)') N
      call get_command_argument( number=3,value=auxString )
      read(auxString,'(f15.10)') radius
      call get_command_argument( number=4,value=auxString )
      read(auxString,'(f15.10)') center(1)
      call get_command_argument( number=5,value=auxString )
      read(auxString,'(f15.10)') center(2)
      call get_command_argument( number=6,value=auxString )
      read(auxString,'(f15.10)') tiltingAngle
      tiltingAngle=tiltingAngle*pi/180.d0
   case default
      write(0,*) 'Incorrect number of arguments. Given: ',nArguments
      write(0,*) 'Expected: 6: number of points, number of sides, radius, center in X and Y, tilting angle(deg)'
      call exit(1)
   end select
   do i=0,nPoints
      call parametricPolygonEquation( theta=dfloat(i)*2.d0*pi/nPoints,r=radius,N=N,theta0=tiltingAngle,x0=center,x=position )
      write(*,*) position(:)
   enddo
end program trajtouPrintPolygon
