program TEST_CRP6D
! Initial declarations
use SYSTEM_MOD, only: initialize_system
use CRP6D_MOD,  only: Crp6d
use DEBUG_MOD,  only: verbose_write
use UNITS_MOD,  only: pi
implicit none
! Variables
type(Crp6d):: thispes
!type(Crp6d):: thisrawpes
character(len=1024):: luaFile
real(kind=8),dimension(6) :: r
real(kind=8),parameter:: ucell=5.44335612578d0
! STEP 1: READ CRP6D INPUT FILES
select case(command_argument_count())
   case(1)
      call get_command_argument(1,luaFile)
   case default
      write(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      write(0,*) "It is only needed one string, which is a config lua file"
      call exit(1)
end select
call initialize_system(trim(luaFile))
call verbose_write('##############################################')
call verbose_write('######### TEST SOME PRE-DEFINED CUTS #########')
call verbose_write('##############################################')
call thispes%initialize()
!CALL thisrawpes%READ(trim(luaFile),'pes')
!CALL thisrawpes%RAWINTERPOL()
!CALL thisrawpes%INTERPOL_NEW_RZGRID(25,50)
! CRP for cartwheel
r=[-1.d0,-1.d0,13.6060281568d0,1.42d0,0.d0,0.d0]
call thispes%PLOT_XYMAP(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z7.2.cart.dat")
r=[-1.d0,-1.d0,13.0391103169d0,1.42d0,0.d0,0.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z6.9.cart.dat")
r=[-1.d0,-1.d0,11.3383567973d0,1.42d0,0.d0,0.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z6.cart.dat")
r=[-1.d0,-1.d0,9.44863066443d0,1.42d0,0.d0,0.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z5.cart.dat")
r=[-1.d0,-1.d0,7.55890453154d0,1.42d0,0.d0,0.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z4.cart.dat")
r=[-1.d0,-1.d0,5.66917839866d0,1.42d0,0.d0,0.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z3.cart.dat")
r=[-1.d0,-1.d0,3.77945226577d0,1.42d0,0.d0,0.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z2.cart.dat")
call thispes%atomiccrp(1)%plot_xyMap(init_xyz=r(1:3),nxpoints=150,nypoints=150,Lx=ucell,Ly=ucell,filename="xycut_Z2.atom.dat")
r=[-1.d0,-1.d0,1.88972613289d0,1.42d0,0.d0,0.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z1.cart.dat")
! crp for helicopter
r=[-1.d0,-1.d0,13.6060281568d0,1.42d0,pi/2.d0,0.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z7.2.hel.dat")
r=[-1.d0,-1.d0,13.0391103169d0,1.42d0,pi/2.d0,0.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z6.9.hel.dat")
r=[-1.d0,-1.d0,11.3383567973d0,1.42d0,pi/2.d0,0.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z6.hel.dat")
r=[-1.d0,-1.d0,9.44863066443d0,1.42d0,pi/2.d0,0.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z5.hel.dat")
r=[-1.d0,-1.d0,7.55890453154d0,1.42d0,pi/2.d0,0.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z4.hel.dat")
r=[-1.d0,-1.d0,5.66917839866d0,1.42d0,pi/2.d0,0.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z3.hel.dat")
r=[-1.d0,-1.d0,3.77945226577d0,1.42d0,pi/2.d0,0.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z2.hel.dat")
r=[-1.d0,-1.d0,1.88972613289d0,1.42d0,pi/2.d0,0.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z1.hel.dat")
! crp for tilted
r=[-1.d0,-1.d0,13.6060281568d0,1.42d0,pi/2.d0,pi/6.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z7.2.tilt.dat")
r=[-1.d0,-1.d0,13.0391103169d0,1.42d0,pi/2.d0,pi/6.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z6.9.tilt.dat")
r=[-1.d0,-1.d0,11.3383567973d0,1.42d0,pi/2.d0,pi/6.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z6.tilt.dat")
r=[-1.d0,-1.d0,9.44863066443d0,1.42d0,pi/2.d0,pi/6.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z5.tilt.dat")
r=[-1.d0,-1.d0,7.55890453154d0,1.42d0,pi/2.d0,pi/6.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z4.tilt.dat")
r=[-1.d0,-1.d0,5.66917839866d0,1.42d0,pi/2.d0,pi/6.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z3.tilt.dat")
r=[-1.d0,-1.d0,3.77945226577d0,1.42d0,pi/2.d0,pi/6.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z2.tilt.dat")
r=[-1.d0,-1.d0,1.88972613289d0,1.42d0,pi/2.d0,pi/6.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z1.tilt.dat")
! crp for tilted2
r=[-1.d0,-1.d0,13.6060281568d0,1.42d0,pi/2.d0,pi/4.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z7.2.tilt2.dat")
r=[-1.d0,-1.d0,13.0391103169d0,1.42d0,pi/2.d0,pi/4.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z6.9.tilt2.dat")
r=[-1.d0,-1.d0,11.3383567973d0,1.42d0,pi/2.d0,pi/4.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z6.tilt2.dat")
r=[-1.d0,-1.d0,9.44863066443d0,1.42d0,pi/2.d0,pi/4.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z5.tilt2.dat")
r=[-1.d0,-1.d0,7.55890453154d0,1.42d0,pi/2.d0,pi/4.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z4.tilt2.dat")
r=[-1.d0,-1.d0,5.66917839866d0,1.42d0,pi/2.d0,pi/4.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z3.tilt2.dat")
r=[-1.d0,-1.d0,3.77945226577d0,1.42d0,pi/2.d0,pi/4.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z2.tilt2.dat")
r=[-1.d0,-1.d0,1.88972613289d0,1.42d0,pi/2.d0,pi/4.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z1.tilt2.dat")
! crp for tilted3
r=[-1.d0,-1.d0,13.6060281568d0,1.42d0,pi/4.d0,pi/4.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z7.2.tilt3.dat")
r=[-1.d0,-1.d0,13.0391103169d0,1.42d0,pi/4.d0,pi/4.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z6.9.tilt3.dat")
r=[-1.d0,-1.d0,11.3383567973d0,1.42d0,pi/4.d0,pi/4.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z6.tilt3.dat")
r=[-1.d0,-1.d0,9.44863066443d0,1.42d0,pi/4.d0,pi/4.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z5.tilt3.dat")
r=[-1.d0,-1.d0,7.55890453154d0,1.42d0,pi/4.d0,pi/4.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z4.tilt3.dat")
r=[-1.d0,-1.d0,5.66917839866d0,1.42d0,pi/4.d0,pi/4.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z3.tilt3.dat")
r=[-1.d0,-1.d0,3.77945226577d0,1.42d0,pi/4.d0,pi/4.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z2.tilt3.dat")
r=[-1.d0,-1.d0,1.88972613289d0,1.42d0,pi/4.d0,pi/4.d0]
call thispes%plot_xyMap(init_point=r,nxpoints=100,nypoints=100,Lx=ucell+2.d0,Ly=ucell+2.d0,filename="xycut_Z1.tilt3.dat")
! some Z cuts
r=[0.d0,0.d0,1.88972613289d0,1.42d0,0.d0,0.d0]
call thispes%plot1d_z( npoints=500,X=r(:),L=20.d0,filename='zCutTopLi.cart.dat' )
r=[0.d0,0.d0,1.88972613289d0,1.42d0,pi/2.d0,0.d0]
call thispes%plot1d_z( npoints=500,X=r(:),L=20.d0,filename='zCutTopLi.hel.dat' )
r=[0.d0,0.d0,1.88972613289d0,1.42d0,pi/2.d0,pi/6.d0]
call thispes%plot1d_z( npoints=500,X=r(:),L=20.d0,filename='zCutTopLi.tilt.dat' )
end program TEST_CRP6D
