! program to test one or more stand alone PES. Just useful for
! quick tests
program patata
   ! initial declarations
   use pes_h2lif001_mod
   implicit none
   ! internal variables
   type(PES_h2lif001) thisPes
   real(kind=8),dimension(6):: geom,derivs
   real(kind=8):: v
   ! Run
   call thisPes%initialize()
   write(*,*) 'Pes test with double initialization'
   write(*,*) '**************************************'
   geom(:)=[0.d0,0.d0,-231.5d0,1.42d0,0.d0,0.d0]
   write(*,*) 'Geom: ',geom(:)
   call thisPes%get_v_and_derivs(geom,v,derivs)
   write(*,*) 'Pot: ',v
   write(*,*) 'Derivs: ',derivs(:)
   write(*,*) '**************************************'
   geom(:)=[0.d0,0.d0,1.5d0,323.42d0,0.d0,0.d0]
   write(*,*) 'Geom: ',geom(:)
   call thisPes%get_v_and_derivs(geom,v,derivs)
   write(*,*) 'Pot: ',v
   write(*,*) 'Derivs: ',derivs(:)
   write(*,*) '**************************************'
   geom(:)=[0.d0,0.d0,4.5d0,1.42d0,0.d0,0.d0]
   write(*,*) 'Geom: ',geom(:)
   call thisPes%get_v_and_derivs(geom,v,derivs)
   write(*,*) 'Pot: ',v
   write(*,*) 'Derivs: ',derivs(:)
end program patata
