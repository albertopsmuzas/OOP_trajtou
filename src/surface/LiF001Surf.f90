!#########################################################
! MODULE SURFACE_MOD
!> @brief
!! Specific implementation of module SURFACE_MOD
!##########################################################
module LiF001SURF_MOD
use SURFACE_MOD
implicit none
!/////////////////////////////////////////////////////////////////////////////////////
! TYPE: LiF001SURF
!> @brief
!! Specific implementation for LiF001 extending Surface Type
!-------------------------------------------------------------------------------------
type,extends(Surface):: LiF001Surf
	contains
	   procedure,public:: initialize => initialize_LiF001
end type LiF001Surf

contains
!###############################################################################
!# SUBROUTINE: INITIALIZE_LiF001 ###############################################
!###############################################################################
!> @brief
!! Specific implementation of initialize
!-------------------------------------------------------------------------------
subroutine initialize_LiF001(surf,filename)
   implicit none
   ! I/O Variables -----------------------------------------------
   class(LiF001Surf),intent(inout):: surf
   character(len=*),intent(in):: filename
   ! Run section --------------------------------------------
   surf%alias='LiF001'
   surf%diff_atoms=2
   surf%symmLabel='p4mm'
   surf%s1=[5.4433561257770959d0,0.d0]
   surf%s2=[0.d0,5.4433561257770959d0]
   allocate( surf%atomType(2) )
   surf%units='au'
   surf%norm_s1=5.4433561257770959d0
   surf%norm_s2=5.4433561257770959d0
   ! Set basis vectors
   surf%surf2cart_mtrx(:,1)=surf%s1(:)
   surf%surf2cart_mtrx(:,2)=surf%s2(:)
   ! Set Lithium specifications
   surf%atomtype(1)%alias='Li'
   surf%atomtype(1)%n=1
   allocate(surf%atomtype(1)%atom(1,3))
   surf%atomtype(1)%atom(1,:)=[0.d0,0.d0,-0.123593758269d0]
   ! Set Fluorinne specifications
   surf%atomtype(2)%alias='F'
   surf%atomtype(2)%n=1
   allocate(surf%atomtype(2)%atom(1,3))
   surf%atomtype(2)%atom(1,:)=[2.72167806289d0,2.72167806289d0,0.d0]
   ! Set conversion matrices
   surf%metricsurf_mtrx=matmul(transpose(surf%surf2cart_mtrx),surf%surf2cart_mtrx)
   call inv_mtrx(2,surf%surf2cart_mtrx,surf%cart2surf_mtrx)
   surf%surfunit2cart_mtrx(:,1)=surf%surf2cart_mtrx(:,1)/surf%norm_s1
   surf%surfunit2cart_mtrx(:,2)=surf%surf2cart_mtrx(:,2)/surf%norm_s2
   call inv_mtrx(2,surf%surfunit2cart_mtrx,surf%cart2surfunit_mtrx)
   surf%recip2cart_mtrx=transpose(surf%cart2surf_mtrx)
   surf%recip2cart_mtrx=2.D0*pi*surf%recip2cart_mtrx
   call inv_mtrx(2,surf%recip2cart_mtrx,surf%cart2recip_mtrx)
   surf%initialized=.true.
   return
end subroutine initialize_LiF001

end module LiF001SURF_MOD
