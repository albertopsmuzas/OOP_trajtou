!###############################################
! MODULE: INITATOM_MOD
!> @brief
!! This module provides routines and onjects to create
!! initial conditions for an atom or list of atoms
!###############################################
MODULE INITATOM_MOD
IMPLICIT NONE
TYPE,EXTENDS(Inicond) :: Inicond_atoms
   LOGICAL :: control_vel, control_posX, control_posY, control_out, control_seed
   REAL*8 :: impact_x, impact_y
   TYPE(Energy) :: E_norm
   TYPE(Velocity) :: vz_angle, vpar_angle
   TYPE(Length) :: init_z ! initial Z value
   TYPE(Mass) :: masss
END TYPE Inicond_atoms
!/////////////////////////////////////////////////////
! TYPE: Atoms
!> @brief
!! Atom subtype dynamics object
!----------------------------------------------------
TYPE,EXTENDS(Dynobject) ::  Atom
   INTEGER :: ireb=0 ! times this atom has changed Pz's direction
   INTEGER :: ixyboun=0 ! times this atoms has changed parallel momentum's direction
END TYPE Atom
!/////////////////////////////////////////////////////
! TYPE: Atom_trajs
!> @brief
!! A list of atom subtype dynamics objects
!----------------------------------------------------
TYPE :: Atom_trajs
   TYPE(Atom),DIMENSION(:),ALLOCATABLE :: atomo
END TYPE Atom_trajs

END MODULE INITATOM_MOD
