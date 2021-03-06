!###########################################################
!# SUBROUTINE: MCTDH_readInputCRP3D
!###########################################################
!> @brief
!! Wrapper subroutine to read a typical CRP3D input file.
!
!> @warning
!! - It will need the static library generated by OOPtrajtou.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Oct/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE MCTDH_readInputCRP3D(crp3d_pes)
   ! Initial declarations
   USE CRP3D_MOD
   IMPLICIT NONE
   ! I/O variables
   TYPE(CRP3D),INTENT(OUT) :: crp3d_pes
   ! Run section
   CALL crp3d_pes%INITIALIZE("INcrp3d.inp")
   RETURN
END SUBROUTINE MCTDH_readInputCRP3D
!###########################################################
!# SUBROUTINE: MCTDH_GetPotCRP3D
!###########################################################
!> @brief
!! Wrapper subroutine to get the CRP3D potential at a given point
!! r in space
!
!> @warning
!! - It will need the static library generated by OOPtrajtou.
!! - It is assumed that the potential is evaluated at a point 'r', which
!!   is given in a.u. normalized surface coordinates, i.e. in atomic units
!!   along the direction of surface vectors (standard way of MCTDH)
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Oct/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE MCTDH_getPotCRP3D(crp3d_pes,r,v)
   ! Initial declarations
   USE CRP3D_MOD
   IMPLICIT NONE
   ! I/O variables
   TYPE(CRP3D),INTENT(IN) :: crp3d_pes
   REAL(KIND=8),DIMENSION(3),INTENT(IN) :: r
   REAL(KIND=8),INTENT(OUT) :: v
   ! Run section
   REAL(KIND=8),DIMENSION(3):: geom,dvdu ! auxiliar geometry, dummy derivatives
   ! Change of coordinates:
   geom(3) = r(3)
   geom(1:2) = crp3d_pes%surf%surfunit2cart(r(1:2))
   call crp3d_pes%GET_V_AND_DERIVS(geom,v,dvdu)
   RETURN
END SUBROUTINE MCTDH_getPotCRP3D
!###########################################################
!# SUBROUTINE: MCTDH_evalPotCRP3D
!###########################################################
!> @brief
!! Calls readInput and getPot subroutines
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Oct/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE crp3d_generic(r,v)
   ! Initial declarations
   USE CRP3D_MOD
   IMPLICIT NONE
   ! I/O variables
   REAL(KIND=8),DIMENSION(3),INTENT(IN) :: r
   REAL(KIND=8),INTENT(IN) :: v
   ! Local variables
   TYPE(CRP3D) :: crp3dPes
   ! Run section
   CALL MCTDH_readInputCRP3D(crp3dPes)
   CALL MCTDH_getPotCRP3D(crp3dPes,r,v)
   RETURN
END SUBROUTINE crp3d_generic
