!###########################################################
!# SUBROUTINE: MCTDH_readInputCRP6D
!###########################################################
!> @brief
!! Wrapper subroutine to read a typical CRP6D input file.
!
!> @warning
!! - It will need the static library generated by OOPtrajtou.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Oct/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE MCTDH_readInputCRP6D(crp6d_pes)
   ! Initial declarations
   USE CRP6D_MOD
   IMPLICIT NONE
   ! I/O variables
   TYPE(CRP6D),INTENT(OUT) :: crp6d_pes
   ! Run section
   ! Change coordinates
   CALL crp6d_pes%INITIALIZE("INcrp6d.inp")
   RETURN
END SUBROUTINE MCTDH_readInputCRP6D
!###########################################################
!# SUBROUTINE: MCTDH_GetPotCRP6D
!###########################################################
!> @brief
!! Wrapper subroutine to get the CRP6D potential at a given point
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
SUBROUTINE MCTDH_getPotCRP6D(crp6d_pes,r,v)
   ! Initial declarations
   USE CRP6D_MOD
   IMPLICIT NONE
   ! I/O variables
   TYPE(CRP6D),INTENT(IN):: crp6d_pes
   REAL(KIND=8),DIMENSION(6),INTENT(IN):: r
   REAL(KIND=8),INTENT(OUT):: v
   ! Local variables
   REAL(KIND=8),DIMENSION(6):: geom,dvdu ! auxiliar geometry, dummy derivatives
   ! Run section
   ! Change of coordinates:
   geom(3:6)=r(3:6)
   geom(1:2) = crp6d_pes%surf%surfunit2cart(r(1:2))
   call crp6d_pes%GET_V_AND_DERIVS(geom,v,dvdu)
   RETURN
END SUBROUTINE MCTDH_getPotCRP6D
!###########################################################
!# SUBROUTINE: MCTDH_evalPotCRP6D
!###########################################################
!> @brief
!! Calls readInput and getPot subroutines
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Oct/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE crp6d_generic(r,v)
   ! Initial declarations
   USE CRP6D_MOD
   IMPLICIT NONE
   ! I/O variables
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: r
   REAL(KIND=8),INTENT(IN) :: v
   ! Local variables
   TYPE(CRP6D) :: crp6dPes
   ! Run section
   CALL MCTDH_readInputCRP6D(crp6dPes)
   CALL MCTDH_getPotCRP6D(crp6dPes,r,v)
   RETURN
END SUBROUTINE crp6d_generic
