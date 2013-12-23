!#########################################################
! RAISON D'ÊTRE:
! - Interpolations for Z variable (XY fixed) 
! UPDATES:
! - Created 19/Dec/2013
! FUNCTIONALITY:
! - 
! IDEAS FOR THE FUTURE:
! - None
!##########################################################
MODULE ZINTERPOL_MOD
! Initial declarations
IMPLICIT NONE
!=====================================================================
! Type Polinom3
!---------------
TYPE, PRIVATE :: Polinom3
   REAL(KIND=8) :: a, b, c, d
END TYPE Polinom3
!=====================================================================
! Type Generic interpolation in Z
!---------------------------------
TYPE, PRIVATE :: Zinterpol
END TYPE Zinterpol
!====================================================================
! Type Cubic Splines
!--------------------
TYPE, EXTENDS(Zinterpol) :: Z_csplines
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: d2vdz ! collection of second derivatives at each point
   TYPE(Polinom3), DIMENSION(:), ALLOCATABLE :: coeff ! coefficients
END TYPE Z_csplines
CONTAINS
!######################################################################
! SUBROUTINE: DSPLIN ##################################################
!######################################################################
! For cubic splines interpolation.
! Si(Xi)=ai*(xi**3)+bi*(xi**2)+ci*xi+di
! Calculates S''i(xi) [ second derivatives at nodes ]
! Splines satisfy continuity , smothness and second derivatives cointinuity.
! There should be two extra conditions to get a  system of compatible
! linear equations. "id1" and "id2" define the kind of conditions that
! are supported by this dspline version. "id1" is the code for "cond1" (first node condition)
! and  "id2" for "cond2" (last node condition).
!============
! id values: |
!============
!-  0: only supported for id1 (cond1). First two intervals have the same cubic spline.
!-  1: S''(x(1))=cond1 if id1=1  or S''(x(n))=cond2 if id2=1
!      n is the number of nodes. There are n-1 cubic splines.
!----------------------------------------------------------------------


END MODULE ZINTERPOL_MOD

