!#############################################################
! MODULE : FOURIER_MOD
!
!> @brief 
!! Allows to create 2d fourier interpolations of data
!
!> @warning
!! - Only for C4v symmetry
!
!> @todo
!! - Add compatibility woth other symmetries
!! - Use planar symmetry notation and not puntual one
!#############################################################
MODULE FOURIER2D_MOD
   USE INTERPOL2D_MOD
IMPLICIT NONE
!////////////////////////////////////////////////////////////
! TYPE: FOURIER2D
!> @brief
!! Contains all information and procedures to perform a fourier symmetry addapted
!! 2D interpolation
!
!> @param coeff - Array of coefficients for @b f function
!> @param coeff_dfdz - Array of coefficients for @b dfdz function
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 06/Feb/2014
!> @version 1.0
!-----------------------------------------------------------
TYPE,EXTENDS(Interpol2d) :: Fourier2d
   PRIVATE
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: coeff
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: coeff_dfdz
CONTAINS
   PROCEDURE,PUBLIC :: SET_COEFF => SET_FOURIER2D_COEFF
   PROCEDURE,PUBLIC :: GET_F_AND_DERIV => GET_F_AND_DFDZ_FOURIER2D
END TYPE Fourier2d
!////////////////////////////////////////////////////////////
CONTAINS
!#############################################################
!# SUBROUTINE: SET_FOURIER_COEFF #############################
!#############################################################
!> @brief
!! Calculates coefficients for a fourier 2d expansion
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 06/Feb/2014
!> @version 1.0
!--------------------------------------------------------------
SUBROUTINE SET_FOURIER2D_COEFF(this,surf)
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   USE SURFACE_MOD
   USE LAPACKCONTROL_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourier2d),INTENT(INOUT) :: this
   TYPE(Surface),INTENT(IN) :: surf
   ! Local variables
   INTEGER :: i,j,n,k ! counter
   INTEGER :: nn, max_n
   REAL*8, DIMENSION(this%n,this%n) :: T, inv_T 
   REAL*8, DIMENSION(2) :: r ! normalized surface coordinates
   CHARACTER(LEN=21),PARAMETER :: routinename="SET_FOURIER2D_COEFF: "
   ! Lapack variables
   INTEGER(KIND=4) :: lwork, info
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: work
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE :: ipiv
	! HEY HO!, LET'S GO! ----------
   ALLOCATE(this%coeff(this%n))
	! Calculate number of terms 
	IF (this%n.EQ.1) THEN
		max_n = 0
	ELSE IF (this%n.EQ.3) THEN
		max_n=1
	ELSE IF (this%n.GT.3) THEN
		! Here, we are supposing that we are using the exact amout of sites that
		! are necessary for a Nth order fourier expansion. See documentation.
		max_n = int((-3.D0+dsqrt(1.D0+8.D0*DFLOAT(this%n)))/2.D0)
	END IF
	! Create Y vector & T matrix
	DO i=1, this%n
		! site data
		r(1)=this%xy(i,1)
		r(2)=this%xy(i,2)
#ifdef DEBUG
		CALL DEBUG_WRITE(routinename,"y(i)",this%f(i))
		CALL DEBUG_WRITE(routinename,"r(1)",r(1))
		CALL DEBUG_WRITE(routinename,"r(2)",r(2))
		CALL DEBUG_WRITE(routinename,"==================================")
		CALL DEBUG_WRITE(routinename,"i=", i)
		CALL DEBUG_WRITE(routinename,"j --->",1)
#endif
		T(i,1) = 1.D0
#ifdef DEBUG
		CALL DEBUG_WRITE(routinename,"g_plus value: ",G_PLUS(surf,0,0,r))
#endif
		DO j=2,this%n
#ifdef DEBUG
			CALL DEBUG_WRITE(routinename,"j --->", j)
#endif
			DO n=1, max_n
				nn=n*(n+1)/2
#ifdef DEBUG
				CALL DEBUG_WRITE(routinename,"     n --->", n)
				CALL DEBUG_WRITE(routinename,"nn", nn)
#endif
				IF (j.EQ.nn+1) THEN
					T(i,j)=G_PLUS(surf,n,0,r)+G_PLUS(surf,0,n,r)
#ifdef DEBUG
					CALL DEBUG_WRITE(routinename,"       Bn, n= ", n)
					CALL DEBUG_WRITE(routinename,"g_plus term: ",G_PLUS(surf,n,0,r)+G_PLUS(surf,0,n,r))
#endif
					EXIT
				ELSE IF (j.EQ.nn+2) THEN
					T(i,j)=G_PLUS(surf,n,n,r)+G_PLUS(surf,-n,n,r)
#ifdef DEBUG
					CALL DEBUG_WRITE(routinename,"       Dn, n=", n)
					CALL DEBUG_WRITE(routinename,"g_plus term:", G_PLUS(surf,n,n,r)+G_PLUS(surf,-n,n,r))
#endif
					EXIT
				ELSE
					DO k=1, n-1 
#ifdef DEBUG
						CALL DEBUG_WRITE(routinename,"       k ---> ", k)
#endif
						IF (j.EQ.nn+2+k) THEN
							T(i,j)=G_PLUS(surf,k,n,r)+G_PLUS(surf,n,k,r)+G_PLUS(surf,-k,n,r)+G_PLUS(surf,n,-k,r)
#ifdef DEBUG
							CALL DEBUG_WRITE(routinename,"         Akn") 
							CALL DEBUG_WRITE(routinename,"k,n ---> ",k,n)
							CALL DEBUG_WRITE(routinename,"g_plus term:",G_PLUS(surf,k,n,r)+G_PLUS(surf,n,k,r)+G_PLUS(surf,-k,n,r)+G_PLUS(surf,n,-k,r))
#endif
							EXIT
						END IF
					END DO
#ifdef DEBUG
				CALL DEBUG_WRITE(routinename,"     Nothing for this n")
#endif            
				END IF
			END DO
		END DO
	END DO
   ! Lapack variables
   lwork=3
   ALLOCATE(work(lwork))
   ALLOCATE(ipiv(3))
   FORALL(i=1:3) ipiv(i)=i
   ! Calculate inverse matrix
   inv_T=T
   CALL DGETRI(3,inv_T,3,ipiv,work,lwork,info)
   CALL LAPACK_CHECK("DGETRI",info)
   ! Get coefficients
   this%coeff = matmul(inv_T,this%f)
   this%coeff_dfdz = matmul(inv_T,this%dfdz)
#ifdef DEBUG
	CALL DEBUG_WRITE(routinename,"Coeff f")
	DO i=1, this%n
		CALL DEBUG_WRITE(routinename,this%coeff(i))
	END DO
	CALL DEBUG_WRITE(routinename,"Coeff dfdz")
	DO i=1, this%n
		CALL DEBUG_WRITE(routinename,this%coeff_dfdz(i))
	END DO
#endif
	RETURN
END SUBROUTINE SET_FOURIER2D_COEFF
!################################################################
! FUNCTION: G_PLUS ##############################################
!################################################################
!> @brief
!! Interesting function to span fourier series
!
!> @warning
! - Input in cartesian coordinates. In C_{4v} symmetry, normalized
!   surface coordinates are equivalent to cartesian ones.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 06/Feb/2014
!> @version 1.0
!
!> @see extra documentation
!----------------------------------------------------------------
REAL*8 FUNCTION  g_plus(surf,h,k,r)
   USE SURFACE_MOD
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables
   REAL*8, DIMENSION(2), INTENT(IN) :: r
   TYPE(Surface), INTENT(IN) :: surf
   INTEGER, INTENT(IN) :: h, k
   ! Local variables
   REAL*8 :: g ! constants
   ! HEY HO!, LET'S GO! ---------------
   g=2.D0*PI/surf%norm_s1 ! 
   g_plus = DCOS(g*(DFLOAT(h)*r(1)+DFLOAT(k)*r(2)))
   RETURN
END FUNCTION g_plus
!################################################################
! FUNCTION: G_MINUS #############################################
!################################################################
!> @brief
!! Interesting function to span fourier series
!
!> @warning
!! - Input in cartesian coordinates. In C_{4v} symmetry, normalized
!!   surface coordinates are equivalent to cartesian ones.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 06/Feb/2014
!> @version 1.0
!
!> @see extra documentation
!----------------------------------------------------------------
REAL*8 FUNCTION  g_minus(surf,h,k,r)
   USE SURFACE_MOD
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables
   REAL*8, DIMENSION(2), INTENT(IN) :: r
   TYPE(Surface), INTENT(IN) :: surf
   INTEGER, INTENT(IN) :: h, k
   ! Local variables
   REAL*8 :: g ! constants
   ! HEY HO!, LET'S GO! ---------------
   g=2.D0*PI/surf%norm_s1
   g_minus = DSIN(g*(DFLOAT(h)*r(1)+DFLOAT(k)*r(2)))
   RETURN
END FUNCTION g_minus
!################################################################
! SUBROUTINE: GET_F_AND_DFDZ_FOURIER2D ##########################
!################################################################
!> @brief
!!  Given a Fourier inteprolation XY calculates the potential for a (x,y,z) 
!!  point and its derivatives in surface directions
!
!> @warning
!! - Only C4v symmetry allowed
!! - We should have already a set of coefficients defined.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 06/Feb/2014
!> @version 1.0
!----------------------------------------------------------------
SUBROUTINE GET_F_AND_DFDZ_FOURIER2D(interpolxy,surf,q,v,dvdu)
   USE SURFACE_MOD
   USE CONSTANTS_MOD
	IMPLICIT NONE
	! I/O variables
	CLASS(Fourier2d),INTENT(IN) :: interpolxy
	TYPE(Surface), INTENT(IN) :: surf
	REAL*8, DIMENSION(3), INTENT(IN) :: q
	REAL*8, DIMENSION(3), INTENT(OUT) :: dvdu
	REAL*8, INTENT(OUT) :: v
	! Local variables
	INTEGER :: i,n,k,nn, max_n ! counters
	REAL*8, DIMENSION(2) :: r
	REAL*8 :: g  ! in Cv4, modul_u1 = modul_u2
	! GABBA , GABBA HEY! ---------------------
	g = 2.D0*PI/surf%norm_s1 ! For this symmetry modul_u1 = modul_u2
	! Storing position
	r(1)=q(1)
	r(2)=q(2)
	! In C4v symmety, cartesian coordinates and normalized surface coordinates are the same.
	! Position has the correct value (inisde IWS cell, surface-normalized coordinates)
	IF (SIZE(interpolxy%coeff).EQ.1) THEN
		max_n = 0
	ELSE IF (SIZE(interpolxy%coeff).EQ.3) THEN
		max_n=1
	ELSE IF (SIZE(interpolxy%coeff).GT.3) THEN
		! Here, we are supposing that we are using the exact amout of sites that
		! are necessaty for a Nth order fourier expansion. See documentation.
		max_n = INT((-3.D0+DSQRT(1.D0+8.D0*DFLOAT(SIZE(interpolxy%coeff))))/2.D0)
	END IF
	v = 0.D0
	dvdu(1) = 0.D0
	dvdu(2) = 0.D0
	dvdu(3) = 0.D0
	v = v+interpolxy%coeff(1)*G_PLUS(surf,0,0,r) ! first term contribution to potential
	dvdu(3) = dvdu(3)+interpolxy%coeff_dfdz(1)*G_PLUS(surf,0,0,r) ! first term contribution to dvdz 
	DO i=2,SIZE(interpolxy%coeff)
		DO n=1, max_n
			nn=n*(n+1)/2
			IF (i.EQ.nn+1) THEN
				dvdu(1) = dvdu(1)+interpolxy%coeff(i)*(-g*DFLOAT(n)*G_MINUS(surf,n,0,r))
				dvdu(2) = dvdu(2)+interpolxy%coeff(i)*(-g*DFLOAT(n)*G_MINUS(surf,0,n,r))
				dvdu(3) = dvdu(3)+interpolxy%coeff_dfdz(i)*(G_PLUS(surf,n,0,r)+G_PLUS(surf,0,n,r))
				v = v+interpolxy%coeff(i)*(G_PLUS(surf,n,0,r)+G_PLUS(surf,0,n,r))
				EXIT
			ELSE IF (i.EQ.nn+2) THEN
				dvdu(1) = dvdu(1)+interpolxy%coeff(i)*(g*DFLOAT(n)*(G_MINUS(surf,-n,n,r)-G_MINUS(surf,n,n,r)))
				dvdu(2) = dvdu(2)+interpolxy%coeff(i)*(-g*DFLOAT(n)*(G_MINUS(surf,-n,n,r)+G_MINUS(surf,n,n,r)))
				dvdu(3) = dvdu(3)+interpolxy%coeff_dfdz(i)*(G_PLUS(surf,n,n,r)+G_PLUS(surf,-n,n,r))
				v = v+interpolxy%coeff(i)*(G_PLUS(surf,n,n,r)+G_PLUS(surf,-n,n,r))
				EXIT
			ELSE
				DO k=1, n-1 
					IF (i.EQ.nn+2+k) THEN
						dvdu(1) = dvdu(1)+interpolxy%coeff(i)*(g*DFLOAT(k)*(G_MINUS(surf,-k,n,r)-G_MINUS(surf,k,n,r))- &
						       g*DFLOAT(n)*(G_MINUS(surf,n,k,r)+G_MINUS(surf,n,-k,r)))
						dvdu(2) = dvdu(2)+interpolxy%coeff(i)*g*(DFLOAT(k)*(G_MINUS(surf,n,-k,r)- &
						       G_MINUS(surf,n,k,r))-DFLOAT(n)*(G_MINUS(surf,-k,n,r)+G_MINUS(surf,k,n,r)))
						dvdu(3) = dvdu(3)+interpolxy%coeff_dfdz(i)*(G_PLUS(surf,k,n,r)+G_PLUS(surf,n,k,r)+ &
						       G_PLUS(surf,-k,n,r)+G_PLUS(surf,n,-k,r))
						v = v+interpolxy%coeff(i)*(G_PLUS(surf,k,n,r)+G_PLUS(surf,n,k,r)+G_PLUS(surf,-k,n,r)+G_PLUS(surf,n,-k,r))
						EXIT
					END IF
				END DO
			END IF
		END DO
	END DO
END SUBROUTINE GET_F_AND_DFDZ_FOURIER2D
END MODULE FOURIER2D_MOD
