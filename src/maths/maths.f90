!#########################################################
! MODULE MATHS
!> @brief
!! Contains some mathematical operations
!##########################################################
MODULE MATHS_MOD
! Initial declarations
implicit none
CONTAINS
!####################################################################
! SUBROUTINE: ORDER #################################################
!####################################################################
!> @brief
!! This subroutine orders an array using the insertion ordering algorithm while
!! a second array is ordered following the same permutations as the prior one.
!
!> @details
!! Let's define some notation:
!! - @b arregl1 has elements @f$v_{i}@f$ 
!! - @b arregl2 has elements @f$w_{i}@f$
!! - If @f$<@f$ and @f$\spadesuit@f$ are order relationships, it holds:
!!   @f$w_{i}\spadesuit w_{i+1} \Leftrightarrow v_{i}<v_{i+1}@f$
!
!> @param[in,out] arreg1 - This array is ordered so that @f$v_{i}<v_{i+1}@f$
!> @param[in,out] arreg2 - This array is ordered so that @f$w_{i}\spadesuit w_{i+1}@f$
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 27/Jan/2014
!> @version 1.0
! -----------------------------------------------------
SUBROUTINE ORDER(ARREG1,ARREG2)
   IMPLICIT NONE
   ! I/O variables ---------------------------
   REAL*8, DIMENSION(:), INTENT(INOUT) :: ARREG1, ARREG2
   ! Local variables -------------------------
   INTEGER :: NELEM
   LOGICAL :: CLAVE
   INTEGER :: I, K, POS
   REAL*8 :: AUX1, AUX2
   ! FIRE IN THE HOLE! ------------------------
   nelem=size(arreg1)
   IF (nelem/=size(arreg2)) THEN
      WRITE(0,*) "ORDER ERR: dimension mismatch in arreg1 and arreg2"
      WRITE(0,*) "arreg1: ", nelem
      WRITE(0,*) "arreg2: ", size(arreg2)
      CALL EXIT(1)
   END IF
   !
   SELECT CASE (nelem)
      CASE(: 1)
         WRITE(0,*) "ORDER ERR: Less than 2 elements inside arrays."
         CALL EXIT(1)
      CASE DEFAULT
         ! do nothing
   END SELECT 
   !
   DO I=2,NELEM
      K = I
      AUX1 = ARREG1(K)
      AUX2 = ARREG2(K)
      CLAVE = .FALSE.
      DO WHILE((K.GT.1).AND.(.NOT.CLAVE))
         IF (ARREG1(K-1).GT.AUX1) THEN
            ARREG1(K) = ARREG1(K-1)
            ARREG2(K) = ARREG2(K-1)
            K = K-1
         ELSE
            CLAVE = .TRUE.
         ENDIF
      END DO
   POS = K
   ARREG1(POS) = AUX1
   ARREG2(POS) = AUX2
   ENDDO
END SUBROUTINE ORDER
!###########################################################
!# SUBROUTINE: ORDER_VECT 
!###########################################################
!> @brief
!! This subroutine orders an order 1 array from low to high values
!> @param[in,out] arreg1 - Vector to order
!!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 31/Jan/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE ORDER_VECT(ARREG1)
	IMPLICIT NONE
	! I/O variables ---------------------------
	REAL(KIND=8),DIMENSION(:),INTENT(INOUT) :: ARREG1
   ! Local variables -------------------------
	INTEGER(KIND=4) :: NELEM
	LOGICAL :: CLAVE
	INTEGER :: I, K, POS
	REAL*8 :: AUX1
	! FIRE IN THE HOLE! ------------------------
   nelem=size(arreg1)
	IF (NELEM.LT.2)  THEN
		WRITE(0,*) "ORDER ERR: Less than 2 elements inside arrays."
		STOP
	END IF
	DO I=2,NELEM
		K = I
		AUX1 = ARREG1(K)
		CLAVE = .FALSE.
		DO WHILE((K.GT.1).AND.(.NOT.CLAVE))
			IF (ARREG1(K-1).GT.AUX1) THEN
				ARREG1(K) = ARREG1(K-1)
				K = K-1
			ELSE
				CLAVE = .TRUE.
			ENDIF
		END DO
	POS = K
	ARREG1(POS) = AUX1
	ENDDO
	RETURN
END SUBROUTINE ORDER_VECT
!####################################################################
! SUBROUTINE: SYMMETRIZE ############################################
!####################################################################
!> @brief
!! This subroutine takes  couples of values @f$(x_{i},F(x_{i}))@f$ and
!! makes them symmetric respect to a given value @f$x_{0}@f$. See details.
!
!> @details
!> - Actually, this routine projects points around a given value @f$x_{0}@f$ so
!! that the final distribution of points is symmetric respect that value. 
!! - This procedure is only done if @f$F(x_{i})>f_{0}@f$, where @f$f_{0}@f$ is a 
!! threshold value given to the routine.
!
!> @param[in] n - Number of initial @f$(x_{i},F(x_{i}))@f$ couples
!> @param[in,out] x,v - Couples @f$(x_{i},F(x_{i}))@f$
!> @param{in} zero - @f$x_{0}@f$ value
!> @param[in] vtop - @f$f_{0}@f$ threshold value
!> @param[out] nsym - Number of final points after symmetrization procedure
!> @param[out] xsym, vsym - Collection of @f$(x_{i},F(x_{i}))@f$ couples after
!!                          symmetrization procedure
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 03/Feb/2014
!> @version 1.0
! -----------------------------------------------------
SUBROUTINE SYMMETRIZE(n,x,v,zero,vtop,nsym,xsym,vsym)
	IMPLICIT NONE
	! I/O variables ------------------------
	INTEGER, INTENT(IN) :: n
	INTEGER, INTENT(OUT) :: nsym
	REAL*8, INTENT(IN) :: vtop, zero
	REAL*8, INTENT(INOUT), DIMENSION(n) :: x, v 
	REAL*8, INTENT(OUT), DIMENSION(:), ALLOCATABLE :: xsym, vsym
	! Internal variables -------------------
	INTEGER :: i,j, k ! Counters
	INTEGER :: npoints, nnonred
	REAL*8, DIMENSION(:), ALLOCATABLE :: auxx, auxv 
	REAL*8, DIMENSION(:), ALLOCATABLE :: nonredx, nonredv 
	INTEGER, DIMENSION(:), ALLOCATABLE :: redundant
	! HEY, HO! LET'S GO! -------------------
	CALL ORDER(x,v) ! Order the initial array
	npoints=0
	! Checking number of points that satisfies F(x) > vtop
	DO i=1,n
		IF (v(i).GE.vtop) THEN
			WRITE(*,*) "SYMMETRIZE: Pair: ",i, x(i), v(i)
			npoints=npoints+1
		END IF
	END DO
	WRITE(*,*) "SYMMETRIZE: Symmetrize subroutine was invoked. Take care of what you are doing."
	WRITE(*,*) "SYMMETRIZE: Found ", npoints, "points to symmetrize. F(X) > ", vtop
	! Check if some of these points are redundant
	ALLOCATE(auxx(1:npoints))
	ALLOCATE(auxv(1:npoints))
	ALLOCATE(redundant(1:npoints))
	DO i=1, npoints
		redundant(i)=0
	END DO
	k=0
	DO i=1,n
		IF (v(i).GE.vtop) THEN
			k=k+1
			auxv(k)=v(i)
			auxx(k)=x(i)
		END IF
	END DO
	k=0
	DO i=1, npoints
		IF (redundant(i).EQ.0) THEN
			DO j=i+1, npoints
				IF (DABS(auxx(i)-zero).EQ.DABS(auxx(j)-zero)) THEN
					k=k+1
					redundant(i)=k
					redundant(j)=k
				END IF
			END DO
		END IF
	END DO
	WRITE(*,*) "SYMMETRIZE: Redundance mapping:"
	DO i=1, npoints
		WRITE(*,*) i, auxx(i), auxv(i), redundant(i)
	END DO
	
!debug	WRITE(*,*) "SYMMETRIZE: ", k, "redundant pairs found."
	! Store non-redundant points
	nnonred=npoints-k*2
	ALLOCATE(nonredx(1:nnonred))
	ALLOCATE(nonredv(1:nnonred))
	k=0
	DO i=1, npoints
		IF (redundant(i).EQ.0) THEN
			k=k+1
			nonredx(k)=auxx(i)
			nonredv(k)=auxv(i)
		END IF
	END DO
	! Adding the new points 
	nsym=n+nnonred
	ALLOCATE(xsym(1:nsym))
	ALLOCATE(vsym(1:nsym))
	k=0
	DO i=1,nnonred 
		IF (nonredx(i)-zero.LE.0) THEN
			k=k+1
			xsym(k)=nonredx(i)+2.D0*DABS(nonredx(i)-zero)
			vsym(k)=nonredv(i)
		ELSE IF (nonredx(i)-zero.GT.0) THEN
			k=k+1
			xsym(k)=nonredx(i)-2.D0*DABS(nonredx(i)-zero)
			vsym(k)=nonredv(i)
		END IF
	END DO
	DO i=nnonred+1,nsym
		xsym(i)=x(i-nnonred)
		vsym(i)=v(i-nnonred)
	END DO
	! The vector is not ordered
	CALL ORDER(xsym,vsym)
	WRITE(*,*) "SYMMETRIZE: New points added."
	RETURN
END SUBROUTINE SYMMETRIZE
!###########################################################################
! SUBROUTINE: TRIDIA #######################################################
!###########################################################################
!> @brief
!! Solves a tridiagonal system. Thomas algorithm.
!! Conserve the system matrix.
!
!> @param[in] n dimensions of the arrays. number of equations
!> @param[in] a(n) - subdiagonal
!> @param[in] b(n) - diagonal
!> @param[in] c(n) - superdiagonal
!> @param[in] d(n) - right part of the equation
!> @param[out] x(n) - answer
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 11/Feb/2014
!> @version 1.0
! -------------------------------------------------------------------
SUBROUTINE TRIDIA(n,a,b,c,d,x)
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: n
        REAL*8,DIMENSION(n),INTENT(IN) :: a,b,c,d
        REAL*8,DIMENSION(n),INTENT(OUT) :: x
        REAL*8,DIMENSION(n) :: cp,dp
        REAL*8 :: m
        INTEGER :: i
! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
	DO i = 2,n
		m = b(i)-cp(i-1)*a(i)
		cp(i) = c(i)/m
		dp(i) = (d(i)-dp(i-1)*a(i))/m
	END DO
! initialize x
	x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
	DO i = n-1, 1, -1
		x(i) = dp(i)-cp(i)*x(i+1)
	END DO
END SUBROUTINE TRIDIA
!#####################################################################
! SUBROUTINE : INV_MTRX ##############################################
!#####################################################################
!> @brief
!! Invert a square matrix with Gauss method. @b mtrx is not destroyed during  the procedure
!
!> @param[in] n - Order of the matrix
!> @param[in] mtrx(n,n) - Initial matrix
!> @param[out] i_mtrx(n,n) - Inverse matrix
!
!> @warning
!! - This algorithm is inefficient for large matrices or sparse matrices. Use
!!   it wisely
!---------------------------------------------------------------------
SUBROUTINE INV_MTRX (n, mtrx, i_mtrx)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n
	REAL(8), DIMENSION(n,n), INTENT(IN):: mtrx
	REAL(8), DIMENSION(n,n), INTENT(OUT):: i_mtrx
	! Local Variables ------------------------
	REAL*8, DIMENSION(n,n) :: B, A
	REAL*8, DIMENSION(n) :: temp
	INTEGER, DIMENSION(n) :: ipvt
	INTEGER, DIMENSION(1) :: imax
	REAL(8) :: c, d
	INTEGER :: i, j, k, m
	! HEY, HO! LET'S GO! ---------------------
	A = mtrx
	B = A
	ipvt = (/ (i, i = 1, n) /)
	DO k = 1,n
		imax = MAXLOC(ABS(b(k:n,k)))
		m = k-1+imax(1)
		IF (m /= k) THEN
			ipvt( (/m,k/) ) = ipvt( (/k,m/) )
			B((/m,k/),:) = B((/k,m/),:)
		END IF
		d = 1/B(k,k)
		temp = B(:,k)
		DO j = 1, n
			c = B(k,j)*d
			B(:,j) = B(:,j)-temp*c
			B(k,j) = c
		END DO
		B(:,k) = temp*(-d)
		B(k,k) = d
	END DO
	A(:,ipvt) = B
	i_mtrx = A
END SUBROUTINE INV_MTRX
!###################################################################
! FUNCTION: cartesianPeriodizer
!###################################################################
!> @brief
!! Given a function f(x), and the cartesian periodizer function g(x,T),
!! it stands that f( g(x,T) ) is the f(x)'s segment between 0<=x<=T repeated
!! indefinitely
!-------------------------------------------------------------------
function cartesianPeriodizer(x,T) result(y)
   ! initial declarations
   implicit none
   ! I/O variables
   real(kind=8),intent(in):: x, T
   ! Dummy function variable
   real(kind=8):: y
   ! Local variables
   real(kind=8),parameter:: pi=dacos(-1.d0)
   ! Run section .....................................
   y=T/2.d0-(T/pi)*datan( 1.d0/dtan(x*pi/T) )
   return
end function cartesianPeriodizer
!###################################################################
! FUNCTION: polarPeriodizer
!###################################################################
!> @brief
!! Defined as cartesian periodizer function g(x,T) with some changes:
!! x=x-x0 and T=2pi/N. Designed to be used with figures with polar
!! symmetry and whole number of lobes.
!> @details
!! - Extracted from E. Chicurel-Uziel/Computer Aided Geometric Design/21(2004) 23-42
!-------------------------------------------------------------------
function polarPeriodizer(theta,N,theta0) result(y)
   ! initial declarations
   implicit none
   ! I/O variables
   real(kind=8),intent(in):: theta,theta0
   integer(kind=4),intent(in):: N
   ! dummy function out variable
   real(kind=8):: y
   ! Local variables
   real(kind=8),parameter:: pi=dacos(-1.d0)
   ! Run section ..............................
   y=pi/dfloat(N)-(2.d0/dfloat(N))*datan( 1.d0/dtan(dfloat(N)*(theta-theta0)/2.d0) )
   return
end function polarPeriodizer
!####################################################################
! FUNCTION: radialPolygonEquation
!####################################################################
!> @brief
!! 
!--------------------------------------------------------------------
function radialPolygonEquation(theta,r,N,theta0) result(y)
   ! initial declarations
   implicit none
   ! I/O variables
   real(kind=8),intent(in):: theta,r,theta0
   integer(kind=4),intent(in):: N
   ! Dummy output variable
   real(kind=8):: y
   ! Local variable
   real(kind=8):: beta
   real(kind=8),parameter:: pi=dacos(-1.d0)
   ! Run section ......................................
   beta=pi*(dfloat(N-2)/dfloat(2*N))
   y=r*datan(beta)/(dsin(polarPeriodizer(theta,N,theta0))+dtan(beta)*dcos(polarPeriodizer(theta,N,theta0)))
   return
end function radialPolygonEquation
!####################################################################
! FUNCTION: parametricPolygonEquation
!####################################################################
!> @brief
!! 
!--------------------------------------------------------------------
subroutine parametricPolygonEquation(theta,r,N,theta0,x0,x)
   ! initial declarations
   implicit none
   ! I/O variables
   real(kind=8),intent(in):: theta,theta0,r
   real(kind=8),dimension(2),intent(in):: x0
   integer(kind=4),intent(in):: N
   real(kind=8),dimension(2),intent(out):: x
   ! Run section ...................................
   x(1)=x0(1)+radialPolygonEquation(theta,r,N,theta0)*dcos(theta)
   x(2)=x0(2)+radialPolygonEquation(theta,r,N,theta0)*dsin(theta)
   return
end subroutine parametricPolygonEquation
!##################################################################
! FUNCTION: checkLoschianOrder
!##################################################################
!> @brief
!! Given an integer number, check the order of the nearest Loschian
!! number. It will be used to give diffraction order for hexagonal
!! lattices.
!> @details
!! - It uses an array of first 64 Loschian numbers. If diffraction is
!!   too high, it may be insufficient and a longer series should be given
!!   and compiled. You can detect this unconviniency if you get too many
!!   peaks with order 63 and delta K changes within this order
!------------------------------------------------------------------
function checkLoschianOrder(num) result(order)
   implicit none
   ! I/O variables
   integer(kind=4),intent(in):: num
   ! Dummy function output variable
   integer(kind=4):: order
   ! Local variables
   integer(kind=4),dimension(0:63),parameter:: loschianNum=[   0,  1,  3,  4,  7,  9, 12, 13, 16, 19, 21, 25, 27, 28, 31, 36,&
                                                              37, 39, 43, 48, 49, 52, 57, 61, 63, 64, 67, 73, 75, 76, 79, 81,&
                                                              84, 91, 93, 97,100,103,108,109,111,112,117,121,124,127,129,133,&
                                                             139,144,147,148,151,156,157,163,169,171,172,175,181,183,189,192]
   integer(kind=4),dimension(0:63):: intDistance
   ! Run section ...............................................................................................................
   intDistance(:)=( loschianNum(:)-num )**2.d0
   order=minloc( array=intDistance(:),dim=1 )
   return
end function checkLoschianOrder
END MODULE MATHS_MOD
