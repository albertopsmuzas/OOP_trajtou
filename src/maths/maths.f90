!#########################################################
! MODULE MATHS
!> @brief
!! Contains some mathematical operations
!##########################################################
MODULE MATHS_MOD
! Initial declarations
IMPLICIT NONE
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
!> @param[in] nelem - Number of elements in @b arreg1 and @b arreg2
!> @param[in,out] arreg1 - This array is ordered so that @f$v_{i}<v_{i+1}@f$
!> @param[in,out] arreg2 - This array is ordered so that @f$w_{i}\spadesuit w_{i+1}@f$
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 27/Jan/2014
!> @version 1.0
! -----------------------------------------------------
SUBROUTINE ORDER(NELEM,ARREG1,ARREG2)
   IMPLICIT NONE
   ! I/O variables ---------------------------
   INTEGER, INTENT(IN) :: NELEM
   REAL*8, DIMENSION(NELEM), INTENT(INOUT) :: ARREG1, ARREG2
   ! Local variables -------------------------
   LOGICAL :: CLAVE
   INTEGER :: I, K, POS
   REAL*8 :: AUX1, AUX2
   ! FIRE IN THE HOLE! ------------------------
   IF (NELEM.LT.2)  THEN
      WRITE(0,*) "ORDER ERR: Less than 2 elements inside arrays."
      STOP
   END IF
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
	CALL ORDER(n,x,v) ! Order the initial array
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
	CALL ORDER(nsym,xsym,vsym)
	WRITE(*,*) "SYMMETRIZE: New points added."
	RETURN
END SUBROUTINE SYMMETRIZE

END MODULE MATHS_MOD