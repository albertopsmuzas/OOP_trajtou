!#########################################################
! MODULE: BULIRSCH_STOER_INTEGRATOR_MOD
!> @brief
!! Contains all subroutines to use Bulirsch-Stoer integrator with
!! time-step control.
!##########################################################
MODULE BULIRSCH_STOER_INTEGRATOR_MOD
! Initial declarations
USE INTEGRATOR_O1_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: BULIRSCH_STOER_INTEGRATOR
!> @brief
!! Type which stores info and subroutines needed for Bulirsch Stoer integration.
!
!> @details
!! - Routines addapted from FORTRAN  77 numerical recipes
!! - strinparam definition: 
!!   1) Kind of extrapolation: Polinomi, Rational (Case sensitive)
!! - realparam not used
!! - intparam not used.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!
!> @see FORTRAN 77 numerical recipes
!----------------------------------------------------------------
TYPE,EXTENDS(Integrator_o1) :: BULIRSCH_STOER_INTEGRATOR
CONTAINS
   PROCEDURE,PUBLIC :: STEP => STEP_BS
   PROCEDURE,PUBLIC :: MMID => MMID_BS
   PROCEDURE,PUBLIC :: POLINOM_EXTRAPOL => PZEXTR_BS
   PROCEDURE,PUBLIC :: RATIONAL_EXTRAPOL => RZEXTR_BS
   PROCEDURE,PUBLIC :: INTEGRATE => INTEGRATE_BS
END TYPE BULIRSCH_STOER_INTEGRATOR
!/////////////////////////////////////////////////////////////////
CONTAINS
!########################################################################
!# SUBROUTINE: INTEGRATE_BS #############################################
!########################################################################
!> @brief
!! Integrates first order differential equations with Bulirsch-Stoer
!! algorithm, with step control.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!------------------------------------------------------------------------
SUBROUTINE INTEGRATE_BS(this,f,dfdx,x,DERIV,switch)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Bulirsch_Stoer_integrator),INTENT(INOUT):: this
   REAL(KIND=8),DIMENSION(:),INTENT(INOUT):: f
   REAL(KIND=8),DIMENSION(:),INTENT(IN):: dfdx
   REAL(KIND=8),INTENT(INOUT):: x
   EXTERNAL DERIV
   LOGICAL,INTENT(INOUT) :: switch
   ! Local variables
   REAL(KIND=8) :: hdid,hnext
   ! Run section
   CALL this%STEP(f,dfdx,x,hdid,hnext,DERIV,switch)
   SELECT CASE(switch)
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         CALL this%SET_TIMESTEP_USED(hdid)
         CALL this%SET_TIMESTEP_NEXT(hnext)
   END SELECT
   RETURN
END SUBROUTINE INTEGRATE_BS
!#########################################################################################
!# SUBROUTINE: MMID_BS ###################################################################
!#########################################################################################
!> @brief 
!! Modified midpoint method. Given a vector of components @f$y_{i}(x)@f$, this subroutine
!! can advance @f$y_{i}(x+H)@f$ by a sequence of n substeps. It is used as an integrator in
!! the more powerful Bulirsch-Stoer technique.
!
!> @details
!! - Adapted from Numerical recipes FORTRAN 77
!! - Adapted for the integration of the equations of motion of an atom
!
!> @param[in] this - We need information stored in this variable
!> @param[in] y(1:6) - array with positions and momenta
!> @param[in] dydx - derivatives of y
!> @param[in] xs - X at which y and dydx are meassured
!> @param[in] htot - Total step size
!> @param[in] nstep - substeps to be used
!> @param[out] yout - output array
!> @param[in,out] switch - Controls if there was a problem calculating the potential at 
!!                         a given value
!
!> @see Fortran 77 numerical recipes
!-----------------------------------------------------------------------------------------
SUBROUTINE MMID_BS(this,y,dydx,xs,htot,nstep,yout,DERIVS,switch)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Bulirsch_Stoer_integrator),INTENT(IN):: this
   INTEGER(KIND=4),INTENT(IN):: nstep
   REAL(KIND=8),DIMENSION(:),INTENT(IN):: y,dydx
   REAL(KIND=8),DIMENSION(:),INTENT(OUT):: yout
   REAL(KIND=8),INTENT(IN):: xs,htot
   LOGICAL,INTENT(INOUT):: switch
   EXTERNAL:: DERIVS ! Subroutine to calculate first derivatives
   ! Local variables
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: ym,yn
   INTEGER(KIND=4):: i,n ! counters
   REAL(KIND=8):: h,h2,swap,x
   ! ROCK THE CASBAH !!! ---------------------
   h=htot/dfloat(nstep) ! Stepsize this trip.
   ALLOCATE(ym(this%getnv()))
   ALLOCATE(yn(this%getnv()))
   ym(:)=y(:)
   yn(:)=y(:)+h*dydx(:)
   x=xs+h
   CALL DERIVS(yn,yout,switch) ! Will use yout for temporary storage of derivatives.
   SELECT CASE(switch)
      CASE(.TRUE.)
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   h2=2.D0*h
   DO n=2,nstep ! General step.
      DO i=1,this%getnv()
         swap=ym(i)+h2*yout(i)
         ym(i)=yn(i)
         yn(i)=swap
      END DO
      x=x+h
      CALL DERIVS(yn,yout,switch)
      SELECT CASE(switch)
         CASE(.TRUE.)
            RETURN
         CASE(.FALSE.)
            ! do nothing
      END SELECT
   END DO
   yout(:)=0.5D0*(ym(:)+yn(:)+h*yout(:))
   RETURN
END SUBROUTINE MMID_BS
!###############################################################################################
!# SUBROUTINE : STEP_BS ########################################################################
!###############################################################################################
!> @brief
!! Bulirsch-Stoer step with monitoring of local truncation error to ensure accuracy and adjust
!! stepsize. 
!
!> @param[in,out] y(:) - @b real(kind=8) @b array. Function to integrate.
!!                       @f$ f(x), f:\mathbb{R}\rightarrow\mathbb{R}^{N} @f$
!> @param[in] dydx(:) - @b real(kind=8) @b array. Function derivatives.
!!                      @f$ \over{df}{dx}(x), \over{df}{dx}:\mathbb{R}\rightarrow\mathbb{R}^{N} @f$
!> @param[in,out] x - @b real(kind=8). Value at which @f$f@f$ and @f$\over{df}{dx}@f$ are evaluated.
!> @param[out] hdid - @b real(kind=8). Time step actually accomplished.
!> @param[out] hnext - @b real(kind=8). Estimated next step size.
!> @param[in,out] switch - @b logical. .TRUE. if there was a problem calculating the potential.
!
!------------------------------------------------------------------------------------------------
SUBROUTINE STEP_BS(this,y,dydx,x,hdid,hnext,DERIV,switch)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Bulirsch_Stoer_integrator),INTENT(IN):: this
	REAL(KIND=8),DIMENSION(:),INTENT(IN):: dydx 
	REAL(KIND=8),DIMENSION(:),INTENT(INOUT):: y
	REAL(KIND=8),INTENT(INOUT):: x
	LOGICAL,INTENT(INOUT):: switch
	REAL(KIND=8),INTENT(OUT):: hdid
	REAL(KIND=8),INTENT(OUT):: hnext
   EXTERNAL :: DERIV
   ! Parameters for this routine
	REAL(KIND=8),PARAMETER:: SAFE1 = 0.25D0 ! Safety factor
	REAL(KIND=8),PARAMETER:: SAFE2 = 0.7D0  ! Safety factor
	REAL(KIND=8),PARAMETER:: SCALMX = 0.1D0 ! 1/SCALMX is the maximum factor by which a stepsize can be increased.
	REAL(KIND=8),PARAMETER:: REDMAX = 1.D-5 ! Maximum factor used when a stepsize is reduced
	REAL(KIND=8),PARAMETER:: REDMIN = 0.7D0 ! The minimum
   INTEGER(KIND=4),PARAMETER:: KMAXX = 8 ! Maxmimum row number used in the extrapolation
   INTEGER(KIND=4),PARAMETER:: IMAX = KMAXX+1 ! Next row number
   CHARACTER(LEN=9),PARAMETER:: routinename = "STEP_BS: "
   ! Local variables
   INTEGER:: nv ! total number of dimensions y, dydx
   INTEGER,DIMENSION(IMAX):: nseq
	REAL(KIND=8),DIMENSION(KMAXX):: err
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: ysav,yseq,yerr
	REAL(KIND=8),DIMENSION(IMAX):: a
	REAL(KIND=8),DIMENSION(KMAXX,KMAXX):: alf
   INTEGER(KIND=4):: i,iq,k,kk,km,kmax,kopt
	REAL(KIND=8):: eps1,epsold,errmax,fact,h,red,scale,work,wrkmin,xest, xnew
	LOGICAL:: first,reduct
   SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
   DATA first/.TRUE./,epsold/-1./
   DATA nseq/2,4,6,8,10,12,14,16,18/
   ! HEY, HO! LET'S GO !!! -----------------------------
   SELECT CASE(size(y)/=size(dydx).OR.size(y)/=this%getnv())
      CASE(.TRUE.)
         WRITE(0,*) "STEP_BS: dimension mismatch between y, dydx and variable nv"
         CALL EXIT(1)
      CASE(.FALSE.)
         ALLOCATE(ysav(this%getnv()))
         ALLOCATE(yseq(this%getnv()))
         ALLOCATE(yerr(this%getnv()))
   END SELECT
   switch = .FALSE.
	IF (this%geterr()/=epsold) THEN !A new tolerance, so reinitialize.
		hnext = -1.D29  ! Impossible values.
		xnew  = -1.D29
		eps1 = SAFE1*this%geterr()
		a(1)=nseq(1)+1  ! Compute work coecients Ak.
		DO k=1,KMAXX
			a(k+1)=a(k)+nseq(k+1)
		END DO
		DO iq=2,KMAXX ! Compute alpha(k,q)
			DO k=1,iq-1
				alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.)*(2*k+1)))
			END DO
		END DO
		epsold = this%geterr()
		DO kopt=2,KMAXX-1 ! Determine optimal row number for convergence
			if(a(kopt+1)>a(kopt)*alf(kopt-1,kopt)) goto 1 
		END DO
1 kmax = kopt
	END IF
	h = this%gettimestep_next()
	DO i=1,this%getnv() ! Save the starting values.
		ysav(i) = y(i)
	END DO
   ! A new stepsize or a new integration: re-establish the order window
   SELECT CASE(h/=hnext .OR. x/=xnew)
      CASE(.TRUE.)
         first=.true.
         kopt=kmax
      CASE(.FALSE.)
         ! Do nothing
   END SELECT
	reduct=.false.
2 DO k=1,kmax   ! Evaluate the sequence of modified midpoint integrations
      xnew=x+h
      CALL this%MMID(ysav,dydx,x,h,nseq(k),yseq,DERIV,switch)
      SELECT CASE(switch)
         CASE(.TRUE.)
            RETURN
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      xest=(h/nseq(k))**2.D0 ! Squared, since error series is even.
      SELECT CASE(this%stringparam(1))
         CASE("Rational")
            CALL this%RATIONAL_EXTRAPOL(k,xest,yseq,y,yerr) ! Perform extrapolation with Rational funcions
         CASE("Polinomi")
            CALL this%POLINOM_EXTRAPOL(k,xest,yseq,y,yerr) ! Perform extrapolation with Polinoms
         CASE DEFAULT
            WRITE(0,*) "STEP_BS ERR: Wrong keyword for extrapolation variable"
            CALL EXIT(1)
      END SELECT
      SELECT CASE(k/=1) ! Compute normalized error estimate 
         CASE(.TRUE.)
            errmax=tiny(0.d0)
            DO i=1,this%getnv()
               errmax=max(errmax,dabs(yerr(i)/this%errscal(i)))
            END DO
            errmax=errmax/this%geterr() ! Scale error relative to tolerance.
            km=k-1
            err(km)=(errmax/SAFE1)**(1./(2*km+1))
         CASE(.FALSE.)
            ! Do nothing
      END SELECT
		IF(k.NE.1.AND.(k.GE.kopt-1.OR.first)) THEN ! In order window.
			IF(errmax.lt.1.) GO TO 4  ! Converged.
			IF (k==kmax .OR. k==kopt+1) THEN ! Check for possible stepsize reduction.
				red=SAFE2/err(km)
				GO TO 3
			ELSE IF (k.EQ.kopt) THEN
				IF(alf(kopt-1,kopt)<err(km)) THEN
					red=1./err(km)
					GO TO 3
				END IF
			ELSE IF(kopt.EQ.kmax) THEN
				IF(alf(km,kmax-1)<err(km)) THEN
					red=alf(km,kmax-1)*SAFE2/err(km)
					GO TO 3
				END IF
			ELSE IF (alf(km,kopt)<err(km)) THEN
				red=alf(km,kopt-1)/err(km)
				GO TO 3
			END IF
		END IF
	END DO
3 red=min(red,REDMIN)       ! Reduce stepsize by at least REDMIN and at
	red=max(red,REDMAX) ! most REDMAX.
	h=h*red
	reduct=.TRUE.
	GO TO 2 ! Try again.
4 x=xnew ! Successful step taken.
	hdid=h
	first=.FALSE.
	wrkmin=1.D35 ! Compute optimal row for convergence and
   DO kk=1,km   ! corresponding stepsize.
      fact= max(err(kk),SCALMX)
      work=fact*a(kk+1)
      SELECT CASE(work<wrkmin)
         CASE(.TRUE.)
            scale=fact
            wrkmin=work
            kopt=kk+1
         CASE(.FALSE.)
            ! Do nothing
      END SELECT
   END DO
   hnext=h/scale
   ! Check for possible order increase, but not if stepsize was just reduced
   SELECT CASE(kopt>=k .AND. kopt/=kmax .AND. .NOT.reduct)
      CASE(.TRUE.)
         fact=max(scale/alf(kopt-1,kopt),SCALMX)
         SELECT CASE(a(kopt+1)*fact<=wrkmin)
            CASE(.TRUE.)
               hnext=h/fact
               kopt=kopt+1
            CASE(.FALSE.)
               ! Do nothing
         END SELECT
      CASE(.FALSE.)
         ! Do nothing
   END SELECT
   RETURN
END SUBROUTINE STEP_BS
!############################################################################################
!# SUBROUTINE : PZEXTR_BS ###################################################################
!############################################################################################
!!> @brief
!! Uses polynomial extrapolation to evaluate nv functions at x = 0 by fitting a polynomial to a
!! sequence of estimates with progressively smaller values x = xest, and corresponding function
!! vectors yest(1:nv).
!
!> @param[in] iest - @b integer(kind=4). Number of midpoints
!> @param[in] xest - @b real(kind=8). Value at which @b yest was evaluated.
!> @param[in] yest(:) - @b real(kind=8) @b array. Functions to integrate: @f$ y_{i}(x) @f$
!> @param[out] yz(:) - @b real(kind=8) @b array. Extrapolated function value.
!> @param[out] dy(:) - @b real(kind=8) @b array. Estimated error.
!
!> @details
!! - Addapted from numerical recipes in FORTRAN 77
!
!> @see Numerical recipes in FORTRAN 77 
!--------------------------------------------------------------------------------------------
SUBROUTINE PZEXTR_BS(this,iest,xest,yest,yz,dy)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Bulirsch_Stoer_integrator),INTENT(IN):: this
   INTEGER(KIND=4),INTENT(IN):: iest
   REAL(KIND=8),INTENT(IN):: xest
   REAL(KIND=8),DIMENSION(:),INTENT(IN):: yest
   REAL(KIND=8),DIMENSION(:),INTENT(OUT):: yz
   REAL(KIND=8),DIMENSION(:),INTENT(OUT):: dy
   ! Local Variables
   INTEGER(KIND=4):: j, k1 ! counters
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: d
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE,SAVE:: x
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE,SAVE:: qcol
   REAL(KIND=8):: delta,f1,f2,q
   ! ROCK THE CASBAH !!!! --------------
   SELECT CASE(allocated(qcol).AND.allocated(x)) ! both should be allocated
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         ALLOCATE(qcol(this%getnv(),iest))
         ALLOCATE(x(iest))
   END SELECT
   ALLOCATE(d(this%getnv()))
   x(iest)=xest !Save current independent variable.
   dy(:)=yest(:)
   yz(:)=yest(:)
   SELECT CASE(iest)
      CASE(1)
         qcol(:,1)=yest(:)
      CASE DEFAULT
         d(:)=yest(:)
         DO k1=1,iest-1
            delta=1.d0/(x(iest-k1)-xest)
            f1=xest*delta
            f2=x(iest-k1)*delta
            DO j=1,this%getnv() ! Propagate tableau 1 diagonal more.
               q=qcol(j,k1)
               qcol(j,k1)=dy(j)
               delta=d(j)-q
               dy(j)=f1*delta
               d(j)=f2*delta
               yz(j)=yz(j)+dy(j)
            END DO
         END DO
         qcol(:,iest)=dy(:)
   END SELECT
   RETURN
END SUBROUTINE PZEXTR_BS
!##################################################################################################
!# SUBROUTINE: RZEXTR_BS ##########################################################################
!################################################################################################## 
!> @brief 
!! - Uses diagonal rational function extrapolation to evaluate functions. Same as PZEXTR_BS but uses
!!   rational functions.
!
!>@details
!! Taken from Numerical recipes in Fortran 77
!> @see pzextr
!--------------------------------------------------------------------------------------------------
SUBROUTINE RZEXTR_BS(this,iest,xest,yest,yz,dy)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Bulirsch_Stoer_integrator),INTENT(IN):: this
   INTEGER(KIND=4),INTENT(IN):: iest
   REAL(KIND=8),INTENT(IN):: xest
   REAL(KIND=8),DIMENSION(:),INTENT(IN):: yest
   REAL(KIND=8),DIMENSION(:),INTENT(OUT):: dy,yz
   ! Local variables 	
   INTEGER(KIND=4):: j,k
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: fx
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE,SAVE:: d
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE,SAVE:: x
   REAL(KIND=8):: b,b1,c,ddy,v,yy
   ! HEY, HO!, LETS GO!! ------------------------
   SELECT CASE(allocated(d).AND.allocated(x))
      CASE(.TRUE.)
         ! Do nothing
      CASE(.FALSE.)
         ALLOCATE(d(this%getnv(),iest))
         ALLOCATE(x(iest))
   END SELECT
   ALLOCATE(fx(this%getnv()))
   x(iest)=xest
   !Save current independent variable.
   SELECT CASE(iest)
      CASE(1)
         yz(:)=yest(:)
         d(:,1)=yest(:)
         dy(:)=yest(:)
      CASE DEFAULT
         DO k=1,iest-1
            fx(k+1)=x(iest-k)/xest
         END DO
         DO j=1,this%getnv()
            ! Evaluate next diagonal in tableau.
            yy=yest(j)
            v=d(j,1)
            c=yy
            d(j,1)=yy
            DO k=2,iest
               b1=fx(k)*v
               b=b1-c
               SELECT CASE(b/=0.)
                  CASE(.TRUE.)
                     b=(c-v)/b
                     ddy=c*b
                     c=b1*b
                  CASE(.FALSE.)
                  ! Care needed to avoid division by 0.
                     ddy=v
               END SELECT 
               SELECT CASE(k/=iest)
                  CASE(.TRUE.)
                     v=d(j,k)
                  CASE(.FALSE.)
                     ! do nothing
               END SELECT
               d(j,k)=ddy
               yy=yy+ddy
            END DO
            dy(j)=ddy
            yz(j)=yy
         END DO
      END SELECT 
   RETURN
END SUBROUTINE RZEXTR_BS
END MODULE BULIRSCH_STOER_INTEGRATOR_MOD 
