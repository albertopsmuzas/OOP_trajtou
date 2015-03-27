!################################################################################################### 
! MODULE: CUBICSPLINES_MOD
!
!> @brief  
!! Module that manages cubic splines interpolations 
!
!##################################################################################################
MODULE CUBICSPLINES_MOD
! Initial declarations
use INTERPOL1D_MOD
use MATHS_MOD
#ifdef DEBUG
   use DEBUG_MOD
#endif
IMPLICIT NONE
!//////////////////////////////////////////////////////////////////////
! TYPE: Cubic Splines
!> @brief
!! Extends type Interpol1d. It stores all information and actions that can be extracted 
!! from a cubic splines interpolation
!
!> @param dv2vdz - Set of @f$\frac{d^{2}F(x_{i})}{dx^{2}}@f$, calculated by DSPLIN
!> @param coeff - Matrix of cubic splines coefficients. Coeff(i,j) stands for the
!!                j'th coefficient of the i'th cubic spline
!> @see dsplin
!------------------------------------------------------------------------
TYPE,EXTENDS(Interpol1d) :: Csplines
   PRIVATE
   LOGICAL :: is_initialized = .FALSE.
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: d2fdx 
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: coeff
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE,PUBLIC :: xmin
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE,PUBLIC :: xroot
   CONTAINS
      ! Public procedures
      PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_CUBIC_SPLINES
      PROCEDURE,PUBLIC :: REINTERPOL => REINTERPOL_CSPLINES
      PROCEDURE,PUBLIC :: getvalue => get_csplines_value
      PROCEDURE,PUBLIC :: getderiv => get_csplines_dfdx_value
      PROCEDURE,PUBLIC :: getderiv2 => get_csplines_df2dx2_value
      PROCEDURE,PUBLIC :: GET_V_AND_DERIVS => GET_V_AND_DERIVS_CSPLINES
      PROCEDURE,PUBLIC :: SET_MINIMUM => SET_MINIMUM_CSPLINES
      PROCEDURE,PUBLIC :: SET_XROOT => SET_XROOT_CSPLINES
      ! Private procedures
      PROCEDURE,PRIVATE :: SET_SECOND_DERIVS => DSPLIN
      PROCEDURE,PRIVATE :: SET_COEFF => SET_CUBIC_SPLINES_COEFF
END TYPE Csplines
!//////////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: REINTERPOL_CSPLINES 
!###########################################################
!> @brief
!! Interpolates an already interpolated CSplines job
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!
!> @see interpol_cubic_splines
!-----------------------------------------------------------
SUBROUTINE REINTERPOL_CSPLINES(this,dz1,id1,dz2,id2)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Csplines),INTENT(INOUT) :: this
   REAL(KIND=8),INTENT(IN) :: dz1,dz2
   INTEGER(KIND=4),INTENT(IN) :: id1,id2
   ! Local variables
   ! Run section
   SELECT CASE(this%is_initialized)
      CASE(.FALSE.)
         WRITE(0,*) "REINTERPOL_CSPLINES ERR: use INTERPOL intead of REINTERPOL, this csplines &
                     job was not prior initialized"
         CALL EXIT(1)
      CASE(.TRUE.)
         ! do nothing
   END SELECT
   DEALLOCATE(this%d2fdx)
   DEALLOCATE(this%coeff)
   this%is_initialized=.FALSE.
   CALL this%INTERPOL(dz1,id1,dz2,id2)
   SELECT CASE(ALLOCATED(this%xmin))
      CASE(.TRUE.)
         DEALLOCATE(this%xmin)
         CALL this%SET_MINIMUM()
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE REINTERPOL_CSPLINES
!###########################################################
!# SUBROUTINE: SET_MINIMUM_CSPLINES 
!###########################################################
!> @brief
!! Locates minimums in the interpolation. Analytical solution
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 26/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_MINIMUM_CSPLINES(this)
   ! Initial declarations  
   IMPLICIT NONE
   ! I/O variables
   CLASS(Csplines),TARGET,INTENT(INOUT)::this
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   REAL(KIND=8),DIMENSION(:),POINTER :: coeff
   REAL(KIND=8),POINTER :: x1,x2
   REAL(KIND=8) :: s
   REAL(KIND=8) :: discriminant
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: aux
   CHARACTER(LEN=22),PARAMETER :: routinename="SET_MINIMUM_CSPLINES: "
   ! Run section
   DO i = 1, this%n-1 ! number of splines
      coeff => this%coeff(i,:)
      x1 => this%x(i)
      x2 => this%x(i+1)
      discriminant=coeff(2)**2.D0-3.D0*coeff(1)*coeff(3)
      SELECT CASE( discriminant > 0.D0)
         CASE(.TRUE.)
            s=(-coeff(2)+dsqrt(discriminant))/(3.D0*coeff(1))
         CASE(.FALSE.)
            CYCLE
      END SELECT
      SELECT CASE(s >= 0.D0 .AND. s <= x2-x1)
         CASE(.TRUE.)
            SELECT CASE(ALLOCATED(this%xmin))
               CASE(.TRUE.)
                  ALLOCATE(aux(size(this%xmin)))
                  aux=this%xmin
                  DEALLOCATE(this%xmin)
                  ALLOCATE(this%xmin(size(aux)+1))
                  this%xmin(size(aux)+1)=s+x1
                  DEALLOCATE(aux)
               CASE(.FALSE.)
                  ALLOCATE(this%xmin(1))
                  this%xmin=s+x1
            END SELECT
         CASE(.FALSE.)
            CYCLE
      END SELECT
   END DO
#ifdef DEBUG
   SELECT CASE(ALLOCATED(this%xmin))
      CASE(.TRUE.)
         CALL DEBUG_WRITE(routinename,"Minimums found: ",size(this%xmin))
         CALL DEBUG_WRITE(routinename,this%xmin)
      CASE(.FALSE.)
         CALL DEBUG_WRITE(routinename,"There are not any minimums")
   END SELECT
#endif
   RETURN
END SUBROUTINE SET_MINIMUM_CSPLINES
!###########################################################
!# SUBROUTINE: SET_XROOT_CSPLINES 
!###########################################################
!> @brief
!! Locates roots in the interpolation. Analytical solution
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 26/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_XROOT_CSPLINES(this)
   ! Initial declarations  
   IMPLICIT NONE
   ! I/O variables
   CLASS(Csplines),TARGET,INTENT(INOUT)::this
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   REAL(KIND=8),DIMENSION(:),POINTER :: coeff
   REAL(KIND=8),POINTER :: x1,x2
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: s
   REAL(KIND=8) :: discriminant
   REAL(KIND=8) :: a1,a2,a3,theta
   REAL(KIND=8):: param_q, param_r
   REAL(KIND=8) :: param_s, param_t
   COMPLEX(KIND=8):: param_sc, param_tc
   INTEGER(KIND=4) :: nold
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: aux
   CHARACTER(LEN=22),PARAMETER :: routinename="SET_MINIMUM_CSPLINES: "
   ! Run section
   DO i = 1, this%n-1 ! number of splines
      coeff => this%coeff(i,:)
      x1 => this%x(i)
      x2 => this%x(i+1)
      a1=coeff(2)/coeff(1)
      a2=coeff(3)/coeff(1)
      a3=coeff(4)/coeff(1)
      param_q=(3.D0*a2-a1**2.D0)/9.D0
      param_r=(9.D0*a1*a2-27.D0*a3-2.D0*a1**3.D0)/54.D0
      discriminant=param_q**3.D0+param_r**2.D0
      SELECT CASE(discriminant > 0.D0) ! one real solution
         CASE(.TRUE.)
            param_s=param_r+dsqrt(discriminant)
            SELECT CASE(param_s >= 0.D0)
               CASE(.TRUE.)
                  param_s=(param_s)**(1.D0/3.D0)
               CASE(.FALSE.)
                  param_s=-(dabs(param_s))**(1.D0/3.D0)
            END SELECT
            param_t=param_r-dsqrt(discriminant)
            SELECT CASE(param_t >= 0.D0)
               CASE(.TRUE.)
                  param_t=(param_t)**(1.D0/3.D0)
               CASE(.FALSE.)
                  param_t=-(dabs(param_t))**(1.D0/3.D0)
            END SELECT
            ALLOCATE(s(1))
            s(1)=param_s+param_t-a1/3.D0
         CASE(.FALSE.) ! i.e. discriminant <= 0.D0
            param_sc=dcmplx(param_r,0.D0)+sqrt(dcmplx(discriminant,0.D0))
            param_sc=param_sc**(1.D0/3.D0)
            param_tc=dcmplx(param_r,0.D0)-sqrt(dcmplx(discriminant,0.D0))
            param_tc=param_tc**(1.D0/3.D0)
            ALLOCATE(s(3))
            s(1)=real(param_sc+param_tc-dcmplx(a1/3.D0,0.D0))
            s(2)=real(-0.5D0*(param_sc+param_tc)+dcmplx(0.D0,0.5D0*dsqrt(3.D0))*(param_sc-param_tc)-dcmplx(a1/3.D0,0.D0))
            s(3)=real(-0.5D0*(param_sc+param_tc)-dcmplx(0.D0,0.5D0*dsqrt(3.D0))*(param_sc-param_tc)-dcmplx(a1/3.D0,0.D0))
      END SELECT
      DO j = 1, size(s) ! loop over roots
         SELECT CASE(s(j) >= 0.D0 .AND. s(j) <= x2-x1)
            CASE(.TRUE.)
               SELECT CASE(allocated(this%xroot))
                  CASE(.TRUE.)
                     nold=size(this%xroot)
                     ALLOCATE(aux(nold))
                     aux=this%xroot
                     DEALLOCATE(this%xroot)
                     ALLOCATE(this%xroot(nold+1))
                     this%xroot(1:nold)=aux
                     this%xroot(nold+1)=s(j)+x1
                     DEALLOCATE(aux)
                  CASE(.FALSE.)
                     ALLOCATE(this%xroot(1))
                     this%xroot(1)=s(j)+x1
             END SELECT
          CASE(.FALSE.)
            CYCLE
         END SELECT
      END DO
      SELECT CASE(allocated(s))
         CASE(.TRUE.)
            DEALLOCATE(s)
         CASE(.FALSE.)
            ! do nothing
      END SELECT
   END DO
#ifdef DEBUG
   SELECT CASE(allocated(this%xroot))
      CASE(.TRUE.)
         CALL DEBUG_WRITE(routinename,"Roots found: ",size(this%xmin))
         CALL DEBUG_WRITE(routinename,this%xroot)
      CASE(.FALSE.)
         CALL DEBUG_WRITE(routinename,"There are not any roots")
   END SELECT
#endif
   RETURN
END SUBROUTINE SET_XROOT_CSPLINES

!######################################################################
! SUBROUTINE: DSPLIN ##################################################
!######################################################################
!> @brief
!> For cubic splines interpolation.
!
!> @detailed   
!! Given a set of couples @f$(x_{i},F(x_{i}))@f$ and a set of two conditions (read parameters section), 
!! this module calculates the set of second derivatives at nodes @f$d^{2}S_{i}(x_{i})\over{dx^{2}}@f$, where: 
!! - @f$F(x)=S_{i}(x)@f$ if @f$x\in[x_{i},x_{i+1}]@f$
!! - @f$S_{i}(x)=A_{i}(x-x_{i})^{3}+B_{i}(x-x_{i})^{2}+C_{i}(x-x_{i})+D_{i}@f$.
!!
!! @b Constrains:
!! - @f$S_{i}(x_{i+1})=S_{i+1}(x_{i+1})@f$
!! - @f$\frac{dS_{i}(x_{i+1})}{dx}=\frac{dS_{i+1}(x_{i+1})}{dx}@f$
!! - @f$\frac{d^{2}S_{i}(x_{i+1})}{dx^{2}}=\frac{d^{2}S_{i+1}(x_{i+1})}{dx^{2}}@f$
!! - There should be two extra conditions to get a  system of compatible
!!   linear equations. @b id1 and @b id2 define the kind of conditions that
!!   are supported by this dspline version. @b id1 is the code for @b cond1 (first node condition)
!!   and @b id2 for @b cond2 (last node condition).
!!
!> @param[in,out] cubicspl - Csplines subtype variable. Contains input and output data
!> @param[in] cond1 - if @b id1=1, @b cond1@f$=\frac{d^{2}S(x_{1})}{dx^{2}}@f$; if @b id1=0, @f$S_{1}(x)=S_{2}(x)@f$
!> @param[in] id1 - Integer parameter that controls the meaning of @b cond1
!> @param[in] cond2 - if @b id2=1, @b cond2@f$=\frac{d^{2}S(x_{N})}{dx^{2}}@f$; if @b id2=0, @f$S_{N-1}(x)=S_{N}(x)@f$
!> @param[in] id2 - Integer parameter that controls the meaning of @b cond2
!
!> @author A.P Muzas - alberto.muzas@uam.es
!> @date    20/Jan/2014
!> @version 1.1
!
!> @see debug_mod 
!----------------------------------------------------------------------
SUBROUTINE DSPLIN(cubicspl,cond1,id1,cond2,id2)
        IMPLICIT NONE
        ! I/O Variables
        CLASS(Csplines),TARGET,INTENT(INOUT) :: cubicspl
        REAL(KIND=8),INTENT(IN) :: cond1 
        INTEGER,INTENT(IN) :: id1
        REAL(KIND=8),INTENT(IN) :: cond2
        INTEGER,INTENT(IN) :: id2 
        ! Internal variables
        REAL(KIND=8),DIMENSION(cubicspl%n),TARGET :: diag, indep
        REAL(KIND=8),DIMENSION(cubicspl%n),TARGET :: supdiag, subdiag
        REAL(KIND=8),DIMENSION(cubicspl%n-1) :: h ! space between nodes
        REAL(KIND=8),DIMENSION(cubicspl%n-2) :: sigma
        REAL(KIND=8),DIMENSION(cubicspl%n-2) :: delta
        INTEGER :: i ! Counter
        CHARACTER(LEN=8), PARAMETER :: routinename = "DSPLIN: "
        !Pointers
        INTEGER(KIND=4),POINTER :: n
        REAL(KIND=8),DIMENSION(:),POINTER :: x,y,d2sdx
        REAL(KIND=8),DIMENSION(:),POINTER :: pointdiag, pointindep, pointsupdiag, pointsubdiag
        REAL*8,DIMENSION(:),POINTER :: pointd2sdx
        ! HEY, HO! LET'S GO!
        !=============================
        n => cubicspl%n
        ALLOCATE(cubicspl%d2fdx(1:n))
        x => cubicspl%x(1:n)
        y => cubicspl%f(1:n)
        d2sdx => cubicspl%d2fdx(1:n)
!----------------------------------------------- auxiliar vectors
        ! Get inter-space vector h
        DO i=1,n-1
                h(i)=x(i+1)-x(i)
        END DO
        ! Get sigma and delta auxiliary vectors
        DO i=1, n-2
                sigma(i)=h(i)+h(i+1)
                delta(i)=(1.D0/h(i))+(1.D0/h(i+1))
        END DO
!--------------------------------------------------------
        ! Set diagonal values
        DO i=1, n-2
                diag(i+1)=2.D0*sigma(i)
        END DO
        ! Set supdiag values
        DO i=2,n-1
                supdiag(i)=h(i)
        END DO
        supdiag(n)=0.D0 ! assumption needed by TRIDIA
        ! Set subdiag values
        subdiag(1)=0.D0 ! assumption needed by TRIDIA
        DO i=1,n-2
                subdiag(i+1)=h(i)
        END DO
        ! Set independent terms
        DO i=2, n-1
	        indep(i)=6.D0*((1.D0/h(i-1))*y(i-1)-delta(i-1)*y(i)+(1.D0/h(i))*y(i+1))
        END DO
!------------------- Set diag, supdiag and subdiag parts that depend on id1 and id2
        IF (id1.EQ.1) THEN
                diag(1)=2.D0
                supdiag(1)=1.D0
                indep(1)=(((y(2)-y(1))/h(1))-cond1)*6D0/h(1)
        ELSE IF (id1.EQ.0) THEN
                diag(1)=0.D0 ! Not needed 
                diag(2)=2.D0*(sigma(1)+h(1))
                subdiag(2)=0.D0 ! Dummy for tridia 
                supdiag(1)=0.D0 ! Not needed
                supdiag(2)=h(2)-h(1)
                ! indep(2) has its usual value
        END IF
!
        IF (id2.EQ.1) THEN
                diag(n)=2.D0
                subdiag(n)=1.D0
                indep(n)=(((y(n-1)-y(n))/h(n-1))+cond2)*6.D0/h(n-1)
        ELSE IF (id2.EQ.0) THEN
                diag(n)=0.D0 ! Not needed
                diag(n-1)=2.D0*(sigma(n-2)+h(n-1))
                supdiag(n-1)=0.D0 ! Dummy for tridia
                subdiag(n)=0.D0 ! Not needed
                subdiag(n-1)=h(n-2)-h(n-1)
                ! indep(n-1) has its usual value
        END IF
!------------------ Set the correct number of dimensions depending on id1 and id2
        IF ((id1.EQ.0).AND.(id2.EQ.0)) THEN
                pointdiag => diag(2:n-1)
                pointsupdiag => supdiag(2:n-1)
                pointsubdiag => subdiag(2:n-1)
                pointindep => indep(2:n-1)
                pointd2sdx => d2sdx(2:n-1)
                CALL TRIDIA(n-2,pointsubdiag,pointdiag,pointsupdiag,pointindep,pointd2sdx)
                d2sdx(1)=2*d2sdx(2)-d2sdx(3)
                d2sdx(n)=2*d2sdx(n-1)-d2sdx(n-2)
        ELSE IF (id1.EQ.0) THEN
                pointdiag => diag(2:n) 
                pointsupdiag => supdiag(2:n)
                pointsubdiag => subdiag(2:n)
                pointindep => indep(2:n)
                pointd2sdx => d2sdx(2:n)
                CALL TRIDIA(n-1,pointsubdiag,pointdiag,pointsupdiag,pointindep,pointd2sdx)
                d2sdx(1)=2*d2sdx(1)-d2sdx(2)
        ELSE IF (id2.EQ.0) THEN
                pointdiag => diag(1:n-1)
                pointsupdiag => supdiag(1:n-1)
                pointsubdiag => subdiag(1:n-1)
                pointindep => indep(1:n-1)
                pointd2sdx => d2sdx(1:n-1)
                CALL TRIDIA(n-1,pointsubdiag,pointdiag,pointsupdiag,pointindep,pointd2sdx)
                d2sdx(n)=2*d2sdx(n-1)-d2sdx(n-2)
        ELSE
                CALL TRIDIA(n,subdiag,diag,supdiag,indep,d2sdx)
        END IF
END SUBROUTINE DSPLIN
!######################################################################
! SUBROUTINE: SET_CUBIC_SPLINES_COEFF #################################
!######################################################################
!> @brief 
!! Calculates cubic splines coefficients a,b,c,d for a given
!! set of Z(i), F(Z(i)), F''(Z(i))
!
!> @warning 
!! - DSPLIN routine should have been executed before, i.e. second derivatives
!!   should be available
!
!> @author A.P Muzas - alberto.muzas@uam.es
!> @date 20/Jan/2014
!> @version 1.0
!
!> @see dsplin 
!----------------------------------------------------------------------
SUBROUTINE SET_CUBIC_SPLINES_COEFF(this)
        IMPLICIT NONE
        ! I/O variable ==============================================
        CLASS(Csplines),TARGET,INTENT(INOUT) :: this
        ! Internal Variables
        INTEGER :: i ! counter
        REAL(KIND=8) :: h ! step between nodes
        REAL(KIND=8) :: x1,x2,y1,y2 ! values of nodes and the function defined there
        CHARACTER(LEN=13), PARAMETER :: routinename = "SET_SPLINES: "
        ! Pointers
        INTEGER(KIND=4),POINTER :: n ! Number of nodes
        REAL(KIND=8),DIMENSION(:),POINTER :: z
        REAL(KIND=8),DIMENSION(:),POINTER :: v
        REAL(KIND=8), POINTER :: a,b,c,d ! pointers for coefficients
        REAL(KIND=8), POINTER :: m1,m2 ! pointers for second derivatives
        ! Run section ========================================
        n => this%n
        z => this%x(1:n)
        v => this%f(1:n)
        !
        SELECT CASE(this%is_initialized)
           CASE(.TRUE.)
               WRITE(0,*) "SET_SPLINES ERR: this splines were initialized before"
               CALL EXIT(1)
           CASE(.FALSE.)
              ! do nothing
        END SELECT
        ALLOCATE (this%coeff(n -1,4)) ! There are N-1 Splines and 4 coefficients for each one of them
        DO i=1, n-1
                a => this%coeff(i,1)
                b => this%coeff(i,2)
                c => this%coeff(i,3)
                d => this%coeff(i,4)
                m1 => this%d2fdx(i)
                m2 => this%d2fdx(i+1)
                y1 = v(i)
                y2 = v(i+1)
                x1 = z(i)
                x2 = z(i+1)
!
                h=x2-x1
                a=(m2-m1)/(6.D0*h)
                b=m1/2.D0
                c=((y2-y1)/h)-((m2+2.D0*m1)/6.D0)*h
                d=y1
        END DO
        this%is_initialized=.TRUE.
        RETURN
END SUBROUTINE SET_CUBIC_SPLINES_COEFF
!######################################################################
! SUBROUTINE: INTERPOL_CUBIC_SPLINES ##################################
!######################################################################
!> @brief
!! This subroutine coordinates a complete cubic spline interpolation.
!
!> @details 
!! - Contidions @b id1 and @b id2 controls the meaning of @b dz1 and @b dz2. See dsplin
!!   routine
!! - This subroutine allocates csplines%d2vdz(:) with the same size as z(:) or v(:).
!! - Obviously, z(:) and v(:) must have the same size.
!
!> @param[in,out] this - Csplines subtype variable that contains all the cubic splines interpolation data
!> @param[in] dz1 - maybe, second derivative at first node
!> @param[in] id1 - Integer control variable that controls the meaning of @ dz1
!> @param[in] dz2 - maybe, second derivative at last node
!> @param[in] id2 - Integer control variable that controls the meaning of @ dz2
!> @param[in] filename - Output file in which cubic splines coefficients will be printed (optional). If it is not
!!                       given or set to "None", no output file will be printed.
!
!> @author A.P. Muzas - alberto.muzas@uam.es
!> @date 03/Feb/2014
!> @version 1.1
!
!> @see interpol_cubic_splines
!----------------------------------------------------------------------
SUBROUTINE INTERPOL_CUBIC_SPLINES(this,dz1,id1,dz2,id2)
        ! Initial declarations
        IMPLICIT NONE
        ! I/O variables ==================================
        CLASS(Csplines),TARGET,INTENT(INOUT) :: this 
        REAL(KIND=8),INTENT(IN) :: dz1,dz2 
        INTEGER(KIND=4),INTENT(IN) :: id1,id2 
        ! Local variables
        CHARACTER(LEN=24),PARAMETER :: routinename="INTERPOL_CUBIC_SPLINES: "
        ! Pointers
        INTEGER(KIND=4),POINTER :: n
        REAL(KIND=8),DIMENSION(:),POINTER :: z, v
        ! Run section =====================================
        ! Allocate pointers
        n => this%n
        z => this%x(1:n)
        v => this%f(1:n)
        CALL this%SET_SECOND_DERIVS(dz1,id1,dz2,id2)
        CALL this%SET_COEFF()
        RETURN
END SUBROUTINE INTERPOL_CUBIC_SPLINES
!###########################################################
!# FUNCTION: get_csplines_value ############################ 
!###########################################################
!> @brief 
!! Gets @f$F(r)@f$ using cubic interpolation
!
!> @param[in] this - Csplines subtype variable that contains all the information needed
!> @param[in] x - Point to evaluate @f$F(x)@f$
!> @param[in] shift - Optional parameter that if it exists and shift=@f$s@f$, we get @f$F(x+s)@f$ instead of @f$F(x)@f$
!
!> @warning
!! - This function doesn't work if we want to evaluate our fuction outside the range
!!   in which the interpolation was defined. Though, there is a confidence value of 1e-6 units around the
!!   extremes of this range.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 28/Jan/2014
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION get_csplines_value(this,x,shift) 
   ! I/O Variables
   CLASS(Csplines),INTENT(IN),TARGET :: this
   REAL(KIND=8),INTENT(IN) :: x
   REAL(KIND=8),OPTIONAL,INTENT(IN) :: shift
   ! Local variables
   INTEGER,POINTER :: n
   REAL(KIND=8),DIMENSION(:),POINTER :: z,v
   REAL(KIND=8):: z1, z2, a, b, c, d
   INTEGER :: i ! Counter
   REAL(KIND=8) :: r
   ! MAY THE FORCE BE WITH YOU --------------------------------------
   n => this%n
   z => this%x(1:n)
   v => this%f(1:n)
   !
   SELECT CASE(present(shift))
      CASE(.TRUE.)
         r=x+shift
      CASE(.FALSE.)
         r=x
   END SELECT
   SELECT CASE(r<this%x(1) .OR. r>this%x(n))
      CASE(.TRUE.)
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: r outside interpolation interval"
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: r = ", r
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: range from ", z(1), "to", z(n)
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   !
   DO i=1,n-1
      z1 = z(i)
      z2 = z(i+1)
      a = this%coeff(i,1)
      b = this%coeff(i,2)
      c = this%coeff(i,3)
      d = this%coeff(i,4)
      SELECT CASE((r<=z2).AND.(r>=z1))
         CASE(.TRUE.)
            EXIT
         CASE(.FALSE.)
            ! do nothing
      END SELECT
   END DO
   ! Now, counter i has the correct label
   get_csplines_value=a*((r-z1)**3.D0)+b*((r-z1)**2.D0)+c*(r-z1)+d
   RETURN
END FUNCTION get_csplines_value
!#################################################################
! FUNCTION: get_csplines_dfdx_value
!> @brief
!! Gets @f$\frac{dF(x)}{dx}@f$ using cubic interpolation
!
!> @param[in] this - Csplines subtype variable that contains all the information needed
!> @param[in] x - Point to evaluate @f$\frac{dF(x)}{dx}@f$
!> @param[in] shift - Optional parameter that if it exists and shift=@f$s@f$,
!!                    we get @f$\frac{dF(x+s)}{dx}@f$ instead of @f$\frac{dF(x)}{dx}@f$
!
!> @warning
!! - This function doesn't work if we want to evaluate our fuction outside the range
!!   in which the interpolation was defined. Though, there is a confidence value of 1e-6 units around the
!!   extremes of this range.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 30/Jan/2014
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION get_csplines_dfdx_value(this,x,shift)
   ! Initial declaration
   IMPLICIT NONE
   ! I/O variables
   CLASS(Csplines),INTENT(IN),TARGET :: this
   REAL(KIND=8),INTENT(IN) :: x
   REAL(KIND=8),OPTIONAL,INTENT(IN) :: shift
   ! Local variables
   INTEGER :: i ! Counter
   REAL(KIND=8) :: r
   ! Pointers 
   INTEGER,POINTER :: n
   REAL(KIND=8):: z1,z2,a,b,c,d
   REAL(KIND=8),DIMENSION(:),POINTER :: z,v
   ! Run section -----------------------------------
   n => this%n
   z => this%x(1:n)
   v => this%f(1:n)
   SELECT CASE(present(shift))
      CASE(.TRUE.)
         r=x+shift
      CASE(.FALSE.)
         r=x
   END SELECT
   SELECT CASE(r<this%x(1) .OR. r>this%x(n))
      CASE(.TRUE.)
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: r outside interpolation interval"
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: r = ", r
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: range from ", z(1), "to", z(n)
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   !
   DO i=1,n-1
      z1 = this%x(i)
      z2 = this%x(i+1)
      a = this%coeff(i,1)
      b = this%coeff(i,2)
      c = this%coeff(i,3)
      SELECT CASE((r<=z2).AND.(r>=z1))
         CASE(.TRUE.)
            EXIT
         CASE(.FALSE.)
            ! do nothing
      END SELECT
   END DO
   ! Now, counter i has the correct label
   get_csplines_dfdx_value=3.D0*(a*((r-z1)**2.D0))+2.D0*b*(r-z1)+c
   RETURN
END FUNCTION get_csplines_dfdx_value

REAL(KIND=8) FUNCTION get_csplines_df2dx2_value(this,x,shift)
   ! Initial declaration
   IMPLICIT NONE
   ! I/O variables
   CLASS(Csplines),INTENT(IN),TARGET :: this
   REAL(KIND=8),INTENT(IN) :: x
   REAL(KIND=8),OPTIONAL,INTENT(IN) :: shift
   ! Local variables
   INTEGER :: i ! Counter
   REAL(KIND=8) :: r
   ! Pointers
   INTEGER,POINTER :: n
   REAL(KIND=8):: z1,z2,a,b,c,d
   REAL(KIND=8),DIMENSION(:),POINTER :: z,v
   ! Run section -----------------------------------
   n => this%n
   z => this%x(1:n)
   v => this%f(1:n)
   SELECT CASE(present(shift))
      CASE(.TRUE.)
         r=x+shift
      CASE(.FALSE.)
         r=x
   END SELECT
   SELECT CASE(r<this%x(1) .OR. r>this%x(n))
      CASE(.TRUE.)
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: r outside interpolation interval"
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: r = ", r
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: range from ", z(1), "to", z(n)
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   !
   DO i=1,n-1
      z1 = this%x(i)
      z2 = this%x(i+1)
      a = this%coeff(i,1)
      b = this%coeff(i,2)
      SELECT CASE((r<=z2).AND.(r>=z1))
         CASE(.TRUE.)
            EXIT
         CASE(.FALSE.)
            ! do nothing
      END SELECT
   END DO
   ! Now, counter i has the correct label
   get_csplines_df2dx2_value=6.D0*(a*(r-z1))+2.D0*b
   RETURN
END FUNCTION get_csplines_df2dx2_value
!###############################################################
! SUBROUTINE: GET_V_AND_DERIVS_CSPLINE
!> @brief
!! Computes F(X) and F'(X) at the same time. Better time performance
!> @details
!! Some guidelines: let's call @f$f(x)=f(r)@f$ where @f$r=x+\delta@f$
!! f(x) is the function shifted and f(r) the old one. Take this into
!! account to use correctly the shift option.
!
!> @param[in]  x  - point to evaluate the function (shifted function) 
!> @param[out] pot   - f(x)
!> @param[out] deriv - @f$\over{df(x)}{dx}@f$
!> @param[in,optional] shift - Turns on the shifted option
!---------------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_CSPLINES(this,x,pot,deriv,shift)
   IMPLICIT NONE
   CLASS(Csplines),TARGET,INTENT(IN) :: this
   REAL(KIND=8),INTENT(IN) :: x
   REAL(KIND=8),INTENT(OUT) :: pot,deriv
   REAL(KIND=8),OPTIONAL,INTENT(IN) :: shift
   ! Local variables
   INTEGER :: i ! Counter
   REAL(KIND=8) ::  r
   ! Pointers 
   INTEGER,POINTER :: n
   REAL(KIND=8) :: z1,z2,a,b,c,d
   REAL(KIND=8),DIMENSION(:),POINTER :: z,v
   ! Run section -----------------------------------
   n => this%n
   z => this%x(1:n)
   v => this%f(1:n)
   SELECT CASE(present(shift))
      CASE(.TRUE.)
         r=x+shift
      CASE(.FALSE.) 
         r=x 
   END SELECT
   SELECT CASE(r<this%x(1).OR. r>this%x(n))
      CASE(.TRUE.)
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: r outside interpolation interval"
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: r = ", r
         WRITE(0,*) "GET_V_AND_DERIVS_CSPLINES ERR: range from ", z(1), "to", z(n)
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   !
   DO i=1,n-1
      z1=this%x(i)
      z2=this%x(i+1)
      a=this%coeff(i,1)
      b=this%coeff(i,2)
      c=this%coeff(i,3)
      d=this%coeff(i,4)
      SELECT CASE((r<=z2).AND.(r>=z1))
         CASE(.TRUE.)
            EXIT
         CASE(.FALSE.)
            ! do nothing
      END SELECT
   END DO
   ! Now, counter i has the correct label
   pot = a*((r-z1)**3.D0)+b*((r-z1)**2.D0)+c*(r-z1)+d
   deriv = 3.D0*(a*((r-z1)**2.D0))+2.D0*b*(r-z1)+c
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_CSPLINES
END MODULE CUBICSPLINES_MOD
