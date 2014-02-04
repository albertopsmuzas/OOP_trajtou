!################################################################################################### 
! MODULE: INTERPOL1D
!
!> @brief  
!! Module that manages different interpolation schemes for one variable
!
!> @details 
!! All types and subprograms intended to create interpolations in 1D should be placed inside this module
!##################################################################################################
MODULE INTERPOL1D_MOD
! Initial declarations
IMPLICIT NONE
!//////////////////////////////////////////////////////////////////////
! TYPE: interpol1d
!> @brief
!! Generic type of one dimensional interpolations
!
!> @details
!! All kinds of interpolations should be declared as extensions of this type
!
!> @param n - Number of data
!> @param x - A set of @f$x_{i}@f$ data
!> @param f - A set of @f$f(x_{i})@f$ data
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 28/Jan/2014
!> @version 1.0
!-----------------------------------------------------------------------
TYPE :: Interpol1d
   PRIVATE
   INTEGER(KIND=4) :: n
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: x
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: f
END TYPE Interpol1d
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
   CONTAINS
      ! Public procedures
      PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_CUBIC_SPLINES
      PROCEDURE,PUBLIC :: getvalue => get_csplines_value
      PROCEDURE,PUBLIC :: getderiv => get_csplines_dfdx_value
      ! Private procedures
      PROCEDURE,PRIVATE :: SET_SECOND_DERIVS => DSPLIN
      PROCEDURE,PRIVATE :: SET_COEFF => SET_CUBIC_SPLINES_COEFF
END TYPE Csplines
!//////////////////////////////////////////////////////////////////////
CONTAINS
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
!> @warning
!! Uses LAPACKCONTROL_MOD and maybe DEBUG_MOD
!
!> @author A.P Muzas - alberto.muzas@uam.es
!> @date    20/Jan/2014
!> @version 1.1
!
!> @see debug_mod, lapackcontrol_mod
!----------------------------------------------------------------------
SUBROUTINE DSPLIN(cubicspl,cond1,id1,cond2,id2)
#ifdef DEBUG
        USE DEBUG_MOD
#endif
        USE LAPACKCONTROL_MOD
        IMPLICIT NONE
        ! I/O Variables
        CLASS(Csplines),TARGET,INTENT(INOUT) :: cubicspl
        REAL(KIND=8),INTENT(IN) :: cond1 
        INTEGER,INTENT(IN) :: id1
        REAL(KIND=8),INTENT(IN) :: cond2
        INTEGER,INTENT(IN) :: id2 
        ! Internal variables
        INTEGER(KIND=4) :: info
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
        ! HEY, HO! LET'S GO!
        !=============================
        n => cubicspl%n
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
                pointdiag => diag(2:n-1) ! should've n-2 dimensions
                pointsupdiag => supdiag(2:n-2) ! should've n-3 dimensions
                pointsubdiag => subdiag(3:n-1) ! should've n-3 dimensions
                pointindep => indep(2:n-1) ! should've n-2 dimensions
                CALL DGTSV(n-2,1,pointsubdiag,pointdiag,pointsupdiag,pointindep,n-2,info)
                CALL LAPACK_CHECK("DGTSV",info)
                !CALL TRIDIA(n-2,pointsubdiag,pointdiag,pointsupdiag,pointindep,pointd2sdx)
                d2sdx(1)=2*d2sdx(2)-d2sdx(3)
                d2sdx(n)=2*d2sdx(n-1)-d2sdx(n-2)
                FORALL(i=2:n-1) d2sdx(i)=indep(i)
        ELSE IF (id1.EQ.0) THEN
                pointdiag => diag(2:n) ! Should've n-1 dimensions
                pointsupdiag => supdiag(2:n-1) ! Should've n-2 dimensions
                pointsubdiag => subdiag(3:n) ! should've n-2 dimensions
                pointindep => indep(2:n) ! should've n-1 dimensions
                CALL DGTSV(n-1,1,pointsubdiag,pointdiag,pointsupdiag,pointindep,n-1,info)
                CALL LAPACK_CHECK("DGTSV",info)
                !CALL TRIDIA(n-1,pointsubdiag,pointdiag,pointsupdiag,pointindep,pointd2sdx)
                d2sdx(1)=2*d2sdx(1)-d2sdx(2)
                FORALL(i=2:n) d2sdx(i)=indep(i)
        ELSE IF (id2.EQ.0) THEN
                pointdiag => diag(1:n-1) ! should've n-1 dimensions
                pointsupdiag => supdiag(1:n-2) ! should've n-2 dimensions
                pointsubdiag => subdiag(2:n-1) ! should've n-2 dimensions
                pointindep => indep(1:n-1) ! should've n-1 dimensions
                CALL DGTSV(n-1,1,pointsubdiag,pointdiag,pointsupdiag,pointindep,n-1,info)
                CALL LAPACK_CHECK("DGTSV",info)
                !CALL TRIDIA(n-1,pointsubdiag,pointdiag,pointsupdiag,pointindep,pointd2sdx)
                d2sdx(n)=2*d2sdx(n-1)-d2sdx(n-2)
                FORALL(i=1:n-1) d2sdx(i)=indep(i)
        ELSE
                pointsubdiag => subdiag(2:n)
                pointsupdiag => supdiag(1:n-1)
                CALL DGTSV(n,1,pointsubdiag,diag,pointsupdiag,indep,n,info)
                CALL LAPACK_CHECK("DGTSV",info)
                d2sdx=indep
                !CALL TRIDIA(n,subdiag,diag,supdiag,indep,d2sdx)
        END IF
#ifdef DEBUG
        CALL VERBOSE_WRITE(routinename,"Second derivatives at nodes calculated -----> Done")
#endif
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
SUBROUTINE SET_CUBIC_SPLINES_COEFF(this,filename)
#ifdef DEBUG
        USE DEBUG_MOD
#endif
        IMPLICIT NONE
        ! I/O variable ==============================================
        CLASS(Csplines),TARGET,INTENT(INOUT) :: this
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: filename
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
        IF (this%is_initialized) THEN
           WRITE(0,*) "SET_SPLINES ERR: this splines were initialized before"
           CALL EXIT(1)
        END IF
        IF (size(z)/=size(v)) THEN
           WRITE(0,*) "SET SPLINES ERR: v and z arrays don't have the same size"
           WRITE(0,*) "SET_SPLINES ERR: Z size is:",size(z)
           WRITE(0,*) "SET SPLINES ERR: V size is:",size(v)
           CALL EXIT(1)
        END IF
        ALLOCATE (this%coeff(n -1,4)) ! There are N-1 Splines and 4 coefficients for each one of them
        IF ((filename.NE."None").AND.(present(filename))) OPEN(11, file=filename, status="replace")
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
!
                IF ((filename.NE."None").AND.(present(filename))) WRITE(11,*) x1, x2, a, b, c, d
        END DO
#ifdef DEBUG
        CALL VERBOSE_WRITE(routinename,"Set coefficients  -----> Done")
        IF ((filename.NE."None").AND.(present(filename))) CALL VERBOSE_WRITE(routinename,"Data stored inside ", filename)
#endif
        IF ((filename.NE."None").AND.(present(filename))) CLOSE(11)
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
SUBROUTINE INTERPOL_CUBIC_SPLINES(this,dz1,id1,dz2,id2,filename)
        ! Initial declarations
#ifdef DEBUG
        USE DEBUG_MOD
#endif
        IMPLICIT NONE
        ! I/O variables ==================================
        CLASS(Csplines),TARGET,INTENT(INOUT) :: this 
        REAL(KIND=8),INTENT(IN) :: dz1,dz2 
        INTEGER(KIND=4),INTENT(IN) :: id1,id2 
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: filename !> file name to store interpolation coefficients
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
        IF (size(z)/=size(v)) THEN
           WRITE(0,*) "INTERPOL_CUBIC_SPLINES ERR: v and z arrays don't have the same size"
           WRITE(0,*) "INTERPOL_CUBIC_SPLINES ERR: Z size is:",size(z)
           WRITE(0,*) "INTERPOL_CUBIC_SPLINES ERR: V size is:",size(v)
           CALL EXIT(1)
        END IF
        !
        CALL this%SET_SECOND_DERIVS(dz1,id1,dz2,id2)
        !
        IF(present(filename)) THEN
           CALL this%SET_COEFF(filename)
        ELSE
           CALL this%SET_COEFF()
        END IF
#ifdef DEBUG
        IF (present(filename)) THEN
           CALL VERBOSE_WRITE(routinename,"1D interpolation stored in",filename)
        END IF
#endif
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
   REAL(KIND=8), POINTER :: z1, z2, a, b, c, d
   INTEGER :: i ! Counter
   REAL(KIND=8) :: margen, r
   ! MAY THE FORCE BE WITH YOU --------------------------------------
   margen = 1.0D-6
   n => this%n
   z => this%x(1:n)
   v => this%f(1:n)
   !
   r=x
   IF(present(shift)) r=r+shift
   IF ((r.lt.(z(1)-margen)).or.(r.gt.(z(n)+margen))) THEN
      WRITE(0,*) "get_csplines_value ERR: r outside the interpolation interval"
      WRITE(0,*) "get_csplines_value ERR: r = ", r
      WRITE(0,*) "get_csplines_value ERR: range from ", z(1)-margen, "to", z(n)+margen
      STOP
   ELSE IF ((r.LE.(z(n)+margen)).AND.(r.GT.z(n))) THEN
      get_csplines_value=z(n)
      RETURN
   ELSE IF ((r.GE.(z(n)+margen)).AND.(r.LT.z(n))) THEN
      get_csplines_value=z(1)
      RETURN
   END IF
   DO i=1,n-1
      z1 => z(i)
      z2 => z(i+1)
      a => this%coeff(i,1)
      b => this%coeff(i,2)
      c => this%coeff(i,3)
      d => this%coeff(i,4)
     IF ((r.le.z2).and.(r.ge.z1)) EXIT
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
        REAL(KIND=8) :: margen, r
        ! Pointers 
        INTEGER,POINTER :: n
        REAL(KIND=8),POINTER :: z1,z2,a,b,c,d
        REAL(KIND=8),DIMENSION(:),POINTER :: z,v
        ! Run section -----------------------------------
        n => this%n
        z => this%x(1:n)
        v => this%f(1:n)
        margen=1D-6
        r=x
        IF(present(shift)) r=r+shift
        !
        IF ((r.lt.this%x(1)).or.(r.gt.this%x(n))) THEN
                WRITE(0,*) "get_csplines_dfdx_value ERR: r outside interpolation interval"
                WRITE(0,*) "get_csplines_dfdx_value ERR: r = ", r
                WRITE(0,*) "get_csplines_dfdx_value ERR: range from ", z(1)-margen, "to", z(n)+margen
                STOP
        ELSE IF ((r.LE.(z(n)+margen)).AND.(r.GT.z(n))) THEN
           get_csplines_dfdx_value=z(n)
           RETURN
        ELSE IF ((r.GE.(z(n)+margen)).AND.(r.LT.z(n))) THEN
           get_csplines_dfdx_value=z(1)
           RETURN
        END IF
        !
        DO i=1,n-1
                z1 => this%x(i)
                z2 => this%x(i+1)
                a => this%coeff(i,1)
                b => this%coeff(i,2)
                c => this%coeff(i,3)
                d => this%coeff(i,4)
                IF ((r.le.z2).and.(r.ge.z1)) EXIT
        END DO
        ! Now, counter i has the correct label
        get_csplines_dfdx_value=3.D0*(a*((r-z1)**2.D0))+2.D0*b*(r-z1)+c
        RETURN
END FUNCTION get_csplines_dfdx_value
END MODULE INTERPOL1D_MOD
