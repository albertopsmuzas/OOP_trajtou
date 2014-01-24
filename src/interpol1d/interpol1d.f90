!######################################################### 
! MODULE: INTERPOL1D
!
!> \brief  Module that manages different interpolation schemes for one variable
!
!> \details All types and subprograms intended to create interpolations in 1D should be placed
!!          inside this module
!
!> \author A.P. Muzas - alberto.muzas@uam.es
!
!> \date  19/Dec/2013
!
!> \version 1.0
!
!##########################################################
MODULE INTERPOL1D_MOD
! Initial declarations
IMPLICIT NONE
!//////////////////////////////////////////////////////////////////////
! TYPE: Zinterpol
TYPE :: interpol1d
END TYPE interpol1d
!//////////////////////////////////////////////////////////////////////
! TYPE: Cubic Splines
TYPE,EXTENDS(interpol1d) :: Csplines
   PRIVATE
   LOGICAL :: is_initialized = .FALSE.
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: d2vdz !> Collection of second derivatives at each point
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: coeff !> Matrix of cubic splines coefficients: 
   !!                                                  - coeff(i,j) stands for the j'th coefficient of the i'th
   !!                                                    cubic spline.
   CONTAINS
      PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_CUBIC_SPLINES
END TYPE Csplines
!//////////////////////////////////////////////////////////////////////
CONTAINS
!######################################################################
! SUBROUTINE: DSPLIN ##################################################
!######################################################################
!> \brief   For cubic splines interpolation.
!
!> \detailed   Calculates S''i(xi) [ second derivatives at nodes ]
!!             Si(Xi)=ai*(xi**3)+bi*(xi**2)+ci*xi+di
!!             Splines satisfy continuity , smothness and second derivatives cointinuity.
!!             There should be two extra conditions to get a  system of compatible
!!             linear equations. "id1" and "id2" define the kind of conditions that
!!             are supported by this dspline version. "id1" is the code for "cond1" (first node condition)
!!             and  "id2" for "cond2" (last node condition).
!!             There are n-1 cubic splines
!
!> \author A.P Muzas - alberto.muzas@uam.es
!
!> \date    20/Jan/2014
!
!> \version 1.1
!----------------------------------------------------------------------
SUBROUTINE DSPLIN(n,x,y,cond1,id1,cond2,id2, d2sdx)
        USE DEBUG_MOD
        USE LAPACKCONTROL_MOD
        IMPLICIT NONE
        ! I/O Variables
        INTEGER,INTENT(IN) :: n !> Number of nodes (number of data)
        REAL(KIND=8),DIMENSION(n),INTENT(IN) :: x !> list of x values
        REAL(KIND=8),DIMENSION(n),INTENT(IN) :: y !> list of y values
        REAL(KIND=8),INTENT(IN) :: cond1 !> if id1=1; this is the value of the first second derivative S''(x(1))=cond1
        !!                                  if id1=0; first two intervals have the same cubic spline. A value is not needed.
        INTEGER,INTENT(IN) :: id1 !> Code that identifies the kind of extra conditions that we are using in the
        !!                            cubic splines procedure. Currently, it supports id1=0,1
        REAL(KIND=8),INTENT(IN) :: cond2 !> if id2=1; this is the value of the last second derivative S''(x(n))=cond2
        !!                                  if id2=0; last two intervals have the same cubic spline. A value is not needed
        INTEGER,INTENT(IN) :: id2 !> Code that identifies the kind of extra conditions that we are using in the 
        !!                           cubic splines procedure. Currently, it supports id2=0,1
        REAL(KIND=8),DIMENSION(n),INTENT(OUT),TARGET :: d2sdx !> Stores second derivatives at nodes.
        ! Internal variables
        INTEGER(KIND=4) :: info
        REAL(KIND=8),DIMENSION(n), TARGET :: diag, indep
        REAL(KIND=8),DIMENSION(n), TARGET :: supdiag, subdiag
        REAL(KIND=8),DIMENSION(n-1) :: h ! space between nodes
        REAL(KIND=8),DIMENSION(n-2) :: sigma
        REAL(KIND=8),DIMENSION(n-2) :: delta
        REAL(KIND=8),POINTER :: pointdiag(:), pointindep(:) , pointsupdiag(:), pointsubdiag(:)
        INTEGER :: i ! Counter
        CHARACTER(LEN=8), PARAMETER :: routinename = "DSPLIN: "
        ! HEY, HO! LET'S GO!
        !=============================
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
! SUBROUTINE: SET_CUBIC_SPLINES #######################################
!######################################################################
!> \brief Calculates cubic splines coefficients a,b,c,d for a given
!!        set of Z(i), F(Z(i)), F''(Z(i))
!
!> \details - csplines%coeff(:) shouldn't be allocated before running this subprogram
!>          - csplines%d2vdz(:) should be allocated before and contain relevant data, i.e.
!!            second derivatives at nodes should've been calculated in a prior calculation
!!            , for example, using DSPLIN
!
!> \author A.P Muzas - alberto.muzas@uam.es
!
!> \date 20/Jan/2014
!
!> \version 1.0
!----------------------------------------------------------------------
SUBROUTINE SET_CUBIC_SPLINES(this,z,v,filename)
        USE DEBUG_MOD
        IMPLICIT NONE
        ! I/O variable ==============================================
        CLASS(Csplines),TARGET,INTENT(INOUT) :: this
        REAL(KIND=8),DIMENSION(:),INTENT(IN) :: z
        REAL(KIND=8),DIMENSION(:),INTENT(IN) :: v
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: filename
        ! Internal Variables
        INTEGER :: i ! counter
        REAL(KIND=8) :: h ! step between nodes
        REAL(KIND=8) :: x1,x2,y1,y2 ! values of nodes and the function defined there
        CHARACTER(LEN=13), PARAMETER :: routinename = "SET_SPLINES: "
        ! Pointers
        INTEGER(KIND=4) :: n ! Number of nodes
        REAL(KIND=8), POINTER :: a,b,c,d ! pointers for coefficients
        REAL(KIND=8), POINTER :: m1,m2 ! pointers for second derivatives
        ! Run section ========================================
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
        n=size(z)
        ALLOCATE (this%coeff(n -1,4)) ! There are N-1 Splines
        ALLOCATE (this%d2vdz(1:n)) ! There are N second derivatives
        IF ((filename.NE."None").AND.(present(filename))) OPEN(11, file=filename, status="replace")
        DO i=1, n-1
                a => this%coeff(i,1)
                b => this%coeff(i,2)
                c => this%coeff(i,3)
                d => this%coeff(i,4)
                m1 => this%d2vdz(i)
                m2 => this%d2vdz(i+1)
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
END SUBROUTINE SET_CUBIC_SPLINES
!######################################################################
! SUBROUTINE: INTERPOL_CUBIC_SPLINES ##################################
!######################################################################
!> \brief This subroutine creates a cubic spline interpolation.
!
!> \details - The boundaries of the cubic splines are the standard ones, i.e. we have
!!            to set the first derivatives at first and last gridpoint.
!!          - This subroutine allocates csplines%d2vdz(:) with the same size as z(:) or v(:).
!!          - Obviously, z(:) and v(:) must have the same size.
!
!> \author A.P. Muzas - alberto.muzas@uam.es
!
!> \date 20/Jan/2014
!
!> \version 1.0
!----------------------------------------------------------------------
SUBROUTINE INTERPOL_CUBIC_SPLINES(this,v,z,dz1,dz2,filename)
        ! Initial declarations
        IMPLICIT NONE
        ! I/O variables ==================================
        CLASS(Csplines),INTENT(INOUT) :: this !> Cubic splines Z_interpolation subtype
        REAL(KIND=8),DIMENSION(:),INTENT(IN) :: v !> Set of values for F(z(i))
        REAL(KIND=8),DIMENSION(:),INTENT(IN) :: z !> Set of Z values
        REAL(KIND=8),INTENT(IN) :: dz1 !> Value of second derivative at first point
        REAL(KIND=8),INTENT(IN) :: dz2 !> Value of second derivative at last point
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: filename !> file name to store interpolation coefficients
        ! Local variables
        INTEGER(KIND=4) :: n
        ! Run section =====================================
        IF (size(z)/=size(v)) THEN
           WRITE(0,*) "INTERPOL_CUBIC_SPLINES ERR: v and z arrays don't have the same size"
           WRITE(0,*) "INTERPOL_CUBIC_SPLINES ERR: Z size is:",size(z)
           WRITE(0,*) "INTERPOL_CUBIC_SPLINES ERR: V size is:",size(v)
           CALL EXIT(1)
        END IF
        n=size(z)
        ALLOCATE(this%d2vdz(n)) ! Allocate second derivatives array
        CALL DSPLIN(n,z,v,dz1,1,dz2,1,this%d2vdz)
        CALL SET_CUBIC_SPLINES(this,z,v,filename)
        RETURN
END SUBROUTINE INTERPOL_CUBIC_SPLINES

END MODULE INTERPOL1D_MOD
