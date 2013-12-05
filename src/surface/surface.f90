MODULE SURFACE_MOD
!=====================================================================================================
! Surface characteristics derived data 
TYPE Surface
PRIVATE
	CHARACTER(LEN=30) :: alias ! name of this surface
	CHARACTER(LEN=30) :: filename ! file that contains all surface information.
	CHARACTER(LEN=10) :: symmetry_alias ! name of the symmetry
	LOGICAL :: initialized=.FALSE.
	INTEGER :: symmetry ! standard id for the symmetry
	REAL(KIND=8) :: norm_s1, norm_s2, norm_s3 ! 2-norm of surface vectors s1&s2. s3 is perpendicular to s1 & s2
	REAL(KIND=8), DIMENSION(3,3) :: surf2cart_mtrx ! surface coord to cartesians
	REAL(KIND=8), DIMENSION(3,3) :: cart2surf_mtrx ! Cartesian to surface coord.
	REAL(KIND=8), DIMENSION(3,3) :: surfunit2cart_mtrx ! Unit surface coord. to cartesians.
	REAL(KIND=8), DIMENSION(3,3) :: cart2surfunit_mtrx ! Cartesian to unit surface coord.
	REAL(KIND=8), DIMENSION(3,3) :: recip2cart_mtrx ! Reciprocal to cartesian coord.
	REAL(KIND=8), DIMENSION(3,3) :: cart2recip_mtrx ! Cartesian to reciprocal coord.
	REAL(KIND=8), DIMENSION(3,3) :: metricsurf_mtrx ! Metric matrix for surface coord.
	INTEGER :: diff_atoms ! Number of different atoms in the surface
	TYPE(Atom_list), DIMENSION(:), ALLOCATABLE :: atomtype ! List of atoms
CONTAINS
	! Initiallize
	PROCEDURE, PUBLIC :: INITIALIZE
	! Operations block
	PROCEDURE, PUBLIC :: surf2cart
	PROCEDURE, PUBLIC :: cart2surf
	PROCEDURE, PUBLIC :: surfunit2cart
	PROCEDURE, PUBLIC :: cart2surfunit
	PROCEDURE, PUBLIC :: recip2cart
	PROCEDURE, PUBLIC :: cart2recip
END TYPE
!====================================================================================================
! Atom_list derived data
TYPE Atom_list
PRIVATE
	CHARACTER(LEN=2) :: alias ! atom name
	INTEGER :: n ! number of atoms in this list
	TYPE(Point_2D), DIMENSION(:), ALLOCATABLE :: atom
END TYPE
! MODULE CONTAINS 
CONTAINS
!###############################################################################
!# SUBROUTINE: INITIALIZE ######################################################
!###############################################################################
! - Initializes surface from file filename
!-------------------------------------------------------------------------------
SUBROUTINE INITIALIZE(surf,filename)
        IMPLICIT NONE
        ! I/O Variables -----------------------------------------------
        CLASS(Surface), INTENT(INOUT) :: surf
        ! Local variables ---------------------------------------------
        INTEGER :: i,j ! Counters
        INTEGER :: control
        CHARACTER(LEN=11), PARAMETER :: routinename = "READ_SURF: "
        TYPE(Quantity), DIMENSION(3,3) :: quant
        TYPE(Quantity) :: quant1,quant2
        ! Run section --------------------------------------------
        OPEN(10,FILE=surf%filename,STATUS="old")
        READ(10,*) ! Should be a dummy line, just with some aclarations
        ! Read surface units
        READ(10,*) surf%units
        ! Set surface coordinates to cartesian change matrix
        DO i=1,3
                READ(10,*)  quant(i,:)%mag
        END DO
        ! Set metadata
        FORALL (i=1:3,j=1:3)
                quant(i,j)%units = surf%units
                quant(i,j)%kind = "length"
        END FORALL
        ! Convert to atomic units
        CALL VERBOSE_WRITE(routinename, "Going to a.u.")
        DO i=1,3
                DO j=1,2 ! We do not want to convert 3rd vector (should be 0 0 1)
                        CALL GO_TO_AU(quant(i,j))
                END DO
        END DO
        ! Store results in a.u.
        FORALL (i=1:3,j=1:3) surf%chg_mtrx(i,j) = quant(i,j)%mag
        CALL VERBOSE_WRITE(routinename,"Basis definition in a.u.: ")
        DO i=1,3
                CALL VERBOSE_WRITE(routinename,surf%chg_mtrx(i,1),surf%chg_mtrx(i,2),surf%chg_mtrx(i,3))
        END DO
        ! Storing actual surface units
        surf%units = "au"
        READ(10,*) surf%diff_atoms
        CALL VERBOSE_WRITE(routinename,"Different atoms found: ", surf%diff_atoms)
        ALLOCATE(surf%atomtype(1:surf%diff_atoms))
        DO i=1, surf%diff_atoms
                READ(10,*) control, surf%atomtype(i)%alias, surf%atomtype(i)%n
                IF (control.NE.i) THEN
                        WRITE(0,*) "READ_SURF ERR: Atom type definitions are not in the correct order"
                        STOP
                END IF
                ALLOCATE(surf%atomtype(i)%atom(1:surf%atomtype(i)%n))
        END DO
        DO i=1, surf%diff_atoms
                DO j=1, surf%atomtype(i)%n
                        READ(10,*) control, surf%atomtype(i)%atom(j)%x, surf%atomtype(i)%atom(j)%y
                         IF (control.NE.i) THEN
                                WRITE(0,*) "READ_SURF ERR: Atom type definitions are not in the correct order"
                                STOP
                        END IF
                END DO
        END DO
        CLOSE(10)
        ! Set metric matrix associated with surface coordinates.
        surf%metric_mtrx=MATMUL(TRANSPOSE(surf%chg_mtrx),surf%chg_mtrx)
        CALL VERBOSE_WRITE(routinename,surf%filename," -----> Done")
        ! Set inverse chg_matrix
        CALL INV_MTRX(3,surf%chg_mtrx,surf%inv_chg_mtrx)
        ! Set modules
        surf%modul_u1=DSQRT(DOT_PRODUCT(surf%chg_mtrx(1:3,1),surf%chg_mtrx(1:3,1)))
        CALL VERBOSE_WRITE(routinename,"Modulus U1: ", surf%modul_u1)
        surf%modul_u2=DSQRT(DOT_PRODUCT(surf%chg_mtrx(1:3,2),surf%chg_mtrx(1:3,2)))
        CALL VERBOSE_WRITE(routinename,"Modulus U2: ", surf%modul_u2)
        surf%modul_u3=DSQRT(DOT_PRODUCT(surf%chg_mtrx(1:3,3),surf%chg_mtrx(1:3,3)))
        CALL VERBOSE_WRITE(routinename,"Modulus U3: ", surf%modul_u3)
        ! Set normalized surface coordinates to cartesian change matrix
        DO i=1, 3
                surf%chg_mtrx_unit(1,i)=surf%chg_mtrx(i,1)/surf%modul_u1
                surf%chg_mtrx_unit(2,i)=surf%chg_mtrx(i,2)/surf%modul_u2
                surf%chg_mtrx_unit(3,i)=surf%chg_mtrx(i,3)/surf%modul_u3
        END DO
        ! Set inverse of the previous matrix 
        CALL INV_MTRX(3,surf%chg_mtrx_unit,surf%inv_chg_mtrx_unit)
        ! Set reciprocal space coordinates to cartesian change matrix
        surf%recip_mtrx=TRANSPOSE(surf%inv_chg_mtrx)
        FORALL(i=1:3,j=1:3) surf%recip_mtrx(i,j)=2*PI*surf%recip_mtrx(i,j)
        ! Set inverse of the previous transformation matrix
        CALL INV_MTRX(3, surf%recip_mtrx, surf%inv_recip_mtrx)
        RETURN
END SUBROUTINE INITIALIZE

END MODULE SURFACE_MOD
