!#########################################################
! MODULE SURFACE_MOD
!
!> @brief 
!! Should contain everything related with periodic 2D surfaces
!##########################################################
MODULE SURFACE_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////////////////////////
! TYPE : Atom_list
!> @brief
!! Auxiliary type data that contains info about a list of atoms of the same type
!
!> @param n - Number of atoms in this list
!> @param alias - Its periodic table symbol
!> @param atom - Matrix of positions of each atom in this list
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 03/Feb/2014
!> @version 1.1
!-------------------------------------------------------------------------------------
TYPE,PRIVATE :: Atom_list
	INTEGER(KIND=4) :: n ! number of atoms in this list
	CHARACTER(LEN=2) :: alias ! atom name, periodic table
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: atom
END TYPE
!/////////////////////////////////////////////////////////////////////////////////////
! TYPE: Surface
!> @brief
!! Class that contains all data needed from a 2D periodic surface
!
!> @param alias - Human-fiendly name for this surface
!> @param filename - Input file that contains all surface information
!> @param symmetry_alias - Name for the symmetry group of this surface
!> @param units -  Units in which distances are stored
!> @param initialized - Controls if this class was initialized or not
!> @param symmetry - Integer that stores standard symmetry group ID
!> @param norm_s1, norm_s2 - Norms of surface vectors s1 & s2
!> @param surf2cart_mtrx - Transformation matrix from surface coordinates to auxiliary cartesian
!> @param cart2surf_mtrx - Transformation matrix from auxiliary cartesian coordinates to surface
!> @param surfunit2cart_mtrx - Transformation matrix from unit surface corrdinates to auxiliary cartesians
!> @param cart2surfunit_mtrx - Transformation matrix from auxiliary cartesians to unit surface coordinates
!> @param recip2cart_mtrx - Transformation matrix from reciprocal lattice coordinates to auxiliary cartesians
!> @param cart2recip_mtrx - Transformation matrix from auxiliary cartesian coordinates to reciprocal lattice
!> @param metricsurf_mtrx - Metric matrix for surface coordinates
!> @param diff_atoms - Number of different types of atoms in the surface
!> @param atomtype - Array of lists of atoms
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 03/Feb/2014
!> @version 1.0
! 
!> @see atom_list
!-------------------------------------------------------------------------------------
TYPE Surface
PRIVATE
	CHARACTER(LEN=30) :: alias
	CHARACTER(LEN=30) :: filename
	CHARACTER(LEN=10) :: symmetry_alias
	LOGICAL :: initialized=.FALSE.
	INTEGER(KIND=4) :: symmetry
	REAL(KIND=8),DIMENSION(2,2) :: surf2cart_mtrx 
	REAL(KIND=8),DIMENSION(2,2) :: cart2surf_mtrx
	REAL(KIND=8),DIMENSION(2,2) :: surfunit2cart_mtrx 
	REAL(KIND=8),DIMENSION(2,2) :: cart2surfunit_mtrx 
	REAL(KIND=8),DIMENSION(2,2) :: recip2cart_mtrx 
	REAL(KIND=8),DIMENSION(2,2) :: cart2recip_mtrx 
	INTEGER(KIND=4) :: diff_atoms 
	TYPE(Atom_list),DIMENSION(:),ALLOCATABLE,PUBLIC :: atomtype 
	REAL(KIND=8),DIMENSION(2,2),PUBLIC :: metricsurf_mtrx 
   CHARACTER(LEN=10),PUBLIC :: units
	REAL(KIND=8),PUBLIC :: norm_s1, norm_s2 
CONTAINS
	! Initiallize
	PROCEDURE, PUBLIC :: INITIALIZE
	! Operations block
	PROCEDURE,PUBLIC :: surf2cart
	PROCEDURE,PUBLIC :: cart2surf
	PROCEDURE,PUBLIC :: surfunit2cart
	PROCEDURE,PUBLIC :: cart2surfunit
	PROCEDURE,PUBLIC :: recip2cart
	PROCEDURE,PUBLIC :: cart2recip
   PROCEDURE,PUBLIC :: project_unitcell
   PROCEDURE,PUBLIC :: project_iwscell
   ! Enquire block
   PROCEDURE, PUBLIC :: is_initialized
END TYPE
! MODULE CONTAINS 
CONTAINS
!###########################################################
!# FUNCTION: is_initialized 
!###########################################################
! - Check if surface type is already initialized
!-----------------------------------------------------------
LOGICAL FUNCTION is_initialized(surf) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface), INTENT(IN) :: surf
   ! Run section
   is_initialized=surf%initialized
   RETURN
END FUNCTION is_initialized
  
!###############################################################################
!# SUBROUTINE: INITIALIZE ######################################################
!###############################################################################
!> @brief
!! Initializes surface from file @b filename
!-------------------------------------------------------------------------------
SUBROUTINE INITIALIZE(surf,filename)
   USE UNITS_MOD
   USE MATHS_MOD 
   USE CONSTANTS_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   ! I/O Variables -----------------------------------------------
   CLASS(Surface), INTENT(INOUT) :: surf
   CHARACTER(LEN=*), INTENT(IN) :: filename
   ! Local variables ---------------------------------------------
   INTEGER :: i,j ! Counters
   INTEGER :: control
   CHARACTER(LEN=17), PARAMETER :: routinename = "INITIALIZE_SURF: "
   TYPE(length) :: len
   REAL(KIND=8), DIMENSION(2,2) :: aux_r
   ! Run section --------------------------------------------
      surf%filename=filename
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Initializing a new surface")
   CALL VERBOSE_WRITE(routinename,"File: ",surf%filename)
#endif
   IF (.NOT.surf%is_initialized()) THEN
      surf%initialized=.FALSE.
      OPEN(10,FILE=surf%filename,STATUS="old")
      READ(10,*) ! dummy line
      READ(10,*) surf%alias
      READ(10,*) surf%units
      ! Read surface basis vectors in cartesian coordinates.
      ! They should be written horizontally
      ! Store results in a.u.
      ! Set surface coordinates to cartesian change matrix and metadata
      DO i=1,2
         READ(10,*)  aux_r(i,:)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Surface vector: ",i)
         CALL VERBOSE_WRITE(routinename,"Units: ",surf%units)
         CALL VERBOSE_WRITE(routinename,aux_r(i,:))
#endif
         DO j = 1, 2
            CALL len%READ(aux_r(i,j),surf%units)
            CALL len%TO_STD()
            aux_r(i,j)=len%getvalue()
         END DO
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Units: a.u.")
         CALL VERBOSE_WRITE(routinename,aux_r(i,:))
#endif
      END DO
      aux_r=transpose(aux_r)
      surf%surf2cart_mtrx=aux_r
#ifdef DEBUG
      CALL DEBUG_WRITE(routinename,"Surf2cart matrix calculated: ")
      DO i=1,2
         CALL DEBUG_WRITE(routinename,surf%surf2cart_mtrx(i,:))
      END DO
#endif
      ! Storing actual surface units (now, everything is in atomic units)
      surf%units = "au"
      ! Read number of different atom types in the surface
      READ(10,*) surf%diff_atoms
#ifdef DEBUG
         CALL DEBUG_WRITE(routinename,"Different atoms found in the basis: ",surf%diff_atoms)
#endif
      ! Set the the atomic basis, i.e., coordinates for non equivalent  atoms inside
      ! the unit cell
      ALLOCATE(surf%atomtype(1:surf%diff_atoms))
      DO i=1, surf%diff_atoms
         READ(10,*) control, surf%atomtype(i)%alias, surf%atomtype(i)%n
         IF (control.NE.i) THEN
            WRITE(0,*) "INITIALIZE_SURF ERR: Atom type definitions are not in the correct order"
            CALL EXIT(1)
         END IF
         ALLOCATE(surf%atomtype(i)%atom(surf%atomtype(i)%n,2))
      END DO
      DO i=1, surf%diff_atoms
            DO j=1, surf%atomtype(i)%n
               READ(10,*) control,surf%atomtype(i)%atom(j,1),surf%atomtype(i)%atom(j,2)
               IF (control.NE.i) THEN
                  WRITE(0,*) "INITIALIZE_SURF ERR: Atom type definitions are not in the correct order"
                  CALL EXIT(1)
               END IF
#ifdef DEBUG
               CALL DEBUG_WRITE(routinename,"Atom type: (surf. coord.)",i)
               CALL DEBUG_WRITE(routinename,surf%atomtype(i)%atom(j,1),surf%atomtype(i)%atom(j,2))
#endif
            END DO
      END DO
      CLOSE(10)
      ! Set metric matrix associated with surface coordinates.
      surf%metricsurf_mtrx=MATMUL(TRANSPOSE(surf%surf2cart_mtrx),surf%surf2cart_mtrx)
#ifdef DEBUG
          CALL VERBOSE_WRITE(routinename,"Surface metric matrix calculated:")
          DO i = 1,2
            CALL DEBUG_WRITE(routinename, surf%metricsurf_mtrx(i,1), surf%metricsurf_mtrx(i,2)) 
          END DO
#endif
      ! Set Matrix: from auxiliar cartesian coordinates to surface coordinates
      CALL INV_MTRX(2,surf%surf2cart_mtrx,surf%cart2surf_mtrx)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Cart2surf matrix calculated:")
         DO i = 1, 2
            CALL DEBUG_WRITE(routinename,surf%cart2surf_mtrx(i,1),surf%cart2surf_mtrx(i,2))
         END DO
#endif
      ! Set primitive vectors norms
      surf%norm_s1=DSQRT(DOT_PRODUCT(surf%surf2cart_mtrx(1:2,1),surf%surf2cart_mtrx(1:2,1)))
      surf%norm_s2=DSQRT(DOT_PRODUCT(surf%surf2cart_mtrx(1:2,2),surf%surf2cart_mtrx(1:2,2)))
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Modulus S1: ", surf%norm_s1)
         CALL VERBOSE_WRITE(routinename,"Modulus S2: ", surf%norm_s2)
#endif
      ! Set Matrix: from normalized surface coordinates to auxiliar cartesian coordinates
      FORALL(i=1:2) surf%surfunit2cart_mtrx(i,1)=surf%surf2cart_mtrx(i,1)/surf%norm_s1
      FORALL(i=1:2) surf%surfunit2cart_mtrx(i,2)=surf%surf2cart_mtrx(i,2)/surf%norm_s2
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Surfunit2cart matrix calculated:")
         DO i = 1, 2
            CALL DEBUG_WRITE(routinename,surf%surfunit2cart_mtrx(i,1),surf%surfunit2cart_mtrx(i,2))
         END DO
#endif
      ! Set Matrix: from auxiliar cartesian coordinates to normalized surface coordinates
      CALL INV_MTRX(2,surf%surfunit2cart_mtrx,surf%cart2surfunit_mtrx)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Cart2surfunit matrix calculated:")
         DO i = 1, 2
            CALL DEBUG_WRITE(routinename,surf%cart2surfunit_mtrx(i,1),surf%cart2surfunit_mtrx(i,2))
         END DO
#endif
      ! Set Matrix: form reciprocal space to auxiliar cartesian coordinates
      surf%recip2cart_mtrx=TRANSPOSE(surf%cart2surf_mtrx)
      !FORALL(i=1:2,j=1:2) surf%recip_mtrx(i,j)=2*PI*surf%recip_mtrx(i,j)
      surf%recip2cart_mtrx=2.D0*PI*surf%recip2cart_mtrx
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Recip2cart matrix calculated:")
         DO i = 1, 2
            CALL DEBUG_WRITE(routinename,surf%recip2cart_mtrx(i,1),surf%recip2cart_mtrx(i,2))
         END DO
#endif
      ! Set Matrix: from auxiliar cartesian coordinates to reciprocal space
      CALL INV_MTRX(2,surf%recip2cart_mtrx,surf%cart2recip_mtrx)
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Cart2recip matrix calculated:")
         DO i = 1, 2
            CALL DEBUG_WRITE(routinename,surf%cart2recip_mtrx(i,1),surf%cart2recip_mtrx(i,2))
         END DO
#endif
      surf%initialized=.TRUE. 
   ELSE
      WRITE(0,*) "INITIALIZE_SURF ERR: surface already initialized"
      CALL EXIT(1)
   END IF
   RETURN
END SUBROUTINE INITIALIZE
!###########################################################
!# FUNCTION: cart2surf 
!###########################################################
! - Goes from auxiliar cartesian coordinates to surface coordinates
!-----------------------------------------------------------
FUNCTION cart2surf(surf,r) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface), INTENT(IN) :: surf
   REAL(KIND=8), DIMENSION(2), INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: cart2surf
   ! Run section
   cart2surf=matmul(surf%cart2surf_mtrx,r)
   RETURN
END FUNCTION cart2surf
!###########################################################
!# FUNCTION: surf2cart 
!###########################################################
! - Goes from surface to auxiliar cartesian coordinates
!-----------------------------------------------------------
FUNCTION surf2cart(surf,r) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface),INTENT(IN) :: surf
   REAL(KIND=8), DIMENSION(2),INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: surf2cart
   ! Run section
   surf2cart=matmul(surf%surf2cart_mtrx,r)
   RETURN
END FUNCTION surf2cart
!###########################################################
!# FUNCTION: surfunit2cart 
!###########################################################
! - Goes from unit surface to auxiliar cartesian coordinates
!-----------------------------------------------------------
FUNCTION surfunit2cart(surf,r) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface),INTENT(IN) :: surf
   REAL(KIND=8), DIMENSION(2),INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: surfunit2cart
   ! Run section
   surfunit2cart=matmul(surf%surfunit2cart_mtrx,r)
   RETURN
END FUNCTION surfunit2cart
!###########################################################
!# FUNCTION: cart2surfunit
!###########################################################
! - Goes from auxiliar cartesian coordinates to surface coordinates
!-----------------------------------------------------------
FUNCTION cart2surfunit(surf,r) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface), INTENT(IN) :: surf
   REAL(KIND=8), DIMENSION(2), INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: cart2surfunit
   ! Run section
   cart2surfunit=matmul(surf%cart2surfunit_mtrx,r)
   RETURN
END FUNCTION cart2surfunit
!###########################################################
!# FUNCTION: cart2recip
!###########################################################
! - Goes from auxiliar cartesian to reciprocal space coordinates
!-----------------------------------------------------------
FUNCTION cart2recip(surf,r) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface),INTENT(IN) :: surf
   REAL(KIND=8), DIMENSION(2),INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: cart2recip
   ! Run section
   cart2recip=matmul(surf%cart2recip_mtrx,r)
   RETURN
END FUNCTION cart2recip
!###########################################################
!# FUNCTION: recip2cart 
!###########################################################
! - Goes from auxiliar cartesian to surface coordinates
!-----------------------------------------------------------
FUNCTION recip2cart(surf,r) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Surface), INTENT(IN) :: surf
   REAL(KIND=8), DIMENSION(2), INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: recip2cart
   ! Run section
   recip2cart=matmul(surf%recip2cart_mtrx,r)
   RETURN
END FUNCTION recip2cart
!################################################################
!# SUBROUTINE: PROJECT_UNITCELL #################################
!################################################################
!> @brief
!! Projects 2D point into the C4v unit cell
!
!> @warning
!! - Input/output in cartesian coordinates (r)
!----------------------------------------------------------------
FUNCTION project_unitcell(surf,r)
	IMPLICIT NONE
	! I/O variables
	CLASS(Surface),INTENT(IN) :: surf
	REAL(KIND=8),DIMENSION(2),INTENT(IN) :: r
	! Local variables
   REAL(KIND=8),DIMENSION(2) :: project_unitcell
	REAL(KIND=8), DIMENSION(2) :: aux
	INTEGER :: i ! counters
	! HEY, HO! LET'S GO !!! ----------------------
   project_unitcell = surf%cart2surf(r)
	FORALL (i=1:2) 
		aux(i)=DFLOAT(INT(project_unitcell(i)))
		project_unitcell(i)=project_unitcell(i)-aux(i)
	END FORALL
   project_unitcell = surf%surf2cart(project_unitcell)
	RETURN
END FUNCTION project_unitcell
!################################################################
! SUBROUTINE: project_iwscell ###################################
!################################################################
!> @brief
!! Projects R into Irreducible WS cell. Needs information from
!! surface main vectors.
!
!> @param[in] surf - Surface specifications
!> @param[in] x - 2D point
!
!> @warning
!! - r is in cartesian coordinates (Input and output)
!! - Only C4v symmetry
!----------------------------------------------------------------
FUNCTION project_iwscell(surf,x)
   USE CONSTANTS_MOD
	IMPLICIT NONE
	! I/O variables
	CLASS(Surface) , INTENT(IN) :: surf
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: x
	! Local variables
	REAL*8,DIMENSION(2) :: r
	REAL*8,DIMENSION(2,2) :: Proj_x, Proj_y, Rot
	REAL*8,DIMENSION(2) :: aux
	REAL*8 :: angle, radius, alpha
	INTEGER :: i ! counters
   REAL(KIND=8),DIMENSION(2) :: project_iwscell
	! HEY, HO! LET'S GO! ------------------
	! Go to surface coordinates
   r = x
	r = surf%cart2surf(r)
	FORALL (i=1:2) 
		aux(i)=DFLOAT(INT(r(i)))
		r(i)=r(i)-aux(i)
	END FORALL
	! Now, r vector is inside the unit cell. Let's define this vector
	! taking as the origin the center of the cell (in surface units is 0.5,0.5):
	FORALL (i=1:2) r(i)=r(i)-0.5D0
	! Calculate angle
	IF (r(1).EQ.0.D0) THEN
		angle = PI/2.D0
	ELSE
		angle = DATAN(DABS(r(2))/DABS(r(1))) ! Only angles from 0 to pi rad.
	END IF
	radius = DSQRT(r(1)**2.D0+r(2)**2.D0) ! Radius
	r(1)=radius*DCOS(angle)
	r(2)=radius*DSIN(angle)
	alpha = PI/2.D0 ! 90 deg.
	! 2D Projectors ========================
	!Rotation
	Rot(1,1)=DCOS(alpha)
	Rot(1,2)=-DSIN(alpha)
	Rot(2,1)=DSIN(alpha)
	Rot(2,2)=DCOS(alpha)
	! Projector X
	Proj_x(1,1)=1.D0
	Proj_x(1,2)=0.D0
	Proj_x(2,1)=0.D0
	Proj_x(2,2)=-1.D0
	! Projector Y
	Proj_y(1,1)=-1.D0
	Proj_y(1,2)=0.D0
	Proj_y(2,1)=0.D0
	Proj_y(2,2)=1.D0
	! Project onto IWS cell
	IF ((angle.LE.(PI/4.D0)).AND.(angle.GE.0.D0)) THEN
		r = MATMUL(Proj_y,r)
		r = MATMUL(Rot, r)
	ELSE IF ((angle.LE.(PI/2.D0)).AND.(angle.GT.(PI/4.D0))) THEN
		r = MATMUL(Proj_y, r)
		r = MATMUL(Proj_x, r)
	ELSE 
		WRITE(*,*) "ERR GET_V_XYZ: Incorrect angle value."
		WRITE(*,*) "angle : ", angle
		STOP
	END IF
	FORALL (i=1:2) r(i)=r(i)+0.5D0
	! Go to cartesian coordinates
   r = surf%surf2cart(r)
   project_iwscell = r
	RETURN
END FUNCTION project_iwscell
END MODULE SURFACE_MOD
