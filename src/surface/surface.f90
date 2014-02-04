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
   CHARACTER(LEN=10) :: units
	LOGICAL :: initialized=.FALSE.
	INTEGER(KIND=4) :: symmetry
	REAL(KIND=8) :: norm_s1, norm_s2 
	REAL(KIND=8),DIMENSION(2,2),PUBLIC :: metricsurf_mtrx 
	REAL(KIND=8),DIMENSION(2,2) :: surf2cart_mtrx 
	REAL(KIND=8),DIMENSION(2,2) :: cart2surf_mtrx
	REAL(KIND=8),DIMENSION(2,2) :: surfunit2cart_mtrx 
	REAL(KIND=8),DIMENSION(2,2) :: cart2surfunit_mtrx 
	REAL(KIND=8),DIMENSION(2,2) :: recip2cart_mtrx 
	REAL(KIND=8),DIMENSION(2,2) :: cart2recip_mtrx 
	INTEGER(KIND=4) :: diff_atoms 
	TYPE(Atom_list),DIMENSION(:),ALLOCATABLE,PUBLIC :: atomtype 
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
! - Initializes surface from file filename
!-------------------------------------------------------------------------------
SUBROUTINE INITIALIZE(surf,alias,filename)
   USE UNITS_MOD
   USE LAPACKCONTROL_MOD
   USE CONSTANTS_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   ! I/O Variables -----------------------------------------------
   CLASS(Surface), INTENT(INOUT) :: surf
   CHARACTER(LEN=*), INTENT(IN) :: filename, alias
   ! Local variables ---------------------------------------------
   INTEGER :: i,j ! Counters
   INTEGER :: control
   CHARACTER(LEN=17), PARAMETER :: routinename = "INITIALIZE_SURF: "
   TYPE(length) :: len
   REAL(KIND=8), DIMENSION(2,2) :: aux_r
   ! Lapack variables
   INTEGER(KIND=4) :: lwork
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: work
   INTEGER(KIND=4) :: info
   INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipiv
   ! Run section --------------------------------------------
      surf%alias=alias
      surf%filename=filename
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Initializing a new surface")
   CALL VERBOSE_WRITE(routinename,"Name: ",surf%alias)
   CALL VERBOSE_WRITE(routinename,"File: ",surf%filename)
#endif
   IF (.NOT.surf%is_initialized()) THEN
      surf%initialized=.FALSE.
      OPEN(10,FILE=surf%filename,STATUS="old")
      READ(10,*) ! Should be a dummy line, just with some aclarations
      ! Read surface units
      READ(10,*) surf%units
      ! Read surface basis vectors in cartesian coordinates.
      ! They should be written horizontally
      ! Store results in a.u.
      ! Set surface coordinates to cartesian change matrix and metadata
      DO i=1,2
            READ(10,*)  aux_r(i,:)
            DO j = 1, 2
               CALL len%READ(aux_r(i,j),surf%units)
               CALL len%TO_STD()
               aux_r(i,j)=len%getvalue()
            END DO
      END DO
      aux_r=transpose(aux_r)
      surf%surf2cart_mtrx=aux_r
#ifdef DEBUG
         CALL DEBUG_WRITE(routinename,"Surf2cart matrix calculated: ")
         DO i=1,2
            CALL DEBUG_WRITE(routinename,surf%surf2cart_mtrx(i,1),surf%surf2cart_mtrx(i,2))
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
      ! Setting lapack work variables
      lwork=3
      ALLOCATE(work(lwork))
      ALLOCATE(ipiv(2))
      FORALL(i=1:2) ipiv(i)=i
      ! Set Matrix: from auxiliar cartesian coordinates to surface coordinates
      surf%cart2surf_mtrx=surf%surf2cart_mtrx
      CALL DGETRI(2,surf%cart2surf_mtrx,2,ipiv,work,lwork,info)
      CALL LAPACK_CHECK("DGETRI",info)
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
      surf%cart2surfunit_mtrx=surf%surfunit2cart_mtrx
      CALL DGETRI(2,surf%cart2surfunit_mtrx,2,ipiv,work,lwork,info)
      CALL LAPACK_CHECK("DGETRI",info)
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
      surf%cart2recip_mtrx=surf%recip2cart_mtrx
      CALL DGETRI(2,surf%cart2recip_mtrx,2,ipiv,work,lwork,info)
      CALL LAPACK_CHECK("DGETRI",info)
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
   REAL(KIND=8), DIMENSION(2), INTENT(INOUT) :: r
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
   REAL(KIND=8), DIMENSION(2) :: r
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
   REAL(KIND=8), DIMENSION(2) :: r
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
   REAL(KIND=8), DIMENSION(2), INTENT(INOUT) :: r
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
   REAL(KIND=8), DIMENSION(2) :: r
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
   REAL(KIND=8), DIMENSION(2), INTENT(INOUT) :: r
   ! Local variables
   REAL(KIND=8),DIMENSION(2) :: recip2cart
   ! Run section
   recip2cart=matmul(surf%recip2cart_mtrx,r)
   RETURN
END FUNCTION recip2cart
END MODULE SURFACE_MOD
