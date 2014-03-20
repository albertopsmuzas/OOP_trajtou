!#######################################################
! MODULE: CRP6D_MOD
!#######################################################
!> @brief
!! Provides all structures and procedures to perform a complete
!! CRP-6D interpolation job
!
!> @warning
!! - Inherits modules CRP3D_MOD, BICSPLINES_MOD
!#######################################################
MODULE CRP6D_MOD
USE CRP3D_MOD
USE BICSPLINES_MOD
USE FOURIER1D_MOD
IMPLICIT NONE
!//////////////////////////////////////////////////
! TYPE: Cut2d
!
!-------------------------------------------------
TYPE Cut2d
   CHARACTER(LEN=30) :: alias
   CHARACTER(LEN=30) :: filename
   REAL(KIND=8) :: x,y,phi,theta
   TYPE(Bicsplines) :: interrz
   CONTAINS
      PROCEDURE,PUBLIC :: READ => READ_CUT2D
END TYPE Cut2d
!//////////////////////////////////////////////////
! TYPE: CRP6D_SITIO
!
!--------------------------------------------------
TYPE CRP6D_SITIO
   INTEGER(KIND=4) :: n2dcuts
   INTEGER(KIND=4) :: nphicuts
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE :: nphipoints
   TYPE(Fourier1d),DIMENSION(:),ALLOCATABLE :: phicut
   TYPE(Cut2d),DIMENSION(:),ALLOCATABLE :: zrcut
END TYPE CRP6D_SITIO

!/////////////////////////////////////////////////
! TYPE: CRP6D
!
!-------------------------------------------------
TYPE,EXTENDS(PES) :: CRP6D
   INTEGER(KIND=4) :: nsites
   INTEGER(KIND=4) :: natomic
   TYPE(CRP6D_Sitio),DIMENSION(:),ALLOCATABLE :: CRP6D_site
   TYPE(CRP3D),DIMENSION(:),ALLOCATABLE :: atomiccrp
   CONTAINS
      PROCEDURE,PUBLIC :: READ => READ_CRP6D
END TYPE CRP6D
CONTAINS
!#####################################################
! SUBROUTINE: READ_CUT2D
!> @brief
!! Set up a Cut2d object from file
!-----------------------------------------------------
SUBROUTINE READ_CUT2D(this,filename)
   ! Initial declarations
   USE MATHS_MOD
   USE UNITS_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(OUT) :: this
   CHARACTER(LEN=30),INTENT(IN) :: filename
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   INTEGER(KIND=4) :: nx,ny
   REAL(KIND=8) :: auxr1,auxr2
   CHARACTER(LEN=10) :: units1,units2,units3
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: z,r
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: f
   TYPE(Length) :: len1,len2
   TYPE(Energy) :: en
   TYPE(Angle) :: angl 
   ! Run section ----------------------
   this%filename = filename
   OPEN (UNIT=10,FILE=filename,STATUS="old",ACTION="read")
   READ(10,*) 
   READ(10,*) this%alias
   READ(10,*) units1,units2,units3
   READ(10,*) auxr1,auxr2
   CALL len1%READ(auxr1,units1)
   CALL len2%READ(auxr2,units1)
   CALL len1%TO_STD()
   CALL len2%TO_STD()
   this%x=len1%getvalue()
   this%y=len2%getvalue()
   READ(10,*) auxr1,auxr2
   CALL angl%READ(auxr1,units3)
   CALL angl%TO_STD()
   this%theta=angl%getvalue()
   CALL angl%READ(auxr2,units3)
   CALL angl%TO_STD()
   this%phi=angl%getvalue()
   READ(10,*) nx,ny
   ALLOCATE(r(nx))
   ALLOCATE(z(ny))
   ALLOCATE(f(nx,ny))
   DO j = 1, ny
      DO i = 1, nx
         READ(10,*) r(i),z(j),f(i,j)
         CALL len1%READ(r(i),units1)
         CALL len2%READ(z(j),units1)
         CALL en%READ(f(i,j),units2)
         CALL len1%TO_STD()
         CALL len2%TO_STD()
         CALL en%TO_STD()
         r(i)=len1%getvalue()
         z(j)=len2%getvalue()
         f(i,j)=en%getvalue()
      END DO
   END DO
   ! order array in a good way:
   DO j = 1, ny
      CALL ORDER(r,f(:,j))
   END DO
   DO i = 1, nx
      CALL ORDER(z,f(i,:))
   END DO
   !
   CALL this%interrz%READGRID(r,z,f)
   CALL this%interrz%SET_COEFF("blabla")
   CLOSE(10)
   RETURN
END SUBROUTINE READ_CUT2D
!###########################################################
!# SUBROUTINE: READ_CRP6D 
!###########################################################
!> @brief
!! Sets up a CRP6D Object
!-----------------------------------------------------------
SUBROUTINE READ_CRP6D(this,filename)
   ! Initial declarations   
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(OUT) :: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   INTEGER(KIND=4) :: natomic ! number of atomic potentials
   INTEGER(KIND=4) :: auxint 
   CHARACTER(LEN=30) :: cut2dfilename, string
   CHARACTER(LEN=12),PARAMETER :: routinename="READ_CRP6D: "
   ! Run section -----------------------
   CALL this%SET_ALIAS("CRP6D PES")
   CALL this%SET_DIMENSIONS(6)
   ! set up molecular crp
   OPEN (UNIT=180,FILE=filename,STATUS="old",ACTION="read")
   READ(180,*) ! dummy line
   READ(180,*) this%natomic
   SELECT CASE(this%natomic)
      CASE(: 0,3 :)
         WRITE(0,*) "READ_CRP6D ERR: Wrong number of atomic potentials. Allowed number: 1 or 2."
         CALL EXIT(1)
      CASE DEFAULT 
         ! do nothing
   END SELECT
   ALLOCATE(this%atomiccrp(this%natomic))
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Setting up: ",this%getalias())
   CALL VERBOSE_WRITE(routinename,"Atomic potentials found: ",this%natomic)
#endif
   DO i = 1, this%natomic
      READ(180,*) string
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Load from file: ",string)
#endif
      CALL  this%atomiccrp(i)%READ(string)
      CALL this%atomiccrp(i)%INTERPOL_Z()
   END DO
   READ(180,*) this%nsites
   SELECT CASE(this%nsites)
      CASE(: 0)
         WRITE(0,*) "READ_CRP6D ERR: Wrong natural number"
         CALL EXIT(0)
      CASE DEFAULT
         ! do nothing
   END SELECT
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"CRP6D sites found: ",this%nsites)
#endif
   ALLOCATE(this%CRP6D_site(this%nsites))
   DO i = 1, this%nsites
      READ(180,*) auxint
      this%CRP6D_site(i)%n2dcuts = auxint
      ALLOCATE(this%CRP6D_site(i)%zrcut(auxint))
      DO j = 1, auxint
         READ(180,*) string
         CALL this%CRP6D_site(i)%zrcut(j)%READ(string)
      END DO
   END DO
   CLOSE(180)
   RETURN
END SUBROUTINE READ_CRP6D
END MODULE CRP6D_MOD
