!#######################################################
! MODULE: CRP6D_MOD
!#######################################################
!> @brief
!! Provides all structures and procedures to perform a complete
!! CRP-6D interpolation job
!
!> @warning
!! - Inherits modules CRP_MOD, BICSPLINES_MOD
!#######################################################
MODULE CRP6D_MOD
USE CRP_MOD
USE BICSPLINES_MOD
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
!/////////////////////////////////////////////////
! TYPE: CRP6D
!
!-------------------------------------------------
TYPE,EXTENDS(PES) :: CRP6D
   INTEGER(KIND=4) :: n
   TYPE(Cut2d),DIMENSION(:),ALLOCATABLE :: corte2d 
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
   CHARACTER(LEN=10) :: units1,units2,units3
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: z,r
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: f
   TYPE(Length) :: len1,len2
   TYPE(Energy) :: en
   ! Run section ----------------------
   this%filename = filename
   OPEN (UNIT=10,FILE=filename,STATUS="old",ACTION="read")
   READ(10,*) 
   READ(10,*) this%alias
   READ(10,*) this%x,this%y
   READ(10,*) this%theta,this%phi
   READ(10,*) units1,units2,units3
   READ(10,*) nx,ny
   ALLOCATE(r(nx))
   ALLOCATE(z(ny))
   ALLOCATE(f(nx,ny))
   DO j = 1, ny
      DO i = 1, nx
         READ(10,*) r(i),z(j),f(i,j)
         CALL len1%READ(r(i),units1)
         CALL len2%READ(z(j),units2)
         CALL en%READ(f(i,j),units3)
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
   INTEGER(KIND=4) :: i ! counters
   CHARACTER(LEN=30) :: cut2dfilename
   CHARACTER(LEN=12),PARAMETER :: routinename="READ_CRP6D: "
   ! Run section -----------------------
   CALL this%SET_ALIAS("CRP6D PES")
   CALL this%SET_DIMENSIONS(6)
   OPEN (UNIT=11,FILE=filename,STATUS="old",ACTION="read")
   READ(11,*) ! dummy line
   READ(11,*) this%n
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Setting up: ",this%getalias())
   CALL VERBOSE_WRITE(routinename,"(r-Z) 2D cuts found: ",this%n)
#endif
   ALLOCATE(this%corte2d(this%n))
   DO i = 1, this%n
      READ(11,*) cut2dfilename
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"2D cut input file found: ",cut2dfilename)
#endif
      CALL this%corte2d(i)%READ(cut2dfilename)
   END DO
   CLOSE(11)
   RETURN
END SUBROUTINE READ_CRP6D
END MODULE CRP6D_MOD
