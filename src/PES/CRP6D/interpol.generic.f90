!#########################################################
! MODULE: INTERPOL_WYCKOFF_GENERIC
!> @brief
!! Provides generic routines to interpol a CRP6D PES 
!##########################################################
MODULE INTERPOL_WYCKOFF_GENERIC_MOD
   USE BICSPLINES_MOD
   USE FOURIER1D_MOD
! Initial declarations
IMPLICIT NONE
!//////////////////////////////////////////////////
! TYPE: Cut2d
!> @brief
!! Contains data for a Z,R interpolation letting X,Y,Theta,Phi
!! fixed
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
! TYPE: WYCKOFFSITIO
!
!> @brief
!! Stores all data related to a single X,Y position
!! of the molecule on the surface.
!
!> @param  n2dcuts - Number of Z,R cuts for this site
!> @param nphicuts - Number of Phi interpolations that will exist
!!                   for different Theta values. 
!> @param letter - Letter that specifies the symmetry subgroup
!!                related to this site. Its meaning changes for
!!                each planar group
!> @brief nphipoints - Array of integer numbers storing the number of points
!!                     in which each phi interpolation is based
!> 
!--------------------------------------------------
TYPE,ABSTRACT :: Wyckoffsitio
   !PRIVATE
   CHARACTER :: id
   INTEGER(KIND=4) :: n2dcuts
   INTEGER(KIND=4) :: nphicuts
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE :: nphipoints
   TYPE(Fourier1d),DIMENSION(:),ALLOCATABLE :: phicut
   TYPE(Cut2d),DIMENSION(:),ALLOCATABLE :: zrcut
   CONTAINS
      PROCEDURE,PUBLIC :: READ => READ_WYCKOFFSITIO
      PROCEDURE,PUBLIC :: SET_ID => SET_ID_WYCKOFFSITIO
END TYPE Wyckoffsitio
!/////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_ID_WYCKOFFSITIO 
!###########################################################
!> @brief
!! Set the correct id (letter) for this wyckoff sitio
!-----------------------------------------------------------
SUBROUTINE SET_ID_WYCKOFFSITIO(this,id)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Wyckoffsitio),INTENT(INOUT) :: this
   CHARACTER,INTENT(IN) :: id
   ! Run section
   this%id=id
   RETURN
END SUBROUTINE SET_ID_WYCKOFFSITIO
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
!# SUBROUTINE: READ_WYCKOFFSITE 
!###########################################################
!> @brief
!! Initialize a Wyckoffsitio from file. Instead of a filename
!! it should be given a unit number because this routine will read
!! a portion of an already opened input file for CRP6D
!
!> @warning
!! - Unit should be opened
!! - Can only read integer numbers with format '(I2)'. In principle,
!!   this should be enough, but can be modified if necessary.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 20/03/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE READ_WYCKOFFSITIO(this,u)
   USE DEBUG_MOD
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Wyckoffsitio),INTENT(INOUT) :: this
   INTEGER(KIND=4),INTENT(IN) :: u
   ! Local variables
   CHARACTER(LEN=30) :: filename
   INTEGER(KIND=4) :: i ! counters
   CHARACTER(LEN=19) :: routinename="READ_WYCKOFFSITIO: "
   ! Run section
   ! Read first line in several steps
   READ(u,'(I2)',advance="no") this%n2dcuts
   ALLOCATE(this%zrcut(this%n2dcuts))
   READ(u,'(I2)',advance="no") this%nphicuts
   ALLOCATE(this%phicut(this%nphicuts))
   ALLOCATE(this%nphipoints(this%nphicuts))
   READ(u,*) this%nphipoints(:)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"n2dcuts: ",this%n2dcuts)
   CALL VERBOSE_WRITE(routinename,"nphicuts: ",this%nphicuts)
   CALL VERBOSE_WRITE(routinename,this%nphipoints(:))
#endif
   DO i = 1, this%n2dcuts
      READ(u,*) filename
      CALL this%zrcut(i)%READ(filename)
   END DO
   RETURN
END SUBROUTINE READ_WYCKOFFSITIO

END MODULE INTERPOL_WYCKOFF_GENERIC_MOD
