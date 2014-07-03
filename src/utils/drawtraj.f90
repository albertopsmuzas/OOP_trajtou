!#########################################################
! MODULE: DRAWTRAJ_MOD
!> @brief
!! Module to analyze OUTtrajxxxxx.out files in order to obtain
!! graphic representations
!
!> @todo
!! - Generalize for CRP6D and other formats
!##########################################################
MODULE DRAWTRAJ_MOD
   USE SURFACE_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: Drawtraj
!> @brief
!! All utils and information to create pgraphic representation of a
!! trajectory
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!----------------------------------------------------------------
TYPE :: Drawtraj
PRIVATE
   CHARACTER(LEN=3) :: format_out
   TYPE(Surface):: surf
   INTEGER(KIND=4) :: npattern
   INTEGER(KIND=4) :: nprojectile
CONTAINS
   ! Initialization block
   PROCEDURE,PUBLIC :: INITIALIZE => INITIALIZE_DRAWTRAJ
   ! Tools block
   PROCEDURE,PUBLIC :: DRAW => DRAW_DRAWTRAJ
END TYPE Drawtraj
CONTAINS
!###########################################################
!# SUBROUTINE: INITIALIZE_DRAWTRAJ 
!###########################################################
!> @brief
!! brief description
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INITIALIZE_DRAWTRAJ(this,format_out,filename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Drawtraj),INTENT(OUT):: this
   CHARACTER(LEN=3),INTENT(IN) :: format_out
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables
   ! Run section
   this%format_out=format_out
   CALL this%surf%INITIALIZE(filename)
   this%npattern=sum(this%surf%atomtype(:)%n)
   RETURN
END SUBROUTINE INITIALIZE_DRAWTRAJ
!###########################################################
!# SUBROUTINE: DRAW_DRAWTRAJ  
!###########################################################
!> @brief
!! Draw a trajectory adding a pattern. The pattern is what is
!! actually moving
!
!> @warning
!! Can be only used with CRP6D dynamics output
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE DRAW_DRAWTRAJ (this,order,dynamicsfilename,outputfilename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Drawtraj),INTENT(IN):: this
   CHARACTER(LEN=*),INTENT(IN) :: dynamicsfilename,outputfilename
   INTEGER(KIND=4),INTENT(IN) :: order
   ! IMPORTANT: units used to read and write
   INTEGER(KIND=4) :: runit=38, wunit=173
   ! Local variables
   TYPE(Surface) :: surf_aux
   INTEGER(KIND=4) :: i  ! counters
   REAL(KIND=8),DIMENSION(15) :: dummy
   REAL(KIND=8),DIMENSION(3) :: init_xcm,xa,xb,xcm
   REAL(KIND=8),DIMENSION(2) :: dx
   REAL(KIND=8) :: natoms
   INTEGER(KIND=4) :: npattern
   INTEGER(KIND=4) :: ioerr
   ! Run section
   SELECT CASE(order)
      CASE(0)
         npattern=this%npattern
      CASE(1:)
         npattern=this%npattern
         DO i = 1, order
            npattern=npattern+8*i*this%npattern
         END DO
      CASE DEFAULT
         WRITE(0,*) "DRAW_DRAWTRAJ ERR: wrong order number"
   END SELECT
   ! Deal with dynamics file
   OPEN (runit,FILE=dynamicsfilename,STATUS="old",ACTION="read")
   OPEN (wunit,FILE=outputfilename,STATUS="replace",ACTION="write")
   READ(runit,*) ! dummy line
   READ(runit,*) ! dummy line
   READ(runit,*) ! dummy line
   READ(runit,*) ! dummy line
   READ(runit,*) ! dummy line
   i=0
   READ(runit,*) dummy(1:6),init_xcm,dummy(7:15),xa,xb 
   WRITE(wunit,*) npattern+2
   WRITE(wunit,*) 
   CALL this%surf%PRINT_PATTERN(wunit,order,this%format_out)
   WRITE(wunit,*) "H ",xa
   WRITE(wunit,*) "H ",xb
   DO 
      i=i+1
      surf_aux=this%surf
      ! read section
      READ(runit,*,iostat=ioerr) dummy(1:6),xcm,dummy(7:15),xa,xb 
      dx(1:2)=xcm(1:2)-init_xcm(1:2)
      ! Check if EOF reached
      SELECT CASE(ioerr/=0)
         CASE(.TRUE.)
            WRITE(*,*) "DRAW_DRAWTRAJ: EOF reached. Job Finished."
            WRITE(*,*) "Input file: ",dynamicsfilename
            WRITE(*,*) "Output file: ",outputfilename
            WRITE(*,*) "IOSTAT code: ",ioerr
            EXIT
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      ! Write section
      WRITE(wunit,*) npattern+2
      WRITE(wunit,*) ! write dummy line
      CALL surf_aux%MOVE_PATTERN(-dx)
      CALL surf_aux%PRINT_PATTERN(wunit,order,this%format_out)
      xa(1:2)=xa(1:2)-dx
      xb(1:2)=xb(1:2)-dx
      WRITE(wunit,*) "H ",xa
      WRITE(wunit,*) "H ",xb
   END DO
   CLOSE(runit)
   CLOSE(wunit)
   RETURN
END SUBROUTINE DRAW_DRAWTRAJ 
END MODULE DRAWTRAJ_MOD
