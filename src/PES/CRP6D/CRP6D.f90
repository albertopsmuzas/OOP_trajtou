!#######################################################
! MODULE: CRP6D_MOD
!#######################################################
!> @brief
!! Provides all structures and procedures to perform a complete
!! CRP-6D interpolation job. This method only works systematically with diatomic
!! molecules
!
!> @warning
!! - Inherits modules CRP3D_MOD, BICSPLINES_MOD
!#######################################################
MODULE CRP6D_MOD
USE CRP3D_MOD
USE INTERPOL_WYCKOFF_P4MM_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////
! TYPE: CRP6D
!
!-------------------------------------------------
TYPE,EXTENDS(PES) :: CRP6D
   INTEGER(KIND=4) :: nsites
   INTEGER(KIND=4) :: natomic
   CLASS(Wyckoffsitio),DIMENSION(:),ALLOCATABLE :: wyckoffsite
   TYPE(CRP3D),DIMENSION(:),ALLOCATABLE :: atomiccrp
   CONTAINS
      PROCEDURE,PUBLIC :: READ => READ_CRP6D
END TYPE CRP6D
CONTAINS
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
   CHARACTER,DIMENSION(:),ALLOCATABLE :: letter
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
   ! Prepare atomic potentials
   DO i = 1, this%natomic
      READ(180,*) string
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Load from file: ",string)
#endif
      CALL  this%atomiccrp(i)%READ(string)
      CALL this%atomiccrp(i)%INTERPOL_Z()
   END DO
   ! Check that all atomic potentials are based on the same surface
   DO i = 1, this%natomic
      SELECT CASE(this%atomiccrp(1)%surf%tellfilename()==this%atomiccrp(i)%surf%tellfilename())
         CASE(.FALSE.)
            WRITE(0,*) "READ_CRP6D ERR: not all atomic potentials have the same surface file"
            CALL EXIT(1)
         CASE(.TRUE.)
            ! do nothing
      END SELECT
   END DO
   ! Read number of wyckoff sites and its letters
   READ(180,'(I2)',advance="no") this%nsites
   ALLOCATE(letter(this%nsites))
   READ(180,*) letter(:)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Wyckoff sites found: ",this%nsites)
   CALL VERBOSE_WRITE(routinename,"Wyckoff letters map:")
   CALL VERBOSE_WRITE(routinename,letter(:))
#endif
! Allocate with the correct type (symmetry) all wyckoff sites
   SELECT CASE(this%atomiccrp(1)%surf%tellsymmlabel())
      CASE("p4mm")
         ALLOCATE(Wyckoffp4mm::this%wyckoffsite(this%nsites))
#ifdef DEBUG
         CALL VERBOSE_WRITE(routinename,"Allocated p4mm Wyckoff sites")
#endif
      CASE DEFAULT
         WRITE(0,*) "READ_CRP6D ERR: surface symmetry is not implemented yet"
         WRITE(0,*) "Good luck!"
         CALL EXIT(1)
   END SELECT
   ! Now that we have the correct type:
   DO i = 1, this%nsites
      CALL this%wyckoffsite(i)%READ(180)
      CALL this%wyckoffsite(i)%SET_ID(letter(i))
   END DO
   CLOSE(180)
   RETURN
END SUBROUTINE READ_CRP6D
END MODULE CRP6D_MOD
