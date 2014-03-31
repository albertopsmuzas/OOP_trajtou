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
USE EXTRAPOL_TO_VACUUM_MOD
USE WYCKOFF_P4MM_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////
! TYPE: CRP6D
!
!-------------------------------------------------
TYPE,EXTENDS(PES) :: CRP6D
   INTEGER(KIND=4) :: nsites
   INTEGER(KIND=4) :: natomic
   LOGICAL :: is_smooth=.FALSE.
   CLASS(Wyckoffsitio),DIMENSION(:),ALLOCATABLE :: wyckoffsite
   TYPE(CRP3D),DIMENSION(:),ALLOCATABLE :: atomiccrp
   TYPE(Vacuumpot) :: farpot
   CONTAINS
      ! Initialization block
      PROCEDURE,PUBLIC :: READ => READ_CRP6D
      ! Tools block
      PROCEDURE,PUBLIC :: SMOOTH => SMOOTH_CRP6D
      PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_CRP6D
      PROCEDURE,PUBLIC :: RAWINTERPOL => RAWINTERPOL_CRP6D
      PROCEDURE,PUBLIC :: EXTRACT_VACUUMSURF => EXTRACT_VACUUMSURF_CRP6D
END TYPE CRP6D
CONTAINS
!###########################################################
!# SUBROUTINE: INTERPOL_CRP6D 
!###########################################################
!> @brief
!! Interpolates CRP6D smooth 2dcuts in R and Z.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INTERPOL_CRP6D(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(INOUT)::this
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   ! Run section
   CALL this%EXTRACT_VACUUMSURF()   
   CALL this%SMOOTH()
   this%is_smooth=.TRUE.
   DO i = 1, this%nsites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         CALL this%wyckoffsite(i)%zrcut(j)%INTERPOL
      END DO
   END DO
   RETURN
END SUBROUTINE INTERPOL_CRP6D
!###########################################################
!# SUBROUTINE: RAWINTERPOL_CRP6D 
!###########################################################
!> @brief
!! Interpolates CRP6D 2dcuts in R and Z. There is not smoothing
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE RAWINTERPOL_CRP6D(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(INOUT)::this
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   ! Run section
   DO i = 1, this%nsites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         CALL this%wyckoffsite(i)%zrcut(j)%INTERPOL
      END DO
   END DO
   RETURN
END SUBROUTINE RAWINTERPOL_CRP6D
!###########################################################
!# SUBROUTINE: FROM_MOLECULAR_TO_ATOMIC 
!###########################################################
!> @brief
!! Go from molecular coordinates x,y,z,r,theta,phi to
!! xa,ya,za,xb,yb,zb
!-----------------------------------------------------------
SUBROUTINE FROM_MOLECULAR_TO_ATOMIC(ma,mb,molcoord,atomcoord)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   REAL(KIND=8),INTENT(IN) :: ma,mb
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: molcoord
   REAL(KIND=8),DIMENSION(6),INTENT(OUT) :: atomcoord
   ! Run section
   atomcoord(1)=molcoord(1)+(mb/(ma+mb))*molcoord(4)*dcos(molcoord(6))*dsin(molcoord(5))
   atomcoord(2)=molcoord(2)+(mb/(ma+mb))*molcoord(4)*dsin(molcoord(6))*dsin(molcoord(5))
   atomcoord(3)=molcoord(3)+(mb/(ma+mb))*molcoord(4)*dcos(molcoord(5))
   atomcoord(4)=molcoord(1)-(ma/(ma+mb))*molcoord(4)*dcos(molcoord(6))*dsin(molcoord(5))
   atomcoord(5)=molcoord(2)-(ma/(ma+mb))*molcoord(4)*dsin(molcoord(6))*dsin(molcoord(5))
   atomcoord(6)=molcoord(3)-(ma/(ma+mb))*molcoord(4)*dcos(molcoord(5))
   RETURN
END SUBROUTINE FROM_MOLECULAR_TO_ATOMIC
!###########################################################
!# SUBROUTINE: SMOOTH_CRP6D
!###########################################################
!> @brief
!! Smooths a an Rz-2dcut of the potential
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 25/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SMOOTH_CRP6D(this)
   ! Initial declarations   
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   ! I/O variables
    CLASS(CRP6D),INTENT(INOUT) :: this
   ! Local variables
   REAL(KIND=8),DIMENSION(6) :: molcoord,atomcoord
   INTEGER(KIND=4) :: nr,nz
   INTEGER(KIND=4) :: i,j,k,l ! counters
   REAL(KIND=8) :: ma,mb
   CHARACTER(LEN=14),PARAMETER :: routinename="SMOOTH_CRP6D: "
   REAL(KIND=8) :: newpot
   ! Run section
   ma=this%atomdat(1)%getmass()
   mb=this%atomdat(2)%getmass()
   DO i = 1, this%nsites ! cycle wyckoff sites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         molcoord(1)=this%wyckoffsite(i)%zrcut(j)%x
         molcoord(2)=this%wyckoffsite(i)%zrcut(j)%y
         molcoord(5)=this%wyckoffsite(i)%zrcut(j)%theta
         molcoord(6)=this%wyckoffsite(i)%zrcut(j)%phi
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizer()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizez()
         DO k = 1, nr
            DO l = 1, nz
               molcoord(3)=this%wyckoffsite(i)%zrcut(j)%getgridvalueZ(l)
               molcoord(4)=this%wyckoffsite(i)%zrcut(j)%getgridvalueR(k)
               CALL FROM_MOLECULAR_TO_ATOMIC(ma,mb,molcoord,atomcoord)
#ifdef DEBUG
               CALL DEBUG_WRITE(routinename,"Molecular coords:")
               CALL DEBUG_WRITE(routinename,molcoord)
               CALL DEBUG_WRITE(routinename,"Atomic coords:")
               CALL DEBUG_WRITE(routinename,atomcoord)
#endif
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)
               SELECT CASE(this%natomic)
                  CASE(1)
                     newpot=newpot-this%atomiccrp(1)%getpot(atomcoord(1:3))-&
                        this%atomiccrp(1)%getpot(atomcoord(4:6))
                     CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  CASE(2)
                     newpot=newpot-this%atomiccrp(1)%getpot(atomcoord(1:3))-&
                        this%atomiccrp(2)%getpot(atomcoord(4:6))
                     CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  CASE DEFAULT
                     WRITE(0,*) "SMOOTH_CRP6D ERR: Something is wrong with number of atomic crp potentials"
                     CALL EXIT(1)
               END SELECT
            END DO
         END DO
      END DO
   END DO
   RETURN
END SUBROUTINE SMOOTH_CRP6D
!###########################################################
!# SUBROUTINE: EXTRACT_VACUUMSURF_CRP6D
!###########################################################
!> @brief
!! Extracts energy at equilibrium distance in the vacuum and surface energy
!! to an Rz-2dcut of the potential. Shifts the vacuum potential as well.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 25/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE EXTRACT_VACUUMSURF_CRP6D(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
    CLASS(CRP6D),INTENT(INOUT) :: this
   ! Local variables
   REAL(KIND=8),DIMENSION(6) :: molcoord,atomcoord
   INTEGER(KIND=4) :: nr,nz
   INTEGER(KIND=4) :: i,j,k,l ! counters
   REAL(KIND=8) :: ma,mb
   REAL(KIND=8) :: newpot
   ! Run section
   DO i = 1, this%nsites ! cycle wyckoff sites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizer()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizez()
         DO k = 1, nr
            DO l = 1, nz
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)
               newpot=newpot-this%farpot%getscalefactor()
               CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
            END DO
         END DO
      END DO
   END DO
   CALL this%farpot%SHIFTPOT()
   RETURN
END SUBROUTINE EXTRACT_VACUUMSURF_CRP6D
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
   USE UNITS_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(OUT) :: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   INTEGER(KIND=4) :: natomic ! number of atomic potentials
   INTEGER(KIND=4) :: auxint 
   REAL(KIND=8) :: auxr
   CHARACTER(LEN=30) :: cut2dfilename, string
   CHARACTER(LEN=12),PARAMETER :: routinename="READ_CRP6D: "
   CHARACTER,DIMENSION(:),ALLOCATABLE :: letter
   REAL(KIND=8),DIMENSION(2) :: masslist
   CHARACTER(LEN=2),DIMENSION(2) :: symbollist
   CHARACTER(LEN=10) :: units
   TYPE(Mass) :: masss
   ! Run section -----------------------
   CALL this%SET_ALIAS("CRP6D PES")
   CALL this%SET_DIMENSIONS(6)
   ! set up molecular crp
   OPEN (UNIT=180,FILE=filename,STATUS="old",ACTION="read")
   READ(180,*) ! dummy line
   READ(180,*) symbollist(1),auxr,units
   CALL masss%READ(auxr,units)
   CALL masss%TO_STD()
   masslist(1)=masss%getvalue()
   READ(180,*) symbollist(2),auxr,units
   CALL masss%READ(auxr,units)
   CALL masss%TO_STD()
   masslist(2)=masss%getvalue()
   CALL this%SET_ATOMS(2,symbollist,masslist)
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
      CALL this%atomiccrp(i)%INTERPOL()
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
   READ(180,*) string
   CALL this%farpot%INITIALIZE(string)
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
