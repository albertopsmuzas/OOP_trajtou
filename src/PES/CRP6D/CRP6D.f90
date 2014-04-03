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
USE FOURIER_P4MM_MOD
USE WYCKOFF_P4MM_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////
! TYPE: CRP6D
!
!-------------------------------------------------
TYPE,EXTENDS(PES) :: CRP6D
   INTEGER(KIND=4) :: nsites
   INTEGER(KIND=4) :: natomic
   LOGICAL :: is_homonucl=.FALSE.
   LOGICAL :: is_smooth=.FALSE.
   CLASS(Wyckoffsitio),DIMENSION(:),ALLOCATABLE :: wyckoffsite
   TYPE(CRP3D),DIMENSION(:),ALLOCATABLE :: atomiccrp
   TYPE(Vacuumpot) :: farpot
   INTEGER(KIND=4),DIMENSION(:,:),ALLOCATABLE :: xyklist
   CONTAINS
      ! Initialization block
      PROCEDURE,PUBLIC :: READ => READ_CRP6D
      ! Get block
      PROCEDURE,PUBLIC :: GET_V_AND_DERIVS => GET_V_AND_DERIVS_CRP6D
      ! Tools block
      PROCEDURE,PUBLIC :: SMOOTH => SMOOTH_CRP6D
      PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_CRP6D
      PROCEDURE,PUBLIC :: RAWINTERPOL => RAWINTERPOL_CRP6D
      PROCEDURE,PUBLIC :: EXTRACT_VACUUMSURF => EXTRACT_VACUUMSURF_CRP6D
      PROCEDURE,PUBLIC :: INTERPOL_NEW_RZGRID => INTERPOL_NEW_RZGRID_CRP6D
END TYPE CRP6D
CONTAINS
!###########################################################
!# SUBROUTINE: INTERPOL_NEW_RZGRID_CRP6D 
!###########################################################
!> @brief
!! Creates a different R-Z grid using the ones stored in files
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INTERPOL_NEW_RZGRID_CRP6D(this,nRpoints,nZpoints)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(INOUT) :: this
   INTEGER(KIND=4),INTENT(IN) :: nrpoints,nzpoints
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   ! Run section
   DO i = 1, this%nsites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         CALL this%wyckoffsite(i)%zrcut(j)%interrz%INTERPOL_NEWGRID(nrpoints,nzpoints)
      END DO
   END DO
   RETURN
END SUBROUTINE INTERPOL_NEW_RZGRID_CRP6D
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_CRP6D 
!###########################################################
!> @brief
!! Gets the potential and the derivatives respect to all DOFs
!
!> @warning
!! - Assumed a Fourier interpolation in XY (symmetry adapted)
!! - Assumed a Fourier interpolation in theta (symmetry adapted)
!! - Assumed a Fourier interpolation in phi (symmetry adapted)
!! - Assumed a bicubic splines interpolation in R-Z (symmetry independent)
!! - Assumed that Z-R grids have the same size
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_CRP6D(this,x,v,dvdu)
   ! Initial declarations   
   USE DEBUG_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(IN):: this
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: x
   REAL(KIND=8),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(6),INTENT(OUT) :: dvdu
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   REAL(KIND=8) :: ma,mb ! masses
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: f ! smooth function and derivs
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: dfdu ! smooth derivatives
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: xy
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: aux1
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: aux2
   REAL(KIND=8),DIMENSION(6) :: atomicx
   REAL(KIND=8) :: va,vb
   REAL(KIND=8),DIMENSION(3) :: dvdu_atomicA
   REAL(KIND=8),DIMENSION(3) :: dvdu_atomicB
   TYPE(Fourierp4mm) :: xyinterpol
   CHARACTER(LEN=24),PARAMETER :: routinename="GET_V_AND_DERIVS_CRP6D: "
   ! Run section
   SELECT CASE(this%is_smooth)
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         WRITE(0,*) "GET_V_AND_DERIVS_CRP6D ERR: smooth the PES first (CALL thispes%SMOOTH()"
         CALl EXIT(1)
   END SELECT
   SELECT CASE(x(3) >= this%wyckoffsite(1)%zrcut(1)%getlastz()) 
      CASE(.TRUE.)
         v=this%farpot%getpot(x(4))
         dvdu(1)=0.D0
         dvdu(2)=0.D0
         dvdu(3)=0.D0
         dvdu(4)=this%farpot%getderiv(x(4))
         dvdu(5)=0.D0
         dvdu(6)=0.D0
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   ALLOCATE(f(5,this%nsites))
   ALLOCATE(xy(this%nsites,2))
   DO i = 1, this%nsites 
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"Wyckoffsite: ",i)
      CALL VERBOSE_WRITE(routinename,"Letter:",this%wyckoffsite(i)%id)
#endif
      CALL this%wyckoffsite(i)%GET_V_AND_DERIVS(x(3:6),f(1,i),f(2:5,i))
      xy(i,1)=this%wyckoffsite(i)%x
      xy(i,2)=this%wyckoffsite(i)%y
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,"At XY: ",(/xy(i,1),xy(i,2)/))
      CALL VERBOSE_WRITE(routinename,"With Z/R/Theta/Phi:")
      CALL VERBOSE_WRITE(routinename,x(3:6))
      CALL VERBOSE_WRITE(routinename,"Pot: ",f(1,i))
      CALL VERBOSE_WRITE(routinename,"dvdz: ",f(2,i))
      CALL VERBOSE_WRITE(routinename,"dvdr: ",f(3,i))
      CALL VERBOSE_WRITE(routinename,"dvdtheta: ",f(4,i))
      CALL VERBOSE_WRITE(routinename,"dvdphi: ",f(5,i))
#endif
   END DO
   ! f(1,:) smooth potential values
   ! f(2,:) smooth dvdz
   ! f(3,:) smooth dvdr
   ! f(4,:) smooth dvdtheta
   ! f(5,:) smooth dvdphi
   WRITE(*,*) "xy "
   DO i = 1, this%nsites
      WRITE(*,*) xy(i,:)
   END DO
   WRITE(*,*) "xyklist"
   DO i = 1, this%nsites
      WRITE(*,*) this%xyklist(i,:)
   END DO
   WRITE(*,*) "f"
   DO i = 1, 5
      WRITE(*,*) f(i,:)
   END DO
   CALL xyinterpol%READ(xy,f,this%xyklist)
   CALL xyinterpol%INTERPOL(this%atomiccrp(1)%surf)
   ALLOCATE(aux1(5))
   ALLOCATE(aux2(5,2))
   CALL xyinterpol%GET_F_AND_DERIVS(this%atomiccrp(1)%surf,x(1:2),aux1,aux2)
   WRITE(*,*) aux1(:)
   DO i = 1, 5
      WRITE(*,*) aux2(i,:)
   END DO
#ifdef DEBUG
   !-------------------------------------
   ! Results for the smooth potential
   !-------------------------------------
   ! v=aux1(1)         ! v
   ! dvdu(1)=aux2(1,1) ! dvdx
   ! dvdu(2)=aux2(1,2) ! dvdy
   ! dvdu(3)=aux1(2)   ! dvdz
   ! dvdu(4)=aux1(3)   ! dvdr
   ! dvdu(5)=aux1(4)   ! dvdtheta
   ! dvdu(6)=aux1(5)   ! dvdphi
   CALL VERBOSE_WRITE(routinename,"Values for smooth potential")
   CALL VERBOSE_WRITE(routinename,"v: ",aux1(1))   ! v
   CALL VERBOSE_WRITE(routinename,"dvdx: ",aux2(1,1)) ! dvdx
   CALL VERBOSE_WRITE(routinename,"dvdy: ",aux2(1,2)) ! dvdy
   CALL VERBOSE_WRITE(routinename,"dvdz: ",aux1(2))   ! dvdz
   CALL VERBOSE_WRITE(routinename,"dvdr: ",aux1(3))   ! dvdr
   CALL VERBOSE_WRITE(routinename,"dvdtheta: ",aux1(4))   ! dvdtheta
   CALL VERBOSE_WRITE(routinename,"dvdphi: ",aux1(5))   ! dvdphi
#endif
   !--------------------------------------
   ! Results for the real potential
   !-------------------------------------
   ma=this%atomdat(1)%getmass()
   mb=this%atomdat(2)%getmass()
   CALL FROM_MOLECULAR_TO_ATOMIC(ma,mb,x,atomicx)
   SELECT CASE(this%natomic)
      CASE(1)
         CALL this%atomiccrp(1)%GET_V_AND_DERIVS(atomicx(1:3),va,dvdu_atomicA)
         CALL this%atomiccrp(1)%GET_V_AND_DERIVS(atomicx(4:6),vb,dvdu_atomicB)
      CASE(2)
         CALL this%atomiccrp(1)%GET_V_AND_DERIVS(atomicx(1:3),va,dvdu_atomicA)
         CALL this%atomiccrp(2)%GET_V_AND_DERIVS(atomicx(4:6),vb,dvdu_atomicB)
      CASE DEFAULT
         WRITE(0,*) "GET_V_AND_DERIVS_CRP6D ERR: wrong number of atomic potentials"
         CALL EXIT(1)
   END SELECT
   v=aux1(1)+va+vb
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Contributions of the atomic potential: ")
   CALL VERBOSE_WRITE(routinename, "Position Atom A: ",atomicx(1:3))
   CALL VERBOSE_WRITE(routinename,"Va: ",va)
   CALL VERBOSE_WRITE(routinename,"dVa/dxa; dVa/dya: ",dvdu_atomicA)
   CALL VERBOSE_WRITE(routinename, "Position Atom B: ",atomicx(4:6))
   CALL VERBOSE_WRITE(routinename,"Vb: ",vb)
   CALL VERBOSE_WRITE(routinename,"dVa/dxa; dVb/dyb: ",dvdu_atomicB)
#endif
   dvdu(1)=aux2(1,1)+dvdu_atomicA(1)+dvdu_atomicB(1)
   dvdu(2)=aux2(1,2)+dvdu_atomicA(2)+dvdu_atomicB(2)
   dvdu(3)=aux1(2)+dvdu_atomicA(3)+dvdu_atomicB(3)
   dvdu(4)=aux1(3)&
      +(mb/(ma+mb))*dvdu_atomicA(1)*dcos(x(6))*dsin(x(5))&
      +(mb/(ma+mb))*dvdu_atomicA(2)*dsin(x(6))*dsin(x(5))&
      +(mb/(ma+mb))*dvdu_atomicA(3)*dcos(x(5))&
      -(ma/(ma+mb))*dvdu_atomicB(1)*dcos(x(6))*dsin(x(5))&
      -(ma/(ma+mb))*dvdu_atomicB(2)*dsin(x(6))*dsin(x(5))&
      -(ma/(ma+mb))*dvdu_atomicB(3)*dcos(x(5))
   dvdu(5)=aux1(4)&
      +(mb/(ma+mb))*x(4)*dvdu_atomicA(1)*dcos(x(6))*dcos(x(5))&
      +(mb/(ma+mb))*x(4)*dvdu_atomicA(2)*dsin(x(6))*dcos(x(5))&
      -(mb/(ma+mb))*x(4)*dvdu_atomicA(3)*dsin(x(5))&
      -(ma/(ma+mb))*x(4)*dvdu_atomicB(1)*dcos(x(6))*dcos(x(5))&
      -(ma/(ma+mb))*x(4)*dvdu_atomicB(2)*dsin(x(6))*dcos(x(5))&
      +(ma/(ma+mb))*x(4)*dvdu_atomicB(3)*dsin(x(5))
   dvdu(6)=aux1(5)&
      -(mb/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicA(1)*dsin(x(6))&
      +(mb/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicA(2)*dcos(x(6))&
      +(ma/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicB(1)*dsin(x(6))&
      -(ma/(ma+mb))*x(4)*dsin(x(5))*dvdu_atomicB(2)*dcos(x(6))
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_CRP6D
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
      CASE(1)
         this%is_homonucl=.TRUE.
      CASE(2)
         this%is_homonucl=.FALSE.
      CASE DEFAULT 
         WRITE(0,*) "READ_CRP6D ERR: Wrong number of atomic potentials. Allowed number: 1 or 2."
         CALL EXIT(1)
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
      CALL this%wyckoffsite(i)%SET_HOMONUCL(this%is_homonucl)
   END DO
   ! Read the final part of the input: kpoints for XY interpolation
   ALLOCATE(this%xyklist(this%nsites,2))
   READ(180,*) ! dummy line
   DO i = 1, this%nsites
      READ(180,*) this%xyklist(i,:)
   END DO
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"xyklist found:")
   DO i = 1, this%nsites
      CALL VERBOSE_WRITE(routinename,this%xyklist(i,:))
   END DO
#endif
   CLOSE(180)
   RETURN
END SUBROUTINE READ_CRP6D
END MODULE CRP6D_MOD
