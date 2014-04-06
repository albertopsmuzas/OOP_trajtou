!#########################################################
! MODULE: INTERPOL_WYCKOFF_GENERIC
!> @brief
!! Provides generic routines to interpol a CRP6D PES 
!##########################################################
MODULE WYCKOFF_GENERIC_MOD
   USE BICSPLINES_MOD
   USE FOURIER1D_MOD
! Initial declarations
IMPLICIT NONE
!//////////////////////////////////////////////////
! TYPE: Fouklist
!> @brief
!! Auxiliar type to get rid with some technical problems
!--------------------------------------------------
TYPE Fouklist
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE :: k
END TYPE
!//////////////////////////////////////////////////
! TYPE: Cut2d
!> @brief
!! Contains data for a Z,R interpolation letting X,Y,Theta,Phi
!! fixed
!-------------------------------------------------
TYPE Cut2d
   PRIVATE
   CHARACTER(LEN=30) :: alias
   CHARACTER(LEN=30) :: filename
   REAL(KIND=8),PUBLIC :: x,y,phi,theta
   TYPE(Bicsplines),PUBLIC :: interrz
   CONTAINS
      ! Initialize block
      PROCEDURE,PUBLIC :: READ => READ_CUT2D
      ! Get block
      PROCEDURE,PUBLIC :: getgridsizer => getgridsizer_cut2d
      PROCEDURE,PUBLIC :: getgridsizez => getgridsizez_cut2d
      PROCEDURE,PUBLIC :: getfirstr => getfirstr_cut2d
      PROCEDURE,PUBLIC :: getfirstz => getfirstz_cut2d
      PROCEDURE,PUBLIC :: getlastz => getlastz_cut2d
      PROCEDURE,PUBLIC :: getlastr => getlastr_cut2d
      PROCEDURE,PUBLIC :: getgridvaluer => getgridvaluer_cut2d
      PROCEDURE,PUBLIC :: getgridvaluez => getgridvaluez_cut2d
      PROCEDURE,PUBLIC :: getalias => getalias_cut2d
      PROCEDURE,PUBLIC :: getpotatgridpoint => getpotatgridpoint_cut2d
      ! Tools block
      PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_CUT2D
      PROCEDURE,PUBLIC :: CHANGEPOT_AT_GRIDPOINT => CHANGEPOT_AT_GRIDPOINT_CUT2D
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
   INTEGER(KIND=4) :: mynumber
   LOGICAL :: is_homonucl=.FALSE.
   REAL(KIND=8) :: x,y
   INTEGER(KIND=4) :: n2dcuts
   INTEGER(KIND=4) :: nphicuts
   REAL(KIND=8) ::rinit,zinit
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE :: nphipoints
   TYPE(Fouklist),DIMENSION(:),ALLOCATABLE :: klistphi
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE :: klisttheta
   TYPE(Cut2d),DIMENSION(:),ALLOCATABLE :: zrcut
   CONTAINS
      ! Initialization block
      PROCEDURE,PUBLIC :: READ => READ_WYCKOFFSITIO
      PROCEDURE,PUBLIC :: SET_ID => SET_ID_WYCKOFFSITIO
      PROCEDURE,PUBLIC :: SET_HOMONUCL => SET_HOMONUCL_WYCKOFFSITIO
      PROCEDURE,PUBLIC :: SET_MYNUMBER => SET_MYNUMBER_WYCKOFFSITIO
      ! Interface procedures
      PROCEDURE,PUBLIC :: GET_V_AND_DERIVS => GET_V_AND_DERIVS_WYCKOFFSITIO
END TYPE Wyckoffsitio
!/////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: SET_MYNUMBER_WYCKOFFSITIO 
!###########################################################
!> @brief
!! Common set function. Sets_mynumber atribute
!-----------------------------------------------------------
SUBROUTINE SET_MYNUMBER_WYCKOFFSITIO(this,mynumber)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Wyckoffsitio),INTENT(INOUT):: this
   INTEGER(KIND=4),INTENT(IN) :: mynumber
   ! Run section
   this%mynumber=mynumber
   RETURN
END SUBROUTINE SET_MYNUMBER_WYCKOFFSITIO
!###########################################################
!# SUBROUTINE: SET_HOMONUCL_WYCKOFFSITIO 
!###########################################################
!> @brief
!! Common set function. Sets homonucl atribute
!-----------------------------------------------------------
SUBROUTINE SET_HOMONUCL_WYCKOFFSITIO(this,is_homonucl)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Wyckoffsitio),INTENT(INOUT):: this
   LOGICAL,INTENT(IN) :: is_homonucl
   ! Run section
   this%is_homonucl=is_homonucl
   RETURN
END SUBROUTINE SET_HOMONUCL_WYCKOFFSITIO
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_WYCKOFFSITIO 
!###########################################################
!> @brief
!! Dummy subroutine, acting as an interface 
!-----------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_WYCKOFFSITIO(this,x,v,dvdu)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Wyckoffsitio),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(4),INTENT(IN) ::x 
   REAL(KIND=8),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(4),INTENT(OUT) :: dvdu
   ! Run section
   WRITE(0,*) "GET_V_AND_DERIVS_WYCKOFFSITIO ERR: if I was invoked, something went wrong with type allocation of wyckoffsitio"
   CALL EXIT(1)
END SUBROUTINE GET_V_AND_DERIVS_WYCKOFFSITIO
!###########################################################
!# FUNCTION: getpotatgridpoint_cut2d 
!###########################################################
!> @brief 
!! Common get function. Gets the potential at a given grid point
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getpotatgridpoint_cut2d(this,i,j) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN) :: this
   INTEGER(KIND=4),INTENT(IN) :: i,j
   ! Local variables
   ! local_vars
   ! Run section
   getpotatgridpoint_cut2d=this%interrz%fgrid(i,j)
   RETURN
END FUNCTION getpotatgridpoint_cut2d
!###########################################################
!# SUBROUTINE: CHANGEPOT_AT_GRIDPOINT_CUT2D 
!###########################################################
!> @brief
!! Changes the value of the potential at a given point of the
!! grid.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 26/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE CHANGEPOT_AT_GRIDPOINT_CUT2D(this,i,j,newpot)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(INOUT):: this
   INTEGER(KIND=4),INTENT(IN) :: i,j
   REAL(KIND=8),INTENT(IN) :: newpot
   ! Run section
   this%interrz%fgrid(i,j)=newpot
   RETURN
END SUBROUTINE CHANGEPOT_AT_GRIDPOINT_CUT2D
!###########################################################
!# FUNCTION: getgridvaluer_cut2d 
!###########################################################
!> @brief 
!! Common get function. Gets ith component of the grid in R
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getgridvaluer_cut2d(this,i) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   INTEGER(KIND=4),INTENT(IN) :: i
   ! Run section
   getgridvaluer_cut2d=this%interrz%x(i)
   RETURN
END FUNCTION getgridvaluer_cut2d
!###########################################################
!# FUNCTION: getgridvaluez_cut2d 
!###########################################################
!> @brief 
!! Common get function. Gets ith component of the grid in Z
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getgridvaluez_cut2d(this,i) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   INTEGER(KIND=4),INTENT(IN) :: i
   ! Run section
   getgridvaluez_cut2d=this%interrz%y(i)
   RETURN
END FUNCTION getgridvaluez_cut2d
!###########################################################
!# FUNCTION: getalias_cut2d 
!###########################################################
!> @brief
!! Common get function. Gets alias of a Cut2d object
!-----------------------------------------------------------
CHARACTER(LEN=30) FUNCTION getalias_cut2d(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Run section
   getalias_cut2d=this%alias
   RETURN
END FUNCTION getalias_cut2d
!###########################################################
!# FUNCTION: getgridsizer_cut2d 
!###########################################################
!> @brief
!! Common get function. Gets R grid size
!-----------------------------------------------------------
INTEGER(KIND=4) FUNCTION getgridsizer_cut2d(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Run section
   getgridsizer_cut2d=size(this%interrz%x)
   RETURN
END FUNCTION getgridsizer_cut2d
!###########################################################
!# FUNCTION: getgridsizez_cut2d 
!###########################################################
!> @brief
!! Common get function. Gets Z grid size
!-----------------------------------------------------------
INTEGER(KIND=4) FUNCTION getgridsizez_cut2d(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Run section
   getgridsizez_cut2d=size(this%interrz%y)
   RETURN
END FUNCTION getgridsizez_cut2d
!###########################################################
!# FUNCTION: getfirstr_cut2d 
!###########################################################
!> @brief 
!! Just common get function to get first R value
!
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getfirstr_cut2d(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Local variables
   
   ! Run section
   getfirstr_cut2d=this%interrz%x(1)
   RETURN
END FUNCTION getfirstr_cut2d
!###########################################################
!# FUNCTION: getlastr_cut2d 
!###########################################################
!> @brief 
!! Just common get function to get last R value
!
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getlastr_cut2d(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Run section
   getlastr_cut2d=this%interrz%x(size(this%interrz%x))
   RETURN
END FUNCTION getlastr_cut2d
!###########################################################
!# FUNCTION: getfirstz_cut2d 
!###########################################################
!> @brief 
!! Just common get function to get first Z value
!
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getfirstz_cut2d(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Run section
   getfirstz_cut2d=this%interrz%y(1)
   RETURN
END FUNCTION getfirstz_cut2d
!###########################################################
!# FUNCTION: getlastz_cut2d 
!###########################################################
!> @brief 
!! Just common get function to get last Z value
!
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getlastz_cut2d(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Run section
   getlastz_cut2d=this%interrz%y(size(this%interrz%y))
   RETURN
END FUNCTION getlastz_cut2d
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
   CALL this%interrz%READ(r,z,f)
   CLOSE(10)
   RETURN
END SUBROUTINE READ_CUT2D
!###########################################################
!# SUBROUTINE: INTERPOl_CUT2D 
!###########################################################
!> @brief
!! Interpolates a RZ-2dcur of the potential
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 25/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INTERPOl_CUT2D(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(INOUT) :: this
   ! Local variables
   ! Run section
   CALL this%interrz%INTERPOL()
   RETURN
END SUBROUTINE INTERPOl_CUT2D
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
   ALLOCATE(this%nphipoints(this%nphicuts))
   ALLOCATE(this%klistphi(this%nphicuts))
   READ(u,*) this%nphipoints(:)
   DO i = 1, this%nphicuts
         ALLOCATE(this%klistphi(i)%k(this%nphipoints(i)))
         READ(u,*) this%klistphi(i)%k(:)
   END DO
   ALLOCATE(this%klisttheta(this%nphicuts))
   READ(u,*) this%klisttheta(:)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"n2dcuts: ",this%n2dcuts)
   CALL VERBOSE_WRITE(routinename,"nphicuts: ",this%nphicuts)
   CALL VERBOSE_WRITE(routinename,this%nphipoints(:))
   DO i = 1, this%nphicuts
         CALL VERBOSE_WRITE(routinename,"Kpoins list for phi interpolations:")
         CALL VERBOSE_WRITE(routinename,this%klistphi(i)%k(:))
   END DO
   CALL VERBOSE_WRITE(routinename,"Kpoints for theta interpolation: ")
   CALL VERBOSE_WRITE(routinename,this%klisttheta)
#endif
   DO i = 1, this%n2dcuts
      READ(u,*) filename
      CALL this%zrcut(i)%READ(filename)
   END DO
   ! All zrcuts should belong to the same XY position (center of mass)
   this%x=this%zrcut(1)%x
   this%y=this%zrcut(1)%y
   RETURN
END SUBROUTINE READ_WYCKOFFSITIO

END MODULE WYCKOFF_GENERIC_MOD
