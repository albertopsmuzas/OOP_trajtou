!#########################################################
! MODULE: INTERPOL_WYCKOFF_GENERIC
!> @brief
!! Provides generic routines to interpol a CRP6D PES 
!##########################################################
MODULE WYCKOFF_GENERIC_MOD
   USE BICSPLINES_MOD
   USE FOURIER1D_MOD
   USE MATHS_MOD
   USE UNITS_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
! Initial declarations
IMPLICIT NONE
!//////////////////////////////////////////////////
! TYPE: Fouklist
!> @brief
!! Auxiliar type to get rid of some technical problems
!--------------------------------------------------
TYPE Fouklist
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE:: k
END TYPE
!//////////////////////////////////////////////////
! TYPE: Cut2d
!> @brief
!! Contains data for a Z,R interpolation letting X,Y,Theta,Phi
!! fixed
!-------------------------------------------------
TYPE Cut2d
   PRIVATE
   CHARACTER(LEN=:),ALLOCATABLE:: alias
   CHARACTER(LEN=:),ALLOCATABLE:: filename
   REAL(KIND=8),PUBLIC:: x,y,phi,theta
   TYPE(Bicsplines),PUBLIC:: interrz
   CONTAINS
      ! Initialize block
      PROCEDURE,PUBLIC:: READ => READ_CUT2D
      ! Get block
      PROCEDURE,PUBLIC:: gettheta => gettheta_cut2d
      PROCEDURE,PUBLIC:: getphi => getphi_cut2d
      PROCEDURE,PUBLIC:: getgridsizer => getgridsizer_cut2d
      PROCEDURE,PUBLIC:: getgridsizez => getgridsizez_cut2d
      PROCEDURE,PUBLIC:: getfirstr => getfirstr_cut2d
      PROCEDURE,PUBLIC:: getfirstz => getfirstz_cut2d
      PROCEDURE,PUBLIC:: getlastz => getlastz_cut2d
      PROCEDURE,PUBLIC:: getlastr => getlastr_cut2d
      PROCEDURE,PUBLIC:: getgridvaluer => getgridvaluer_cut2d
      PROCEDURE,PUBLIC:: getgridvaluez => getgridvaluez_cut2d
      PROCEDURE,PUBLIC:: getalias => getalias_cut2d
      PROCEDURE,PUBLIC:: getpotatgridpoint => getpotatgridpoint_cut2d
      ! Tools block
      PROCEDURE,PUBLIC:: PRINT_INPUT => PRINT_INPUT_CUT2D
      PROCEDURE,PUBLIC:: INTERPOL => INTERPOL_CUT2D
      PROCEDURE,PUBLIC:: CHANGEPOT_AT_GRIDPOINT => CHANGEPOT_AT_GRIDPOINT_CUT2D
END TYPE Cut2d
!////////////////////////////////////////////////////////////////////////////////////////////////////
! TYPE: WYCKOFFSITIO
!
!> @brief
!! Stores all data related to a single X,Y position
!! of the molecule on the surface.
!
!> @param id - Letter that specifies the symmetry subgroup
!!             related to this site. Its meaning changes for
!!             each wallpaper group
!> @param mynumber - Number that identifies this wyckoff site
!> @param x,y - Position in XY plane of this wyckoffsite
!> @param is_homonuclear - True if the molecule is homonuclear
!> @param  n2dcuts - Number of Z,R cuts for this site
!> @param nphicuts - Number of Phi interpolations that will exist
!!                   for different Theta values. 
!
!> @brief nphipoints - Array of integer numbers storing the number of points
!!                     in which each phi interpolation is based
!> 
!---------------------------------------------------------------------------------------------------
TYPE,ABSTRACT:: Wyckoffsitio
   !PRIVATE
   CHARACTER:: id
   INTEGER(KIND=4):: mynumber
   LOGICAL:: is_homonucl=.FALSE.
   REAL(KIND=8):: x,y
   INTEGER(KIND=4):: n2dcuts
   INTEGER(KIND=4):: nphicuts
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE:: nphipoints
   TYPE(Cut2d),DIMENSION(:),ALLOCATABLE:: zrcut
   CONTAINS
      ! Initialization block
      PROCEDURE,PUBLIC,NON_OVERRIDABLE:: INITIALIZE => INITIALIZE_WYCKOFFSITIO
      PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_ID => SET_ID_WYCKOFFSITIO
      PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_HOMONUCL => SET_HOMONUCL_WYCKOFFSITIO
      PROCEDURE,PUBLIC,NON_OVERRIDABLE:: SET_MYNUMBER => SET_MYNUMBER_WYCKOFFSITIO
      ! Interface procedures
      PROCEDURE(GET_V_AND_DERIVS_WYCKOFFSITIO),PUBLIC,DEFERRED:: GET_V_AND_DERIVS
END TYPE Wyckoffsitio

ABSTRACT INTERFACE
   !###########################################################
   !# SUBROUTINE: GET_V_AND_DERIVS_WYCKOFFSITIO 
   !###########################################################
   !-----------------------------------------------------------
   SUBROUTINE GET_V_AND_DERIVS_WYCKOFFSITIO(this,x,v,dvdu)
      IMPORT Wyckoffsitio
      CLASS(Wyckoffsitio),INTENT(IN) :: this
      REAL(KIND=8),DIMENSION(4),INTENT(IN) ::x 
      REAL(KIND=8),INTENT(OUT) :: v
      REAL(KIND=8),DIMENSION(4),INTENT(OUT) :: dvdu
   END SUBROUTINE GET_V_AND_DERIVS_WYCKOFFSITIO
END INTERFACE
!//////////////////////////////////////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: PRINT_INPUT_CUT2D 
!###########################################################
!> @brief
!! Prints an input file from data stored in RAM
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE PRINT_INPUT_CUT2D(this,filename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables
   INTEGER(KIND=4) :: nr,nz
   INTEGER(KIND=4) :: i,j ! counters
   ! Run section
   nr=this%getgridsizeR()
   nz=this%getgridsizeZ()
   OPEN (10,FILE=filename,STATUS="replace",ACTION="write")
   WRITE(10,*) "# FILE generated by PRINT_CUT2D_INPUT CRP6D"
   WRITE(10,*) this%alias
   WRITE(10,*) "au au rad"
   WRITE(10,*) this%x,this%y
   WRITE(10,*) this%theta,this%phi
   WRITE(10,*) nr,nz
   DO j = 1, nz
      DO i = 1, nr
         WRITE(10,*) this%getgridvalueR(i),this%getgridvalueZ(j),this%getpotatgridpoint(i,j)
      END DO
   END DO
   CLOSE(10)
   RETURN
END SUBROUTINE PRINT_INPUT_CUT2D
!###########################################################
!# FUNCTION: gettheta_cut2d 
!###########################################################
!> @brief 
!! Common get functions. Gets theta atribute
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION gettheta_cut2d(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Run section
   gettheta_cut2d=this%theta
   RETURN
END FUNCTION gettheta_cut2d
!###########################################################
!# FUNCTION: getphi_cut2d 
!###########################################################
!> @brief 
!! Common get functions. Gets phi atribute
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION getphi_cut2d(this) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(IN):: this
   ! Run section
   getphi_cut2d=this%phi
   RETURN
END FUNCTION getphi_cut2d
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
   CLASS(Wyckoffsitio),INTENT(INOUT):: this
   CHARACTER,INTENT(IN):: id
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
   IMPLICIT NONE
   ! I/O variables
   CLASS(Cut2d),INTENT(OUT):: this
   CHARACTER(LEN=*),INTENT(IN):: filename
   ! Local variables
   INTEGER(KIND=4):: i,j ! counters
   INTEGER(KIND=4):: nx,ny
   REAL(KIND=8):: auxr1,auxr2
   CHARACTER(LEN=10):: units1,units2,units3
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: z,r
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: f
   TYPE(Length):: len1,len2
   TYPE(Energy):: en
   TYPE(Angle):: angl 
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
!# SUBROUTINE: INITIALIZE_WYCKOFFSITIO 
!###########################################################
!> @brief
!! Initialize a Wyckoffsitio from file.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 20/03/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INITIALIZE_WYCKOFFSITIO(this,nphipoints,filenames,letter,is_homonucl,mynumber)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(Wyckoffsitio),INTENT(INOUT):: this
   INTEGER(KIND=4),DIMENSION(:),INTENT(IN):: nphipoints
   CHARACTER(LEN=*),DIMENSION(:),INTENT(IN):: filenames
   CHARACTER,INTENT(IN):: letter
   LOGICAL,INTENT(IN):: is_homonucl
   INTEGER(KIND=4),INTENT(IN):: mynumber
   ! Local variables
   INTEGER(KIND=4):: i ! counters
   ! Parameters 
   CHARACTER(LEN=*),PARAMETER:: routinename="INITIALIZE_WYCKOFFSITIO: "
   ! Run section
   this%mynumber=mynumber
   this%is_homonucl=is_homonucl
   this%id=letter
   this%n2dcuts=sum(nphipoints(:))
   this%nphicuts=size(nphipoints(:))
   ALLOCATE(this%zrcut(this%n2dcuts))
   ALLOCATE(this%nphipoints(this%nphicuts))
   this%nphipoints(:)=nphipoints(:)
   SELECT CASE(this%n2dcuts==size(filenames(:)))
      CASE(.true.)
         ! do nothing
      CASE(.false.)
         WRITE(0,*) "INITIALIZE_WYCKOFFSITIO ERR: mismatch between number of n2dcuts and number of files to open"
         CALL EXIT(1)
   END SELECT
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"n2dcuts: ",this%n2dcuts)
   CALL VERBOSE_WRITE(routinename,"nphicuts: ",this%nphicuts)
   CALL VERBOSE_WRITE(routinename,'theta structure: ',this%nphipoints(:))
   DO i = 1, this%n2dcuts
      CALL VERBOSE_WRITE(routinename,trim(filenames(i)))
   END DO
#endif
   DO i = 1, this%n2dcuts
      CALL this%zrcut(i)%READ(trim(filenames(i)))
   END DO
   ! All zrcuts should belong to the same XY position (center of mass)
   this%x=this%zrcut(1)%x
   this%y=this%zrcut(1)%y
   ! DEBUGGING PART
   RETURN
END SUBROUTINE INITIALIZE_WYCKOFFSITIO

END MODULE WYCKOFF_GENERIC_MOD
