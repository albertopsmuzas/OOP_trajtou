!#########################################################
! MODULE: HLiF001_WS_MOD
!> @brief
!! CRP3D specific implementation for H/LiF001
!##########################################################
module PES_LIF001_MOD
! Initial declarations
use SYSTEM_MOD
use LiF001SURF_MOD
use PES_MOD
use CUBICSPLINES_MOD
use LINK_FOURIER2D_MOD
use LINK_FUNCTION1D_MOD
use MATHS_MOD, only: order_vect,order,symmetrize
implicit none
!///////////////////////////////////////////////////////////////////////////////
! TYPE & SUBTYPES: Symmetric point
!------------------------------------------------------------------------------
type :: Symmpoint
private
   character(len=:),allocatable:: filename
   character(len=:),allocatable:: alias
   integer(kind=4):: n
   real(kind=8):: x
   real(kind=8):: y
   character(len=10):: units_z
   character(len=10):: units_v
   real(kind=8),dimension(:),allocatable:: z
   real(kind=8),dimension(:),allocatable:: v
   real(KIND=8):: dz1
   real(KIND=8):: dz2
   type(Csplines),public:: interz
   contains
      procedure,public:: PLOT_DATA => PLOT_DATA_SYMMPOINT
      procedure,public:: PLOT => PLOT_INTERPOL_SYMMPOINT
end type Symmpoint

type,extends(Symmpoint) :: Pair_pot
   private
   integer(kind=4):: id
   real(kind=8):: vasint
	real(kind=8):: rumpling
   contains
      procedure,public:: read => read_standard_PAIRPOT
      procedure,public:: get_v_and_derivs => get_v_and_derivs_PAIRPOT
end type Pair_pot

type,extends(Symmpoint) :: Sitio
   private
	real(kind=8),dimension(:),allocatable:: dvdx,dvdy,dvdz
   contains
      procedure,public:: read => READ_STANDARD_SITIO
end type Sitio
!////////////////////////////////////////////////////////////////////////////////
! TYPE: PES_LIF001
!------------------------------------------------------------------------------
type,extends(PES) :: PES_LIF001
   integer(kind=4):: max_order
   type(Pair_pot),dimension(:),allocatable:: all_pairpots
   type(Sitio),dimension(:),allocatable:: all_sites
   integer(kind=4),dimension(:,:),allocatable:: klist
   class(Function1d),allocatable:: dampfunc
   contains
      ! Initialization block
      procedure,public:: READ => READ_PES_LIF001
      procedure,public:: INITIALIZE => INITIALIZE_PES_LIF001
      ! Get block 
      procedure,public:: GET_V_AND_DERIVS => GET_V_AND_DERIVS_PES_LIF001
      procedure,public:: GET_V_AND_DERIVS_CORRECTION => GET_V_AND_DERIVS_CORRECTION_PES_LIF001
      procedure,public:: GET_REPUL_CORRECTIONS => GET_REPUL_CORRECTIONS_PES_LIF001
      procedure,public:: getpot => getpot_crp3d
      ! Enquire block
      procedure,public:: is_allowed => is_allowed_PES_LIF001
      ! Tools block
      procedure,public:: EXTRACT_VASINT => EXTRACT_VASINT_PES_LIF001
      procedure,public:: SMOOTH => SMOOTH_PES_LIF001
      procedure,public:: INTERPOL => INTERPOL_Z_PES_LIF001
      ! Plot tools
      procedure,public:: PLOT_XYMAP => PLOT_XYMAP_PES_LIF001
      procedure,public:: PLOT_DIRECTION1D => PLOT_DIRECTION1D_PES_LIF001
      procedure,public:: PLOT_SITIOS => PLOT_SITIOS_PES_LIF001
      procedure,public:: PLOT_PAIRPOTS => PLOT_PAIRPOTS_PES_LIF001
      procedure,public:: PLOT_Z => PLOT_Z_PES_LIF001
end type PES_LIF001
!///////////////////////////////////////////////////////////////////////////
contains
!###########################################################
!# SUBROUTINE: INITIALIZE_PES_LIF001
!###########################################################
SUBROUTINE INITIALIZE_PES_LIF001(this,filename,tablename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_LIF001),INTENT(OUT):: this
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: filename,tablename
   ! Local variables
   CHARACTER(LEN=:),ALLOCATABLE:: auxstring
   ! Run section
   SELECT CASE(allocated(system_inputfile) .or. .not.present(filename))
      CASE(.TRUE.)
         auxstring=trim(system_inputfile)
      CASE(.FALSE.)
         auxstring=trim(filename)
   END SELECT
   SELECT CASE(present(tablename))
      CASE(.TRUE.) ! present tablename
         CALL this%READ(filename=trim(auxstring),tablename=trim(tablename))
      CASE(.FALSE.) ! not present tablename
         CALL this%READ(filename=trim(auxstring),tablename='pes')
   END SELECT
   CALL this%INTERPOL()
   RETURN
END SUBROUTINE INITIALIZE_PES_LIF001
!###########################################################
!# SUBROUTINE: GET_REPUL_CORRECTIONS_PES_LIF001
!###########################################################
SUBROUTINE GET_REPUL_CORRECTIONS_PES_LIF001(this,P,v,dvdz,dvdx,dvdy)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_LIF001),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(3),INTENT(IN) :: P
   REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT):: dvdx,dvdy,dvdz ! corrections to the derivatives
   ! Local variables
   INTEGER(KIND=4) :: npairpots
   INTEGER(KIND=4) :: l,k ! counters
   REAL(KIND=8) :: aux1,aux2,aux3,aux4
   ! Run section
   npairpots=size(this%all_pairpots)
   FORALL(l=1:npairpots)
      v(l)=0.D0
      dvdz(l)=0.D0
      dvdx(l)=0.D0
      dvdy(l)=0.D0
   END FORALL
   DO l = 1, npairpots
      DO k = 0, this%max_order
         CALL INTERACTION_AENV(k,P,this%all_pairpots(l),this%dampfunc,aux1,aux2,aux3,aux4)
         v(l)=v(l)+aux1
         dvdz(l)=dvdz(l)+aux2
         dvdx(l)=dvdx(l)+aux3
         dvdy(l)=dvdy(l)+aux4
      END DO
   END DO
   RETURN
END SUBROUTINE GET_REPUL_CORRECTIONS_PES_LIF001
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_PAIRPOT 
!###########################################################
SUBROUTINE GET_V_AND_DERIVS_PAIRPOT(this,x,v,dvdu)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Pair_pot),INTENT(IN):: this
   REAL(KIND=8),INTENT(IN) :: x
   REAL(KIND=8),INTENT(OUT) :: v,dvdu
   ! Run section
   SELECT CASE(X>this%z(this%n)-this%rumpling)
      CASE(.TRUE.)
         v=0.D0
         dvdu=0.D0
      CASE(.FALSE.)
         CALL this%interz%GET_V_AND_DERIVS(x,v,dvdu,this%rumpling)
   END SELECT
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_PAIRPOT
!###########################################################
!# SUBROUTINE: READ_STANDARD_PAIRPOT #######################
!###########################################################
SUBROUTINE READ_STANDARD_PAIRPOT(pairpot,filename)
      ! Initial declarations
   IMPLICIT NONE
   ! I/O variables ------------------------------
   CLASS(Pair_pot),INTENT(INOUT) :: pairpot
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables ----------------------------
   INTEGER :: i
   CHARACTER(LEN=14), PARAMETER :: routinename = "READ_PAIRPOT: "
   ! Run section ---------------------------------
   pairpot%filename=filename
   OPEN(10,FILE=pairpot%filename,STATUS='OLD')
   READ(10,*) ! dummy line
   READ(10,*) ! dummy line
   READ(10,*) pairpot%alias
   READ(10,*) pairpot%vasint
   READ(10,*) pairpot%dz1
   READ(10,*) pairpot%dz2
   READ(10,*) pairpot%id,pairpot%rumpling
   READ(10,*) pairpot%n
   ALLOCATE(pairpot%z(1:pairpot%n))
   ALLOCATE(pairpot%v(1:pairpot%n))
   DO i=1, pairpot%n
      READ(10,*) pairpot%z(i), pairpot%v(i)
   END DO
   CLOSE(10)
   CALL pairpot%interz%READ(pairpot%z,pairpot%v)
#ifdef DEBUG
      CALL VERBOSE_WRITE(routinename,pairpot%filename,pairpot%alias)
#endif
   RETURN
END SUBROUTINE READ_STANDARD_PAIRPOT
!###########################################################
!# SUBROUTINE: READ_STANDARD_SITIO #########################
!###########################################################
SUBROUTINE READ_STANDARD_SITIO(site,filename)
   IMPLICIT NONE
   ! I/O variables
   CLASS(Sitio),INTENT(INOUT):: site
   CHARACTER(LEN=*),INTENT(IN):: filename
   ! Local variables
   INTEGER:: i ! counter
   CHARACTER(LEN=*),PARAMETER:: routinename = "READ_SITIO: "
   !
   site%filename=filename
   OPEN (10,file=site%filename,status="old")
   READ(10,*) ! dummy line
   READ(10,*) ! dummy line
   READ(10,*) site%alias
   READ(10,*) site%x, site%y
   READ(10,*) site%n
   READ(10,*) site%dz1
   READ(10,*) site%dz2
   ALLOCATE(site%z(1:site%n))
   ALLOCATE(site%v(1:site%n))
   DO i=1, site%n
      READ(10,*) site%z(i), site%v(i)
   END DO
   CLOSE(10)
   CALL site%interz%READ(site%z,site%v)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,site%filename,site%alias)
#endif
   RETURN
END SUBROUTINE READ_STANDARD_SITIO
!###########################################################
!# SUBROUTINE: READ_PES_LIF001
!###########################################################
SUBROUTINE READ_PES_LIF001(this,filename,tablename)
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_LIF001),INTENT(OUT) :: this
   CHARACTER(LEN=*),INTENT(IN) :: filename
   CHARACTER(LEN=*),INTENT(IN):: tablename
   ! Local variables
   INTEGER(KIND=4) :: n_pairpots,n_sites
   CHARACTER(LEN=1024),DIMENSION(:),ALLOCATABLE:: files_pairpots,files_sites
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: param
   INTEGER(KIND=4):: ierr
   INTEGER(KIND=4) :: i ! counter
   ! Lua-related variables
   TYPE(flu_State):: conf ! Lua state
   INTEGER(KIND=4):: pes_table,pairpot_table,sitio_table,dampfunc_table,param_table,fourier_table ! tables
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE:: subtables
   ! Auxiliar (dummy) variables
   INTEGER(KIND=4):: auxint
   CHARACTER(LEN=1024):: auxstring
   ! Parameters
   CHARACTER(LEN=*),PARAMETER :: routinename="READ_PES_LIF001: "
   ! HEY HO!, LET'S GO!! ------------------
   ! Open Lua file
   CALL OPEN_CONFIG_FILE(L=conf,filename=filename,ErrCode=ierr)
   SELECT CASE(ierr)
      CASE(0)
         ! do nothing
      CASE DEFAULT
         WRITE(0,*) "READ_PES_LIF001 ERR: error reading Lua config file: ",filename
         CALL EXIT(1)
   END SELECT
   ! Open PES table
   CALL AOT_TABLE_OPEN(L=conf,thandle=pes_table,key=tablename)
   ! Set pestype (kind)
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='kind',val=auxstring)
   CALL this%SET_PESTYPE(trim(auxstring))
   SELECT CASE(trim(auxstring))
      CASE('PES_LIF001')
         ! do nothing
      CASE DEFAULT
         WRITE(0,*) "READ_PES_LIF001 ERR: wrong type of PES. Expected: PES_LIF001. Encountered: "//trim(auxstring)
         CALL EXIT(1)
   END SELECT
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'Type of PES: '//trim(auxstring))
#endif
   ! Set alias (name)
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='name',val=auxstring)
   CALL this%SET_ALIAS(trim(auxstring))
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'PES Name: '//trim(auxstring))
#endif
   ! Set dimensions
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='dimensions',val=auxint)
   CALL this%SET_DIMENSIONS(auxint)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'PES dimensions: ',auxint)
#endif
   ! Set max environment
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='maxEnvironment',val=auxint)
   this%max_order=auxint
   ! Set pair potentials
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=pairpot_table,key='pairPotentials')
   n_pairpots=aot_table_length(L=conf,thandle=pairpot_table)
   ALLOCATE(files_pairpots(n_pairpots))
   ALLOCATE(this%all_pairpots(n_pairpots))
   DO i = 1, n_pairpots
      CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pairpot_table,pos=i,val=files_pairpots(i))
      CALL this%all_pairpots(i)%READ(trim(files_pairpots(i)))
   END DO
   CALL AOT_TABLE_CLOSE(L=conf,thandle=pairpot_table)
   ! Set damping function
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=dampfunc_table,key='dampFunction')
   CALL AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=dampfunc_table,key='kind',val=auxstring)
   SELECT CASE(trim(auxstring))
      CASE("Logistic")
         ALLOCATE(Logistic_func::this%dampfunc)
         ALLOCATE(param(2))
         ! open param table
         CALL AOT_TABLE_OPEN(L=conf,parent=dampfunc_table,thandle=param_table,key='param')
         auxint=aot_table_length(L=conf,thandle=param_table)
         SELECT CASE(auxint/=2) 
            CASE(.TRUE.)
               WRITE(0,*) "READ_PES_LIF001 ERR: wrong number of parameters in pes.dampFunc.param table"
               CALL EXIT(1)
            CASE(.FALSE.)
               ! do nothing
         END SELECT
         DO i = 1, 2
            CALL AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=param_table,pos=i,val=param(i))
         END DO
         CALL this%dampfunc%READ(param)
         CALL AOT_TABLE_CLOSE(L=conf,thandle=param_table)
      CASE("None")
         ALLOCATE(One_func::this%dampfunc)
      CASE DEFAULT
         WRITE(0,*) "READ_PES_LIF001 ERR: dampfunction keyword is not implemented"
         WRITE(0,*) "Implemented ones: Logistic, None"
         WRITE(0,*) "Case sensitive"
         CALL EXIT(1)
   END SELECT
   CALL AOT_TABLE_CLOSE(L=conf,thandle=dampfunc_table)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'Damping function used: '//trim(auxstring))
   SELECT CASE(allocated(param))
      CASE(.TRUE.)
         CALL VERBOSE_WRITE(routinename,'Damping function parameters: ',param(:))
      CASE(.FALSE.)
         ! do nothing
   END SELECT
#endif
   ! Set sitios
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=sitio_table,key='sitios')
   n_sites=aot_table_length(L=conf,thandle=sitio_table)
   ALLOCATE(files_sites(n_sites))
   ALLOCATE(this%all_sites(n_sites))
   DO i = 1, n_sites
      CALl AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=sitio_table,pos=i,val=files_sites(i))
      CALL this%all_sites(i)%READ(trim(files_sites(i))) 
   END DO
   CALL AOT_TABLE_CLOSE(L=conf,thandle=sitio_table)
   ! Read fourier Kpoints
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=fourier_table,key='fourierKpoints')
   auxint=aot_table_length(L=conf,thandle=fourier_table)
   SELECT CASE(auxint/=n_sites)
      CASE(.TRUE.)
         WRITE(0,*) "READ_PES_LIF001 ERR: dimension mismatch between fourierKpoints and number of sitios"
         WRITE(0,*) 'Number of Fourier Kpoints: ',auxint
         WRITE(0,*) 'Number of sitios: ',n_sites
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   ALLOCATE(this%klist(n_sites,2))
   ALLOCATE(subtables(n_sites))
   DO i = 1, n_sites
      CALL AOT_TABLE_OPEN(L=conf,parent=fourier_table,thandle=subtables(i),pos=i)
      CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=subtables(i),pos=1,val=auxint)
      this%klist(i,1)=auxint
      CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=subtables(i),pos=2,val=auxint)
      this%klist(i,2)=auxint
      CALL AOT_TABLE_CLOSE(L=conf,thandle=subtables(i))
   END DO
   CALL AOT_TABLE_CLOSE(L=conf,thandle=fourier_table)
   CALL AOT_TABLE_CLOSE(L=conf,thandle=pes_table)
   CALL CLOSE_CONFIG(conf)
   ! VERBOSE PRINT 
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Maximum environmental order: ",this%max_order)
   CALL VERBOSE_WRITE(routinename,"Number of pair potentials: ",n_pairpots)
   CALL VERBOSE_WRITE(routinename,"Pair potentials input files:")
   DO i = 1, n_pairpots
      CALL VERBOSE_WRITE(routinename,trim(files_pairpots(i)))
   END DO
   CALL VERBOSE_WRITE(routinename,"Number of sitios: ",n_sites)
   CALL VERBOSE_WRITE(routinename,"Sitios input files:")
   DO i = 1, n_sites
      CALL VERBOSE_WRITE(routinename,trim(files_sites(i)))
   END DO
   CALL VERBOSE_WRITE(routinename,"List of Kpoints for Fourier interpolation: ")
   DO i = 1, n_sites
      CALL VERBOSE_WRITE(routinename,this%klist(i,:))
   END DO
#endif
   RETURN
END SUBROUTINE READ_PES_LIF001
!#######################################################################
!# SUBROUTINE: EXTRACT_VASINT_PES_LIF001 ######################################
!#######################################################################
SUBROUTINE EXTRACT_VASINT_PES_LIF001(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_LIF001),INTENT(INOUT) :: this
   ! Local variables
   INTEGER(KIND=4) :: npairpots, nsites
   INTEGER(KIND=4) :: i,j ! counters
   REAL(KIND=8) :: control_vasint
   CHARACTER(LEN=22) :: routinename="EXTRACT_VASINT_PES_LIF001: "
   ! Run section ------------------------
   npairpots=size(this%all_pairpots)
   control_vasint=this%all_pairpots(1)%vasint
   DO i = 1, npairpots
      IF (this%all_pairpots(1)%vasint/=control_vasint) THEN
         WRITE(0,*) "EXTRACT_VASINT_PES_LIF001 ERR: Incoherences in vasint values found"
         WRITE(0,*) "EXTRACT_VASINT_PES_LIF001 ERR: vasint's value at pairpot",1,control_vasint
         WRITE(0,*) "EXTRACT_VASINT_PES_LIF001 ERR: vasint's value at pairpot",i,control_vasint
         CALL EXIT(1)
      END IF
      DO j = 1, this%all_pairpots(i)%n
         this%all_pairpots(i)%interz%f(j)=this%all_pairpots(i)%interz%f(j)-this%all_pairpots(i)%vasint
      END DO
#ifdef DEBUG
      CALL DEBUG_WRITE(routinename,"Vasint extracted from pair potential ",i)
#endif
   END DO
   nsites=size(this%all_sites)
   DO i = 1, nsites
      DO j = 1, this%all_sites(i)%n
         this%all_sites(i)%interz%f(j)=this%all_sites(i)%interz%f(j)-this%all_pairpots(1)%vasint
      END DO
#ifdef DEBUG
      CALL DEBUG_WRITE(routinename,"Vasint extracted from pair site ",i)
#endif
   END DO
   RETURN
END SUBROUTINE EXTRACT_VASINT_PES_LIF001
!############################################################
!# SUBROUTINE: SMOOTH_PES_LIF001 ############################
!############################################################
SUBROUTINE SMOOTH_PES_LIF001(this)
   ! Initial declaraitons
   IMPLICIT NONE
   CLASS(PES_LIF001),INTENT(INOUT):: this
   ! Local variables
   REAL(KIND=8),DIMENSION(3):: A
   INTEGER(KIND=4):: i,j ! counters
   INTEGER(KIND=4):: npairpots,nsites
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: v,dvdzr,dummy
   CHARACTER(LEN=*),PARAMETER:: routinename="SMOOTH_PES_LIF001: "
   ! Run section ----------
   nsites = size(this%all_sites)
   npairpots = size(this%all_pairpots)
   ALLOCATE(v(npairpots))
   ALLOCATE(dvdzr(npairpots))
   ALLOCATE(dummy(npairpots))
   DO i = 1, nsites ! loop over sites
      A(1) = this%all_sites(i)%x
      A(2) = this%all_sites(i)%y
      DO j = 1, this%all_sites(i)%n ! loop over pairs v,z
         A(3)=this%all_sites(i)%z(j)
         CALL this%GET_REPUL_CORRECTIONS(A,v,dvdzr,dummy,dummy)
         this%all_sites(i)%interz%f(j)=this%all_sites(i)%interz%f(j)-sum(v)
         IF (j.EQ.1) THEN
            this%all_sites(i)%dz1=this%all_sites(i)%dz1-sum(dvdzr) ! correct first derivative
         ELSE IF (j.EQ.this%all_sites(i)%n) THEN
            this%all_sites(i)%dz2=this%all_sites(i)%dz2-sum(dvdzr) ! correct first derivative
         END IF
      END DO
   END DO
   RETURN
END SUBROUTINE SMOOTH_PES_LIF001
!############################################################
!# SUBROUTINE: INTERPOL_Z_PES_LIF001 ########################
!############################################################
SUBROUTINE INTERPOL_Z_PES_LIF001(this)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_LIF001),INTENT(INOUT) :: this
   ! Local variables
   INTEGER(KIND=4) :: nsites,npairpots
   REAL(KIND=8) :: dz1,dz2
   INTEGER(KIND=4) :: i ! counters
   ! Run secton------------------------
   nsites=size(this%all_sites)
   npairpots=size(this%all_pairpots)
   CALL this%EXTRACT_VASINT()
   DO i = 1, npairpots ! loop pairpots
      dz1=this%all_pairpots(i)%dz1
      dz2=this%all_pairpots(i)%dz2
      CALL this%all_pairpots(i)%interz%INTERPOL(dz1,1,dz2,1)
   END DO
   CALL this%SMOOTH()
   DO i = 1, nsites
      dz1=this%all_sites(i)%dz1
      dz2=this%all_sites(i)%dz2
      CALL this%all_sites(i)%interz%INTERPOL(dz1,1,dz2,1) 
   END DO
   RETURN
END SUBROUTINE INTERPOL_Z_PES_LIF001
!##################################################################
!# SUBROUTINE: INTERACTION_AP #####################################
!##################################################################
SUBROUTINE INTERACTION_AP(A,P,pairpot,dampfunc,interac,dvdz_corr,dvdx_corr,dvdy_corr)
   IMPLICIT NONE
   ! I/O VAriables --------------------------------------------
   REAL(KIND=8),DIMENSION(3),INTENT(IN):: A, P
   TYPE(Pair_pot),INTENT(IN):: pairpot
   CLASS(Function1d),INTENT(IN):: dampfunc
   REAL(KIND=8),INTENT(OUT):: interac,dvdz_corr,dvdx_corr,dvdy_corr
   ! Local variables ------------------------------------------
   REAL(KIND=8):: r ! distance
   REAL(KIND=8):: v,pre
   REAL(KIND=8):: aux ! dv/dr
   CHARACTER(LEN=*),PARAMETER:: routinename = "INTERACTION_AP: "
   ! GABBA, GABBA HEY! ----------------------------------------
   ! Find the distance between A and P, in a.u.
   r=dsqrt((A(1)-P(1))**2.D0+(A(2)-P(2))**2.D0+(A(3)-P(3))**2.D0)
   CALL pairpot%GET_V_AND_DERIVS(r,v,aux)
   interac=v*dampfunc%getvalue(r)
   pre=aux*dampfunc%getvalue(r)+v*dampfunc%getderiv(r)
   dvdz_corr=pre*(A(3)-P(3))/r
   dvdx_corr=pre*(A(1)-P(1))/r
   dvdy_corr=pre*(A(2)-P(2))/r
   RETURN
END SUBROUTINE INTERACTION_AP
!##################################################################
!# SUBROUTINE: INTERACTION_AENV ###################################
!##################################################################
SUBROUTINE INTERACTION_AENV(n,A,pairpot,dampfunc,interac,dvdz_term,dvdx_term,dvdy_term)
   IMPLICIT NONE
   ! I/O VAriables ------------------------------------------
   INTEGER,INTENT(IN) :: n
   REAL(KIND=8),DIMENSION(3),INTENT(IN) :: A
   TYPE(Pair_pot),INTENT(IN) :: pairpot
   CLASS(Function1d),INTENT(IN) :: dampfunc
   REAL(KIND=8),INTENT(OUT) :: interac, dvdz_term, dvdx_term, dvdy_term
   ! Local variables ----------------------------------------
   REAL(KIND=8),DIMENSION(3) :: P
   REAL(KIND=8),DIMENSION(3) :: ghost_A ! A in cartesians, but inside unitcell
   REAL(KIND=8),DIMENSION(2) :: aux
   REAL(KIND=8) :: dummy1, dummy2, dummy3, dummy4
   REAL(KIND=8) :: atomx, atomy
   INTEGER :: pairid
   INTEGER :: i, k ! Counters
   CHARACTER(LEN=18), PARAMETER :: routinename = "INTERACTION_AENV: "
   ! SUSY IS A HEADBANGER !!!! -------------------
   ! Defining some aliases to make the program simpler:
   pairid = pairpot%id
   interac=0.D0
   dvdz_term=0.D0
   dvdx_term=0.D0
   dvdy_term=0.D0
   ! ghost A definition
   ghost_A(1:2)=system_surface%project_unitcell(A(1:2))
   ghost_A(3)=A(3)

   SELECT CASE(n)
      CASE(0)
         DO i=1, system_surface%atomtype(pairid)%n
            P(:)=system_surface%atomtype(pairid)%atom(i,:)
            CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac=interac+dummy1
            dvdz_term=dvdz_term+dummy2
            dvdx_term=dvdx_term+dummy3
            dvdy_term=dvdy_term+dummy4
         END DO
         RETURN

      CASE(1 :)
         DO i=1, system_surface%atomtype(pairid)%n
            atomx=system_surface%atomtype(pairid)%atom(i,1)
            atomy=system_surface%atomtype(pairid)%atom(i,2)
            P(3)=system_surface%atomtype(pairid)%atom(i,3)
            DO k= -n,n
               aux(1)=dfloat(n)
               aux(2)=dfloat(k)
               aux=system_surface%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac=interac+dummy1
               dvdz_term=dvdz_term+dummy2
               dvdx_term=dvdx_term+dummy3
               dvdy_term=dvdy_term+dummy4
               !
               aux(1)=dfloat(-n)
               aux(2)=dfloat(k)
               aux=system_surface%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
            END DO
            DO k= -n+1, n-1
               aux(1)=dfloat(k)
               aux(2)=dfloat(n)
               aux=system_surface%surf2cart(aux)
               P(1) =atomx+aux(1)
               P(2) =atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
               !
               aux(1)=dfloat(k)
               aux(2)=dfloat(-n)
               aux=system_surface%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               CALL INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
            END DO
         END DO
         RETURN

      CASE DEFAULT
         WRITE(0,*) "INTERACTION_AENV ERR: Wrong environment order."
         CALL EXIT(1)
   END SELECT
END SUBROUTINE INTERACTION_AENV
!############################################################
!# SUBROUTINE: GET_V_AND_DERIVS_PES_LIF001 ##################
!############################################################
!> @brief
!! Subroutine that calculates the 3D potential for a point A and
!! its derivatives in cartesian coordinates.
!
!> @param[in] this - PES_LIF001 PES
!> @param[in] X - Point in space to calculate the potential and it's derivatives. Cartesian's
!> @param[out] v - Value of the potential at X
!> @param[out] dvdu - derivatives, cartesian coordinates 
!
!> @warning
!! - All sitios should have been interpolated in Z direction previously. It
!!   is optinal to extract "vasint".
!! - Corrugation, MUST have been extracted from sitios previously (smoothed) 
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 06/Feb/2014
!> @version 1.0
!------------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_PES_LIF001(this,X,v,dvdu,errCode)
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_LIF001),TARGET,INTENT(IN):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN):: X
   REAL(KIND=8),INTENT(OUT):: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT):: dvdu
   integer(kind=1),optional,intent(out):: errCode
   ! Local variables
   INTEGER(KIND=4):: nsites,npairpots
   CLASS(Fourier2d),ALLOCATABLE:: interpolxy
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: pot,dvdz,dvdx,dvdy
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: potarr
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: f,derivarr ! arguments to the xy interpolation
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE:: xy ! arguments to the xy interpolation
   INTEGER :: i ! counters
   ! Pointers
	REAL(KIND=8),POINTER:: zmax
   CHARACTER(LEN=*),PARAMETER:: routinename="GET_V_AND_DERIVS_PES_LIF001: "
   zmax => this%all_sites(1)%z(this%all_sites(1)%n)
   npairpots = size(this%all_pairpots)
   nsites = size(this%all_sites)
   ! GABBA, GABBA HEY! ----------------------
   SELECT CASE(X(3)>zmax)
      CASE(.TRUE.)
         dvdu = (/0.D0,0.D0,0.D0/)
         v = 0.D0
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   ! Initialization section
   ALLOCATE(pot(npairpots))
   ALLOCATE(dvdz(npairpots))
   ALLOCATE(dvdx(npairpots))
   ALLOCATE(dvdy(npairpots))
   CALL this%GET_REPUL_CORRECTIONS(X,pot,dvdz,dvdx,dvdy)
   v=sum(pot)
   dvdu(1)=sum(dvdx)
   dvdu(2)=sum(dvdy)
   dvdu(3)=sum(dvdz)
   ! Now, we have all the repulsive interaction and corrections to the derivarives
   ! stored in v(:) and dvdu(:) respectively.
   ! Let's get v and derivatives from xy interpolation of the corrugationless function
   ! f(1,i) ==> v values interpolated for Z at site "i"
   ! f(2,i) ==> dvdz values interpolated for Z at site "i"
   ALLOCATE(f(2,nsites))
   ALLOCATE(xy(nsites,2))
   ! potarr(1) ==> v interpolated at X,Y for a given Z
   ! potarr(2) ==> dvdz interpolated at X,Y for a given Z
   ALLOCATE(potarr(2))
   ! derivarr(1,1:2) ==> dvdx, dvdy interpolated at X,Y
   ! derivarr(2,1:2) ==> d(dvdz)/dx, d(dvdz)/dy interpolated at X,Y. This is not needed
   ALLOCATE(derivarr(2,2))
   DO i=1,nsites
      xy(i,1)=this%all_sites(i)%x
      xy(i,2)=this%all_sites(i)%y
      CALL this%all_sites(i)%interz%GET_V_AND_DERIVS(X(3),f(1,i),f(2,i))
   END DO
   CALL this%SET_FOURIER_SYMMETRY(interpolxy)
   CALL interpolxy%READ(xy,f,this%klist)
   CALL interpolxy%INTERPOL(system_surface)
   CALL interpolxy%GET_F_AND_DERIVS(system_surface,X,potarr,derivarr)
   ! Corrections from the smoothing procedure
   v = v + potarr(1)
   dvdu(1)=dvdu(1)+derivarr(1,1)
   dvdu(2)=dvdu(2)+derivarr(1,2)
   dvdu(3)=dvdu(3)+potarr(2)
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_PES_LIF001
!############################################################
!# SUBROUTINE: GET_V_AND_DERIVS_CORRECTION_PES_LIF001 ############
!############################################################
!> @brief
!! Subroutine that calculates the correction to the 3D PES for a point A and
!! its derivatives in cartesian coordinates.
!
!> @param[in] this - PES_LIF001 PES
!> @param[in] X - Point in space to calculate the potential and it's derivatives. Cartesian's
!> @param[out] v - Value of the potential at X
!> @param[out] dvdu - derivatives, cartesian coordinates 
!
!> @warning
!! - All sitios should have been interpolated in Z direction previously. It
!!   is optinal to extract "vasint".
!! - Corrugation, MUST have been extracted from sitios previously (smoothed) 
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 06/Feb/2014
!> @version 1.0
!------------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_CORRECTION_PES_LIF001(this,X,v,dvdu)
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_LIF001),TARGET,INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(3), INTENT(IN) :: X
   REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(3),INTENT(OUT) :: dvdu
   ! Local variables
   INTEGER(KIND=4) :: npairpots
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: pot,dvdz,dvdx,dvdy
   INTEGER :: i ! counters
   ! Pointers
	REAL(KIND=8), POINTER :: zmax
   CHARACTER(LEN=24),PARAMETER :: routinename="GET_V_AND_DERIVS_PES_LIF001: "
   zmax => this%all_sites(1)%interz%x(this%all_sites(1)%n)
   npairpots = size(this%all_pairpots)
   SELECT CASE(size(v)==npairpots+1)
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         WRITE(*,*) "GET_V_AND_DERIVS_CORRECTION_PES_LIF001 ERR: wrong number of dimensions array v"
         CALL EXIT(1)
   END SELECT
   ! GABBA, GABBA HEY! ----------------------
   SELECT CASE(X(3)>zmax)
      CASE(.TRUE.)
         dvdu = (/0.D0,0.D0,0.D0/)
         FORALL(i=1:npairpots+1) v(i) = 0.D0
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   ! Initializing variables
   ALLOCATE(pot(npairpots))
   ALLOCATE(dvdz(npairpots))
   ALLOCATE(dvdx(npairpots))
   ALLOCATE(dvdy(npairpots))
   FORALL(i=1:npairpots+1) v(i) = 0.D0
   ! Compute
   CALL this%GET_REPUL_CORRECTIONS(X,pot,dvdz,dvdx,dvdy)
   v(1)=sum(pot)
   FORALL(i=2:npairpots+1) v(i)=pot(i-1)
   dvdu(1)=sum(dvdx)
   dvdu(2)=sum(dvdy)
   dvdu(3)=sum(dvdz)
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_CORRECTION_PES_LIF001
!############################################################
!# FUNCTION: getpot_crp3d ###################################
!############################################################
!> @brief
!! Subroutine that calculates the 3D potential for a point A
!
!> @param[in] this - PES_LIF001 PES
!> @param[in] X - Point in space to calculate the potential. Cartesian's
!
!> @warning
!! - All sitios should have been interpolated in Z direction previously. It
!!   is optinal to extract "vasint".
!! - Corrugation, MUST have been extracted from sitios previously (smoothed) 
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 06/Feb/2014
!> @version 1.0
!------------------------------------------------------------
REAL(KIND=8) FUNCTION getpot_crp3d(this,X)
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_LIF001),TARGET,INTENT(IN) :: this
	REAL(KIND=8),DIMENSION(3), INTENT(IN) :: X
   ! Local variables
   CLASS(Fourier2d),ALLOCATABLE:: interpolxy
   REAL(KIND=8):: v
   INTEGER(KIND=4) :: nsites,npairpots
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: pot,dummy
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: potarr
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: f,derivarr ! arguments to the xy interpolation
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: xy ! arguments to the xy interpolation
   INTEGER :: i ! counters
   ! Pointers
   REAL(KIND=8), POINTER :: zmax
   CHARACTER(LEN=24),PARAMETER :: routinename="GET_V_AND_DERIVS_PES_LIF001: "
   zmax => this%all_sites(1)%interz%x(this%all_sites(1)%n)
   npairpots = size(this%all_pairpots)
   nsites = size(this%all_sites)
   ! GABBA, GABBA HEY! ----------------------
   ALLOCATE(pot(npairpots))
   ALLOCATE(dummy(npairpots))
   SELECT CASE(X(3)>zmax)
      CASE(.TRUE.)
         getpot_crp3d=0.D0
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   !
   CALL this%GET_REPUL_CORRECTIONS(X,pot,dummy,dummy,dummy)
   v=sum(pot)
   ! Now, we have all the repulsive interaction and corrections to the derivarives
   ! stored in v(:) and dvdu(:) respectively.
   ! Let's get v and derivatives from xy interpolation of the corrugationless function
   ALLOCATE(f(2,nsites))
   ALLOCATE(xy(nsites,2))
   ALLOCATE(potarr(2))
   ALLOCATE(derivarr(2,2))
   DO i=1,nsites
      xy(i,1)=this%all_sites(i)%x
      xy(i,2)=this%all_sites(i)%y
      CALL this%all_sites(i)%interz%GET_V_AND_DERIVS(X(3),f(1,i),f(2,i))
   END DO
   CALL this%SET_FOURIER_SYMMETRY(interpolxy)
   CALL interpolxy%READ(xy,f,this%klist)
   CALL interpolxy%INTERPOL(system_surface)
   CALL interpolxy%GET_F_AND_DERIVS(system_surface,X,potarr,derivarr)
   ! Corrections from the smoothing procedure
   getpot_crp3d=v+potarr(1)
   RETURN
END FUNCTION getpot_crp3d
!######################################################################
! SUBROUTINE: PLOT_DATA_SYMMPOINT #####################################
!######################################################################
SUBROUTINE PLOT_DATA_SYMMPOINT(this,filename)
   IMPLICIT NONE
   ! I/O Variables ----------------------
   CLASS(Symmpoint), INTENT(IN) :: this
   CHARACTER(LEN=*), INTENT(IN) :: filename
   ! Local variables --------------------
   INTEGER :: i ! Counter
   ! GABBA GABBA HEY!! ------------------
   OPEN(11, FILE=filename, STATUS="replace")
   DO i=1,this%n
      WRITE(11,*) this%z(i),this%v(i)
   END DO
   CLOSE (11)
   WRITE(*,*) "PLOT_DATA_SYMMPOINT: ",this%alias,filename," file created"
END SUBROUTINE PLOT_DATA_SYMMPOINT
!#######################################################################
! SUBROUTINE: PLOT_XYMAP_PES_LIF001
!#######################################################################
SUBROUTINE PLOT_XYMAP_PES_LIF001(this,filename,init_xyz,nxpoints,nypoints,Lx,Ly)
   IMPLICIT NONE
   CLASS(PES_LIF001),INTENT(IN) :: this
   REAL*8,DIMENSION(3),INTENT(IN) :: init_xyz ! Initial position to start the scan (in a.u.)
   INTEGER,INTENT(IN) :: nxpoints, nypoints ! number of points in XY plane
   CHARACTER(LEN=*),INTENT(IN) :: filename ! filename
   REAL*8,INTENT(IN) :: Lx ! Length of X axis 
   REAL*8,INTENT(IN) :: Ly ! Length of X axis 
   ! Local variables
   REAL*8 :: xmin, ymin, xmax, ymax, z
   REAL*8, DIMENSION(3) :: r, dvdu
   REAL*8 :: xdelta, ydelta
   INTEGER :: xinpoints, nxdelta
   INTEGER :: yinpoints, nydelta
   INTEGER :: i, j ! counters
   REAL*8 :: v ! potential
	! GABBA, GABBA HEY! ---------
   xmin = init_xyz(1)
   ymin = init_xyz(2)
   z = init_xyz(3)
   xmax = init_xyz(1)+Lx
   ymax = init_xyz(2)+Ly
   ! For X, grid parameters
   xinpoints=nxpoints-2
   nxdelta=nxpoints-1
   xdelta=Lx/DFLOAT(nxdelta)
   ! For Y, grid parameters
   yinpoints=nypoints-2
   nydelta=nypoints-1
   ydelta=(ymax-ymin)/DFLOAT(nydelta)
   ! Let's go! 
   ! 1st XY point
   OPEN(11,file=filename,status="replace")
   r(1) = xmin
   r(2) = ymin
   r(3) = z
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v
   END DO
   r(2) = ymax
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v
   ! inpoints in XY
   DO i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      r(3) = z
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v
      DO j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         CALL this%GET_V_AND_DERIVS(r,v,dvdu)
         WRITE(11,*) r(1), r(2), v
      END DO
      r(2) = ymax
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v
   END DO
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   r(3) = z
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1), r(2), v
   END DO
   r(2) = ymax
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1), r(2), v
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_XYMAP_PES_LIF001
!#######################################################################
! SUBROUTINE: PLOT_DIRECTION1D_PES_LIF001 ###################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. To define 
!! the direction, the angle alpha is given. 
!
!> @param[in] this - PES_LIF001 PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] angle - Angle between the surface vector S1 and the direction of the
!!                    cut. It should be given in degrees.
!> @param[in] z - Distance to the surface. All points have the same z.
!> @param[in] L - Length of the graphic
!
!> @warning
!! - The graph starts always at 0,0
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 09/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT_DIRECTION1D_PES_LIF001(this,filename,npoints,angle,z,L)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_LIF001),INTENT(IN) :: this
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*), INTENT(IN) :: filename
   REAL*8, INTENT(IN) :: z, angle
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL*8 :: delta,L,v,s,alpha
   REAL*8 :: xmax, xmin, ymax, ymin 
   REAL*8, DIMENSION(3) :: r, dvdu
   INTEGER :: i ! Counter
   CHARACTER(LEN=24), PARAMETER :: routinename = "PLOT_DIRECTION1D_PES_LIF001: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT_DIRECTION1D_PES_LIF001 ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   ! Change alpha to radians
   alpha = angle * PI / 180.D0
   !
   xmin = 0.D0
   ymin = 0.D0
   xmax = L*DCOS(alpha)
   ymax = L*DSIN(alpha)
   r(3) = z
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=L/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(1) = xmin
   r(2) = ymin
   s = DSQRT(r(1)**2.D0+r(2)**2.D0)
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) s, v , dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   ! cycle for inpoints
   DO i=1, inpoints
      r(1)=xmin+(DFLOAT(i)*delta)*DCOS(alpha)
      r(2)=ymin+(DFLOAT(i)*delta)*DSIN(alpha)
      s = DSQRT(r(1)**2.D0+r(2)**2.D0)
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) s, v , dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   END DO
   ! Final value
   r(1) = xmax
   r(2) = ymax
   s = DSQRT(r(1)**2.D0+r(2)**2.D0)
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) s, v, dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT_DIRECTION1D_PES_LIF001
!#######################################################################
!# SUBROUTINE: PLOT_Z_PES_LIF001 #######################################
!#######################################################################
SUBROUTINE PLOT_Z_PES_LIF001(this,npoints,xyz,L,filename)
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(PES_LIF001),INTENT(IN):: this
   INTEGER,INTENT(IN):: npoints
   CHARACTER(LEN=*),INTENT(IN):: filename
   REAL(KIND=8),DIMENSION(3),INTENT(IN):: xyz
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8):: delta,L
   REAL(KIND=8):: zmax, zmin 
   REAL(KIND=8),DIMENSION(3):: r, dvdu
   REAL(KIND=8):: v
   INTEGER:: i ! Counter
   CHARACTER(LEN=*),PARAMETER:: routinename = "PLOT_DIRECTION1D_PES_LIF001: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT_Z_PES_LIF001 ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   r(1)=xyz(1)
   r(2)=xyz(2)
   zmin=xyz(3)
   zmax=xyz(3)+L
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=L/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(3) = zmin
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)
   ! cycle for inpoints
   DO i=1, inpoints
      r(3)=zmin+(DFLOAT(i)*delta)
      CALL this%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(3),v,dvdu(:)
   END DO
   ! Final value
   r(3) = zmax
   CALL this%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT_Z_PES_LIF001
!###########################################################
!# SUBROUTINE: PLOT_INTERPOL_SYMMPOINT 
!###########################################################
SUBROUTINE PLOT_INTERPOL_SYMMPOINT(this,npoints,filename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Symmpoint),INTENT(IN):: this
   INTEGER(KIND=4),INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Run section
   CALL this%interz%PLOT(npoints,filename)
   RETURN
END SUBROUTINE PLOT_INTERPOL_SYMMPOINT
!###########################################################
!# SUBROUTINE: PLOT_PAIRPOTS_PES_LIF001
!###########################################################
SUBROUTINE PLOT_PAIRPOTS_PES_LIF001(this,npoints)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_LIF001),INTENT(IN)::this
   INTEGER(KIND=4),INTENT(IN) :: npoints
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   INTEGER(KIND=4) :: npairpots
   CHARACTER(LEN=12) :: stringbase
   CHARACTER(LEN=13) :: filename
   ! Run section
   npairpots=size(this%all_pairpots)
   WRITE(stringbase,'(A12)') "-pairpot.dat"
   DO i = 1, npairpots
      WRITE(filename,'(I1,A12)') i,stringbase
      CALL this%all_pairpots(i)%PLOT(npoints,filename)
   END DO
   RETURN
END SUBROUTINE PLOT_PAIRPOTS_PES_LIF001
!###########################################################
!# SUBROUTINE: PLOT_SITIOS_PES_LIF001
!###########################################################
SUBROUTINE PLOT_SITIOS_PES_LIF001(this,npoints)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_LIF001),INTENT(IN)::this
   INTEGER(KIND=4),INTENT(IN) :: npoints
   ! Local variables
   INTEGER(KIND=4) :: i ! counters
   INTEGER(KIND=4) :: nsitios
   CHARACTER(LEN=10) :: stringbase
   CHARACTER(LEN=11) :: filename
   ! Run section
   nsitios=size(this%all_sites)
   WRITE(stringbase,'(A10)') "-sitio.dat"
   DO i = 1, nsitios
      WRITE(filename,'(I1,A10)') i,stringbase
      CALL this%all_sites(i)%PLOT(npoints,filename)
   END DO
   RETURN
END SUBROUTINE PLOT_SITIOS_PES_LIF001
!###########################################################
!# FUNCTION: is_allowed_PES_LIF001
!###########################################################
LOGICAL FUNCTION is_allowed_PES_LIF001(this,x)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(PES_LIF001),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8) :: xmin,xmax
   ! Run section
   xmin=this%all_sites(1)%z(1)
   xmax=this%all_sites(1)%z(this%all_sites(1)%n)
   SELECT CASE(size(x)/=3)
      CASE(.TRUE.)
         WRITE(0,*) "is_allowed_PES_LIF001 ERR: array doesn't have 3 dimensions: 3"
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE( x(3)<xmin )
      CASE(.TRUE.)
         is_allowed_PES_LIF001=.FALSE.
      CASE(.FALSE.)
         is_allowed_PES_LIF001=.TRUE.
   END SELECT
   RETURN
END FUNCTION is_allowed_PES_LIF001
END MODULE PES_LIF001_MOD
