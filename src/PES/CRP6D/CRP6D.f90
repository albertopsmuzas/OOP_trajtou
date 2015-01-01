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
   USE UNITS_MOD
   USE CRP3D_MOD
   USE SYSTEM_MOD
   USE EXTRAPOL_TO_VACUUM_MOD
   USE FOURIER_P4MM_MOD
   USE WYCKOFF_P4MM_MOD
   USE AOTUS_MODULE, ONLY: flu_State, OPEN_CONFIG_FILE, CLOSE_CONFIG, AOT_GET_VAL
   USE AOT_TABLE_MODULE, ONLY: AOT_TABLE_OPEN, AOT_TABLE_CLOSE, AOT_TABLE_LENGTH, AOT_TABLE_GET_VAL
   USE LINK_FUNCTION1D_MOD
#if DEBUG
   USE DEBUG_MOD
#endif
IMPLICIT NONE
!/////////////////////////////////////////////////
! TYPE: CRP6D
!
!-------------------------------------------------
TYPE,EXTENDS(PES) :: CRP6D
   INTEGER(KIND=4):: nsites
   INTEGER(KIND=4):: natomic
   LOGICAL:: is_interpolated=.FALSE.
   LOGICAL:: is_homonucl=.FALSE.
   LOGICAL:: is_smooth=.FALSE.
   LOGICAL:: is_shifted=.FALSE.
   LOGICAL:: is_resized=.FALSE.
   REAL(KIND=8) :: zvacuum
   INTEGER(KIND=4),DIMENSION(2) :: grid=[0,0]
   CLASS(Wyckoffsitio),DIMENSION(:),ALLOCATABLE:: wyckoffsite
   TYPE(CRP3D),DIMENSION(:),ALLOCATABLE:: atomiccrp
   TYPE(Vacuumpot):: farpot
   CLASS(Function1d),ALLOCATABLE:: dumpfunc
   CHARACTER(LEN=30):: extrapol2vac_flag
   INTEGER(KIND=4),DIMENSION(:,:),ALLOCATABLE:: xyklist
   CONTAINS
      ! Initialization block
      PROCEDURE,PUBLIC :: READ => READ_CRP6D
      PROCEDURE,PUBLIC :: INITIALIZE => INITIALIZE_CRP6D
      ! Set block
      PROCEDURE,PUBLIC :: SET_SMOOTH => SET_SMOOTH_CRP6D
      ! Get block
      PROCEDURE,PUBLIC :: GET_V_AND_DERIVS => GET_V_AND_DERIVS_CRP6D
      PROCEDURE,PUBLIC :: GET_V_AND_DERIVS_PURE => GET_V_AND_DERIVS_PURE_CRP6D
      PROCEDURE,PUBLIC :: GET_V_AND_DERIVS_SMOOTH => GET_V_AND_DERIVS_SMOOTH_CRP6D
      PROCEDURE,PUBLIC :: GET_ATOMICPOT_AND_DERIVS => GET_ATOMICPOT_AND_DERIVS_CRP6D
      ! Tools block
      PROCEDURE,PUBLIC :: SMOOTH => SMOOTH_CRP6D
      PROCEDURE,PUBLIC :: SMOOTH_EXTRA => SMOOTH_EXTRA_CRP6D
      PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_CRP6D
      PROCEDURE,PUBLIC :: RAWINTERPOL => RAWINTERPOL_CRP6D
      PROCEDURE,PUBLIC :: EXTRACT_VACUUMSURF => EXTRACT_VACUUMSURF_CRP6D
      PROCEDURE,PUBLIC :: ADD_VACUUMSURF => ADD_VACUUMSURF_CRP6D
      PROCEDURE,PUBLIC :: INTERPOL_NEW_RZGRID => INTERPOL_NEW_RZGRID_CRP6D
      PROCEDURE,PUBLIC :: CHEAT_CARTWHEEL_ONTOP => CHEAT_CARTWHEEL_ONTOP_CRP6D
      ! Plot tools
      PROCEDURE,PUBLIC :: PLOT1D_THETA => PLOT1D_THETA_CRP6D
      PROCEDURE,PUBLIC :: PLOT1D_THETA_SMOOTH => PLOT1D_THETA_SMOOTH_CRP6D
      PROCEDURE,PUBLIC :: PLOT1D_ATOMIC_INTERAC_THETA => PLOT1D_ATOMIC_INTERAC_THETA_CRP6D
      PROCEDURE,PUBLIC :: PLOT1D_PHI => PLOT1D_PHI_CRP6D
      PROCEDURE,PUBLIC :: PLOT1D_PHI_SMOOTH => PLOT1D_PHI_SMOOTH_CRP6D
      PROCEDURE,PUBLIC :: PLOT1D_ATOMIC_INTERAC_PHI => PLOT1D_ATOMIC_INTERAC_PHI_CRP6D
      PROCEDURE,PUBLIC :: PLOT1D_R => PLOT1D_R_CRP6D
      PROCEDURE,PUBLIC :: PLOT1D_R_SMOOTH => PLOT1D_R_SMOOTH_CRP6D
      PROCEDURE,PUBLIC :: PLOT1D_Z => PLOT1D_Z_CRP6D
      PROCEDURE,PUBLIC :: PLOT1D_Z_SMOOTH => PLOT1D_Z_SMOOTH_CRP6D
      PROCEDURE,PUBLIC :: PLOT_XYMAP => PLOT_XYMAP_CRP6D
      PROCEDURE,PUBLIC :: PLOT_RZMAP => PLOT_RZMAP_CRP6D
      PROCEDURE,PUBLIC :: PLOT_ATOMIC_INTERAC_RZ => PLOT_ATOMIC_INTERAC_RZ_CRP6D
      ! Enquire block
      PROCEDURE,PUBLIC :: is_allowed => is_allowed_CRP6D
END TYPE CRP6D
CONTAINS
!###########################################################
!# SUBROUTINE: INITIALIZE_CRP6D 
!###########################################################
!> @brief
!! Specific implementation of initialize PES from file
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INITIALIZE_CRP6D(this,filename,tablename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(OUT)::this
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: filename,tablename
   ! Local variables
   CHARACTER(LEN=:),ALLOCATABLE:: auxstring
   ! Run section
   SELECT CASE(allocated(system_inputfile) .or. .not.present(filename))
      CASE(.TRUE.)
         auxstring=system_inputfile
      CASE(.FALSE.)
         auxstring=filename
   END SELECT
   SELECT CASE(present(tablename))
      CASE(.TRUE.)
         CALL this%READ(filename=auxstring,tablename=tablename)
      CASE(.FALSE.)
         CALL this%READ(filename=auxstring,tablename='pes')
   END SELECT
   CALL this%INTERPOL()
   SELECT CASE(this%is_resized)
      CASE(.TRUE.)
         CALL this%INTERPOL_NEW_RZGRID(this%grid(1),this%grid(2))
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE INITIALIZE_CRP6D
!###########################################################
!# SUBROUTINE: GET_ATOMICPOT_AND_DERIVS_CRP6D 
!###########################################################
!> @brief
!! Get atomic potentials and derivatives
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_ATOMICPOT_AND_DERIVS_CRP6D(this,molecx,atomicx,v,dvdu)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(IN):: this
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: molecx
   REAL(KIND=8),DIMENSION(2),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(6),INTENT(OUT) :: dvdu
   REAL(KIND=8),DIMENSION(6),INTENT(OUT) :: atomicx
   ! Local variables
   REAL(KIND=8) :: ma,mb
   REAL(KIND=8) :: vcorra,vcorrb
   ! Run section
   ma=system_mass(1)
   mb=system_mass(2)
   CALL FROM_MOLECULAR_TO_ATOMIC(ma,mb,molecx,atomicx)
   SELECT CASE(this%natomic)
      CASE(1)
         CALL this%atomiccrp(1)%GET_V_AND_DERIVS(atomicx(1:3),v(1),dvdu(1:3))
         CALL this%atomiccrp(1)%GET_V_AND_DERIVS(atomicx(4:6),v(2),dvdu(4:6))
      CASE(2)
         CALL this%atomiccrp(1)%GET_V_AND_DERIVS(atomicx(1:3),v(1),dvdu(1:3))
         CALL this%atomiccrp(2)%GET_V_AND_DERIVS(atomicx(4:6),v(2),dvdu(4:6))
      CASE DEFAULT
         WRITE(0,*) "GET_ATOMICPOT_AND_DERIVS_CRP6D ERR: wrong number of atomic potentials"
         CALL EXIT(1)
   END SELECT
   vcorra=this%dumpfunc%getvalue(atomicx(3))
   vcorrb=this%dumpfunc%getvalue(atomicx(6))
   dvdu(1)=dvdu(1)*vcorra
   dvdu(2)=dvdu(2)*vcorra
   dvdu(3)=dvdu(3)*vcorra+v(1)*this%dumpfunc%getderiv(atomicx(3))
   dvdu(4)=dvdu(4)*vcorrb
   dvdu(5)=dvdu(5)*vcorrb
   dvdu(6)=dvdu(6)*vcorrb+v(2)*this%dumpfunc%getderiv(atomicx(6))
   v(1)=v(1)*vcorra
   v(2)=v(2)*vcorrb
   RETURN
END SUBROUTINE GET_ATOMICPOT_AND_DERIVS_CRP6D
!###########################################################
!# SUBROUTINE: SET_SMOOTH_CRP6D
!###########################################################
!> @brief
!! Common set function. Sets is_smooth atribute
!-----------------------------------------------------------
SUBROUTINE SET_SMOOTH_CRP6D(this,is_smooth)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(INOUT):: this
   LOGICAL,INTENT(IN) :: is_smooth
   ! Run section
   this%is_smooth=is_smooth
   RETURN
END SUBROUTINE SET_SMOOTH_CRP6D
!###########################################################
!# SUBROUTINE: CHEAT_CARTWHEEL_ONTOP_CRP6D 
!###########################################################
!> @brief
!! This subroutine changes data from a raw cut2d input by the
!! value interpolated from atomic potentials. Only data too close
!! to the plane @f$ \pi: z=\over{r}{2}+rumpling@f$
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE CHEAT_CARTWHEEL_ONTOP_CRP6D(this,wyckoff,cut2d,toptype,dmax)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O 
   CLASS(CRP6D),INTENT(INOUT) :: this
   INTEGER(KIND=4),INTENT(IN) :: wyckoff,cut2d,toptype
   REAL(KIND=8),INTENT(IN) :: dmax
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   REAL(KIND=8) :: r,z,rump,interac
   REAL(KIND=8) :: x,y ! x,y position
   REAL(KIND=8) :: m1, m2 ! masses
   REAL(KIND=8) :: dist ! distance to the plane
   INTEGER(KIND=4) :: nr,nz
   REAL(KIND=8),DIMENSION(3) :: pos1,pos2
   ! Run section
   nr=this%wyckoffsite(wyckoff)%zrcut(cut2d)%getgridsizeR()
   nz=this%wyckoffsite(wyckoff)%zrcut(cut2d)%getgridsizeZ()
   x=this%wyckoffsite(wyckoff)%x
   y=this%wyckoffsite(wyckoff)%y
   m1=system_mass(1)
   m2=system_mass(2)
   rump=this%atomiccrp(1)%getrumpling(toptype)
   DO i = 1, nr
      DO j = 1, nz
        r=this%wyckoffsite(wyckoff)%zrcut(cut2d)%getgridvalueR(i)
        z=this%wyckoffsite(wyckoff)%zrcut(cut2d)%getgridvalueZ(j)
        dist=(2.D0/sqrt(5.D0))*(z+rump-r/2.D0)
        pos1=(/x,y,z+(m2/(m1+m2))*r/)
        pos2=(/x,y,z-(m1/(m1+m2))*r/)
        SELECT CASE(dist <= dmax) 
           CASE(.TRUE.)
              SELECT CASE(this%is_homonucl)
                 CASE(.TRUE.)
                    interac=this%atomiccrp(1)%getpot(pos1)+this%atomiccrp(1)%getpot(pos2)!+this%farpot%getpot(r)
                 CASE(.FALSE.)
                    interac=this%atomiccrp(1)%getpot(pos1)+this%atomiccrp(2)%getpot(pos2)!+this%farpot%getpot(r)
              END SELECT
              CALL this%wyckoffsite(wyckoff)%zrcut(cut2d)%CHANGEPOT_AT_GRIDPOINT(i,j,interac)
           CASE(.FALSE.)
              ! do nothing
        END SELECT
      END DO
   END DO
   RETURN
END SUBROUTINE CHEAT_CARTWHEEL_ONTOP_CRP6D
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
!# SUBROUTINE: GET_V_AND_DERIVS_PURE_CRP6D 
!###########################################################
!> @brief
!! Gets the potential and derivatives respect to all DOFs if we
!! are in the pure CRP6D region. Otherwise, there'd be errors and
!! unexpected behaviors.
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
SUBROUTINE GET_V_AND_DERIVS_PURE_CRP6D(this,x,v,dvdu)
   ! Initial declarations   
   USE DEBUG_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),TARGET,INTENT(IN):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   REAL(KIND=8),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: dvdu
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   REAL(KIND=8) :: ma,mb
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: f ! smooth function and derivs
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: dfdu ! smooth derivatives
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: xy
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: aux1
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: aux2
   REAL(KIND=8),DIMENSION(6) :: atomicx
   REAL(KIND=8),DIMENSION(2) :: atomic_v
   REAL(KIND=8),DIMENSION(6) :: atomic_dvdu
   REAL(KIND=8),DIMENSION(3) :: dvdu_atomicA,dvdu_atomicB
   TYPE(Fourierp4mm) :: xyinterpol
   CHARACTER(LEN=24),PARAMETER :: routinename="GET_V_AND_DERIVS_CRP6D: "
   ! Run section
   SELECT CASE(this%is_smooth)
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         WRITE(0,*) "GET_V_AND_DERIVS_CRP6D ERR: smooth the PES first (CALL thispes%SMOOTH())"
         CALl EXIT(1)
   END SELECT
   ALLOCATE(f(5,this%nsites))
   ALLOCATE(xy(this%nsites,2))
   DO i = 1, this%nsites 
#ifdef DEBUG
      CALL VERBOSE_SEPARATOR1()
      CALL VERBOSE_WRITE(routinename,"WYCKOFF SITE: ",i)
      CALL VERBOSE_WRITE(routinename,"Letter: ",this%wyckoffsite(i)%id)
#endif
      CALL this%wyckoffsite(i)%GET_V_AND_DERIVS(x(3:6),f(1,i),f(2:5,i))
      xy(i,1)=this%wyckoffsite(i)%x
      xy(i,2)=this%wyckoffsite(i)%y
   END DO
   ! f(1,:) smooth potential values
   ! f(2,:) smooth dvdz
   ! f(3,:) smooth dvdr
   ! f(4,:) smooth dvdtheta
   ! f(5,:) smooth dvdphi
   CALL xyinterpol%READ(xy,f,this%xyklist)
   CALL xyinterpol%INTERPOL(this%atomiccrp(1)%surf)
   ALLOCATE(aux1(5))
   ALLOCATE(aux2(5,2))
   CALL xyinterpol%GET_F_AND_DERIVS(this%atomiccrp(1)%surf,x(1:2),aux1,aux2)
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
   CALL VERBOSE_WRITE(routinename,"List of XY: ")
   DO i = 1,this%nsites 
     CALL VERBOSE_WRITE(routinename,xy(i,:))
   END DO
   CALL VERBOSE_WRITE(routinename,"Smooth V at each wyckoff site:")
   CALL VERBOSE_WRITE(routinename,f(1,:))
   CALL VERBOSE_WRITE(routinename,"Smooth dv/dz at each wyckoff site:")
   CALL VERBOSE_WRITE(routinename,f(2,:))
   CALL VERBOSE_WRITE(routinename,"Smooth dv/dr at each wyckoff site:")
   CALL VERBOSE_WRITE(routinename,f(3,:))
   CALL VERBOSE_WRITE(routinename,"Smooth dv/dtheta at each wyckoff site:")
   CALL VERBOSE_WRITE(routinename,f(4,:))
   CALL VERBOSE_WRITE(routinename,"Smooth dv/dphi at each wyckoff site:")
   CALL VERBOSE_WRITE(routinename,f(5,:))
   CALL VERBOSE_WRITE(routinename,"Smooth interpolated values:")
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
   CALL this%GET_ATOMICPOT_AND_DERIVS(x,atomicx,atomic_v,atomic_dvdu)
   dvdu_atomicA=atomic_dvdu(1:3)
   dvdu_atomicB=atomic_dvdu(4:6)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Contributions of the atomic potential: ")
   CALL VERBOSE_WRITE(routinename, "Position Atom A: ",atomicx(1:3))
   CALL VERBOSE_WRITE(routinename,"Va: ",atomic_v(1))
   CALL VERBOSE_WRITE(routinename,"dVa/dxa; dVa/dya: ",dvdu_atomicA)
   CALL VERBOSE_WRITE(routinename, "Position Atom B: ",atomicx(4:6))
   CALL VERBOSE_WRITE(routinename,"Vb: ",atomic_v(2))
   CALL VERBOSE_WRITE(routinename,"dVa/dxa; dVb/dyb: ",dvdu_atomicB)
#endif
   v=aux1(1)+sum(atomic_v)

   ma=system_mass(1)
   mb=system_mass(2)
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
END SUBROUTINE GET_V_AND_DERIVS_PURE_CRP6D
!###########################################################
!# SUBROUTINE: GET_V AND DERIVS_CRP6D 
!###########################################################
!> @brief
!! Gets potential for any configuration. Discriminates between 3 regions:
!! - Pure CRP6D region: zmin to zcrp
!! - Extrapolation region zcrp to zvacuum
!! - Vacuum region: further than zvacuum
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_CRP6D(this,X,v,dvdu)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),TARGET,INTENT(IN):: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   REAL(KIND=8),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: dvdu
   ! Local variables
   REAL(KIND=8) :: zcrp, zvac ! last CRP6D z value and Z infinity
   REAL(KIND=8) :: vzcrp, vzvac ! potentials at zcrp and zvac
   REAL(KIND=8),DIMENSION(6) :: dvducrp ! derivatives at zcrp
   REAL(KIND=8),DIMENSION(6) :: dvduvac ! derivatives at vacuum
   REAL(KIND=8) :: alpha,beta,gama ! parameters
   CLASS(Function1d),ALLOCATABLE:: extrapolfunc
   INTEGER(KIND=4) :: i !counter
   ! Local Parameter
   REAL(KIND=8),PARAMETER :: zero=0.D0 ! what we will consider zero (a.u.)
   REAL(KIND=8),PARAMETER :: dz=0.5D0 ! 0.25 Angstroems approx
   REAL(KIND=8) :: dummy
   ! Run section
   zcrp=this%wyckoffsite(1)%zrcut(1)%getlastZ()
   zvac=this%zvacuum
   ! Check if we are in the pure CRP6D region
   SELECT CASE(x(3)<= zcrp) !easy
      CASE(.TRUE.)
         CALL this%GET_V_AND_DERIVS_PURE(x,v,dvdu)
         RETURN
      CASE(.FALSE.)
         ! do nothing, next switch
   END SELECT
   ! Check if we are in the extrapolation region
   SELECT CASE(x(3)>zcrp .AND. x(3)<zvac)
      CASE(.TRUE.) ! uff
         ! Set potential and derivs
         vzvac=this%farpot%getpot(x(4))
         CALL this%GET_V_AND_DERIVS_PURE([x(1),x(2),zcrp,x(4),x(5),x(6)],vzcrp,dvducrp)
         dvduvac(1:3)=zero
         dvduvac(4)=this%farpot%getderiv(x(4))
         dvduvac(5:6)=zero
         ! Check kind of extrapolation
         SELECT CASE(this%extrapol2vac_flag)
            CASE("Xexponential")
               ALLOCATE(Xexponential_func::extrapolfunc)
               ! Extrapol potential
               beta=-1.D0/zvac
               alpha=(vzcrp-vzvac)/(zcrp*dexp(beta*zcrp)-zvac*dexp(beta*zvac))
               gama=vzvac-alpha*zvac*dexp(beta*zvac)
               CALL extrapolfunc%READ([alpha,beta])
               v=extrapolfunc%getvalue(x(3))+gama
               dvdu(3)=extrapolfunc%getderiv(x(3))
               ! Extrapol derivatives
               DO i = 1, 6
                  SELECT CASE(i)
                     CASE(3)
                        ! Skip dvdz
                     CASE DEFAULT
                        beta=-1.D0/zvac
                        alpha=(dvducrp(i)-dvduvac(i))/(zcrp*dexp(beta*zcrp)-zvac*dexp(beta*zvac))
                        gama=dvduvac(i)-alpha*zvac*dexp(beta*zvac)
                        CALL extrapolfunc%READ([alpha,beta])
                        dvdu(i)=extrapolfunc%getvalue(x(3))+gama
                  END SELECT
               END DO
               DEALLOCATE(extrapolfunc)
               RETURN

            CASE("Linear")
               ALLOCATE(Linear_func::extrapolfunc)
               ! Extrapol potential
               beta=(vzcrp*zvac-vzvac*zcrp)/(zvac-zcrp)
               alpha=(vzvac-beta)/zvac
               CALL extrapolfunc%READ([alpha,beta])
               v=extrapolfunc%getvalue(x(3))
               dvdu(3)=extrapolfunc%getderiv(x(3))
               ! Extrapol derivatives
               DO i = 1, 6
                  SELECT CASE(i)
                     CASE(3)
                        ! do nothing
                     CASE DEFAULT
                        beta=(dvducrp(i)*zvac-dvduvac(i)*zcrp)/(zvac-zcrp)
                        alpha=(dvduvac(i)-beta)/zvac
                        CALL extrapolfunc%READ([alpha,beta])
                        dvdu(i)=extrapolfunc%getvalue(x(3))
                  END SELECT
               END DO
               DEALLOCATE(extrapolfunc)
               RETURN

            CASE DEFAULT
               WRITE(0,*) "GET_V_AND_DERIVS_CRP6D ERR: type of extrapolation function isn't implemented yet"
               WRITE(0,*) "Implemented ones: Linear, Xexponential, None"
               CALL EXIT(1)
         END SELECT
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   ! Check if we are in the Vacuum region
    SELECT CASE(x(3)>=zvac) !easy
       CASE(.TRUE.)
          v=this%farpot%getpot(x(4))
          dvdu(1:3)=0.D0
          dvdu(4)=this%farpot%getderiv(x(4))
          dvdu(5:6)=0.D0
       CASE(.FALSE.) ! this's the last switch!
          WRITE(0,*) "GET_V_AND_DERIVS_CRP6D ERR: unclassificable point. Where are we?"
          WRITE(0,*) "Asking potential at Z: ", x(3)
          WRITE(0,*) "Zvacuum: ",zvac
          WRITE(0,*) "Zcrp: ",zcrp
          CALL EXIT(1)
    END SELECT
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_CRP6D
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_SMOOTH_CRP6D 
!###########################################################
!> @brief
!! Gets the SMOOTH potential and the derivatives respect to all DOFs
!
!> @details
!! This routine doesn't check PES smoothness, so it can be used to get
!! interpolated values for the raw PES
!> @warning
!! - Assumed a 2D Fourier interpolation in XY (symmetry adapted)
!! - Assumed a 1D Fourier interpolation in theta (symmetry adapted)
!! - Assumed a 1D Fourier interpolation in phi (symmetry adapted)
!! - Assumed a bicubic splines interpolation in R-Z (symmetry independent)
!! - Assumed that Z-R grids have the same size. This can be ensured using
!!   using subroutine INTERPOL_NEWGRID
!
!> @see
!! interpol_newgrid_crp6d
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_SMOOTH_CRP6D(this,x,v,dvdu)
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
   TYPE(Fourierp4mm) :: xyinterpol
   CHARACTER(LEN=24),PARAMETER :: routinename="GET_V_AND_DERIVS_CRP6D: "
   ! Run section
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
      CALL VERBOSE_SEPARATOR1()
      CALL VERBOSE_WRITE(routinename,"WYCKOFFSITE: ",i)
      CALL VERBOSE_WRITE(routinename,"Letter:",this%wyckoffsite(i)%id)
#endif
      CALL this%wyckoffsite(i)%GET_V_AND_DERIVS(x(3:6),f(1,i),f(2:5,i))
      xy(i,1)=this%wyckoffsite(i)%x
      xy(i,2)=this%wyckoffsite(i)%y
   END DO
   ! f(1,:) smooth potential values
   ! f(2,:) smooth dvdz
   ! f(3,:) smooth dvdr
   ! f(4,:) smooth dvdtheta
   ! f(5,:) smooth dvdphi
   CALL xyinterpol%READ(xy,f,this%xyklist)
   CALL xyinterpol%INTERPOL(this%atomiccrp(1)%surf)
   ALLOCATE(aux1(5))
   ALLOCATE(aux2(5,2))
   CALL xyinterpol%GET_F_AND_DERIVS(this%atomiccrp(1)%surf,x(1:2),aux1,aux2)
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
   CALL VERBOSE_WRITE(routinename,"List of XY: ")
   DO i = 1,this%nsites 
     CALL VERBOSE_WRITE(routinename,xy(i,:))
   END DO
   CALL VERBOSE_WRITE(routinename,"Pots at each wyckoff site:")
   CALL VERBOSE_WRITE(routinename,f(1,:))
   CALL VERBOSE_WRITE(routinename,"dv/dz at each wyckoff site:")
   CALL VERBOSE_WRITE(routinename,f(2,:))
   CALL VERBOSE_WRITE(routinename,"dv/dr at each wyckoff site:")
   CALL VERBOSE_WRITE(routinename,f(3,:))
   CALL VERBOSE_WRITE(routinename,"dv/dtheta at each wyckoff site:")
   CALL VERBOSE_WRITE(routinename,f(4,:))
   CALL VERBOSE_WRITE(routinename,"dv/dphi at each wyckoff site:")
   CALL VERBOSE_WRITE(routinename,f(5,:))
#endif
   v=aux1(1)
   dvdu(1)=aux2(1,1)
   dvdu(2)=aux2(1,2)
   dvdu(3)=aux1(2)
   dvdu(4)=aux1(3)
   dvdu(5)=aux1(4)
   dvdu(6)=aux1(5)
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_SMOOTH_CRP6D
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
   !CALL this%SMOOTH_EXTRA()
   DO i = 1, this%nsites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         CALL this%wyckoffsite(i)%zrcut(j)%INTERPOL()
      END DO
   END DO
   this%is_interpolated=.TRUE.
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
   CALL this%EXTRACT_VACUUMSURF()   
   DO i = 1, this%nsites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         CALL this%wyckoffsite(i)%zrcut(j)%INTERPOL()
      END DO
   END DO
   this%is_interpolated=.TRUE.
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
!# SUBROUTINE: FROM_ATOMIC_TO_MOLECULAR 
!###########################################################
!> @brief
!! Goes from atomic coordinates xa,ya,za,xb,yb,zb to molecular
!! coordinates x,y,z,r,theta,phi.
!> @details
!! - We have enforced @f$\theta \in [0,\pi]@f$ and @f$\phi \in [0,2\pi)@f$
!-----------------------------------------------------------
SUBROUTINE FROM_ATOMIC_TO_MOLECULAR(ma,mb,atomcoord,molcoord)
   ! Initial declarations
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables
   REAL(KIND=8),INTENT(IN) :: ma,mb
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: atomcoord
   REAL(KIND=8),DIMENSION(6),INTENT(OUT) :: molcoord
   ! Local variables
   ! Run section
   molcoord(1)=(1.D0/(ma+mb))*(atomcoord(1)*ma+atomcoord(4)*mb)
   molcoord(2)=(1.D0/(ma+mb))*(atomcoord(2)*ma+atomcoord(5)*mb)
   molcoord(3)=(1.D0/(ma+mb))*(atomcoord(3)*ma+atomcoord(6)*mb)
   molcoord(4)=dsqrt((atomcoord(1)-atomcoord(4))**2.D0+&
      (atomcoord(2)-atomcoord(5))**2.D0+(atomcoord(3)-atomcoord(6))**2.D0)
   molcoord(5)=dacos((atomcoord(3)-atomcoord(6))/molcoord(4))
   SELECT CASE(atomcoord(1)<atomcoord(4)) ! II or III Quadrant
      CASE(.TRUE.)
         molcoord(6)=PI-datan((atomcoord(2)-atomcoord(5))/(atomcoord(1)-atomcoord(4)))
         RETURN
      CASE(.FALSE.)
         ! do nothing   
   END SELECT
   SELECT CASE(atomcoord(1)>atomcoord(4) .AND. atomcoord(2)>=atomcoord(5)) ! I Quadrant
      CASE(.TRUE.)
         molcoord(6)=datan((atomcoord(2)-atomcoord(5))/(atomcoord(1)-atomcoord(4)))
         RETURN
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(atomcoord(1)>atomcoord(4) .AND. atomcoord(2)<atomcoord(5)) ! IV Quadrant
      CASE(.TRUE.)
         molcoord(6)=2.D0*PI+datan((atomcoord(2)-atomcoord(5))/(atomcoord(1)-atomcoord(4)))
         RETURN
      CASE(.FALSE.)
         !do nothing
   END SELECT
   SELECT CASE(atomcoord(1)==atomcoord(4) .AND. atomcoord(2)>atomcoord(5))
      CASE(.TRUE.)
         molcoord(6)=PI/2.D0
      CASE(.FALSE.)
         !do nothing
   END SELECT
   SELECT CASE(atomcoord(1)==atomcoord(4) .AND. atomcoord(2)<atomcoord(5))
      CASE(.TRUE.)
         molcoord(6)=3.D0*PI/2.D0
      CASE(.FALSE.)
         !do nothing
   END SELECT
   SELECT CASE(atomcoord(1)==atomcoord(4) .AND. atomcoord(2)==atomcoord(5)) ! cartwheel, cannot be defined
      CASE(.TRUE.)
         molcoord(6)=0.D0 
      CASE(.FALSE.)
         WRITE(*,*) "FROM_ATOMIC_TO_MOLECULAR ERR: Geometry was not taken into account, check code"
         CALL EXIT(1)
   END SELECT
   RETURN
END SUBROUTINE FROM_ATOMIC_TO_MOLECULAR
!###########################################################
!# SUBROUTINE: SMOOTH_CRP6D
!###########################################################
!> @brief
!! Smooths a RZ-2dcut of the potential using atomic potentials
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
   REAL(KIND=8),DIMENSION(6) :: molcoord,atomcoord,dvdu
   INTEGER(KIND=4) :: nr,nz
   INTEGER(KIND=4) :: i,j,k,l ! counters
   CHARACTER(LEN=14),PARAMETER :: routinename="SMOOTH_CRP6D: "
   REAL(KIND=8) :: newpot
   REAL(KIND=8),DIMENSION(2) :: atomic_v
   ! Run section
   DO i = 1, this%nsites ! cycle wyckoff sites
      molcoord(1)=this%wyckoffsite(i)%x
      molcoord(2)=this%wyckoffsite(i)%y
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         molcoord(5)=this%wyckoffsite(i)%zrcut(j)%theta
         molcoord(6)=this%wyckoffsite(i)%zrcut(j)%phi
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizeR()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizeZ()
         DO k = 1, nr
            DO l = 1, nz
               molcoord(3)=this%wyckoffsite(i)%zrcut(j)%getgridvalueZ(l)
               molcoord(4)=this%wyckoffsite(i)%zrcut(j)%getgridvalueR(k)
               CALL this%GET_ATOMICPOT_AND_DERIVS(molcoord,atomcoord,atomic_v,dvdu)
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)-sum(atomic_v)
               CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
            END DO
         END DO
      END DO
   END DO
   this%is_smooth=.TRUE.
   RETURN
END SUBROUTINE SMOOTH_CRP6D
!###########################################################
!# SUBROUTINE: ROUGH_CRP6D
!###########################################################
!> @brief
!! ROUGHs a RZ-2dcut of the potential using atomic potentials
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE ROUGH_CRP6D(this)
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
   CHARACTER(LEN=13),PARAMETER :: routinename="ROUGH_CRP6D: "
   REAL(KIND=8) :: newpot
   ! Run section
   ma=system_mass(1)
   mb=system_mass(2)
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
                     newpot=newpot+this%atomiccrp(1)%getpot(atomcoord(1:3))+&
                        this%atomiccrp(1)%getpot(atomcoord(4:6))
                     CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  CASE(2)
                     newpot=newpot+this%atomiccrp(1)%getpot(atomcoord(1:3))+&
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
   this%is_smooth=.FALSE.
   RETURN
END SUBROUTINE ROUGH_CRP6D
!###########################################################
!# SUBROUTINE: SMOOTH_EXTRA_CRP6D
!###########################################################
!> @brief
!! Smooths a RZ-2dcut of the potential using atomic potentials as well
!! as the vacuum potential
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SMOOTH_EXTRA_CRP6D(this)
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
   CHARACTER(LEN=20),PARAMETER :: routinename="SMOOTH_EXTRA_CRP6D: "
   REAL(KIND=8) :: newpot
   ! Run section
   ma=system_mass(1)
   mb=system_mass(2)
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
                     newpot=newpot-this%atomiccrp(1)%getpot(atomcoord(1:3))&
                        -this%atomiccrp(1)%getpot(atomcoord(4:6))-this%farpot%getpot(molcoord(4))&
                        -0.8*dexp(-molcoord(4))
                     CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  CASE(2)
                     newpot=newpot-this%atomiccrp(1)%getpot(atomcoord(1:3))-&
                        this%atomiccrp(2)%getpot(atomcoord(4:6))-this%farpot%getpot(molcoord(4))
                     CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
                  CASE DEFAULT
                     WRITE(0,*) "SMOOTH_EXTRA_CRP6D ERR: Something is wrong with number of atomic crp potentials"
                     CALL EXIT(1)
               END SELECT
            END DO
         END DO
      END DO
   END DO
   this%is_smooth=.TRUE.
   RETURN
END SUBROUTINE SMOOTH_EXTRA_CRP6D
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
   REAL(KIND=8) :: newpot
   CHARACTER(LEN=26),PARAMETER :: routinename="EXTRACT_VACUUMSURF_CRP6D: "
   ! Run section
   DO i = 1, this%nsites ! cycle wyckoff sites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizeR()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizeZ()
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
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Potential shifted: ",this%farpot%getscalefactor())
#endif
   this%is_shifted=.TRUE.
   RETURN
END SUBROUTINE EXTRACT_VACUUMSURF_CRP6D
!###########################################################
!# SUBROUTINE: ADD_VACUUMSURF_CRP6D
!###########################################################
!> @brief
!! Adds energy at equilibrium distance in the vacuum and surface energy
!! to an Rz-2dcut of the potential. Shifts the vacuum potential as well.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date May/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE ADD_VACUUMSURF_CRP6D(this)
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
   REAL(KIND=8) :: newpot
   CHARACTER(LEN=22),PARAMETER :: routinename="ADD_VACUUMSURF_CRP6D: "
   ! Run section
   DO i = 1, this%nsites ! cycle wyckoff sites
      DO j = 1, this%wyckoffsite(i)%n2dcuts
         nr=this%wyckoffsite(i)%zrcut(j)%getgridsizeR()
         nz=this%wyckoffsite(i)%zrcut(j)%getgridsizeZ()
         DO k = 1, nr
            DO l = 1, nz
               newpot=this%wyckoffsite(i)%zrcut(j)%getpotatgridpoint(k,l)
               newpot=newpot+this%farpot%getscalefactor()
               CALL this%wyckoffsite(i)%zrcut(j)%CHANGEPOT_AT_GRIDPOINT(k,l,newpot)
            END DO
         END DO
      END DO
   END DO
   CALL this%farpot%SHIFTPOT_UNDO()
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Potential re-shifted: ",this%farpot%getscalefactor())
#endif
   this%is_shifted=.FALSE.
   RETURN
END SUBROUTINE ADD_VACUUMSURF_CRP6D
!###########################################################
!# SUBROUTINE: READ_CRP6D 
!###########################################################
!> @brief
!! Sets up a CRP6D Object
!-----------------------------------------------------------
SUBROUTINE READ_CRP6D(this,filename,tablename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(OUT):: this
   CHARACTER(LEN=*),INTENT(IN):: filename,tablename
   ! IMPORTANT: unit used to read
   INTEGER(KIND=4),PARAMETER:: runit=180
   ! Local variables
   INTEGER(KIND=4):: i,j ! counters
   INTEGER(KIND=4):: natomic ! number of atomic potentials
   REAL(KIND=8):: auxr
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: param_dump
   CHARACTER(LEN=30):: cut2dfilename, string, surfacefile
   CHARACTER,DIMENSION(:),ALLOCATABLE:: letter
   REAL(KIND=8),DIMENSION(2):: masslist
   CHARACTER(LEN=2),DIMENSION(2):: symbollist
   CHARACTER(LEN=10):: units
   CHARACTER(LEN=6):: resize
   TYPE(Length):: len
   ! Lua specifications
   TYPE(flu_State):: conf
   INTEGER(KIND=4):: ierr
   INTEGER(KIND=4):: pes_table,crp3d_table,vacfunc_table
   ! Auxiliar, dummy variables
   INTEGER(KIND=4):: auxint
   CHARACTER(LEN=1024):: auxstring
   ! Parameters
   CHARACTER(LEN=*),PARAMETER:: routinename="READ_CRP6D: "
   ! Run section -----------------------
   ! Open Lua config file
   CALL OPEN_CONFIG_FILE(L=conf,ErrCode=ierr,filename=filename)
   ! Open PES table
   CALL AOT_TABLE_OPEN(L=conf,thandle=pes_table,key=tablename)
   ! get pes.name
   CALL AOT_GET_VAL(L=conf,thandle=pes_table,key='name',val=auxstring)
   CALL this%SET_ALIAS(trim(auxstring))
   ! get pes.kind
   CALL AOT_GET_VAL(L=conf,thandle=pes_table,key='kind',val=auxstring)
   CALL this%SET_PESTYPE(trim(auxstring))
   ! get pes.dimensions
   CALL AOT_GET_VAL(L=conf,thandle=pes_table,key='dimensions',val=auxint)
   CALL this%SET_DIMENSIONS(auxint)
   ! deal with surface
   CALL AOT_GET_VAL(L=conf,thandle=pes_table,key='surface',val=auxstring)
   CALL this%surf%INITIALIZE(trim(auxstring))
   ! get crp3d subpes
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=crp3d_table,key='crp3dPes')
   this%natomic=aot_table_length(L=conf,thandle=crp3d_table)
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
   DO i = 1, this%natomic
      CALL AOT_TABLE_GET_VALUE(L=conf,ErrCode=ierr,thandle=crp3d_table,pos=i,val=auxstring)
      CALL this%atomiccrp(i)%INITIALIZE(filename=filename,tablename=trim(auxstring))
   END DO
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Setting up: ",this%getalias())
   CALL VERBOSE_WRITE(routinename,"Atomic potentials found: ",this%natomic)
   CALL VERBOSE_WRITE(routinename,"Load from file: ",string)
#endif
   DO i = 1, this%natomic
      SELECT CASE(this%atomiccrp(1)%surf%tellfilename()==this%atomiccrp(i)%surf%tellfilename())
         CASE(.FALSE.)
            WRITE(0,*) "READ_CRP6D ERR: not all atomic potentials have the same surface file"
            CALL EXIT(1)
         CASE(.TRUE.)
            ! do nothing
      END SELECT
   END DO
   CALL AOT_TABLE_CLOSE(L=conf,thandle=crp3d_table)
   ! Read vacuum potential file -----------------------
   CALL AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=vacfunc_table,key='vacuumFunction')
   READ(runit,*) string
   CALL this%farpot%INITIALIZE(string)
   ! Read dumping function ------------------------
   READ(runit,*) string
   SELECT CASE(string)
      CASE("Logistic")
         ALLOCATE(Logistic_func::this%dumpfunc)
         ALLOCATE(param_dump(2))
         READ(runit,*) param_dump
         CALL this%dumpfunc%READ(param_dump)
         DEALLOCATE(param_dump)
      CASE("None") 
         ALLOCATE(One_func::this%dumpfunc)
         READ(runit,*) ! dummy line
      CASE DEFAULT
         WRITE(0,*) "READ_CRP6D ERR: Keyword for dumping function needed"
         WRITE(*,*) "Currently implemented: Logistic, None"
         CALL EXIT(1)
   END SELECT
   ! Read extrapolation function
   READ(runit,*) string
   SELECT CASE(string)
      CASE("Xexponential")
         this%extrapol2vac_flag=string
         READ(runit,*) auxr,units
         CALL len%READ(auxr,units)
         CALL len%TO_STD()
         this%zvacuum=len%getvalue()
      CASE("Linear")
         this%extrapol2vac_flag=string
         READ(runit,*) auxr,units
         CALL len%READ(auxr,units)
         CALL len%TO_STD()
         this%zvacuum=len%getvalue()
      CASE("None")
         READ(runit,*) ! dummy line 
         this%zvacuum=0.D0
      CASE DEFAULT
         WRITE(0,*) "READ_CRP6D ERR: Keyword for extrapolation function needed"
         WRITE(*,*) "Currently implemented: None, Logistic, Exponential, Xexponential"
         WRITE(*,*) "Case sensitive."
         CALL EXIT(1)
   END SELECT
   ! Read zvacuum
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"Z vacuum: ",this%zvacuum)
#endif
   ! Read Resize options
   READ(runit,'(A6,1X,L1)',advance="no") resize,this%is_resized
   SELECT CASE( resize=="RESIZE" .AND. this%is_resized )
      CASE(.TRUE.)
         READ(runit,*) this%grid(:)
      CASE(.FALSE.)
         READ(runit,*) ! just ignore the rest of the line
   END SELECT
   ! Read number of wyckoff sites and its letters ----------------
   READ(runit,'(I2)',advance="no") this%nsites
   ALLOCATE(letter(this%nsites))
   READ(runit,*) letter(:)
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
      CALL this%wyckoffsite(i)%READ(runit)
      CALL this%wyckoffsite(i)%SET_ID(letter(i))
      CALL this%wyckoffsite(i)%SET_HOMONUCL(this%is_homonucl)
      CALL this%wyckoffsite(i)%SET_MYNUMBER(i)
   END DO
   ! Read the final part of the input: kpoints for XY interpolation
   ALLOCATE(this%xyklist(this%nsites,2))
   READ(runit,*) ! dummy line
   DO i = 1, this%nsites
      READ(runit,*) this%xyklist(i,:)
   END DO
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"xyklist found:")
   DO i = 1, this%nsites
      CALL VERBOSE_WRITE(routinename,this%xyklist(i,:))
   END DO
#endif
   CLOSE(runit)
   RETURN
END SUBROUTINE READ_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_PHI_CRP6D #######################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0. Initial PHI value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 09/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_PHI_CRP6D(thispes,npoints,X,filename)
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER,INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin, ymax, ymin 
   REAL(KIND=8),DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=18),PARAMETER :: routinename = "PLOT1D_PHI_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_PHI_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:5)=x(1:5)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(6)=xmin
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(6),v,dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(6)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(6),v,dvdu(:)
   END DO
   ! Final value
   r(6) = xmax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(6),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_PHI_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_ATOMIC_INTERAC_PHI_CRP6D #######################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the atomic interaction
!! PES
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0. Initial PHI value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_ATOMIC_INTERAC_PHI_CRP6D(thispes,npoints,X,filename)
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER,INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta
   REAL(KIND=8),DIMENSION(2) :: v
   REAL(KIND=8) :: xmax, xmin, ymax, ymin 
   REAL(KIND=8),DIMENSION(6) :: r, dvdu, ratom
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=18),PARAMETER :: routinename = "PLOT1D_PHI_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_PHI_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:5)=x(1:5)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(6)=xmin
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
   WRITE(11,*) r(6),sum(v),v(1),v(2),ratom(:),dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(6)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
      WRITE(11,*) r(6),sum(v),v(1),v(2),ratom(:),dvdu(:)  
   END DO
   ! Final value
   r(6) = xmax
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
   WRITE(11,*) r(6),sum(v),v(1),v(2),ratom(:),dvdu(:)  
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_ATOMIC_INTERAC_PHI_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_PHI_SMOOTH_CRP6D #######################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the smooth PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0. Initial PHI value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 09/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_PHI_SMOOTH_CRP6D(thispes,npoints,X,filename)
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER,INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin, ymax, ymin 
   REAL(KIND=8),DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=25),PARAMETER :: routinename = "PLOT1D_PHI_SMOOTH_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_PHI_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:5)=x(1:5)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(6)=xmin
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(6),v,dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(6)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(6),v,dvdu(:)
   END DO
   ! Final value
   r(6) = xmax
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(6),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_PHI_SMOOTH_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_THETA_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial THETA value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_THETA_CRP6D(thispes,npoints,X,filename)
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin, ymax, ymin 
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=20),PARAMETER :: routinename = "PLOT1D_THETA_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_THETA_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:4)=x(1:4)
   r(6)=x(6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(5)=xmin
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(5),v,dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(5)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(5),v,dvdu(:)
   END DO
   ! Final value
   r(5) = xmax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(5),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_THETA_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_ATOMIC_INTERAC_THETA_CRP6D ########################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the atomic PES 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial THETA value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_ATOMIC_INTERAC_THETA_CRP6D(thispes,npoints,X,filename)
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta
   REAL(KIND=8) :: xmax, xmin, ymax, ymin 
   REAL(KIND=8), DIMENSION(6) :: r,dvdu,ratom
   REAL(KIND=8),DIMENSION(2) :: v
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=20),PARAMETER :: routinename = "PLOT1D_THETA_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_THETA_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:4)=x(1:4)
   r(6)=x(6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(5)=xmin
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
   WRITE(11,*) r(5),sum(v),v(:),ratom(:),dvdu(:)  

   ! cycle for inpoints
   DO i=1, inpoints
      r(5)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
      WRITE(11,*) r(5),sum(v),v(:),ratom(:),dvdu(:)
   END DO
   ! Final value
   r(5) = xmax
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,ratom,v,dvdu)
   WRITE(11,*) r(5),sum(v),v(:),ratom(:),dvdu(:)  
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_ATOMIC_INTERAC_THETA_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_THETA_SMOOTH_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the smooth PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial THETA value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_THETA_SMOOTH_CRP6D(thispes,npoints,X,filename)
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin, ymax, ymin 
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=27),PARAMETER :: routinename = "PLOT1D_THETA_SMOOTH_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_THETA_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = 0.D0
   xmax = 2.D0*PI
   r(1:4)=x(1:4)
   r(6)=x(6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=2.D0*PI/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(5)=xmin
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(5),v,dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(5)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(5),v,dvdu(:)
   END DO
   ! Final value
   r(5) = xmax
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(5),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_THETA_SMOOTH_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_R_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial R value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_R_CRP6D(thispes,npoints,X,filename)
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin, ymax, ymin 
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=16),PARAMETER :: routinename = "PLOT1D_R_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_R_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = thispes%wyckoffsite(1)%zrcut(1)%getfirstr()
   xmax = thispes%wyckoffsite(1)%zrcut(1)%getlastr()
   r(1:3)=x(1:3)
   r(5:6)=x(5:6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=(xmax-xmin)/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(4)=xmin
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),v,dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(4)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(4),v,dvdu(:)
   END DO
   ! Final value
   r(4) = xmax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_R_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_R_SMOOTH_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the smooth PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial R value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_R_SMOOTH_CRP6D(thispes,npoints,X,filename)
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin, ymax, ymin 
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=23),PARAMETER :: routinename = "PLOT1D_R_SMOOTH_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_R_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = thispes%wyckoffsite(1)%zrcut(1)%getfirstr()
   xmax = thispes%wyckoffsite(1)%zrcut(1)%getlastr()
   r(1:3)=x(1:3)
   r(5:6)=x(5:6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=(xmax-xmin)/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(4)=xmin
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(4),v,dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(4)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(4),v,dvdu(:)
   END DO
   ! Final value
   r(4) = xmax
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(4),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_R_SMOOTH_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_Z_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!> @param[in] L - Length to plot from minimim Z
!
!> @warning
!! - The graph starts always at 0,0, Initial Z value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_Z_CRP6D(thispes,npoints,X,L,filename)
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   REAL(KIND=8),INTENT(IN) :: L
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin, ymax, ymin 
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=16),PARAMETER :: routinename = "PLOT1D_Z_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   IF (npoints.lt.2) THEN
      WRITE(0,*) "PLOT1D_Z_CRP6D ERR: Less than 2 points"
      CALL EXIT(1)
   END IF
   !
   xmin = thispes%wyckoffsite(1)%zrcut(1)%getfirstz()
   xmax = xmin+L
   r(1:2)=x(1:2)
   r(4:6)=x(4:6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=(xmax-xmin)/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(3)=xmin
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(3)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(3),v,dvdu(:)
   END DO
   ! Final value
   r(3) = xmax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_Z_CRP6D
!#######################################################################
! SUBROUTINE: PLOT1D_Z_SMOOTH_CRP6D #########################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the smooth PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!> @param[in] x - array with X,Y,Z,R,THETA,PHI values
!
!> @warning
!! - The graph starts always at 0,0, Initial Z value in X array is ignored
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT1D_Z_SMOOTH_CRP6D(thispes,npoints,X,filename)
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables -------------------------------
   CLASS(CRP6D),INTENT(IN) :: thispes
   INTEGER, INTENT(IN) :: npoints
   CHARACTER(LEN=*),INTENT(IN) :: filename
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: X
   ! Local variables -----------------------------
   INTEGER :: inpoints, ndelta
   REAL(KIND=8) :: delta,v
   REAL(KIND=8) :: xmax, xmin, ymax, ymin 
   REAL(KIND=8), DIMENSION(6) :: r, dvdu
   INTEGER(KIND=4) :: i ! Counter
   CHARACTER(LEN=23),PARAMETER :: routinename = "PLOT1D_Z_SMOOTH_CRP6D: "
   ! HE HO ! LET'S GO ----------------------------
   SELECT CASE(npoints)
      CASE(: 1)
         WRITE(0,*) "PLOT1D_Z_CRP6D ERR: Less than 2 points"
         CALL EXIT(1)
      CASE DEFAULT
         ! do nothing
   END SELECT
   !
   xmin = thispes%wyckoffsite(1)%zrcut(1)%getfirstz()
   xmax = thispes%wyckoffsite(1)%zrcut(1)%getlastz()
   r(1:2)=x(1:2)
   r(4:6)=x(4:6)
   !
   inpoints=npoints-2
   ndelta=npoints-1
   delta=(xmax-xmin)/dfloat(ndelta)
   !
   OPEN(11,file=filename,status="replace")
   ! Initial value
   r(3)=xmin
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)  
   ! cycle for inpoints
   DO i=1, inpoints
      r(3)=xmin+DFLOAT(i)*delta
      CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      WRITE(11,*) r(3),v,dvdu(:)
   END DO
   ! Final value
   r(3) = xmax
   CALL thispes%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   WRITE(11,*) r(3),v,dvdu(:)
   WRITE(*,*) routinename, "file created ",filename
   CLOSE(11)
END SUBROUTINE PLOT1D_Z_SMOOTH_CRP6D
!#######################################################################
! SUBROUTINE: PLOT_XYMAP_CRP6D
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (X,Y) of the PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the file to print the output
!> @param[in] init_point - Initial position to start the scan (a.u. and radians). 6 DIMENSIONAL
!> @param[in] nxpoints - Number of points in X axis (auxiliar cartesian coordinates)
!> @param[in] nypoints - Number of points in Y axis (auxiliar cartesian coordinates)
!> @param[in] Lx - Length of X axis (a.u.)
!> @param[in] Ly - Length of Y axis (a.u.)
!
!> @warning
!! - Z,R,THETA,PHI parameters are taken from @b init_point
!
!> @author A.S. Muzas
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT_XYMAP_CRP6D(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(CRP6D),INTENT(IN) :: thispes
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: init_point ! Initial position to start the scan (in a.u. and radians)
   INTEGER,INTENT(IN) :: nxpoints, nypoints ! number of points in XY plane
   CHARACTER(LEN=*),INTENT(IN) :: filename ! filename
   REAL(KIND=8),INTENT(IN) :: Lx ! Length of X axis 
   REAL(KIND=8),INTENT(IN) :: Ly ! Length of X axis 
   ! Local variables
   REAL(KIND=8) :: xmin, ymin, xmax, ymax
   REAL(KIND=8),DIMENSION(6) :: r,dvdu
   REAL(KIND=8) :: xdelta, ydelta
   INTEGER :: xinpoints, nxdelta
   INTEGER :: yinpoints, nydelta
   INTEGER :: i, j ! counters
   REAL(KIND=8) :: v ! potential
   ! GABBA, GABBA HEY! ---------
   xmin = init_point(1)
   ymin = init_point(2)
   xmax = init_point(1)+Lx
   ymax = init_point(2)+Ly
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
   r(3:6)=init_point(3:6)
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
   END DO
   r(2) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   ! inpoints in XY
   DO i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      r(3:6) = init_point(3:6)
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
      DO j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
         WRITE(11,*) r(1:2),v,dvdu(:)
      END DO
      r(2) = ymax
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
   END DO
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   r(3:6) = init_point(3:6)
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   DO i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(1:2),v,dvdu(:)
   END DO
   r(2) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(1:2),v,dvdu(:)
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_XYMAP_CRP6D
!#######################################################################
! SUBROUTINE: PLOT_RZMAP_CRP6D
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (R,Z) of the PES. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the file to print the output
!> @param[in] init_point - Initial position to start the scan (a.u. and radians). 6 DIMENSIONAL
!> @param[in] nxpoints - Number of points in R axis
!> @param[in] nypoints - Number of points in Z axis
!> @param[in] Lx - Length of R axis (a.u.)
!> @param[in] Ly - Length of Z axis (a.u.)
!
!> @warning
!! - X,Y,THETA,PHI parameters are taken from @b init_point
!
!> @author A.S. Muzas
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT_RZMAP_CRP6D(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(CRP6D),INTENT(IN) :: thispes
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: init_point 
   INTEGER,INTENT(IN) :: nxpoints, nypoints 
   CHARACTER(LEN=*),INTENT(IN) :: filename 
   REAL(KIND=8),INTENT(IN) :: Lx
   REAL(KIND=8),INTENT(IN) :: Ly
   ! Local variables
   REAL(KIND=8) :: xmin, ymin, xmax, ymax
   REAL(KIND=8),DIMENSION(6) :: r,dvdu
   REAL(KIND=8) :: xdelta, ydelta
   INTEGER :: xinpoints, nxdelta
   INTEGER :: yinpoints, nydelta
   INTEGER :: i, j ! counters
   REAL(KIND=8) :: v ! potential
   ! GABBA, GABBA HEY! ---------
   xmin = init_point(4)
   ymin = init_point(3)
   xmax = init_point(4)+Lx
   ymax = init_point(3)+Ly
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
   r(4) = xmin
   r(3) = ymin
   r(1:2)=init_point(1:2)
   r(5:6)=init_point(5:6)
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),r(3),v,dvdu(:)
   DO i =1, yinpoints
      r(3) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(4),r(3),v,dvdu(:)
   END DO
   r(3) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),r(3),v,dvdu(:)
   ! inpoints in XY
   DO i = 1, xinpoints
      r(4) = xmin+DFLOAT(i)*xdelta
      r(3) = ymin
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(4),r(3),v,dvdu(:)
      DO j = 1, yinpoints
         r(3) = ymin + DFLOAT(j)*ydelta
         CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
         WRITE(11,*) r(4),r(3),v,dvdu(:)
      END DO
      r(3) = ymax
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(4),r(3),v,dvdu(:)
   END DO
   ! Last point in XY plane
   r(4) = xmax
   r(3) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),r(3),v,dvdu(:)
   DO i =1, yinpoints
      r(3) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
      WRITE(11,*) r(4),r(3),v,dvdu(:)
   END DO
   r(3) = ymax
   CALL thispes%GET_V_AND_DERIVS(r,v,dvdu)
   WRITE(11,*) r(4),r(3),v,dvdu(:)
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_RZMAP_CRP6D
!#######################################################################
! SUBROUTINE: PLOT_ATOMIC_INTERAC_RZ_CRP6D
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (R,Z) of the sub of atomic potentials. 
!
!> @param[in] thispes - CRP6D PES used
!> @param[in] filename - Name of the file to print the output
!> @param[in] init_point - Initial position to start the scan (a.u. and radians). 6 DIMENSIONAL
!> @param[in] nxpoints - Number of points in R axis
!> @param[in] nypoints - Number of points in Z axis
!> @param[in] Lx - Length of R axis (a.u.)
!> @param[in] Ly - Length of Z axis (a.u.)
!
!> @warning
!! - X,Y,THETA,PHI parameters are taken from @b init_point
!
!> @author A.S. Muzas
!> @date Apr/2014
!> @version 1.0
!----------------------------------------------------------------------
SUBROUTINE PLOT_ATOMIC_INTERAC_RZ_CRP6D(thispes,init_point,nxpoints,nypoints,Lx,Ly,filename)
   IMPLICIT NONE
   CLASS(CRP6D),INTENT(IN) :: thispes
   REAL(KIND=8),DIMENSION(6),INTENT(IN) :: init_point 
   INTEGER,INTENT(IN) :: nxpoints, nypoints 
   CHARACTER(LEN=*),INTENT(IN) :: filename 
   REAL(KIND=8),INTENT(IN) :: Lx
   REAL(KIND=8),INTENT(IN) :: Ly
   ! Local variables
   REAL(KIND=8) :: xmin, ymin, xmax, ymax
   REAL(KIND=8),DIMENSION(6) :: r
   REAL(KIND=8) :: xdelta, ydelta
   INTEGER :: xinpoints, nxdelta
   INTEGER :: yinpoints, nydelta
   REAL(KIND=8),DIMENSION(2) :: v
   REAL(KIND=8),DIMENSION(6) :: atomicx
   REAL(KIND=8),DIMENSION(6) :: dvdu
   INTEGER :: i, j ! counters
   ! GABBA, GABBA HEY! ---------
   xmin = init_point(4)
   ymin = init_point(3)
   xmax = init_point(4)+Lx
   ymax = init_point(3)+Ly
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
   r(4) = xmin
   r(3) = ymin
   r(1:2)=init_point(1:2)
   r(5:6)=init_point(5:6)
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
   WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   DO i =1, yinpoints
      r(3) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
      WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   END DO
   r(3) = ymax
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
   WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   ! inpoints in XY
   DO i = 1, xinpoints
      r(4) = xmin+DFLOAT(i)*xdelta
      r(3) = ymin
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
      WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
      DO j = 1, yinpoints
         r(3) = ymin + DFLOAT(j)*ydelta
         CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
         WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
      END DO
      r(3) = ymax
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
      WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   END DO
   ! Last point in XY plane
   r(4) = xmax
   r(3) = ymax
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
   WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   DO i =1, yinpoints
      r(3) = ymin + DFLOAT(i)*ydelta
      CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
      WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   END DO
   r(3) = ymax
   CALL thispes%GET_ATOMICPOT_AND_DERIVS(r,atomicx,v,dvdu)
   WRITE(11,*) r(4),r(3),sum(v),v(1),v(2)
   CLOSE(11)
   RETURN
END SUBROUTINE PLOT_ATOMIC_INTERAC_RZ_CRP6D
!###########################################################
!# FUNCTION: is_allowed_CRP6D 
!###########################################################
!> @brief 
!! Enquires if the CRP6D potential can be calculated at @b X
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
LOGICAL FUNCTION is_allowed_CRP6D(this,x) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(CRP6D),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(:),INTENT(IN) :: x
   ! Local variables
   REAL(KIND=8) :: zmin, rmin, rmax
   ! Run section
   zmin=this%wyckoffsite(1)%zrcut(1)%getfirstZ()
   rmin=this%wyckoffsite(1)%zrcut(1)%getfirstR()
   rmax=this%wyckoffsite(1)%zrcut(1)%getlastR()
   SELECT CASE(size(x)/=6)
      CASE(.TRUE.)
         WRITE(0,*) "is_allowed_CRP6D ERR: wrong number of dimensions"
         CALL EXIT(1)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   SELECT CASE(x(3)<zmin .OR. x(4)<rmin .OR. x(4)>rmax)
      CASE(.TRUE.)
         is_allowed_CRP6D=.FALSE.
      CASE(.FALSE.)
         is_allowed_CRP6D=.TRUE.
   END SELECT
   RETURN
END FUNCTION is_allowed_CRP6D

END MODULE CRP6D_MOD
