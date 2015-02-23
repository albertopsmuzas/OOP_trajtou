!#########################################################
! MODULE: WYCKOFF_P4MM
!> @brief
!! Provides tools to interpolate through wyckoff sites belonging to
!! p4mm wallpaper symmetry
!##########################################################
MODULE WYCKOFF_P4MM_MOD
! Initial declarations
USE WYCKOFF_GENERIC_MOD
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////
! TYPE: Wyckoffp4mm
!> @brief
!! Subclass of Wyckoffsitio for generic p4mm symmetry.
!! in p4mm symmetry does not matter
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 20/03/2014 
!> @version 1.0
!----------------------------------------------------------------
TYPE,EXTENDS(Wyckoffsitio) :: Wyckoffp4mm
   CONTAINS
      PROCEDURE,PUBLIC :: GET_V_AND_DERIVS => GET_V_AND_DERIVS_WYCKOFFP4MM
END TYPE Wyckoffp4mm
!/////////////////////////////////////////////////////////////////
CONTAINS
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_WYCKOFFP4MM 
!###########################################################
!> @brief
!! Symmetry adapted interpolation for p4mm wallpaper group.
!
!> @param[in] x - Array which stands for Z,r,theta,phi
!> @param[out] v - potential at X
!> @param[out] dvdu - Array with derivatives: dvdz, dvdr, dvdtheta dvdphi
!
!> @warning
!! - Only stands for a,b,c, and f p4mm Wyckoff sites.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014
!> @version 1.2
!-----------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_WYCKOFFP4MM(this,x,v,dvdu)
   ! Initial declarations   
#ifdef DEBUG
   USE DEBUG_MOD
   USE UNITS_MOD
#endif
   USE FOURIER1D_2_MOD
   USE FOURIER1D_4MM_MOD
   USE FOURIER1D_M_MOD
   USE FOURIER1D_M45_MOD
   USE FOURIER1D_MM2_MOD
   USE FOURIER1D_M45M1352_MOD
   USE FOURIER1D_E_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(Wyckoffp4mm),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(4),INTENT(IN) ::x 
   REAL(KIND=8),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(4),INTENT(OUT) :: dvdu
   ! Local variables
   INTEGER(KIND=4) :: i,j,h ! counters
   REAL(KIND=8) :: z,r,theta,phi
   CLASS(Fourier1d),DIMENSION(:),ALLOCATABLE :: phicut
   CLASS(Fourier1d),ALLOCATABLE :: thetacut
   REAL(KIND=8) :: aux_theta
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: aux
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: f
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: dfdz
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: dfdr
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: dfdphi
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: philist
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: thetalist
   CHARACTER(LEN=2) :: theta_irrep, phi_irrep
   REAL(KIND=8) :: phi_shift
   CHARACTER(LEN=30),PARAMETER :: routinename="GET_V_AND_DERIVS_WYCKOFFP4MM: "
#ifdef DEBUG
   TYPE(Angle),DIMENSION(:),ALLOCATABLE :: beta
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: philistdeg
   CHARACTER(LEN=18) :: filename
#endif
   ! Run section
   z=x(1)
   r=x(2)
   theta=x(3)
   phi=x(4)
   ! CONDITIONS: edit this part to include/edit symmetries
   SELECT CASE(this%id)
      CASE("a" : "b")
         ALLOCATE(Fourier1d_4mm::phicut(this%nphicuts))
         phi_irrep="A1"
         phi_shift=0.D0
         SELECT CASE(this%is_homonucl)
            CASE(.TRUE.)
               ALLOCATE(Fourier1d_mm2::thetacut)
               theta_irrep="A1"
            CASE(.FALSE.)
               ALLOCATE(Fourier1d_m::thetacut)
               theta_irrep="Ap"
         END SELECT
         
      CASE("c")
         ALLOCATE(Fourier1d_mm2::phicut(this%nphicuts))
         phi_irrep="A1"
         phi_shift=0.D0
         SELECT CASE(this%is_homonucl)
            CASE(.TRUE.)
               ALLOCATE(Fourier1d_mm2::thetacut)
               theta_irrep="A1"
            CASE(.FALSE.)
               ALLOCATE(Fourier1d_m::thetacut)
               theta_irrep="Ap"
         END SELECT

      CASE("f")
         ALLOCATE(Fourier1d_m45::phicut(this%nphicuts))
         !phi_irrep="Ap" ! should be decided later. There are special symmetries depending upon theta
         phi_shift=0.D0
         SELECT CASE(this%is_homonucl)
            CASE(.TRUE.)
               ALLOCATE(Fourier1d_2::thetacut)
               theta_irrep="A"
            CASE(.FALSE.)
               ALLOCATE(Fourier1d_e::thetacut)
               theta_irrep="A"
         END SELECT

      CASE DEFAULT
         WRITE(0,*) "GET_V_AND_DERIVS_WYCKOFFP4MM: Unexpected error with Wyckoff id"
         CALL EXIT(1)
   END SELECT
   ! PHI INTERPOLATION ----------------------------------------------------
   h=0 ! initialize h
   DO i = 1, this%nphicuts ! loop over specific phi cuts (each one for a different theta value) 
      ALLOCATE(f(this%nphipoints(i)))
      ALLOCATE(dfdr(this%nphipoints(i)))
      ALLOCATE(dfdz(this%nphipoints(i)))
      ALLOCATE(philist(this%nphipoints(i)))
      ALLOCATE(aux(2,this%nphipoints(i)))
      DO j = 1, this%nphipoints(i) ! loop over number of zrcuts inside
         h=h+1 ! numbering of zrcuts
         f(j)=this%zrcut(h)%interrz%getvalue((/r,z/)) ! storing potential at this site
         dfdr(j)=this%zrcut(h)%interrz%getderivx((/r,z/)) ! storing d/dr at this site
         dfdz(j)=this%zrcut(h)%interrz%getderivy((/r,z/)) ! storing d/dz at this site
         philist(j)=this%zrcut(h)%phi
         aux_theta=this%zrcut(h)%theta
      END DO
      SELECT CASE(this%id=="f" .AND. this%is_homonucl) ! "f" has special simmetry if theta=90 and we've a homonuclear molecule
         CASE(.TRUE.)
            SELECT CASE(aux_theta>=dacos(0.D0)-1.D-6 .AND. aux_theta<=dacos(0.D0)+1.D-6) ! check if theta is pi/2
               CASE(.TRUE.)
                  phi_irrep="A1" ! expanded symmetry m45 -> m45m1352
               CASE(.FALSE.)
                  phi_irrep="Ap"
            END SELECT
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      aux(1,:)=dfdz(:)
      aux(2,:)=dfdr(:)
      CALL phicut(i)%READ(philist,f)
      CALL phicut(i)%ADD_MOREFUNCS(aux)
      CALL phicut(i)%SET_IRREP(phi_irrep)
      CALL phicut(i)%SET_SHIFT(phi_shift)
      CALL phicut(i)%INTERPOL()
#ifdef DEBUG
      CALL DEBUG_WRITE(routinename,"NEW PHICUT")
      CALL DEBUG_WRITE(routinename,"For theta: ",this%zrcut(h)%theta)
      CALL DEBUG_WRITE(routinename,"Is homonuclear: ",this%is_homonucl)
      ALLOCATE(beta(size(philist)))
      ALLOCATE(philistdeg(size(philist)))
      DO j = 1, size(philist)
         CALL beta(j)%READ(philist(j),"rad")
         CALL beta(j)%TO_DEG()
         philistdeg(j)=beta(j)%getvalue()
      END DO
      CALL DEBUG_WRITE(routinename,"At Phi: (deg)",philistdeg)
      CALL DEBUG_WRITE(routinename,"At Phi: (rad)",philist)
      CALL DEBUG_WRITE(routinename,"Irrep: ",phi_irrep)
      CALL DEBUG_WRITE(routinename,"Klist: ",phicut(i)%getklist())
      CALL DEBUG_WRITE(routinename,"f: ",f)
      CALL DEBUG_WRITE(routinename,"dfdz: ",aux(1,:))
      CALL DEBUG_WRITE(routinename,"dfdr: ",aux(2,:))
      SELECT CASE(get_debugmode())
         CASE(.TRUE.)
            WRITE(filename,'(I1,A1,I1,A1,A14)') this%mynumber,"-",i,"-","wyckoffphi.raw" 
            CALL phicut(i)%PLOTDATA(filename)
            WRITE(filename,'(I1,A1,I1,A1,A14)') this%mynumber,"-",i,"-","wyckoffphi.cyc" 
            CALL phicut(i)%PLOTCYCLIC_ALL(300,filename)
         CASE(.FALSE.)
            ! do nothing
      END SELECT
      DEALLOCATE(beta)
      DEALLOCATE(philistdeg)
#endif
      DEALLOCATE(philist)
      DEALLOCATE(f)
      DEALLOCATE(dfdr)
      DEALLOCATE(dfdz)
      DEALLOCATE(aux)
   END DO
   ! THETA INTERPOLATION --------------------------
   ALLOCATE(thetalist(this%nphicuts))
   ALLOCATE(f(this%nphicuts))
   ALLOCATE(dfdz(this%nphicuts))
   ALLOCATE(dfdr(this%nphicuts))
   ALLOCATE(dfdphi(this%nphicuts))
   h=0 ! reboot h
   DO i = 1, this%nphicuts
      DO j = 1, this%nphipoints(i)
         h=h+1 
         thetalist(i)=this%zrcut(h)%theta 
      END DO
      ALLOCATE(aux(2,3))
      CALL phicut(i)%GET_ALLFUNCS_AND_DERIVS(phi,aux(1,:),aux(2,:))
      f(i)=aux(1,1)
      dfdz(i)=aux(1,2)
      dfdr(i)=aux(1,3)
      dfdphi(i)=aux(2,1)
      DEALLOCATE(aux)
   END DO
   ALLOCATE(aux(3,this%nphicuts))
   aux(1,:)=dfdz(:)
   aux(2,:)=dfdr(:)
   aux(3,:)=dfdphi(:)
   CALL thetacut%READ(thetalist,f)
   CALL thetacut%ADD_MOREFUNCS(aux)
   CALL thetacut%SET_IRREP(theta_irrep)
   CALL thetacut%INTERPOL()
#ifdef DEBUG
   CALL DEBUG_WRITE(routinename,"NEW THETACUT")
   ALLOCATE(beta(size(thetalist)))
   ALLOCATE(philistdeg(size(thetalist)))
   DO j = 1, size(thetalist)
      CALL beta(j)%READ(thetalist(j),"rad")
      CALL beta(j)%TO_DEG()
      philistdeg(j)=beta(j)%getvalue()
   END DO
   CALL DEBUG_WRITE(routinename,"At Theta: (deg) ",philistdeg)
   CALL DEBUG_WRITE(routinename,"At Theta: (rad) ",thetalist)
   CALL DEBUG_WRITE(routinename,"Klist:          ",thetacut%getklist())
   CALL DEBUG_WRITE(routinename,"f:              ",f)
   CALL DEBUG_WRITE(routinename,"dfdz:           ",dfdz)
   CALL DEBUG_WRITE(routinename,"dfdr:           ",dfdr)
   CALL DEBUG_WRITE(routinename,"dfdphi:         ",dfdphi)
   SELECT CASE(get_debugmode())
      CASE(.TRUE.)
         WRITE(filename,'(I1,A1,A16)') this%mynumber,"-","wyckofftheta.raw" 
         CALL thetacut%PLOTDATA(filename)
         WRITE(filename,'(I1,A1,A16)') this%mynumber,"-","wyckofftheta.cyc" 
         CALL thetacut%PLOTCYCLIC(300,filename)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   DEALLOCATE(beta)
   DEALLOCATE(philistdeg)
#endif
   DEALLOCATE(f)
   DEALLOCATE(dfdr)
   DEALLOCATE(dfdz)
   DEALLOCATE(dfdphi)
   DEALLOCATE(aux)
   DEALLOCATE(thetalist)
   ALLOCATE(aux(2,4))
   CALL thetacut%GET_ALLFUNCS_AND_DERIVS(theta,aux(1,:),aux(2,:))
   v=aux(1,1) ! value of the potential
   dvdu(1)=aux(1,2) ! dvdz
   dvdu(2)=aux(1,3) ! dvdr
   dvdu(3)=aux(2,1) ! dvdtheta
   dvdu(4)=aux(1,4) ! dvdphi
   DEALLOCATE(aux)
   RETURN
END SUBROUTINE GET_V_AND_DERIVS_WYCKOFFP4MM
END MODULE WYCKOFF_P4MM_MOD
