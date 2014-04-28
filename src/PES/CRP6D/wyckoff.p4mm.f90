!#########################################################
! MODULE: WYCKOFF_P4MM
!> @brief
!! Provides tools to interpolate through wyckoff sites belonging to
!! p4mm wallpaper symmetry
!##########################################################
MODULE WYCKOFF_P4MM_MOD
! Initial declarations
USE WYCKOFF_GENERIC_MOD
USE FOURIER1D_MOD
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
!! Interpolates all phicuts and thetacuts for a given value of r and z.
!! Gets the potential and derivatives for z,r,theta,phi
!
!> @param[in] x - Array which stands for Z,r,theta,phi
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_V_AND_DERIVS_WYCKOFFP4MM(this,x,v,dvdu)
   ! Initial declarations   
#ifdef DEBUG
   USE DEBUG_MOD
   USE UNITS_MOD
#endif
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(Wyckoffp4mm),INTENT(IN) :: this
   REAL(KIND=8),DIMENSION(4),INTENT(IN) ::x 
   REAL(KIND=8),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(4),INTENT(OUT) :: dvdu
   ! Local variables
   INTEGER(KIND=4) :: i,j,h ! counters
   REAL(KIND=8) :: dummy
   REAL(KIND=8) :: period_phi,period_theta
   REAL(KIND=8) :: z,r,theta,phi
   TYPE(Fourier1d),DIMENSION(:),ALLOCATABLE :: phicut
   TYPE(Fourier1d) :: thetacut
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: aux
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: f
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: dfdz
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: dfdr
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: dfdphi
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: philist
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: thetalist
   LOGICAL :: phi_is_even,theta_is_even
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
   h=0
   SELECT CASE(this%id)
      !----------------------------------------------------
      CASE("a" : "c") 
         ! Intrinsic periodicities
         SELECT CASE(this%id)
            CASE("a" : "b")
               period_phi=PI/2.D0
               phi_is_even=.TRUE.
               theta_is_even=.TRUE. ! due to inversion center
            CASE("c")
               period_phi=PI
               phi_is_even=.TRUE.
               theta_is_even=.TRUE. ! due to inversion center
            CASE DEFAULT
               WRITE(0,*) "GET_V_AND_DERIVS_WYCKOFFP4MM: Unexpected error with Wyckoff id"
               CALL EXIT(1)
         END SELECT
         period_theta=PI ! this stands always for a~c sites (due to inversion center)
         ! Prepare phi interpolation
         ALLOCATE(phicut(this%nphicuts))
         DO i = 1, this%nphicuts
            ALLOCATE(f(this%nphipoints(i)))
            ALLOCATE(dfdr(this%nphipoints(i)))
            ALLOCATE(dfdz(this%nphipoints(i)))
            ALLOCATE(philist(this%nphipoints(i)))
            ALLOCATE(aux(2,this%nphipoints(i)))
            DO j = 1, this%nphipoints(i) ! loop over number of zrcuts inside
               h=h+1 ! numbering of zrcuts
               f(j)=this%zrcut(h)%interrz%getvalue((/r,z/))
               dfdr(j)=this%zrcut(h)%interrz%getderivx((/r,z/))
               dfdz(j)=this%zrcut(h)%interrz%getderivy((/r,z/))
               philist(j)=this%zrcut(h)%phi
            END DO
            aux(1,:)=dfdz(:)
            aux(2,:)=dfdr(:)
            CALL phicut(i)%READ(philist,f)
            CALL phicut(i)%ADD_MOREFUNCS(aux)
            CALL phicut(i)%READ_EXTRA(period_phi,this%klistphi(i)%k,phi_is_even)
            CALL phicut(i)%INTERPOL()
#ifdef DEBUG
            CALL VERBOSE_WRITE(routinename,"New phicut:")
            CALL VERBOSE_WRITE(routinename,"At Phi: (deg)")
            ALLOCATE(beta(size(philist)))
            ALLOCATE(philistdeg(size(philist)))
            DO j = 1, size(philist)
               CALL beta(j)%READ(philist(j),"rad")
               CALL beta(j)%TO_DEG()
               philistdeg(j)=beta(j)%getvalue()
            END DO
            CALl VERBOSE_WRITE(routinename,philistdeg)
            CALL VERBOSE_WRITE(routinename,"At Phi: (rad)")
            CALl VERBOSE_WRITE(routinename,philist)
            CALL VERBOSE_WRITE(routinename,"f:")
            CALl VERBOSE_WRITE(routinename,f)
            CALL VERBOSE_WRITE(routinename,"dfdz:")
            CALL VERBOSE_WRITE(routinename,aux(1,:))
            CALL VERBOSE_WRITE(routinename,"dfdr:")
            CALL VERBOSE_WRITE(routinename,aux(2,:))
            WRITE(filename,'(I1,A1,I1,A1,A14)') this%mynumber,"-",i,"-","wyckoffphi.dat" 
            CALL phicut(i)%PLOT(100,filename)
            WRITE(filename,'(I1,A1,I1,A1,A14)') this%mynumber,"-",i,"-","wyckoffphi.raw" 
            CALL phicut(i)%PLOTDATA(filename)
            WRITE(filename,'(I1,A1,I1,A1,A14)') this%mynumber,"-",i,"-","wyckoffphi.cyc" 
            CALL phicut(i)%PLOTCYCLIC(300,filename)
            DEALLOCATE(beta)
            DEALLOCATE(philistdeg)
#endif
            DEALLOCATE(philist)
            DEALLOCATE(f)
            DEALLOCATE(dfdr)
            DEALLOCATE(dfdz)
            DEALLOCATE(aux)
         END DO
      !----------------------------------------------------
      CASE("f") ! in this case we should interpol respect to phi-pi/4
         phi_is_even=.TRUE. ! in phi-pi/4 is even
         theta_is_even=.FALSE.
         ! In this case, homonuclear case increases periodicity
         SELECT CASE(this%is_homonucl)
            CASE(.TRUE.)
               period_phi=PI
               period_theta=PI
            CASE(.FALSE.)
               period_phi=2.D0*PI
               period_theta=2.D0*PI
         END SELECT
         ALLOCATE(phicut(this%nphicuts))
         DO i = 1, this%nphicuts
            ALLOCATE(f(this%nphipoints(i)))
            ALLOCATE(dfdr(this%nphipoints(i)))
            ALLOCATE(dfdz(this%nphipoints(i)))
            ALLOCATE(philist(this%nphipoints(i)))
            ALLOCATE(aux(2,this%nphipoints(i)))
            DO j = 1, this%nphipoints(i) ! loop over number of zrcuts inside
               h=h+1 ! numbering of zrcuts
               f(j)=this%zrcut(h)%interrz%getvalue((/r,z/))
               dfdr(j)=this%zrcut(h)%interrz%getderivx((/r,z/))
               dfdz(j)=this%zrcut(h)%interrz%getderivy((/r,z/))
               philist(j)=this%zrcut(h)%phi-PI/4.D0
            END DO
            aux(1,:)=dfdz(:)
            aux(2,:)=dfdr(:)
            CALL phicut(i)%READ(philist,f)
            CALL phicut(i)%ADD_MOREFUNCS(aux)
            CALL phicut(i)%READ_EXTRA(period_phi,this%klistphi(i)%k,phi_is_even)
            CALL phicut(i)%INTERPOL()
#ifdef DEBUG
            CALL VERBOSE_WRITE(routinename,"New phicut:")
            CALL VERBOSE_WRITE(routinename,"At Phi: (deg)")
            ALLOCATE(beta(size(philist)))
            ALLOCATE(philistdeg(size(philist)))
            DO j = 1, size(philist)
               CALL beta(j)%READ(philist(j),"rad")
               CALL beta(j)%TO_DEG()
               philistdeg(j)=beta(j)%getvalue()
            END DO
            CALl VERBOSE_WRITE(routinename,philistdeg)
            CALL VERBOSE_WRITE(routinename,"At Phi: (rad)")
            CALl VERBOSE_WRITE(routinename,philist)
            CALL VERBOSE_WRITE(routinename,"f:")
            CALl VERBOSE_WRITE(routinename,f)
            CALL VERBOSE_WRITE(routinename,"dfdz:")
            CALL VERBOSE_WRITE(routinename,aux(1,:))
            CALL VERBOSE_WRITE(routinename,"dfdr:")
            CALL VERBOSE_WRITE(routinename,aux(2,:))
            WRITE(filename,'(I1,A1,I1,A1,A14)') this%mynumber,"-",i,"-","wyckoffphi.dat" 
            CALL phicut(i)%PLOT(100,filename)
            WRITE(filename,'(I1,A1,I1,A1,A14)') this%mynumber,"-",i,"-","wyckoffphi.raw" 
            CALL phicut(i)%PLOTDATA(filename)
            WRITE(filename,'(I1,A1,I1,A1,A14)') this%mynumber,"-",i,"-","wyckoffphi.cyc" 
            CALL phicut(i)%PLOTCYCLIC(300,filename)
            DEALLOCATE(beta)
            DEALLOCATE(philistdeg)
#endif
            DEALLOCATE(philist)
            DEALLOCATE(f)
            DEALLOCATE(dfdr)
            DEALLOCATE(dfdz)
            DEALLOCATE(aux)
         END DO
      !--------------------------------------------------------
      CASE DEFAULT
         WRITE(0,*) "INTERPOL_WYCKOFFP4MM ERR: Wrong wyckoff letter for p4mm symmetry or it is not yet implemented"
         CALL EXIT(1)
   END SELECT
   ! each phicut is related to a theta value
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
      SELECT CASE(this%id)
         CASE("a" : "c")
            CALL phicut(i)%GET_ALLFUNCS_AND_DERIVS(phi,aux(1,:),aux(2,:))
         CASE("f")
            ! correction to phi
            CALL phicut(i)%GET_ALLFUNCS_AND_DERIVS(phi-PI/4.D0,aux(1,:),aux(2,:))
         CASE DEFAULT
            WRITE(0,*) "INTERPOL_WYCKOFFP4MM ERR: Wrong wyckoff letter for p4mm symmetry or it is not yet implemented"
            CALL EXIT(1)
      END SELECT
      f(i)=aux(1,1)
      dfdz(i)=aux(1,2)
      dfdr(i)=aux(1,3)
      dfdphi(i)=aux(2,1)
      DEALLOCATE(aux)
   END DO
   CALL thetacut%READ(thetalist,f)
   CALL thetacut%READ_EXTRA(period_theta,this%klisttheta,theta_is_even)
   ALLOCATE(aux(3,this%nphicuts))
   aux(1,:)=dfdz(:)
   aux(2,:)=dfdr(:)
   aux(3,:)=dfdphi(:)
   CALL thetacut%ADD_MOREFUNCS(aux)
   CALL thetacut%INTERPOL()
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"New thetacut:")
   CALL VERBOSE_WRITE(routinename,"At Theta: (deg)")
   ALLOCATE(beta(size(thetalist)))
   ALLOCATE(philistdeg(size(thetalist)))
   DO j = 1, size(thetalist)
      CALL beta(j)%READ(thetalist(j),"rad")
      CALL beta(j)%TO_DEG()
      philistdeg(j)=beta(j)%getvalue()
   END DO
   CALl VERBOSE_WRITE(routinename,philistdeg)
   CALL VERBOSE_WRITE(routinename,"At Theta: (rad)")
   CALl VERBOSE_WRITE(routinename,thetalist)
   CALL VERBOSE_WRITE(routinename,"f:")
   CALl VERBOSE_WRITE(routinename,f)
   CALL VERBOSE_WRITE(routinename,"dfdz:")
   CALL VERBOSE_WRITE(routinename,dfdz)
   CALL VERBOSE_WRITE(routinename,"dfdr:")
   CALL VERBOSE_WRITE(routinename,dfdr)
   CALL VERBOSE_WRITE(routinename,"dfdphi:")
   CALL VERBOSE_WRITE(routinename,dfdphi)
   SELECT CASE(get_debugmode())
      CASE(.TRUE.)
         WRITE(filename,'(I1,A1,A16)') this%mynumber,"-","wyckofftheta.dat" 
         CALL thetacut%PLOT(100,filename)
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
