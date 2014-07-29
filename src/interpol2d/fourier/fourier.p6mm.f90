!#########################################################
! MODULE: FOURIER_P6MM_MOD
!> @brief
!! Provides tools to genererate a symmetry adapted fourier
!! interpolation
!##########################################################
MODULE FOURIER_P6MM_MOD
!use other modules?
USE FOURIER2D_MOD
IMPLICIT NONE
TYPE,EXTENDS(Fourier2d) :: FourierP6MM
   PRIVATE
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: coeff
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE :: termmap
   CONTAINS
      PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_FOURIERP6MM
      PROCEDURE,PUBLIC :: GET_F_AND_DERIVS => GET_F_AND_DERIVS_FOURIERP6MM
      PROCEDURE,PUBLIC :: SET_TERMMAP => SET_TERMMAP_FOURIERP6MM
END TYPE FourierP6MM
! variables and types, body
CONTAINS
!###########################################################
!# FUNCTION: termfouP6MM
!###########################################################
!
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfouP6MM(id,surf,k,r) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: id
   TYPE(Surface),INTENT(IN) :: surf
   INTEGER(KIND=4),DIMENSION(2),INTENT(IN) :: k
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8) :: g
    INTEGER(KIND=4) ::i ! counters
   ! Run section
   g=2.D0*PI/surf%norm_s1
   SELECT CASE(id)
      CASE(0)
         termfouP6MM=1.D0
      CASE(1)
         termfouP6MM=2.D0*dcos(g*dfloat(k(1))*r(1))*dcos(g*dfloat(k(1))*r(2)/dsqrt(3.D0))&
            +dcos(2.D0*g*dfloat(k(1))*r(2)/dsqrt(3.D0))
      CASE(2)
         termfouP6MM=2.D0*dcos(g*dfloat(k(1))*r(1)/2.D0)*dcos(g*dfloat(k(1))*r(2)*dsqrt(3.D0)/2.D0)&
            +dcos(g*dfloat(k(1))*r(1))
      CASE(3)
         termfouP6MM=dcos(g*dfloat(k(1))*r(1))*dcos(g*dfloat(k(1)-2*k(2))*r(2)/dsqrt(3.D0))&
            +dcos(g*dfloat(k(2))*r(1))*dcos(g*dfloat(k(1)-2*k(2))*r(2)/dsqrt(3.D0))&
            +dcos(g*dfloat(k(1)-k(2))*r(1))*dcos(g*dfloat(k(1)+k(2))*r(2)/dsqrt(3.D0))
      CASE DEFAULT
         WRITE(0,*) "termfouP6MM ERR: Incorrect fourier term id: ", id
         CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfouP6MM
!###########################################################
!# FUNCTION: termfouP6MM_dx 
!###########################################################
!> @brief 
!! ! type brief explanation
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfouP6MM_dx(id,surf,k,r) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: id
   TYPE(Surface),INTENT(IN) :: surf
   INTEGER(KIND=4),DIMENSION(2),INTENT(IN) :: k
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8):: g
   INTEGER(KIND=4) :: i ! counter
   ! Run section
   g=2.D0*PI/surf%norm_s1
   SELECT CASE(id)
      CASE(0)
         termfouP6MM_dx=0.D0
      CASE(1)
         termfouP6MM_dx=-2.D0*g*dfloat(k(1))*dsin(g*dfloat(k(1))*r(1))*dcos(g*dfloat(k(1))*r(2)/dsqrt(3.D0))
      CASE(2)
         termfouP6MM_dx=-g*dfloat(k(1))*dsin(g*dfloat(k(1))*r(1)/2.D0)*dcos(g*dfloat(k(1))*r(2)*dsqrt(3.D0)/2.D0)&
            -g*dfloat(k(1))*dsin(g*dfloat(k(1))*r(1))
      CASE(3)
         termfouP6MM_dx=-g*dfloat(k(1))*dsin(g*dfloat(k(1))*r(1))*dcos(g*dfloat(k(1)-2*k(2))*r(2)/dsqrt(3.D0))&
            -g*dfloat(k(2))*dsin(g*dfloat(k(2))*r(1))*dcos(g*dfloat(k(1)-2*k(2))*r(2)/dsqrt(3.D0))&
            -g*dfloat(k(1)-k(2))*dsin(g*dfloat(k(1)-k(2))*r(1))*dcos(g*dfloat(k(1)+k(2))*r(2)/dsqrt(3.D0))
      CASE DEFAULT
         WRITE(0,*) "termfouP6MM_dx ERR: Incorrect fourier term id: ", id
         CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfouP6MM_dx
!###########################################################
!# FUNCTION: termfouP6MM_dy
!###########################################################
!> @brief 
!! ! type brief explanation
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfouP6MM_dy(id,surf,k,r) 
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: id
   TYPE(Surface),INTENT(IN) :: surf
   INTEGER(KIND=4),DIMENSION(2),INTENT(IN) :: k
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8):: g,a
   INTEGER(KIND=4) :: i ! counters
   ! Run section
   g=2.D0*PI/surf%norm_s1
   SELECT CASE(id)
      CASE(0)
         termfouP6MM_dy=0.D0
      CASE(1)
         a=2.D0*g*dfloat(k(1))/dsqrt(3.D0)
         termfouP6MM_dy=-a*dcos(g*dfloat(k(1))*r(1))*dsin(g*dfloat(k(1))*r(2)/dsqrt(3.D0))&
            -a*dsin(2.D0*g*dfloat(k(1))*r(2)/dsqrt(3.D0))
      CASE(2)
         a=g*dfloat(k(1))/dsqrt(3.D0)
         termfouP6MM_dy=-a*dcos(g*dfloat(k(1))*r(1)/2.D0)*dsin(g*dfloat(k(1))*r(2)*dsqrt(3.D0)/2.D0)
      CASE(3)
         a=g/dsqrt(3.D0)
         termfouP6MM_dy=-a*dfloat(k(1)-2*k(2))*dcos(g*dfloat(k(1))*r(1))*dsin(g*dfloat(k(1)-2*k(2))*r(2)/dsqrt(3.D0))&
            -a*dfloat(k(1)-2*k(2))*dcos(g*dfloat(k(2))*r(1))*dsin(g*dfloat(k(1)-2*k(2))*r(2)/dsqrt(3.D0))&
            -a*dfloat(k(1)+k(2))*dcos(g*dfloat(k(1)-k(2))*r(1))*dsin(g*dfloat(k(1)+k(2))*r(2)/dsqrt(3.D0))
      CASE DEFAULT
         WRITE(0,*) "termfouP6MM ERR: Incorrect fourier term id: ", id
         CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfouP6MM_dy
!###########################################################
!# SUBROUTINE: INTERPOL_FOURIERP6MM 
!###########################################################
!> @brief
!! Interpols with a symmetry adapted fourier series for P6MM
!! wallpaper group
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 28/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INTERPOL_FOURIERP6MM(this,surf,filename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(FourierP6MM),INTENT(INOUT) :: this
   TYPE(Surface),INTENT(IN) :: surf
   CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
   ! Local variables
   INTEGER(KIND=4) :: i,j !counters
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: tmtrx, inv_tmtrx
   ! Run section
   ALLOCATE(tmtrx(this%n,this%n))
   ALLOCATE(inv_tmtrx(this%n,this%n))
   ALLOCATE(this%coeff(this%n,this%nfunc))
   CALL this%SET_TERMMAP()
   DO i = 1, this%n ! loop over points
      DO j = 1, this%n ! loop over terms (one for each k point)
         tmtrx(i,j)=termfouP6MM(this%termmap(j),surf,this%klist(j,:),this%xy(i,:)) 
      END DO
   END DO
   CALL INV_MTRX(this%n,tmtrx,inv_tmtrx)
   DO i = 1,this%nfunc ! looop over functions
      this%coeff(:,i)=matmul(inv_tmtrx,this%f(i,:))
   END DO
   SELECT CASE(present(filename)) ! Check if we want to print coefficients
      CASE(.TRUE.)
         OPEN (10,FILE=filename,STATUS="replace",ACTION="write")
         DO i = 1,this%nfunc
            WRITE(10,*) this%coeff(:,i)
         END DO
         CLOSE(10)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE INTERPOL_FOURIERP6MM
!###########################################################
!# SUBROUTINE: SET_TERMMAP_FOURIERP6MM 
!###########################################################
!> @brief
!! Sets the map that identifies the kind of terms that should
!! be used during the fourier interpolation. These terms depend
!! on the kpoints chosen.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 28/Mar/2014 
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE SET_TERMMAP_FOURIERP6MM(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(FourierP6MM),INTENT(INOUT)::this
   ! Local variables
   INTEGER(KIND=4) :: i !counters
   CHARACTER(LEN=25),PARAMETER :: routinename="SET_TERMMAP_FOURIERP6MM: "
   ! Run section
   ALLOCATE(this%termmap(this%n))
   DO i = 1, this%n
      SELECT CASE(this%klist(i,1)==0 .AND. this%klist(i,2)==0)
         CASE(.TRUE.)
            this%termmap(i)=0
            CYCLE
         CASE DEFAULT
            ! do nothing
      END SELECT
      !
      SELECT CASE(this%klist(i,1)>0 .AND. this%klist(i,2)==0)
         CASE(.TRUE.)
            this%termmap(i)=1
            CYCLE
         CASE DEFAULT
            ! do nothing
      END SELECT
      SELECT CASE(this%klist(i,1)>0 .AND. mod(this%klist(i,1),2)==0 .AND. this%klist(i,2)==this%klist(i,1)/2 )
         CASE(.TRUE.)
            this%termmap(i)=2
            CYCLE
         CASE DEFAULT
            ! do nothing
      END SELECT
      !
      SELECT CASE(this%klist(i,1)>0 .AND. this%klist(i,2)>0 .AND. this%klist(i,2)<this%klist(i,1)/2 )
         CASE(.TRUE.)
            this%termmap(i)=3
            CYCLE
         CASE DEFAULT
            WRITE(0,*) "INTERPOL_FOURIERP6MM ERR: Kpoint list has a bad item at position ",i
            CALL EXIT(1)
      END SELECT
      !
   END DO
   RETURN
END SUBROUTINE SET_TERMMAP_FOURIERP6MM
!###########################################################
!# SUBROUTINE: GET_F_AND_DERIVS_FOURIERP6MM 
!###########################################################
!> @brief
!! Gets all values for functions at a given point r
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_F_AND_DERIVS_FOURIERP6MM(this,surf,r,v,dvdu)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(FourierP6MM),INTENT(IN):: this
   TYPE(Surface),INTENT(IN):: surf
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: r
   REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: v
   REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: dvdu
   ! Local variables
   INTEGER(KIND=4) :: i,j ! counters
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: terms,terms_dx,terms_dy
   ! Run section
   ALLOCATE(terms(this%n))
   ALLOCATE(terms_dx(this%n))
   ALLOCATE(terms_dy(this%n))
   SELECT CASE(size(v) == this%nfunc .AND. size(dvdu(:,1)) == this%nfunc .AND.&
      size(dvdu(1,:)) == 2)
      CASE(.FALSE.)
         WRITE(0,*) "GET_F_AND_DERIVS_FOURIERP6MM ERR: size mismatch in v and stored values of f"
         WRITE(0,*) "GET_F_AND_DERIVS_FOURIERP6MM ERR: nfunc: ", this%nfunc
         WRITE(0,*) "GET_F_AND_DERIVS_FOURIERP6MM ERR: size(v): ",size(v)
         WRITE(0,*) "GET_F_AND_DERIVS_FOURIERP6MM ERR: size(vdvdu): ",size(dvdu(:,1))
         CALL EXIT(1)
      CASE(.TRUE.)
         ! do nothing
   END SELECT
   DO i = 1, this%nfunc
      DO j = 1, this%n
         terms(j)=termfouP6MM(this%termmap(j),surf,this%klist(j,:),r)
         terms_dx(j)=termfouP6MM_dx(this%termmap(j),surf,this%klist(j,:),r)
         terms_dy(j)=termfouP6MM_dy(this%termmap(j),surf,this%klist(j,:),r)
      END DO
      v(i)=dot_product(terms,this%coeff(:,i))
      dvdu(i,1)=dot_product(terms_dx,this%coeff(:,i))
      dvdu(i,2)=dot_product(terms_dy,this%coeff(:,i))
   END DO
   DEALLOCATE(terms)
   DEALLOCATE(terms_dx)
   DEALLOCATE(terms_dy)
   RETURN
END SUBROUTINE GET_F_AND_DERIVS_FOURIERP6MM
END MODULE FOURIER_P6MM_MOD 
