!#########################################################
! MODULE: FOURIER_P4MM_MOD
!> @brief
!! Provides tools to genererate a symmetry adapted fourier
!! interpolation
!##########################################################
MODULE FOURIER_P4MM_MOD
!use other modules?
USE FOURIER2D_MOD
IMPLICIT NONE
TYPE,EXTENDS(Fourier2d) :: Fourierp4mm
   PRIVATE
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: coeff
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE :: termmap
   CONTAINS
      PROCEDURE,PUBLIC :: INTERPOL => INTERPOL_FOURIERP4MM
      PROCEDURE,PUBLIC :: GET_F_AND_DERIVS => GET_F_AND_DERIVS_FOURIERP4MM
      PROCEDURE,PUBLIC :: SET_TERMMAP => SET_TERMMAP_FOURIERP4MM
END TYPE Fourierp4mm
! variables and types, body
CONTAINS
!###########################################################
!# FUNCTION: termfoup4mm
!###########################################################
!
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfoup4mm(id,surf,k,r) 
   ! Initial declarations   
   USE SURFACE_MOD
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: id
   TYPE(Surface),INTENT(IN) :: surf
   INTEGER(KIND=4),DIMENSION(2),INTENT(IN) :: k
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8) :: g
   ! Run section
   g=2.D0*PI/surf%norm_s1
   SELECT CASE(id)
      CASE(0)
         termfoup4mm=1.D0
      CASE(1)
         termfoup4mm=dcos(g*dfloat(k(1))*r(1))+dcos(g*dfloat(k(1))*r(2))
      CASE(2)
         termfoup4mm=dcos(g*dfloat(k(1))*r(1))*dcos(g*dfloat(k(1))*r(2))
      CASE(3)
         termfoup4mm=dcos(g*dfloat(k(1))*r(1))*dcos(g*dfloat(k(2))*r(2))+dcos(g*dfloat(k(2))*r(1))*dcos(g*dfloat(k(1))*r(2))
      CASE DEFAULT
         WRITE(0,*) "termfourier ERR: Incorrect fourier term id: ", id
         CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfoup4mm
!###########################################################
!# FUNCTION: termfoup4mm_dx 
!###########################################################
!> @brief 
!! ! type brief explanation
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfoup4mm_dx(id,surf,k,r) 
   ! Initial declarations   
   USE SURFACE_MOD
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: id
   TYPE(Surface),INTENT(IN) :: surf
   INTEGER(KIND=4),DIMENSION(2),INTENT(IN) :: k
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8):: g
   ! Run section
   g=2.D0*PI/surf%norm_s1
   SELECT CASE(id)
      CASE(0)
         termfoup4mm_dx=0.D0
      CASE(1)
         termfoup4mm_dx=-g*k(1)*dsin(g*dfloat(k(1))*r(1))
      CASE(2)
         termfoup4mm_dx=-g*k(1)*dsin(g*dfloat(k(1))*r(1))*dcos(g*dfloat(k(1))*r(2))
      CASE(3)
         termfoup4mm_dx=-g*k(1)*dsin(g*dfloat(k(1))*r(1))*dcos(g*dfloat(k(2))*r(2))-&
            g*k(2)*dsin(g*dfloat(k(2))*r(1))*dcos(g*dfloat(k(1))*r(2))
      CASE DEFAULT
         WRITE(0,*) "termfourier ERR: Incorrect fourier term id: ", id
         CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfoup4mm_dx
!###########################################################
!# FUNCTION: termfoup4mm_dy
!###########################################################
!> @brief 
!! ! type brief explanation
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
REAL(KIND=8) FUNCTION termfoup4mm_dy(id,surf,k,r) 
   ! Initial declarations   
   USE SURFACE_MOD
   USE CONSTANTS_MOD
   IMPLICIT NONE
   ! I/O variables
   INTEGER(KIND=4),INTENT(IN) :: id
   TYPE(Surface),INTENT(IN) :: surf
   INTEGER(KIND=4),DIMENSION(2),INTENT(IN) :: k
   REAL(KIND=8),DIMENSION(2),INTENT(IN) :: r
   ! Local variables
   REAL(KIND=8):: g
   ! Run section
   g=2.D0*PI/surf%norm_s1
   SELECT CASE(id)
      CASE(0)
         termfoup4mm_dy=0.D0
      CASE(1)
         termfoup4mm_dy=-g*k(1)*dsin(g*dfloat(k(1))*r(2))
      CASE(2)
         termfoup4mm_dy=-g*k(1)*dcos(g*dfloat(k(1))*r(1))*dsin(g*dfloat(k(1))*r(2))
      CASE(3)
         termfoup4mm_dy=-g*k(2)*dcos(g*dfloat(k(1))*r(1))*dsin(g*dfloat(k(2))*r(2))-&
            g*k(1)*dcos(g*dfloat(k(2))*r(1))*dsin(g*dfloat(k(1))*r(2))
      CASE DEFAULT
         WRITE(0,*) "termfourier ERR: Incorrect fourier term id: ", id
         CALL EXIT(1)
   END SELECT
   RETURN
END FUNCTION termfoup4mm_dy
!###########################################################
!# SUBROUTINE: INTERPOL_FOURIERP4MM 
!###########################################################
!> @brief
!! Interpols with a symmetry adapted fourier series for p4mm
!! wallpaper group
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 28/Mar/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INTERPOL_FOURIERP4MM(this,surf,filename)
   ! Initial declarations   
   USE SURFACE_MOD
   USE MATHS_MOD
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourierp4mm),INTENT(INOUT) :: this
   TYPE(Surface),INTENT(IN) :: surf
   CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
   ! Local variables
   INTEGER(KIND=4) :: i,j !counters
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: tmtrx, inv_tmtrx
   ! Run section
   ALLOCATE(tmtrx(this%n,this%n))
   ALLOCATE(inv_tmtrx(this%n,this%n))
   ALLOCATE(this%coeff(this%n,size(this%f(:,1))))
   CALL this%SET_TERMMAP()
   DO i = 1, this%n ! loop over points
      DO j = 1, this%n ! loop over terms (one for each k point)
         tmtrx(i,j)=termfoup4mm(this%termmap(j),surf,this%klist(j,:),this%xy(i,:)) 
      END DO
   END DO
   CALL INV_MTRX(this%n,tmtrx,inv_tmtrx)
   DO i = 1, size(this%f(:,1))
      this%coeff(i,:)=matmul(inv_tmtrx,this%f(i,:))
   END DO
   SELECT CASE(present(filename)) ! Check if we want to print coefficients
      CASE(.TRUE.)
         OPEN (10,FILE=filename,STATUS="replace",ACTION="write")
         DO i = 1, size(this%f(:,1))
            WRITE(10,*) this%coeff(i,:)
         END DO
         CLOSE(10)
      CASE(.FALSE.)
         ! do nothing
   END SELECT
   RETURN
END SUBROUTINE INTERPOL_FOURIERP4MM
!###########################################################
!# SUBROUTINE: SET_TERMMAP_FOURIERP4MM 
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
SUBROUTINE SET_TERMMAP_FOURIERP4MM(this)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourierp4mm),INTENT(INOUT)::this
   ! Local variables
   INTEGER(KIND=4) :: i !counters
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
      !
      SELECT CASE(this%klist(i,1)>0 .AND. this%klist(i,2)==this%klist(i,1))
         CASE(.TRUE.)
            this%termmap(i)=2
            CYCLE
         CASE DEFAULT
            ! do nothing
      END SELECT
      !
      SELECT CASE(this%klist(i,2)<this%klist(i,1) .AND. this%klist(i,2)>0)
         CASE(.TRUE.)
            this%termmap(i)=3
            CYCLE
         CASE DEFAULT
            WRITE(0,*) "INTERPOL_FOURIER4MM ERR: Kpoint list has a bad item at position ",i
            CALL EXIT(1)
      END SELECT
   END DO
   RETURN
END SUBROUTINE SET_TERMMAP_FOURIERP4MM
!###########################################################
!# SUBROUTINE: GET_F_AND_DERIVS_FOURIERP4MM 
!###########################################################
!> @brief
!! Gets all values for functions at a given point r
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date ! type a date
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE GET_F_AND_DERIVS_FOURIERP4MM(this,surf,r,v,dvdu)
   USE SURFACE_MOD
   USE CONSTANTS_MOD
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourierp4mm),INTENT(IN):: this
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
   DO i = 1, size(this%f(:,1))
      DO j = 1, this%n
         terms(j)=termfoup4mm(this%termmap(j),surf,this%klist(j,:),r)
         terms_dx(j)=termfoup4mm_dx(this%termmap(j),surf,this%klist(j,:),r)
         terms_dy(j)=termfoup4mm_dy(this%termmap(j),surf,this%klist(j,:),r)
      END DO
      v(i)=dot_product(terms,this%coeff(:,i))
      dvdu(i,1)=dot_product(terms_dx,this%coeff(:,i))
      dvdu(i,2)=dot_product(terms_dy,this%coeff(:,i))
   END DO
   RETURN
END SUBROUTINE GET_F_AND_DERIVS_FOURIERP4MM
END MODULE FOURIER_P4MM_MOD 
