!#########################################################
! MODULE: FOURIER_P4MM_MOD
!> @brief
!! Provides tools to genererate a symmetry adapted fourier
!! interpolation
!##########################################################
MODULE FOURIER_P4MM_MOD
use FOURIER2D_MOD
use UNITS_MOD, only: pi
use MATHS_MOD, only: INV_MTRX
#ifdef DEBUG
use DEBUG_MOD, only: VERBOSE_WRITE, DEBUG_WRITE
#endif
IMPLICIT NONE

TYPE,EXTENDS(Fourier2d):: Fourierp4mm
   PRIVATE
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: coeff
   CONTAINS
      PROCEDURE,PUBLIC:: INTERPOL => INTERPOL_FOURIERP4MM
      PROCEDURE,PUBLIC:: GET_F_AND_DERIVS => GET_F_AND_DERIVS_FOURIERP4MM
END TYPE Fourierp4mm
! variables and types, body
CONTAINS
!###########################################################
!# FUNCTION: termfoup4mm
!###########################################################
!
!-----------------------------------------------------------
function termfoup4mm(surf,k,parity,irrep,r) result(answer)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(surface),intent(in):: surf
   integer(kind=4),dimension(2),intent(in):: k
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   real(kind=8),dimension(2),intent(in):: r
   ! Dummy out variable
   real(kind=8):: answer
   ! Local variables
   real(kind=8):: g
   ! Paramters
   character(len=*),parameter:: routinename='termfoup4mm: '
   ! Run section -------------------------------------------
   g=2.D0*PI/surf%norm_s1
   select case( irrep )
   case('A1')
      select case( parity )
      case('+')
         answer=dcos( g*dfloat(k(1))*r(1) )*dcos( g*dfloat(k(2))*r(2) )+&
                dcos( g*dfloat(k(2))*r(1) )*dcos( g*dfloat(k(1))*r(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select
   case default
      write(0,*) 'ERR '//routinename//' irrep not implemented: "'//irrep//'"'
      call exit(1)
   end select
   return
end function termfoup4mm
!###########################################################
!# FUNCTION: termfoup4mm_dx 
!###########################################################
!> @brief 
!! ! type brief explanation
!-----------------------------------------------------------
function termfoup4mm_dx(surf,k,parity,irrep,r) result(answer)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(surface),intent(in):: surf
   integer(kind=4),dimension(2),intent(in):: k
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   real(kind=8),dimension(2),intent(in):: r
   ! Dummy out variable
   real(kind=8):: answer
   ! Local variables
   real(kind=8):: g
   ! Paramters
   character(len=*),parameter:: routinename='termfoup4mm: '
   ! Run section -------------------------------------------
   g=2.D0*PI/surf%norm_s1
   select case( irrep )
   case('A1')
      select case( parity )
      case('+')
         answer=-g*( dfloat(k(1))*dsin( g*dfloat(k(1))*r(1) )*dcos( g*dfloat(k(2))*r(2) )+&
                     dfloat(k(2))*dsin( g*dfloat(k(2))*r(1) )*dcos( g*dfloat(k(1))*r(2) ) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select
   case default
      write(0,*) 'ERR '//routinename//' irrep not implemented: "'//irrep//'"'
      call exit(1)
   end select
   return
end function termfoup4mm_dx
!###########################################################
!# FUNCTION: termfoup4mm_dy
!###########################################################
!> @brief 
!! ! type brief explanation
!-----------------------------------------------------------
function termfoup4mm_dy(surf,k,parity,irrep,r) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(surface),intent(in):: surf
   integer(kind=4),dimension(2),intent(in):: k
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   real(kind=8),dimension(2),intent(in):: r
   ! Dummy out variable
   real(kind=8):: answer
   ! Local variables
   real(kind=8):: g
   ! Paramters
   character(len=*),parameter:: routinename='termfoup4mm: '
   ! Run section -------------------------------------------
   g=2.D0*PI/surf%norm_s1
   select case( irrep )
   case('A1')
      select case( parity )
      case('+')
         answer=-g*( dfloat(k(2))*dcos( g*dfloat(k(1))*r(1) )*dsin( g*dfloat(k(2))*r(2) )+&
                     dfloat(k(1))*dcos( g*dfloat(k(2))*r(1) )*dsin( g*dfloat(k(1))*r(2) ) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select
   case default
      write(0,*) 'ERR '//routinename//' irrep not implemented: "'//irrep//'"'
      call exit(1)
   end select
   return
end function termfoup4mm_dy
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
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourierp4mm),INTENT(INOUT) :: this
   class(Surface),INTENT(IN) :: surf
   CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
   ! Local variables
   INTEGER(KIND=4) :: i,j !counters
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: tmtrx, inv_tmtrx
   ! Run section
   ALLOCATE(tmtrx(this%n,this%n))
   ALLOCATE(inv_tmtrx(this%n,this%n))
   ALLOCATE(this%coeff(this%n,this%nfunc))
   DO i = 1, this%n ! loop over points
      DO j = 1, this%n ! loop over terms (one for each k point)
         tmtrx(i,j)=termfoup4mm( surf=surf,k=this%kList(j,:),parity=this%parityList(j),irrep=this%irrepList(j),r=this%xy(i,:) )
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
END SUBROUTINE INTERPOL_FOURIERP4MM
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
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   class(Fourierp4mm),INTENT(IN):: this
   class(Surface),INTENT(IN):: surf
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
   SELECT CASE( size(v) == this%nfunc .AND. size(dvdu(:,1)) == this%nfunc .AND.&
                size(dvdu(1,:)) == 2)
      CASE(.FALSE.)
         WRITE(0,*) "GET_F_AND_DERIVS_FOURIERP4MM ERR: size mismatch in v and stored values of f"
         WRITE(0,*) "GET_F_AND_DERIVS_FOURIERP4MM ERR: nfunc: ", this%nfunc
         WRITE(0,*) "GET_F_AND_DERIVS_FOURIERP4MM ERR: size(v): ",size(v)
         WRITE(0,*) "GET_F_AND_DERIVS_FOURIERP4MM ERR: size(vdvdu): ",size(dvdu(:,1))
         CALL EXIT(1)
      CASE(.TRUE.)
         ! do nothing
   END SELECT
   DO i = 1, this%nfunc
      DO j = 1, this%n
         terms(j)   =termfoup4mm( surf=surf,k=this%kList(j,:),parity=this%parityList(j),irrep=this%irrepList(j),r=r )
         terms_dx(j)=termfoup4mm( surf=surf,k=this%kList(j,:),parity=this%parityList(j),irrep=this%irrepList(j),r=r )
         terms_dy(j)=termfoup4mm( surf=surf,k=this%kList(j,:),parity=this%parityList(j),irrep=this%irrepList(j),r=r )
      END DO
      v(i)=dot_product(terms,this%coeff(:,i))
      dvdu(i,1)=dot_product(terms_dx,this%coeff(:,i))
      dvdu(i,2)=dot_product(terms_dy,this%coeff(:,i))
   END DO
   DEALLOCATE(terms)
   DEALLOCATE(terms_dx)
   DEALLOCATE(terms_dy)
   RETURN
END SUBROUTINE GET_F_AND_DERIVS_FOURIERP4MM
END MODULE FOURIER_P4MM_MOD 
