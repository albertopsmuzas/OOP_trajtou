!#########################################################
! MODULE: FOURIER_P4MM_MOD
!> @brief
!! Provides tools to genererate a symmetry adapted fourier
!! interpolation
!##########################################################
MODULE FOURIER_P4MM_MOD
use SYSTEM_MOD, only: system_surface
use FOURIER2D_MOD, only: Fourier2d,TermCalculator2d
use UNITS_MOD, only: pi
use MATHS_MOD, only: INV_MTRX
#ifdef DEBUG
use DEBUG_MOD, only: VERBOSE_WRITE, DEBUG_WRITE
#endif
IMPLICIT NONE
!////////////////////////////////////////////////////////////////
! TYPE: TermCalculator2d_p4mm
!> @brief
!! Type extension of TermCalculator for p4mm symmetry
!----------------------------------------------------------------º
type,extends(TermCalculator2d):: TermCalculator2d_p4mm
   contains
      procedure,public:: getValue  => termFoup4mm
      procedure,public:: getDeriv1 => termFoup4mm_dx
      procedure,public:: getDeriv2 => termFoup4mm_dy
end type TermCalculator2d_p4mm
!/////////////////////////////////////////////////////////
! TYPE: Fourierp4mm
!> @brief
!! Extends Fourier2d interpolation for p4mm symmetry
!---------------------------------------------------------
TYPE,EXTENDS(Fourier2d):: Fourierp4mm
   PRIVATE
   REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: coeff
   CONTAINS
      PROCEDURE,PUBLIC:: INTERPOL => INTERPOL_FOURIERP4MM
      PROCEDURE,PUBLIC:: GET_F_AND_DERIVS => GET_F_AND_DERIVS_FOURIERP4MM
      procedure,public:: initializeTerms => initializeTerms_FOURIERP4MM
END TYPE Fourierp4mm
! variables and types, body
CONTAINS
!###########################################################
!# FUNC: initializeTerms_FOURIERP4MM
!###########################################################
!
!-----------------------------------------------------------
subroutine initializeTerms_FOURIERP4MM(this)
   implicit none
   class(Fourierp4mm),intent(inout):: this
   allocate( TermCalculator2d_p4mm::this%term )
end subroutine initializeTerms_FOURIERP4MM
!###########################################################
!# FUNCTION: termfoup4mm
!###########################################################
!
!-----------------------------------------------------------
function termFoup4mm(this,k,parity,irrep,x) result(answer)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(TermCalculator2d_p4mm),intent(in):: this
   integer(kind=4),dimension(2),intent(in):: k
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   real(kind=8),dimension(2),intent(in):: x
   ! Dummy out variable
   real(kind=8):: answer
   ! Local variables
   real(kind=8):: g
   ! Paramters
   character(len=*),parameter:: routinename='termfoup4mm: '
   ! Run section -------------------------------------------
   g=2.D0*PI/system_surface%norm_s1
   select case( irrep )
   case('A1')
      select case( parity )
      case('+')
         answer=dcos( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) )+&
                dcos( g*dfloat(k(2))*x(1) )*dcos( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('A2')
      select case( parity )
      case('+')
         answer=-dsin( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) )+&
                 dsin( g*dfloat(k(2))*x(1) )*dsin( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('B1')
      select case( parity )
      case('+')
         answer=dcos( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) )-&
                dcos( g*dfloat(k(2))*x(1) )*dcos( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('B2')
      select case( parity )
      case('+')
         answer=dsin( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) )-&
                dsin( g*dfloat(k(2))*x(1) )*dsin( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('EE','E')
      select case( parity )
      case('-')
         answer=dsin( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) )+&
                dcos( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case default
      write(0,*) 'ERR '//routinename//' irrep not implemented: "'//irrep//'"'
      call exit(1)
   end select
   return
end function termFoup4mm
!###########################################################
!# FUNCTION: termfoup4mm_dx 
!###########################################################
!> @brief 
!! ! type brief explanation
!-----------------------------------------------------------
function termfoup4mm_dx(this,k,parity,irrep,x) result(answer)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(TermCalculator2d_p4mm),intent(in):: this
   integer(kind=4),dimension(2),intent(in):: k
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   real(kind=8),dimension(2),intent(in):: x
   ! Dummy out variable
   real(kind=8):: answer
   ! Local variables
   real(kind=8):: g
   ! Paramters
   character(len=*),parameter:: routinename='termfoup4mm_dx: '
   ! Run section -------------------------------------------
   g=2.D0*PI/system_surface%norm_s1
   select case( irrep )
   case('A1')
      select case( parity )
      case('+')
         answer=-g*( dfloat(k(1))*dsin( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) )+&
                     dfloat(k(2))*dsin( g*dfloat(k(2))*x(1) )*dcos( g*dfloat(k(1))*x(2) ) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('A2')
      select case( parity )
      case('+')
         answer=-g*dfloat(k(1))*dcos( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) )+&
                 g*dfloat(k(2))*dcos( g*dfloat(k(2))*x(1) )*dsin( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('B1')
      select case( parity )
      case('+')
         answer=-g*dfloat(k(1))*dsin( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) )+&
                 g*dfloat(k(2))*dsin( g*dfloat(k(2))*x(1) )*dcos( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('B2')
      select case( parity )
      case('+')
         answer=g*dfloat(k(1))*dcos( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) )-&
                g*dfloat(k(2))*dcos( g*dfloat(k(2))*x(1) )*dsin( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('EE','E')
      select case( parity )
      case('-')
         answer=g*dfloat(k(1))*( dcos( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) )-&
                                 dsin( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) ) )
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
function termfoup4mm_dy(this,k,parity,irrep,x) result(answer)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(TermCalculator2d_p4mm),intent(in):: this
   integer(kind=4),dimension(2),intent(in):: k
   character(len=1),intent(in):: parity
   character(len=2),intent(in):: irrep
   real(kind=8),dimension(2),intent(in):: x
   ! Dummy out variable
   real(kind=8):: answer
   ! Local variables
   real(kind=8):: g
   ! Paramters
   character(len=*),parameter:: routinename='termfoup4mm_dy: '
   ! Run section -------------------------------------------
   g=2.D0*PI/system_surface%norm_s1
   select case( irrep )
   case('A1')
      select case( parity )
      case('+')
         answer=-g*( dfloat(k(2))*dcos( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) )+&
                     dfloat(k(1))*dcos( g*dfloat(k(2))*x(1) )*dsin( g*dfloat(k(1))*x(2) ) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('A2')
      select case( parity )
      case('+')
         answer=-g*dfloat(k(2))*dsin( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) )+&
                 g*dfloat(k(1))*dsin( g*dfloat(k(2))*x(1) )*dcos( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('B1')
      select case( parity )
      case('+')
         answer=-g*dfloat(k(2))*dcos( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) )+&
                 g*dfloat(k(1))*dcos( g*dfloat(k(2))*x(1) )*dsin( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('B2')
      select case( parity )
      case('+')
         answer=g*dfloat(k(2))*dsin( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) )-&
                g*dfloat(k(1))*dsin( g*dfloat(k(2))*x(1) )*dcos( g*dfloat(k(1))*x(2) )
      case default
         write(0,*) 'ERR '//routinename//' parity and irrep selection is not compatible'
         call exit(1)
      end select

   case('EE','E')
      select case( parity )
      case('-')
         answer=g*dfloat(k(2))*( -dsin( g*dfloat(k(1))*x(1) )*dsin( g*dfloat(k(2))*x(2) )+&
                                  dcos( g*dfloat(k(1))*x(1) )*dcos( g*dfloat(k(2))*x(2) ) )
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
SUBROUTINE INTERPOL_FOURIERP4MM(this,filename)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   CLASS(Fourierp4mm),INTENT(INOUT) :: this
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
         tmtrx(i,j)=this%term%getValue( k=this%kList(j,:),parity=this%parityList(j),irrep=this%irrepList(j),x=this%xy(i,:) )
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
SUBROUTINE GET_F_AND_DERIVS_FOURIERP4MM(this,r,v,dvdu)
   ! Initial declarations   
   IMPLICIT NONE
   ! I/O variables
   class(Fourierp4mm),INTENT(IN):: this
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
         terms(j)   =this%term%getValue ( k=this%kList(j,:),parity=this%parityList(j),irrep=this%irrepList(j),x=r )
         terms_dx(j)=this%term%getDeriv1( k=this%kList(j,:),parity=this%parityList(j),irrep=this%irrepList(j),x=r )
         terms_dy(j)=this%term%getDeriv2( k=this%kList(j,:),parity=this%parityList(j),irrep=this%irrepList(j),x=r )
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
