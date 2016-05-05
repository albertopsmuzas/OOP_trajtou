!#########################################################
! MODULE: CRP3D_MOD
!> @brief
!! Provides all those routines and types needed to define a CRP3D PES
!
!> @warning
!! - This module uses PES_MOD, SURFACE_MOD and INTERPOL1D_MOD
!
!> @see pes_mod, surface_mod, interpol1d_mod
!##########################################################
module CRP3D_MOD
   use PES_MOD
   use SYSTEM_MOD
   use CUBICSPLINES_MOD
   use LINK_FOURIER2D_MOD
   use LINK_FUNCTION1D_MOD
   use AOTUS_MODULE, only: flu_State, OPEN_CONFIG_FILE, CLOSE_CONFIG, AOT_GET_VAL
   use AOT_TABLE_MODULE, only: AOT_TABLE_OPEN, AOT_TABLE_CLOSE, AOT_TABLE_LENGTH, AOT_TABLE_GET_VAL
   use MATHS_MOD, only: ORDER_VECT, ORDER, SYMMETRIZE
   use UNITS_MOD, only: pi
#ifdef DEBUG
   use DEBUG_MOD, only: VERBOSE_WRITE, DEBUG_WRITE
#endif
! Initial declarations
implicit none
!/////////////////////////////////////////////////////////////////
! TYPE: CRP3D_details
!> @brief
!! Type used in jobs to generate good inputs for a CRP3D PES
!
!> @param vasint - potential when the atom is far from the surface (in the vacuum)
!> @param dfin - arbitrary boundary: 1st derivative of all potentials at zgrid(nzgrid)
!> @param nrumpling - number of different rumplings defined
!> @param rumpling - storage of rumplings
!> @param zeropos - position of ith rumpling in the grid
!> @param nzgrid - dimension of zgrid
!> @param zgrid - grid proposed to generate input files
!> @param first, last - first and last points in the grid

!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jul/2014
!> @version 1.0
!----------------------------------------------------------------
type :: CRP3D_details
private
   real(kind=8):: vasint 
   real(kind=8):: dfin 
   integer(kind=4):: nrumpling 
   real(kind=8),dimension(:),allocatable:: rumpling 
   integer(kind=4),dimension(:),allocatable :: zeropos 
   integer(kind=4) :: nzgrid 
   real(kind=8),dimension(:),allocatable :: zgrid 
   real(kind=8):: first, last 
   contains
      procedure,public :: READ => READ_CRP3D_details
end type CRP3D_details
!///////////////////////////////////////////////////////////////////////////////
! TYPE: Symmetric point
! ---------------------
!
!> @brief 
!! For CRP3D PES. Symmetric point in XY plane. Composed by a scan in Z @f$(Z_{i},F(Z_{i}))@f$
!! variable, perpendicular to the surface.
!
!> @details
!! - Should be chosen inside the irreducible Wegner-Seitz cell of the surface
!> @param filename - File that contains initial information for this point
!> @param alias - Human-friendly alias for the Symmetric point
!> @param n - Number of points in the scan
!> @param x,y - Position in the XY plane
!> @param units_z,units_v - Units in which v(:) and z(:) are stored. Should be compatible
!!                          with units_mod module
!> @param dz1,dz2 - Second derivatives at Z(1) and Z(n)
!> @param interz - One dimensional interpolation for Z variable
!
!> @author A.P. Muzas
!> @version 1.0
!> @date 21/Jan/2014
!------------------------------------------------------------------------------
type :: SymmPoint
   ! private part
   character(len=:),allocatable,private:: filename
   character(len=:),allocatable,private:: alias
   character(len=10),private:: units_z
   character(len=10),private:: units_v
   real(kind=8),dimension(:),allocatable,private:: v
   real(kind=8),private:: dz1
   real(kind=8),private:: dz2
   ! public part
   integer(kind=4),public:: n
   real(kind=8),public:: x
   real(kind=8),public:: y
   real(kind=8),dimension(:),allocatable,public:: z
   type(Csplines),public:: interz 
   contains
      procedure,public:: initializeRaw => initializeRaw_SymmPoint
      procedure,public:: printSymmetrizedRawInput => printSymmetrizedRawInput_SymmPoint
      procedure,public:: PLOT_DATA => PLOT_DATA_SYMMPOINT
      procedure,public:: PLOT => PLOT_INTERPOL_SYMMPOINT
end type Symmpoint
!///////////////////////////////////////////////////////////////////////
! SUBTYPE: Pair potential  
! -----------------------
!> @brief 
!! For CRP3D PES. Symmetric point used to create a pair potential to extract corrugation.
!
!> @details
!! - It is needed one pair potentials for each different kind of atom in the surface
!! - Usually, the pair potentian is built on a top site of the surface.
!
!> @param id - Integer number to identify the pair potential
!> @param vasint - Potential far from the surface
!> @param rumpling - Rumpling associated with this pair potential. As pair potentials
!!                      are built on top sites and are related to different kinds of atoms
!!                      in the surface, the position of the coulomb repulsive singularity
!!                      may not be situated at Z=0. This rumpling is the shift so that
!!                      the singularity is plazed at zero. Should be a negative number o 0.D0.
!
!> @author A.P. Muzas
!> @version 1.0
!> @date 21/Jan/2014
!------------------------------------------------------------------------------
type,extends(SymmPoint) :: Pair_pot
   private
   integer(kind=4):: id 
   real(kind=8):: vasint
	real(kind=8):: rumpling
   contains
      ! Initialization block
      procedure,public:: READ => READ_STANDARD_PAIRPOT
      ! Tools block
      procedure,public:: GET_V_AND_DERIVS => GET_V_AND_DERIVS_PAIRPOT
end type Pair_pot
!/////////////////////////////////////////////////////////////////////
! SUBTYPE: Sitio
! --------------
!> @brief 
!! For CRP3D PES. Symmetric point inside the irreducible Wigner-Seitz cell of the surface.
!! Used just to build up the PES.
!
!> @param dvdx,dvdy,dvdz - Collection of derivatives for each @f$V(Z_{i})@f$
!
!> @author A.P. Muzas
!> @version 1.0
!> @date 21/Jan/2014
!------------------------------------------------------------------------------
type,extends(Symmpoint) :: Sitio
   private
	real(kind=8),dimension(:),allocatable:: dvdx,dvdy,dvdz 
   contains
      procedure,public:: READ => READ_STANDARD_SITIO
end type Sitio
!////////////////////////////////////////////////////////////////////////////////
! SUBTYPE: CRP3D
! ------------
!> @brief 
!! Contains all information that can be extracted or needed during a CRP3D interpolation
!
!> @details
!! - This is a subtype PES. In principle, should be compatible with all available 
!!   dynamics jobs.
!
!> @param max_order - Integer number. Order of the environment used to extract corrugation
!> @param all_pairpots,all_sitios - Compendium of all pair potentials and sitios used
!!                                  in the interpolation procedure
!> @param surf - Surface parameters in which the CRP3D surface is defined
!
!> @warning
!! - Use this type variable and related procedures once you've obtained a reliable set of input
!!   files to read. In order to create such a set, use Newinput variables and member procedures
!!
!> @author A.P. Muzas
!> @version 1.1
!> @date 04/Feb/2014
!
!> @see newinput
!------------------------------------------------------------------------------
type,extends(PES) :: CRP3D
   integer(kind=4):: max_order
   type(Pair_pot),dimension(:),allocatable:: all_pairpots
   type(Sitio),dimension(:),allocatable:: all_sites
   integer(kind=4),dimension(:,:),allocatable:: kList
   character(len=1),dimension(:),allocatable:: parityList
   character(len=2),dimension(:),allocatable:: irrepList
   class(Function1d),allocatable:: dampfunc
   contains
      ! Initialization block
      procedure,public:: READ => READ_CRP3D
      procedure,public:: INITIALIZE => INITIALIZE_CRP3D
      ! Get block 
      procedure,public:: GET_V_AND_DERIVS => GET_V_AND_DERIVS_CRP3D
      procedure,public:: GET_V_AND_DERIVS_CORRECTION => GET_V_AND_DERIVS_CORRECTION_CRP3D
      procedure,public:: GET_V_AND_DERIVS_SMOOTH => GET_V_AND_DERIVS_SMOOTH_CRP3D
      procedure,public:: GET_REPUL_CORRECTIONS => GET_REPUL_CORRECTIONS_CRP3D
      procedure,public:: getpot => getpot_crp3d
      procedure,public:: getrumpling => getrumpling_CRP3D
      ! Enquire block
      procedure,public:: is_allowed => is_allowed_CRP3D
      ! Tools block
      procedure,public:: EXTRACT_VASINT => EXTRACT_VASINT_CRP3D
      procedure,public:: SMOOTH => SMOOTH_CRP3D
      procedure,public:: INTERPOL => INTERPOL_Z_CRP3D
      procedure,public:: RAWINTERPOL => RAWINTERPOL_Z_CRP3D
      procedure,public:: SET_FOURIER_SYMMETRY => SET_FOURIER_SYMMETRY_CRP3D
      ! Plot tools
      procedure,public:: PLOT_XYMAP => PLOT_XYMAP_CRP3D
      procedure,public:: PLOT_XYMAP_SMOOTH => PLOT_XYMAP_CRP3D
      procedure,public:: PLOT_XYMAP_CORRECTION => PLOT_XYMAP_CORRECTION_CRP3D
      procedure,public:: PLOT_DIRECTION1D => PLOT_DIRECTION1D_CRP3D
      procedure,public:: PLOT_DIRECTION1D_SMOOTH => PLOT_DIRECTION1D_SMOOTH_CRP3D
      procedure,public:: PLOT_DIRECTION1D_CORRECTION => PLOT_DIRECTION1D_CORRECTION_CRP3D
      procedure,public:: PLOT_SITIOS => PLOT_SITIOS_CRP3D
      procedure,public:: PLOT_PAIRPOTS => PLOT_PAIRPOTS_CRP3D
      procedure,public:: PLOT_Z => PLOT_Z_CRP3D
      procedure,public:: plot_z_smooth => plot_z_smooth_CRP3D
end type CRP3D
!///////////////////////////////////////////////////////////////////////////
contains
!###########################################################
!# SUBROUTINE: READ_CRP3D_details 
!###########################################################
!> @brief
!! Reads from file interesting details to generate new inputs for a CRP3D PES.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jul/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine READ_CRP3D_details(this,filename)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP3D_details),intent(out):: this
   character(len=*),intent(in) :: filename
   ! WARNING: unit used to read
   integer(kind=4) :: runit=233
   ! Local variables
   integer(kind=4) :: i,k ! counters
   type(Energy):: en
   type(Length):: len
   real(kind=8) :: raux,raux2
   character(len=10) :: units
   character(len=20),parameter :: routinename="READ_CRP3D_details: "
   character(len=4) :: control
   logical:: exists_zero
   ! Run section
   open (runit,file=filename,status="old",action="read")
   read(runit,*) ! dummy line
   read(runit,*) raux,units
   call en%READ(raux,units)
   call en%TO_STD()
   this%vasint=en%getvalue()
   read(runit,*) this%nrumpling
   allocate(this%rumpling(this%nrumpling))
   this%dfin = 0.D0
   do i=1,this%nrumpling
      read(runit,*) raux,units
      call len%READ(raux,units)
      call len%TO_STD()
      this%rumpling(i)=len%getvalue()
   end do
#ifdef DEBUG
   call VERBOSE_WRITE(routinename,"Rumplings: ",this%rumpling(:))
#endif
   read(runit,*) this%nzgrid,control,raux,raux2,units
   allocate(this%zgrid(this%nzgrid))
   call len%READ(raux,units)
   call len%TO_STD()
   this%first=len%getvalue()
   call len%READ(raux2,units)
   call len%TO_STD()
   this%last=len%getvalue()
   select case(control)
      case("MANU") ! manual grid input
         do i = 1, this%nzgrid
            read(runit,*) raux
            call len%READ(raux,units)
            call len%TO_STD()
            this%zgrid(i)=len%getvalue()
         end do
         call ORDER_VECT(this%zgrid)
         allocate(this%zeropos(this%nrumpling))
         forall(i=1:this%nrumpling) this%zeropos(i)=0
         do i=1,this%nzgrid
            do k=1,this%nrumpling
               select case(this%zgrid(i)==this%rumpling(k))
                  case(.true.)
                     this%zeropos(k)=i
                  case(.false.)
                     ! do nothing
               end select
            end do
         end do
         exists_zero=.true.
         do k=1,this%nrumpling
            select case(this%zeropos(k))
               case(0)
                  exists_zero=.false.
               case default
                  ! do nothing
            end select
         end do
         select case(exists_zero)
            case(.false.)
               write(0,*) "READ_CRP3D_details ERR: Manual grid does not have a rumpling point in the grid"
               call EXIT(1)
            case(.true.)
               ! do nothing
         end select
      !
      case default
         write(0,*) "READ_CRP3D_details ERR: wrong grid control keyword"
         write(0,*) "You geave: ", control
         write(0,*) "Implemented ones: MANU"
         write(0,*) "Warning: case-sensitive"
         call EXIT(1)
   end select

   close(runit)
   return
end subroutine READ_CRP3D_details
!###########################################################
!# SUBROUTINE: SET_FOURIER_SYMMETRY_CRP3D
!###########################################################
!> @brief
!! Given a Fourier2d allocatable variable, allocates it with the correct surface
!! symmetry
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jul/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine SET_FOURIER_SYMMETRY_CRP3D(this,interpolxy)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP3D),intent(in):: this
   class(Fourier2d),allocatable,intent(inout):: interpolxy
   ! Run section
   select case(system_surface%getsymmlabel())
      case("p4mm")
         allocate(Fourierp4mm::interpolxy)
      case default
         write(0,*) "SET_FOURIER_SYMMETRY_CRP3D ERR: Incorrect surface symmlabel"
         write(0,*) "Used: ",system_surface%getsymmlabel()
         write(0,*) "Implemented ones: p4mm, p6mm"
         call EXIT(1)
   end select
   return
end subroutine SET_FOURIER_SYMMETRY_CRP3D
!###########################################################
!# SUBROUTINE: INITIALIZE_CRP3D 
!###########################################################
!> @brief
!! Specific implementation of initialize PES from input file
!
!> @param[in] filename - char(len=*),optional: input file. If system_inputfile
!!                       is allocated, it is ignored.
!> @param[in] tablename - char(len=*),optional: name of Lua table where the PES is
!!                        stored. By default, it has the value 'pes'.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014 
!> @version 1.0
!-----------------------------------------------------------
subroutine INITIALIZE_CRP3D(this,filename,tablename)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP3D),intent(out):: this
   character(len=*),optional,intent(in):: filename,tablename
   ! Local variables
   character(len=:),allocatable:: auxstring
   ! Run section
   select case(allocated(system_inputfile) .or. .not.present(filename))
      case(.true.)
         auxstring=trim(system_inputfile)
      case(.false.)
         auxstring=trim(filename)
   end select
   select case(present(tablename))
      case(.true.) ! present tablename
         call this%READ(filename=trim(auxstring),tablename=trim(tablename))
      case(.false.) ! not present tablename
         call this%READ(filename=trim(auxstring),tablename='pes')
   end select
   call this%INTERPOL()
   return
end subroutine INITIALIZE_CRP3D
!###########################################################
!# SUBROUTINE: GET_REPUL_CORRECTIONS_CRP3D 
!###########################################################
!> @brief
!! Collects interactions 
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine GET_REPUL_CORRECTIONS_CRP3D(this,P,v,dvdz,dvdx,dvdy)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP3D),intent(in) :: this
   real(kind=8),dimension(3),intent(in) :: P
   real(kind=8),dimension(:),intent(out) :: v
   real(kind=8),dimension(:),intent(out):: dvdx,dvdy,dvdz ! corrections to the derivatives
   ! Local variables
   integer(kind=4) :: npairpots
   integer(kind=4) :: l,k ! counters
   real(kind=8) :: aux1,aux2,aux3,aux4
   ! Run section
   npairpots=size(this%all_pairpots)
   forall(l=1:npairpots)
      v(l)=0.D0
      dvdz(l)=0.D0
      dvdx(l)=0.D0
      dvdy(l)=0.D0
   end forall
   do l = 1, npairpots
      do k = 0, this%max_order
         call INTERACTION_AENV(k,P,this%all_pairpots(l),this%dampfunc,aux1,aux2,aux3,aux4)
         v(l)=v(l)+aux1
         dvdz(l)=dvdz(l)+aux2
         dvdx(l)=dvdx(l)+aux3
         dvdy(l)=dvdy(l)+aux4
      end do
   end do
   return
end subroutine GET_REPUL_CORRECTIONS_CRP3D
!###########################################################
!# SUBROUTINE: GET_V_AND_DERIVS_PAIRPOT 
!###########################################################
!> @brief
!! Gets the pairpot potential shifted to 0 (maximum interaction is at 0)
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine GET_V_AND_DERIVS_PAIRPOT(this,x,v,dvdu)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(Pair_pot),intent(in):: this
   real(kind=8),intent(in) :: x
   real(kind=8),intent(out) :: v,dvdu
   ! Run section
   select case(X>this%z(this%n)-this%rumpling)
      case(.true.)
         v=0.D0
         dvdu=0.D0
      case(.false.)
         call this%interz%GET_V_AND_DERIVS(x,v,dvdu,this%rumpling)
   end select
   return
end subroutine GET_V_AND_DERIVS_PAIRPOT
!###########################################################
!# FUNCTION: getrumpling_CRP3D 
!###########################################################
!> @brief 
!! Common get function. Gets rumpling value of parit potential "toptype"
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date MAy/2014
!> @version 1.0
!-----------------------------------------------------------
real(kind=8) function getrumpling_CRP3D(this,toptype) 
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP3D),intent(in):: this
   integer(kind=4),intent(in) :: toptype
   ! Run section
   getrumpling_CRP3D=this%all_pairpots(toptype)%rumpling
   return
end function getrumpling_CRP3D
!######################################################################
! SUBROUTINE: initializeRaw_SymmPoint #################################
!######################################################################
!> @brief 
!! Initializes all data contained in a generic Symmpoint type variable 
!! from a "raw" input file.
!
!> @details
!! - Data is stored in a.u. automatically despite the units in input file.
!
!> @param[in,out] symmraw - A Symmpoint type variable. All its members will be loaded
!!                          from the input file
!> @param[in] filename - File that contains the "raw" input file
!
!> @author A.P. Muzas
!> @version 1.0
!> @date 22/Jan/2014
!
!> @warning
!! - This is the @b format needed for the "raw" input file:
!!    -# line 1: @b real(kind=8),@b real(kind=8),@b character(len=10); X location, Y location, length units
!!    -# line 2: @b integer(kind=4); Number of points (N) contained in the file
!!    -# line 3: @b character(len=10),@b character(len=10); unit codes for Z(length) and potential(energy)
!!    -# lines 4~4+N: @b real(kind=8),@b real(kind=8); cuples @f$(Z_{i},V(Z_{i}))@f$. They should
!!                    be consistent with the previous choice of units. There is no need of a specific ordering
!!                    of the data
!!
!! - This way of initializing a generic SymmPoint sets dz2=0.d0
!
!> @see units_mod, debug_mod, maths_mod
!---------------------------------------------------------------------- 
subroutine initializeRaw_SymmPoint(this,fileName)
   implicit none
   ! I/O Variables ----------------
   class(SymmPoint),intent(inout):: this
   character(len=*),intent(in):: fileName
   ! Local variables
   integer:: i ! Counter
   real(kind=8):: aux_r1, aux_r2
   real(kind=8),dimension(:),allocatable :: aux_v1,aux_v2
   character(len=20),parameter:: routinename = "READ_SYMMPOINT_RAW: "
   type(Length):: x,y,z
   type(Energy):: v
   character(len=10):: units1,units2
   ! FIRE IN THE HOLE!
   this%filename=filename
   open(10,file=this%filename,status="OLD")
      read(10,*) aux_r1,aux_r2,units1
      call x%READ(aux_r1,units1)
      call y%READ(aux_r2,units1)
      this%x = x%getvalue()
      this%y = y%getvalue()
      ! ------
      read(10,*) this%n
      read(10,*) units1, units2
      allocate(this%v(this%n))
      allocate(this%z(this%n))
      allocate(aux_v1(this%n))
      allocate(aux_v2(this%n))
      do i=1,this%n
         read(10,*) aux_r1,aux_r2
         call z%READ(aux_r1,units1)
         call v%READ(aux_r2,units2)
         call z%TO_STD() ! go to standard units (a.u.)
         call v%TO_STD() ! go to standard units (a.u.)
         aux_v1(i) = z%getvalue()
         aux_v2(i) = v%getvalue()
      end do
   close(10)
   call ORDER(aux_v1,aux_v2) ! Should order Z(i) and V(i)
   this%z = aux_v1
   this%v = aux_v2
   this%units_z = z%getunits()
   this%units_v = v%getunits()
   ! Raw interpolation
   call this%interZ%READ(x=this%z, f=this%v)
   call this%interZ%INTERPOL(dz1=0.d0,id1=0,dz2=0.d0,id2=1)
   this%dz1=this%interZ%getderiv2(x=this%z(1))
   this%dz2=0.d0
#ifdef DEBUG
   ! Write status message
   call VERBOSE_WRITE(routinename, symmraw%filename, "Done")
#endif
   return
end subroutine initializeRaw_SymmPoint
!###########################################################
!# SUBROUTINE: READ_STANDARD_PAIRPOT #######################
!###########################################################
!> @brief
!! Initializes all data for a Pair_pot subtype variable from
!! an standard input file, i.e. one written in a.u.
!
!> @param[in,out] pairpot - Pair_pot subtype variable to be initialized
!> @param[in] filename - String with the name of the standard pairpot
!!                        file to be loaded
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 24/Jan/2014
!> @version 1.0
!
!> @warning
!! - There is a routine to create standard input files from "raw" ones.
!! - This is the @b format needed for a standard pair potential input:
!!    -# line 1: Dummy line, just for some comments
!!    -# line 2: Dummy line, just for some comments
!!    -# line 3: @b character(len=30); human-frienfly alias for this pair potential
!!    -# line 4: @b real(kind=8); Energy at long dinstance from the surface
!!    -# line 5: @b real(kind=8); @f$\partial^{2}{V(Z_{1})}\over{\partial{Z}^{2}}@f$
!!    -# line 6: @b real(kind=8); @f$\partial^{2}{V(Z_{N})}\over{\partial{Z}^{2}}@f$
!!    -# line 7: @b integer, @b real(kind=8); integer ID number for the pair potential, Z value
!!               in which @f$V(Z)=\pm\infty@f$
!!    -# line 8: @b integer; Number of points (N)
!!    -# line 9~9+N: @b real(kind=8), @b real(kind=8); couples @f$(Z_{i},V(Z_{i})), Z_{i}<Z_{i+1}@f$ in a.u.
!
!> @see 
!! debug_mod
!-----------------------------------------------------------
subroutine READ_STANDARD_PAIRPOT(pairpot,filename)
      ! Initial declarations
   implicit none
   ! I/O variables ------------------------------
   class(Pair_pot),intent(inout) :: pairpot
   character(len=*),intent(in) :: filename
   ! Local variables ----------------------------
   integer :: i
   character(len=14), parameter :: routinename = "READ_PAIRPOT: "
   character(len=1024):: auxString
   ! Run section ---------------------------------
   pairpot%filename=filename
   open(10,file=pairpot%filename,status='OLD')
   read(10,*) ! dummy line
   read(10,*) ! dummy line
   read(10,*) auxString
   pairpot%alias=trim(auxString)
   read(10,*) pairpot%vasint
   read(10,*) pairpot%dz1
   read(10,*) pairpot%dz2
   read(10,*) pairpot%id,pairpot%rumpling
   read(10,*) pairpot%n
   allocate(pairpot%z(1:pairpot%n))
   allocate(pairpot%v(1:pairpot%n))
   do i=1, pairpot%n
      read(10,*) pairpot%z(i), pairpot%v(i)
   end do
   close(10)
   call pairpot%interz%READ(pairpot%z,pairpot%v)
#ifdef DEBUG
      call VERBOSE_WRITE(routinename,pairpot%filename,pairpot%alias)
#endif
   return
end subroutine READ_STANDARD_PAIRPOT
!###########################################################
!# SUBROUTINE: READ_STANDARD_SITIO #########################
!###########################################################
!> @brief
!! Initializes all data for a Sitio subtype variable from
!! an standard input file, i.e. one written in a.u.
!
!> @param[in,out] site - Sitio subtype variable to be initialized
!> @param[in] filename - String with the name of the standard site
!!                        file to be loaded
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 24/Jan/2014
!> @version 1.0
!
!> @warning
!! - There is a routine to create standard input files from "raw" ones.
!! - This is the @b format needed for a standard Sitio input:
!!    -# line 1: Dummy line, just for some comments
!!    -# line 2: Dummy line, just for some comments
!!    -# line 3; @b character(len=30); human_friendly alias for this site
!!    -# line 4: @b real(kind=8),@b real(kind=8); X, Y location in a.u.
!!    -# line 5: @b integer; Number of points (N)
!!    -# line 6: @b real(kind=8); @f$\partial{V(Z_{1})}\over{\partial{Z}}@f$
!!    -# line 7: @b real(kind=8); @f$\partial{V(Z_{N})}\over{\partial{Z}}@f$
!!    -# line 8~8+N: @b real(kind=8),@b real(kind=8); couples @f$(Z_{i},V(Z_{i})), Z_{i}<Z_{i+1}@f$ in a.u.
!
!-----------------------------------------------------------
subroutine READ_STANDARD_SITIO(site,filename)
   implicit none
   ! I/O variables
   class(Sitio),intent(inout):: site
   character(len=*),intent(in):: filename
   ! Local variables
   integer:: i ! counter
   character(len=*),parameter:: routinename = "READ_SITIO: "
   character(len=1024):: auxString
   !
   site%filename=filename
   open (10,file=site%filename,status="old")
   read(10,*) ! dummy line
   read(10,*) ! dummy line
   read(10,*) auxString
   site%alias=trim(auxString)
   read(10,*) site%x, site%y
   read(10,*) site%n
   read(10,*) site%dz1
   read(10,*) site%dz2
   allocate(site%z(1:site%n))
   allocate(site%v(1:site%n))
   do i=1, site%n
      read(10,*) site%z(i), site%v(i)
   end do
   close(10)
   call site%interz%READ(site%z,site%v)
#ifdef DEBUG
   call VERBOSE_WRITE(routinename,site%filename,site%alias)
#endif
   return
end subroutine READ_STANDARD_SITIO
!####################################################################
! SUBROUTINE:  GET_SYMMETRIZED_RAW_INPUT ############################
!####################################################################
!> @brief
!! Tool that adds new points to a "raw" input (can be other kind of input)
!! so that the function is symmetric respect to @b zero. These changes
!! are made only for those couples that satisfy @f$V(z_{i})> f_{0}@f$. Changes
!! are stored in a file called @b filename, which is a standard "raw" input.
!
!> @param[in,out] symmraw - Symmpoint variable to be symmetrized
!> @param[in,out] zero - Length quantity, contains @f$x_{0}@f$. After execution, in a.u.
!> @param[in,out] vtop - Energy quantity, contains @f$f_{0}@f$. After execution, in a.u.
!> @param[in] filename - Name of file to generate
!
!> @warning
!! - This subroutine should be used carefully. The only aim is to add new points
!!   to a specific symmpoint instead of performing new DFT calculations. 
!!   If @b vtop is chosen wisely (high enough) we do not
!!   commit any error projecting points in the sorroundings of @b zero. Repulsive 
!!   potentials in the vecinity of a singularity are symmetric.
!! - @b Symmraw, should have been read before applying this routine.
!
!> @author A.s. Muzas - alberto.muzas@uam.es
!> @date 03/Feb/2014
!> @version 1.0
!
!> @see symmetrize, maths_mod, units_mod
!--------------------------------------------------------------------
subroutine printSymmetrizedRawInput_SymmPoint(this,zero,vtop,filename)
	implicit none
	! I/O variable ---------------------------
	class(SymmPoint),intent(inout):: this
   type(Length),intent(inout):: zero
	type(Energy),intent(inout):: vtop
	character(len=*),intent(in):: fileName
	! Local variable -------------------------
	real(kind=8),dimension(:),allocatable:: x,v
	integer(kind=4):: n
	integer(kind=4):: i ! Counter
   character(len=*),parameter:: routineName="printSymmetrizedRawInput_Symmpoint: "
	! GABBA GABBA HEY! ------------------------
	call zero%TO_STD()
	call vtop%TO_STD()
	call symmetrize(this%n,this%z,this%v,zero%getvalue(),vtop%getvalue(),n,x,v)
	deallocate(this%z)
	deallocate(this%v)
	this%n = n
	allocate(this%z(this%n))
	allocate(this%v(this%n))
	do i=1, this%n
		this%z(i)=x(i)
		this%v(i)=v(i)
	end do
	! Store data in filename ----------------
	open(11,file=fileName,status="replace")
	write(11,*) this%x, this%y, ' au    <----(X,Y) location in a.u.'
	write(11,*) this%n , " SYMM job with V(z) > (a.u.) ", vtop%getvalue()
	write(11,*) "au         au       <---- everything in a.u."
	do i=1, this%n
		write(11,*) this%z(i), this%v(i)
	end do
	close(11)
#ifdef DEBUG
   call verbose_write(routinename,"Symmetrized input generated: ",filename)
#endif
	return
end subroutine printSymmetrizedRawInput_SymmPoint
!###########################################################
!# SUBROUTINE: READ_CRP3D
!###########################################################
!> @brief
!! Read CRP3D input file from Lua config file.
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine READ_CRP3D(this,filename,tablename)
   implicit none
   ! I/O variables
   class(CRP3D),intent(out) :: this 
   character(len=*),intent(in) :: filename
   character(len=*),intent(in):: tablename
   ! Local variables
   integer(kind=4) :: n_pairpots,n_sites
   character(len=1024),dimension(:),allocatable:: files_pairpots,files_sites
   real(kind=8),dimension(:),allocatable:: param
   integer(kind=4):: ierr
   integer(kind=4) :: i ! counter
   ! Lua-related variables
   type(flu_State):: conf ! Lua state
   integer(kind=4):: pes_table,pairpot_table,sitio_table,dampfunc_table,param_table,fourier_table ! tables
   integer(kind=4),dimension(:),allocatable:: subtables
   ! Auxiliar (dummy) variables
   integer(kind=4):: auxint
   character(len=1024):: auxstring
   ! Parameters
   character(len=*),parameter :: routinename="READ_CRP3D: "
   ! HEY HO!, LET'S GO!! ------------------
   ! Open Lua file
   call OPEN_CONFIG_FILE(L=conf,filename=filename,ErrCode=ierr)
   select case(ierr)
      case(0)
         ! do nothing
      case default
         write(0,*) "READ_CRP3D ERR: error reading Lua config file: ",filename
         call EXIT(1)
   end select
   ! Open PES table
   call AOT_TABLE_OPEN(L=conf,thandle=pes_table,key=tablename)
   ! Set pestype (kind)
   call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='kind',val=auxstring)
   call this%SET_PESTYPE(trim(auxstring))
   select case(trim(auxstring))
      case('CRP3D')
         ! do nothing
      case default
         write(0,*) "READ_CRP3D ERR: wrong type of PES. Expected: CRP3D. Encountered: "//trim(auxstring)
         call EXIT(1)
   end select
#ifdef DEBUG
   call VERBOSE_WRITE(routinename,'Type of PES: '//trim(auxstring))
#endif
   ! Set alias (name)
   call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='name',val=auxstring)
   call this%SET_ALIAS(trim(auxstring))
#ifdef DEBUG
   call VERBOSE_WRITE(routinename,'PES Name: '//trim(auxstring))
#endif
   ! Set dimensions
   call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='dimensions',val=auxint)
   call this%SET_DIMENSIONS(auxint)
#ifdef DEBUG
   call VERBOSE_WRITE(routinename,'PES dimensions: ',auxint)
#endif
   ! Set max environment
   call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pes_table,key='maxEnvironment',val=auxint)
   this%max_order=auxint
   ! Set pair potentials
   call AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=pairpot_table,key='pairPotentials')
   n_pairpots=aot_table_length(L=conf,thandle=pairpot_table)
   allocate(files_pairpots(n_pairpots))
   allocate(this%all_pairpots(n_pairpots))
   do i = 1, n_pairpots
      call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=pairpot_table,pos=i,val=files_pairpots(i))
      call this%all_pairpots(i)%READ(trim(files_pairpots(i)))
   end do
   call AOT_TABLE_CLOSE(L=conf,thandle=pairpot_table)
   ! Set damping function
   call AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=dampfunc_table,key='dampFunction')
   call AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=dampfunc_table,key='kind',val=auxstring)
   select case(trim(auxstring))
      case("Logistic")
         allocate(Logistic_func::this%dampfunc)
         allocate(param(2))
         ! open param table
         call AOT_TABLE_OPEN(L=conf,parent=dampfunc_table,thandle=param_table,key='param')
         auxint=aot_table_length(L=conf,thandle=param_table)
         select case(auxint/=2) 
            case(.true.)
               write(0,*) "READ_CRP3D ERR: wrong number of parameters in pes.dampFunc.param table"
               call EXIT(1)
            case(.false.)
               ! do nothing
         end select
         do i = 1, 2
            call AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=param_table,pos=i,val=param(i))
         end do
         call this%dampfunc%READ(param)
         call AOT_TABLE_CLOSE(L=conf,thandle=param_table)
      case("fullCRP")
         allocate(One_func::this%dampfunc)
      case("fullRaw")
         allocate(Zero_func::this%dampfunc)
      case default
         write(0,*) "READ_CRP3D ERR: dampfunction keyword is not implemented"
         write(0,*) "Implemented ones: Logistic, fullCRP, fullRaw"
         write(0,*) "Case sensitive"
         call EXIT(1)
   end select
   call AOT_TABLE_CLOSE(L=conf,thandle=dampfunc_table)
#ifdef DEBUG
   call VERBOSE_WRITE(routinename,'Damping function used: '//trim(auxstring))
   select case(allocated(param))
      case(.true.)
         call VERBOSE_WRITE(routinename,'Damping function parameters: ',param(:))
      case(.false.)
         ! do nothing
   end select
#endif
   ! Set sitios
   call AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=sitio_table,key='sitios')
   n_sites=aot_table_length(L=conf,thandle=sitio_table)
   allocate(files_sites(n_sites))
   allocate(this%all_sites(n_sites))
   do i = 1, n_sites
      call AOT_TABLE_GET_VAL(L=conf,ErrCode=ierr,thandle=sitio_table,pos=i,val=files_sites(i))
      call this%all_sites(i)%READ(trim(files_sites(i))) 
   end do
   call AOT_TABLE_CLOSE(L=conf,thandle=sitio_table)
   ! Read fourier Kpoints
   call AOT_TABLE_OPEN(L=conf,parent=pes_table,thandle=fourier_table,key='fourierKpoints')
   auxint=aot_table_length(L=conf,thandle=fourier_table)
   allocate( this%irrepList(auxInt) )
   allocate( this%parityList(auxInt) )
   this%irrepList(:)='A1'
   this%parityList(:)='+'
   select case(auxint/=n_sites)
      case(.true.)
         write(0,*) "READ_CRP3D ERR: dimension mismatch between fourierKpoints and number of sitios"
         write(0,*) 'Number of Fourier Kpoints: ',auxint
         write(0,*) 'Number of sitios: ',n_sites
         call EXIT(1)
      case(.false.)
         ! do nothing
   end select
   allocate(this%klist(n_sites,2))
   allocate(subtables(n_sites))
   do i = 1, n_sites
      call AOT_TABLE_OPEN(L=conf,parent=fourier_table,thandle=subtables(i),pos=i)
      call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=subtables(i),pos=1,val=auxint)
      this%klist(i,1)=auxint
      call AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=subtables(i),pos=2,val=auxint)
      this%klist(i,2)=auxint
      call AOT_TABLE_CLOSE(L=conf,thandle=subtables(i))
   end do
   call AOT_TABLE_CLOSE(L=conf,thandle=fourier_table)
   call AOT_TABLE_CLOSE(L=conf,thandle=pes_table)
   call CLOSE_CONFIG(conf)
   ! VERBOSE PRINT 
#ifdef DEBUG
   call VERBOSE_WRITE(routinename,"Maximum environmental order: ",this%max_order)
   call VERBOSE_WRITE(routinename,"Number of pair potentials: ",n_pairpots)
   call VERBOSE_WRITE(routinename,"Pair potentials input files:")
   do i = 1, n_pairpots
      call VERBOSE_WRITE(routinename,trim(files_pairpots(i)))
   end do
   call VERBOSE_WRITE(routinename,"Number of sitios: ",n_sites)
   call VERBOSE_WRITE(routinename,"Sitios input files:")
   do i = 1, n_sites
      call VERBOSE_WRITE(routinename,trim(files_sites(i)))
   end do
   call VERBOSE_WRITE(routinename,"List of Kpoints for Fourier interpolation: ")
   do i = 1, n_sites
      call VERBOSE_WRITE(routinename,this%klist(i,:))
   end do
#endif
   return
end subroutine READ_CRP3D
!#######################################################################
!# SUBROUTINE: EXTRACT_VASINT_CRP3D ######################################
!#######################################################################
!> @brief
!! Extracts the potential at long distances from every value stored in pairpots
!! and sites that belongs to @b this. Basically, this procedure sets
!! our zero potential value at the vacuum.
!
!> @param[in,out] this - A CRP3D subtype variable
!
!> @warning 
!! - The value of the potential at long distance should be the same in all
!!   pair potentials, read from files previously
!! - At least one pair potential should've defined in the CRP3D variable
!! - Obviously, @ this should contain reliable data
!! - Only modifies things inside @b interz
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 04/Feb/2014
!> @version 1.0
!------------------------------------------------------------------------
subroutine EXTRACT_VASINT_CRP3D(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(CRP3D),intent(inout) :: this
   ! Local variables
   integer(kind=4) :: npairpots, nsites
   integer(kind=4) :: i,j ! counters
   real(kind=8) :: control_vasint
   character(len=22) :: routinename="EXTRACT_VASINT_CRP3D: "
   ! Run section ------------------------
   npairpots=size(this%all_pairpots)
   control_vasint=this%all_pairpots(1)%vasint
   do i = 1, npairpots
      if (this%all_pairpots(1)%vasint/=control_vasint) then
         write(0,*) "EXTRACT_VASINT_CRP3D ERR: Incoherences in vasint values found"
         write(0,*) "EXTRACT_VASINT_CRP3D ERR: vasint's value at pairpot",1,control_vasint
         write(0,*) "EXTRACT_VASINT_CRP3D ERR: vasint's value at pairpot",i,control_vasint
         call EXIT(1)
      end if
      do j = 1, this%all_pairpots(i)%n
         this%all_pairpots(i)%interz%f(j)=this%all_pairpots(i)%interz%f(j)-this%all_pairpots(i)%vasint
      end do
#ifdef DEBUG
      call DEBUG_WRITE(routinename,"Vasint extracted from pair potential ",i)
#endif
   end do
   nsites=size(this%all_sites)
   do i = 1, nsites
      do j = 1, this%all_sites(i)%n
         this%all_sites(i)%interz%f(j)=this%all_sites(i)%interz%f(j)-this%all_pairpots(1)%vasint
      end do
#ifdef DEBUG
      call DEBUG_WRITE(routinename,"Vasint extracted from pair site ",i)
#endif
   end do
   return
end subroutine EXTRACT_VASINT_CRP3D
!############################################################
! SUBROUTINE: SMOOTH_CRP3D ####################################
!############################################################
!> @brief
!! Substracts to all sitio potentials the pair-potential interactions 
!! with all atoms of the surface inside an environment of order @b
!! max order
! 
!> @param[in,out] this - CRP3D PES variable to be smoothed (only it's sites)
!
!> @warning
!! - Vasint should've  been extracted before
!! - Only afects values stored in @b interz variable
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 05/Feb/2014
!> @version 1.0
!
!> @see extract_vasint_CRP3D, initialize_CRP3D_pes
!-----------------------------------------------------------
subroutine SMOOTH_CRP3D(this)
   ! Initial declaraitons
   implicit none
   class(CRP3D),intent(inout):: this
   ! Local variables
   real(kind=8),dimension(3):: A
   integer(kind=4):: i,j ! counters
   integer(kind=4):: npairpots,nsites
   real(kind=8),dimension(:),allocatable:: v,dvdzr,dummy
   character(len=*),parameter:: routinename="SMOOTH_CRP3D: "
   ! Run section ----------
   nsites = size(this%all_sites)
   npairpots = size(this%all_pairpots)
   allocate(v(npairpots))
   allocate(dvdzr(npairpots))
   allocate(dummy(npairpots))
   do i = 1, nsites ! loop over sites
      A(1) = this%all_sites(i)%x
      A(2) = this%all_sites(i)%y
      do j = 1, this%all_sites(i)%n ! loop over pairs v,z
         A(3)=this%all_sites(i)%z(j)
         call this%GET_REPUL_CORRECTIONS(A,v,dvdzr,dummy,dummy)
         this%all_sites(i)%interz%f(j)=this%all_sites(i)%interz%f(j)-sum(v)
         if (j.eq.1) then
            this%all_sites(i)%dz1=this%all_sites(i)%dz1-sum(dvdzr) ! correct first derivative
         else if (j.eq.this%all_sites(i)%n) then
            this%all_sites(i)%dz2=this%all_sites(i)%dz2-sum(dvdzr) ! correct first derivative
         end if
      end do
#ifdef DEBUG
      call VERBOSE_WRITE(routinename,"Site smoothed: ",i)
#endif
   end do
   return
end subroutine SMOOTH_CRP3D
!############################################################
! SUBROUTINE: INTERPOL_Z_CRP3D ################################
!############################################################
!> @brief
!! Interpolates all sitios and pairpot potentials that belong
!! to an specific CRP3D PES. Sitio potentials are smoothed
!
!> @param[in,out] this - CRP3D PES to be interpolated 
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 05/Feb/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine INTERPOL_Z_CRP3D(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(CRP3D),intent(inout) :: this
   ! Local variables
   integer(kind=4) :: nsites,npairpots
   real(kind=8) :: dz1,dz2
   integer(kind=4) :: i ! counters
   ! Run secton------------------------
   nsites=size(this%all_sites)
   npairpots=size(this%all_pairpots)
   call this%EXTRACT_VASINT()
   do i = 1, npairpots ! loop pairpots
      dz1=this%all_pairpots(i)%dz1
      dz2=this%all_pairpots(i)%dz2
      call this%all_pairpots(i)%interz%INTERPOL(dz1,1,dz2,1)
   end do
   call this%SMOOTH()
   do i = 1, nsites
      dz1=this%all_sites(i)%dz1
      dz2=this%all_sites(i)%dz2
      call this%all_sites(i)%interz%INTERPOL(dz1,1,dz2,1) 
   end do
   return
end subroutine INTERPOL_Z_CRP3D
!############################################################
! SUBROUTINE: RAWINTERPOL_Z_CRP3D ################################
!############################################################
!> @brief
!! Interpolates all sitios and pairpot potentials that belong
!! to an specific CRP3D PES. Sitio potentials aren't smoothed
!
!> @param[in,out] this - CRP3D PES to be interpolated 
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine RAWINTERPOL_Z_CRP3D(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(CRP3D) :: this
   ! Local variables
   integer(kind=4) :: nsites,npairpots
   real(kind=8) :: dz1,dz2
   integer(kind=4) :: i ! counters
   ! Run secton------------------------
   nsites=size(this%all_sites)
   npairpots=size(this%all_pairpots)
   call this%EXTRACT_VASINT()
   do i = 1, npairpots ! loop pairpots
      dz1=this%all_pairpots(i)%dz1
      dz2=this%all_pairpots(i)%dz2
      call this%all_pairpots(i)%interz%INTERPOL(dz1,1,dz2,1)
   end do
   do i = 1, nsites
      dz1=this%all_sites(i)%dz1
      dz2=this%all_sites(i)%dz2
      call this%all_sites(i)%interz%INTERPOL(dz1,1,dz2,1) 
   end do
   return
end subroutine RAWINTERPOL_Z_CRP3D
!##################################################################
!# SUBROUTINE: INTERACTION_AP #####################################
!##################################################################
!> @brief
!! Calculates the interaction between two given 3D Points.
!! The potential is extracted from @b pairpot, which should've been interpolated before
!
!> @param[in] A, P - Pair of 3D points that are interacting through a @b pairpot @b potential
!> @param[in] pairpot - Pair potential. Source of @f$V(r)@f$, where @b r is the distance
!!                      between @b A and @b P. @f$V(r)@f$ is just a shifted version of
!!                      @f$V(z)@f$.
!> @param[out] interac - Actual interaction between @b A and @b P
!> @param[out] dvdz_corr - Corrections to first derivatives: 
!!                         @f$\frac {\partial V(r)}{\partial z}@f$
!> @param[out] dvdx_corr - Corrections to first derivatives: 
!!                         @f$\frac {\partial V(r)}{\partial x}@f$
!> @param[out] dvdy_corr - Corrections to first derivatives: 
!!                         @f$\frac{\partial V(r)}{\partial y}@f$
!
!> @warning
!! - Sitio and Pairpot should have been shifted so that the interaction
!!   at large distance from the surface is 0 (extracted vasint)
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 04/Feb/2014
!> @version 1.0
!
!> @see extra documentation
!------------------------------------------------------------------
subroutine INTERACTION_AP(A,P,pairpot,dampfunc,interac,dvdz_corr,dvdx_corr,dvdy_corr)
   implicit none
   ! I/O VAriables --------------------------------------------
   real(kind=8),dimension(3),intent(in):: A, P
   type(Pair_pot),intent(in):: pairpot
   class(Function1d),intent(in):: dampfunc
   real(kind=8),intent(out):: interac,dvdz_corr,dvdx_corr,dvdy_corr
   ! Local variables ------------------------------------------
   real(kind=8):: r ! distance
   real(kind=8):: v,pre
   real(kind=8):: aux ! dv/dr
   character(len=*),parameter:: routinename = "INTERACTION_AP: "
   ! GABBA, GABBA HEY! ----------------------------------------
   ! Find the distance between A and P, in a.u.
   r=dsqrt((A(1)-P(1))**2.D0+(A(2)-P(2))**2.D0+(A(3)-P(3))**2.D0)
   call pairpot%GET_V_AND_DERIVS(r,v,aux)
   interac=v*dampfunc%getvalue(r)
   pre=aux*dampfunc%getvalue(r)+v*dampfunc%getderiv(r)
   dvdz_corr=pre*(A(3)-P(3))/r
   dvdx_corr=pre*(A(1)-P(1))/r
   dvdy_corr=pre*(A(2)-P(2))/r
   return
end subroutine INTERACTION_AP
!##################################################################
!# SUBROUTINE: INTERACTION_AENV ###################################
!##################################################################
!> @brief
!! Calculates the complete interaction between a given point @b A with an
!! environment of order @b n. Depending on the surface, it can use a hexagonal
!! environment or an octagonal one. If a new surface symmetry is implemented, this routine should be
!! edited as well.
!
!> @param[in] n - Environment order
!> @param[in] A - Point in space, cartesian coordinates
!> @param[in] pairpot - Pair potential which defines interaction potential
!> @param[out] interac - Total interaction with the environment defined
!> @param[out] dvdx_term - 
!> @param[out] dvdy_term - 
!> @param[out] dvdz_term - 
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Feb/2014;Jun/2014
!> @version 2.0
!------------------------------------------------------------------
subroutine INTERACTION_AENV(n,A,pairpot,dampfunc,interac,dvdz_term,dvdx_term,dvdy_term)
   ! Initial declarations   
   implicit none
   ! I/O variables
   integer,intent(in):: n
   real(kind=8),dimension(3),intent(in):: A
   type(Pair_pot),intent(in):: pairpot
   class(Function1d),intent(in):: dampfunc
   real(kind=8),intent(out):: interac, dvdz_term, dvdx_term, dvdy_term
   ! Run section
   select case(system_surface%getsymmlabel())
      case("p4mm")
         call INTERACTION_AENV_OCTA(n,A,pairpot,dampfunc,interac,dvdz_term,dvdx_term,dvdy_term)
      case("p6mm")
         call INTERACTION_AENV_HEXA(n,A,pairpot,dampfunc,interac,dvdz_term,dvdx_term,dvdy_term)
      case default
         write(0,*) "INTERACTION_AENV ERR: wrong surface symmlabel or it's not been implemented yet"
         write(0,*) "Surface reads: ",system_surface%getsymmlabel()
         write(0,*) "Implemented ones: p4mm, p6mm"
         call EXIT(1)      
   end select
   return
end subroutine INTERACTION_AENV
subroutine INTERACTION_AENV_OCTA(n,A,pairpot,dampfunc,interac,dvdz_term,dvdx_term,dvdy_term)
   implicit none
   ! I/O VAriables ------------------------------------------
   integer,intent(in) :: n
   real(kind=8),dimension(3),intent(in) :: A
   type(Pair_pot),intent(in) :: pairpot
   class(Function1d),intent(in) :: dampfunc
   real(kind=8),intent(out) :: interac, dvdz_term, dvdx_term, dvdy_term
   ! Local variables ----------------------------------------
   real(kind=8),dimension(3) :: P
   real(kind=8),dimension(3) :: ghost_A ! A in cartesians, but inside unitcell
   real(kind=8),dimension(2) :: aux
   real(kind=8) :: dummy1, dummy2, dummy3, dummy4
   real(kind=8) :: atomx, atomy
   integer :: pairid
   integer :: i, k ! Counters
   character(len=18), parameter :: routinename = "INTERACTION_AENV: "
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

   select case(n)
      case(0)
         do i=1, system_surface%atomtype(pairid)%n
            P(:)=system_surface%atomtype(pairid)%atom(i,:)
            call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac=interac+dummy1
            dvdz_term=dvdz_term+dummy2
            dvdx_term=dvdx_term+dummy3
            dvdy_term=dvdy_term+dummy4
         end do
         return

      case(1 :)
         do i=1, system_surface%atomtype(pairid)%n
            atomx=system_surface%atomtype(pairid)%atom(i,1)
            atomy=system_surface%atomtype(pairid)%atom(i,2)
            P(3)=system_surface%atomtype(pairid)%atom(i,3)
            do k= -n,n
               aux(1)=dfloat(n)
               aux(2)=dfloat(k)
               aux=system_surface%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
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
               call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
            end do
            do k= -n+1, n-1
               aux(1)=dfloat(k)
               aux(2)=dfloat(n)
               aux=system_surface%surf2cart(aux)
               P(1) =atomx+aux(1)
               P(2) =atomy+aux(2)
               call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
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
               call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
            end do
         end do
         return

      case default
         write(0,*) "INTERACTION_AENV ERR: Wrong environment order."
         call EXIT(1)
   end select
end subroutine INTERACTION_AENV_OCTA
subroutine INTERACTION_AENV_HEXA(n,A,pairpot,dampfunc,interac,dvdz_term,dvdx_term,dvdy_term)
   implicit none
   ! I/O VAriables ------------------------------------------
   integer,intent(in) :: n
   real(kind=8),dimension(3),intent(in) :: A
   type(Pair_pot),intent(in) :: pairpot
   class(Function1d),intent(in) :: dampfunc
   real(kind=8),intent(out) :: interac, dvdz_term, dvdx_term, dvdy_term
   ! Local variables ----------------------------------------
   real(kind=8),dimension(3) :: P
   real(kind=8),dimension(3) :: ghost_A ! A in cartesians, but inside unitcell
   real(kind=8),dimension(2) :: aux
   real(kind=8) :: dummy1, dummy2, dummy3, dummy4
   real(kind=8) :: atomx, atomy
   integer :: pairid
   integer :: i, k ! Counters
   character(len=18), parameter :: routinename = "INTERACTION_AENV: "
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
   select case(n)
      case(0)
         do i=1, system_surface%atomtype(pairid)%n
            P(:)=system_surface%atomtype(pairid)%atom(i,:)
            call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac=interac+dummy1
            dvdz_term=dvdz_term+dummy2
            dvdx_term=dvdx_term+dummy3
            dvdy_term=dvdy_term+dummy4
         end do
         return
      case(1)
         do i=1, system_surface%atomtype(pairid)%n
            atomx=system_surface%atomtype(pairid)%atom(i,1)
            atomy=system_surface%atomtype(pairid)%atom(i,2)
            P(3)=system_surface%atomtype(pairid)%atom(i,3)
            k=0
            !
            aux(1)=dfloat(n)
            aux(2)=dfloat(-k)
            aux=system_surface%surf2cart(aux)
            P(1) =atomx+aux(1)
            P(2) =atomy+aux(2)
            call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
            !
            aux(1)=dfloat(-n)
            aux(2)=dfloat(k)
            aux=system_surface%surf2cart(aux)
            P(1)=atomx+aux(1)
            P(2)=atomy+aux(2)
            call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
            !
            aux(1)=dfloat(-k)
            aux(2)=dfloat(n)
            aux=system_surface%surf2cart(aux)
            P(1)=atomx+aux(1)
            P(2)=atomy+aux(2)
            call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
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
            call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
            !
            aux(1)=dfloat(-n)
            aux(2)=dfloat(n)
            aux=system_surface%surf2cart(aux)
            P(1)=atomx+aux(1)
            P(2)=atomy+aux(2)
            call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
            !
            aux(1)=dfloat(n)
            aux(2)=dfloat(-n)
            aux=system_surface%surf2cart(aux)
            P(1)=atomx+aux(1)
            P(2)=atomy+aux(2)
            call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
         end do
         return
      case(2 :)
         do i=1, system_surface%atomtype(pairid)%n
            atomx=system_surface%atomtype(pairid)%atom(i,1)
            atomy=system_surface%atomtype(pairid)%atom(i,2)
            P(3)=system_surface%atomtype(pairid)%atom(i,3)
            do k= 1,n-1
               aux(1)=dfloat(k)
               aux(2)=dfloat(n-k)
               aux=system_surface%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac=interac+dummy1
               dvdz_term=dvdz_term+dummy2
               dvdx_term=dvdx_term+dummy3
               dvdy_term=dvdy_term+dummy4
               !
               aux(1)=dfloat(-k)
               aux(2)=dfloat(k-n)
               aux=system_surface%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
            end do
            do k= 0, n-1
               aux(1)=dfloat(n)
               aux(2)=dfloat(-k)
               aux=system_surface%surf2cart(aux)
               P(1) =atomx+aux(1)
               P(2) =atomy+aux(2)
               call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
               !
               aux(1)=dfloat(-n)
               aux(2)=dfloat(k)
               aux=system_surface%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
               !
               aux(1)=dfloat(-k)
               aux(2)=dfloat(n)
               aux=system_surface%surf2cart(aux)
               P(1)=atomx+aux(1)
               P(2)=atomy+aux(2)
               call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
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
               call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
               interac = interac+dummy1
               dvdz_term = dvdz_term+dummy2
               dvdx_term = dvdx_term+dummy3
               dvdy_term = dvdy_term+dummy4
            end do
            !
            aux(1)=dfloat(-n)
            aux(2)=dfloat(n)
            aux=system_surface%surf2cart(aux)
            P(1)=atomx+aux(1)
            P(2)=atomy+aux(2)
            call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
            !
            aux(1)=dfloat(n)
            aux(2)=dfloat(-n)
            aux=system_surface%surf2cart(aux)
            P(1)=atomx+aux(1)
            P(2)=atomy+aux(2)
            call INTERACTION_AP(ghost_A,P,pairpot,dampfunc,dummy1,dummy2,dummy3,dummy4)
            interac = interac+dummy1
            dvdz_term = dvdz_term+dummy2
            dvdx_term = dvdx_term+dummy3
            dvdy_term = dvdy_term+dummy4
         end do
         return

      case default
         write(0,*) "INTERACTION_AENV ERR: Wrong environment order."
         call EXIT(1)
   end select
end subroutine INTERACTION_AENV_HEXA
!############################################################
!# SUBROUTINE: GET_V_AND_DERIVS_CRP3D #########################
!############################################################
!> @brief
!! Subroutine that calculates the 3D potential for a point A and
!! its derivatives in cartesian coordinates.
!
!> @param[in] this - CRP3D PES
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
subroutine GET_V_AND_DERIVS_CRP3D(this,X,v,dvdu,errCode)
   implicit none
   ! I/O variables
   class(CRP3D),target,intent(in):: this
   real(kind=8),dimension(:),intent(in):: X
   real(kind=8),intent(out):: v
   real(kind=8),dimension(:),intent(out):: dvdu
   integer(kind=1),optional,intent(out):: errCode
   ! Local variables
   integer(kind=4):: nsites,npairpots
   class(Fourier2d),allocatable:: interpolxy
   real(kind=8),dimension(:),allocatable:: pot,dvdz,dvdx,dvdy
   real(kind=8),dimension(:),allocatable:: potarr
   real(kind=8),dimension(:,:),allocatable:: f,derivarr ! arguments to the xy interpolation
   real(kind=8),dimension(:,:),allocatable:: xy ! arguments to the xy interpolation
   integer :: i ! counters
   ! Pointers
	real(kind=8),pointer:: zmax
   character(len=*),parameter:: routinename="GET_V_AND_DERIVS_CRP3D: "
   zmax => this%all_sites(1)%z(this%all_sites(1)%n)
   npairpots = size(this%all_pairpots)
   nsites = size(this%all_sites)
   ! GABBA, GABBA HEY! ----------------------
   select case(X(3)>zmax)
      case(.true.)
         dvdu = (/0.D0,0.D0,0.D0/)
         v = 0.D0
         return
      case(.false.)
         ! do nothing
   end select
   ! Initialization section
   allocate(pot(npairpots))
   allocate(dvdz(npairpots))
   allocate(dvdx(npairpots))
   allocate(dvdy(npairpots))
   call this%GET_REPUL_CORRECTIONS(X,pot,dvdz,dvdx,dvdy)
   v=sum(pot)
   dvdu(1)=sum(dvdx)
   dvdu(2)=sum(dvdy)
   dvdu(3)=sum(dvdz)
   ! Now, we have all the repulsive interaction and corrections to the derivarives
   ! stored in v(:) and dvdu(:) respectively.
   ! Let's get v and derivatives from xy interpolation of the corrugationless function
   ! f(1,i) ==> v values interpolated for Z at site "i"
   ! f(2,i) ==> dvdz values interpolated for Z at site "i"
   allocate(f(2,nsites))
   allocate(xy(nsites,2))
   ! potarr(1) ==> v interpolated at X,Y for a given Z
   ! potarr(2) ==> dvdz interpolated at X,Y for a given Z
   allocate(potarr(2))
   ! derivarr(1,1:2) ==> dvdx, dvdy interpolated at X,Y
   ! derivarr(2,1:2) ==> d(dvdz)/dx, d(dvdz)/dy interpolated at X,Y. This is not needed
   allocate(derivarr(2,2))
   do i=1,nsites
      xy(i,1)=this%all_sites(i)%x
      xy(i,2)=this%all_sites(i)%y
      call this%all_sites(i)%interz%GET_V_AND_DERIVS(X(3),f(1,i),f(2,i))
   end do
   call this%SET_FOURIER_SYMMETRY(interpolxy)
   call interpolxy%READ( xy=xy,f=f,kList=this%klist,irrepList=this%irrepList(:),parityList=this%parityList )
   call interpolXY%initializeTerms()
   call interpolxy%INTERPOL()
   call interpolxy%GET_F_AND_DERIVS(X,potarr,derivarr)
   ! Corrections from the smoothing procedure
   v = v + potarr(1)
   dvdu(1)=dvdu(1)+derivarr(1,1)
   dvdu(2)=dvdu(2)+derivarr(1,2)
   dvdu(3)=dvdu(3)+potarr(2)
   call interpolXY%cleanAll()
   return
end subroutine GET_V_AND_DERIVS_CRP3D
!############################################################
!# SUBROUTINE: GET_V_AND_DERIVS_CORRECTION_CRP3D ############
!############################################################
!> @brief
!! Subroutine that calculates the correction to the 3D PES for a point A and
!! its derivatives in cartesian coordinates.
!
!> @param[in] this - CRP3D PES
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
subroutine GET_V_AND_DERIVS_CORRECTION_CRP3D(this,X,v,dvdu)
   implicit none
   ! I/O variables
   class(CRP3D),target,intent(in) :: this
   real(kind=8),dimension(3), intent(in) :: X
   real(kind=8),dimension(:),intent(out) :: v
   real(kind=8),dimension(3),intent(out) :: dvdu
   ! Local variables
   integer(kind=4) :: npairpots
   real(kind=8),dimension(:),allocatable :: pot,dvdz,dvdx,dvdy
   integer :: i ! counters
   ! Pointers
	real(kind=8), pointer :: zmax
   character(len=24),parameter :: routinename="GET_V_AND_DERIVS_CRP3D: "
   zmax => this%all_sites(1)%interz%x(this%all_sites(1)%n)
   npairpots = size(this%all_pairpots)
   select case(size(v)==npairpots+1)
      case(.true.)
         ! do nothing
      case(.false.)
         write(*,*) "GET_V_AND_DERIVS_CORRECTION_CRP3D ERR: wrong number of dimensions array v"
         call EXIT(1)
   end select
   ! GABBA, GABBA HEY! ----------------------
   select case(X(3)>zmax)
      case(.true.)
         dvdu = (/0.D0,0.D0,0.D0/)
         forall(i=1:npairpots+1) v(i) = 0.D0
         return
      case(.false.)
         ! do nothing
   end select
   ! Initializing variables
   allocate(pot(npairpots))
   allocate(dvdz(npairpots))
   allocate(dvdx(npairpots))
   allocate(dvdy(npairpots))
   forall(i=1:npairpots+1) v(i) = 0.D0
   ! Compute
   call this%GET_REPUL_CORRECTIONS(X,pot,dvdz,dvdx,dvdy)
   v(1)=sum(pot)
   forall(i=2:npairpots+1) v(i)=pot(i-1)
   dvdu(1)=sum(dvdx)
   dvdu(2)=sum(dvdy)
   dvdu(3)=sum(dvdz)
   return
end subroutine GET_V_AND_DERIVS_CORRECTION_CRP3D
!############################################################
!# SUBROUTINE: GET_V_AND_DERIVS_SMOOTH_CRP3D #########################
!############################################################
!> @brief
!! Subroutine that calculates the 3D potential for a point A and
!! its derivatives in cartesian coordinates without any correction to
!! the smooth potential
!
!> @param[in] this - CRP3D PES
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
!> @date Jun/2014
!> @version 1.0
!------------------------------------------------------------
subroutine GET_V_AND_DERIVS_SMOOTH_CRP3D(this,X,v,dvdu)
   implicit none
   ! I/O variables
   class(CRP3D),target,intent(in) :: this
   real(kind=8),dimension(3), intent(in) :: X
   real(kind=8),intent(out) :: v
   real(kind=8),dimension(3),intent(out) :: dvdu
   ! Local variables
   class(Fourier2d),allocatable:: interpolxy
   integer(kind=4) :: nsites,npairpots
   real(kind=8),dimension(3) :: deriv
   real(kind=8),dimension(:),allocatable :: potarr
   real(kind=8),dimension(:,:),allocatable :: f,derivarr ! arguments to the xy interpolation
   real(kind=8),dimension(:,:),allocatable :: xy ! arguments to the xy interpolation
   integer :: i ! counters
   ! Pointers
	real(kind=8), pointer :: zmax
   character(len=24),parameter :: routinename="GET_V_AND_DERIVS_CRP3D: "
   zmax => this%all_sites(1)%interz%x(this%all_sites(1)%n)
   npairpots = size(this%all_pairpots)
   nsites = size(this%all_sites)
   ! GABBA, GABBA HEY! ----------------------
   if (X(3).gt.zmax) then
      dvdu = 0.D0
      v = 0.D0
      return
   end if
   !
   v = 0.D0 ! Initialization value
   forall(i=1:3) 
      dvdu(i) = 0.D0 ! Initialization value
      deriv(i) = 0.D0 ! Initialization value
   end forall 
   ! Let's get v and derivatives from xy interpolation of the corrugationless function
   ! f(1,i) ==> v values interpolated for Z at site "i"
   ! f(2,i) ==> dvdz values interpolated for Z at site "i"
   allocate(f(2,nsites))
   allocate(xy(nsites,2))
   ! potarr(1) ==> v interpolated at X,Y for a given Z
   ! potarr(2) ==> dvdz interpolated at X,Y for a given Z
   allocate(potarr(2))
   ! derivarr(1,1:2) ==> dvdx, dvdy interpolated at X,Y
   ! derivarr(2,1:2) ==> d(dvdz)/dx, d(dvdz)/dy interpolated at X,Y. It is not needed
   allocate(derivarr(2,2))
   do i=1,nsites
      xy(i,1)=this%all_sites(i)%x
      xy(i,2)=this%all_sites(i)%y
      call this%all_sites(i)%interz%GET_V_AND_DERIVS(X(3),f(1,i),f(2,i))
   end do
   call this%SET_FOURIER_SYMMETRY(interpolxy)
   call interpolxy%READ( xy=xy,f=f,kList=this%klist,irrepList=this%irrepList(:),parityList=this%parityList )
   call interpolXY%initializeTerms()
   call interpolxy%INTERPOL()
   call interpolxy%GET_F_AND_DERIVS(X,potarr,derivarr)
   ! Corrections from the smoothing procedure
   v = v + potarr(1)
   dvdu(1)=dvdu(1)+derivarr(1,1)
   dvdu(2)=dvdu(2)+derivarr(1,2)
   dvdu(3)=dvdu(3)+potarr(2)
   return
end subroutine GET_V_AND_DERIVS_SMOOTH_CRP3D
!############################################################
!# FUNCTION: getpot_crp3d ###################################
!############################################################
!> @brief
!! Subroutine that calculates the 3D potential for a point A
!
!> @param[in] this - CRP3D PES
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
real(kind=8) function getpot_crp3d(this,X)
   implicit none
   ! I/O variables
   class(CRP3D),target,intent(in) :: this
	real(kind=8),dimension(3), intent(in) :: X
   ! Local variables
   class(Fourier2d),allocatable:: interpolxy
   real(kind=8):: v
   integer(kind=4) :: nsites,npairpots
   real(kind=8),dimension(:),allocatable :: pot,dummy
   real(kind=8),dimension(:),allocatable :: potarr
   real(kind=8),dimension(:,:),allocatable :: f,derivarr ! arguments to the xy interpolation
   real(kind=8),dimension(:,:),allocatable :: xy ! arguments to the xy interpolation
   integer :: i ! counters
   ! Pointers
   real(kind=8), pointer :: zmax
   character(len=24),parameter :: routinename="GET_V_AND_DERIVS_CRP3D: "
   zmax => this%all_sites(1)%interz%x(this%all_sites(1)%n)
   npairpots = size(this%all_pairpots)
   nsites = size(this%all_sites)
   ! GABBA, GABBA HEY! ----------------------
   allocate(pot(npairpots))
   allocate(dummy(npairpots))
   select case(X(3)>zmax)
      case(.true.)
         getpot_crp3d=0.D0
         return
      case(.false.)
         ! do nothing
   end select
   !
   call this%GET_REPUL_CORRECTIONS(X,pot,dummy,dummy,dummy)
   v=sum(pot)
   ! Now, we have all the repulsive interaction and corrections to the derivarives
   ! stored in v(:) and dvdu(:) respectively.
   ! Let's get v and derivatives from xy interpolation of the corrugationless function
   allocate(f(2,nsites))
   allocate(xy(nsites,2))
   allocate(potarr(2))
   allocate(derivarr(2,2))
   do i=1,nsites
      xy(i,1)=this%all_sites(i)%x
      xy(i,2)=this%all_sites(i)%y
      call this%all_sites(i)%interz%GET_V_AND_DERIVS(X(3),f(1,i),f(2,i))
   end do
   call this%SET_FOURIER_SYMMETRY(interpolxy)
   call interpolxy%READ( xy=xy,f=f,kList=this%klist,irrepList=this%irrepList(:),parityList=this%parityList )
   call interpolXY%initializeTerms()
   call interpolxy%INTERPOL()
   call interpolxy%GET_F_AND_DERIVS(X,potarr,derivarr)
   ! Corrections from the smoothing procedure
   getpot_crp3d=v+potarr(1)
   return
end function getpot_crp3d
!######################################################################
! SUBROUTINE: PLOT_DATA_SYMMPOINT #####################################
!######################################################################
!> @brief
!! Creates a data file called "filename" with the data z(i) and v(i) 
!! for a Symmpoint type variable
!----------------------------------------------------------------------
subroutine PLOT_DATA_SYMMPOINT(this,filename)
   implicit none
   ! I/O Variables ----------------------
   class(Symmpoint), intent(in) :: this
   character(len=*), intent(in) :: filename
   ! Local variables --------------------
   integer :: i ! Counter
   ! GABBA GABBA HEY!! ------------------
   open(11, file=filename, status="replace")
   do i=1,this%n
      write(11,*) this%z(i),this%v(i)
   end do
   close (11)
   write(*,*) "PLOT_DATA_SYMMPOINT: ",this%alias,filename," file created"
end subroutine PLOT_DATA_SYMMPOINT
!#######################################################################
! SUBROUTINE: PLOT_XYMAP_CRP3D
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (X,Y) of the PES. 
!
!> @param[in] this - CRP3D PES used
!> @param[in] filename - Name of the file to print the output
!> @param[in] init_xyz - Initial position to start the scan (a.u.)
!> @param[in] nxpoints - Number of points in X axis (auxiliar cartesian coordinates)
!> @param[in] nypoints - Number of points in Y axis (auxiliar cartesian coordinates)
!> @param[in] Lx - Length of X axis (a.u.)
!> @param[in] Ly - Length of Y axis (a.u.)
!
!> @warning
!! - Z parameter is taken from @b init_xyz and is constant
!
!> @author A.S. Muzas
!> @date 08/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
subroutine PLOT_XYMAP_CRP3D(this,filename,init_xyz,nxpoints,nypoints,Lx,Ly)
   implicit none
   class(CRP3D),intent(in) :: this
   real*8,dimension(3),intent(in) :: init_xyz ! Initial position to start the scan (in a.u.)
   integer,intent(in) :: nxpoints, nypoints ! number of points in XY plane
   character(len=*),intent(in) :: filename ! filename
   real*8,intent(in) :: Lx ! Length of X axis 
   real*8,intent(in) :: Ly ! Length of X axis 
   ! Local variables
   real*8 :: xmin, ymin, xmax, ymax, z
   real*8, dimension(3) :: r, dvdu
   real*8 :: xdelta, ydelta
   integer :: xinpoints, nxdelta
   integer :: yinpoints, nydelta
   integer :: i, j ! counters
   real*8 :: v ! potential
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
   open(11,file=filename,status="replace")
   r(1) = xmin
   r(2) = ymin
   r(3) = z
   call this%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(1),r(2),v,dvdu(:)
   do i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      call this%GET_V_AND_DERIVS(r,v,dvdu)
      write(11,*) r(1), r(2), v,dvdu(:)
   end do
   r(2) = ymax
   call this%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(1), r(2), v,dvdu(:)
   ! inpoints in XY
   do i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      r(3) = z
      call this%GET_V_AND_DERIVS(r,v,dvdu)
      write(11,*) r(1), r(2), v,dvdu(:)
      do j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         call this%GET_V_AND_DERIVS(r,v,dvdu)
         write(11,*) r(1), r(2), v,dvdu(:)
      end do
      r(2) = ymax
      call this%GET_V_AND_DERIVS(r,v,dvdu)
      write(11,*) r(1), r(2), v,dvdu(:)
   end do
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   r(3) = z
   call this%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(1), r(2), v,dvdu(:)
   do i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      call this%GET_V_AND_DERIVS(r,v,dvdu)
      write(11,*) r(1), r(2), v,dvdu(:)
   end do
   r(2) = ymax
   call this%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(1), r(2), v,dvdu(:)
   close(11)
   return
end subroutine PLOT_XYMAP_CRP3D
!#######################################################################
! SUBROUTINE: PLOT_XYMAP_CORRECTION_CRP3D
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (X,Y) of the corrections
!! to the PES
!
!> @param[in] this - CRP3D PES used
!> @param[in] filename - Name of the file to print the output
!> @param[in] init_xyz - Initial position to start the scan (a.u.)
!> @param[in] nxpoints - Number of points in X axis (auxiliar cartesian coordinates)
!> @param[in] nypoints - Number of points in Y axis (auxiliar cartesian coordinates)
!> @param[in] Lx - Length of X axis (a.u.)
!> @param[in] Ly - Length of Y axis (a.u.)
!
!> @warning
!! - Z parameter is taken from @b init_xyz and is constant
!
!> @author A.S. Muzas
!> @date 08/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
subroutine PLOT_XYMAP_CORRECTION_CRP3D(this,filename,init_xyz,nxpoints,nypoints,Lx,Ly)
   implicit none
   class(CRP3D),intent(in) :: this
   real*8,dimension(3),intent(in) :: init_xyz ! Initial position to start the scan (in a.u.)
   integer,intent(in) :: nxpoints, nypoints ! number of points in XY plane
   character(len=*),intent(in) :: filename ! filename
   real*8,intent(in) :: Lx ! Length of X axis 
   real*8,intent(in) :: Ly ! Length of X axis 
   ! Local variables
   real*8 :: xmin, ymin, xmax, ymax, z
   real*8, dimension(3) :: r, dvdu
   real*8 :: xdelta, ydelta
   integer :: xinpoints, nxdelta, npairpots
   integer :: yinpoints, nydelta
   integer :: i, j ! counters
   real*8,dimension(:),allocatable :: v ! potential
	! GABBA, GABBA HEY! ---------
   npairpots=size(this%all_pairpots)
   allocate(v(npairpots+1))
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
   open(11,file=filename,status="replace")
   r(1) = xmin
   r(2) = ymin
   r(3) = z
   call this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
   write(11,*) r(1), r(2), v(:)
   do i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      call this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
      write(11,*) r(1), r(2), v(:)
   end do
   r(2) = ymax
   call this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
   write(11,*) r(1), r(2), v(:)
   ! inpoints in XY
   do i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      r(3) = z
      call this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
      write(11,*) r(1), r(2), v(:)
      do j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         call this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
         write(11,*) r(1), r(2), v(:)
      end do
      r(2) = ymax
      call this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
      write(11,*) r(1), r(2), v(:)
   end do
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   r(3) = z
   call this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
   write(11,*) r(1), r(2), v(:)
   do i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      call this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
      write(11,*) r(1), r(2), v(:)
   end do
   r(2) = ymax
   call this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
   write(11,*) r(1), r(2), v(:)
   close(11)
   return
end subroutine PLOT_XYMAP_CORRECTION_CRP3D
!#######################################################################
! SUBROUTINE: PLOT_XYMAP_SMOOTH_CRP3D
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 2D cut (X,Y) of the PES without
!! corrections to the potential
!
!> @param[in] this - CRP3D PES used
!> @param[in] filename - Name of the file to print the output
!> @param[in] init_xyz - Initial position to start the scan (a.u.)
!> @param[in] nxpoints - Number of points in X axis (auxiliar cartesian coordinates)
!> @param[in] nypoints - Number of points in Y axis (auxiliar cartesian coordinates)
!> @param[in] Lx - Length of X axis (a.u.)
!> @param[in] Ly - Length of Y axis (a.u.)
!
!> @warning
!! - Z parameter is taken from @b init_xyz and is constant
!
!> @author A.S. Muzas
!> @date 08/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
subroutine PLOT_XYMAP_SMOOTH_CRP3D(this,filename,init_xyz,nxpoints,nypoints,Lx,Ly)
   implicit none
   class(CRP3D),intent(in) :: this
   real*8,dimension(3),intent(in) :: init_xyz ! Initial position to start the scan (in a.u.)
   integer,intent(in) :: nxpoints, nypoints ! number of points in XY plane
   character(len=*),intent(in) :: filename ! filename
   real*8,intent(in) :: Lx ! Length of X axis 
   real*8,intent(in) :: Ly ! Length of X axis 
   ! Local variables
   real*8 :: xmin, ymin, xmax, ymax, z
   real*8, dimension(3) :: r, dvdu
   real*8 :: xdelta, ydelta
   integer :: xinpoints, nxdelta
   integer :: yinpoints, nydelta
   integer :: i, j ! counters
   real*8 :: v ! potential
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
   open(11,file=filename,status="replace")
   r(1) = xmin
   r(2) = ymin
   r(3) = z
   call this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(1), r(2), v
   do i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      call this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      write(11,*) r(1), r(2), v
   end do
   r(2) = ymax
   call this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(1), r(2), v
   ! inpoints in XY
   do i = 1, xinpoints
      r(1) = xmin+DFLOAT(i)*xdelta
      r(2) = ymin
      r(3) = z
      call this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      write(11,*) r(1), r(2), v
      do j = 1, yinpoints
         r(2) = ymin + DFLOAT(j)*ydelta
         call this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
         write(11,*) r(1), r(2), v
      end do
      r(2) = ymax
      call this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      write(11,*) r(1), r(2), v
   end do
   ! Last point in XY plane
   r(1) = xmax
   r(2) = ymax
   r(3) = z
   call this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(1), r(2), v
   do i =1, yinpoints
      r(2) = ymin + DFLOAT(i)*ydelta
      call this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      write(11,*) r(1), r(2), v
   end do
   r(2) = ymax
   call this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(1), r(2), v
   close(11)
   return
end subroutine PLOT_XYMAP_SMOOTH_CRP3D
!#######################################################################
! SUBROUTINE: PLOT_DIRECTION1D_CRP3D ###################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. To define 
!! the direction, the angle alpha is given. 
!
!> @param[in] this - CRP3D PES used
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
subroutine PLOT_DIRECTION1D_CRP3D(this,filename,npoints,angle,z,L)
   implicit none
   ! I/O variables -------------------------------
   class(CRP3D),intent(in) :: this
   integer, intent(in) :: npoints
   character(len=*), intent(in) :: filename
   real*8, intent(in) :: z, angle
   ! Local variables -----------------------------
   integer :: inpoints, ndelta
   real*8 :: delta,L,v,s,alpha
   real*8 :: xmax, xmin, ymax, ymin 
   real*8, dimension(3) :: r, dvdu
   integer :: i ! Counter
   character(len=24), parameter :: routinename = "PLOT_DIRECTION1D_CRP3D: "
   ! HE HO ! LET'S GO ----------------------------
   if (npoints.lt.2) then
      write(0,*) "PLOT_DIRECTION1D_CRP3D ERR: Less than 2 points"
      call EXIT(1)
   end if
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
   open(11,file=filename,status="replace")
   ! Initial value
   r(1) = xmin
   r(2) = ymin
   s = DSQRT(r(1)**2.D0+r(2)**2.D0)
   call this%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) s, v , dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   ! cycle for inpoints
   do i=1, inpoints
      r(1)=xmin+(DFLOAT(i)*delta)*DCOS(alpha)
      r(2)=ymin+(DFLOAT(i)*delta)*DSIN(alpha)
      s = DSQRT(r(1)**2.D0+r(2)**2.D0)
      call this%GET_V_AND_DERIVS(r,v,dvdu)
      write(11,*) s, v , dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   end do
   ! Final value
   r(1) = xmax
   r(2) = ymax
   s = DSQRT(r(1)**2.D0+r(2)**2.D0)
   call this%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) s, v, dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   write(*,*) routinename, "file created ",filename
   close(11)
end subroutine PLOT_DIRECTION1D_CRP3D
!#######################################################################
! SUBROUTINE: PLOT_DIRECTION1D_CORRECTION_CRP3D ###################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the corrections to the PES. To define 
!! the direction, the angle alpha is given. 
!
!> @param[in] this - CRP3D PES used
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
subroutine PLOT_DIRECTION1D_CORRECTION_CRP3D(this,filename,npoints,angle,z,L)
   implicit none
   ! I/O variables -------------------------------
   class(CRP3D),intent(in) :: this
   integer, intent(in) :: npoints
   character(len=*), intent(in) :: filename
   real*8, intent(in) :: z, angle
   ! Local variables -----------------------------
   integer :: inpoints, ndelta,npairpots
   real*8 :: delta,L,s,alpha
   real*8 :: xmax, xmin, ymax, ymin 
   real*8, dimension(3) :: r, dvdu
   real(kind=8),dimension(:),allocatable :: v
   integer :: i ! Counter
   character(len=24), parameter :: routinename = "PLOT_DIRECTION1D_CRP3D: "
   ! HE HO ! LET'S GO ----------------------------
   if (npoints.lt.2) then
      write(0,*) "PLOT_DIRECTION1D_CRP3D ERR: Less than 2 points"
      call EXIT(1)
   end if
   npairpots=size(this%all_pairpots)
   allocate(v(npairpots+1))
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
   open(11,file=filename,status="replace")
   ! Initial value
   r(1) = xmin
   r(2) = ymin
   s = DSQRT(r(1)**2.D0+r(2)**2.D0)
   call this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
   write(11,*) s, v(:)
   ! cycle for inpoints
   do i=1, inpoints
      r(1)=xmin+(DFLOAT(i)*delta)*DCOS(alpha)
      r(2)=ymin+(DFLOAT(i)*delta)*DSIN(alpha)
      s = DSQRT(r(1)**2.D0+r(2)**2.D0)
      call this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
      write(11,*) s, v(:)
   end do
   ! Final value
   r(1) = xmax
   r(2) = ymax
   s = DSQRT(r(1)**2.D0+r(2)**2.D0)
   call this%GET_V_AND_DERIVS_CORRECTION(r,v,dvdu)
   write(11,*) s, v(:)
   write(*,*) routinename, "file created ",filename
   close(11)
end subroutine PLOT_DIRECTION1D_CORRECTION_CRP3D
!#######################################################################
! SUBROUTINE: PLOT_Z_CRP3D #############################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut along z direction
!
!> @param[in] this - CRP3D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!!                    cut. It should be given in degrees.
!> @param[in] xyz - Initial point
!> @param[in] L - Length of the graphic
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 09/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
subroutine PLOT_Z_CRP3D(this,npoints,xyz,L,filename)
   implicit none
   ! I/O variables -------------------------------
   class(CRP3D),intent(in):: this
   integer,intent(in):: npoints
   character(len=*),intent(in):: filename
   real(kind=8),dimension(3),intent(in):: xyz
   ! Local variables -----------------------------
   integer :: inpoints, ndelta
   real(kind=8):: delta,L
   real(kind=8):: zmax, zmin 
   real(kind=8),dimension(3):: r, dvdu
   real(kind=8):: v
   integer:: i ! Counter
   character(len=*),parameter:: routinename = "PLOT_DIRECTION1D_CRP3D: "
   ! HE HO ! LET'S GO ----------------------------
   if (npoints.lt.2) then
      write(0,*) "PLOT_Z_CRP3D ERR: Less than 2 points"
      call EXIT(1)
   end if
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
   open(11,file=filename,status="replace")
   ! Initial value
   r(3) = zmin
   call this%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(3),v,dvdu(:)
   ! cycle for inpoints
   do i=1, inpoints
      r(3)=zmin+(DFLOAT(i)*delta)
      call this%GET_V_AND_DERIVS(r,v,dvdu)
      write(11,*) r(3),v,dvdu(:)
   end do
   ! Final value
   r(3) = zmax
   call this%GET_V_AND_DERIVS(r,v,dvdu)
   write(11,*) r(3),v,dvdu(:)
   write(*,*) routinename, "file created ",filename
   close(11)
end subroutine PLOT_Z_CRP3D
!#######################################################################
! SUBROUTINE: PLOT_Z_CRP3D #############################################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut along z direction
!
!> @param[in] this - CRP3D PES used
!> @param[in] filename - Name of the output file
!> @param[in] npoints - Number of points in the graphic. npoints>=2
!!                    cut. It should be given in degrees.
!> @param[in] xyz - Initial point
!> @param[in] L - Length of the graphic
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date 09/Feb/2014
!> @version 1.0
!----------------------------------------------------------------------
subroutine plot_z_smooth_CRP3D(this,npoints,xyz,L,filename)
   implicit none
   ! I/O variables -------------------------------
   class(CRP3D),intent(in):: this
   integer,intent(in):: npoints
   character(len=*),intent(in):: filename
   real(kind=8),dimension(3),intent(in):: xyz
   ! Local variables -----------------------------
   integer :: inpoints, ndelta
   real(kind=8):: delta,L
   real(kind=8):: zmax, zmin
   real(kind=8),dimension(3):: r, dvdu
   real(kind=8):: v
   integer:: i ! Counter
   character(len=*),parameter:: routinename = "PLOT_DIRECTION1D_CRP3D: "
   ! HE HO ! LET'S GO ----------------------------
   if (npoints.lt.2) then
      write(0,*) "PLOT_Z_CRP3D ERR: Less than 2 points"
      call EXIT(1)
   end if
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
   open(11,file=filename,status="replace")
   ! Initial value
   r(3) = zmin
   call this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(3),v,dvdu(:)
   ! cycle for inpoints
   do i=1, inpoints
      r(3)=zmin+(DFLOAT(i)*delta)
      call this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      write(11,*) r(3),v,dvdu(:)
   end do
   ! Final value
   r(3) = zmax
   call this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) r(3),v,dvdu(:)
   write(*,*) routinename, "file created ",filename
   close(11)
end subroutine plot_z_smooth_CRP3D
!#######################################################################
! SUBROUTINE: PLOT_DIRECTION1D_SMOOTH_CRP3D ############################
!#######################################################################
!> @brief
!! Creates a file with name "filename" with a 1D cut of the PES. To define 
!! the direction, the angle alpha is given. There's not any correction to
!! PES values.
!
!> @param[in] this - CRP3D PES used
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
subroutine PLOT_DIRECTION1D_SMOOTH_CRP3D(this,filename,npoints,angle,z,L)
   implicit none
   ! I/O variables -------------------------------
   class(CRP3D),intent(in) :: this
   integer, intent(in) :: npoints
   character(len=*), intent(in) :: filename
   real*8, intent(in) :: z, angle
   ! Local variables -----------------------------
   integer :: inpoints, ndelta
   real*8 :: delta,L,v,s,alpha
   real*8 :: xmax, xmin, ymax, ymin 
   real*8, dimension(3) :: r, dvdu
   integer :: i ! Counter
   character(len=24), parameter :: routinename = "PLOT_DIRECTION1D_CRP3D: "
   ! HE HO ! LET'S GO ----------------------------
   if (npoints.lt.2) then
      write(0,*) "PLOT_DIRECTION1D_CRP3D ERR: Less than 2 points"
      call EXIT(1)
   end if
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
   open(11,file=filename,status="replace")
   ! Initial value
   r(1) = xmin
   r(2) = ymin
   s = DSQRT(r(1)**2.D0+r(2)**2.D0)
   call this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) s, v , dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   ! cycle for inpoints
   do i=1, inpoints
      r(1)=xmin+(DFLOAT(i)*delta)*DCOS(alpha)
      r(2)=ymin+(DFLOAT(i)*delta)*DSIN(alpha)
      s = DSQRT(r(1)**2.D0+r(2)**2.D0)
      call this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
      write(11,*) s, v , dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   end do
   ! Final value
   r(1) = xmax
   r(2) = ymax
   s = DSQRT(r(1)**2.D0+r(2)**2.D0)
   call this%GET_V_AND_DERIVS_SMOOTH(r,v,dvdu)
   write(11,*) s, v, dvdu(1), dvdu(2), DCOS(alpha)*dvdu(1)+DSIN(alpha)*dvdu(2)
   write(*,*) routinename, "file created ",filename
   close(11)
end subroutine PLOT_DIRECTION1D_SMOOTH_CRP3D
!###########################################################
!# SUBROUTINE: PLOT_INTERPOL_SYMMPOINT 
!###########################################################
!> @brief
!! Plots a Z scan along a symmpoint once it's been interpolated
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine PLOT_INTERPOL_SYMMPOINT(this,npoints,filename)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(Symmpoint),intent(in):: this
   integer(kind=4),intent(in) :: npoints
   character(len=*),intent(in) :: filename
   ! Run section
   call this%interz%PLOT(npoints,filename)
   return
end subroutine PLOT_INTERPOL_SYMMPOINT
!###########################################################
!# SUBROUTINE: PLOT_PAIRPOTS_CRP3D 
!###########################################################
!> @brief
!! Plots all pairpot potentials
!
!> @warning
!! - It only works if there are less than 9 pairpotentials defined
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine PLOT_PAIRPOTS_CRP3D(this,npoints)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP3D),intent(in)::this
   integer(kind=4),intent(in) :: npoints
   ! Local variables
   integer(kind=4) :: i ! counters
   integer(kind=4) :: npairpots
   character(len=12) :: stringbase
   character(len=13) :: filename
   ! Run section
   npairpots=size(this%all_pairpots)
   write(stringbase,'(A12)') "-pairpot.dat"
   do i = 1, npairpots
      write(filename,'(I1,A12)') i,stringbase
      call this%all_pairpots(i)%PLOT(npoints,filename)
   end do
   return
end subroutine PLOT_PAIRPOTS_CRP3D
!###########################################################
!# SUBROUTINE: PLOT_SITIOS_CRP3D 
!###########################################################
!> @brief
!! Plots all sitio potentials
!
!> @warning
!! - It only works if there are less than 9 sitio potentials defined
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Jun/2014
!> @version 1.0
!-----------------------------------------------------------
subroutine PLOT_SITIOS_CRP3D(this,npoints)
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP3D),intent(in)::this
   integer(kind=4),intent(in) :: npoints
   ! Local variables
   integer(kind=4) :: i ! counters
   integer(kind=4) :: nsitios
   character(len=10) :: stringbase
   character(len=11) :: filename
   ! Run section
   nsitios=size(this%all_sites)
   write(stringbase,'(A10)') "-sitio.dat"
   do i = 1, nsitios
      write(filename,'(I1,A10)') i,stringbase
      call this%all_sites(i)%PLOT(npoints,filename)
   end do
   return
end subroutine PLOT_SITIOS_CRP3D
!###########################################################
!# FUNCTION: is_allowed_CRP3D
!###########################################################
!> @brief
!! Determines if the potential can be calculated in this point
!-----------------------------------------------------------
logical function is_allowed_CRP3D(this,x) 
   ! Initial declarations   
   implicit none
   ! I/O variables
   class(CRP3D),intent(in) :: this
   real(kind=8),dimension(:),intent(in) :: x
   ! Local variables
   real(kind=8) :: xmin,xmax
   ! Run section
   xmin=this%all_sites(1)%z(1)
   xmax=this%all_sites(1)%z(this%all_sites(1)%n)
   select case(size(x)/=3)
      case(.true.)
         write(0,*) "is_allowed_CRP3D ERR: array doesn't have 3 dimensions: 3"
         call EXIT(1)
      case(.false.)
         ! do nothing
   end select
   select case( x(3)<xmin )
      case(.true.)
         is_allowed_CRP3D=.false.
      case(.false.)
         is_allowed_CRP3D=.true.
   end select
   return
end function is_allowed_CRP3D
end module CRP3D_MOD
