!######################################################################
! MODULE: INPUTGEN_MOD
!
!> @brief
!! This module provides tools to generate "standard" inputs from a set 
!! of "raw" inputs
!
!######################################################################
module NEWINPUT3D_MOD
! Initial declarations
#ifdef DEBUG
   use DEBUG_MOD
#endif
use AOTUS_MODULE, only: flu_State, open_config_file, close_config, aot_get_val
use AOT_TABLE_MODULE, only: aot_table_open, aot_table_close, aot_table_length, aot_table_get_val
use UNITS_MOD, only: Length, Energy
use MATHS_MOD, only: order_vect
use CRP3D_MOD, only: SymmPoint
implicit none
!/////////////////////////////////////////////////////////////////////
! TYPE: NewInput3d
!
!> @brief
!! Stores all data needed to generate a new set of standard input files
!
!> @param vasint - Potential when the atom is far from the surface (in the vacuum)
!> @param nRumpling - Number of different rumplings defined
!> @param nFiles - Number of total input files
!> @param rumpling - Storage of rumplings
!> @param zeroPos - Position of ith rumpling in the grid
!> @param nzgrid - Number of points in zgrid 
!> @param zgrid - Grid proposed to generate input files
!
!> @author A.S. Muzas
!> @date 30/Jan/2014
!> @version 1.0
!---------------------------------------------------------------------
type NewInput3d
   ! private section
   real(kind=8),private:: vasint
   real(kind=8),private:: dfin
   integer(kind=4),private:: nRumpling
   integer(kind=4),private:: nFiles
   real(kind=8),dimension(:),allocatable,private:: rumpling
   integer,dimension(:),allocatable,private:: zeropos
   integer(kind=4),private:: nzgrid
   real(kind=8),dimension(:),allocatable,private:: zgrid
   character(len=50),dimension(:),allocatable,private:: inputSitioList
   character(len=50),dimension(:),allocatable,private:: inputPairpotList
   character(len=50),dimension(:),allocatable,private:: outputSitioList
   character(len=50),dimension(:),allocatable,private:: outputPairpotList
   type(SymmPoint),dimension(:),allocatable,private:: protoSitios
   type(SymmPoint),dimension(:),allocatable,private:: protoPairpots
   ! public section
   contains
      procedure,public:: initialize => initialize_NewInput3d
      procedure,public:: printStatus => printStatus_NewInput3d
      procedure,public:: printPairpot => printPairpot_NewInput3d
      procedure,public:: printSitio => printSitio_NewInput3d
end type NewInput3d
!/////////////////////////////////////////////////////////////////////
contains
!##################################################################
! SUBROUTINE: READ_NEWINPUT
!
!> @brief
!! Read information from file "fileName" to create an standard input
!! generation job
!
!> @param this - Newinput type variable
!> @param fileName - Lua input file to be read
!
!> @author A.S. Muzas
!> @date 30/Jan/2014, 02/May/2016
!> @version 1.0, 2.0
!
!> @warning
!! - In this routine Newinput%dfin is set to zero, which is a good value for potentials far from the surface
!! - There is no need to have an ordered grid, this subroutine will order it.
!-----------------------------------------------------------------
subroutine initialize_NewInput3d(this,fileName)
   implicit none
   ! I/O Variables -----------------
   class(NewInput3d),intent(out):: this
   character(len=*),intent(in):: fileName
   ! Local variables
   type(flu_State):: conf ! Lua file manager
   integer(kind=4):: iErr ! code for Lua errors
   integer(kind=4):: tab_newInput,tab_magnitude,tab_rumpl,tab_gridInfo,tab_grid ! Lua table handlers
   integer(kind=4):: tab_files, tab_subFile ! Lua table handlers
   real(kind=8):: auxReal
   character(len=1012):: auxString,auxString2
   integer(kind=4):: i,k ! counters
	type(Length) :: dist
   type(Energy) :: en
	logical :: exist_zero
   integer :: nbefore, nafter !
   real*8, dimension(:), allocatable :: aux1, aux2
	character(len=15), parameter :: routineName = "READ_NEWINPUT: "
   ! Run section -------------------------------------- -----------------
   ! Open Lua conf. file
   call open_config_file(L=conf,errCode=iErr,fileName=trim(fileName))
   ! Open new Input table
   call aot_table_open(L=conf,thandle=tab_newInput,key='newInput3d')
   ! Manage vasint
   call aot_table_open(L=conf,parent=tab_newInput,thandle=tab_magnitude,key='vacuumPot')
   call aot_get_val(L=conf,errCode=iErr,thandle=tab_magnitude,pos=1,val=auxReal)
   call aot_get_val(L=conf,errCode=iErr,thandle=tab_magnitude,pos=2,val=auxString)
   call en%read(mag=auxReal,units=trim(auxString))
   call en%to_std()
   this%vasint=en%getValue()
   call aot_table_close(L=conf,thandle=tab_magnitude)
   ! Manage rumpling list
   call aot_table_open(L=conf,parent=tab_newInput,thandle=tab_rumpl,key='rumplingList')
   call aot_get_val(L=conf,errCode=iErr,thandle=tab_rumpl,key='units',val=auxString)
   this%nRumpling=aot_table_length(L=conf,thandle=tab_rumpl)-1
   allocate(this%rumpling(this%nRumpling))
   do i=1,this%nRumpling
      call aot_get_val(L=conf,errCode=iErr,thandle=tab_rumpl,pos=i,val=auxReal)
      call dist%read(mag=auxReal,units=trim(auxString))
      call dist%to_std()
      this%rumpling(i)=dist%getValue()
   enddo
   call aot_table_close(L=conf,thandle=tab_rumpl)
   ! Manage Zgrid
   call aot_table_open(L=conf,parent=tab_newInput,thandle=tab_gridInfo,key='gridInfo')
   call aot_get_val(L=conf,errCode=iErr,thandle=tab_gridInfo,key='nPoints',val=this%nzgrid)
   allocate(this%zgrid(this%nzgrid))
   call aot_get_val(L=conf,errCode=iErr,thandle=tab_gridInfo,key='units',val=auxString)
   call aot_get_val(L=conf,errCode=iErr,thandle=tab_gridInfo,key='kind',val=auxString2)
   select case( trim(auxString2) )
   case('Manual','manual','MANUAL')
      call aot_table_open(L=conf,parent=tab_gridInfo,thandle=tab_grid,key='points')
      select case( this%nzgrid==aot_table_length(L=conf,thandle=tab_grid) )
      case(.true.)
         ! OK, continue
      case(.false.)
         write(0,*) 'NEWINPUT3D ERR: nPoints option is different than given manual grid length'
         call exit(1)
      end select
      do i=1,this%nzgrid
         call aot_get_val(L=conf,errCode=iErr,thandle=tab_grid,pos=i,val=auxReal)
         call dist%read(mag=auxReal,units=trim(auxString))
         call dist%to_std()
         this%zgrid(i)=dist%getValue()
      enddo
      call aot_table_close(L=conf,thandle=tab_grid)
      call order_vect(this%zgrid)
   case default
      write(0,*) 'INITIALIZE INPUTGEN ERR: wrong kind of grid. The chosen one is: '//trim(auxString2)
      write(0,*) 'Implemented ones: Manual'
      call exit(1)
   end select
   call aot_table_close(L=conf,thandle=tab_gridInfo)
   ! Check that the manual grid contains rumplings
   allocate(this%zeropos(this%nRumpling))
   this%zeropos(:)=0
	do i=1,this%nzgrid
		do k=1, this%nrumpling
		   select case( this%zgrid(i)==this%rumpling(k) )
		   case(.true.)
				this%zeropos(k)=i
		   case(.false.)
		      ! cycle
		   end select
		enddo
	enddo
	exist_zero=.true.
	do k=1, this%nrumpling
		select case( this%zeropos(k)==0 )
		case(.true.)
			exist_zero=.false.
		case(.false.)
		   ! cycle
		end select
	end do
	! Close the program if there was not any 0.D0 point defined inside the grid
	select case( exist_zero )
	case(.true.)
	   ! continue
	case(.false.)
		write(0,*) "READ_NEWINPUT ERR: Manual grid does not have a rumpling point in the grid"
		call exit(1)
	end select
	! Manage Sitio files
   call aot_table_open(L=conf,parent=tab_newInput,thandle=tab_files,key='sitioFiles')
   this%nFiles=aot_table_length(L=conf, thandle=tab_files)
   allocate(this%inputSitioList(this%nFiles))
   allocate(this%outputSitioList(this%nFiles))
   allocate(this%protoSitios(this%nFiles))
   do i=1,this%nFiles
      call aot_table_open(L=conf, parent=tab_files, thandle=tab_subFile, pos=i)
      call aot_get_val(L=conf,errCode=iErr,thandle=tab_subFile,key='inp',val=this%inputSitioList(i))
      call aot_get_val(L=conf,errCode=iErr,thandle=tab_subFile,key='out',val=this%outputSitioList(i))
      call this%protoSitios(i)%initializeRaw( fileName=trim(this%inputSitioList(i)) )
      call this%printSitio(symmRaw=this%protoSitios(i), fileName=trim(this%outputSitioList(i)))
      call aot_table_close(L=conf,thandle=tab_subFile)
   enddo
   call aot_table_close(L=conf, thandle=tab_files)
   ! Manage pairpot files
   call aot_table_open(L=conf, parent=tab_newInput, thandle=tab_files, key='pairpotFiles')
   select case( this%nRumpling == aot_table_length(L=conf, thandle=tab_files) )
   case(.true.)
      ! continue
   case(.false.)
      write(0,*) 'INITIALIZE NEWINPUT3D ERR: bad number of pairpot files and rumplings'
      write(0,*) 'Number of rumplings defined: ',this%nRumpling
      write(0,*) 'Number of pairpot defined: ', aot_table_length(L=conf, thandle=tab_files)
      call exit(1)
   end select
   allocate(this%inputPairpotList(this%nRumpling))
   allocate(this%outputPairpotList(this%nRumpling))
   allocate(this%protoPairpots(this%nRumpling))
   do i=1,this%nRumpling
      call aot_table_open(L=conf, parent=tab_files, thandle=tab_subFile, pos=i)
      call aot_get_val(L=conf,errCode=iErr,thandle=tab_subFile,key='inp',val=this%inputPairpotList(i))
      call aot_get_val(L=conf,errCode=iErr,thandle=tab_subFile,key='out',val=this%outputPairpotList(i))
      call this%protoPairpots(i)%initializeRaw( fileName=trim(this%inputPairpotList(i)) )
      call this%printPairpot(symmRaw=this%protoPairpots(i),rumplingID=i,fileName=trim(this%outputPairpotList(i)))
      call aot_table_close(L=conf,thandle=tab_subFile)
   enddo
   call aot_table_close(L=conf, thandle=tab_files)
   ! Initialize symmPoints

   ! Close new input table
   call aot_table_close(L=conf,thandle=tab_newInput)
   ! Close Lua conf. file
   call close_config(L=conf)
   call this%printStatus()
	return
end subroutine initialize_NewInput3d
!##########################################################################
!# SUBROUTINE: printStatus_NewInput3d
!##########################################################################
!> @brief
!! This subroutine prints atributes of class NewInput3d
!> @author A.S. Muzas
!> @date May/2016
!> @version 1.0
!--------------------------------------------------------------------------
subroutine printStatus_NewInput3d(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(NewInput3d),intent(in):: this
   ! Local variables
   character(:),parameter:: header='NEWINPUT3D status: '
   integer(kind=4):: i ! counters
   ! Hey Oh, let's go
   print *, header//'Vasint is (au): ',this%vasint
   print *, header//'Number of rumplings: ',this%nRumpling
   do i=1,this%nRumpling
      print *, header//'Rumpling value (bohr): ',this%rumpling(i),' Zero pos: ',this%zeroPos(i)
   enddo
   do i=1,this%nRumpling
      print *, header//'Pairpot: input: '//trim(this%inputPairpotList(i))//' | output: '//trim(this%outputPairpotList(i))
   enddo
   do i=1,this%nFiles
      print *, header//'Sitio: input: '//trim(this%inputSitioList(i))//' | output: '//trim(this%outputSitioList(i))
   enddo
   return
end subroutine printStatus_NewInput3d
!########################################################################
!# SUBROUTINE: GEN_INPUT_PAIRPOT ########################################
!########################################################################
!> @brief
!! This routine creates an standard Pairpot input for future calculations.
!
!> @details
!! - The purpose is to generate a set of CRP pairpots inputs from "raw" ones.
!
!> @param[in] symmRaw - Symmpoint variable whose information will be used in the input generation
!> @param[in] interpol - Newinput variable that contains information about the grid and stuff
!> @param[in] rump_control - Integer(>=0) to select which rumpling correction we want to use 
!> @param[in] fileName - Name of the output file
!
!> @warning
!! - The Symmpoint variable should have been initialized before and
!!   interpolated in Z variable
!! - The Newinput variable should have been initialized before
!! - The output file is written in a.u.
!
!> @author A.M. Muzas - alberto.muzas@uam.es
!> @date 03/Feb/2014
!> @version 1.0
!
!> @see symmRaw, newinput
!--------------------------------------------------------------------
subroutine printPairpot_NewInput3d(this,symmRaw,rumplingID,fileName)
   ! Initial declarations.
   implicit none
   ! I/O Variables ------------------------------
   class(NewInput3d),intent(in):: this
   class(SymmPoint),intent(inout):: symmRaw
	integer,intent(in):: rumplingID ! 0 if there is not rumpling
   character(len=*),intent(in):: fileName
   ! Local variables ----------------------------
	character(len=:),parameter:: routineName = "PRINTPAIRPOT_NEWINPUT3D: "
	integer:: i ! Counter
	real(kind=8):: dz1, dz2
   ! HEY HO! LET'S GO! --------------------------
#ifdef DEBUG
	call debug_write(routineName,"Alias: ",symmRaw%alias)
	call debug_write(routineName,"From file: ",symmRaw%fiNename)
	call debug_write(routineName,"Output file: ",fileName)
	call verbose_write(routineName,"Using rumpling correction: ",rumplingID)
#endif
	! Case we want rumpling correction
	select case( rumplingID > 0 )
	case(.true.)
	   ! continue
	case(.false.)
		write(0,*) 'PRINTPAIRPOT_NEWINPUT3D ERR: Wrong "rumplingID" value'
		call exit(1)
	end select
	! Manage output ------------------------------------
	open(unit=11,file=fileName,status="replace")
	write(11,*) "# This input file was generated by printPairpot_NewInput3d"
	write(11,*) "# Do not modify anything. Everything is in a.u."
	write(11,*) this%vasint, '     VASINT'
	dz1=symmRaw%interz%getderiv(this%zgrid(this%zeropos(rumplingID)))
	dz2=0.d0
	write(11,*) dz1,'     DZ1'
	write(11,*) dz2,'     DZ2'
	write(11,*) rumplingID, this%rumpling(rumplingID), '   ID, RUMPLING'
	write(11,*) this%nzgrid-this%zeropos(rumplingID)+1
	do i=this%zeropos(rumplingID),this%nzgrid
		! We can have problems if the grid goes higher than the values defined for symmRaw
		select case( this%zgrid(i) <= symmRaw%z(symmRaw%n) )
		case(.true.)
			write(11,*) this%zgrid(i), symmRaw%interz%getvalue(this%zgrid(i))
	   case(.false.)
			write(11,*) this%zgrid(i),this%vasint ! Remaining values setted to vasint
		end select
	end do
	close(unit=11)
#ifdef DEBUG
	call verbose_write(routineName,"Pairpot input created successfully")
#endif
	return
end subroutine printPairpot_NewInput3d
!########################################################################
!# SUBROUTINE: GEN_INPUT_SITIO ##########################################
!########################################################################
!> @brief
!! This routine creates an standard Sitio input for future calculations.
!
!> @details
!! - The purpose is to generate a set of CRP Sitio inputs from "raw" ones.
!
!> @param[in] interpol - Newinput variable that contains information about the grid and stuff
!> @param[in] symmRaw - Symmpoint variable whose information will be used in the input generation
!> @param[in] fileName - Name of the output file
!
!> @warning
!! - The Symmpoint variable should have been initialized before and
!!   interpolated in Z variable
!! - The Newinput variable should have been initialized before
!! - The output file is written in a.u.
!
!> @author A.M. Muzas - alberto.muzas@uam.es
!> @date 03/Feb/2014
!> @version 1.0
!
!> @see symmRaw, newinput
!--------------------------------------------------------------------
subroutine printSitio_NewInput3d(this,symmRaw,fileName)
	implicit none
	! I/O variables --------------------
	class(NewInput3d),intent(in):: this
	class(SymmPoint),intent(inout):: symmRaw
	character(len=*),intent(in):: fileName
	! Local variables -------------------
	character(len=:),parameter:: routineName = "GEN_INPUT_SITIO: "
	integer(kind=4):: i ! Counter
	real(kind=8):: dz1,dz2
	! FIRE IN THE HOLE ! ----------------
#ifdef DEBUG
	call debug_write(routineName,"Alias: ",symmRaw%alias)
	call debug_write(routineName,"From file: ",symmRaw%fileName)
	call debug_write(routineName,"Output file: ",fileName)
#endif
	! Prepare input file
	open(11,file=fileName,status="replace")
	write(11,*) "# Input file generated by GEN_INPUT_SITIO"
	write(11,*) "# Do not modify anything. Everything is in a.u."
	write(11,*) symmRaw%x, symmRaw%y, ' <----(X,Y) location in a.u.'
	write(11,*) this%nzgrid, '   NZ'
	dz1=symmRaw%interz%getderiv(this%zgrid(1))
	dz2=this%dfin
	write(11,*) dz1, '    DZ1'
	write(11,*) dz2, '    DZ2'
	do i=1, this%nzgrid
		! We can have problems if the grid goes higher than the values defined for sitio
		select case( this%zgrid(i) <= symmRaw%z(symmRaw%n) )
		case(.true.)
			write(11,*) this%zgrid(i),symmRaw%interz%getvalue(this%zgrid(i))
		case(.false.)
			write(11,*) this%zgrid(i),this%vasint ! Remaining values setted to vasint
		end select
	end do
	close(unit=11)
#ifdef DEBUG
	call verbose_write(routineName,"Sitio input created successfully")
#endif
end subroutine printSitio_NewInput3d

end module NEWINPUT3D_MOD
