!######################################################################
! MODULE: INPUTGEN_MOD
!
!> @brief
!! This module provides tools to generate "standard" inputs from a set 
!! of "raw" inputs
!
!######################################################################
module INPUTGEN_MOD
! Initial declarations
#ifdef DEBUG
   use DEBUG_MOD
#endif
use AOTUS_MODULE, only: flu_State, open_config_file, close_config, aot_get_val
use AOT_TABLE_MODULE, only: aot_table_open, aot_table_close, aot_table_length, aot_table_get_val
use UNITS_MOD
use MATHS_MOD
use CRP3D_MOD
implicit none
!/////////////////////////////////////////////////////////////////////
! TYPE: NewInput
!
!> @brief
!! Stores all data needed to generate a new set of standard input files
!
!> @param vasint - Potential when the atom is far from the surface (in the vacuum)
!> @param dfin - Arbitrary boundary: 1st derivative of all potentials at zgrid(nzgrid)
!> @param nrumpling - Number of different rumplings defined 
!> @param rumpling - Storage of rumplings
!> @param zeropos - Position of ith rumpling in the grid
!> @param nzgrid - Number of points in zgrid 
!> @param zgrid - Grid proposed to generate input files
!
!> @author A.S. Muzas
!> @date 30/Jan/2014
!> @version 1.0
!---------------------------------------------------------------------
type NewInput
   real(kind=8),private:: vasint
   real(kind=8),private:: dfin
   integer(kind=4),private:: nRumpling
   real(kind=8),dimension(:),allocatable,private:: rumpling
   integer,dimension(:),allocatable,private:: zeropos
   integer(kind=4),private:: nzgrid
   real(kind=8),dimension(:),allocatable,private:: zgrid
   contains
      procedure,public:: initialize => initialize_NewInput
      procedure,public:: printStatus => printStatus_NewInput
!      procedure,public:: GEN_PAIRPOT => GEN_INPUT_PAIRPOT
!      procedure,public:: GEN_SITIO => GEN_INPUT_SITIO
end type NewInput
!/////////////////////////////////////////////////////////////////////
contains
!##################################################################
! SUBROUTINE: READ_NEWINPUT
!
!> @brief
!! Read information from file "filename" to create an standard input 
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
subroutine initialize_NewInput(this,filename)
   implicit none
   ! I/O Variables -----------------
   class(NewInput),intent(out):: this
   character(len=*),intent(in):: fileName
   ! Local variables
   type(flu_State):: conf ! Lua file manager
   integer(kind=4):: iErr ! code for Lua errors
   integer(kind=4):: tab_newInput,tab_magnitude,tab_rumpl,tab_gridInfo,tab_grid ! Lua table handlers
   real(kind=8):: auxReal
   character(len=1012):: auxString,auxString2
   integer(kind=4):: i,k ! counters
	type(Length) :: dist
   type(Energy) :: en
	logical :: exist_zero
   integer :: nbefore, nafter !
   real*8, dimension(:), allocatable :: aux1, aux2
	character(len=15), parameter :: routinename = "READ_NEWINPUT: "
   ! Run section -------------------------------------- -----------------
   this%dfin=0.d0
   ! Open Lua conf. file
   call open_config_file(L=conf,errCode=iErr,fileName=trim(fileName))
   ! Open new Input table
   call aot_table_open(L=conf,thandle=tab_newInput,key='newInput')
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
   select case( trim(auxString) )
   case('Manual','manual','MANUAL')
      call aot_table_open(L=conf,parent=tab_gridInfo,thandle=tab_grid,key='points')
      do i=1,this%nzgrid
         call aot_get_val(L=conf,errCode=iErr,thandle=tab_grid,pos=i,val=auxReal)
         call dist%read(mag=auxReal,units=trim(auxString))
         call dist%to_std()
         this%zgrid(i)=dist%getValue()
      enddo
      call aot_table_close(L=conf,thandle=tab_grid)
   end select
   call aot_table_close(L=conf,thandle=tab_gridInfo)



   ! Close new input table
   call aot_table_close(L=conf,thandle=tab_newInput)



   ! Close Lua conf. file
   call close_config(L=conf)
   call this%printStatus()


!   ! Manual grid input (only option available for the moment)
!	if (control.eq."MANU") then
!      allocate(interpol%zgrid(1:interpol%nzgrid))
!      do i=1, interpol%nzgrid
!			read(11,*) aux_r
!         call long%READ(aux_r,units)
!         call long%TO_STD()
!         interpol%zgrid(i) = long%getvalue()
!      end do
!		! Order the grid from low to high values
!		call ORDER_VECT(interpol%zgrid)
!		! Check that the zgrid contains the rumplings defined
!		allocate(interpol%zeropos(1:interpol%nrumpling))
!		do i=1, interpol%nrumpling
!			interpol%zeropos(i)=0 ! Initialize default values for zeropos
!		end do
!		do i=1,interpol%nzgrid
!			do k=1, interpol%nrumpling
!				if (interpol%zgrid(i).eq.interpol%rumpling(k)) then
!					interpol%zeropos(k)=i
!				end if
!			end do
!		end do
!		exist_zero=.true.
!		do k=1, interpol%nrumpling
!			if (interpol%zeropos(k).eq.0) then
!				exist_zero=.false.
!			end if
!		end do
!		! Close the program if there was not any 0.D0 point defined inside the grid
!		if(.not.exist_zero) then
!			write(0,*) "READ_NEWINPUT ERR: Manual grid does not have a rumpling point in the grid"
!			stop
!		end if
!	! Should define two grids: one from the lowest value to 0 and the other from 0
!	! to the highest value in the grid.
!        else
!                write(0,*) 'READ_NEWINPUT ERR: Variable "control" is not correct'
!                write(0,*) 'READ_NEWINPUT ERR: It can only be "MANU"'
!        end if
!        close(11)
!	return
end subroutine initialize_NewInput
!##########################################################################
!# SUBROUTINE: printStatus_NewInput
!##########################################################################
!> @brief
!! This subroutine prints atributes of class NewInput
!> @author A.S. Muzas
!> @date May/2016
!> @version 1.0
!--------------------------------------------------------------------------
subroutine printStatus_NewInput(this)
   ! Initial declarations
   implicit none
   ! I/O variables
   class(NewInput),intent(in):: this
   ! Local variables
   character(:),parameter:: header='NEWINPUT status: '
   integer(kind=4):: i ! counters
   ! Hey Oh, let's go
   print *, header//'Vasint is (au): ',this%vasint
   print *, header//'Number of rumplings: ',this%nRumpling
   do i=1,this%nRumpling
      print *, header//'Rumpling value (bohr): ',this%rumpling(i)
   enddo
   return
end subroutine printStatus_NewInput
!########################################################################
!# SUBROUTINE: GEN_INPUT_PAIRPOT ########################################
!########################################################################
!> @brief
!! This routine creates an standard Pairpot input for future calculations.
!
!> @details
!! - The purpose is to generate a set of CRP pairpots inputs from "raw" ones.
!
!> @param[in] symmraw - Symmpoint variable whose information will be used in the input generation
!> @param[in] interpol - Newinput variable that contains information about the grid and stuff
!> @param[in] rump_control - Integer(>=0) to select which rumpling correction we want to use 
!> @param[in] filename - Name of the output file 
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
!> @see symmraw, newinput
!--------------------------------------------------------------------
!subroutine GEN_INPUT_PAIRPOT(interpol,symmraw,rump_control,filename)
!   implicit none
!   ! I/O Variables ------------------------------
!   class(Newinput),intent(in)  :: interpol
!   class(Symmpoint),intent(inout) :: symmraw
!	integer,intent(in) :: rump_control ! 0 if there is not rumpling
!   character(len=*), intent(in) :: filename
!	character(len=20), parameter :: routinename = "GEN_INPUT_PAIRPOT: "
!        ! Local variables ----------------------------
!	integer :: i ! Counter
!	real*8 :: dz1, dz2
!        ! HEY HO! LET'S GO! --------------------------
!#ifdef DEBUG
!	call DEBUG_WRITE(routinename,"Alias: ", symmraw%alias)
!	call DEBUG_WRITE(routinename,"From file: ", symmraw%filename)
!	call DEBUG_WRITE(routinename,"Output file: ", filename)
!	call VERBOSE_WRITE(routinename,"Using rumpling correction: ", rump_control)
!#endif
!	! Case we want rumpling correction
!	if (rump_control.le.0) then
!		write(0,*) 'MAKE_INPUT_PAIRPOT ERR: Wrong "rump_control" value'
!		stop
!	end if
!	! Manage output ------------------------------------
!	open(11,file=filename,status="replace")
!	write(11,*) "# This input file was generated by GEN_INPUT_PAIRPOT"
!	write(11,*) "# Do not modify anything. Everything is in a.u."
!	write(11,*) interpol%vasint, '     VASINT'
!	dz1=symmraw%interz%getderiv(interpol%zgrid(interpol%zeropos(rump_control)))
!	dz2=interpol%dfin
!	write(11,*) dz1,'     DZ1'
!	write(11,*) dz2,'     DZ2'
!	write(11,*) rump_control, interpol%rumpling(rump_control), '   ID, RUMPLING'
!	write(11,*) interpol%nzgrid-interpol%zeropos(rump_control)+1
!	do i=interpol%zeropos(rump_control), interpol%nzgrid
!		! We can have problems if the grid goes higher than the values defined for symmraw
!		if (interpol%zgrid(i).le.symmraw%z(symmraw%n)) then
!			write(11,*) interpol%zgrid(i), symmraw%interz%getvalue(interpol%zgrid(i))
!		else if (interpol%zgrid(i).gt.symmraw%z(symmraw%n)) then
!			write(11,*) interpol%zgrid(i), interpol%vasint ! Remaining values setted to vasint
!		else
!			write(0,*) "GEN_INPUT_PAIRPOT ERR: Take a look at what you are doing!"
!			stop
!
!		end if
!	end do
!	close(11)
!#ifdef DEBUG
!	call VERBOSE_WRITE(routinename,"Pairpot input created successfully")
!#endif
!	return
!end subroutine GEN_INPUT_PAIRPOT
!!########################################################################
!!# SUBROUTINE: GEN_INPUT_SITIO ##########################################
!!########################################################################
!!> @brief
!!! This routine creates an standard Sitio input for future calculations.
!!
!!> @details
!!! - The purpose is to generate a set of CRP Sitio inputs from "raw" ones.
!!
!!> @param[in] interpol - Newinput variable that contains information about the grid and stuff
!!> @param[in] symmraw - Symmpoint variable whose information will be used in the input generation
!!> @param[in] filename - Name of the output file
!!
!!> @warning
!!! - The Symmpoint variable should have been initialized before and
!!!   interpolated in Z variable
!!! - The Newinput variable should have been initialized before
!!! - The output file is written in a.u.
!!
!!> @author A.M. Muzas - alberto.muzas@uam.es
!!> @date 03/Feb/2014
!!> @version 1.0
!!
!!> @see symmraw, newinput
!!--------------------------------------------------------------------
!subroutine GEN_INPUT_SITIO(interpol,symmraw,filename)
!	implicit none
!	! I/O variables --------------------
!	class(Newinput), intent(in), target :: interpol
!	class(Symmpoint),intent(inout), target :: symmraw
!	character(len=*), intent(in) :: filename
!	character(len=17), parameter :: routinename = "GEN_INPUT_SITIO: "
!	! Local variables -------------------
!	integer :: i ! Counter
!	real*8 :: dz1, dz2
!	! FIRE IN THE HOLE ! ----------------
!#ifdef DEBUG
!	call DEBUG_WRITE(routinename,"Alias: ", symmraw%alias)
!	call DEBUG_WRITE(routinename,"From file: ", symmraw%filename)
!	call DEBUG_WRITE(routinename,"Output file: ", filename)
!#endif
!	! Prepare input file
!	open(11,file=filename,status="replace")
!	write(11,*) "# Input file generated by GEN_INPUT_SITIO"
!	write(11,*) "# Do not modify anything. Everything is in a.u."
!	write(11,*) symmraw%x, symmraw%y, ' <----(X,Y) location in a.u.'
!	write(11,*) interpol%nzgrid, '   NZ'
!	dz1=symmraw%interz%getderiv(interpol%zgrid(1))
!	dz2=interpol%dfin
!	write(11,*) dz1, '    DZ1'
!	write(11,*) dz2, '    DZ2'
!	do i=1, interpol%nzgrid
!		! We can have problems if the grid goes higher than the values defined for sitio
!		if (interpol%zgrid(i).le.symmraw%z(symmraw%n)) then
!			write(11,*) interpol%zgrid(i), symmraw%interz%getvalue(interpol%zgrid(i))
!		else if (interpol%zgrid(i).gt.symmraw%z(symmraw%n)) then
!			write(11,*) interpol%zgrid(i), interpol%vasint ! Remaining values setted to vasint
!		else
!			write(0,*) "MAKE_INPUT_SITIO ERR: Take a look at what you are doing!"
!			stop
!		end if
!	end do
!	close(11)
!#ifdef DEBUG
!	call VERBOSE_WRITE(routinename,"Sitio input created successfully")
!#endif
!end subroutine GEN_INPUT_SITIO
end module INPUTGEN_MOD
