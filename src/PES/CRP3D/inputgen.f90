!######################################################################
! MODULE: INPUTGEN_MOD
!
!> @brief
!! This module provides tools to generate "standard" inputs from a set 
!! of "raw" inputs
!
!######################################################################
MODULE INPUTGEN_MOD
! Initial declarations
IMPLICIT NONE
!/////////////////////////////////////////////////////////////////////
! TYPE: Newinput
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
TYPE Newinput
   PRIVATE
   REAL(KIND=8) :: vasint 
   REAL(KIND=8) :: dfin 
   INTEGER :: nrumpling
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: rumpling 
   INTEGER,DIMENSION(:),ALLOCATABLE :: zeropos
   INTEGER :: nzgrid 
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: zgrid 
   CONTAINS
      PROCEDURE,PUBLIC :: READ => READ_NEWINPUT
      PROCEDURE,PUBLIC :: GEN_PAIRPOT => GEN_INPUT_PAIRPOT
      PROCEDURE,PUBLIC :: GEN_SITIO => GEN_INPUT_SITIO
END TYPE Newinput
!/////////////////////////////////////////////////////////////////////
CONTAINS
!##################################################################
! SUBROUTINE: READ_NEWINPUT
!
!> @brief
!! Read information from file "filename" to create an standard input 
!! generation job
!
!> @param interpol - Newinput type variable
!> @param filename - Input file to be read
!
!> @author A.S. Muzas
!> @date 30/Jan/2014
!> @version 1.0
!
!> @warning
!! - In this routine Newinput%dfin is set to zero, which is a good value for potentials far from the surface
!! - Input file should have the following structure:
!!    -# line 1: @b real(kind=8), @b character(len=10) ; energy at long distance from the surface, units
!!    -# line 2: @b integer(kind=4); (N) number of different rumplings defined (usually, number of different atoms in the surface)
!!    -# lines 3~3+N: @b real(kind=8), @b character(len=10); list of rumplings and units
!!    -# line 3+N+1: @b integer(kind=4), @b character(len=4), @b real(kind=8), @b real(kind=8), @b character(len=10); number of
!!                   points in the grid (M), keyword: MANU (for manual grid, only option available now), first and last values 
!!                   in the grid, units
!!    -# lines 3N+2~3N+2+M: real(kind=8); grid value. There is no need to have an ordered grid, this subroutine will order it.
!-----------------------------------------------------------------
SUBROUTINE READ_NEWINPUT(interpol,filename)
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   USE UNITS_MOD
   USE MATHS_MOD
   IMPLICIT NONE
   ! I/O Variables -----------------
   CLASS(Newinput),INTENT(OUT) :: interpol
   CHARACTER(LEN=*),INTENT(IN) :: filename
   ! Local variables
   INTEGER :: i,k ! counter
	TYPE(Length) :: long, long2
   TYPE(Energy) :: vasint
   REAL(KIND=8) :: aux_r, aux_r2
	CHARACTER(LEN=10) :: units
   CHARACTER(LEN=4) :: control
	LOGICAL :: exist_zero
   INTEGER :: nbefore, nafter !
   REAL*8, DIMENSION(:), ALLOCATABLE :: aux1, aux2
	CHARACTER(LEN=15), PARAMETER :: routinename = "READ_NEWINPUT: "
   ! Run section -------------------------------------- -----------------
	CALL VERBOSE_WRITE(routinename,"New interpolation scheme defined")
   OPEN(11,FILE=filename,STATUS="old")
	!----
   READ(11,*) aux_r, units
   CALL vasint%READ(aux_r,units)
   CALL vasint%TO_STD()
   interpol%vasint = vasint%getvalue()
	!------
   READ(11,*) interpol%nrumpling
	!-----
   ALLOCATE(interpol%rumpling(1:interpol%nrumpling))
	interpol%dfin = 0.D0
	DO i=1,interpol%nrumpling
      READ(11,*) aux_r, units
      CALL long%READ(aux_r,units)
      CALL long%TO_STD()
		interpol%rumpling(i) = long%getvalue()
	END DO
   READ(11,*) interpol%nzgrid,control,aux_r,aux_r2,units
   CALL long%READ(aux_r,units)
   CALL long2%READ(aux_r2,units)
   CALL long%TO_STD()
   CALL long2%TO_STD()
	interpol%first = long%getvalue()
	interpol%last = long2%getvalue()
   ! Manual grid input (only option available for the moment)
	IF (control.EQ."MANU") THEN
      ALLOCATE(interpol%zgrid(1:interpol%nzgrid))
      DO i=1, interpol%nzgrid
			READ(11,*) aux_r
         CALL long%READ(aux_r,units)
         CALL long%TO_STD()
         interpol%zgrid(i) = long%getvalue()
      END DO
		! Order the grid from low to high values
		CALL ORDER_VECT(interpol%zgrid)
		! Check that the zgrid contains the rumplings defined
		ALLOCATE(interpol%zeropos(1:interpol%nrumpling))
		DO i=1, interpol%nrumpling
			interpol%zeropos(i)=0 ! Initialize default values for zeropos
		END DO
		DO i=1,interpol%nzgrid
			DO k=1, interpol%nrumpling
				IF (interpol%zgrid(i).EQ.interpol%rumpling(k)) THEN
					interpol%zeropos(k)=i
				END IF
			END DO
		END DO
		exist_zero=.TRUE.
		DO k=1, interpol%nrumpling
			IF (interpol%zeropos(k).EQ.0) THEN
				exist_zero=.FALSE.
			END IF
		END DO
		! Close the program if there was not any 0.D0 point defined inside the grid
		IF(.NOT.exist_zero) THEN
			WRITE(0,*) "READ_NEWINPUT ERR: Manual grid does not have a rumpling point in the grid"
			STOP
		END IF
	! Should define two grids: one from the lowest value to 0 and the other from 0
	! to the highest value in the grid. 
        ELSE
                WRITE(0,*) 'READ_NEWINPUT ERR: Variable "control" is not correct'
                WRITE(0,*) 'READ_NEWINPUT ERR: It can only be "MANU"'
        END IF
        CLOSE(11)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,"New input job read from file",filename)
#endif
	RETURN
END SUBROUTINE READ_NEWINPUT
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
SUBROUTINE GEN_INPUT_PAIRPOT(interpol,symmraw,rump_control,filename)
   USE CRP_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   IMPLICIT NONE
   ! I/O Variables ------------------------------
   CLASS(Newinput),INTENT(IN)  :: interpol
   CLASS(Symmpoint),INTENT(INOUT) :: symmraw
	INTEGER,INTENT(IN) :: rump_control ! 0 if there is not rumpling
   CHARACTER(LEN=*), INTENT(IN) :: filename
	CHARACTER(LEN=20), PARAMETER :: routinename = "GEN_INPUT_PAIRPOT: "
        ! Local variables ----------------------------
	INTEGER :: i ! Counter
	REAL*8 :: dz1, dz2
        ! HEY HO! LET'S GO! --------------------------
#ifdef DEBUG
	CALL DEBUG_WRITE(routinename,"Alias: ", symmraw%alias)
	CALL DEBUG_WRITE(routinename,"From file: ", symmraw%filename)
	CALL DEBUG_WRITE(routinename,"Output file: ", filename)
	CALL VERBOSE_WRITE(routinename,"Using rumpling correction: ", rump_control)
#endif
	! Case we want rumpling correction
	IF (rump_control.LE.0) THEN
		WRITE(0,*) 'MAKE_INPUT_PAIRPOT ERR: Wrong "rump_control" value'
		STOP
	END IF
	! Manage output ------------------------------------
	OPEN(11,FILE=filename,STATUS="replace")
	WRITE(11,*) "# This input file was generated by GEN_INPUT_PAIRPOT"
	WRITE(11,*) "# Do not modify anything. Everything is in a.u."
	WRITE(11,*) interpol%vasint, '     VASINT'
	dz1=symmraw%interz%getderiv(interpol%zgrid(interpol%zeropos(rump_control)))
	dz2=interpol%dfin
	WRITE(11,*) dz1,'     DZ1'
	WRITE(11,*) dz2,'     DZ2'
	WRITE(11,*) rump_control, interpol%rumpling(rump_control), '   ID, RUMPLING'
	WRITE(11,*) interpol%nzgrid-interpol%zeropos(rump_control)+1
	DO i=interpol%zeropos(rump_control), interpol%nzgrid
		! We can have problems if the grid goes higher than the values defined for symmraw
		IF (interpol%zgrid(i).LE.symmraw%z(symmraw%n)) THEN
			WRITE(11,*) interpol%zgrid(i), symmraw%interz%getvalue(interpol%zgrid(i))
		ELSE IF (interpol%zgrid(i).GT.symmraw%z(symmraw%n)) THEN
			WRITE(11,*) interpol%zgrid(i), interpol%vasint ! Remaining values setted to vasint
		ELSE
			WRITE(0,*) "GEN_INPUT_PAIRPOT ERR: Take a look at what you are doing!"
			STOP
			
		END IF
	END DO
	CLOSE(11)
#ifdef DEBUG
	CALL VERBOSE_WRITE(routinename,"Pairpot input created successfully")
#endif
	RETURN
END SUBROUTINE GEN_INPUT_PAIRPOT
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
!> @param[in] symmraw - Symmpoint variable whose information will be used in the input generation
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
SUBROUTINE GEN_INPUT_SITIO(interpol,symmraw,filename)
#ifdef DEBUG
   USE DEBUG_MOD
#endif
   USE CRP_MOD
	IMPLICIT NONE
	! I/O variables --------------------
	CLASS(Newinput), INTENT(IN), TARGET :: interpol
	CLASS(Symmpoint),INTENT(INOUT), TARGET :: symmraw
	CHARACTER(LEN=*), INTENT(IN) :: filename
	CHARACTER(LEN=17), PARAMETER :: routinename = "GEN_INPUT_SITIO: "
	! Local variables -------------------
	INTEGER :: i ! Counter
	REAL*8 :: dz1, dz2
	! FIRE IN THE HOLE ! ----------------
#ifdef DEBUG
	CALL DEBUG_WRITE(routinename,"Alias: ", symmraw%alias)
	CALL DEBUG_WRITE(routinename,"From file: ", symmraw%filename)
	CALL DEBUG_WRITE(routinename,"Output file: ", filename)
#endif
	! Prepare input file
	OPEN(11,FILE=filename,STATUS="replace")
	WRITE(11,*) "# Input file generated by GEN_INPUT_SITIO"
	WRITE(11,*) "# Do not modify anything. Everything is in a.u."
	WRITE(11,*) symmraw%x, symmraw%y, ' <----(X,Y) location in a.u.'
	WRITE(11,*) interpol%nzgrid, '   NZ'
	dz1=symmraw%interz%getderiv(interpol%zgrid(1))
	dz2=interpol%dfin
	WRITE(11,*) dz1, '    DZ1'
	WRITE(11,*) dz2, '    DZ2'
	DO i=1, interpol%nzgrid
		! We can have problems if the grid goes higher than the values defined for sitio
		IF (interpol%zgrid(i).LE.symmraw%z(symmraw%n)) THEN
			WRITE(11,*) interpol%zgrid(i), symmraw%interz%getvalue(interpol%zgrid(i))
		ELSE IF (interpol%zgrid(i).GT.symmraw%z(symmraw%n)) THEN
			WRITE(11,*) interpol%zgrid(i), interpol%vasint ! Remaining values setted to vasint
		ELSE
			WRITE(0,*) "MAKE_INPUT_SITIO ERR: Take a look at what you are doing!"
			STOP
		END IF
	END DO
	CLOSE(11)
#ifdef DEBUG
	CALL VERBOSE_WRITE(routinename,"Sitio input created successfully")
#endif
END SUBROUTINE GEN_INPUT_SITIO
END MODULE INPUTGEN_MOD
