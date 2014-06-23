!###############################################
! MODULE: INITATOM_MOD
!> @brief
!! This module provides routines and onjects to create
!! initial conditions for an atom or list of atoms
!###############################################
MODULE INITATOM_MOD
USE INICOND_MOD
USE UNITS_MOD
IMPLICIT NONE
TYPE,EXTENDS(Inicond) :: Initatom
   LOGICAL :: control_vel, control_posX, control_posY, control_out, control_seed
   REAL*8 :: impact_x, impact_y
   TYPE(Energy) :: E_norm
   TYPE(Angle) :: vz_angle, vpar_angle
   TYPE(Length) :: init_z ! initial Z value
   TYPE(Mass) :: masss
   CONTAINS
      PROCEDURE,PUBLIC :: READ => READ_INICOND_ATOM
      PROCEDURE,PUBLIC :: GENERATE_TRAJS => GENERATE_TRAJS_ATOMS
END TYPE Initatom
!/////////////////////////////////////////////////////
! TYPE: Atoms
!> @brief
!! Atom subtype dynamics object
!----------------------------------------------------
TYPE,EXTENDS(Dynobject) ::  Atom
   INTEGER :: ireb=0 ! times this atom has changed Pz's direction
   INTEGER :: ixyboun=0 ! times this atoms has changed parallel momentum's direction
END TYPE Atom
!/////////////////////////////////////////////////////
! TYPE: Atom_trajs
!> @brief
!! A list of atom subtype dynamics objects
!----------------------------------------------------
TYPE :: Atom_trajs
   TYPE(Atom),DIMENSION(:),ALLOCATABLE :: atomo
END TYPE Atom_trajs
!////////////////////////////////////////////////////
CONTAINS
!##################################################################################
!# SUBROUTINE: READ_INICOND_ATOM ##################################################
!##################################################################################
!> @brief
!! Reads from an input file enough data to generate initial conditions
!! for a batch of atoms. Adapted to atom/surface systems
!
!> @param[out] inicondat - Initial conditions for a batch of atoms
!> @param[in] filename - Name of the input file
!
!> @warning
!! - Input file syntax:
!!    -# line 1: character(len=30); human friendly alias 
!!    -# line 2: character(len=10); kind of input. "Atoms" is the only available label
!!    -# line 3: integer(kind=4); initial trajectory
!!    -# line 4: integer(kind=4); final trajectory
!!    -# line 5: real(kind=8),character(len=10); mass, units
!!    -# line 6: real(kind=8),character(len=10); initial perpendicular energy, units
!!    -# line 7: real(kind=8),character(len=10); angle of velocity respect to surface plane, units
!!    -# line 8: real(kind=8),character(len=10); angle of velocity respect to surface vector S1, units
!!    -# line 9: real(kind=8),character(len=10); initial Z, units
!!    -# line 10: logical,logical; random seed generation in X and Y?
!!    -# line 11: real(kind=8); initial X parameter in surface coordinates 1<= x >= 0
!!    -# line 12: real(kind=8); initial Y parameter in surface coordinates 1<= x >= 0
!!    -# line 13: logical; Create output file?
!!    -# line 14: character(len=30); name of the output file
!!    -# line 15: logical; read seed from seed.inp?
!
!> @author A.S. Muzas alberto.muzas@uam.es
!> @date 10/Feb/2014
!> @version 1.0
!------------------------------------------------------------------------------
SUBROUTINE READ_INICOND_ATOM(inicondat,filename)
   USE UNITS_MOD
   USE CONSTANTS_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
	IMPLICIT NONE
	! I/O variables
	CLASS(Initatom), INTENT(OUT) :: inicondat
	CHARACTER(LEN=*), INTENT(IN) :: filename
	! Local variables
   REAL(KIND=8) :: aux
   CHARACTER(LEN=10) :: units
	INTEGER :: i ! counters
	INTEGER :: size_seed, clock
	CHARACTER(LEN=23), PARAMETER :: routinename = "DEFINE_INICOND_SCHEME: "
	! YIPEE KI YAY !! -------
	inicondat%input_file = filename
	OPEN(11,FILE=filename,STATUS="old")
	READ(11,*) inicondat%alias
	READ(11,*) inicondat%kind
	READ(11,*) inicondat%nstart
	READ(11,*) inicondat%ntraj
#ifdef DEBUG
	CALL VERBOSE_WRITE(routinename,"Alias: ",inicondat%alias)
	CALL VERBOSE_WRITE(routinename,"Kind: ",inicondat%kind)
	CALL VERBOSE_WRITE(routinename,"Initial traj: ",inicondat%nstart)
	CALL VERBOSE_WRITE(routinename,"Final traj: ",inicondat%ntraj)
#endif
	IF (inicondat%kind.EQ."Atoms") THEN
      READ(11,*) aux,units
      CALL inicondat%masss%READ(aux,units)
		CALL inicondat%masss%TO_STD()
		
      READ(11,*) aux,units
		CALL inicondat%E_norm%READ(aux,units)
		CALL inicondat%E_norm%TO_STD()
		
      READ(11,*) aux,units
		CALL inicondat%vz_angle%READ(aux,units)
		CALL inicondat%vz_angle%TO_STD()
		IF (inicondat%vz_angle%getvalue().GT.90.D0*PI/180.D0) THEN
			WRITE(0,*) "DEFINE_INICOND_SCHEME ERR: incidence angle greater than 90.0 deg (1.57079632679D0 rad)"
			WRITE(0,*) "Incidence angle (rad) : ",inicondat%vz_angle%getvalue()
			CALL EXIT(1)
		ELSE IF (inicondat%vz_angle%getvalue().LT.0.D0) THEN
			WRITE(0,*) "DEFINE_INICOND_SCHEME ERR: incidence angle lower than 0.0 deg."
			WRITE(0,*) "Incidence angle (rad) : ",inicondat%vz_angle%getvalue()
			CALL EXIT(1)
		END IF
		
      READ(11,*) aux,units
		CALL inicondat%vpar_angle%READ(aux,units)
		CALL inicondat%vpar_angle%TO_STD()
		IF (inicondat%vpar_angle%getvalue().GT.90.D0*PI/180.D0) THEN
			WRITE(0,*) "DEFINE_INICOND_SCHEME ERR: parallele angle greater than 90.0 deg (1.57079632679D0 rad)"
			WRITE(0,*) "Incidence angle (rad) : ",inicondat%vpar_angle%getvalue()
			STOP
		ELSE IF (inicondat%vz_angle%getvalue().LT.0.D0) THEN
			WRITE(0,*) "DEFINE_INICOND_SCHEME ERR: parallele angle lower than 0.0 deg."
			WRITE(0,*) "Incidence angle (rad) : ",inicondat%vpar_angle%getvalue()
			STOP
		END IF
		
      READ(11,*) aux,units
		CALL inicondat%init_z%READ(aux,units)
		CALL inicondat%init_z%TO_STD()
#ifdef DEBUG
		CALL VERBOSE_WRITE(routinename,"Mass in au: ",inicondat%masss%getvalue())
		CALL VERBOSE_WRITE(routinename,"Initial normal energy in au:",inicondat%E_norm%getvalue())
		CALL VERBOSE_WRITE(routinename,"Incidence angle respect to surface plane in radians: ",inicondat%vz_angle%getvalue())
		CALL VERBOSE_WRITE(routinename,"Angle between trajectory and S1 vector in radians",inicondat%vpar_angle%getvalue())
		CALL VERBOSE_WRITE(routinename,"Initial Z in au: ",inicondat%init_z%getvalue())
#endif
		READ(11,*) inicondat%control_posX, inicondat%control_posY 
		READ(11,*) inicondat%impact_x
		READ(11,*) inicondat%impact_y
#ifdef DEBUG
		CALL VERBOSE_WRITE(routinename,"Random X impact parameter?: ",inicondat%control_posX)
		CALL VERBOSE_WRITE(routinename,"Random Y impact parameter?: ",inicondat%control_posY)
		CALL VERBOSE_WRITE(routinename,"Impact param X: ",inicondat%impact_x)
		CALL VERBOSE_WRITE(routinename,"Impact param Y: ",inicondat%impact_x)
#endif
		IF(((inicondat%impact_x.GT.1.D0).OR.(inicondat%impact_x.LT.0.D0)).AND.(.NOT.inicondat%control_posX)) THEN
			WRITE(0,*) "DEFINE_INICOND_SCHEME ERR: X impact parameter outside range 0-1"
			CALL EXIT(1)
		ELSE IF(((inicondat%impact_y.GT.1.D0).OR.(inicondat%impact_y.LT.0.D0)).AND.(.NOT.inicondat%control_posy))THEN
			WRITE(0,*) "DEFINE_INICOND_SCHEME ERR: Y impact parameter outside range 0-1"
			CALL EXIT(1)
		END IF
		READ(11,*) inicondat%control_out
		READ(11,*) inicondat%output_file
		READ(11,*) inicondat%control_seed
#ifdef DEBUG
		CALL VERBOSE_WRITE(routinename,"Output files?: ",inicondat%control_out)
		CALL VERBOSE_WRITE(routinename,"Output file name: ",inicondat%output_file)
		CALL VERBOSE_WRITE(routinename,"Seed read from file?: ",inicondat%control_seed)
#endif
		IF (inicondat%control_seed.EQV..TRUE.) THEN
			CALL RANDOM_SEED(SIZE=size_seed)
#ifdef DEBUG
			CALL VERBOSE_WRITE(routinename,"Default size for seed array: ",size_seed)
#endif
			ALLOCATE(inicondat%seed(1:size_seed))
			OPEN(12,FILE="seed.inp",STATUS="old")
			READ(12,*) inicondat%seed
			CLOSE(12)
		ELSE
			CALL RANDOM_SEED(SIZE=size_seed)
			ALLOCATE(inicondat%seed(1:size_seed))
			CALL SYSTEM_CLOCK(COUNT=clock)
			inicondat%seed = clock+ 37*(/ (i - 1, i = 1, size_seed) /)
			CALL RANDOM_SEED(PUT=inicondat%seed)
#ifdef DEBUG
			CALL VERBOSE_WRITE(routinename,"Seed generated from CPU time: ")
			CALL VERBOSE_WRITE(routinename,"CPU time: ", clock)
#endif
			IF (inicondat%control_out.EQV..TRUE.) THEN
				OPEN(12,FILE="seed.inp",STATUS="replace")
				WRITE(12,*) inicondat%seed
				CLOSE(12)
			END IF
		END IF
#ifdef DEBUG
		DO i=1,size_seed
			CALL VERBOSE_WRITE(routinename,inicondat%seed(i))
		END DO
#endif
	ELSE
		WRITE(0,*) "DEFINE_INICOND_SCHEME ERR: Wrong kind of initial conditions"
		WRITE(0,*) "Available kinds: Atoms"
		WRITE(0,*) "You wrote: ", inicondat%kind
		CALL EXIT(1)
	END IF
	CLOSE(11)
	RETURN
END SUBROUTINE READ_INICOND_ATOM
!####################################################################
!# SUBROUTINE: GENERATE_TRAJS_ATOMS #################################
!####################################################################
!> @brief
!! Creates trajectories of atoms. All data stored as atomic units.
!
!> @param[in] inicondat - Initial conditions for atoms to be used
!> @param[in] surf - Surface used
!> @param[out] final_trajs - Trajectories with initial conditions allocated
!
!> @warning
!! - Auxiliar cartesian coordinates.
!! - For atom/surface systems
!! - @b final_trajs should not be allocated before this routine is invoked
!
!> @author A.S. Muzas alberto.muzas@uam.es
!> @date 10/Feb/2014
!> @version 1.0
!--------------------------------------------------------------------
SUBROUTINE GENERATE_TRAJS_ATOMS(inicondat,thispes,final_trajs)
   USE UNITS_MOD
   USE CONSTANTS_MOD
   USE CRP3D_MOD
#ifdef DEBUG
   USE DEBUG_MOD
#endif
	IMPLICIT NONE
	! I/O variables
	CLASS(Initatom),TARGET,INTENT(IN) :: inicondat
	TYPE(Atom_trajs),INTENT(OUT) :: final_trajs
	TYPE(CRP3D),INTENT(IN) :: thispes
	! Local variables
	INTEGER :: i, k,j ! counters
	CHARACTER(LEN=22), PARAMETER :: routinename = "GENERATE_TRAJS_ATOMS: "
	TYPE(Atom_trajs) :: trajs
	REAL*8,DIMENSION(2) :: proj_iws_r
	REAL*8 :: delta, alpha, Enorm, masa
	! YIPPIEE KI YAY !! -------------------
#ifdef DEBUG
	CALL VERBOSE_WRITE(routinename,"New set of trajectories")
	CALL VERBOSE_WRITE(routinename,"Allocating trajs: ", inicondat%ntraj)
#endif
	ALLOCATE(trajs%atomo(1:inicondat%ntraj))
	ALLOCATE(final_trajs%atomo(1:inicondat%ntraj-inicondat%nstart+1))
	IF((inicondat%control_posX).AND.(.NOT.inicondat%control_posY)) THEN
		DO i=1, inicondat%ntraj
			CALL trajs%atomo(i)%INITIALIZE(i,3)
			trajs%atomo(i)%r(3) = inicondat%init_z%getvalue()
			CALL RANDOM_NUMBER(trajs%atomo(i)%r(1))
			trajs%atomo(i)%r(2) = inicondat%impact_y
		END DO
	ELSE IF ((.NOT.inicondat%control_posX).AND.(inicondat%control_posY)) THEN
		DO i=1, inicondat%ntraj
			CALL trajs%atomo(i)%INITIALIZE(i,3)
			trajs%atomo(i)%r(3) = inicondat%init_z%getvalue()
			trajs%atomo(i)%r(1) = inicondat%impact_x
			CALL RANDOM_NUMBER(trajs%atomo(i)%r(2))
		END DO
	ELSE IF ((inicondat%control_posX).AND.(inicondat%control_posY)) THEN
		DO i=1, inicondat%ntraj
			CALL trajs%atomo(i)%INITIALIZE(i,3)
			trajs%atomo(i)%r(3) = inicondat%init_z%getvalue()
			CALL RANDOM_NUMBER(trajs%atomo(i)%r(1:2))
		END DO
	ELSE
		DO i=1, inicondat%ntraj
			CALL trajs%atomo(i)%INITIALIZE(i,3)
			trajs%atomo(i)%r(3) = inicondat%init_z%getvalue()
			trajs%atomo(i)%r(1) = inicondat%impact_x
			trajs%atomo(i)%r(2) = inicondat%impact_y
		END DO
	END IF
	DO i=1,inicondat%ntraj
		! Change to cartesian coordinates (impact parameters are in surface coordinates)
		trajs%atomo(i)%r(1:2) = thispes%surf%surf2cart(trajs%atomo(i)%r(1:2))
		! projectin into IWS cell leads to errors in the dynamics (wrong sampling)
		IF (thispes%surf%units/="au") THEN
         WRITE(0,*) "GENERATE_TRAJS_ATOMS ERR: Surface should be in atomic units!"
         CALL EXIT(1)
      END IF   
		delta = inicondat%vpar_angle%getvalue()
		alpha = inicondat%vz_angle%getvalue()
		masa = inicondat%masss%getvalue()
		Enorm = inicondat%E_norm%getvalue()
		trajs%atomo(i)%p(3) = -DSQRT(2.D0*masa*Enorm)  ! Z momentum (m*v_z), negative (pointing uppon the surface)
		trajs%atomo(i)%p(1) = DCOS(delta)*DSQRT(2.D0*masa*Enorm/(DTAN(alpha)**2.D0))
		trajs%atomo(i)%p(2) = DSIN(delta)*DSQRT(2.D0*masa*Enorm/(DTAN(alpha)**2.D0))
		trajs%atomo(i)%E = Enorm/(DSIN(alpha)**2.D0)
		! Setting initial values
		trajs%atomo(i)%init_r(1) = trajs%atomo(i)%r(1)
		trajs%atomo(i)%init_r(2) = trajs%atomo(i)%r(2)
		trajs%atomo(i)%init_r(3) = trajs%atomo(i)%r(3)
		trajs%atomo(i)%init_p(1) = trajs%atomo(i)%p(1)
		trajs%atomo(i)%init_p(2) = trajs%atomo(i)%p(2)
		trajs%atomo(i)%init_p(3) = trajs%atomo(i)%p(3)
	END DO
#ifdef DEBUG
	CALL VERBOSE_WRITE(routinename,"Using only user-defined trajectories starting from nstart: ")
	CALL VERBOSE_WRITE(routinename,"Size of final_trajs: ",SIZE(final_trajs%atomo))
	CALL VERBOSE_WRITE(routinename,"Size of trajs: ",SIZE(trajs%atomo))
#endif
	DO i=1 , inicondat%ntraj-inicondat%nstart+1
		k = i+inicondat%nstart-1
		CALL final_trajs%atomo(i)%INITIALIZE(k,3)
		final_trajs%atomo(i)%r(1) = trajs%atomo(k)%r(1)
		final_trajs%atomo(i)%r(2) = trajs%atomo(k)%r(2)
		final_trajs%atomo(i)%r(3) = trajs%atomo(k)%r(3)
		final_trajs%atomo(i)%p(1) = trajs%atomo(k)%p(1)
		final_trajs%atomo(i)%p(2) = trajs%atomo(k)%p(2)
		final_trajs%atomo(i)%p(3) = trajs%atomo(k)%p(3)
		final_trajs%atomo(i)%E = trajs%atomo(k)%E
		final_trajs%atomo(i)%init_r(1) = trajs%atomo(k)%init_r(1)
		final_trajs%atomo(i)%init_r(2) = trajs%atomo(k)%init_r(2)
		final_trajs%atomo(i)%init_r(3) = trajs%atomo(k)%init_r(3)
		final_trajs%atomo(i)%init_p(1) = trajs%atomo(k)%init_p(1)
		final_trajs%atomo(i)%init_p(2) = trajs%atomo(k)%init_p(2)
		final_trajs%atomo(i)%init_p(3) = trajs%atomo(k)%init_p(3)
	END DO
	DEALLOCATE(trajs%atomo)
   ! Print if the option was given
	IF(inicondat%control_out.EQV..TRUE.) THEN
		OPEN(11,FILE=inicondat%output_file,STATUS="replace")
		WRITE(11,*) "# FILE CREATED BY : GENERATE_TRAJS_ATOMS ===============================================================" 
		WRITE(11,*) "# Format:   traj_num      X,Y,Z (a.u.)      Px,Py,Pz(a.u.)   X,Y(Projected IWS, a.u.)"
		WRITE(11,*) "# Initial total Energy (a.u.) / (eV) : ", final_trajs%atomo(1)%E, " /  ", final_trajs%atomo(1)%E*au2ev
		WRITE(11,*) "# Perpendicular Energy (a.u.) / (eV) : ", inicondat%E_norm%getvalue(), " / ", inicondat%E_norm%getvalue()*au2ev
		WRITE(11,*) "# MASS (a.u.) / proton_mass : ", inicondat%masss%getvalue()," / ", inicondat%masss%getvalue()/pmass2au
		WRITE(11,*) "# Incidence angle (deg): ", inicondat%vz_angle%getvalue()*180.D0/PI
		WRITE(11,*) "# Parallel velocity direction (deg): ", inicondat%vpar_angle%getvalue()*180.D0/PI
		IF ((inicondat%control_posX).AND.(.NOT.inicondat%control_posY)) THEN
			WRITE(11,*) "# Random X impact parameter"
			WRITE(11,*) "# Seed used: ", inicondat%seed
		ELSE IF ((.NOT.inicondat%control_posX).AND.(inicondat%control_posY)) THEN
			WRITE(11,*) "# Random Y impact parameter"
			WRITE(11,*) "# Seed used: ", inicondat%seed
		ELSE IF ((inicondat%control_posX).AND.(inicondat%control_posY)) THEN
			WRITE(11,*) "# Random X and Y impact parameters"
			WRITE(11,*) "# Seed used: ", inicondat%seed
		ELSE
			WRITE(11,*) "# X, Y values are not random numbers"
		END IF
		WRITE(11,*) "#======================================================================================================="
		DO i=1, inicondat%ntraj-inicondat%nstart+1
			FORALL(j=1:2) proj_iws_r(j) = final_trajs%atomo(i)%init_r(j)
         proj_iws_r = thispes%surf%project_iwscell(proj_iws_r)
			WRITE(11,'(1X,I10,3(3F20.7))') final_trajs%atomo(i)%id, final_trajs%atomo(i)%init_r, final_trajs%atomo(i)%init_p, proj_iws_r
		END DO
		CLOSE(11)
#ifdef DEBUG
		CALL VERBOSE_WRITE(routinename,"Outputfile generated: ",inicondat%output_file)
#endif
	END IF
	RETURN
END SUBROUTINE GENERATE_TRAJS_ATOMS
END MODULE INITATOM_MOD
