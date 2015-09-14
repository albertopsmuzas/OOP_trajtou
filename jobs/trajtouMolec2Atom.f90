!##################################################
! PROGRAM: molec2atom
!> @brief
!! Program to go from molecular coordinates: X,Y,Z,R,Theta,Phi
!! to atomic coordinates X1,Y1,Z1,X2,Y2,Z2
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014 
!> @version 1.0
!##################################################
PROGRAM molec2atom
! Initial declarations
use DEBUG_MOD, only: VERBOSE_WRITE
use SYSTEM_MOD, only: INITIALIZE_SYSTEM,from_molecular_to_atomic_phaseSpace
IMPLICIT NONE
! variables
CHARACTER(LEN=1024):: auxstring,luafile
REAL(KIND=8),DIMENSION(12):: molCoordPhaseSpace
! Run yeah, tun
SELECT CASE(command_argument_count())
   CASE(13)
      CALL GET_COMMAND_ARGUMENT(1,luafile)
      CALL GET_COMMAND_ARGUMENT(2,auxstring)
      READ(auxstring,*) molCoordPhaseSpace(1)
      CALL GET_COMMAND_ARGUMENT(3,auxstring)
      READ(auxstring,*) molCoordPhaseSpace(2)
      CALL GET_COMMAND_ARGUMENT(4,auxstring)
      READ(auxstring,*) molCoordPhaseSpace(3)
      CALL GET_COMMAND_ARGUMENT(5,auxstring)
      READ(auxstring,*) molCoordPhaseSpace(4)
      CALL GET_COMMAND_ARGUMENT(6,auxstring)
      READ(auxstring,*) molCoordPhaseSpace(5)
      CALL GET_COMMAND_ARGUMENT(7,auxstring)
      READ(auxstring,*) molCoordPhaseSpace(6)
      CALL GET_COMMAND_ARGUMENT(8,auxstring)
      READ(auxstring,*) molCoordPhaseSpace(7)
      CALL GET_COMMAND_ARGUMENT(9,auxstring)
      READ(auxstring,*) molCoordPhaseSpace(8)
      CALL GET_COMMAND_ARGUMENT(10,auxstring)
      READ(auxstring,*) molCoordPhaseSpace(9)
      CALL GET_COMMAND_ARGUMENT(11,auxstring)
      READ(auxstring,*) molCoordPhaseSpace(10)
      CALL GET_COMMAND_ARGUMENT(12,auxstring)
      READ(auxstring,*) molCoordPhaseSpace(11)
      CALL GET_COMMAND_ARGUMENT(13,auxstring)
      READ(auxstring,*) molCoordPhaseSpace(12)
   CASE DEFAULT
      WRITE(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      CALL EXIT(1)
END SELECT
CALL INITIALIZE_SYSTEM(trim(luafile))
CALL VERBOSE_WRITE("***********************************************")
CALL VERBOSE_WRITE("** FROM MOLECULAR COORDINATES TO ATOMIC ONES **")
CALL VERBOSE_WRITE("***********************************************")
WRITE(*,*) from_molecular_to_atomic_phaseSpace(molCoordPhaseSpace(:))
END PROGRAM molec2atom
