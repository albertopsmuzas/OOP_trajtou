!##################################################
! PROGRAM: molec2atom
!> @brief
!! Program to go from atomic coordinates X1,Y1,Z1,X2,Y2,Z2
!! fo molecular coordinates: X,Y,Z,R,Theta,Phi
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Apr/2014 
!> @version 1.0
!##################################################
PROGRAM atom2molec
! Initial declarations
use DEBUG_MOD, only: VERBOSE_WRITE
use SYSTEM_MOD, only: INITIALIZE_SYSTEM,system_mass
use CRP6D_MOD, only: FROM_ATOMIC_TO_MOLECULAR
IMPLICIT NONE
! variables
CHARACTER(LEN=1024):: auxstring,luafile
REAL(KIND=8),DIMENSION(6):: molec_coord, atomic_coord
! Run yeah, tun
SELECT CASE(command_argument_count())
   CASE(7)
      CALL GET_COMMAND_ARGUMENT(1,luafile)
      CALL GET_COMMAND_ARGUMENT(2,auxstring)
      READ(auxstring,*) atomic_coord(1)
      CALL GET_COMMAND_ARGUMENT(3,auxstring)
      READ(auxstring,*) atomic_coord(2)
      CALL GET_COMMAND_ARGUMENT(4,auxstring)
      READ(auxstring,*) atomic_coord(3)
      CALL GET_COMMAND_ARGUMENT(5,auxstring)
      READ(auxstring,*) atomic_coord(4)
      CALL GET_COMMAND_ARGUMENT(6,auxstring)
      READ(auxstring,*) atomic_coord(5)
      CALL GET_COMMAND_ARGUMENT(7,auxstring)
      READ(auxstring,*) atomic_coord(6)
   CASE DEFAULT
      WRITE(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      CALL EXIT(1)
END SELECT
CALL INITIALIZE_SYSTEM(trim(luafile))
CALL VERBOSE_WRITE("***********************************************")
CALL VERBOSE_WRITE("** FROM ATOMIC COORDINATES TO MOLECULER ONES **")
CALL VERBOSE_WRITE("***********************************************")
CALL FROM_ATOMIC_TO_MOLECULAR(system_mass(1),system_mass(2),atomic_coord,molec_coord)
WRITE(*,*) molec_coord(:)
END PROGRAM atom2molec
