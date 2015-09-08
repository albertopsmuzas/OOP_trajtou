!##################################################
! PROGRAM: GETPOT_CRP6D
!> @brief
!! Test if system part of Lua conf file is read correctly
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!##################################################
PROGRAM GETPOT_CRP6D
! Initial declarations
USE DEBUG_MOD
USE CRP6D_MOD
USE SYSTEM_MOD
! use some modules?
IMPLICIT NONE
! Variables
REAL(KIND=4),DIMENSION(2):: timearr
REAL(KIND=8),DIMENSION(6):: x,dvdx
REAL(KIND=8):: v
TYPE(CRP6D):: thispes
REAL(KIND=4):: timer
CHARACTER(LEN=1024):: luafile
CHARACTER(LEN=1024):: auxstring
integer(kind=4):: i,loop
! GABBA GABBA HEY! ===============================
!
! STEP 1: INITIALIZE SYSTEM VIA LUA CONFIG FILE
SELECT CASE(command_argument_count())
   CASE(8)
      CALL GET_COMMAND_ARGUMENT(1,luafile)
      CALL GET_COMMAND_ARGUMENT(2,auxstring)
      READ(auxstring,*) x(1)
      CALL GET_COMMAND_ARGUMENT(3,auxstring)
      READ(auxstring,*) x(2)
      CALL GET_COMMAND_ARGUMENT(4,auxstring)
      READ(auxstring,*) x(3)
      CALL GET_COMMAND_ARGUMENT(5,auxstring)
      READ(auxstring,*) x(4)
      CALL GET_COMMAND_ARGUMENT(6,auxstring)
      READ(auxstring,*) x(5)
      CALL GET_COMMAND_ARGUMENT(7,auxstring)
      READ(auxstring,*) x(6)
      call get_command_argument(8,auxString)
      read(auxString,*) loop
   CASE DEFAULT
      WRITE(0,*) "ERR: Bad number of arguments: ",command_argument_count()
      WRITE(0,*) "Seven arguments needed: lua config file, 6 real numbers"
      CALL EXIT(1)
END SELECT
CALL ETIME(timearr,timer)
CALL INITIALIZE_SYSTEM(trim(luafile))
CALL thispes%INITIALIZE()
do i=1,loop
   CALL thispes%GET_V_AND_DERIVS(x,v,dvdx)
   CALL VERBOSE_WRITE('Potential calculated at:',x(:))
   WRITE(*,*) "Potential (a.u.): ", v
   WRITE(*,*) "Derivatives (X,Y,Z,r,theta,phi) (a.u./rad.): ", dvdx(:)
enddo
CALL ETIME(timearr,timer)
! STEP 2: GET VALUES
!CALL VERBOSE_WRITE('Potential calculated at:',x(:))
!WRITE(*,*) "Potential (a.u.): ", v
!WRITE(*,*) "Derivatives (X,Y,Z,r,theta,phi) (a.u./rad.): ", dvdx(:)

! STEP 3: PRINT TIME 
CALL VERBOSE_WRITE("****************** RUN TIME ***************************")
CALL VERBOSE_WRITE('',"User time: ",real(timearr(1),kind=8))
CALL VERBOSE_WRITE('',"System time: ",real(timearr(2),kind=8))
CALL VERBOSE_WRITE('',"Total time: ",real(timer,kind=8))
CALL VERBOSE_WRITE("******************************************************")
END PROGRAM GETPOT_CRP6D
