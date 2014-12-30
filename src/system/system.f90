!#########################################################
! MODULE: SYSTEM_MOD
!> @brief
!! Public module which compiles interesting variables to 
!! keep during runtime. Debug variables are separated in the
!! DEBUG_MOD
!##########################################################
MODULE SYSTEM_MOD
   USE UNITS_MOD
   USE AOTUS_MODULE, ONLY: flu_State, OPEN_CONFIG_FILE, CLOSE_CONFIG, AOT_GET_VAL
   USE AOT_TABLE_MODULE, ONLY: AOT_TABLE_OPEN, AOT_TABLE_CLOSE, AOT_TABLE_LENGTH
#ifdef DEBUG
   USE DEBUG_MOD
#endif
! Initial declarations
IMPLICIT NONE
! Global variables
REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: system_mass
CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE:: system_atomsymbols
CHARACTER(LEN=:),ALLOCATABLE:: system_inputfile
CHARACTER(LEN=:),ALLOCATABLE:: system_pespath
INTEGER(KIND=4):: system_dimensions
INTEGER(KIND=4):: system_natoms
CONTAINS
!###########################################################
!# SUBROUTINE: INITIALIZE_SYSTEM
!###########################################################
!> @brief
!! Load system parameters
!
!> @author A.S. Muzas - alberto.muzas@uam.es
!> @date Dec/2014
!> @version 1.0
!-----------------------------------------------------------
SUBROUTINE INITIALIZE_SYSTEM(filename)
   ! Initial declarations
   IMPLICIT NONE
   ! I/O variables
   CHARACTER(LEN=*),INTENT(IN):: filename
   ! Local variables
   TYPE(flu_State):: conf ! Lua file
   INTEGER(KIND=4):: ierr
   INTEGER(KIND=4):: sys_table,sym_table,mass_table
   INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE:: subtable
   CHARACTER(LEN=*),PARAMETER:: routinename="INITIALIZE_SYSTEM: "
   TYPE(Mass):: masa
   REAL(KIND=8):: numero
   CHARACTER(LEN=10):: units
   INTEGER(KIND=4):: i ! counters
   CHARACTER(LEN=1024):: string
   ! Run section
   ALLOCATE(system_inputfile,source=filename)
   CALL OPEN_CONFIG_FILE(L=conf,filename=system_inputfile,ErrCode=ierr)
   SELECT CASE(ierr)
      CASE(0)
         ! do nothing
      CASE DEFAULT
         WRITE(0,*) "INITIALIZE_SYSTEM: error opening Lua config file"
         CALL EXIT(1)
   END SELECT
   ! Open table system
   CALL AOT_TABLE_OPEN(L=conf,thandle=sys_table,key='system')
   ! Get system dimensiopns
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=sys_table,key='dimensions',val=system_dimensions)
   ! Get number of atoms, masses and atomic symbols
   CALL AOT_TABLE_OPEN(L=conf,parent=sys_table,thandle=sym_table,key='symbols')
   CALL AOT_TABLE_OPEN(L=conf,parent=sys_table,thandle=mass_table,key='masses')
   system_natoms=aot_table_length(L=conf,thandle=sym_table)
   SELECT CASE(system_natoms==aot_table_length(L=conf,thandle=mass_table) .AND. system_natoms/=0)
      CASE(.TRUE.)
         ! do nothing
      CASE(.FALSE.)
         WRITE(0,*) "INITIALIZE_SYSTEM ERR: dimensions mismatch (or zero) between tables system.masses and system.symbols"
         CALL EXIT(1)
   END SELECT
   ALLOCATE(system_mass(system_natoms))
   ALLOCATE(system_atomsymbols(system_natoms))
   ALLOCATE(subtable(system_natoms))
   DO i = 1, system_natoms
      CALL AOT_TABLE_OPEN(L=conf,parent=mass_table,thandle=subtable(i),pos=i)
      CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=subtable(i),pos=1,val=numero)
      CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=subtable(i),pos=2,val=units)
      CALL masa%READ(numero,trim(units))
      CALL masa%TO_STD()
      system_mass(i)=masa%getvalue()
      CALL AOT_TABLE_CLOSE(L=conf,thandle=subtable(i))
      CALl AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=sym_table,pos=i,val=system_atomsymbols(i))
   END DO
   CALL AOT_TABLE_CLOSE(L=conf,thandle=mass_table)
   CALL AOT_TABLE_CLOSE(L=conf,thandle=sym_table)
   CALL AOT_GET_VAL(L=conf,ErrCode=ierr,thandle=sys_table,key='pesPath',val=string)
   ALLOCATE(CHARACTER(LEN=len(trim(string))):: system_pespath)
   system_pespath=trim(string)
   CALL CLOSE_CONFIG(conf)
#ifdef DEBUG
   CALL VERBOSE_WRITE(routinename,'config file: ',system_inputfile)
   CALL VERBOSE_WRITE(routinename,'default PES path: ',system_pespath)
   CALL VERBOSE_WRITE(routinename,'Natoms: ',system_natoms)
   CALl VERBOSE_WRITE(routinename,'dimensions: ',system_dimensions)
   CALL VERBOSE_WRITE(routinename,'masses(au): ',system_mass(:))
   CALL VERBOSE_WRITE(routinename,'atomic symbols: ',system_atomsymbols(:))
#endif
   RETURN
END SUBROUTINE INITIALIZE_SYSTEM
END MODULE SYSTEM_MOD
