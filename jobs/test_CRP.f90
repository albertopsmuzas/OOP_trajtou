PROGRAM CRP_TEST
! Initial declarations
USE DEBUG_MOD
USE CRP_PES
IMPLICIT NONE
! Variables
TYPE(CRP) :: crp_pes
CHARACTER(LEN=30), DIMENSION(:), ALLOCATABLE :: site_aliases, site_filenames, pp_aliases, pp_filenames
CHARACTER(LEN=30) :: surf_filename, surf_alias, filename
INTEGER(KIND=4) :: max_order
INTEGER(KIND=4) :: i ! counters
! Run section===================================================================================0
! STEP 0: HELLO!
WRITE(*,*) "***************************************" 
WRITE(*,*) "TEST_PES program executed"
WRITE(*,*) "***************************************" 
WRITE(*,*) "Set max_order: (INTEGER)"
READ(*,*) max_order
CALL SET_DEBUG_MODE(.FALSE.)
CALL SET_VERBOSE_MODE(.FALSE.)
! STEP 1: INITIALIZE FILENAMES AND ALIASES -------------------------------------------------------
npairpots=2
nsites=6
ALLOCATE (site_aliases(1:nsites))
ALLOCATE (site_filenames(1:nsites))
ALLOCATE (pp_aliases(1:npairpots))
ALLOCATE (pp_filenames(1:npairpots))
site_aliases = [CHARACTER(LEN=30) :: "Top_Li","Top_F","Hollow","Bridge","Half_hollow-F","Half_hollow-Li"]
site_filenames = [CHARACTER(LEN=30) :: "Sitio1.dat","Sitio2.dat","Sitio3.dat","Sitio4.dat","Sitio5.dat","Sitio6.dat"]
pp_aliases = [CHARACTER(LEN=30) :: "Li_repul","F_repul"]
pp_filenames = [CHARACTER(LEN=30) :: "intrep1.dat", "intrep2.dat"]
! STEP 2: INITIALIZE CRP PES:
CALL crp_pes%INITIALIZE(npairpots,pp_filenames,pp_aliases,nsites,site_filenames,site_aliases,surf_filename,surf_alias,max_order)
! STEP 3: DO Z INTERPOLATION EXTRACTING VASINT AND SMOOTHING SITES
CALL crp_pes%EXTRACT_VASINT()
CALL crp_pes%SMOOTH()
CALL crp_pes%INTERPOL_Z()
! STEP 4: PLOT SOME GRAPHS
DO i = 1, nsites
   WRITE(filename,'(I2)') i
   filename=trim(filename)
   filename="Site_data_"//filename//".dat"
   CALL crp_pes%all_sites(i)%PLOT_DATA(filename)
END DO
DO i = 1, npairpots
   WRITE(filename,'(I2)') i
   filename=trim(filename)
   filename="Pairpot_data_"//filename//".dat"
   CALL crp_pes%all_sites(i)%PLOT_DATA(filename)
END DO
END PROGRAM CRP_TEST
