PROGRAM ANALYSIS
USE DEBUG_MOD
USE DIFFRACTIONCRP6D_MOD
TYPE(Allowed_peaksCRP6D) :: this
CALL SET_VERBOSE_MODE(.TRUE.)
CALL this%INITIALIZE("INsurface.inp","INinicond6d.inp")
CALL this%SETUP()
CALL this%ASSIGN_PEAKS()
CALL this%PRINT_LABMOMENTA_AND_ANGLES()
END PROGRAM ANALYSIS
