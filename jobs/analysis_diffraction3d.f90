PROGRAM ANALYSIS
USE DEBUG_MOD
USE DIFFRACTION_MOD
TYPE(Allowed_peaks) :: this
CALL this%INITIALIZE()
CALL this%SETUP()
CALL this%ASSIGN_PEAKS()
END PROGRAM ANALYSIS
