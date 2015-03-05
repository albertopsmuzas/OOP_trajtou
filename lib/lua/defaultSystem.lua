--////////////////////////////////////////////////
--       DEFAULT SYSTEM INFORMATION INPUT
--////////////////////////////////////////////////
--
-- Initializes all system variables with their
-- default values
--
--***********************************************
local SYSTEM={
	debugMode=false,    								-- Controls debug mode (more tedious than verbosemode)
	verboseMode=false,  								-- Controls verbosity
	surface='None',                     		-- Name of surface input file, if any.
	symbols={},         								-- Table of atomic symbols
	masses={},          								-- Table of magnitudes (mass units)
	surfaceOscillator={
    wx=0.0,wy=0.0,wz=0.0,
    mass={0.0,'au'},
  },
  binning={type='Morse',dissociationEnergy=0.0,equilibriumDistance=0.0,width=0.0}, -- Morse is the only one implemented
	pesPath=os.getenv('OOPTRAJTOUPES')..'/'	-- Path to search for data input files
}
return SYSTEM
