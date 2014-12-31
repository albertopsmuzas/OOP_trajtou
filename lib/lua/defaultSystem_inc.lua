--////////////////////////////////////////////////
--       DEFAULT SYSTEM INFORMATION INPUT
--////////////////////////////////////////////////
--
-- Initializes all system variables with their
-- default values
--
--***********************************************
system={
	debugMode=false,    						-- Controls debug mode (more tedious than verbosemode)
	verboseMode=false,  						-- Controls verbosity
	surface='None',                     -- Name of surface input file, if any.
	symbols={},         						-- Table of atomic symbols
	masses={},          						-- Table of magnitudes (mass units)
	pesPath=os.getenv('OOPTRAJTOUPES')	-- Path to search for data input files
}
