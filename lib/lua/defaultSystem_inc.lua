--////////////////////////////////////////////////
--       DEFAULT SYSTEM INFORMATION INPUT
--////////////////////////////////////////////////
--
-- Initializes all system variables with their
-- default values
--
--***********************************************
system={
	dimensions=0,		  						-- Dimensions fo the system
	debugMode=false,    						-- Controls debug mode (more tedious than verbosemode)
	verboseMode=false,  						-- Controls verbosity
	symbols={},         						-- Table of atomic symbols
	masses={},          						-- Table of magnitudes (mass units)
	pesPath=os.getenv('OOPTRAJTOUPES')	-- Path to search for PES files
}
