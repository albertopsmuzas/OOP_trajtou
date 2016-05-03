--/////////////////////////////////////////////////////
--              DEFAULT NEWINPUT INPUT
--/////////////////////////////////////////////////////
--
-- Must include all parameters and info needed 
-- to create new 3D input files: sitios, pairpots, etc.
--
--*****************************************************
-- New input block
local newInput={
  vacuumPot={0.0,'au'},
  rumplingList={units='au'},
  gridInfo={nPoints=0,units='au',kind='Manual',
    points={},
  }, 
}
-- Output newInput
return newInput