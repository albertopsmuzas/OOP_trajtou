--/////////////////////////////////////////////////////
--              DEFAULT NEWINPUT INPUT
--/////////////////////////////////////////////////////
--
-- Must include all parameters and info needed 
-- to create new 3D input files: sitios, pairpots, etc.
--
--*****************************************************
-- New input block
local newInput3d={
  vacuumPot={0.0,'au'},
  rumplingList={units='au'},
  gridInfo={nPoints=0,units='au',kind='Manual',points={}}, 
  pairpotFiles={}, -- should be given with the same order as in rumplingList
  -- example: 
  -- pairpotFiles={ {inp='fooRaw.dat',out='foo.dat'}, {inp='foo2Raw.dat',out='foo2.dat'}, etc. }
  sitioFiles={},
  -- example: 
  -- sitioFiles={ {inp='fooRaw.dat',out='foo.dat'}, {inp='foo2Raw.dat',out='foo2.dat'}, etc. }
}
-- Output newInput
return newInput3d