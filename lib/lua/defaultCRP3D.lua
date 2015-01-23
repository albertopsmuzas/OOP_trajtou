--///////////////////////////////////////////////
--           DEFAULT CRP3D PES INPUT
--///////////////////////////////////////////////
local CRP3D={
      kind='CRP3D',                          -- CRP3D kind of PES
      name='default',                        -- PES files will be at $OOPTRAJTOUPES/'name'
      dimensions=3,                          -- Number of dimensions
      surfaceInput='INsurface.inp',          -- Surface input filename
      maxEnvironment=0,                      -- Max environment involved in corrugation extraction
      dampFunction={kind='None',param={}},   -- Damp function used to soften PES at far distances
      pairPotentials={},                     -- List of file names
      sitios={},										--	List of file names
      fourierKpoints={},                     -- List of integer points, e.g.: {0,0}, {0,1} ...
}
return CRP3D
