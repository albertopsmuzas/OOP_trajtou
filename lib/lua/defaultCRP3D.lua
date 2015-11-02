--///////////////////////////////////////////////
--           DEFAULT CRP3D PES INPUT
--///////////////////////////////////////////////
local CRP3D={
      kind='CRP3D',                                            -- CRP3D kind of PES
      name='default',                                          -- PES files will be at $OOPTRAJTOUPES/'name'
      dimensions=3,                                            -- Number of dimensions
      surfaceInput='INsurface.inp',                            -- Surface input filename
      maxEnvironment=0,                                        -- Max environment involved in corrugation extraction
      dampFunction={kind='None',param={}},                     -- Damp function used to soften PES at far distances
      pairPotentials={},                                       -- List of file names
      sitios={},										                  --	List of file names
      fourierKpoints={},                                       -- List of integer 2D points, e.g.: {0,0}, {0,1} ...
		-- resize={r=integer,z=integer},                         -- Optional. Resizes ZR grid so that derivatives are smoother
		dampFunction={kind='None',param={'list of parameters'}}, -- Pair potentials may be multiplied by a control function. 
}
return CRP3D
