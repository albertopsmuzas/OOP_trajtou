--///////////////////////////////////////////////
--           DEFAULT CRP6D PES INPUT
--///////////////////////////////////////////////
pes={
		kind='CRP6D',                            -- CRP3D kind of PES
		name='default',                          -- PES files will be at $OOPTRAJTOUPES/'name'
		dimensions=6,                            -- Number of dimensions
		crp3dPes={},                             -- CRP3D tables needed.
		surfaceInput='surfacefile',              -- Surface input filename
		resize={r=0,z=0},                               -- If we want a new denser grid for ZR-cuts
		dampFunction={kind='None',param={}},              -- Damp funtion used to soften CRP3D additions & substractions at far distances
		extrapolFunction={kind='None',upToZ={0.e0,'au'}}, -- Function used to extrapolate to the vacuum. For Z>upToZ, we have vacuum potential
		vacuumFunction={kind='Numerical',source='None'},  -- Source for vacuum potential
		fourierKpoints={},                                -- List of integer points, e.g.: {0,0}, {0,1} ...
}
pes.wyckoffSite={}
--************************************************************************
-- EXAMPLE OF HOW TO DEFINE A pes.wyckoffSite
--
-- pes.wyckoffSite[1]={
-- 	kind='wyckoffSymbol',
-- 	name='some string',
--		{theta=0.e0,files={'some list of files'}},
-- }
