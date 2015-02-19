--/////////////////////////////////////////////////////
--              DEFAULT DYNAMICS INPUT
--/////////////////////////////////////////////////////
--
-- Must include all parameters and info needed for a 
-- dynamics job (3D and 6D)
--
--*****************************************************

-- dynamics block
local dynamics={
	kind="Atoms or Molecules",
	extrapolation='Rational',
	scaling='Equal',
	timeStep={0.e0,'au'},
	maxTime={0.e0,'au'},
	stopAtZ={-100.e0,'au'},
	stopAtZ_dZ={0.e0,'au'},
	follow={},
	outTrajConditions={
		reflection={},
		adsorption={},
		absorption={},
		dissociation={},
		trappedAfter=1000,
	},
}
return dynamics
