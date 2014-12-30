--/////////////////////////////////////////////////////
--              DEFAULT DYNAMICS INPUT
--/////////////////////////////////////////////////////
--
-- Must include all parameters and info needed for a 
-- dynamics job.
--
--*****************************************************

-- dynamics block
dynamics={
	kind="default",            -- Avail: Atoms, Molec
	precision=0.e0,
	scaling="Equal",           -- Avail: Equal, Smart
	extrapolation="Rational",  -- Avail: Rational, Polinomi
	timeStep={0.e0,'au'},
	maxTime={0.e0,'au'},
	stopAtZ={-100.e0,'au'},
	stopAtZ_dZ={0.e0,'au'},
	integrator={},
	follow={},
}

-- outTrajConditions block
outTrajConditions={
	reflec={},
	adsorp={},
	absorp={},
}
