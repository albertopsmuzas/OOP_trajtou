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
	precision=0.e0,
	scaling="Equal",           -- Avail: Equal, Smart
	extrapolation="Rational",  -- Avail: Rational, Polinomi
	timeStep={0.e0,'au'},
	maxTime={0.e0,'au'},
	stopAtZ={-100.e0,'au'},
	stopAtZ_dZ={0.e0,'au'},
	integrator={
		dt={0.0,'au'},
	},
	follow={},
	outTrajConditions={
		reflec={},
		adsorp={},
		absorp={},
	},
}
return dynamics
