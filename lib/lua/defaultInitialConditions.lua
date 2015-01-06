--///////////////////////////////////////////////////////
--         INITIAL CONDITIONS GENERAL INPUT
--//////////////////////////////////////////////////////
local initialConditions={
	-- For 3D initial conditions
	kind='Atoms or Molecules',
	trajList={from=0,to=0};
	Enormal={3.e0,'ev'},
	incidenceAngle={0.e0,'rad'},
	directionAngle={0.e0,'rad'},
	initialZ={0.e0,'au'},
	randomXY={X={true,0.e0},Y={true,0.e0}},
	outputFile={'OUTinicond.out'},
	seedRead=true,
	-- For 6D initial conditions
	vacuumFunction={kind='Numerical',source='filename'},
	vibrationalState={v=0,J=0},
	internalEnergy={0.e0,'au'},
	integrator={
		timeStep={0.e0,'au'},
		precision={0.e0,'au'},
		extrapolation='Rational',
	},
}
return initialConditions
