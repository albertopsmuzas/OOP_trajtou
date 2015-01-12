--///////////////////////////////////////////////////////
--         INITIAL CONDITIONS GENERAL INPUT
--//////////////////////////////////////////////////////
local initialConditions={
	-- For 3D initial conditions
	kind='Atoms or Molecules',
	trajList={from=0,to=0};
	Enormal={0,'ev'},
	incidenceAngle={0.e0,'rad'},
	directionAngle={0.e0,'rad'},
	initialZ={0.e0,'au'},
	randomXY={X={true,0.e0},Y={true,0.e0}},
	outputFile={'OUTinicond.out'},
	seedRead=true,
	-- For 6D initial conditions
	vibrationalFunction={kind='Numerical',source='filename'},
	internalState={v=0,J=0},
	internalEnergy={0.0,'au'},
	timeStep={0.0,'au'},
	extrapolation='Rational',
}
return initialConditions
