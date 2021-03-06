--/////////////////////////////////////////////////////////////
--            INITIAL CONDITIONS GENERAL INPUT
--/////////////////////////////////////////////////////////////
-- There should be an entry for every posible input
-- parameter. Values writen here act as default initialization
-- values
-- ************************************************************
local initialConditions={
	-- For 3D-atom initial conditions
	-- ******************************
	kind='Atoms or Molecules',
	trajList={from=0,to=0};
	Enormal={0,'ev'},
	incidenceAngle={0.e0,'rad'},
	directionAngle={0.e0,'rad'},
	initialZ={0.e0,'au'},
	randomXY={X={true,0.e0},Y={true,0.e0}},
	outputFile={'OUTinicond.out'},
	seedRead=true,

	-- For 6D-atom + surf conditions
	-- *****************************
	sufaceTemperature={0.0,'Kelvin'},

	-- For 6D-molecule initial conditions
	-- **********************************
	vibrationalFunction={kind='Numerical',source='filename'},
	internalState={v=0,J=0,mJ='Average'},
	internalEnergy={0.0,'au'},
	timeStep={0.0,'au'},
	pecision=1e-6,
	extrapolation='Rational',
}
return initialConditions
