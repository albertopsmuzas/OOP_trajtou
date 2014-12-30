--{{ include default specifications
require("defaultDynamics_inc")
--}}

PES.name="INcrp3d.inp"
initialConditions.name="INinicond3d.inp"

--{{ dynamics block
	dynamics.kind="Atom"
	dynamics.precision=1.e-6
	dynamics.timeStep:Read(.5e0,"angst")
	dynamics.maxTime:Read(15.e0,"ps")
--}} end

--{{ outTrajConditions block
	outTrajConditions.reflec:Read(5.e0,"angst")
	outTrajConditions.adsorp:Read(3.e0,"angst")
	outTrajConditions.absorp:Read(-0.5e0,"angst")
--}} end
