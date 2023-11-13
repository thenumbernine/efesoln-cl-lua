#!/usr/bin/env luajit
local ffi = require 'ffi'
local class = require 'ext.class'
local math = require 'ext.math'
local table = require 'ext.table'
local path = require 'ext.path'
local template = require 'template'
local vec3d = require 'vec-ffi.vec3d'
local gl = require 'gl'
local CLEnv = require 'cl.obj.env'
local clnumber = require 'cl.obj.number'
-- parameters:

-- helper for indexing symmetric matrices
local function sym(i,j)
	if i <= j then return i..j else return j..i end
end

-- constants

local c = 299792458			-- m/s 
local G = 6.67384e-11		-- m^3 / (kg s^2)

-- spherical body, no charge, no velocity

local SphericalBody = class()

function SphericalBody:init(args)
	self.radius = args.radius
	self.mass = args.mass
	self.volume = 4/3 * math.pi * self.radius * self.radius * self.radius	-- m^3
	self.density = self.mass / self.volume	-- 1/m^2

	self.useMatter = true
	self.useVel = false
	self.useEM = false

	self.init = template([[
	real3 const x = getX(i);
	real const r = real3_len(x);
	
	real const rho0 = <?=self.density?>;
	real const radius = <?=self.radius?>;
	real const radius3 = radius * radius * radius;
	real const mass = <?=self.mass?>;
	real const r2 = r * r;
	TPrim->rho = r < radius ? rho0 : 0;
	TPrim->P = r < radius ? (rho0 * (
		(sqrt(1. - 2. * mass * r2 / radius3) - sqrt(1. - 2. * mass / radius))
		/ (3. * sqrt(1. - 2. * mass / radius) - sqrt(1. - 2. * mass * r2 / radius3))
	)) : 0;

]], {self=self})
end

local EMRing = class()

function EMRing:init(args)
	self.radius = args.radius
	
	self.useMatter = false
	self.useVel = false
	self.useEM = true
	
	self.init = template([[
	real const radius = <?=self.radius?>;
	
	real3 const x = getX(i);
	
	real const polar_rSq = x.x*x.x + x.y*x.y;
	real const polar_r = sqrt(polar_rSq);		//r in polar coordinates	
	real const dr = polar_r - radius;			//difference from polar radius to torus big radius
	real const r = sqrt(x.z*x.z + dr*dr);		//r in torus radial coordinates
	real const theta = atan2(x.z, dr);		//angle around the small radius
	real const phi = atan2(x.x, x.y);			//angle around the big radius

	//F^uv_;v = -4 π J^u
	// means that the divergence of the EM is the 4-current 
	//the divergence of the exterior of the 4-potential is the 4-current
	//so if the 4-current is a Dirac delta function along the line in space where there is current
	// then the EM tensor is going to be an inverse falloff around it

	//4-current
	//t is current density
	//i is charge density = current density * drift velocity
	//stressEnergyPrims.chargeDensity
	//stressEnergyPrims.A

	/*
	point on the surface:
		r * cos(phi) * cos(theta)
		r * sin(phi) * cos(theta)
		r * sin(theta)
	*/

	TPrim->E.x = -x.y / polar_rSq;
	TPrim->E.y = x.x / polar_rSq;
	TPrim->E.z = 0;

	TPrim->B.x = cos(theta) / r * cos(phi);
	TPrim->B.y = cos(theta) / r * sin(phi);
	TPrim->B.z = -sin(theta) / r;

]], {self=self})
end

local EMLine = class()

function EMLine:init(args)
	self.radius = args.radius
	
	self.useMatter = false
	self.useVel = false
	self.useEM = true
	
	-- from cpu-test/grem.lua
	local s_in_m = 1 / c
	local kg_in_m = G / c^2
	local ke = 8.9875517873681764e+9
	local C_in_m = math.sqrt(ke * G) / c^2	-- m
	local N_in_m = kg_in_m / s_in_m^2	-- m^0
	local V_in_m = N_in_m / C_in_m	-- m^0
	local Ohm_in_m = kg_in_m / (s_in_m * C_in_m^2)	-- m^0
	local mu0 = 4 * math.pi	-- m^0
	local eps0 = 1 / mu0	-- m^0
	local e_const = 6.2415093414e+18
	local e_in_m = C_in_m / e_const
	local h_in_m = 1.61622938e-35
	local alpha_const = (e_in_m/h_in_m)^2
	-- source: http://hyperphysics.phy-astr.gsu.edu/hbase/electric/resis.html
	local wire_resistivities = table{	-- at 20' Celsius, in ohm m
		aluminum = 2.65e-8,
		copper = 1.724e-8,
		iron = 9.71e-8,
		nichrome = 1e-6,
		gold = 2.24e-8,
		silver = 1.59e-8,
		platinum = 1.06e-7,
		tungsten = 5.65e-8,
	}:map(function(v) return v * Ohm_in_m end)	-- ohm m => m^0
	local in_in_m = .0254 
	local wire_diameters = table{	-- starts in inches
		electrical_range = .1019,
		household_circuit = .0808,
		switch_leads = .0640,
	}:map(function(v) return v * in_in_m end)	-- in => m 
	local wire_radius = .5 * wire_diameters.electrical_range	-- m
	local wire_cross_section_area = math.pi * wire_radius^2	-- m^2
	local wire_length = 12 * in_in_m	-- m
	local wire_resistivity = wire_resistivities.gold
	local wire_resistance = wire_resistivity * wire_length / wire_cross_section_area	-- m^0
	local battery_voltage = 1.5 * V_in_m	-- m^0
	local wire_current = battery_voltage / wire_resistance	-- amps = C / s = m / m = m^0, likewise volts = m^0, ohms = m^0, so amps = volts / ohms = m^0
	local wire_charge_density = 0	-- C / m^3 = m^-2
	local wire_charge_density_per_length = wire_charge_density * wire_cross_section_area	-- m^-2 * m^2 = m^0
	local wire_surface_charge_density = 0
	-- https://physics.stackexchange.com/questions/291779/electric-field-outside-wire-with-stationary-current?rq=1 
	local rEr = wire_surface_charge_density * wire_radius / math.sqrt(eps0)	-- m^0 * m / m = m^0
	local Ez = wire_current * wire_resistivity * math.sqrt(eps0)	-- m^0
	-- http://www.ifi.unicamp.br/~assis/Found-Phys-V29-p729-753(1999).pdf
	local rBt = math.sqrt(mu0) * wire_current / (2 * math.pi)	-- m^-1
	local rBz = 0

	self.init = template([[
	real const radius = <?=self.radius?>;
	
	real3 const x = getX(i);
	
	real const r2Sq = x.x*x.x + x.y*x.y;
	real const r2 = sqrt(r2Sq);		//r in polar coordinates	

	real const Er = <?=clnumber(rEr)?> / r2;
	real const Ez = <?=clnumber(Ez)?>;

	TPrim->E.x = x.x/r2 * Er;
	TPrim->E.y = x.y/r2 * Er;
	TPrim->E.z = Ez;

	real const Bt = <?=clnumber(rBt)?> / r2;
	real const Bz = <?=clnumber(Bz)?>;

	TPrim->B.x = -x.y/r2 * Bt;
	TPrim->B.y = x.x/r2 * Bt;
	TPrim->B.z = Bz;

]], {
		self = self,
		clnumber = clnumber,
		rEr = rEr,
		Ez = Ez,
		rBt = rBt,
		Bz = Bz,
	})
end


-- body parameters:

local bodies = {
	vacuum = {
		radius = 1,
	},
	Earth = SphericalBody{
		radius = 6.37101e+6,	-- m
		mass = 5.9736e+24 * G / c / c,	-- m
	},
	Sun = SphericalBody{
		radius = 6.960e+8,	-- m
		mass = 1.9891e+30 * G / c / c,	-- m
	},
	['EM ring'] = EMRing{
		radius = 2,
	},
	['EM line'] = EMLine{
		radius = 2,
	},
	['EM constant'] = {
		useEM = true,
		radius = 2,
		init = [[
	TPrim->E = _real3(1,0,0);
	TPrim->B = _real3(0,1,0);
]],
	},
}

-- initial conditions:

local EFESolver = CLEnv:subclass()

EFESolver.useFourPotential = false

EFESolver.initConds = table{
	{flat = ''},
	{['stellar Schwarzschild'] = [[
	real const radius = <?=solver.body.radius?>;
	real const mass = <?=solver.body.mass?>;
	real const density = <?=solver.body.density?>;
	
	real const matterRadius = (real)min(r, radius);
	real const volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
	real const m = density * volumeOfMatterRadius;	// m^3		
	
	/*
	g_ti = beta_i = 0
	g_tt = -alpha^2 + beta^2 = -alpha^2 = -1 + Rs/r <=> alpha = sqrt(1 - Rs/r)
	g_ij = gamma_ij = delta_ij + x^i x^j / r^2 2M/(r - 2M)		<- but x is upper, and you can't lower it without specifying gamma_ij
	 ... which might be why the contravariant spatial metrics of spherical and cartesian look so similar 
	*/
	/*
	I'm going by MTW box 23.2 eqn 6 d/dt (proper time) = sqrt(1 - R/r) for r > R
		and ( 3/2 sqrt(1 - 2 M / R) - 1/2 sqrt(1 - 2 M r^2 / R^3) ) for r < R
		for M = total mass 
		and R = planet radius 
	*/
	gPrim->alpha = r > radius 
		? sqrt(1 - 2*mass/r)
		: (1.5 * sqrt(1 - 2*mass/radius) - .5 * sqrt(1 - 2*mass*r*r/(radius*radius*radius)));
<?
for i=0,sDim-1 do
?>	gPrim->betaU.s<?=i?> = 0;
<?	for j=i,sDim-1 do
?>	gPrim->gammaLL.s<?=i..j?> = <?if i==j then?>1. + <?end?>x.s<?=i?>/r * x.s<?=j?>/r * 2*m/(r - 2*m);
<?	end
end
?>
	/*
	dr^2's coefficient
	spherical: 1/(1 - 2M/r) = 1/((r - 2M)/r) = r/(r - 2M)
	spherical contravariant: 1 - 2M/r
	cartesian contravariant: delta_ij - x/r y/r 2M/r
	hmm, contravariant terms of cartesian vs spherical look more similar than covariant terms do

	in the OV metric, dr^2's coefficient is exp(2 Lambda) = 1/(1 - 2 m(r) / r) where m(r) is the enclosing mass
	so the contravariant coefficient would be exp(-2 Lambda) = 1 - 2 m(r) / r
	I'm going to do the lazy thing and guess this converts to delta^ij - 2 m(r) x^i x^j / r^3
	*/

#if 0	//rotating about a distance
	/*
	now if we are going to rotate this
	at a distance of L and at an angular frequency of omega
	(not considering relativistic Thomas precession just yet)
	
	this might be a mess, but I'm (1) calculating the change in time as if I were in a frame rotating by L exp(i omega t)
	then (2) re-centering the frame at L exp(i omega t) ... so I can use the original coordinate system
	*/
	real const dr_alpha = r > radius 
		? (mass / (r * r * sqrt(1. - 2. * mass / r))) 
		: (mass * r / (radius * radius * radius * sqrt(1. - 2. * mass * r * r / (radius * radius * radius))));
	real const dr_m = r > radius ? 0 : (4. * M_PI * r * r * density);
	MetricPrims & dt_metricPrims = dt_metricPrimGrid(index);
	real const L = 149.6e+9;	//distance from earth to sun, in m 
	//real omega = 0; //no rotation
	//real omega = 2. * M_PI / (60. * 60. * 24. * 365.25) / c;	//one revolution per year in m^-1 
	//real omega = 1;	//angular velocity of the speed of light
	real const omega = c;	//I'm trying to find a difference ...
	real const t = 0;	//where the position should be.  t=0 means the body is moved by [L, 0], and its derivatives are along [0, L omega] 
	Vector<real,2> dt_xHat(L * omega * sin(omega * t), -L * omega * cos(omega * t));
	dt_metricPrims.alpha = dr_alpha * (xi(0)/r * dt_xHat(0) + xi(1)/r * dt_xHat(1));
	for (int i = 0; i < sDim; ++i) {
		dt_metricPrims.betaU(i) = 0;
	}
	for (int i = 0; i < sDim; ++i) {
		for (int j = 0; j < sDim; ++j) {
			real sum = 0;
			for (int k = 0; k < 2; ++k) {
				//gamma_ij = f/g
				//so d/dxHat^k gamma_ij = 
				real dxHat_k_of_gamma_ij = 
				// f' / g
				(
					((i==k)*xi(j) + xi(i)*(j==k)) * 2.*m + xi(i)*xi(j) * 2.*dr_m * xi(k)/r
				) / (r * r * (r - 2 * m))
				// - f g' / g^2
				- (xi(i) * xi(j) * 2 * m) * ( (xi(k) - 2 * dr_m * xi(k)) * r + 2 * xi(k) * (r - 2 * m) )
				/ (r * r * r * r * (r - 2 * m) * (r - 2 * m));
				sum += dxHat_k_of_gamma_ij * dt_xHat(k);
			}
			dt_metricPrims.gammaLL(i,j) = sum;
		}
	}
#endif

#if 0	//work that beta
	/*
	so if we can get the gamma^ij beta_j,t components to equal the Gamma^t_tt components ...
	 voila, gravity goes away.
	I'm approximating this as beta^i_,t ... but it really is beta^j_,t gamma_ij + beta^j gamma_ij,t
	 ... which is the same as what I've got, but I'm setting gamma_ij,t to zero
	*/
	//expanding ...
	MetricPrims& dt_metricPrims = dt_metricPrimGrid(index);
	TensorLsub betaL;
	for (int i = 0; i < sDim; ++i) {
		//negate all gravity by throttling the change in space/time coupling of the metric
		real dm_dr = 0;
		betaL(i) = -(2*m * (r - 2*m) + 2 * dm_dr * r * (2*m - r)) / (2 * r * r * r) * xi(i)/r;
	}
	TensorSUsub gammaUU = inverse(metricPrims.gammaLL);
	for (int i = 0; i < sDim; ++i) {
		real sum = 0;
		for (int j = 0; j < sDim; ++j) {
			sum += gammaUU(i,j) * betaL(j);
		}
		dt_metricPrims.betaU(i) = sum;
	}
#endif


]]},
	-- looks like an error
	{['stellar Kerr-Newman'] = [[
	real const radius = <?=solver.body.radius?>;
	real const mass = <?=solver.body.mass?>;
	real const density = <?=solver.body.density?>;
	
	real const angularVelocity = 2. * M_PI / (60. * 60. * 24.) / c;	//angular velocity, in m^-1
	real const inertia = 2. / 5. * mass * radius * radius;	//moment of inertia about a sphere, in m^3
	real const angularMomentum = inertia * angularVelocity;	//angular momentum in m^2
	real const a = angularMomentum / mass;	//m
			
	//real r is the solution of (x*x + y*y) / (r*r + a*a) + z*z / (r*r) = 1 
	// r^4 - (x^2 + y^2 + z^2 - a^2) r^2 - a^2 z^2 = 0
	real const RSq_minus_aSq = real3_lenSq(x) - a*a;
	//so we have two solutions ... which do we use? 
	//from gnuplot it looks like the two of these are the same ...
	r = sqrt((RSq_minus_aSq + sqrt(RSq_minus_aSq * RSq_minus_aSq + 4.*a*a*x.z*x.z)) / 2.);	//use the positive root

	//should I use the Kerr-Schild 'r' coordinate?
	//well, if 'm' is the mass enclosed within the coordinate
	// and that determines 'a', the angular momentum per mass within the coordinate (should it?)
	// then we would have a circular definition
	//real R = real3_len(x);
	real const matterRadius = min(r, radius);
	real const volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
	real const m = density * volumeOfMatterRadius;	// m^3

	real const Q = 0;	//charge
	real const H = (r*m - Q*Q/2.)/(r*r + a*a*x.z*x.z/(r*r));

	//3.4.33 through 3.4.35 of Alcubierre "Introduction to 3+1 Numerical Relativity"
	
	/*TODO fix this for the metric within the star
	 in other news, this is an unsolved problem!
	https://arxiv.org/pdf/1503.02172.pdf section 3.11
	https://arxiv.org/pdf/1410.2130.pdf section 4.2 last paragraph
	*/
	//metricPrims.alpha = 1./sqrt(1. + 2*H);
	gPrim->alpha = sqrt(1. - 2*H/(1+2*H) );
	
	real3 const l = _real3( (r*x.x + a*x.y)/(r*r + a*a), (r*x.y - a*x.x)/(r*r + a*a), x.z/r );
<?
for i=0,sDim-1 do
?>	gPrim->betaU.s<?=i?> = 2. * H * l.s<?=i?> / (1. + 2. * H);
<?
	for j=i,sDim-1 do
?>	gPrims->gammaLL.s<?=i..j?> = <?if i==j then?>1. + <?end?>2. * H * l.s<?=i?> * l.s<?=j?>;
<?	end
end
?>
]]},
}:map(function(kv)
	local k,v = next(kv)
	return {name=k, code=v}
end)

EFESolver.updateMethods = {'Newton', 'ConjRes', 'GMRes', 'JFNK'}

function EFESolver:init(args)
	self.app = args.app
	local config = args.config
	
	self.stDim = 4 		-- spacetime dim
	self.sDim = 3		-- space dim
	
	self.body = bodies[config.body]

	-- CLEnv:init calls CLEnv:getTypeCode()
	-- which (in EFESolver:getTypeCode) includes efe.h
	-- which depends on the body
	-- so if the body is changed.
	--  the ffi.cdef needs to be updated
	--  and self.code needs to be updated
	--  (and all subsequently influenced kernels / buffers need to be recompiled/allocated)
	EFESolver.super.init(self, {
		-- TODO rename to 'gridSize' ?
		size = {config.size, config.size, config.size},
		verbose = true,
	})

	-- parameters:

	self.xmin = vec3d(-1,-1,-1) * self.body.radius * config.bodyRadii
	self.xmax = vec3d(1,1,1) * self.body.radius * config.bodyRadii

	-- append efe.cl to the environment code
	self.code = self.code .. '\n' 
		.. self:compileTemplates(path'efe.cl':read())

	self.updateLambda = config.updateLambda

	self.initCond = self.initConds:find(nil, function(initCond)
		return initCond.name == config.initCond
	end) or 1

	-- what do we want to converge
	-- upon changing these, regenerate the gradientDescent.cl kernels
	self.convergeAlpha = true
	self.convergeBeta = false
	self.convergeGamma = false	-- TODO option for converging a scalar gamma vs a matrix gamma

	local function makeDiv(field)
		return template([[
	real div = 0.;
	<? for i=0,sDim-1 do ?>{
		int4 iL = i;
		iL.s<?=i?> = max(i.s<?=i?> - 1, 0);
		int const indexL = indexForInt4(iL);
		global <?=TPrim_t?> const * const TPrim_prev = TPrims + indexL;
		
		int4 iR = i;
		iR.s<?=i?> = min(i.s<?=i?> + 1, size.s<?=i?> - 1);
		int const indexR = indexForInt4(iR);
		global <?=TPrim_t?> const * const TPrim_next = TPrims + indexR;
		
		div += (TPrim_next-><?=field?>.s<?=i?> - TPrim_prev-><?=field?>.s<?=i?>) * .5 * inv_dx.s<?=i?>;
	}<? end ?>
	
	texCLBuf[index] = div;
]], {
	field = field,
	sDim = self.sDim,
	TPrim_t = self.TPrim_t,
})
	end

	-- this needs to be updated every time self.body changes
	-- once buffers are initialized, make displayVars
	--converts solver buffers to float[]
	self.displayVars = table()
	:append(
		table.map({
			{[{'TPrims'}] = table()
				:append(self.body.useMatter and {
					{['rho (g/cm^3)'] = 'texCLBuf[index] = TPrims[index].rho * c * c / G / 1000;'},
					{['P (m)'] = 'texCLBuf[index] = TPrims[index].P;'},
					{['eInt (m)'] = 'texCLBuf[index] = TPrims[index].eInt;'},
				} or nil)
				:append(self.body.useMatter and self.body.useVel and {
					{['v (m/s)'] = 'texCLBuf[index] = TPrims[index].v * c;'},
				} or nil)
				:append(self.body.useEM and {
					{['|E|'] = 'texCLBuf[index] = real3_len(TPrims[index].E);'},
					{['div E'] = makeDiv'E'},
					{['|B|'] = 'texCLBuf[index] = real3_len(TPrims[index].B);'},
					{['div B'] = makeDiv'B'},
				} or nil)
			},
			{[{'gPrims'}] = {
				{['alpha-1'] = 'texCLBuf[index] = gPrims[index].alpha - 1.;'},
				{['|beta|'] = 'texCLBuf[index] = real3_len(gPrims[index].betaU);'},
				{['det|gamma|-1'] = 'texCLBuf[index] = sym3_det(gPrims[index].gammaLL) - 1.;'},
			}},
			{[{'GammaULLs'}] = {
				{['numerical gravity'] = [[
	real3 const x = getX(i);
	real const r = real3_len(x);
	global tensor_4sym4 const * const GammaULL = GammaULLs + index;
	texCLBuf[index] = (0.
		+ GammaULL->s1.s00 * x.s0 / r
		+ GammaULL->s2.s00 * x.s1 / r
		+ GammaULL->s3.s00 * x.s2 / r) * c * c;
]]},
			}},	
			{[{}] = self.body.density and {
				{['analytical gravity'] = [[
	real3 const x = getX(i);
	real const r = real3_len(x);
	real const matterRadius = min(r, (real)<?=solver.body.radius?>);
	real const volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
	real const m = <?=solver.body.density?> * volumeOfMatterRadius;	// m^3
	real const dm_dr = 0;
	texCLBuf[index] = (2*m * (r - 2*m) + 2 * dm_dr * r * (2*m - r)) / (2 * r * r * r)
		* c * c;	//+9 at earth surface, without matter derivatives
]]},		
			} or {}},
			{[{'EFEs'}] = {
				{['EFE_tt (g/cm^3)'] = 'texCLBuf[index] = EFEs[index].s00 / (8. * M_PI) * c * c / G / 1000.;'},
				{['|EFE_ti|*c'] = [[
	global sym4 const * const EFE = EFEs + index;	
	texCLBuf[index] = sqrt(0.
<? for i=0,sDim-1 do ?>
		+ EFE->s0<?=i+1?> * EFE->s0<?=i+1?>
<? end ?>) * c;
]]},
				{['|EFE_ij|'] = [[
	global sym4 const * const EFE = EFEs + index;
	texCLBuf[index] = sqrt(0.
<? for i=0,sDim-1 do
	for j=0,sDim-1 do
	?>	+ EFE->s<?=sym(i+1,j+1)?> * EFE->s<?=sym(i+1,j+1)?>
<?	end
end ?>);
]]},
			}},
			{[{'gLLs', 'gUUs', 'GammaULLs'}] = {
				{['|Einstein_ab|'] = [[
	sym4 const EinsteinLL = calc_EinsteinLL(gLLs, gUUs, GammaULLs);
	texCLBuf[index] = sqrt(sym4_dot(EinsteinLL, EinsteinLL));
]]},
				{['Einstein_tt (g/cm^3)'] = [[
	sym4 const EinsteinLL = calc_EinsteinLL(gLLs, gUUs, GammaULLs);
	texCLBuf[index] = EinsteinLL.s00 / (8. * M_PI) * c * c / G / 1000.;
]]},
				{['|Einstein_ti|*c'] = [[
	sym4 const EinsteinLL = calc_EinsteinLL(gLLs, gUUs, GammaULLs);
	texCLBuf[index] = sqrt(0.
<? for i=0,sDim-1 do ?>
		+ EinsteinLL.s0<?=i+1?> * EinsteinLL.s0<?=i+1?>
<? end ?>) * c;
]]},		
				{['|Einstein_ij| (g/cm^3)'] = [[
	sym4 const EinsteinLL = calc_EinsteinLL(gLLs, gUUs, GammaULLs);
	texCLBuf[index] = sqrt(0.
<? for i=0,sDim-1 do
	for j=0,sDim-1 do
	?>	+ EinsteinLL.s<?=sym(i+1,j+1)?> * EinsteinLL.s<?=sym(i+1,j+1)?>
<?	end
end ?>) / (8. * M_PI) * c * c / G / 1000.;
]]},
		
			}},
		}, function(kv)
			local bufs, funcs = next(kv)
			return table.map(funcs, function(kv)
				local k,v = next(kv)
				return {name=k, argsIn=bufs, body=v}
			end)
		end):unpack()
	)
	self.displayVarNames = table.map(self.displayVars, function(displayVar) return displayVar.name end)

	self.updateMethod = table.find(self.updateMethods, config.solver) or 1
	self.useLineSearch = not not config.useLineSearch

	self:initBuffers()	
	self:refreshKernels()


	-- EFE: G_ab(g_ab) = 8 π T_ab(g_ab)
	-- consider x = α, β^i, γ_ij
	-- b = 8 π T_ab (and ignore the fact that it is based on x as well)
	-- linearize: G x = b
	
	local CLConjResSolver = require 'solver.cl.conjres':subclass()

	-- cache buffers
	function CLConjResSolver:newBuffer(name)
		self.conjResBuffers = self.conjResBuffers or {}
		if not self.conjResBuffers[name] then
			self.conjResBuffers[name] = CLConjResSolver.super.newBuffer(self, name)
		end
		return self.conjResBuffers[name]
	end

	self.conjResSolver = CLConjResSolver{
		env = self,
		A = function(y,x)
			-- treat 'x' as the gPrims
			-- change any kernels bound to gPrims to x instead
			self.calc_gLLs_and_gUUs.obj:setArg(2, x)
		
			-- TODO the EFE's don't need to be updated in this call
			-- they only need to be updated for ...
			-- * the Newton gradient descent solver
			-- * the display kernels
			-- so I say move calc_EFEs() to updateNewton, and code the EFE's directly into the display kernels that use them
			--  (because more often than not we're not watching the EFE constraints themselves)
			self:updateAux()
			
			-- bind our output to 'y'
			self.calc_EinsteinLLs.obj:setArg(0, y)
			self.calc_EinsteinLLs()
			
			-- fix the kernel arg state changes
			self.calc_gLLs_and_gUUs.obj:setArg(2, self.gPrims.obj)
		end,
		-- TODO if alpha, beta, or gamma are disabled then this can be a rectangular solver 
		x = self.gPrims,
		b = self._8piTLLs,
		type = 'real',
		size = self.base.volume * ffi.sizeof'gPrim_t' / ffi.sizeof'real',
		errorCallback = function(residual, iter)
			io.stderr:write('residual='..tostring(residual)..', iter='..tostring(iter)..'\n')
			if not math.isfinite(residual) then
				print("got a non-finite residual! "..tostring(residual))
				return
			end
		end,
		maxiter = self.base.volume * 10,
	}

	local CLGMResSolver = require 'solver.cl.gmres':subclass()

	-- cache buffers
	function CLGMResSolver:newBuffer(name)
		self.gmresBuffers = self.gmresBuffers or {}
		if not self.gmresBuffers[name] then
			self.gmresBuffers[name] = CLGMResSolver.super.newBuffer(self, name)
		end
		return self.gmresBuffers[name]
	end

	self.gmresSolver = CLGMResSolver{
		env = self,
		A = function(y,x)
			-- treat 'x' as the gPrims
			-- change any kernels bound to gPrims to x instead
			self.calc_gLLs_and_gUUs.obj:setArg(2, x)
		
			-- TODO the EFE's don't need to be updated in this call
			-- they only need to be updated for ...
			-- * the Newton gradient descent solver
			-- * the display kernels
			-- so I say move calc_EFEs() to updateNewton, and code the EFE's directly into the display kernels that use them
			--  (because more often than not we're not watching the EFE constraints themselves)
			self:updateAux()
			
			-- bind our output to 'y'
			self.calc_EinsteinLLs.obj:setArg(0, y)
			self.calc_EinsteinLLs()
			
			-- fix the kernel arg state changes
			self.calc_gLLs_and_gUUs.obj:setArg(2, self.gPrims.obj)
		end,
		-- TODO if alpha, beta, or gamma are disabled then this can be a rectangular solver 
		x = self.gPrims,
		b = self._8piTLLs,
		type = 'real',
		size = self.base.volume * ffi.sizeof'gPrim_t' / ffi.sizeof'real',
		errorCallback = function(residual, iter)
			io.stderr:write('residual='..tostring(residual)..', iter='..tostring(iter)..'\n')
			if not math.isfinite(residual) then
				print("got a non-finite residual! "..tostring(residual))
				return
			end
		end,
		maxiter = self.base.volume * 10,
	}
	-- I'm going to use the gmresSolver.dot
	-- for the norms of my Newton descent
	-- also TODO reuse the same dot with conjres, gmres, and jfnk
	
	local CLJFNKSolver = require 'solver.cl.jfnk':subclass()

	function CLJFNKSolver:newBuffer(name)
assert(type(name)=='string')	
		self.jfnkBuffers = self.jfnkBuffers or {}
		if not self.jfnkBuffers[name] then
			self.jfnkBuffers[name] = CLJFNKSolver.super.newBuffer(self, name)
		end
		return self.jfnkBuffers[name]
	end

	local gmresResidualFile -- = io.open('gmres_err.txt', 'w')
	local jfnkResidualFile -- = io.open('jfnk_err.txt', 'w')
	local jfnkIter
	self.jfnkSolver = CLJFNKSolver{
		env = self,
		f = function(y,x)
			
			-- solve for zero EFE_ab = G_ab - 8 π T_ab
			self.calc_gLLs_and_gUUs.obj:setArg(2, x)	-- input arg
			self.calc_EFEs.obj:setArg(0, y)
		
			-- calc_EFEs
			self:updateAux()
	
			-- fix the kernel arg state changes
			self.calc_gLLs_and_gUUs.obj:setArg(2, self.gPrims.obj)
			self.calc_EFEs.obj:setArg(0, self.EFEs.obj)
		end,
		x = self.gPrims,
		type = 'real',
		size = self.base.volume * ffi.sizeof'gPrim_t' / ffi.sizeof'real',
		errorCallback = function(residual, iter)
			--io.stderr:write('jfnk residual='..tostring(residual)..', iter='..tostring(iter)..'\n')
			if not math.isfinite(residual) then
				print("JFNK got a non-finite error! "..tostring(residual))
				return
			end
			jfnkIter = iter
			if jfnkResidualFile then
				jfnkResidualFile:write(iter,'\t',residual,'\n')
				jfnkResidualFile:flush()
			end
			if gmresResidualFile then
				gmresResidualFile:write'\n'
				gmresResidualFile:flush()
			end
		end,
		gmres = {
			MInv = function(y,x)
				-- the EFE is on the order of G/c^4 ~ 7e-47 so scale it up
				self.jfnkSolver.args.scale(y, x, c^4 / G)
				-- C++ version also scales by 1e-3
				-- and, for converging alpha only, scales by another c^2
			end,
			errorCallback = function(residual, iter)
				--io.stderr:write('gmres residual='..tostring(residual)..', iter='..tostring(iter)..'\n')
				if not math.isfinite(residual) then 
					print("GMRES got a non-finite error! "..tostring(residual))
					return
				end
				if gmresResidualFile then
					gmresResidualFile:write(jfnkIter,'\t',iter,'\t',residual,'\n')
					gmresResidualFile:flush()
				end
			end,
		},
		maxiter = 1,--self.base.volume * 10,
		jfnkEpsilon = 1e-6,
		epsilon = 1e-48,	-- efe error is G/c^4 ~ 6.67e-47
	}

	self:resetState()
end

--[[
this is called by CLEnv:init to get the initial set of cdef'd code
it is called right after self.real is set to float or double
therefore if self.real changes, a whole lot of other stuff will have to change too

Aside from that, in efesoln, some structs are populated based on the following:
	solver.body
	solver.useFourPotential
If either of these change, the structs have to change as well.
Luckily all those are in TPrim_t
so I'll just template out the name.
--]]
function EFESolver:getTypeCode()
	self.TPrim_t = 'TPrim_'..tostring(self.body.useMatter)
				..'_'..tostring(self.body.useVel)
				..'_'..tostring(self.body.useEM)
				..'_'..tostring(self.body.useFourPotential)
				..'_t'
	
	-- update this every time body changes
	return EFESolver.super.getTypeCode(self)..'\n'
	..template(path'efe.h':read(), {
		solver = self,
		TPrim_t = self.TPrim_t,
	})
end

function EFESolver:initBuffers()
	self.TPrims = self:buffer{name='TPrims', type=self.TPrim_t}
	self.gPrims = self:buffer{name='gPrims', type='gPrim_t'}
	self.gLLs = self:buffer{name='gLLs', type='sym4'}
	self.gUUs = self:buffer{name='gUUs', type='sym4'}
	self.GammaULLs = self:buffer{name='GammaULLs', type='tensor_4sym4'}
	
	-- used by updateNewton:
	self.EFEs = self:buffer{name='EFEs', type='sym4'}
	self.dPhi_dgPrims = self:buffer{name='dPhi_dgPrims', type='gPrim_t'}
	
	-- used by updateNewton's line trace
	 self.gPrimsCopy = self:buffer{name='gPrimsCopy', type='gPrim_t'}

	-- used by updateConjRes and updateGMRes:
	self._8piTLLs = self:buffer{name='_8piTLLs', type='sym4'}

	self.tex = require 'gl.tex3d'{
		width = tonumber(self.base.size.x),
		height = tonumber(self.base.size.y),
		depth = tonumber(self.base.size.z),
		internalFormat = gl.GL_RGBA32F,
		format = gl.GL_RGBA,
		type = gl.GL_FLOAT,
		minFilter = gl.GL_NEAREST,
		magFilter = gl.GL_LINEAR,
		wrap = {s=gl.GL_REPEAT, t=gl.GL_REPEAT, r=gl.GL_REPEAT},
	}

	-- TODO finishme
	self:clalloc(self.base.volume * ffi.sizeof'real', 'reduceBuf', 'real')
	self:clalloc(self.base.volume * ffi.sizeof'real' / self.base.localSize1d.x, 'reduceSwapBuf', 'real')
	self.reduceResultPtr = ffi.new('real[1]', 0)

	-- used for downloading visualization data
	self.texCLBuf = self:buffer{name='texCLBuf', type='float'} 

	if self.useGLSharing then
		self.texCLMem = require 'cl.imagegl'{context=self.ctx, tex=self.tex, write=true}
	else
		self.texCPUBuf = ffi.new('float[?]', self.base.volume)
	end
end

function EFESolver:compileTemplates(code)
	return template(code, {
		clnumber = clnumber,
		sym = sym,
		sDim = self.base.sDim,
		stDim = self.base.stDim,
		size = self.base.size,
		stDim = self.stDim,
		sDim = self.sDim,
		xmin = self.xmin,
		xmax = self.xmax,
		solver = self,
		c = c,
		G = G,
		TPrim_t = self.TPrim_t,
	})
end

function EFESolver:refreshKernels()
	-- create code
	print'preprocessing code...'

	local code = self:compileTemplates(table{
		path'calcVars.cl':read(),
		path'gradientDescent.cl':read(),
	}:concat'\n')

	-- keep all these kernels in one program.  what's the advantage?  less compiling I guess.
	local program = self:program{code=code} 

	-- init
	self.init_TPrims = program:kernel{name='init_TPrims', argsOut={self.TPrims}}
	-- compute values for EFE
	self.calc_gLLs_and_gUUs = program:kernel{name='calc_gLLs_and_gUUs', argsOut={self.gLLs, self.gUUs}, argsIn={self.gPrims}}
	self.calc_GammaULLs = program:kernel{name='calc_GammaULLs', argsOut={self.GammaULLs}, argsIn={self.gLLs, self.gUUs}}
	if self.useFourPotential then
		self.solveAL = program:kernel{name='solveAL', argsOut={self.TPrims}}
	end

	-- used by updateNewton and updateJFNK:
	self.calc_EFEs = program:kernel{name='calc_EFEs', argsOut={self.EFEs}, argsIn={self.TPrims, self.gLLs, self.gUUs, self.GammaULLs}}
	
	-- used by updateNewton:
	self.calc_dPhi_dgPrims = program:kernel{name='calc_dPhi_dgPrims', argsOut={self.dPhi_dgPrims}, argsIn={self.TPrims, self.gPrims, self.gLLs, self.gUUs, self.GammaULLs, self.EFEs}}
	self.update_gPrims = program:kernel{name='update_gPrims', argsOut={self.gPrims}, argsIn={self.dPhi_dgPrims}}

	-- used by updateConjRes and updateGMRes
	self.calc_EinsteinLLs = program:kernel{
		name = 'calc_EinsteinLLs',
		-- don't provide an actual buffer here
		-- the conjResSolver will provide its own
		argsOut = {{name='EinsteinLLs', type='sym4', obj=true}},
		argsIn = {self.gLLs, self.gUUs, self.GammaULLs},
	}
	self.calc_8piTLLs = program:kernel{name='calc_8piTLLs', argsOut={self._8piTLLs}, argsIn={self.TPrims, self.gLLs}}

	print'compiling code...'
	program:compile()
	
	print'done compiling code!'

	self:refreshInitCond()

	self.displayVar = 1
	self:refreshDisplayVarKernel()
end

function EFESolver:refreshInitCond()
	local initCond = self.initConds[self.initCond]
	self.init_gPrims = self:kernel{
		argsOut = {self.gPrims},
		body = [[
	real3 const x = getX(i);
	real const r = real3_len(x);

	global gPrim_t * const gPrim = gPrims + index;

	//init to flat by default
	*gPrim = (gPrim_t){
		.alpha = 1,
		.betaU = real3_zero,
		.gammaLL = sym3_ident,
	};

]]..self:compileTemplates(initCond.code),
	}
end

function EFESolver:refreshDisplayVarKernel()
	local displayVar = self.displayVars[self.displayVar]
	self.updateDisplayVarKernel = self:kernel(table(
		displayVar, {
			name = 'display_'..tostring(displayVar):sub(10),
			body = template(displayVar.body, {
				sDim = self.sDim,
				sym = sym,
				solver = self,
			}),
			argsIn = displayVar.argsIn and table.map(displayVar.argsIn, function(arg) 
				return self[arg]
			end) or nil,
			argsOut = {self.texCLBuf},
		}
	))
	self.updateDisplayVarKernel:compile()
	self:updateTex()
end

function EFESolver:resetState()
	self.init_gPrims()
	self.init_TPrims()

	self:updateAux()
	self:updateTex()

	self.iteration = 0
end

function EFESolver:update()
	local updateMethod = self.updateMethods[self.updateMethod]
	if updateMethod == 'Newton' then
		self:updateNewton()
	elseif updateMethod == 'ConjRes' then
		self:updateConjRes()
	elseif updateMethod == 'GMRes' then
		self:updateGMRes()
	elseif updateMethod == 'JFNK' then
		self:updateJFNK()
	end
	
	self.iteration = self.iteration + 1
end

function EFESolver:updateNewton()
	--[[
	iteration:

	∂g_ab/∂t = -∂Φ/∂g_ab
	for Φ = 1/2 Σ_ab (G_ab - 8 π T_ab)^2

	two approaches:
	1) do this per g_ab, so you only need to allocate as big as you would for solving the constraints themselves
	2) do this for g_ab as a whole, which would mean x4^2 symmetric = x10 allocation, but would take less kernel passes
	I'll try for 2 and hope I have enough memory
	--]]

	-- here's the newton update method
	self.calc_dPhi_dgPrims()

	-- now that we have ∂Φ/∂g_ab
	-- trace along g_ab - λ * ∂Φ/∂g_ab 
	-- to find what λ gives us minimal residual

	if self.useLineSearch then	-- do bisect line search
		-- store a backup.  TODO env:copy, and have ConjGrad:copy reference it.	
		self.conjResSolver.args.copy(self.gPrimsCopy, self.gPrims)

		local lineSearchMaxIter = 100
		local lambdaPtr = ffi.new'real[1]'
		local function residualAtLambda(lambda)
			self.conjResSolver.args.copy(self.gPrims, self.gPrimsCopy)
			lambdaPtr[0] = lambda	
			self.update_gPrims.obj:setArg(2, lambdaPtr)
			self.update_gPrims()
			self:updateAux()	-- calcs from gPrims on down to EFE
			local residual = self.conjResSolver.args.dot(self.EFEs, self.EFEs) 
print(string.format('lambda=%.16e residual=%.16e', lambda, residual))
			return residual
		end
		local function bisect(lambdaL, lambdaR)
			local residualL = residualAtLambda(lambdaL)
			local residualR = residualAtLambda(lambdaR)
			for i=0,lineSearchMaxIter do
				local lambdaMid = .5 * (lambdaL + lambdaR)
				local residualMid = residualAtLambda(lambdaMid)
				if residualMid > residualL and residualMid > residualR then break end
				if residualMid < residualL and residualMid < residualR then
					if residualL <= residualR then
						lambdaR, residualR  = lambdaMid, residualMid
					else
						lambdaL, residualL = lambdaMid, residualMid
					end
				elseif residualMid < residualL then
					lambdaL, residualL = lambdaMid, residualMid
				else
					lambdaR, residualR = lambdaMid, residualMid
				end
			end
			if residualL < residualR then
				return lambdaL, residualL
			else
				return lambdaR, residualR
			end
		end
		
		local lambdaFwd, residualFwd = bisect(0, self.updateLambda)
print(string.format('fwd lambda=%.16e residual=%.16e', lambdaFwd, residualFwd))
		local lambdaRev, residualRev = bisect(0, -self.updateLambda)	
print(string.format('rev lambda=%.16e residual=%.16e', lambdaRev, residualRev))

		self.conjResSolver.args.copy(self.gPrims, self.gPrimsCopy)
		local lambda, residual
		if residualFwd < residualRev then
			lambda = lambdaFwd
			residual = residualFwd
		else
			lambda = lambdaRev
			residual = residualRev
		end
print(string.format('using lambda=%.16e residual=%.16e', lambda, residualRev))
		lambdaPtr[0] = lambda
		self.update_gPrims.obj:setArg(2, lambdaPtr)
		self.update_gPrims()
	else	-- no line search
		-- update gPrims from dPhi/dg_ab 
		self.update_gPrims.obj:setArg(2, ffi.new('real[1]', self.updateLambda))
		self.update_gPrims()
	end

	--[[
	then there's the krylov solver treat-it-as-a-linear-system method
	G_ab(gPrims) = 8 π T_ab(also gPrims, but let's pretend not)
	A x = y
	solve using Jacobi method ... means isolating the diagonal terms
	solve using conjugate gradient / residual / bicgstab / gmres ...
	how about conj grad? ... in OpenCL ... 
	--]]

	self:updateAux()
print('residual', self.conjResSolver.args.dot(self.EFEs, self.EFEs))
	self:updateTex()
end

function EFESolver:updateConjRes()
	self.calc_8piTLLs()		-- update the b vector
	self.conjResSolver()	-- solve x in A x = b
end

function EFESolver:updateGMRes()
	self.calc_8piTLLs()		-- update the b vector
	self.gmresSolver()	-- solve x in A x = b
end

-- minimize alpha, beta, gamma
function EFESolver:updateJFNK()
	self.calc_EFEs()
	self.jfnkSolver()
end

function EFESolver:updateAux()
	-- every time gPrims changes, update these:
	self.calc_gLLs_and_gUUs()
	self.calc_GammaULLs()

	if self.useFourPotential then
		-- if we're using charge and current densities then here we'll have to recalculate A_u from J_u
		-- using Jacobi iteration
		for i=1,20 do
			self.solveAL()
		end
	end

	self.calc_EFEs()
end

function EFESolver:updateTex()
	self.updateDisplayVarKernel()
	
	-- TODO run a reduce on the display var stuff
	-- get the min and max
	-- then rescale the data according

-- notice soem of my results might not survive the double->float cast
-- because they exist in 1e-40 and what not

	-- now copy from cl buffer to gl buffer
	self.cmds[1]:enqueueReadBuffer{buffer=self.texCLBuf.obj, block=true, size=ffi.sizeof'float' * self.base.volume, ptr=self.texCPUBuf}

	local min, max = self.texCPUBuf[0], self.texCPUBuf[0]
	for i=1,self.base.volume-1 do
		local x = self.texCPUBuf[i]
		min = math.min(min, x) 
		max = math.max(max, x) 
	end
	self.app.minValue = min
	self.app.maxValue = max
	for i=1,self.base.volume-1 do
		self.texCPUBuf[i] = (self.texCPUBuf[i] - min) / (max - min)
	end

	self.tex:bind(0)
	for z=0,self.tex.depth-1 do
		gl.glTexSubImage3D(gl.GL_TEXTURE_3D, 0, 0, 0, z, self.tex.width, self.tex.height, 1, gl.GL_RED, gl.GL_FLOAT, self.texCPUBuf + self.tex.width * self.tex.height * z)
	end
	self.tex:unbind(0)
end

return EFESolver
