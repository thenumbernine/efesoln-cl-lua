#!/usr/bin/env luajit
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local vec3d = require 'ffi.vec.vec3d'
local gl = require 'ffi.OpenGL'
local CLEnv = require 'clenv'

-- parameters:

local config = require 'config'

-- also in hydro-cl

local function clnumber(x)
	local s = tostring(tonumber(x))
	if s:find'e' then return s end
	if not s:find'%.' then s = s .. '.' end
	return s
end

-- end also in hydro-cl

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
	real3 x = getX(i);
	real r = real3_len(x);
	
	const real rho0 = <?=self.density?>;
	const real radius = <?=self.radius?>;
	const real radius3 = radius * radius * radius;
	const real mass = <?=self.mass?>;
	real r2 = r * r;
	TPrim->rho = r < radius ? rho0 : 0;
	TPrim->P = r < radius ? (rho0 * (
		(sqrt(1. - 2. * mass * r2 / radius3) - sqrt(1. - 2. * mass / radius))
		/ (3. * sqrt(1. - 2. * mass / radius) - sqrt(1. - 2. * mass * r2 / radius3))
	)) : 0;

]], {self=self})
end

local EMFieldBody = class()

function EMFieldBody:init(args)
	self.radius = args.radius
	
	self.useMatter = false
	self.useVel = false
	self.useEM = true
	
	self.init = template([[
	const real radius = <?=self.radius?>;
	
	real3 x = getX(i);
	
	real polar_rSq = x.x*x.x + x.y*x.y;
	real polar_r = sqrt(polar_rSq);		//r in polar coordinates	
	real dr = polar_r - radius;			//difference from polar radius to torus big radius
	real r = sqrt(x.z*x.z + dr*dr);		//r in torus radial coordinates
	real theta = atan2(x.z, dr);		//angle around the small radius
	real phi = atan2(x.x, x.y);			//angle around the big radius

	//F^uv_;v = -4 pi J^u
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

-- body parameters:

local bodies = {
	earth = SphericalBody{
		radius = 6.37101e+6,	-- m
		mass = 5.9736e+24 * G / c / c,	-- m
	},
	sun = SphericalBody{
		radius = 6.960e+8,	-- m
		mass = 1.9891e+30 * G / c / c,	-- m
	},
	['EM Field'] = EMFieldBody{
		radius = 2,
	},
}

-- initial conditions:

local EFESolver = class(CLEnv)

EFESolver.initConds = table{
	{flat = ''},
	{['stellar Schwarzschild'] = [[
	const real radius = <?=solver.body.radius?>;
	const real mass = <?=solver.body.mass?>;
	const real density = <?=solver.body.density?>;
	
	real matterRadius = (real)min(r, radius);
	real volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
	real m = density * volumeOfMatterRadius;	// m^3		
	
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
for i=0,subDim-1 do
?>	gPrim->betaU.s<?=i?> = 0;
<?	for j=i,subDim-1 do
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
	real dr_alpha = r > radius 
		? (mass / (r * r * sqrt(1. - 2. * mass / r))) 
		: (mass * r / (radius * radius * radius * sqrt(1. - 2. * mass * r * r / (radius * radius * radius))));
	real dr_m = r > radius ? 0 : (4. * M_PI * r * r * density);
	MetricPrims& dt_metricPrims = dt_metricPrimGrid(index);
	real L = 149.6e+9;	//distance from earth to sun, in m 
	//real omega = 0; //no rotation
	//real omega = 2. * M_PI / (60. * 60. * 24. * 365.25) / c;	//one revolution per year in m^-1 
	//real omega = 1;	//angular velocity of the speed of light
	real omega = c;	//I'm trying to find a difference ...
	real t = 0;	//where the position should be.  t=0 means the body is moved by [L, 0], and its derivatives are along [0, L omega] 
	Vector<real,2> dt_xHat(L * omega * sin(omega * t), -L * omega * cos(omega * t));
	dt_metricPrims.alpha = dr_alpha * (xi(0)/r * dt_xHat(0) + xi(1)/r * dt_xHat(1));
	for (int i = 0; i < subDim; ++i) {
		dt_metricPrims.betaU(i) = 0;
	}
	for (int i = 0; i < subDim; ++i) {
		for (int j = 0; j < subDim; ++j) {
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
	for (int i = 0; i < subDim; ++i) {
		//negate all gravity by throttling the change in space/time coupling of the metric
		real dm_dr = 0;
		betaL(i) = -(2*m * (r - 2*m) + 2 * dm_dr * r * (2*m - r)) / (2 * r * r * r) * xi(i)/r;
	}
	TensorSUsub gammaUU = inverse(metricPrims.gammaLL);
	for (int i = 0; i < subDim; ++i) {
		real sum = 0;
		for (int j = 0; j < subDim; ++j) {
			sum += gammaUU(i,j) * betaL(j);
		}
		dt_metricPrims.betaU(i) = sum;
	}
#endif


]]},
	{['stellar Kerr-Newman'] = [[
	const real radius = <?=solver.body.radius?>;
	const real mass = <?=solver.body.mass?>;
	
	real angularVelocity = 2. * M_PI / (60. * 60. * 24.) / c;	//angular velocity, in m^-1
	real inertia = 2. / 5. * mass * radius * radius;	//moment of inertia about a sphere, in m^3
	real angularMomentum = inertia * angularVelocity;	//angular momentum in m^2
	real a = angularMomentum / mass;	//m
			
	//real r is the solution of (x*x + y*y) / (r*r + a*a) + z*z / (r*r) = 1 
	// r^4 - (x^2 + y^2 + z^2 - a^2) r^2 - a^2 z^2 = 0
	real RSq_minus_aSq = x*x + y*y + z*z - a*a;
	//so we have two solutions ... which do we use? 
	//from gnuplot it looks like the two of these are the same ...
	real r = sqrt((RSq_minus_aSq + sqrt(RSq_minus_aSq * RSq_minus_aSq + 4.*a*a*z*z)) / 2.);	//use the positive root

	//should I use the Kerr-Schild 'r' coordinate?
	//well, if 'm' is the mass enclosed within the coordinate
	// and that determines 'a', the angular momentum per mass within the coordinate (should it?)
	// then we would have a circular definition
	//real R = sqrt(x*x + y*y + z*z); 
	real matterRadius = std::min<real>(r, radius);
	real volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
	real m = density * volumeOfMatterRadius;	// m^3

	real Q = 0;	//charge
	real H = (r*m - Q*Q/2.)/(r*r + a*a*z*z/(r*r));

	//3.4.33 through 3.4.35 of Alcubierre "Introduction to 3+1 Numerical Relativity"
	
	/*TODO fix this for the metric within the star
	 in other news, this is an unsolved problem!
	https://arxiv.org/pdf/1503.02172.pdf section 3.11
	https://arxiv.org/pdf/1410.2130.pdf section 4.2 last paragraph
	*/
	//metricPrims.alpha = 1./sqrt(1. + 2*H);
	metricPrims.alpha = sqrt(1. - 2*H/(1+2*H) );
	
	real3 l = _real3( (r*x + a*y)/(r*r + a*a), (r*y - a*x)/(r*r + a*a), z/r );
<?
for i=0,subDim-1 do
?>	gPrim->betaU.s<?=i?> = 2. * H * l.s<?=i?> / (1. + 2. * H);
<?
	for j=i,subDim-1 do
?>	gPrims->gammaLL.s<?=i..j?> = <?if i==j then?>1. + <?end?>2. * H * l.s<?=i?> * l.s<?=j?>;
<?	end
end
?>
]]},
}:map(function(kv)
	local k,v = next(kv)
	return {name=k, code=v}
end)

function EFESolver:init(args)
	self.app = args.app
	self.config = args.config
	
	CLEnv.init(self, {
		-- TODO rename to 'gridSize' ?
		size = {self.config.size, self.config.size, self.config.size},
		gridDim = 3,
	})
	
	self.dim = 4 		-- spacetime dim
	self.subDim = 3		-- space dim

	self.updateAlpha = ffi.new('float[1]', self.config.updateAlpha)
	
	self.body = bodies[self.config.body]

	self.initCondPtr = ffi.new('int[1]', 
		self.initConds:find(nil, function(initCond)
			return initCond.name == self.config.initCond
		end) or 1)

	-- what do we want to converge
	-- upon changing these, regenerate the gradientDescent.cl kernels
	self.convergeAlpha = ffi.new('bool[1]', true)
	self.convergeBeta = ffi.new('bool[1]', false)
	self.convergeGamma = ffi.new('bool[1]', false)	-- TODO option for converging a scalar gamma vs a matrix gamma

	-- parameters:

	self.xmin = vec3d(-1,-1,-1) * self.body.radius * self.config.bodyRadii
	self.xmax = vec3d(1,1,1) * self.body.radius * self.config.bodyRadii

	-- update this every time body changes
	self.typeCode = template(file['efe.h'], {
		solver = self,
	})

	-- luajit the types so I can see the sizeof (I hope OpenCL agrees with padding)

	ffi.cdef(self.typeCode)

	local function makeDiv(field)
		return template([[
	real div = 0.;
	<? for i=0,gridDim-1 do ?>{
		int4 iL = i;
		iL.s<?=i?> = max(i.s<?=i?> - 1, 0);
		int indexL = indexForInt4(iL);
		global const TPrim_t* TPrim_prev = TPrims + indexL;
		
		int4 iR = i;
		iR.s<?=i?> = min(i.s<?=i?> + 1, size.s<?=i?> - 1);
		int indexR = indexForInt4(iR);
		global const TPrim_t* TPrim_next = TPrims + indexR;
		
		div += (TPrim_next-><?=field?>.s<?=i?> - TPrim_prev-><?=field?>.s<?=i?>) * .5 * inv_dx.s<?=i?>;
	}<? end ?>
	
	texCLBuf[index] = div;
]], {
	field = field,
	gridDim = self.gridDim,
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
	real3 x = getX(i);
	real r = real3_len(x);
	global const tensor_4sym4* GammaULL = GammaULLs + index;
	texCLBuf[index] = (0.
		+ GammaULL->s1.s00 * x.s0 / r
		+ GammaULL->s2.s00 * x.s1 / r
		+ GammaULL->s3.s00 * x.s2 / r) * c * c;
]]},
			}},	
			{[{}] = self.body.density and {
				{['analytical gravity'] = [[
	real3 x = getX(i);
	real r = real3_len(x);
	real matterRadius = min(r, (real)<?=solver.body.radius?>);
	real volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
	real m = <?=solver.body.density?> * volumeOfMatterRadius;	// m^3
	real dm_dr = 0;
	texCLBuf[index] = (2*m * (r - 2*m) + 2 * dm_dr * r * (2*m - r)) / (2 * r * r * r)
		* c * c;	//+9 at earth surface, without matter derivatives
]]},		
			} or {}},
			{[{'EFEs'}] = {
				{['EFE_tt (g/cm^3)'] = 'texCLBuf[index] = EFEs[index].s00 / (8. * M_PI) * c * c / G / 1000.;'},
				{['|EFE_ti|'] = [[
	global const sym4* EFE = EFEs + index;	
	texCLBuf[index] = sqrt(0.
<? for i=0,subDim-1 do ?>
		+ EFE->s0<?=i+1?> * EFE->s0<?=i+1?>
<? end ?>) * c;
]]},
				{['|EFE_ij|'] = [[
	global const sym4* EFE = EFEs + index;
	texCLBuf[index] = sqrt(0.
	<? for i=0,subDim-1 do
	for j=0,subDim-1 do
	?>	+ EFE->s<?=sym(i+1,j+1)?> * EFE->s<?=sym(i+1,j+1)?>
	<?	end
	end ?>);
]]},
			}},
			{[{'gLLs', 'gUUs', 'GammaULLs'}] = {
				{['|Einstein_ab|'] = [[
	sym4 EinsteinLL = calc_EinsteinLL(gLLs, gUUs, GammaULLs);
	texCLBuf[index] = sqrt(sym4_dot(EinsteinLL, EinsteinLL));
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

	self:initBuffers()	
	self:refreshKernels()
	self:resetState()
end

function EFESolver:initBuffers()
	self.TPrims = self:makeBuffer{name='TPrims', type='TPrim_t'}
	self.gPrims = self:makeBuffer{name='gPrims', type='gPrim_t'} 
--	self.gPrimsCopy = self:makeBuffer{name='gPrimsCopy', type='gPrim_t'}
	self.gLLs = self:makeBuffer{name='gLLs', type='sym4'}
	self.gUUs = self:makeBuffer{name='gUUs', type='sym4'}
	self.GammaULLs = self:makeBuffer{name='GammaULLs', type='tensor_4sym4'}
	self.EFEs = self:makeBuffer{name='EFEs', type='sym4'}
	self.dPhi_dgPrims = self:makeBuffer{name='dPhi_dgPrims', type='gPrim_t'}
	
	self.tex = require 'gl.tex3d'{
		width = tonumber(self.size.x),
		height = tonumber(self.size.y),
		depth = tonumber(self.size.z),
		internalFormat = gl.GL_RGBA32F,
		format = gl.GL_RGBA,
		type = gl.GL_FLOAT,
		minFilter = gl.GL_NEAREST,
		magFilter = gl.GL_LINEAR,
		wrap = {s=gl.GL_REPEAT, t=gl.GL_REPEAT, r=gl.GL_REPEAT},
	}

	-- TODO finishme
	self:clalloc(self.volume * ffi.sizeof'real', 'reduceBuf', 'real')
	self:clalloc(self.volume * ffi.sizeof'real' / self.localSize1d, 'reduceSwapBuf', 'real')
	self.reduceResultPtr = ffi.new('real[1]', 0)

	-- used for downloading visualization data
	self.texCLBuf = self:makeBuffer{name='texCLBuf', type='float'} 

	if self.useGLSharing then
		self.texCLMem = require 'cl.imagegl'{context=self.ctx, tex=self.tex, write=true}
	else
		self.texCPUBuf = ffi.new('float[?]', self.volume)
	end
end

function EFESolver:compileTemplates(code)
	return template(code, {
		clnumber = clnumber,
		sym = sym,
		gridDim = self.gridDim,
		dim = self.dim,
		subDim = self.subDim,
		size = self.size,
		xmin = self.xmin,
		xmax = self.xmax,
		solver = self,
		c = c,
		G = G,
	})
end

function EFESolver:refreshKernels()
	-- create code
	print'preprocessing code...'

	local code = self:compileTemplates(table{
		self.typeCode,
		file['efe.cl'],
		file['calcVars.cl'],
		file['gradientDescent.cl'],
	}:concat'\n')

	-- keep all these kernels in one program.  what's the advantage?  less compiling I guess.
	local program = self:makeProgram{code=code} 

	-- init
	self.init_TPrims = program:kernel{name='init_TPrims', argsOut={self.TPrims}}
	-- compute values for EFE
	self.calc_gLLs_and_gUUs = program:kernel{name='calc_gLLs_and_gUUs', argsOut={self.gLLs, self.gUUs}, argsIn={self.gPrims}}
	self.calc_GammaULLs = program:kernel{name='calc_GammaULLs', argsOut={self.GammaULLs}, argsIn={self.gLLs, self.gUUs}}
	self.calc_EFEs = program:kernel{name='calc_EFEs', argsOut={self.EFEs}, argsIn={self.gPrims, self.TPrims, self.gLLs, self.gUUs, self.GammaULLs}}

	self.calc_dPhi_dgPrims = program:kernel{name='calc_dPhi_dgPrims', argsOut={self.dPhi_dgPrims}, argsIn={self.TPrims, self.gPrims, self.gLLs, self.gUUs, self.GammaULLs, self.EFEs}}
	
	self.update_gPrims = program:kernel{name='update_gPrims', argsOut={self.gPrims}, argsIn={self.dPhi_dgPrims}}

	print'compiling code...'
	program:compile()
	
	print'done compiling code!'

	self:refreshInitCond()

	self.displayVarPtr = ffi.new('int[1]', 0)
	self:refreshDisplayVarKernel()
end

function EFESolver:refreshInitCond()
	local initCond = self.initConds[self.initCondPtr[0]]
	self.init_gPrims = self:kernel{
		argsOut = {self.gPrims},
		header = self:compileTemplates(table{
			self.typeCode,
			file['efe.cl'],
		}:concat'\n'),
		body = [[
	real3 x = getX(i);
	real r = real3_len(x);

	global gPrim_t* gPrim = gPrims + index;

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
	local displayVar = self.displayVars[self.displayVarPtr[0]+1]
	self.updateDisplayVarKernel = self:kernel(table(
		displayVar, {
			name = 'display_'..tostring(displayVar):sub(10),
			header = self:compileTemplates(table{
				self.typeCode,
				file['efe.cl'],
			}:concat'\n'),
			body = template(displayVar.body, {
				subDim = self.subDim,
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

	-- every time gPrims changes, update these:
	self.calc_gLLs_and_gUUs()
	self.calc_GammaULLs()
	self.calc_EFEs()

	self:updateTex()

	self.iteration = 0
end

function EFESolver:update()
	--[[
	iteration:

	dg_ab/dt = -dPhi/dg_ab
	for Phi = 1/2 Sum_ab (G_ab - 8 pi T_ab)^2

	two approaches:
	1) do this per g_ab, so you only need to allocate as big as you would for solving the constraints themselves
	2) do this for g_ab as a whole, which would mean x4^2 symmetric = x10 allocation, but would take less kernel passes
	I'll try for 2 and hope I have enough memory
	--]]

	self.calc_dPhi_dgPrims()

	-- TODO now that we have dPhi/dg_ab
	-- trace along g_ab - alpha * dPhi/dg_ab
	-- to find what alpha gives us minimal error
--	self.cmds:enqueueCopyBuffer{src=self.gPrims.buf, dst=self.gPrimsCopy.buf, size=ffi.sizeof'gPrim_t' * self.volume}

	--[[ update g_ab and refresh gPrims from this
	self:clcall(update_gLLs)
	self:clcall(calc_gPrims_from_gLLs)
	--]]
	-- [[ or update gPrims directly from dPhi/dg_ab 
	self.update_gPrims.kernel:setArg(2, self.updateAlpha)
	self.update_gPrims()
	--]]

	-- update gPrim aux values 
	self.calc_gLLs_and_gUUs()
	self.calc_GammaULLs()
	self.calc_EFEs()
	
	-- and update the display buffer
	self:updateTex()
	
	self.iteration = self.iteration + 1
end

function EFESolver:updateTex()
	self.updateDisplayVarKernel()
	
	-- TODO run a reduce on the display var stuff
	-- get the min and max
	-- then rescale the data according

-- notice soem of my results might not survive the double->float cast
-- because they exist in 1e-40 and what not

	-- now copy from cl buffer to gl buffer
	self.cmds:enqueueReadBuffer{buffer=self.texCLBuf.buf, block=true, size=ffi.sizeof'float' * self.volume, ptr=self.texCPUBuf}

	local min, max = self.texCPUBuf[0], self.texCPUBuf[0]
	for i=1,self.volume-1 do
		local x = self.texCPUBuf[i]
		min = math.min(min, x) 
		max = math.max(max, x) 
	end
	self.app.minValue = min
	self.app.maxValue = max
	for i=1,self.volume-1 do
		self.texCPUBuf[i] = (self.texCPUBuf[i] - min) / (max - min)
	end

	self.tex:bind(0)
	for z=0,self.tex.depth-1 do
		gl.glTexSubImage3D(gl.GL_TEXTURE_3D, 0, 0, 0, z, self.tex.width, self.tex.height, 1, gl.GL_RED, gl.GL_FLOAT, self.texCPUBuf + self.tex.width * self.tex.height * z)
	end
	self.tex:unbind(0)
end

return EFESolver
