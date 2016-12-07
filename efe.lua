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
	//spherical body:
	
	TPrim->rho = r < <?=self.radius?> ? <?=self.density?> : 0;
	//TODO init this with the hydrostatic term of the schwarzschild equation of structure 

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
}

-- initial conditions:

local initConds = {
	flat = {
		code = '',
	},
	stellar_schwarzschild = {
		code = [[
	real radius = <?=body.radius?>;
	real mass = <?=body.mass?>;

	real matterRadius = (real)min(r, radius);
	real volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
	real m = <?=body.density?> * volumeOfMatterRadius;	// m^3		

	gPrim->alpha = r > radius 
		? sqrt(1 - 2*mass/r)
		: (1.5 * sqrt(1 - 2*mass/radius) - .5 * sqrt(1 - 2*mass*r*r/(radius*radius*radius)));

	<? for i=0,subDim-1 do ?>
		gPrim->betaU.s<?=i?> = 0;
		<? for j=i,subDim-1 do ?>
			gPrim->gammaLL.s<?=i?><?=j?> = 
				x.s<?=i?>/r * x.s<?=j?>/r * 2*m/(r - 2*m)
				<? if i == j then ?> + 1. <? end ?>;
		<? end ?>
	<? end ?>

]],
	},
}


local EFESolver = class(CLEnv)

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

	self.initCond = initConds[self.config.initCond]
	self.initCond.code = template(self.initCond.code, {
		subDim = self.subDim,
		body = self.body,
	})

	-- what do we want to converge
	-- TODO make these checkboxes on the GUI
	-- except that means recompiling the gradient descent kernels if they're checked
	self.convergeAlpha = true
	self.convergeBeta = false
	self.convergeGamma = false	-- TODO option for converging a scalar gamma vs a matrix gamma

	-- parameters:

	self.xmin = vec3d(-1,-1,-1) * self.body.radius * self.config.bodyRadii
	self.xmax = vec3d(1,1,1) * self.body.radius * self.config.bodyRadii

	print'generating code...'

	self.typeCode = template(file['efe.h'], {
		solver = self,
	})

	-- luajit the types so I can see the sizeof (I hope OpenCL agrees with padding)

	ffi.cdef(self.typeCode)
	

	local solver = self 
	
	self:initBuffers()
	
	-- once buffers are initialized, make displayVars
	--converts solver buffers to float[]
	self.displayVars = {
		{
			name = 'rho (g/cm^3)',
			argsIn = {'TPrims'},
			body = 'texCLBuf[index] = TPrims[index].rho * c * c / G / 1000;',
		},
		{
			name = 'alpha-1',
			argsIn = {'gPrims'},
			body = [[
	texCLBuf[index] = gPrims[index].alpha - 1.;
]],
		},
		{
			name = '|beta|',
			argsIn = {'gPrims'},
			body = 'texCLBuf[index] = real3_len(gPrims[index].betaU);'
		},
		{
			name = 'det|gamma|-1',
			argsIn = {'gPrims'},
			body = 'texCLBuf[index] = sym3_det(gPrims[index].gammaLL) - 1.;',
		},
		{
			name = 'numerical gravity',
			argsIn = {'GammaULLs'},
			body = [[
	real3 x = getX(i);
	real r = real3_len(x);
	global const tensor_4sym4* GammaULL = GammaULLs + index;
	texCLBuf[index] = (0.
		+ GammaULL->s1.s00 * x.s0 / r
		+ GammaULL->s2.s00 * x.s1 / r
		+ GammaULL->s3.s00 * x.s2 / r) * c * c;
]],
		},
		{
			name = 'analytical gravity',
			body = [[
	real3 x = getX(i);
	real r = real3_len(x);
	real matterRadius = min(r, (real)<?=solver.body.radius?>);
	real volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
	real m = <?=solver.body.density?> * volumeOfMatterRadius;	// m^3
	real dm_dr = 0;
	texCLBuf[index] = (2*m * (r - 2*m) + 2 * dm_dr * r * (2*m - r)) / (2 * r * r * r)
		* c * c;	//+9 at earth surface, without matter derivatives
]],
		},
		{
			name = 'EFE_tt (g/cm^3)',
			argsIn = {'EFEs'},
			body = 'texCLBuf[index] = EFEs[index].s00 / (8. * M_PI) * c * c / G / 1000.;',
		},
		{
			name = '|EFE_ti|',
			argsIn = {'EFEs'},
			body = [[
	global const sym4* EFE = EFEs + index;	
	texCLBuf[index] = sqrt(0.
<? for i=0,subDim-1 do ?>
		+ EFE->s0<?=i+1?> * EFE->s0<?=i+1?>
<? end ?>) * c;
]],
		},
		{
			name = '|EFE_ij|',
			argsIn = {'EFEs'},
			body = [[
	global const sym4* EFE = EFEs + index;
	texCLBuf[index] = sqrt(0.
<? for i=0,subDim-1 do
	for j=0,subDim-1 do ?>
		+ EFE->s<?=sym(i+1,j+1)?> * EFE->s<?=sym(i+1,j+1)?>
<?	end
end ?>);
]],
		},
		{
			name = '|Einstein_ab|',
			argsIn = {'gLLs', 'gUUs', 'GammaULLs'},
			body = [[
	sym4 EinsteinLL = calc_EinsteinLL(gLLs, gUUs, GammaULLs);
	texCLBuf[index] = sqrt(sym4_dot(EinsteinLL, EinsteinLL));
]],
		},
	}
	self.displayVarNames = table.map(self.displayVars, function(displayVar) return displayVar.name end)
	
	self:initKernels()
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

function EFESolver:initKernels()
	-- create code
	print'preprocessing code...'

	local code = self:compileTemplates(table{
		self.typeCode,
		file['efe.cl'],
		file['calcVars.cl'],
		file['gradientDescent.cl'],
	}:concat'\n')
	
	print'compiling code...'
	self.program = require 'cl.program'{context=self.ctx, devices={self.device}, code=code}
	
	print'done compiling code!'

	-- keep all these kernels in one program.  what's the advantage?  less compiling I guess.
	local program = self:makeProgram{code=code} 

	-- init
	self.init_gPrims = program:kernel{name='init_gPrims', argsOut={self.gPrims}}
	self.init_TPrims = program:kernel{name='init_TPrims', argsOut={self.TPrims}}
	-- compute values for EFE
	self.calc_gLLs_and_gUUs = program:kernel{name='calc_gLLs_and_gUUs', argsOut={self.gLLs, self.gUUs}, argsIn={self.gPrims}}
	self.calc_GammaULLs = program:kernel{name='calc_GammaULLs', argsOut={self.GammaULLs}, argsIn={self.gLLs, self.gUUs}}
	self.calc_EFEs = program:kernel{name='calc_EFEs', argsOut={self.EFEs}, argsIn={self.gPrims, self.TPrims, self.gLLs, self.gUUs, self.GammaULLs}}

	self.calc_dPhi_dgPrims = program:kernel{name='calc_dPhi_dgPrims', argsOut={self.dPhi_dgPrims}, argsIn={self.TPrims, self.gPrims, self.gLLs, self.gUUs, self.GammaULLs, self.EFEs}}
	
	self.update_gPrims = program:kernel{name='update_gPrims', argsOut={self.gPrims}, argsIn={self.dPhi_dgPrims}}

	program:compile()

	self.displayVarPtr = ffi.new('int[1]', 0)
	self:refreshDisplayVarKernel()
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
