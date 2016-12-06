#!/usr/bin/env luajit
local ffi = require 'ffi'
require 'ext'
local template = require 'template'
local vec3d = require 'ffi.vec.vec3d'
local gl = require 'ffi.OpenGL'

-- parameters:

local config = require 'config'

-- also in hydro-cl

local vec3sz = require 'ffi.vec.create_ffi'(3,'size_t','sz')

local function clnumber(x)
	local s = tostring(tonumber(x))
	if s:find'e' then return s end
	if not s:find'%.' then s = s .. '.' end
	return s
end

-- end also in hydro-cl

local function get64bit(list)
	local best = list:map(function(item)
		local exts = item:getExtensions():lower():trim()
		return {item=item, fp64=exts:match'cl_%w+_fp64'}
	end):sort(function(a,b)
		return (a.fp64 and 1 or 0) > (b.fp64 and 1 or 0)
	end)[1]
	return best.item, best.fp64
end

-- helper for indexing symmetric matrices
local function sym(i,j)
	if i <= j then return i..j else return j..i end
end

--[[ generate the constraint error functions
--	this is going slow - requires some symmath optimizations (described in symmath/diffgeom.lua)

local function genCode()
	local symmath = require 'symmath'
	local Tensor = symmath.Tensor
	local var = symmath.var

	local t,x,y,z = symmath.vars('t', 'x', 'y', 'z')
	local coords = table{t,x,y,z}
	Tensor.coords{
		{variables=coords},
	}

	local gLLvars = table()
	local gUUvars = table()
	for a=1,4 do
		for b=a,4 do
			gLLvars[a..b] = var('g_{'..coords[a].name..coords[b].name..'}', coords)
			gUUvars[a..b] = var('g^{'..coords[a].name..coords[b].name..'}', coords)
		end
	end
	local gLL = Tensor('_ab', function(a,b)
		return gLLvars[sym(a,b)]
	end)
	local gUU = Tensor('^ab', function(a,b)
		return gUUvars[sym(a,b)]
	end)

	symmath.tostring = require 'symmath.tostring.MultiLine'
	local diffgeom = require 'symmath.diffgeom'(gLL, gUU)
	
	local MathJax = require 'symmath.tostring.MathJax'
	symmath.tostring = MathJax
	
	local output = table()
	diffgeom:print(function(s)
		output:insert(tostring(s)..'<br>\n')
	end)

	file['efe_gen.html'] = 
		MathJax.header
		..output:concat'\n'
		..MathJax.footer
end

genCode()
os.exit()
--]]

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


local EFESolver = class()

function EFESolver:init(args)
	self.app = args.app
	self.config = args.config
	
	self.platform = get64bit(require 'cl.platform'.getAll())
	self.device, self.fp64 = get64bit(self.platform:getDevices{gpu=true})
	
	local exts = string.split(string.trim(self.device:getExtensions()):lower(),'%s+')
	self.useGLSharing = exts:find(nil, function(ext) 
		return ext:match'cl_%w+_gl_sharing' 
	end)
self.useGLSharing = false 	-- for now
	self.ctx = require 'cl.context'{
		platform = self.platform, 
		device = self.device,
		glSharing = self.useGLSharing,
	}
	self.cmds = require 'cl.commandqueue'{context=self.ctx, device=self.device}

	-- initialize types

	self.real = self.fp64 and 'double' or 'float'

	self.dim = 4 		-- spacetime dim
	self.subDim = 3	-- space tim
	self.gridDim = 3	-- grid dim
	self.size = vec3sz(self.config.size, self.config.size, self.config.size)
	self.volume = tonumber(self.size:volume())

	
	-- https://stackoverflow.com/questions/15912668/ideal-global-local-work-group-sizes-opencl
	-- product of all local sizes must be <= max workgroup size
	local maxWorkGroupSize = tonumber(self.device:getInfo'CL_DEVICE_MAX_WORK_GROUP_SIZE')
	print('maxWorkGroupSize',maxWorkGroupSize)

	-- for volumes
	local localSize1d = math.min(maxWorkGroupSize, self.volume)

	-- for boundaries
	local localSizeX = math.min(tonumber(self.size.x), 2^math.ceil(math.log(maxWorkGroupSize,2)/2))
	local localSizeY = maxWorkGroupSize / localSizeX
	local localSize2d = table{localSizeX, localSizeY}

	--	localSize = gridDim < 3 and vec3sz(16,16,16) or vec3sz(4,4,4)
	-- TODO better than constraining by math.min(self.size),
	-- look at which sizes have the most room, and double them accordingly, until all of maxWorkGroupSize is taken up
	local localSize = vec3sz(1,1,1)
	local rest = maxWorkGroupSize
	local localSizeX = math.min(tonumber(self.size.x), 2^math.ceil(math.log(rest,2)/self.gridDim))
	localSize.x = localSizeX
	if self.gridDim > 1 then
		rest = rest / localSizeX
		if self.gridDim == 2 then
			localSize.y = math.min(tonumber(self.size.y), rest)
		elseif self.gridDim == 3 then
			local localSizeY = math.min(tonumber(self.size.y), 2^math.ceil(math.log(math.sqrt(rest),2)))
			localSize.y = localSizeY
			localSize.z = math.min(tonumber(self.size.z), rest / localSizeY)
		end
	end

	print('localSize1d',localSize1d)
	print('localSize2d',localSize2d:unpack())
	print('localSize3d',localSize:unpack())
	self.localSize = localSize

	
	self.body = bodies[self.config.body]

	self.initCond = initConds[self.config.initCond]
	self.initCond.code = template(self.initCond.code, {
		subDim = self.subDim,
		body = self.body,
	})

	-- parameters:

	self.xmin = vec3d(-1,-1,-1) * self.body.radius * self.config.bodyRadii
	self.xmax = vec3d(1,1,1) * self.body.radius * self.config.bodyRadii

	print'generating code...'

	self.typeCode = template(file['efe.h'], {
		solver = self,
	})

	-- luajit the types so I can see the sizeof (I hope OpenCL agrees with padding)

	-- boilerplate
	ffi.cdef(template([[
	typedef union {
		<?=real?> s[2];
		struct { <?=real?> s0, s1; };
		struct { <?=real?> x, y; };
	} <?=real?>2;

	//for real4 I'm using x,y,z,w to match OpenCL
	//...though for my own use I am storing t,x,y,z
	typedef union {
		<?=real?> s[4];
		struct { <?=real?> s0, s1, s2, s3; };
		struct { <?=real?> x, y, z, w; };	
	} <?=real?>4;


	]], {real=self.real}))

	ffi.cdef(self.typeCode)
	
	self.totalGPUSize = 0

	local solver = self 
	self.MetaBuffer = class()
	function self.MetaBuffer:init(args)
		self.name = args.name
		self.type = args.type
		self.buf = solver:clalloc(solver.volume * ffi.sizeof(self.type), name, self.type)
	end
	function self.MetaBuffer:toCPU()
		local cpuMem = ffi.new(self.type..'[?]', solver.volume)
		solver.cmds:enqueueReadBuffer{buffer=self.buf, block=true, size=ffi.sizeof(self.type) * self.volume, ptr=cpuMem}
		return cpuMem
	end

	self:initBuffers()
	self:initKernels()
	self:resetState()

end

function EFESolver:clalloc(size, name, ctype)
	self.totalGPUSize = self.totalGPUSize + tonumber(size)
	print((name and (name..' ') or '')..'allocating '..tonumber(size)..' bytes of type '..ctype..' with size '..ffi.sizeof(ctype)..', total '..self.totalGPUSize)
	return self.ctx:buffer{rw=true, size=size} 
end

function EFESolver:initBuffers()
	self.TPrims = self.MetaBuffer{name='TPrims', type='TPrim_t'}
	self.gPrims = self.MetaBuffer{name='gPrims', type='gPrim_t'} 
	self.gLLs = self.MetaBuffer{name='gLLs', type='sym4'}
	self.gUUs = self.MetaBuffer{name='gUUs', type='sym4'}
	self.GammaULLs = self.MetaBuffer{name='GammaULLs', type='tensor_4sym4'}
	self.EFEs = self.MetaBuffer{name='EFEs', type='sym4'}
	self.dPhi_dgLLs = self.MetaBuffer{name='dPhi_dgLL', type='sym4'}
	
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
	
	-- used for downloading visualization data
	self.texCLBuf =self.MetaBuffer{name='calcDisplayVarBuf', type='float'} 

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
		body = self.body,
		initCond = self.initCond,
		c = c,
		G = G,
		updateAlpha = self.config.updateAlpha,
	})
end

function EFESolver:clcall(kernel, ...)
	if select('#', ...) then
		kernel:setArgs(...)
	end
	self.cmds:enqueueNDRangeKernel{kernel=kernel, dim=self.gridDim, globalSize=self.size:ptr(), localSize=self.localSize:ptr()}
end

function EFESolver:initKernels()
	-- create code
	print'preprocessing code...'

	local code = self:compileTemplates(table{
		self.typeCode,
		file['efe.cl'],
		file['calcVars.cl'],
		file['gradientDescent.cl'],
		file['calcOutputVars.cl'],
	}:concat'\n')
	
	print'compiling code...'
	self.program = require 'cl.program'{context=self.ctx, devices={self.device}, code=code}
	
	print'done compiling code!'

	-- init
	self.init_gPrims = self.program:kernel('init_gPrims', self.gPrims.buf)
	self.init_TPrims = self.program:kernel('init_TPrims', self.TPrims.buf)
	-- compute values for EFE
	self.calc_gLLs_and_gUUs = self.program:kernel('calc_gLLs_and_gUUs', self.gLLs.buf, self.gUUs.buf, self.gPrims.buf)
	self.calc_GammaULLs = self.program:kernel('calc_GammaULLs', self.GammaULLs.buf, self.gLLs.buf, self.gUUs.buf)
	self.calc_EFEs = self.program:kernel('calc_EFEs', self.EFEs.buf, self.gPrims.buf, self.TPrims.buf, self.gLLs.buf, self.gUUs.buf, self.GammaULLs.buf)

	self.calc_dPhi_dgLLs = self.program:kernel('calc_dPhi_dgLLs', self.dPhi_dgLLs.buf, self.TPrims.buf, self.gLLs.buf, self.gUUs.buf, self.GammaULLs.buf, self.EFEs.buf)
	--local update_gLLs = self.program:kernel('update_gLLs', gLLs.buf, dPhi_dgLLs.buf)
	--local calc_gPrims_from_gLLs = self.program:kernel('calc_gPrims_from_gLLs', gPrims.buf, gLLs.buf)
	self.update_gPrims = self.program:kernel('update_gPrims', self.gPrims.buf, self.dPhi_dgLLs.buf)

	self.updateDisplayVarKernel = self.program:kernel('updateDisplayVar', self.texCLBuf.buf)
	
	-- alpha-1
	--self.updateDisplayVarKernel:setArg(1, self.gPrims.buf)
	-- numericalGravity
	self.updateDisplayVarKernel:setArg(1, self.GammaULLs.buf)
end

function EFESolver:resetState()
	self:clcall(self.init_gPrims)
	self:clcall(self.init_TPrims)

	-- every time gPrims changes, update these:
	self:clcall(self.calc_gLLs_and_gUUs)
	self:clcall(self.calc_GammaULLs)
	self:clcall(self.calc_EFEs)

	self:updateTex()
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

	self:clcall(self.calc_dPhi_dgLLs)

	--[[ update g_ab and refresh gPrims from this
	self:clcall(update_gLLs)
	self:clcall(calc_gPrims_from_gLLs)
	--]]
	-- [[ or update gPrims directly from dPhi/dg_ab 
	self:clcall(self.update_gPrims)
	--]]

	-- update gPrim aux values 
	self:clcall(self.calc_gLLs_and_gUUs)
	self:clcall(self.calc_GammaULLs)
	self:clcall(self.calc_EFEs)
	
	-- and update the display buffer
	self:updateTex()
end

function EFESolver:updateTex()
	self:clcall(self.updateDisplayVarKernel)
	
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
	print('alpha-1 min',min,'max',max)
	self.app.minValue = min
	self.app.maxValue = max
	for i=1,self.volume-1 do
		self.texCPUBuf[i] = (self.texCPUBuf[i] - min) / (max - min)
	end
	
	for z=0,self.tex.depth-1 do
		gl.glTexSubImage3D(gl.GL_TEXTURE_3D, 0, 0, 0, z, self.tex.width, self.tex.height, 1, gl.GL_RED, gl.GL_FLOAT, self.texCPUBuf + self.tex.width * self.tex.height * z)
	end
end

--[[ 
	then do some calculations
	EFE
	numericalGravity
	analyticalGravity

	then output it all
		i.x i.y i.z
		rho
		det(gamma_ij)-1
		alpha-1
		gravity
		analyticalGravity
		EFE_tt
		|EFE_ti|
		|EFE_ij|
		|G_ab|
--]]

function EFESolver:updateAuxBuffers()

	-- compute values for output
	self:clcall((program:kernel('calc_detGammas', tmp.buf, gPrims.buf)))
	local detGammas = tmp:toCPU()

	self:clcall((program:kernel('calc_numericalGravity', tmp.buf, GammaULLs.buf)))
	local numericalGravity = tmp:toCPU()

	self:clcall((program:kernel('calc_analyticalGravity', tmp.buf)))
	local analyticalGravity = tmp:toCPU()

	self:clcall((program:kernel('calc_norm_EFE_tts', tmp.buf, EFEs.buf)))
	local norm_EFE_tts = tmp:toCPU()

	self:clcall((program:kernel('calc_norm_EFE_tis', tmp.buf, EFEs.buf)))
	local norm_EFE_tis = tmp:toCPU()

	self:clcall((program:kernel('calc_norm_EFE_ijs', tmp.buf, EFEs.buf)))
	local norm_EFE_ijs = tmp:toCPU()

	self:clcall((program:kernel('calc_norm_EinsteinLLs', tmp.buf, gLLs.buf, gUUs.buf, GammaULLs.buf)))
	local norm_EinsteinLLs = tmp:toCPU()

	print'copying to cpu...'
	local gPrimsCPU = gPrims:toCPU()
	local TPrimsCPU = TPrims:toCPU()
	local dPhi_dgLLsCPU = dPhi_dgLLs:toCPU() 

	print'outputting...'

	local displayCols = {
		{ix = function(index,i,j,k) return i end},
		{iy = function(index,i,j,k) return j end},
		{iz = function(index,i,j,k) return k end},
		{['rho(g/cm^3)'] = function(index) return TPrimsCPU[index].rho * c * c / G / 1000 end},
		{['det-1'] = function(index) return -1+detGammas[index] end},
		{['alpha-1'] = function(index) return -1+gPrimsCPU[index].alpha end},
		{gravity = function(index) return numericalGravity[index] end},
		{analyticalGravity = function(index) return analyticalGravity[index] end},
	--[[ debugging
		{['EFE_tt(g/cm^3)'] = function(index) return norm_EFE_tts[index] end},
		{['|EFE_ti|(g/cm^3)'] = function(index) return norm_EFE_tis[index] end},
		{['|EFE_ij|(g/cm^3)'] = function(index) return norm_EFE_ijs[index] end},
		{['|G_ab|'] = function(index) return norm_EinsteinLLs[index] end},
	--]]
	-- [[ debugging the gradient descent
		{['dPhi/dg_tt'] = function(index) return dPhi_dgLLsCPU[index].s00 end}, 
		{['dPhi/dg_tx'] = function(index) return dPhi_dgLLsCPU[index].s01 end}, 
	--]]
	}

	local file = assert(io.open(config.outputFilename, 'w'))
	do
		file:write'#'
		local sep = ''
		for _,col in ipairs(displayCols) do
			file:write(sep..next(col))
			sep = '\t'
		end
	end
	file:write'\n'
	file:flush()

	local index = 0
	for i=0,tonumber(self.size.x-1) do
		for j=0,tonumber(self.size.y-1) do
			for k=0,tonumber(self.size.z-1) do
				local sep = ''
				for _,col in ipairs(displayCols) do
					local name, f = next(col)
					local x = f(index,i,j,k)
					--print(index,i,j,k,name,x)
					assert(x, "failed for col "..next(col))
					file:write(sep..('%.16e'):format(x))
					sep = '\t'
				end
				index = index + 1
				file:write'\n'
				file:flush()
			end
		end
	end
	file:close()

	print'done!'
end

return EFESolver
