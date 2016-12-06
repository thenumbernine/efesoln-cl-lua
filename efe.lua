#!/usr/bin/env luajit
local ffi = require 'ffi'
require 'ext'
local template = require 'template'
local vec3d = require 'ffi.vec.vec3d'

-- parameters:

local config = {
	updateAlpha = 1,
	
	body = 'earth',
	
	bodyRadii = 2,
	
	initCond = 'flat',
	--initCond = 'stellar_schwarzschild',
	
	maxiter = 1,
	outputFilename = 'out.txt',
}

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

local platform = get64bit(require 'cl.platform'.getAll())
local device, fp64 = get64bit(platform:getDevices{gpu=true})

local ctx = require 'cl.context'{platform=platform, device=device}
local cmds = require 'cl.commandqueue'{context=ctx, device=device}

-- initialize types

local real = fp64 and 'double' or 'float'

local dim = 4 		-- spacetime dim
local subDim = 3	-- space tim
local gridDim = 3	-- grid dim
local size = vec3sz(64,64,64)
local volume = size:volume()

-- helper for indexing symmetric matrices
local function sym(i,j)
	if i <= j then return i..j else return j..i end
end

-- https://stackoverflow.com/questions/15912668/ideal-global-local-work-group-sizes-opencl
-- product of all local sizes must be <= max workgroup size
local maxWorkGroupSize = tonumber(device:getInfo'CL_DEVICE_MAX_WORK_GROUP_SIZE')
print('maxWorkGroupSize',maxWorkGroupSize)

-- for volumes
local offset = vec3sz(0,0,0)
local localSize1d = math.min(maxWorkGroupSize, tonumber(size:volume()))

-- for boundaries
local localSizeX = math.min(tonumber(size.x), 2^math.ceil(math.log(maxWorkGroupSize,2)/2))
local localSizeY = maxWorkGroupSize / localSizeX
local localSize2d = table{localSizeX, localSizeY}

--	localSize = gridDim < 3 and vec3sz(16,16,16) or vec3sz(4,4,4)
-- TODO better than constraining by math.min(size),
-- look at which sizes have the most room, and double them accordingly, until all of maxWorkGroupSize is taken up
local localSize = vec3sz(1,1,1)
local rest = maxWorkGroupSize
local localSizeX = math.min(tonumber(size.x), 2^math.ceil(math.log(rest,2)/gridDim))
localSize.x = localSizeX
if gridDim > 1 then
	rest = rest / localSizeX
	if gridDim == 2 then
		localSize.y = math.min(tonumber(size.y), rest)
	elseif gridDim == 3 then
		local localSizeY = math.min(tonumber(size.y), 2^math.ceil(math.log(math.sqrt(rest),2)))
		localSize.y = localSizeY
		localSize.z = math.min(tonumber(size.z), rest / localSizeY)
	end
end

print('localSize1d',localSize1d)
print('localSize2d',localSize2d:unpack())
print('localSize3d',localSize:unpack())

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

local body = bodies[config.body]

-- initial conditions:

local initConds = {
	flat = {code = ''},
	stellar_schwarzschild = {
		code = template([[
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

]], 	{
			subDim = subDim,
			body = body,
		})
	},
}

local initCond = initConds[config.initCond]

-- parameters:

local xmin = vec3d(-1,-1,-1) * body.radius * config.bodyRadii
local xmax = vec3d(1,1,1) * body.radius * config.bodyRadii

print'generating code...'

local typeCode = template(file['efe.h'], {
	real = real,
	body = body,
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


]], {real=real}))

ffi.cdef(typeCode)

-- allocate buffers

local totalGPUSize = 0
local function clalloc(size, name, ctype)
	totalGPUSize = totalGPUSize + tonumber(size)
	print((name and (name..' ') or '')..'allocating '..tonumber(size)..' bytes of type '..ctype..' with size '..ffi.sizeof(ctype)..', total '..totalGPUSize)
	return ctx:buffer{rw=true, size=size} 
end

local MetaBuffer = class()
function MetaBuffer:init(args)
	self.name = args.name
	self.type = args.type
	self.buf = clalloc(volume * ffi.sizeof(self.type), name, self.type)
end
function MetaBuffer:toCPU()
	local cpuMem = ffi.new(self.type..'[?]', volume)
	cmds:enqueueReadBuffer{buffer=self.buf, block=true, size=ffi.sizeof(self.type) * volume, ptr=cpuMem}
	return cpuMem
end

local TPrims = MetaBuffer{name='TPrims', type='TPrim_t'}
local gPrims = MetaBuffer{name='gPrims', type='gPrim_t'} 
local gLLs = MetaBuffer{name='gLLs', type='sym4'}
local gUUs = MetaBuffer{name='gUUs', type='sym4'}
local GammaULLs = MetaBuffer{name='GammaULLs', type='tensor_4sym4'}
local EFEs = MetaBuffer{name='EFEs', type='sym4'}

local function compileTemplates(code)
	return template(code, {
		clnumber = clnumber,
		sym = sym,
		gridDim = gridDim,
		dim = dim,
		subDim = subDim,
		size = size,
		xmin = xmin,
		xmax = xmax,
		body = body,
		initCond = initCond,
		c = c,
		G = G,
		updateAlpha = config.updateAlpha,
	})
end

local function clcall(kernel, ...)
	if select('#', ...) then
		kernel:setArgs(...)
	end
	cmds:enqueueNDRangeKernel{kernel=kernel, dim=gridDim, globalSize=size:ptr(), localSize=localSize:ptr()}
end

-- create code

local code = compileTemplates(table{
	typeCode,
	file['efe.cl'],
	file['calcVars.cl'],
}:concat'\n')
local program = require 'cl.program'{context=ctx, devices={device}, code=code}

-- init
local init_gPrims = program:kernel('init_gPrims', gPrims.buf)
local init_TPrims = program:kernel('init_TPrims', TPrims.buf)
-- compute values for EFE
local calc_gLLs_and_gUUs = program:kernel('calc_gLLs_and_gUUs', gLLs.buf, gUUs.buf, gPrims.buf)
local calc_GammaULLs = program:kernel('calc_GammaULLs', GammaULLs.buf, gLLs.buf, gUUs.buf)
local calc_EFEs = program:kernel('calc_EFEs', EFEs.buf, gPrims.buf, TPrims.buf, gLLs.buf, gUUs.buf, GammaULLs.buf)

-- run the kernels

print'executing...'

clcall(init_gPrims)
clcall(init_TPrims)
--[[
iteration:

dg_ab/dt = -dPhi/dg_ab
for Phi = 1/2 Sum_ab (G_ab - 8 pi T_ab)^2

two approaches:
1) do this per g_ab, so you only need to allocate as big as you would for solving the constraints themselves
2) do this for g_ab as a whole, which would mean x4^2 symmetric = x10 allocation, but would take less kernel passes
I'll try for 2 and hope I have enough memory
--]]

local dPhi_dgLLs = MetaBuffer{name='dPhi_dgLL', type='sym4'}
if config.maxiter > 0 then 
	print'compiling gradient descent code...'
	local code = compileTemplates(table{
		typeCode,
		file['efe.cl'],
		file['gradientDescent.cl'],
	}:concat'\n')
	local program = require 'cl.program'{context=ctx, devices={device}, code=code}
	
	local calc_dPhi_dgLLs = program:kernel('calc_dPhi_dgLLs', dPhi_dgLLs.buf, TPrims.buf, gLLs.buf, gUUs.buf, GammaULLs.buf, EFEs.buf)
	--local update_gLLs = program:kernel('update_gLLs', gLLs.buf, dPhi_dgLLs.buf)
	--local calc_gPrims_from_gLLs = program:kernel('calc_gPrims_from_gLLs', gPrims.buf, gLLs.buf)
	local update_gPrims = program:kernel('update_gPrims', gPrims.buf, dPhi_dgLLs.buf)

	print'executing gradient descent...'
	for i=1,config.maxiter do
		clcall(calc_gLLs_and_gUUs)
		clcall(calc_GammaULLs)
		clcall(calc_EFEs)

		clcall(calc_dPhi_dgLLs)
	
		--[[ update g_ab and refresh gPrims from this
		clcall(update_gLLs)
		clcall(calc_gPrims_from_gLLs)
		--]]
		-- [[ or update gPrims directly from dPhi/dg_ab 
		clcall(update_gPrims)
		--]]
	end

	print'...done!'
end

clcall(calc_gLLs_and_gUUs)
clcall(calc_GammaULLs)
clcall(calc_EFEs)

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

-- too many allocations or something and luajit is giving me nils when accessing primitive arrays

print'calculating aux values...'

local code = compileTemplates(table{
	typeCode,
	file['efe.cl'],
	file['calcOutputVars.cl'],
}:concat'\n')
local program = require 'cl.program'{context=ctx, devices={device}, code=code}

local tmp = MetaBuffer{name='tmp', type=real} 

-- compute values for output
clcall((program:kernel('calc_detGammas', tmp.buf, gPrims.buf)))
local detGammas = tmp:toCPU()

clcall((program:kernel('calc_numericalGravity', tmp.buf, GammaULLs.buf)))
local numericalGravity = tmp:toCPU()

clcall((program:kernel('calc_analyticalGravity', tmp.buf)))
local analyticalGravity = tmp:toCPU()

clcall((program:kernel('calc_norm_EFE_tts', tmp.buf, EFEs.buf)))
local norm_EFE_tts = tmp:toCPU()

clcall((program:kernel('calc_norm_EFE_tis', tmp.buf, EFEs.buf)))
local norm_EFE_tis = tmp:toCPU()

clcall((program:kernel('calc_norm_EFE_ijs', tmp.buf, EFEs.buf)))
local norm_EFE_ijs = tmp:toCPU()

clcall((program:kernel('calc_norm_EinsteinLLs', tmp.buf, gLLs.buf, gUUs.buf, GammaULLs.buf)))
local norm_EinsteinLLs = tmp:toCPU()

print'copying to cpu...'
local gPrimsCPU = gPrims:toCPU()
local TPrimsCPU = TPrims:toCPU()
local dPhi_dgLLsCPU = dPhi_dgLLs:toCPU() 

print'outputting...'

local cols = {
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
	for _,col in ipairs(cols) do
		file:write(sep..next(col))
		sep = '\t'
	end
end
file:write'\n'
file:flush()

local index = 0
for i=0,tonumber(size.x-1) do
	for j=0,tonumber(size.y-1) do
		for k=0,tonumber(size.z-1) do
			local sep = ''
			for _,col in ipairs(cols) do
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
