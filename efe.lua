#!/usr/bin/env luajit
local ffi = require 'ffi'
require 'ext'
local template = require 'template'
local vec3d = require 'ffi.vec.vec3d'

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
local localSize2d = {localSizeX, localSizeY}

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
print('localSize2d',localSize2d)
print('localSize3d',localSize)

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

-- body parameters:

local body = {
	radius = 6.37101e+6,	-- m
	mass = 5.9736e+24 * G / c / c,	-- m
}
body.volume = 4/3 * math.pi * body.radius * body.radius * body.radius	-- m^3
body.density = body.mass / body.volume	-- 1/m^2

body.init = template([[
		
	//spherical body:
	
	TPrim->rho = r < <?=body.radius?> ? <?=body.density?> : 0;
	//TODO init this with the hydrostatic term of the schwarzschild equation of structure 

]], {
	body = body,
})

-- initial conditions:
--[=[ flat 
local initCond = {code = ''}
--]=]
-- [=[ stellar schwarzschild 
local initCond = {
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

]], {
	subDim = subDim,
	body = body,
})
}
--]=]

-- end body parameters:

local bodyRadii = 2

local xmin = vec3d(-1,-1,-1) * body.radius * bodyRadii
local xmax = vec3d(1,1,1) * body.radius * bodyRadii

print'generating code...'

local typeCode = template(file['efe.h'], {
	real = real,
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

ffi.cdef(table{
	typeCode,
}:concat'\n')

-- header includes macros as well as ffi-cdef code

local headerCode = table{
	typeCode,
	file['efe.cl'],
}:concat'\n'

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
	})
end

function compileMetaKernels(mks)
	local code = compileTemplates(
		table{headerCode}
			:append(table.map(mks, function(mk) return mk.code end))
			:concat'\n'
	)

	local program = require 'cl.program'{context=ctx, devices={device}, code=code}

	for _,mk in ipairs(mks) do
		mk.program = program
		mk.kernel = program:kernel(mk.name, mk.argBuffers:unpack())
	end

	return program
end

local function clcall(kernel, ...)
	if select('#', ...) then
		kernel:setArgs(...)
	end
	cmds:enqueueNDRangeKernel{kernel=kernel, dim=gridDim, globalSize=size:ptr(), localSize=localSize:ptr()}
end

local MetaKernel = class()

function MetaKernel:init(args)
	self.name = args.name
	self.argsOut = args.argsOut
	self.argsIn = args.argsIn
	self.argBuffers = table()
		:append(self.argsOut)
		:append(self.argsIn)
		:map(function(arg) return arg.buf end)
	self.code = template([[
kernel void <?=self.name?>(
<?
local sep = ''
for _,arg in ipairs(self.argsOut or {}) do ?>
	<?=sep?>global <?=arg.type?>* <?=arg.name?>
<?
	sep = ', '
end
for _,arg in ipairs(self.argsIn or {}) do ?>
	<?=sep?>global const <?=arg.type?>* <?=arg.name?>
<?
	sep = ', '
end
?>
) {
	INIT_KERNEL();
	<?=args.code?>
}
]], {self=self, args=args})
end

function MetaKernel:__call(...)
	-- if we get a call request when we have no kernel/program, make sure to get one 
	if not self.kernel then
		compileMetaKernels{self}
	end

	clcall(self.kernel, ...)
end

-- create code

local program = compileMetaKernels({})
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
clcall(calc_gLLs_and_gUUs)
clcall(calc_GammaULLs)
clcall(calc_EFEs)

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

local alpha = 1
local update_dgLLs = MetaKernel{
	name = 'update_dgLLs',
	argsOut = {gLLs},
	argsIn = {dPhi_dgLLs},
	code = template([[
	gLLs[index] = sym4_sub(gLLs[index], sym4_scale(dPhi_dgLLs[index], -<?=alpha?>));
]], {alpha=alpha}),
}

local calc_dPhi_dgLL = MetaKernel{
	name = 'calc_dPhi_dgLL',
	argsIn = {EFEs, TPrims, gLLs, gUUs, GammaULLs},
	argsOut = {dPhi_dgLLs},
	code = [[
	
	global const sym4* gLL = gLLs + index;
	global const sym4* gUU = gUUs + index;
	global const tensor_4sym4* GammaULL = GammaULLs + index;

	//this is also in the Ricci computation, but should I be storing it?  is it too big?
	tensor_44sym4 dGammaLULL;
	dGammaLULL.s0 = tensor_4sym4_zero;
	<? for i=0,gridDim-1 do ?>{
		int4 iL = i;
		iL.s<?=i?> = max(i.s<?=i?> - 1, 0);
		int indexL = indexForInt4(iL);
		global const tensor_4sym4* GammaULL_prev = GammaULLs + indexL;
		
		int4 iR = i;
		iR.s<?=i?> = min(i.s<?=i?> + 1, size.s<?=i?> - 1);
		int indexR = indexForInt4(iR);
		global const tensor_4sym4* GammaULL_next = GammaULLs + indexR;
		
		dGammaLULL.s<?=i+1?> = tensor_4sym4_scale(
			tensor_4sym4_sub(*GammaULL_next, *GammaULL_prev),
			.5 * inv_dx.s<?=i?> );
	}<? end ?>

	tensor_44sym4 RiemannULLL = (tensor_44sym4){
<? for a=0,dim-1 do ?>
		.s<?=a?> = (tensor_4sym4){
	<? for b=0,dim-1 do ?>
			.s<?=b?> = (sym4){
		<? for c=0,dim-1 do ?>
			<? for d=c,dim-1 do ?>
				.s<?=c?><?=d?> = 
					dGammaLULL.s<?=c?>.s<?=a?>.s<?=sym(b,d)?>
					- dGammaLULL.s<?=d?>.s<?=a?>.s<?=sym(b,c)?> 
				<? for e=0,dim-1 do ?>
					+ GammaULL->s<?=a?>.s<?=sym(e,c)?> * GammaULL->s<?=e?>.s<?=sym(b,d)?>
					- GammaULL->s<?=a?>.s<?=sym(e,d)?> * GammaULL->s<?=e?>.s<?=sym(b,c)?>
				<? end ?>,
			<? end ?>
		<? end ?>
			},
	<? end ?>
		},
<? end ?>
	};

	real GammaUUL[4][4][4] = {
<? for a=0,dim-1 do ?>
	{<? for b=0,dim-1 do ?>
		{<? for c=0,dim-1 do ?>
			0.
			<? for d=0,dim-1 do ?>
			+ GammaULL->s<?=a?>.s<?=sym(d,c)?> * gUU->s-><?=sym(d,b)?>
			<? end ?>,
		<? end ?>},
	<? end ?>},
<? end ?>
	};

	tensor_44sym4 RiemannULLL = (tensor_44sym4){
	};

	//dR_ab/dg_pq
	tensor_sym4sym4 dRLL_dgLL = (tensor_sym4sym4){
<? for e=0,dim-1 do ?>
	<? for f=e,dim-1 do ?>
		.s<?=e?><?=f?> = (sym4){
		<? for a=0,dim-1 do ?>
			<? for b=a,dim-1 do ?>
				.s<?=a?><?=b?> = 0.
				<? for c=0,dim-1 do ?>
					+ GammaUUL[<?=e?>][<?=c?>][<?=c?>] * GammaULL.s<?=f?>.s<?=sym(b,a)?>
					- GammaUUL[<?=e?>][<?=c?>][<?=b?>] * GammaULL.s<?=f?>.s<?=sym(c,a)?>
					- gUU->s<?=sym(c,e)?> * RiemannULLL.s<?=f?>.s<?=a?>.s<?=sym(c,b)?>
				<? end ?>,
			<? end ?>
		<? end ?>
		},
	<? end ?>
<? end ?>
	};

	//dG_ab/dg_pq = tensor_sym4sym4.s[pq].s[ab]
	tensor_sym4sym4 dEinsteinLL_dgLL = (tensor_sym4sym4){
<? for e=0,dim-1 do ?>
	<? for f=e,dim-1 do ?>
		.s<?=e?><?=f?> = (sym4){
		<? for a=0,dim-1 do ?>
			<? for b=a,dim-1 do ?>
			.s<?=a?><?=b?> = dRLL_dgLL.s<?=e?><?=f?>.s<?=a?><?=b?>
				- .5 * (0.
				<? if e==a and f==b then ?> + Gaussian <? end ?>
				+ gLL->s<?=a?><?=b?> * (
					-RUU.s<?=e?><?=f?>
				<? for c=0,dim-1 do ?>
					<? for d=0,dim-1 do ?>
						gUU->s<?=sym(c,d)?> * dRLL_dgLL.s<?=e?><?=f?>.s<?=sym(c,d)?>
					<? end ?>
				<? end ?>
				)),
			<? end ?>
		<? end ?>},
	<? end ?>
<? end ?>
	};

	global const TPrim_t* TPrim = TPrims + index;

	real4 EU = (real4)(0 <?for i=0,2 do ?>, TPrim->E.s<?=i?> <? end ?>); 
	real4 EL = sym4_real4_mul(*gLL, EU);
	real ESq = dot(EL, EU);
	
	real4 BU = (real4)(0 <?for i=0,2 do ?>, TPrim->B.s<?=i?> <? end ?>); 
	real4 BL = sym4_real4_mul(*gLL, BU);
	real BSq = dot(BL, BU);

	real sqrt_det_g = sqrt(sym4_det(*gLL));
	real3 SL = real3_scale(real3_cross(TPrim->E, TPrim->B), sqrt_det_g);

	tensor_sym4sym4 d_8piTLL_dgLL = (tensor_sym4sym4){
<? for e=0,dim-1 do ?>
	<? for f=e,dim-1 do ?>
		.s<?=e?><?=f?> = (sym4){
			//dTEM_00/dg_<?=p?><?=q?>
			.s00 = <?=(e==0 or f==0) and '0' or 'TPrim.E.s'..(e-1) * TPrim.E.s'..(f-1)' ?>
					+ <?=(e==0 or f==0) and '0' or 'TPrim.B.s'..(e-1) * TPrim.B.s'..(f-1)' ?>,
		<? for i=0,subDim-1 do ?>
			.s0<?=i+1?> = -SL.s<?=i?> * gUU->s<?=e?><?=f?>,
			<? for j=i,subDim-1 do ?>
			.s<?=i+1?><?=j+1?> = 0.
				<? if e==i+1 and f==j+1 then ?> + ESq + BSq <? end ?>
				+ gLL->s<?=i?><?=j?> * (
					<?=(e==0 or f==0) and '0' or 'TPrim.E.s'..(e-1) * TPrim.E.s'..(f-1)' ?>
					+ <?=(e==0 or f==0) and '0' or 'TPrim.B.s'..(e-1) * TPrim.B.s'..(f-1)' ?>
				)
				- 2. * (0.
					<? if e==i+1 then ?>
					+ EU.s<?=f?> * EL.s<?=j+1?> + BU.s<?=f?> * BL.s<?=j+1?>
					<? end ?>
					<? if e==j+1 then ?>
					+ EU.s<?=f?> * EL.s<?=i+1?> + BU.s<?=f?> * BL.s<?=i+1?>
					<? end ?>
				),
			<? end ?>
		<? end ?>},
	<? end ?>
<? end ?>
	};

	global const sym4* EFE = EFEs + index;	// G_ab - 8 pi T_ab
	<? for p=0,dim-1 do ?>
		<? for q=p,dim-1 do ?>
	dPhi_dgLLs[index].s<?=p?><?=q?> = 0.
			<? for a=0,dim-1 do ?>
				<? for b=0,dim-1 do ?>
					+ EFE->s<?=sym(a,b)?> * (dEinsteinLL_dgLL.s<?=p?><?=q?>.s<?=sym(a,b)?> - d_8piTLL_dgLL.s<?=p?><?=q?>.s<?=sym(a,b)?>)
				<? end ?>
			<? end ?>;
		<? end ?>
	<? end ?>
]],
}


local maxiter = 0
for i=1,maxiter do
	calc_dPhi_dgLL()
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

-- too many allocations or something and luajit is giving me nils when accessing primitive arrays

print'calculating aux values...'

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
	{['EFE_tt(g/cm^3)'] = function(index) return norm_EFE_tts[index] end},
	{['|EFE_ti|(g/cm^3)'] = function(index) return norm_EFE_tis[index] end},
	{['|EFE_ij|(g/cm^3)'] = function(index) return norm_EFE_ijs[index] end},
	{['|G_ab|'] = function(index) return norm_EinsteinLLs[index] end},
}

local file = assert(io.open('out.txt', 'w'))
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
