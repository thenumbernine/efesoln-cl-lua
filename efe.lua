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

-- generate the constraint error functions

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

-- constants

local c = 299792458			-- m/s 
local G = 6.67384e-11		-- m^3 / (kg s^2)

-- body parameters:

local body = {
	radius = 6.37101e+6,	-- m
	mass = 5.9736e+24 * G / c / c,	-- m
}
body.volume = 4/3 * math.pi * body.radius^3	-- m^3
body.density = body.mass / body.volume	-- 1/m^2

body.init = [[
		
	//spherical body:
	
	TPrim->rho = r < <?=body.radius?> ? <?=body.density?> : 0;
	//TODO init this with the hydrostatic term of the schwarzschild equation of structure 

]]

-- initial conditions:
-- [=[ flat 
local initCond = {code = ''}
--]=]
-- [=[ stellar schwarzschild 
local initCond = {
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
}
--]=]

-- end body parameters:

local bodyRadii = 2

local xmin = vec3d(-bodyRadii*body.radius,-bodyRadii*body.radius,-bodyRadii*body.radius)
local xmax = vec3d(bodyRadii*body.radius,bodyRadii*body.radius,bodyRadii*body.radius)

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

local function clalloc(size) 
	return ctx:buffer{rw=true, size=size} 
end

local MetaBuffer = class()
function MetaBuffer:init(args)
	self.name = args.name
	self.type = args.type
	self.buf = clalloc(volume * ffi.sizeof(self.type))
end
function MetaBuffer:toCPU()
	if not self.cpuMem then
		self.cpuMem = ffi.new(self.type..'[?]', volume)
	end
	cmds:enqueueReadBuffer{buffer=self.buf, block=true, size=ffi.sizeof(self.type) * volume, ptr=self.cpuMem}
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
end

local MetaKernel = class()

function MetaKernel:init(args)
	self.name = args.name
	self.argsOut = args.argsOut
	self.argsIn = args.argsIn
	self.argBuffers = table():append(self.argsOut):append(self.argsIn):map(function(arg) return arg.buf end)
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
	
	if select('#', ...) then
		self.kernel:setArgs(...)
	end
	cmds:enqueueNDRangeKernel{kernel=self.kernel, dim=gridDim, globalSize=size:ptr(), localSize=localSize:ptr()}
end

function MetaKernel:toCPU()
	for _,mb in ipairs(self.argsOut) do
		mb:toCPU()
	end
end

local init_gPrim = MetaKernel{
	name = 'init_gPrim',
	argsOut = {gPrims},
	code = [[
	real3 x = getX(i);
	real r = real3_len(x);

	global gPrim_t* gPrim = gPrims + index;

	//init to flat by default
	*gPrim = (gPrim_t){
		.alpha = 1,
		.betaU = _real3(0,0,0),
		.gammaLL = sym3_ident,
	};
]] .. initCond.code,
}

local init_TPrim = MetaKernel{
	name = 'init_TPrim',
	argsOut = {TPrims},
	code = [[
	real3 x = getX(i);
	real r = real3_len(x);
	
	global TPrim_t* TPrim = TPrims + index;
	
	*TPrim = (TPrim_t){
		.rho = 0,
		.eInt = 0,
		.P = 0,
		.v = _real3(0,0,0),
		.E = _real3(0,0,0),
		.B = _real3(0,0,0),
	};

]]..body.init,
}

local calc_gLLs_and_gUUs = MetaKernel{
	name = 'calc_gLLs_and_gUUs',
	argsOut = {gLLs, gUUs},
	argsIn = {gPrims},
	code = template([[
	global const gPrim_t* gPrim = gPrims + index;

	real alphaSq = gPrim->alpha * gPrim->alpha;
	real3 betaL = sym3_real3_mul(gPrim->gammaLL, gPrim->betaU);
	real betaSq = real3_dot(betaL, gPrim->betaU);

	global sym4* gLL = gLLs + index;
	gLL->s00 = -alphaSq;
	<? for i=0,subDim-1 do ?>{
		gLL->s0<?=i+1?> = betaL.s<?=i?>;
		<? for j=i,subDim-1 do ?>{
			gLL->s<?=i+1?><?=j+1?> = gPrim->gammaLL.s<?=i?><?=j?>;
		}<? end ?>
	}<? end ?>

	real det_gammaUU = sym3_det(gPrim->gammaLL);
	sym3 gammaUU = sym3_inv(det_gammaUU, gPrim->gammaLL);

	global sym4* gUU = gUUs + index;
	gUUs->s00 = -1./alphaSq;
	<? for i=0,subDim-1 do ?>{
		gUUs->s0<?=i+1?> = gPrim->betaU.s<?=i?> / alphaSq;
		<? for j=i,subDim-1 do ?>{
			gUUs->s<?=i+1?><?=j+1?> = gammaUU.s<?=i?><?=j?> - gPrim->betaU.s<?=i?> * gPrim->betaU.s<?=j?> / alphaSq;
		}<? end ?>
	}<? end ?>
]], {
		subDim = subDim,
	}),
}

local calc_GammaULLs = MetaKernel{
	name = 'calc_GammaULLs',
	argsOut = {GammaULLs},
	argsIn = {gLLs, gUUs},
	code = template([[
	
	//here's where the finite difference stuff comes in ...
	tensor_4sym4 dgLLL;
	dgLLL.s0 = sym4_zero;
	<? for i=0,gridDim-1 do ?>{
		int4 iL = i;
		iL.s<?=i?> = min(i.s<?=i?> + 1, size.s<?=i?> - 1);
		int indexL = indexForInt4(iL);
		global const sym4* gLL_prev = gLLs + indexL;
		
		int4 iR = i;
		iR.s<?=i?> = max(i.s<?=i?> - 1, 0);
		int indexR = indexForInt4(iR);
		global const sym4* gLL_next = gLLs + indexR;
		
		dgLLL.s<?=i+1?> = sym4_scale(
			sym4_sub(*gLL_next, *gLL_prev),
			1. / (2. * dx.s<?=i?> ));
	}<? end ?>

	tensor_4sym4 GammaLLL = (tensor_4sym4){
	<? for a=0,dim-1 do ?>
		<? for b=0,dim-1 do ?>
			<? for c=b,dim-1 do ?>
				.s<?=a?>.s<?=b?><?=c?> = .5 * (
					dgLLL.s<?=c?>.s<?=sym(a,b)?>
					+ dgLLL.s<?=b?>.s<?=sym(a,c)?>
					- dgLLL.s<?=a?>.s<?=sym(b,c)?>),
			<? end ?>
		<? end ?>
	<? end ?>};

	global const sym4* gUU = gUUs + index;
	global tensor_4sym4* GammaULL = GammaULLs + index;
	<? for a=0,dim-1 do ?>
		<? for b=0,dim-1 do ?>
			<? for c=b,dim-1 do ?>
				GammaULL->s<?=a?>.s<?=b?><?=c?> = 0.
				<? for d=0,dim-1 do ?>
					+ gUU->s<?=sym(a,d)?> * GammaLLL.s<?=d?>.s<?=b?><?=c?>
				<? end ?>;
			<? end ?>
		<? end ?>
	<? end ?>
]], {
		sym = sym,
		dim = dim,
		gridDim = gridDim,
	}),
}

local calc_EFE_constraint = MetaKernel{
	name = 'calc_EFE_constraint',
	argsOut = {EFEs},
	argsIn = {gPrims, TPrims, gLLs, gUUs, GammaULLs},
	code = template([[
	sym4 EinsteinLL = calc_EinsteinLL(gLLs, gUUs, GammaULLs);
	sym4 _8piTLL = calc_8piTLL(gLLs+index, TPrims+index);
	EFEs[index] = sym4_sub(EinsteinLL, _8piTLL);
]], {
	
}),
}

local metaKernels = table{
	init_gPrim,
	init_TPrim,
	calc_gLLs_and_gUUs,
	calc_GammaULLs,
	calc_EFE_constraint,
}

-- create code

compileMetaKernels(metaKernels)

-- run the kernels

print'executing...'

init_gPrim()
init_TPrim()
calc_gLLs_and_gUUs()
calc_GammaULLs()
calc_EFE_constraint()

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

print'calculating aux values...'
	
local detGammas = MetaBuffer{name='detGammas', type=real} 
MetaKernel{
	name = 'calc_detGammas',
	argsOut = {detGammas},
	argsIn = {gPrims},
	code = 'detGammas[index] = sym3_det(gPrims[index].gammaLL);',
}()
	
local numericalGravity = MetaBuffer{name='numericalGravity', type=real}
MetaKernel{
	name = 'calc_numericalGravity',
	argsOut = {numericalGravity},
	argsIn = {GammaULLs},
	code = [[
	real3 x = getX(i);
	real r = real3_len(x);
	global const tensor_4sym4* GammaULL = GammaULLs + index;
	numericalGravity[index] = (0.
		+ GammaULL->s1.s00 * x.x / r
		+ GammaULL->s2.s00 * x.y / r
		+ GammaULL->s3.s00 * x.z / r) * c * c;
]],
}()

local analyticalGravity = MetaBuffer{name='analyticalGravity', type=real}
MetaKernel{
	name = 'calc_analyticalGravity',
	argsOut = {analyticalGravity},
	code = [[
	real3 x = getX(i);
	real r = real3_len(x);
	real matterRadius = min(r, (real)<?=body.radius?>);
	real volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
	real m = <?=body.density?> * volumeOfMatterRadius;	// m^3
	real dm_dr = 0;
	analyticalGravity[index] = (2*m * (r - 2*m) + 2 * dm_dr * r * (2*m - r)) / (2 * r * r * r)
		* c * c;	//+9 at earth surface, without matter derivatives
]],
}()

local EinsteinLLs = MetaBuffer{name='EinsteinLLs', type='sym4'}
MetaKernel{
	name = 'calc_EinsteinLLs',
	argsIn = {gLLs, gUUs, GammaULLs},
	argsOut = {EinsteinLLs},
	code = [[
	EinsteinLLs[index] = calc_EinsteinLL(gLLs+index, gUUs+index, GammaULLs+index);
]],
}()

--I used to use a struct of its own, but LuaJIT kept returning nil when accessing it
-- so now I'm trying real4 ...
local EFE_and_Einstein_norms = MetaBuffer{
	name = 'EFE_and_Einstein_norms',
	type = 'real4',
}
MetaKernel{
	name = 'calc_EFE_and_Einstein_norms',
	argsOut = {EFE_and_Einstein_norms},
	argsIn = {EFEs, EinsteinLLs},
	code = [[
	global const sym4* EFE = EFEs + index;	
	global const sym4* EinsteinLL = EinsteinLLs + index;
	EFE_and_Einstein_norms[index] = (real4){
		//EFE_tt
		EFE->s00 / (8. * M_PI) * c * c / G / 1000.,
		
		//EFE_ti
		sqrt(0.
<? for i=0,subDim-1 do ?>
			+ EFE->s0<?=i+1?> * EFE->s0<?=i+1?>
<? end ?>) * c,
		
		//EFE_ij
		sqrt(0.
<? for i=0,subDim-1 do
	for j=0,subDim-1 do ?>
			+ EFE->s<?=sym(i+1,j+1)?> * EFE->s<?=sym(i+1,j+1)?>
<?	end
end ?>),
		
		//G_ab
		sqrt(0.
<? for a=0,dim-1 do
	for b=0,dim-1 do ?>
			+ EinsteinLL->s<?=sym(a,b)?> * EinsteinLL->s<?=sym(a,b)?>
<?	end
end ?>),
	};
]],
}()

print'copying to cpu...'
gPrims:toCPU()
TPrims:toCPU()
detGammas:toCPU()
numericalGravity:toCPU()
analyticalGravity:toCPU()
EinsteinLLs:toCPU()
EFE_and_Einstein_norms:toCPU()

print'outputting...'

local cols = {
	{ix = function(index,i,j,k) return i end},
	{iy = function(index,i,j,k) return j end},
	{iz = function(index,i,j,k) return k end},
	{rho = function(index) return TPrims.cpuMem[index].rho end},
	{['det-1'] = function(index) return detGammas.cpuMem[index]-1 end},
	{['alpha-1'] = function(index) return gPrims.cpuMem[index].alpha-1 end},
	{gravity = function(index) return numericalGravity.cpuMem[index] end},
	{analyticalGravity = function(index) return analyticalGravity.cpuMem[index] end},
	{['EFE_tt(g/cm^3)'] = function(index) return EFE_and_Einstein_norms.cpuMem[index].s0 end},
	{['EFE_ti(g/cm^3)'] = function(index) return EFE_and_Einstein_norms.cpuMem[index].s1 end},
	{['EFE_ij(g/cm^3)'] = function(index) return EFE_and_Einstein_norms.cpuMem[index].s2 end},
	{['G_ab'] = function(index) return EFE_and_Einstein_norms.cpuMem[index].s3 end},
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
