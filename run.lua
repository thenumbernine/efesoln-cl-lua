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

local c = 299792458			-- m/s 
local G = 6.67384e-11		-- m^3 / (kg s^2)

-- body parameters:

local body = {
	radius = 6.37101e+6,	-- m
	mass = 5.9736e+24 * G / c / c,	-- m
}
body.volume = 4/3 * math.pi * body.radius^3	-- m^3
body.density = body.mass / body.volume	-- 1/m^2

body.init = template([[
		
	//spherical body:
	
	TPrim->rho = r < <?=body.radius?> ? <?=body.density?> : 0;
	//TODO init this with the hydrostatic term of the schwarzschild equation of structure 

]], {body=body})

-- initial conditions:
-- [=[ flat 
local initCond = {code = ''}
--]=]
-- [=[ stellar schwarzschild 
local initCond = {
	code = 	template([[
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
}),
}
--]=]

-- end body parameters:

local bodyRadii = 2

local xmin = vec3d(-bodyRadii*body.radius,-bodyRadii*body.radius,-bodyRadii*body.radius)
local xmax = vec3d(bodyRadii*body.radius,bodyRadii*body.radius,bodyRadii*body.radius)

print'generating code...'

local typeCode = table{
	template([[
typedef <?=real?> real;
]], {real=real}),
	[[
typedef union {
	real s[3];
	struct { real s0, s1, s2; };
	struct { real x, y, z; };
} real3;

typedef union {
	real s[6];
	struct { real s00, s01, s02, s11, s12, s22; };	//useful for templated code
	struct { real xx, xy, xz, yy, yz, zz; };
} sym3;

typedef union {
	real s[10];
	struct { real s00, s01, s02, s03, s11, s12, s13, s22, s23, s33; };	//useful for templated code
	struct { real tt, tx, ty, tz, xx, xy, xz, yy, yz, zz; };
} sym4;

typedef union {
	sym4 s[3];
	struct { sym4 s0, s1, s2; };
} tensor_3sym4;


typedef union {
	sym4 s[4];
	struct { sym4 s0, s1, s2, s3; };
} tensor_4sym4;

typedef struct {
	real alpha;
	real3 betaU;
	sym3 gammaLL;
} gPrim_t;

typedef struct {
	//source terms:
	real rho;	//matter density
	real P;		//pressure ... due to matter.  TODO what about magnetic pressure?
	real eInt;	//specific internal energy

	real3 v;	//3-vel (upper, spatial)

//this needs to be lienar solved for ... but it's an easy problem (at least when geometry is flat)
//	real chargeDensity;
//	TensorUsub currentDensity;	//TODO how does this relate to matter density?

//in the mean time ...
	real3 E, B;	//upper, spatial

} TPrim_t;
]]
}:concat'\n'

-- luajit the types so I can see the sizeof (I hope OpenCL agrees with padding)

ffi.cdef(typeCode)

-- header includes macros as well as ffi-cdef code

local header = table{
	typeCode,
	[[

#define _real3(a,b,c) (real3){.s={a,b,c}}

static inline real real3_dot(real3 a, real3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline real real3_lenSq(real3 a) {
	return real3_dot(a,a);
}

static inline real real3_len(real3 a) {
	return sqrt(real3_lenSq(a));
}

static inline real3 real3_scale(real3 a, real s) {
	return _real3(a.x * s, a.y * s, a.z * s);
}

static inline real3 real3_add(real3 a, real3 b) {
	return _real3(a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline real3 real3_sub(real3 a, real3 b) {
	return _real3(a.x - b.x, a.y - b.y, a.z - b.z);
}

static inline real sym3_det(sym3 m) {
	return m.xx * m.yy * m.zz
		+ m.xy * m.yz * m.xz
		+ m.xz * m.xy * m.yz
		- m.xz * m.yy * m.xz
		- m.yz * m.yz * m.xx
		- m.zz * m.xy * m.xy;
}

static inline sym3 sym3_inv(real d, sym3 m) {
	return (sym3){
		.xx = (m.yy * m.zz - m.yz * m.yz) / d,
		.xy = (m.xz * m.yz - m.xy * m.zz) / d,
		.xz = (m.xy * m.yz - m.xz * m.yy) / d,
		.yy = (m.xx * m.zz - m.xz * m.xz) / d,
		.yz = (m.xz * m.xy - m.xx * m.yz) / d,
		.zz = (m.xx * m.yy - m.xy * m.xy) / d,
	};
}

static inline real3 sym3_real3_mul(sym3 m, real3 v) {
	return _real3(
		m.xx * v.x + m.xy * v.y + m.xz * v.z,
		m.xy * v.y + m.yy * v.y + m.yz * v.z,
		m.xz * v.z + m.yz * v.y + m.zz * v.z);
}
]],
	template([[
constant const int dim = <?=dim?>;	
constant const int subDim = <?=subDim?>;
constant const int gridDim = <?=gridDim?>;

constant const int4 size = (int4)(<?=clnumber(size.x)?>, <?=clnumber(size.y)?>, <?=clnumber(size.z)?>, 0);
constant const int4 stepsize = (int4)(1, <?=size.x?>, <?=size.x * size.y?>, <?=size.x * size.y * size.z?>);

#define INIT_KERNEL() \
	int4 i = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0); \
	if (i.x >= size.x || i.y >= size.y || i.z >= size.z) return; \
	int index = i.x + size.x * (i.y + size.y * i.z);

constant real3 xmin = _real3(<?=xmin.x?>, <?=xmin.y?>, <?=xmin.z?>);
constant real3 xmax = _real3(<?=xmax.x?>, <?=xmax.y?>, <?=xmax.z?>);
constant real3 dx = _real3(
	<?=tonumber(xmax.x - xmin.x) / tonumber(size.x)?>,
	<?=tonumber(xmax.x - xmin.x) / tonumber(size.x)?>,
	<?=tonumber(xmax.x - xmin.x) / tonumber(size.x)?>);

#define getX(i) _real3( \
	xmin.x + ((real)i.x + .5)/(real)size.x * (xmax.x - xmin.x),	\
	xmin.y + ((real)i.y + .5)/(real)size.y * (xmax.y - xmin.y),	\
	xmin.z + ((real)i.z + .5)/(real)size.z * (xmax.z - xmin.z));

]], {
	clnumber = clnumber,
	gridDim = gridDim,
	dim = dim,
	subDim = subDim,
	size = size,
	xmin = xmin,
	xmax = xmax,
}),
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

local TPrims = MetaBuffer{name='TPrims', type='TPrim_t'}
local gPrims = MetaBuffer{name='gPrims', type='gPrim_t'} 
local gLLs = MetaBuffer{name='gLLs', type='sym4'}
local gUUs = MetaBuffer{name='gUUs', type='sym4'}
local GammaULLs = MetaBuffer{name='GammaULLs', type='tensor_4sym4'}


local function clcall(kernel)
	cmds:enqueueNDRangeKernel{kernel=kernel, dim=gridDim, globalSize=size:ptr(), localSize=localSize:ptr()}
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
	if select('#', ...) then
		self.kernel:setArgs(...)
	end
	clcall(self.kernel)
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
		.gammaLL = (sym3){.xx=1, .yy=1, .zz=1, .xy=0, .xz=0, .yz=0},
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

]] .. body.init,
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
	tensor_3sym4 dgLLL3;
	tensor_4sym4 dgLLL;
	<? for i=0,gridDim-1 do ?>{
		global const sym4* gLL_prev = gLLs + index - stepsize.s<?=i?>;
		global const sym4* gLL_next = gLLs + index + stepsize.s<?=i?>;
		
		<? for a=0,dim-1 do ?>{
			<? for b=a,dim-1 do?>{
				real dg_iab = (gLL_next->s<?=a?><?=b?> - gLL_prev->s<?=a?><?=b?>) / (2. * dx.s<?=i?>);
				dgLLL3.s<?=i?>.s<?=a?><?=b?> = dg_iab;
				dgLLL.s<?=i+1?>.s<?=a?><?=b?> = dg_iab;
			}<? end ?>
		}<? end ?>
	}<? end ?>
	<? for a=0,dim-1 do ?>{
		<? for b=a,dim-1 do?>{
			dgLLL.s0.s<?=a?><?=b?> = 0;
		}<? end ?>
	}<? end ?>

	tensor_4sym4 GammaLLL;
	<? for a=0,dim-1 do ?>{
		<? for b=0,dim-1 do ?>{
			<? for c=b,dim-1 do ?>{
				GammaLLL.s<?=a?>.s<?=b?><?=c?> = .5 * (
					dgLLL.s<?=c?>.s<?=sym(a,b)?>
					+ dgLLL.s<?=b?>.s<?=sym(a,c)?>
					- dgLLL.s<?=a?>.s<?=sym(b,c)?>);
			}<? end ?>
		}<? end ?>
	}<? end ?>

	global const sym4* gUU = gUUs + index;
	tensor_4sym4 GammaULL;
	<? for a=0,dim-1 do ?>{
		<? for b=0,dim-1 do ?>{
			<? for c=b,dim-1 do ?>{
				real sum = 0;
				<? for d=0,dim-1 do ?>{
					sum += gUU->s<?=sym(a,d)?> * GammaLLL.s<?=d?>.s<?=b?><?=c?>;
				}<? end ?>
				GammaULL.s<?=a?>.s<?=b?><?=c?> = sum;
			}<? end ?>
		}<? end ?>
	}<? end ?>
]], {
		sym = sym,
		dim = dim,
		gridDim = gridDim,
	}),
}

local calc_EFE_constraint = MetaKernel{
	name = 'calc_EFE_constraint',
	
}

local metaKernels = table{
	init_gPrim,
	init_TPrim,
	calc_gLLs_and_gUUs,
	calc_GammaULLs,
}

-- create code

local code = table{header}
	:append(metaKernels:map(function(metaKernel) return metaKernel.code end))
	:concat'\n'

local program = require 'cl.program'{context=ctx, devices={device}, code=code}

-- get kernels

for _,metaKernel in ipairs(metaKernels) do
	metaKernel.kernel = program:kernel(metaKernel.name, metaKernel.argBuffers:unpack())
end

-- run the kernels

print'executing...'

init_gPrim()
init_TPrim()
calc_gLLs_and_gUUs()
calc_GammaULLs()
-- TODO solver here

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

print'done!'
