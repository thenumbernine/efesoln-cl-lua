#!/usr/bin/env luajit
local ffi = require 'ffi'
local class = require 'ext.class'
local string = require 'ext.string'	-- concat
local timer = require 'ext.timer'
local math = require 'ext.math'
local table = require 'ext.table'
local path = require 'ext.path'
local template = require 'template'
local struct = require 'struct'
local vec3d = require 'vec-ffi.vec3d'
local gl = require 'gl'
local CLEnv = require 'cl.obj.env'
local clnumber = require 'cl.obj.number'

local exec = require 'make.exec'
local makeTargets = require 'make.targets'

-- only write if the contents have changed, to not update the write date
local function writeIfChanged(filename, data)
	local srcpath = path(filename)
	local srcdata = srcpath:read()
	if srcdata ~= data then
--[[ if you want to diff what's been changed ...
		path'cache/tmp':write(data)
		exec(table{'diff', ('%q'):format(filename), 'cache/tmp'}:concat' ', false)
		path'cache/tmp':remove()
--]]
		assert(srcpath:write(data))
	end
end

-- TODO move this to cl?
-- and rename this? 
-- or just use the original CLProgram name?
local CLProgram = require 'cl.obj.program'
local CLClangProgram = CLProgram:subclass()

function CLClangProgram:init(args)
	if args.cacheFileName then
		self.cacheFileName = args.cacheFileName
		path'cache':mkdir()
		self.cacheFileCL = 'cache/'..self.cacheFileName..'.cl'
		self.cacheFileBC = 'cache/'..self.cacheFileName..'.bc'
		self.cacheFileSPV = 'cache/'..self.cacheFileName..'.spv'
	end

	-- ok I don't want super to hit the args.code condition
	-- but I want to store the code for later, so
	self.code = args.code
	args.code = nil
	
	CLClangProgram.super.init(self, args)
end

local function clangCompile(dst, src, buildOptions)
	exec(table{
		'clang',
		buildOptions or '',
		'-v',
		-- residual and numerical gravity giving nans ... 
		-- problem starts from GammaULLs getting nans.
		-- that comes from partial_xU_of_gLL getting nans
		-- that comes from return-struct-by-value producing corrupted values
		--'--target=spirv64-unknown-unknown',	
		--'--target=spirv',	--unsupported
		'--target=spir-unknown-unknown',
		'-emit-llvm',
		'-c',
		--'-O0',	-- -O0 makes some code break ... smh
		--'-O3',
		'-o', ('%q'):format(path(dst):fixpathsep()),
		('%q'):format(path(src):fixpathsep()),
	}:concat' ')
end

-- TODO would be nice to do this once and let the caller / subclass modify it
function CLClangProgram:setupBuildTargets(args)
	self.buildTargets = makeTargets{
		verbose = true,
		{
			srcs = {self.cacheFileCL},
			dsts = {self.cacheFileBC},
			rule = function(rule)
				assert(#rule.dsts == 1)
				assert(#rule.srcs == 1)
				clangCompile(rule.dsts[1], rule.srcs[1], args and arg.buildOption or nil)
			end,
		},
		{
			srcs = {self.cacheFileBC},
			dsts = {self.cacheFileSPV},
			rule = function()
				exec(table{
					'llvm-spirv',
					args and args.linkOptions or '',
					('%q'):format(self.cacheFileBC),
					'-o', ('%q'):format(self.cacheFileSPV),
				}:concat' ')
			end,
		},
	}
end

function CLClangProgram:build(args)
	local code = self:getCode()
	-- only write (and invalidate) when necessary
	writeIfChanged(self.cacheFileCL, code)
	self:setupBuildTargets(args)
	self.buildTargets:run(self.cacheFileBC)
end

function CLClangProgram:link(args)
	self:setupBuildTargets(args)
	self.buildTargets:run(self.cacheFileSPV)
	self.IL = assert(path(self.cacheFileSPV):read())
end

function CLClangProgram:compile(args)
	if self.cacheFileCL or self.cacheFileBC or self.cacheFileSPV then
		self:build(args)
		self:link(args)
		args = table(args):setmetatable(nil)
		args.verbose = true
	end
	local results = CLClangProgram.super.compile(self, args)
	assert(self.obj, "there must have been an error in your error handler")	-- otherwise it would have thrown an error
	do--if self.obj then	-- did compile
		print((self.cacheFileName and self.cacheFileName..' ' or '')..'log:')
		-- TODO log per device ...
		print(string.trim(self.obj:getLog(self.env.devices[1])))
	end
	return results
end

-- helper for indexing symmetric matrices
local function sym(i,j)
	if i <= j then return i..j else return j..i end
end

-- constants

local c = 299792458			-- m/s
local G = 6.67384e-11		-- m^3 / (kg s^2)

-- parameters:

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
	real3 const x = getX(i);						//[m]
	real const r = real3_len(x);					//[m]

	real const rho0 = <?=self.density?>;			//[1/m^2]
	real const radius = <?=self.radius?>;			//[m]
	real const radius3 = radius * radius * radius;	//[m^3]
	real const mass = <?=self.mass?>;				//[m]
	real const r2 = r * r;							//[m^2]
	TPrim->rho = r < radius ? rho0 : 0;				//[1/m^2]

	// equation of structure from 1973 MTW Gravitation, box 23.2 eqn 5...
	TPrim->P = r < radius ? (rho0 * (
		(sqrt(1. - 2. * mass * r2 / radius3) - sqrt(1. - 2. * mass / radius))
		/ (3. * sqrt(1. - 2. * mass / radius) - sqrt(1. - 2. * mass * r2 / radius3))
	)) : 0;											//[1/m^2]
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
		-- rho = 5515 kg/m^3 (correct avg overall is 5513.4 kg/m^3 ... density based on PREM is 13088 in center to 1020 at radius
		-- https://physics.stackexchange.com/questions/184032/what-is-the-pressure-at-the-center-of-the-earth-or-a-neutron-star
		-- ... P should reach 3.65e+11 Pa = kg / (m s^2) at the center of the Earth ...
		-- I'm just getting 1.7 ... which is half ...
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

EFESolver.Program = CLClangProgram

-- finite difference

-- source: https://en.wikipedia.org/wiki/Finite_difference_coefficient
-- derivCoeffs[derivative][order] = {coeffs...}
-- separating out the denom is improving my numerical accuracy.  without doing so, bssnok-fd-num, cartesian, minkowski, RK4 diverges.
local derivCoeffs = {
	-- centered, antisymmetric 1st deriv coefficients
	{
		[2] = {1, denom=2},
		[4] = {8, -1, denom=12},
		[6] = {45, -9, 1, denom=60},
		[8] = {672, -168, 32, -3, denom=840},
		[10] = {5/6, -5/21, 5/84, -5/504, 1/1260, denom=1},
	},
	-- centered, symmetric 2nd deriv coefficients
	{
		[2] = {[0] = -2, 1, denom=1},
		[4] = {[0] = -30, 16, -1, denom=12},
		[6] = {[0] = -490, 270, -27, 2, denom=180},
		[8] = {[0] = -205/72, 8/5, -1/5, 8/315, -1/560, denom=1},
	},
	-- centered, antisymmetric 3rd deriv coefficients
	{
		[2] = {-2, 1, denom=2},
		[4] = {-13, 8, -1, denom=8},
		[6] = {-488, 338, -72, 7, denom=240},
	},
	-- centered, symmetric 4th deriv coefficients
	{
		[2] = {[0] = 6, -4, 1, denom=1},
		[4] = {[0] = 56, -39, 12, -1, denom=6},
		--[6] = {[0] = 2730, -1952, 676, -96, 7},
	}
}

-- bake in denominator
local table = require 'ext.table'
for deriv, coeffsPerOrder in ipairs(derivCoeffs) do
	for order, coeffs in pairs(coeffsPerOrder) do
		local denom = coeffs.denom
		coeffs.denom = nil
		local keys = table.keys(coeffs)
		for _,i in ipairs(keys) do
			coeffs[i] = coeffs[i] / denom
		end
	end
end

-- 1st deriv
function EFESolver:finiteDifference(args)
	return template([[<?
local range = require 'ext.range'
local srcType = args.srcType
local getValue = args.getValue or function(args)
	return args.bufferName.."["..args.index.."]"
end
-- for if you want to getValue based on "i" instead of "index" ...
-- hmm this is ugly
local dstType = srcType == "" and "real4" or "real4x"..srcType
local resultName = args.resultName
local getBoundary = args.getBoundary or function(args)
	return "real"..srcType.."_zero"
end
-- derivCoeffs[1] have [0]=0, so they start at 1, and have implied antisymmetry (for d[-i], use -d[i])
local d1coeffs = assert(derivCoeffs[1][order])
?>	<?=dstType?> const <?=resultName?> = (<?=dstType?>){
		.s0 = real<?=srcType?>_zero,
<?
for i=0,sDim-1 do
?>		.s<?=i+1?> = real<?=srcType?>_add<?=#d1coeffs*2?>(
<?	for offset_i,coeff in ipairs(d1coeffs) do
?>			real<?=srcType?>_real_mul(
<?			-- setup rhs index
			args.i = "i + (int4)("
				..range(0,3):mapi(function(ii) return ii==i and offset_i or 0 end):concat', '
				..")"
			args.index = "index + stepsize.s"..i.." * "..offset_i
?>				(i.s<?=i?> + <?=offset_i?> >= size.s<?=i?>) ? <?=getBoundary(args)?> : <?=getValue(args)?>,
				<?=coeff?> * inv_dx.s<?=i?>
			),
			real<?=srcType?>_real_mul(
<?			-- setup lhs index
			args.i = "i - (int4)("
				..range(0,3):mapi(function(ii) return ii==i and offset_i or 0 end):concat', '
				..")"
			args.index = "index - stepsize.s"..i.." * "..offset_i
?>				(i.s<?=i?> - <?=offset_i?> < 0) ? <?=getBoundary(args)?> : <?=getValue(args)?>,
				<?=-coeff?> * inv_dx.s<?=i?>
			)<?=offset_i == #d1coeffs and "" or ","?>
<?	end
?>		),
<?
end
?>	};
]], {
		args = args,
		derivCoeffs = derivCoeffs,
		sDim = self.sDim,
		order = self.diffOrder,
	})
end

-- 2nd deriv
function EFESolver:finiteDifference2(args)
	return template([[<?
local range = require 'ext.range'
local srcType = args.srcType
local getValue = args.getValue or function(args)
	return args.bufferName.."["..args.index.."]"
end
local dstType = "real4s4x"..srcType
local resultName = args.resultName
local getBoundary = args.getBoundary or function(args)
	return "real"..srcType.."_zero"
end
local d1coeffs = assert(derivCoeffs[1][order])
local d2coeffs = assert(derivCoeffs[2][order], "couldn't find d2 coeffs for order "..order)
?>	<?=dstType?> <?=resultName?> = <?=dstType?>_zero;
	<? for i=0,sDim-1 do ?>{
		<? for j=i,sDim-1 do ?>{
<? if i == j then -- 2nd-deriv kernel
?>
			<? for k,coeff in pairs(d2coeffs) do
				args.i = "i + (int4)("
					..range(0,3):mapi(function(ii) return ii==i and -k or 0 end):concat', '
					..")"
				args.index = "index - stepsize.s"..i.." * "..k
			?>{
				real<?=srcType?> const yL = (i.s<?=i?> - <?=k?> < 0) ? <?=getBoundary(args)?> : <?=getValue(args)?>;
<?				args.i = "i + (int4)("
					..range(0,3):mapi(function(ii) return ii==i and k or 0 end):concat', '
					..")"
				args.index = "index + stepsize.s"..i.." * "..k
?>				real<?=srcType?> const yR = (i.s<?=i?> + <?=k?> >= size.s<?=i?>) ? <?=getBoundary(args)?> : <?=getValue(args)?>;
				<?=resultName?>.s<?=i+1?><?=j+1?> = real<?=srcType?>_add(
					<?=resultName?>.s<?=i+1?><?=j+1?>,
					real<?=srcType?>_real_mul(
						real<?=srcType?>_add(yR, yL),
						<?=coeff?> * inv_dx.s<?=i?> * inv_dx.s<?=j?>
					)
				);
			}<? end ?>

<? else	-- two 1st-deriv kernels
?>
			<? for k,coeff_k in pairs(d1coeffs) do ?>{
				<? for l,coeff_l in pairs(d1coeffs) do
					args.i = "i + (int4)("
						..range(0,3):mapi(function(ii) return ii==i and -k or (ii==j and -l or 0) end):concat', '
						..")"
					args.index = "index - stepsize.s"..i.." * "..k.." - stepsize.s"..j.." * "..l
				?>{
					real<?=srcType?> const yLL = (i.s<?=i?> - <?=k?> < 0 || i.s<?=j?> - <?=l?> < 0) ? <?=getBoundary(args)?> : <?=getValue(args)?>;
<?					
					args.i = "i + (int4)("
						..range(0,3):mapi(function(ii) return ii==i and -k or (ii==j and l or 0) end):concat', '
						..")"
					args.index = "index - stepsize.s"..i.." * "..k.." + stepsize.s"..j.." * "..l
?>					real<?=srcType?> const yLR = (i.s<?=i?> - <?=k?> < 0 || i.s<?=j?> + <?=l?> >= size.s<?=j?>) ? <?=getBoundary(args)?> : <?=getValue(args)?>;
<?					
					args.i = "i + (int4)("
						..range(0,3):mapi(function(ii) return ii==i and k or (ii==j and -l or 0) end):concat', '
						..")"
					args.index = "index + stepsize.s"..i.." * "..k.." - stepsize.s"..j.." * "..l
?>					real<?=srcType?> const yRL = (i.s<?=i?> + <?=k?> >= size.s<?=i?> || i.s<?=j?> - <?=l?> < 0) ? <?=getBoundary(args)?> : <?=getValue(args)?>;
<?					
					args.i = "i + (int4)("
						..range(0,3):mapi(function(ii) return ii==i and k or (ii==j and l or 0) end):concat', '
						..")"
					args.index = "index + stepsize.s"..i.." * "..k.." + stepsize.s"..j.." * "..l
?>					real<?=srcType?> const yRR = (i.s<?=i?> + <?=k?> >= size.s<?=i?> || i.s<?=j?> + <?=l?> >= size.s<?=j?>) ? <?=getBoundary(args)?> : <?=getValue(args)?>;
					
					<?=resultName?>.s<?=i+1?><?=j+1?> = real<?=srcType?>_add(
						<?=resultName?>.s<?=i+1?><?=j+1?>,
						real<?=srcType?>_real_mul(
							real<?=srcType?>_sub(
								real<?=srcType?>_add(yRR, yLL),
								real<?=srcType?>_add(yLR, yRL)
							),
							<?=coeff_k * coeff_l?> * inv_dx.s<?=i?> * inv_dx.s<?=j?>
						)
					);
				}<? end ?>
			}<? end ?>
<? end ?>
		}<? end ?>
	}<? end ?>
<?
]], {
		args = args,
		derivCoeffs = derivCoeffs,
		sDim = self.sDim,
		order = self.diffOrder,
	})
end



EFESolver.useFourPotential = false

EFESolver.initConds = table{
	{flat = 'calc_gPrim_flat'},
	{['stellar Schwarzschild'] = 'calc_gPrim_stellar_Schwarzschild'},
	-- looks like an error
	{['stellar Kerr-Newman'] = 'calc_gPrim_stellar_Kerr_Newman'},
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
	self.diffOrder = config.diffOrder	-- finite-difference order

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

	-- [[ I would put this in the type code, but it requires the type code to already be cdef'd
	-- that means it has to be inserted into the kernels' codes, or appended to the self.code
	do
		-- ffi.cdef already has stdint.h's uint8_t etc defined
		-- but OpenCL doesn't
		-- so ...
		self.code = self.code .. [[
typedef uchar uint8_t;
typedef char int8_t;
]]

		--[[
		-- if I use 'struct' then I only get one listing of either sij or tt..zz
		-- TODO change 'struct' to also have 'union', and make these as nested union+struct
		-- that'd mean one extra layer of metatypes but meh
		local real4s4x4s4_mt, real4s4x4s4_code = struct{
			unionField = 's',
			unionType = 'real4s4',
			fields = table.append(range(0,3):mapi(function(i)
				return range(i,3):mapi(function(j)
					return {'s'..i..j = 'real4s4'}
				end)
			end):unpack())
		}
		--]]

		self.gPrim_mt, self.gPrim_code = struct{
			name = 'gPrim_t',
			fields = {
				{alpha = 'real'},
				{betaU = 'real3'},
				{gammaLL = 'real3s3'},
			},
			-- real s[] as union access
			unionType = 'real',
			unionField = 's',
		}
--print(self.gPrim_code)
		self.code = self.code..'\n'..self.gPrim_code

		local TPrim_fields = table()
		--source terms:
		if self.body.useMatter then
			TPrim_fields:append{
				{rho = 'real'},		-- [1/m^2] matter density
				{P = 'real'},		-- [1/m^2] pressure ... due to matter.  TODO what about magnetic pressure?
				{eInt = 'real'},	-- [1] specific internal energy
			}
			if self.body.useVel then
				TPrim_fields:insert{v = 'real3'}	-- [1] 3-vel (upper, spatial)
			end
		end
		if self.body.useEM then
			if self.useFourPotential then
				TPrim_fields:append{
					{JL = 'real4'},	--4-current: rho, j  use to solve ...
					{AL = 'real4'},	--4-potential: phi, A
				}
			else
				--this needs to be lienar solved for ... but it's an easy problem (at least when geometry is flat)
				--TPrim_fields:append{
				--	{chargeDensity = 'real'},
				--	{currentDensity = 'TensorUsub'},	//TODO how does this relate to matter density?
				--}
				--in the mean time ...
				TPrim_fields:append{
					{E = 'real3'},
					{B = 'real3'},	--upper, spatial
				}
			end
		end

		self.TPrim_t = 'TPrim_'..tostring(self.body.useMatter)
			..'_'..tostring(self.body.useVel)
			..'_'..tostring(self.body.useEM)
			..'_'..tostring(self.useFourPotential)
			..'_t'

		self.TPrim_mt, self.TPrim_code = struct{
			name = self.TPrim_t,
			fields = TPrim_fields,
			-- real s[] as union access
			unionType = 'real',
			unionField = 's',
		}
--print(self.TPrim_code)
		self.code = self.code..'\n'..self.TPrim_code
	end
	--]]

	-- while we're here, create ffi.metatypes for all structs
	-- TODO use the struct-lua project for this for automatic string serialization
	-- TODO use vec-ffi for real3
	local ffi = require 'ffi'
	ffi.metatype('real3', {
		__tostring = function(x)
			return '{'..x.s0..', '..x.s1..', '..x.s2..'}'
		end,
		__concat = string.concat,
	})
	ffi.metatype('real3s3', {
		__tostring = function(x)
			return '{'
			..x.s00..', '..x.s01..', '..x.s02..', '
			..x.s11..', '..x.s12..', '
			..x.s22..'}'
		end,
		__concat = string.concat,
	})
	ffi.metatype('real4s4', {
		__tostring = function(x)
			return '{'
			..x.s00..', '..x.s01..', '..x.s02..', '..x.s03..', '
			..x.s11..', '..x.s12..', '..x.s13..', '
			..x.s22..', '..x.s23..', '
			..x.s33..'}'
		end,
		__concat = string.concat,
	})
	ffi.metatype('real4x4s4', {
		__tostring = function(x)
			return '{'..x.s0..', '..x.s1..', '..x.s2..', '..x.s3..'}'
		end,
		__concat = string.concat,
	})

	--[[
	self:checkStructSizes{
		'real3',
		'real3s3',
		'real4s4',	-- if its aligned between cl and cpu ... then why does the .sij and .s[ij] access produce different values?
		'real4x4s4',
		'real4x4x4s4',
		'real4s4x4s4',
		'gPrim_t',
		self.TPrim_t,
	}
	os.exit()
	--]]

	-- parameters:

	self.xmin = vec3d(-1,-1,-1) * self.body.radius * config.bodyRadii
	self.xmax = vec3d(1,1,1) * self.body.radius * config.bodyRadii

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
		<?=TPrim_t?> TPrim_prev;
		if (i.s<?=i?> > 0) {
			int4 iL = i;
			--iL.s<?=i?>;
			int const indexL = indexForInt4(iL);
			TPrim_prev = TPrims[indexL];
		} else {
			// boundary condition:
			TPrim_prev = TPrims[index];
		}

		<?=TPrim_t?> TPrim_next;
		if (i.s<?=i?> < size.s<?=i?> - 1) {
			int4 iR = i;
			++iR.s<?=i?>;
			int const indexR = indexForInt4(iR);
			TPrim_next = TPrims[indexR];
		} else {
			// boundary condition
			TPrim_next = TPrims[index];
		}

		div += (TPrim_next.<?=field?>.s<?=i?> - TPrim_prev.<?=field?>.s<?=i?>) * .5 * inv_dx.s<?=i?>;
	}<? end ?>

	texCLBuf[index] = div;
]], {
	field = field,
	sDim = self.sDim,
	TPrim_t = self.TPrim_t,
})
	end

	--[[
	List of predefined display vars
	Once buffers are initialized, make displayVars
	converts solver buffers to float[]
	--]]
	self.displayVars = table()
	:append{
		{['alpha-1'] = 'texCLBuf[index] = gPrims[index].alpha - 1.;'},
		{['|beta|'] = 'texCLBuf[index] = real3_len(gPrims[index].betaU);'},
		{['det|gamma|-1'] = 'texCLBuf[index] = real3s3_det(gPrims[index].gammaLL) - 1.;'},
	-- u'^i = -Γ^i_ab u^a u^b
	-- for weak-field, (u^i)^2 ≈ 0, u^t ≈ 1
	-- u'^i = -Γ^i_tt
	-- a^i = c^2 u'^i = -c^2 Γ^i_tt
	-- |a^i| = c^2 |Γ^i_tt|
		{['|Gamma^i_tt|'] = [[
texCLBuf[index] = real3_len(real4x4s4_i00(GammaULLs[index]));
]]},
		{['numerical gravity'] = [[
texCLBuf[index] = real3_len(real4x4s4_i00(GammaULLs[index])) * c * c;
]]},
		{['norm|EFE_ab|'] = 'texCLBuf[index] = real4s4_norm(EFEs[index]);'},
		{['EFE_tt (kg/m^3)'] = 'texCLBuf[index] = EFEs[index].s00 / (8. * M_PI) * c * c / G;'},
		{['|EFE_ti|*c'] = [[texCLBuf[index] = real3_len(real4s4_i0(EFEs[index])) * c;]]},
		{['det|EFE_ij| (kg/m s^2))'] = [[texCLBuf[index] = real3s3_det(real4s4_ij(EFEs[index])) / (8. * M_PI) * c * c * c * c / G;]]},
		{['norm|EFE_ij| (kg/m s^2))'] = [[texCLBuf[index] = real3s3_norm(real4s4_ij(EFEs[index])) / (8. * M_PI) * c * c * c * c / G;]]},
		{['|Einstein_ab|'] = [[
real4s4 const EinsteinLL = calc_EinsteinLL(i, gPrims, gLLs, gUUs, GammaULLs);
texCLBuf[index] = real4s4_norm(EinsteinLL);
]]},
		{['Einstein_tt (kg/m^3)'] = [[
real4s4 const EinsteinLL = calc_EinsteinLL(i, gPrims, gLLs, gUUs, GammaULLs);
texCLBuf[index] = EinsteinLL.s00 / (8. * M_PI) * c * c / G;
]]},
		{['|Einstein_ti|*c'] = [[
real4s4 const EinsteinLL = calc_EinsteinLL(i, gPrims, gLLs, gUUs, GammaULLs);
texCLBuf[index] = real3_len(real4s4_i0(EinsteinLL)) * c;
]]},
		{['det|Einstein_ij| (kg/(m s^2))'] = [[
real4s4 const EinsteinLL = calc_EinsteinLL(i, gPrims, gLLs, gUUs, GammaULLs);
texCLBuf[index] = real3s3_det(real4s4_ij(EinsteinLL)) / (8. * M_PI) * c * c * c * c / G;
]]},
		{['Gaussian'] = [[
real4s4 const RicciLL = calc_RicciLL(i, gPrims, gLLs, gUUs, GammaULLs);
texCLBuf[index] = real4s4_dot(RicciLL, gUUs[index]);
]]},
		{['partial_gLL_of_Phi'] = [[
real4s4 const partial_gLL_of_Phi = calc_partial_gLL_of_Phi(i, TPrims, gPrims, gLLs, gUUs, GammaULLs, EFEs);
texCLBuf[index] = real4s4_normSq(partial_gLL_of_Phi);
]]},
		{['partial_gPrim_of_Phi.alpha'] = [[
gPrim_t const partial_gPrim_of_Phi = calc_partial_gPrim_of_Phi(i, TPrims, gPrims, gLLs, gUUs, GammaULLs, EFEs);
texCLBuf[index] = partial_gPrim_of_Phi.alpha;
]]},
		{['norm|partial_gPrim_of_Phi.beta|'] = [[
gPrim_t const partial_gPrim_of_Phi = calc_partial_gPrim_of_Phi(i, TPrims, gPrims, gLLs, gUUs, GammaULLs, EFEs);
texCLBuf[index] = real3_norm(partial_gPrim_of_Phi.betaU);
]]},
		{['norm|partial_gPrim_of_Phi.gamma|'] = [[
gPrim_t const partial_gPrim_of_Phi = calc_partial_gPrim_of_Phi(i, TPrims, gPrims, gLLs, gUUs, GammaULLs, EFEs);
texCLBuf[index] = real3s3_norm(partial_gPrim_of_Phi.gammaLL);
]]},
--[=[
--[[
testing how well gradient-descent works for finite-difference stuff
Phi(x) = 1/2 |x|^2
dPhi(x)/dx^i = x^i
--]]
		{['test f(x)'] = [[
real3 const x = getX(i);
texCLBuf[index] = .5 * real3_lenSq(x);
]]},
		{['test |∇f(x)|'] = [[
real3 const x = getX(i);
texCLBuf[index] = real3_len(x);
]]},
		{['test |D[f(x)]|'] = [[
real3 const x = getX(i);
<?=solver:finiteDifference{
	getValue = function(args) return ".5 * real3_lenSq(getX("..args.i.."))" end,
	getBoundary = function(args) return ".5 * real3_lenSq(getX("..args.i.."))" end,
	srcType = "",
	resultName = "dPhi_dx",
}?>
texCLBuf[index] = real3_len(real4_to_real3(dPhi_dx));
]]},
		{['test |D[f(x)] - ∇f(x)|'] = [[
real3 const x = getX(i);
<?=solver:finiteDifference{
	getValue = function(args) return ".5 * real3_lenSq(getX("..args.i.."))" end,
	getBoundary = function(args) return ".5 * real3_lenSq(getX("..args.i.."))" end,
	srcType = "",
	resultName = "dPhi_dx",
}?>
texCLBuf[index] = real3_len(real3_sub(real4_to_real3(dPhi_dx), getX(i)));
]]},
--]=]
	}
	-- body-specific:
	:append(self.body.density and {
		{['analytical gravity'] = [[
real3 const x = getX(i);
real const r = real3_len(x);
real const matterRadius = min(r, (real)<?=solver.body.radius?>);
real const volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
real const m = <?=solver.body.density?> * volumeOfMatterRadius;	// m^3
real const dm_dr = 0;
texCLBuf[index] = (2*m * (r - 2*m) + 2 * dm_dr * r * (2*m - r)) / (2 * r * r * r) * c * c;	//+9 at earth surface, without matter derivatives
]]},
	} or nil)
	:append(self.body.useMatter and {
		{['rho (kg/m^3)'] = 'texCLBuf[index] = TPrims[index].rho * c * c / G;'},
		--[P] in 1/m^2
		--[P c^2 / G] in 1/m^2 * kg/m = kg/m^3
		--[P c^4 / G] in 1/m^2 * kg/m * m^2/s^2 = kg/(m s^2)
		{['P (kg/(m s^2))'] = 'texCLBuf[index] = TPrims[index].P * c * c * c * c / G;'},
		{['eInt (m^2/s^2)'] = 'texCLBuf[index] = TPrims[index].eInt * c * c;'},
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
	:mapi(function(kv)
		local k,v = next(kv)
		return {
			name = k,
			argsIn = bufs,
			body = v,
		}
	end)
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
	-- update this every time body changes
	local efe_h = template(path'efe.h':read(), {
		solver = self,
	})

--[=[ ok super defines real2 real4 real8
-- but turns out OpenCL's float2 float4 float8 don't have [] operator ...
-- and I need that so ...
	return EFESolver.super.getTypeCode(self)..'\n'
		..efe_h
--]=]
-- [=[ ... so just setup the 16/64 extensions and skip the real2/4/8
	return template([[
<? if real == 'double' then ?>
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
<? elseif real == 'half' then ?>
#pragma OPENCL EXTENSION cl_khr_fp16 : enable
<? end ?>
typedef <?=real?> real;
]],	{
		real = self.real,
	})..efe_h
--]=]
end

function EFESolver:initBuffers()
	self.TPrims = self:buffer{name='TPrims', type=self.TPrim_t}
	self.gPrims = self:buffer{name='gPrims', type='gPrim_t'}
	self.gLLs = self:buffer{name='gLLs', type='real4s4'}
	self.gUUs = self:buffer{name='gUUs', type='real4s4'}
	self.GammaULLs = self:buffer{name='GammaULLs', type='real4x4s4'}

	-- used by updateNewton:
	self.EFEs = self:buffer{name='EFEs', type='real4s4'}	-- 10 reals per size
	self.partial_gPrim_of_Phis = self:buffer{name='partial_gPrim_of_Phis', type='gPrim_t'}

	-- used by updateNewton's line trace
	 self.gPrimsCopy = self:buffer{name='gPrimsCopy', type='gPrim_t'}

	-- used by norms of EFEs
	 self.tmpBuf = self:buffer{name='tmpBuf', type='real4s4'}

	-- used by updateConjRes and updateGMRes:
	self._8piTLLs = self:buffer{name='_8piTLLs', type='real4s4'}

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

function EFESolver:template(code)
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
		derivCoeffs = derivCoeffs,
		c = c,
		G = G,
		TPrim_t = self.TPrim_t,
	})
end

function EFESolver:refreshKernels()
	-- create code
	print'preprocessing code...'

	-- append efe.cl to the environment code
	self.efeCode = self:template(table{
			path'efe.cl':read(),
			path'calcVars.cl':read(),
			path'gradientDescent.cl':read(),
		}:concat'\n')

	path'cache':mkdir()

	writeIfChanged('cache/efe.funcs.h', self:template(assert(path'efe.funcs.h':read())))
	writeIfChanged('cache/calcVars.h', self:template(assert(path'calcVars.h':read())))

	local efeProgram = self:program{
		cacheFileName = 'efe',	-- produces cache/efe.bc and cache/efe.spv
		code = self.efeCode,
	}

-- why does clCreateProgramWithIL take 45 seconds, while clCreateProgramWithSource takes 10 ... and clCreateProgramWithBinary takes no time at all.
timer('compiling code...', function()
	-- builds efe.cl, efe.bc, efe.spv
	efeProgram:compile()
print'efeProgram'
print'CL_PROGRAM_KERNEL_NAMES:'
print(require 'ext.tolua'(efeProgram.obj:getInfo'CL_PROGRAM_KERNEL_NAMES'))
end)

	-- init
	self.init_TPrims = efeProgram:kernel{
		name = 'init_TPrims',
		argsOut = {
			self.TPrims,
		},
	}

	-- compute values for EFE
	self.calc_gLLs_and_gUUs = efeProgram:kernel{
		name = 'calc_gLLs_and_gUUs',
		argsOut = {
			self.gLLs,
			self.gUUs,
		},
		argsIn = {
			self.gPrims,
		},
	}

	self.calc_GammaULLs = efeProgram:kernel{
		name = 'calc_GammaULLs',
		argsOut = {
			self.GammaULLs,
		},
		argsIn = {
			self.gPrims,
			self.gLLs,
			self.gUUs,
		},
	}

	if self.useFourPotential then
		self.solveAL = efeProgram:kernel{
			name = 'solveAL',
			argsOut = {
				self.TPrims,
			},
		}
	end

	-- used by updateNewton and updateJFNK:
	self.calc_EFEs = efeProgram:kernel{
		name = 'calc_EFEs',
		argsOut = {
			self.EFEs,
		},
		argsIn = {
			self.TPrims,
			self.gPrims,
			self.gLLs,
			self.gUUs,
			self.GammaULLs,
		},
	}

	-- used by updateNewton:
	self.calc_partial_gPrim_of_Phis_kernel = efeProgram:kernel{
		name = 'calc_partial_gPrim_of_Phis_kernel',
		argsOut = {
			self.partial_gPrim_of_Phis,
		},
		argsIn = {
			self.TPrims,
			self.gPrims,
			self.gLLs,
			self.gUUs,
			self.GammaULLs,
			self.EFEs,
		},
	}

	self.update_gPrims = efeProgram:kernel{
		name = 'update_gPrims',
		argsOut = {
			self.gPrims,
		},
		argsIn = {
			self.partial_gPrim_of_Phis,
		},
	}

	-- used by updateConjRes and updateGMRes
	self.calc_EinsteinLLs = efeProgram:kernel{
		name = 'calc_EinsteinLLs',
		argsOut = {
			-- don't provide an actual buffer here
			-- the conjResSolver will provide its own
			{name='EinsteinLLs', type='real4s4', obj=true},
		},
		argsIn = {
			self.gPrims,
			self.gLLs,
			self.gUUs,
			self.GammaULLs,
		},
	}
	self.calc_8piTLLs = efeProgram:kernel{
		name = 'calc_8piTLLs',
		argsOut = {self._8piTLLs},
		argsIn = {self.TPrims, self.gLLs},
	}

	self.init_gPrims = efeProgram:kernel{
		name = 'init_gPrims',
		argsOut = {
			self.gPrims,
		},
		argsIn = {
			{type='int'},
		},
	}


	print'done compiling code!'

	self.displayVarIndex = 1
	self.displayCode = self.displayVars[self.displayVarIndex].body
	self:refreshDisplayKernel()	-- rebuild the displayKernel
	self:refreshDisplayVar()	-- set displayVarIndex of the displayKernel, and set solver.displayCode
end

function EFESolver:refreshDisplayKernel()
	timer('refreshDisplayKernel', function()
		-- if we got bad display code then don't crash the whole app
		xpcall(function()

-- [=[ get rid of multiple constants 
self.code = table{
	self:getTypeCode(),
	self.gPrim_code,
	self.TPrim_code,
	self:template([[
//macro for the index
#define globalInt4()	(int4)((int)get_global_id(0), (int)get_global_id(1), (int)get_global_id(2), 0)

//macros for arbitrary sizes
#define indexForInt4ForSize(i, sx, sy, sz) (i.x + sx * (i.y + sy * i.z))
#define initKernelForSize(sx, sy, sz) \
	int4 i = globalInt4(); \
	if (i.x >= sx || i.y >= sy || i.z >= sz) return; \
	int index = indexForInt4ForSize(i, sx, sy, sz);


// THIS DIFFERS FROM THE ENV CODE
// ONLY COMPILE THAT ONCE, EVERYONE ELSE USE EXTERNS TO ACCESS IT
extern constant const int dim;
extern constant const int4 size;
extern constant const int4 stepsize;

//macros for the base domain
#define indexForInt4(i)	indexForInt4ForSize(i, size.x, size.y, size.z)
#define initKernel()	initKernelForSize(size.x, size.y, size.z)
]], {
		size = self.base.size,
	}),
}:concat'\n'
--]=]			
			local displayProgram = self:program{
				-- TODO can I use cl/obj/kernel.lua's codegen + link to other cl programs?
				cacheFileName = 'display',
				code = self:template[[
#include "efe.funcs.h"
#include "calcVars.h"

kernel void display(
	global float * const texCLBuf,
	int const displayVarIndex,
	global <?=TPrim_t?> const * const TPrims,
	global gPrim_t const * const gPrims,
	global real4s4 const * const gLLs,
	global real4s4 const * const gUUs,
	global real4x4s4 const * const GammaULLs,
	global real4s4 const * const EFEs
) {
	initKernel();

	if (displayVarIndex == 0) {
<?=solver.displayCode?>
<? for i,var in ipairs(solver.displayVars) do
?>	} else if (displayVarIndex == <?=i?>) {
<?=var.body?>
<? end
?>	}
}
]],
			}
			--[[ can't just link multiple bc files into a spv file ...
			displayProgram:compile{
				linkOptions = ('%q'):format('cache/efe.bc'),
			}
			--]]
			-- [[ maybe i have to do this? 
			-- TODO sort out how to specify this in the CLClangProgram class
			displayProgram.cacheFileCL = 'cache/display.cl'
			displayProgram.cacheFileBC = 'cache/display.bc'	-- cl -> bc 
			displayProgram.cacheFileOutBC = 'cache/display-out.bc'	-- bc + other bc's -> merged bc 
			displayProgram.cacheFileSPV = 'cache/display.spv'	-- merged bc -> spv 
			--displayProgram:build()
			local displayCode = self:template(displayProgram:getCode())
			writeIfChanged(displayProgram.cacheFileCL, displayCode)
			makeTargets{
				verbose = true,
				{
					srcs = {displayProgram.cacheFileCL},
					dsts = {displayProgram.cacheFileBC},
					rule = function(rule)
						-- cache/display.cl => cache/display.bc
						assert(#rule.dsts == 1)
						assert(#rule.srcs == 1)
						clangCompile(rule.dsts[1], rule.srcs[1], args and arg.buildOption or nil)
					end,
				},
				{
					srcs = {'cache/efe.bc', displayProgram.cacheFileBC},
					dsts = {displayProgram.cacheFileOutBC},
					rule = function()
						-- cache/efe.bc + cache/display.bc => cache/display-out.bc
						exec(table{
							'llvm-link',
							('%q'):format'cache/efe.bc',
							('%q'):format(displayProgram.cacheFileBC),
							'-o', ('%q'):format(displayProgram.cacheFileOutBC),
						}:concat' ')
					end,
				},
				{
					srcs = {displayProgram.cacheFileOutBC},
					dsts = {displayProgram.cacheFileSPV},
					rule = function()
						-- cache/display-out.bc => cache/display.spv
						--displayProgram:link()
						exec(table{
							'llvm-spirv',
							('%q'):format(displayProgram.cacheFileOutBC),
							'-o',
							('%q'):format(displayProgram.cacheFileSPV),
						}:concat' ')
					end,
				},
			}:run(displayProgram.cacheFileSPV)
			displayProgram.IL = assert(path(displayProgram.cacheFileSPV):read())
			CLClangProgram.super.compile(displayProgram)
			--]]
print'displayProgram'
print'CL_PROGRAM_KERNEL_NAMES:'
print(require 'ext.tolua'(displayProgram.obj:getInfo'CL_PROGRAM_KERNEL_NAMES'))
			self.displayKernel = displayProgram:kernel{
				name = 'display',
				setArgs = {
					self.texCLBuf,
					{type='int', name='displayVarIndex'},
					assert(self.TPrims),
					assert(self.gPrims),
					assert(self.gLLs),
					assert(self.gUUs),
					assert(self.GammaULLs),
					assert(self.EFEs),
				},
			}

			self:updateTex()
		end, function(err)
			print(err..'\n'..debug.traceback())
		end)
	end)
end

function EFESolver:refreshDisplayVar()
	local displayVar = self.displayVars[self.displayVarIndex]
	-- update the displayCode for the gui, 
	--  but don't rebuild , intead just change the displayVarIndex
	--  only rebuild when changing custom code
	self.displayCode = displayVar.body
	
	local displayVarIndex = ffi.new('int[1]') 
	displayVarIndex[0] = self.displayVarIndex
	self.displayKernel.obj:setArg(1, displayVarIndex)
	self:updateTex()
end

function EFESolver:resetState()
--print'init_gPrims'
	local initCondIndex = ffi.new('int[1]') 
	initCondIndex[0] = self.initCond
	self.init_gPrims.obj:setArg(1, initCondIndex)
	self.init_gPrims()	-- initialize gPrims
--self:printbuf'gPrims'

--print'init_TPrims'
	self.init_TPrims()	-- initialize TPrims
--self:printbuf'TPrims'

	self:updateAux()	-- calc gLLs, gUUs, GammaLLs, EFEs
--self:printbuf'gLLs'
--self:printbuf'gUUs'
--self:printbuf'GammaULLs'
--self:printbuf'EFEs'

	self:updateTex()
--print('residual', self:calcBufferNorm())

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

--[[
TODO here - gauge conditions - maybe?
first attempt -- harmonic slicing condition:
d/dt α = (∂_t - L_β) α = -α^2 K^ij γ_ij
α_,t - β^k α_,k = -α^2 K^ij γ_ij

but K_ij isn't a state variable (or is it?)
K_ab = -⟂ ∇_(a n_b)
K_ab = -⟂ n_(b;a)
K_ab = -⟂ (n_(b,a) - Γ^c_ba n_c)
K_ab = -γ_a^u γ_b^v (n_(v,u) - Γ^c_vu n_c)
... for γ_ab = g_ab + n_a n_b
for distinct time coordinate: n_a = [-α, 0]
so for distinct time coordinate:
	K_ab = -γ_a^u γ_b^v (n_(v,u) - Γ^c_vu n_c)
	K_ij = -α Γ^t_ij ... the 4D connection along the time coordinate
	... if we keep going ...
	K_ij =
	= -1/2 α g^tu (g_ui,j + g_uj,i - g_ij,u)
	= -1/2 α g^tu (g_ui,j + g_uj,i - g_ij,u)
	= -1/2 α (g^tt (g_ti,j + g_tj,i - g_ij,t) + g^tk (g_ki,j + g_kj,i - g_ij,k))
	= -1/2 α (-1/α^2 (β_i,j + β_j,i - γ_ij,t) + β^k / α^2 (γ_ki,j + γ_kj,i - γ_ij,k))
	= -1/2 (β^k γ_ki,j + β^k γ_kj,i - β^k γ_ij,k - β_i,j - β_j,i + γ_ij,t) / α
	... using β_i,j = (β^k γ_ki)_,j = β^k_,j γ_ki + β^k γ_ki,j
	= -1/2 (β^k γ_ki,j + β^k γ_kj,i - β^k γ_ij,k - β^k_,j γ_ki - β^k γ_ki,j - β^k_,i γ_kj - β^k γ_kj,i + γ_ij,t) / α
	= -1/2 (- β^k γ_ij,k - β^k_,j γ_ki - β^k_,i γ_kj + γ_ij,t) / α

substitute back into harmonic slicing condition:
α_,t - β^k α_,k = -α^2 K^ij γ_ij
α_,t - β^k α_,k = α^3 Γ^t_ij γ_ij
	... as a constraint ...
Φ = 1/2 |α_,t - β^k α_,k - α^3 Γ^t_ij γ_ij|^2

	... then update?
	∂g_ab/∂t = -∂Φ/∂g_ab
--]]

	self.iteration = self.iteration + 1
end

-- calc norm of self.EFEs
function EFESolver:calcBufferNorm(srcbuf)
	srcbuf = srcbuf or self.EFEs
	--[[ hmmm....
	self.jfnkSolver.args.scale(
		self.tmpBuf,
		srcbuf,
		1 / G)
	return self.conjResSolver.args.dot(self.tmpBuf, self.tmpBuf)
	--]]
	-- [[
	local cpu = srcbuf:toCPU()
	local sum = 0
	local m = ffi.sizeof(srcbuf.type) / ffi.sizeof'real'
	assert(m == 10)	-- real4s4
	local volume = tonumber(self.base.size:volume())
	for i=0,volume-1 do
		local fptr = ffi.cast('real*', cpu[i].s)
		for j=0,m-1 do
--print(i,j, cpu[i].s[j])
--			sum = sum + (fptr[j] / G)^2
			sum = sum + (fptr[j])^2
		end
	end
	return math.sqrt(sum) / volume
	--]]
end

function EFESolver:updateNewton()
	--[[
	iteration:

	∂/∂t(g_ab) = -∂Φ/∂g_ab
	for Φ = 1/2 Sum_ab (G_ab - 8 π T_ab)^2

	two approaches:
	1) do this per g_ab, so you only need to allocate as big as you would for solving the constraints themselves
	2) do this for g_ab as a whole, which would mean x4^2 symmetric = x10 allocation, but would take less kernel passes
	I'll try for 2 and hope I have enough memory
	--]]

	-- here's the newton update method
print'calc_partial_gPrim_of_Phis_kernel'
	self.calc_partial_gPrim_of_Phis_kernel()
--self:printbuf'partial_gPrim_of_Phis'
print('partial_gPrim_of_Phis norm '..self:calcBufferNorm(self.partial_gPrim_of_Phis))
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
			local residual = self:calcBufferNorm()
print(('lambda=%.16e residual=%.16e'):format(lambda, residual))
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
print(('fwd lambda=%.16e residual=%.16e'):format(lambdaFwd, residualFwd))
		local lambdaRev, residualRev = bisect(0, -self.updateLambda)
print(('rev lambda=%.16e residual=%.16e'):format(lambdaRev, residualRev))

		self.conjResSolver.args.copy(self.gPrims, self.gPrimsCopy)
		local lambda, residual
		if residualFwd < residualRev then
			lambda = lambdaFwd
			residual = residualFwd
		else
			lambda = lambdaRev
			residual = residualRev
		end
print(('using lambda=%.16e residual=%.16e'):format(lambda, residualRev))
		lambdaPtr[0] = lambda
		self.update_gPrims.obj:setArg(2, lambdaPtr)
		self.update_gPrims()
	else	-- no line search
		-- update gPrims from ∂Φ/∂g_ab
print('self.update_gPrims')
print('lambda', self.updateLambda)
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
print('residual', self:calcBufferNorm())
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

--print'calc_EFEs'
	self.calc_EFEs()
--self:printbuf'EFEs'
end

function EFESolver:updateTex()
	self.displayKernel()

	-- TODO run a reduce on the display var stuff
	-- get the min and max
	-- then rescale the data according

-- notice some of my results might not survive the double->float cast
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

function EFESolver:printbuf(name)
	print(tostring(name)..':')
	local buf = assert(self[name])
	local m = ffi.sizeof(buf.type) / ffi.sizeof'real'
	local cpu = buf:toCPU()
	for i=0,tonumber(self.base.size:volume()-1) do
		local z = i
		local x = tonumber(z % self.base.size.x)
		z = tonumber((z - x) / self.base.size.x)
		local y = tonumber(z % self.base.size.y)
		z = tonumber((z - y) / self.base.size.y)
		io.write(name..'('..x..', '..y..', '..z..') = '..tostring(cpu[i]))
		io.write(' ... reals:')
		local fptr = ffi.cast('real*', cpu[i].s)
		for j=0,m-1 do
			io.write(' ', fptr[j])
			if not math.isfinite(fptr[j]) then
				error("found a nan!")
			end
		end
		print()
	end
	print()
end

return EFESolver
