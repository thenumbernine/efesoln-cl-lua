//

<?
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

--[[
args:
	bufferName
	bufferType
	resultName
	resultType
--]]
local function finiteDifference(args)
	local order
	local bufferName = args.bufferName
	local bufferType = args.bufferType
	local resultName = args.resultName
	local resultType = args.resultType
	local order = 2	-- TODO put in config.lua
	local coeffs = assert(derivCoeffs[1][order])
?>	<?=resultType?> <?=resultName?> = <?=resultType?>_zero;
	<? for i=0,sDim-1 do ?>{
		<? for j,coeff in pairs(coeffs) do ?>{
			<?=bufferType?> const yL = (i.s<?=i?> - <?=j?> < 0)
				? <?=bufferType?>_zero	// lhs boundary condition
				: <?=bufferName?>[index - stepsize.s<?=i?> * <?=j?>];

			<?=bufferType?> const yR = (i.s<?=i?> + <?=j?> >= size.s<?=i?>)
				? <?=bufferType?>_zero	// rhs boundary condition
				: <?=bufferName?>[index + stepsize.s<?=i?> * <?=j?>];

			<?=resultName?>.s<?=i+1?> = <?=bufferType?>_add(
				<?=resultName?>.s<?=i+1?>,
				<?=bufferType?>_real_mul(
					<?=bufferType?>_sub(yR, yL),
					<?=coeff?> * inv_dx.s<?=i?>
				)
			);
		}<? end ?>
	}<? end ?>
<?
end
?>

// init

kernel void init_TPrims(
	global <?=TPrim_t?> * const TPrims
) {
	initKernel();

	global <?=TPrim_t?> * const TPrim = TPrims + index;

	*TPrim = (<?=TPrim_t?>){
<? if solver.body.useMatter then ?>
		.rho = 0,
		.eInt = 0,
		.P = 0,
<? 	if solver.body.useVel then ?>
		.v = real3_zero,
<? 	end
end
if solver.body.useEM then ?>
		.E = real3_zero,
		.B = real3_zero,
<? end ?>
	};

<?=solver.body.init and solver.body.init or ''?>

<? if solver.useFourPotential then ?>
	TPrim->AL = real4_real_mul(TPrim->JU, -1);
<? end ?>
}

// compute buffers to compute EFE

kernel void calc_gLLs_and_gUUs(
	global real4s4 * const gLLs,
	global real4s4 * const gUUs,
	global gPrim_t const * const gPrims
) {
	initKernel();

	gPrim_t const gPrim = gPrims[index];
	real const alpha = gPrim.alpha;
	real3 const betaU = gPrim.betaU;
	real3s3 const gammaLL = gPrim.gammaLL;

	real const alphaSq = alpha * alpha;
	real const invAlphaSq = 1. / alphaSq;
	real3 const betaL = real3s3_real3_mul(gammaLL, betaU);
	real const betaSq = real3_dot(betaL, betaU);

	global real4s4 * const gLL = gLLs + index;
	gLL->s00 = -alphaSq + betaSq;
	<? for i=0,sDim-1 do ?>
	gLL->s0<?=i+1?> = betaL.s<?=i?>;
		<? for j=i,sDim-1 do ?>
	gLL->s<?=i+1?><?=j+1?> = gammaLL.s<?=i?><?=j?>;
		<? end ?>
	<? end ?>

	real const det_gammaLL = real3s3_det(gammaLL);
	real3s3 const gammaUU = real3s3_inv(det_gammaLL, gammaLL);

	global real4s4 * const gUU = gUUs + index;
	gUU->s00 = -invAlphaSq;
	<? for i=0,sDim-1 do ?>
	gUU->s0<?=i+1?> = betaU.s<?=i?> * invAlphaSq;
		<? for j=i,sDim-1 do ?>
	gUU->s<?=i+1?><?=j+1?> = gammaUU.s<?=i?><?=j?> - betaU.s<?=i?> * betaU.s<?=j?> * invAlphaSq;
		<? end ?>
	<? end ?>
}

constant real4s4 const gLL_flat = (real4s4){
	.tt = -1, .tx = 0, .ty = 0, .tz = 0,
	.xx = 1, .xy = 0, .xz = 0,
	.yy = 1, .yz = 0,
	.zz = 1,
};

kernel void calc_GammaULLs(
	global real4x4s4 * const GammaULLs,
	global real4s4 const * const gLLs,
	global real4s4 const * const gUUs
) {
	initKernel();

	//g_ab,c := dgLLL.c.ab
	//here's where the finite difference stuff comes in ...
	//TODO modular finite-difference kernels & boundary conditions
#if 0	// auto-gen ... breaking
<?finiteDifference{
	bufferName = "gLLs",
	bufferType = "real4s4",
	resultName = "dgLLL",
	resultType = "real4x4s4",
}?>
#else
	real4x4s4 dgLLL;
	dgLLL.s0 = real4s4_zero;
	<? for i=0,sDim-1 do ?>{
		real4s4 const gLL_prev = (i.s<?=i?> - 1 < 0)
			? gLL_flat	// boundary condition:
			: gLLs[index - stepsize.s<?=i?>];
		
		real4s4 const gLL_next = (i.s<?=i?> + 1 >= size.s<?=i?> - 1)
			? gLL_flat	//boundary condition:
			: gLLs[index + stepsize.s<?=i?>];

		dgLLL.s<?=i+1?> = real4s4_real_mul(real4s4_sub(gLL_next, gLL_prev), .5 * inv_dx.s<?=i?>);
	}<? end ?>
#endif

	//Γ_abc := GammaLLL.a.bc
	//Γ_abc = 1/2 (g_ab,c + g_ac,b - g_bc,a)
	real4x4s4 const GammaLLL = (real4x4s4){
	<? for a=0,stDim-1 do ?>
		<? for b=0,stDim-1 do ?>
			<? for c=b,stDim-1 do ?>
		.s<?=a?>.s<?=b?><?=c?> = .5 * (
			dgLLL.s<?=c?>.s<?=sym(a,b)?>
			+ dgLLL.s<?=b?>.s<?=sym(a,c)?>
			- dgLLL.s<?=a?>.s<?=sym(b,c)?>),
			<? end ?>
		<? end ?>
	<? end ?>};

	//Γ^a_bc = GammaULL.a.bc
	//Γ^a_bc = g^ad Γ_dbc
	real4s4 const gUU = gUUs[index];
	global real4x4s4 * const GammaULL = GammaULLs + index;
	<? for a=0,stDim-1 do ?>
		<? for b=0,stDim-1 do ?>
			<? for c=b,stDim-1 do ?>
	GammaULL->s<?=a?>.s<?=b?><?=c?> = 0.
				<? for d=0,stDim-1 do ?>
		+ gUU.s<?=sym(a,d)?> * GammaLLL.s<?=d?>.s<?=b?><?=c?>
				<? end ?>;
			<? end ?>
		<? end ?>
	<? end ?>
}

/*
J_a = (rho, j_i)

flat space:

F_uv^,v = 4 pi J_u
F_uv = A_v,u - A_u,v
A_v,u^v - A_u,v^v = 4 pi J^u

use the gauge A_v,^v = 0

A_u,v^v = -4 pi J_u

curved space:
A_v;u^v - A_u;v^v + R^u_v A^v = 4 pi J^u

use gauge A^u_;u = 0

-A_u;v^v + R^u_v A^v = 4 pi J^u
A_a;u^u - R_a^u A_u = -4 pi J_a

to enforce the gauge, A^u_;u = 0
we need to subtract the potential gradient component of A
...or don't use the gauge :-p

D A_a = -J_a
A_a = -D^-1 J_a for some D...

what is D?

A_v;u^v - A_u;v^v + R_u^v A_v = 4 pi J_u
= g^vw (A_v;uw - A_u;vw + R_uv A_w) = 4 pi J_u
= g^vw (A_v;u;w - A_u;v;w + R_uv A_w) = 4 pi J_u
= g^vw (
	(A_v,u - Gamma^r_vu A_r)_;w
	- (A_u,v - Gamma^r_uv A_r)_;w
	+ R_uv A_w) = 4 pi J_u
= g^vw (
	(A_v,u - Gamma^r_vu A_r)_,w
	- (A_s,u - Gamma^r_su A_r) Gamma^s_vw
	- (A_v,s - Gamma^r_vs A_r) Gamma^s_uw
	- (A_u,v - Gamma^r_uv A_r)_,w
	+ (A_u,s - Gamma^r_us A_r) Gamma^s_vw
	+ (A_s,v - Gamma^r_sv A_r) Gamma^s_uw
	+ R_uv A_w) = 4 pi J_u
= g^vw (
	A_v,uw
	- A_u,vw
	- Gamma^s_vw A_s,u
	+ Gamma^s_uw A_s,v
	- Gamma^s_uw A_v,s
	+ Gamma^s_vw A_u,s
	+ R_uv A_w) = 4 pi J_u


or how about I enforce the A^u_;u = 0 constraint as well?

so every iteration that we converge J_a for -1/(4pi) (g^vw D_v D_w delta^u_a  - R_a^u) A_u = J_a
we also constrain A^u_;u = 0, which means divergence-free,
which means subtract out the potential (in curved space)
... how?
*/
<? if solver.useFourPotential then ?>
kernel void solveAL(
	global <?=TPrim_t?> * const TPrims
) {
	initKernel();

	real4 skewSum = real4_zero;
	<? for i=0,sDim-1 do ?>{

		<?=TPrim_t?> TPrim_prev;
		if (i.s<?=i?> > 0) {
			int4 iL = i;
			--iL.s<?=i?>;
			int const indexL = indexForInt4(iL);
			TPrim_prev = TPrims[indexL];
		} else {
			// boundary condition
			TPrim_prev = TPrims[index];
		}

		<?=TPrim_t?> TPrim_next;
		if (i.s<?=i?> < size.s<?=i?> - 1) {
			int4 iR = i;
			++iR.s<?=i?>;
			int indexR = indexForInt4(iR);
			TPrim_next = TPrims[indexR];
		} else {
			// boundary condition
			TPrim_next = TPrims[index];
		}

		skewSum = real4_add(
			skewSum,
			real4_real_mul(
				real4_add(TPrim_prev.JU, TPrim_next.JU),
				inv_dx.s<?=i?> * inv_dx.s<?=i?>)
			);
	}<? end ?>

	real const diag = -2. * (0
<? for i=0,solver.stDim-1 do ?>
		+ 1. / (dx<?=i?> * dx<?=i?>)
<? end ?>
	);

	global <?=TPrim_t?> * const TPrim = TPrims + index;

	TPrim->AL = real4_sub(TPrim->AL, skewSum) / diag;
}
<? end ?>

//used by linearized solvers of G x = 8 pi T
// where G is the (non)linear function G_ab
// T is T_ab (which is also a function of x but don't tell)
// and x is gPrims
kernel void calc_EinsteinLLs(
	global real4s4 * const EinsteinLLs,
	global real4s4 const * const gLLs,
	global real4s4 const * const gUUs,
	global real4x4s4 const * const GammaULLs
) {
	initKernel();
	EinsteinLLs[index] = calc_EinsteinLL(gLLs, gUUs, GammaULLs);
}

kernel void calc_8piTLLs(
	global real4s4 * const _8piTLLs,
	global <?=TPrim_t?> const * const TPrims,
	global real4s4 const * const gLLs
) {
	initKernel();
	_8piTLLs[index] = calc_8piTLL(gLLs[index], TPrims[index]);
}

kernel void calc_EFEs(
	global real4s4 * const EFEs,
	global <?=TPrim_t?> const * const TPrims,
	global real4s4 const * const gLLs,
	global real4s4 const * const gUUs,
	global real4x4s4 const * const GammaULLs
) {
	initKernel();
	<?=TPrim_t?> const TPrim = TPrims[index];
	real4s4 const gLL = gLLs[index];
	real4s4 const EinsteinLL = calc_EinsteinLL(gLLs, gUUs, GammaULLs);
	real4s4 const _8piTLL = calc_8piTLL(gLL, TPrim);
	// EFEs(x) = G_ab(x) - 8 π T_ab(x)
	EFEs[index] = real4s4_sub(EinsteinLL, _8piTLL);
}
