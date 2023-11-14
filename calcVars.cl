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
	global sym4 * const gLLs,
	global sym4 * const gUUs,
	global gPrim_t const * const gPrims
) {
	initKernel();

	global gPrim_t const * const gPrim = gPrims + index;

	real alphaSq = gPrim->alpha * gPrim->alpha;
	real3 betaL = sym3_real3_mul(gPrim->gammaLL, gPrim->betaU);
	real betaSq = real3_dot(betaL, gPrim->betaU);

	global sym4 * const gLL = gLLs + index;
	gLL->s00 = -alphaSq + betaSq;
	<? for i=0,sDim-1 do ?>
	gLL->s0<?=i+1?> = betaL.s<?=i?>;
		<? for j=i,sDim-1 do ?>
	gLL->s<?=i+1?><?=j+1?> = gPrim->gammaLL.s<?=i?><?=j?>;
		<? end ?>
	<? end ?>

	real det_gammaLL = sym3_det(gPrim->gammaLL);
	sym3 gammaUU = sym3_inv(det_gammaLL, gPrim->gammaLL);

	global sym4 * const gUU = gUUs + index;
	gUU->s00 = -1./alphaSq;
	<? for i=0,sDim-1 do ?>
	gUU->s0<?=i+1?> = gPrim->betaU.s<?=i?> / alphaSq;
		<? for j=i,sDim-1 do ?>
	gUU->s<?=i+1?><?=j+1?> = gammaUU.s<?=i?><?=j?> - gPrim->betaU.s<?=i?> * gPrim->betaU.s<?=j?> / alphaSq;
		<? end ?>
	<? end ?>
}

constant sym4 const gLL_flat = (sym4){
	.tt = -1, .tx = 0, .ty = 0, .tz = 0,
	.xx = 1, .xy = 0, .xz = 0,
	.yy = 1, .yz = 0,
	.zz = 1,
};

kernel void calc_GammaULLs(
	global real4x4s4 * const GammaULLs,
	global sym4 const * const gLLs,
	global sym4 const * const gUUs
) {
	initKernel();

	//g_ab,c = dgLLL.c.ab
	//here's where the finite difference stuff comes in ...
	//TODO modular finite-difference kernels & boundary conditions
	real4x4s4 dgLLL;
	dgLLL.s0 = sym4_zero;
	<? for i=0,sDim-1 do ?>{
		sym4 gLL_prev;
		if (i.s<?=i?> > 0) {
			int4 iL = i;
			--iL.s<?=i?>;
			int const indexL = indexForInt4(iL);
			gLL_prev = gLLs[indexL];
		} else {
			// boundary condition:
			gLL_prev = gLL_flat;
		}
		
		sym4 gLL_next;
		if (i.s<?=i?> < size.s<?=i?> - 1) {
			int4 iR = i;
			++iR.s<?=i?>;
			int const indexR = indexForInt4(iR);
			gLL_next = gLLs[indexR];
		} else {
			//boundary condition:
			gLL_next = gLL_flat;
		}
		
		dgLLL.s<?=i+1?> = (sym4){
		<? for a=0,stDim-1 do ?>
			<? for b=a,stDim-1 do ?>
			.s<?=a..b?> = (gLL_next.s<?=a..b?> - gLL_prev.s<?=a..b?>) * .5 * inv_dx.s<?=i?>,
			<? end ?>
		<? end ?>
		};
	}<? end ?>

	//Γ_abc = GammaLLL.a.bc
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
	sym4 const gUU = gUUs[index];
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
		int4 iL = i;
		iL.s<?=i?> = max(i.s<?=i?> - 1, 0);
		int const indexL = indexForInt4(iL);
		global <?=TPrim_t?> * const TPrim_prev = TPrims + indexL;

		skewSum = real4_add(
			skewSum,
			real4_real_mul(
				TPrim_prev->JU,
				inv_dx.s<?=i?> * inv_dx.s<?=i?>)
			);
		
		int4 iR = i;
		iR.s<?=i?> = min(i.s<?=i?> + 1, size.s<?=i?> - 1);
		int indexR = indexForInt4(iR);
		global <?=TPrim_t?> * const TPrim_next = TPrims + indexR;

		skewSum = real4_add(
			skewSum,
			real4_real_mul(
				TPrim_next->JU,
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
	global sym4 * const EinsteinLLs,
	global sym4 const * const gLLs,
	global sym4 const * const gUUs,
	global real4x4s4 const * const GammaULLs
) {
	initKernel();
	EinsteinLLs[index] = calc_EinsteinLL(gLLs, gUUs, GammaULLs);
}

kernel void calc_8piTLLs(
	global sym4 * const _8piTLLs,
	global <?=TPrim_t?> const * const TPrims,
	global sym4 const * const gLLs
) {
	initKernel();
	_8piTLLs[index] = calc_8piTLL(gLLs[index], TPrims[index]);
}

kernel void calc_EFEs(
	global sym4 * const EFEs,
	global <?=TPrim_t?> const * const TPrims,
	global sym4 const * const gLLs,
	global sym4 const * const gUUs,
	global real4x4s4 const * const GammaULLs
) {
	initKernel();
	<?=TPrim_t?> const TPrim = TPrims[index];
	sym4 const gLL = gLLs[index];
	sym4 const EinsteinLL = sym4_zero;//calc_EinsteinLL(gLLs, gUUs, GammaULLs);
	sym4 const _8piTLL = calc_8piTLL(gLL, TPrim);	// getting nans
	// EFEs(x) = G_ab(x) - 8 π T_ab(x)
	EFEs[index] = sym4_sub(EinsteinLL, _8piTLL);
}
