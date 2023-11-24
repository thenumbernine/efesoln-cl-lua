// different initial conditions / boundary conditions of gPrims
// this is the only code that uses solver.body

#define new_gPrim_flat() ((gPrim_t){\
	.alpha = 1,\
	.betaU = new_real3_zero(),\
	.gammaLL = real3s3_ident,\
})

constant gPrim_t const gPrim_flat = new_gPrim_flat();

gPrim_t calc_gPrim_flat(real3 const x) {
	return gPrim_flat;
}

gPrim_t calc_gPrim_stellar_Schwarzschild(real3 const x) {
	gPrim_t gPrim = calc_gPrim_flat(x);

	real const r = real3_len(x);
	real const radius = <?=solver.body.radius?>;
	real const mass = <?=solver.body.mass?>;
	real const density = <?=solver.body.density?>;

	real const matterRadius = (real)min(r, radius);
	real const volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
	real const m = density * volumeOfMatterRadius;	// m^3

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
	gPrim.alpha = r > radius
		? sqrt(1 - 2*mass/r)
		: (1.5 * sqrt(1 - 2*mass/radius) - .5 * sqrt(1 - 2*mass*r*r/(radius*radius*radius)));

	for (int i = 0; i < sDim; ++i) {
		gPrim.betaU.s[i] = 0;
		for (int j = i; j < sDim; ++j) {
			gPrim.gammaLL.s[sym3[i][j]] = (i == j ? 1. : 0.) + x.s[i]/r * x.s[j]/r * 2*m/(r - 2*m);
		}
	}

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
	real const dr_alpha = r > radius
		? (mass / (r * r * sqrt(1. - 2. * mass / r)))
		: (mass * r / (radius * radius * radius * sqrt(1. - 2. * mass * r * r / (radius * radius * radius))));
	real const dr_m = r > radius ? 0 : (4. * M_PI * r * r * density);
	MetricPrims & dt_metricPrims = dt_metricPrimGrid(index);
	real const L = 149.6e+9;	//distance from earth to sun, in m
	//real omega = 0; //no rotation
	//real omega = 2. * M_PI / (60. * 60. * 24. * 365.25) / c;	//one revolution per year in m^-1
	//real omega = 1;	//angular velocity of the speed of light
	real const omega = c;	//I'm trying to find a difference ...
	real const t = 0;	//where the position should be.  t=0 means the body is moved by [L, 0], and its derivatives are along [0, L omega]
	Vector<real,2> dt_xHat(L * omega * sin(omega * t), -L * omega * cos(omega * t));
	dt_metricPrims.alpha = dr_alpha * (xi(0)/r * dt_xHat(0) + xi(1)/r * dt_xHat(1));
	for (int i = 0; i < sDim; ++i) {
		dt_metricPrims.betaU(i) = 0;
	}
	for (int i = 0; i < sDim; ++i) {
		for (int j = 0; j < sDim; ++j) {
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
	for (int i = 0; i < sDim; ++i) {
		//negate all gravity by throttling the change in space/time coupling of the metric
		real dm_dr = 0;
		betaL(i) = -(2*m * (r - 2*m) + 2 * dm_dr * r * (2*m - r)) / (2 * r * r * r) * xi(i)/r;
	}
	TensorSUsub gammaUU = inverse(metricPrims.gammaLL);
	for (int i = 0; i < sDim; ++i) {
		real sum = 0;
		for (int j = 0; j < sDim; ++j) {
			sum += gammaUU(i,j) * betaL(j);
		}
		dt_metricPrims.betaU(i) = sum;
	}
#endif

	return gPrim;
}

gPrim_t calc_gPrim_stellar_Kerr_Newman(real3 const x) {
	gPrim_t gPrim = calc_gPrim_flat(x);

	real const radius = <?=solver.body.radius?>;
	real const mass = <?=solver.body.mass?>;
	real const density = <?=solver.body.density?>;

	real const angularVelocity = 2. * M_PI / (60. * 60. * 24.) / c;	//angular velocity, in m^-1
	real const inertia = 2. / 5. * mass * radius * radius;	//moment of inertia about a sphere, in m^3
	real const angularMomentum = inertia * angularVelocity;	//angular momentum in m^2
	real const a = angularMomentum / mass;	//m

	//real r is the solution of (x*x + y*y) / (r*r + a*a) + z*z / (r*r) = 1
	// r^4 - (x^2 + y^2 + z^2 - a^2) r^2 - a^2 z^2 = 0
	real const RSq_minus_aSq = real3_lenSq(x) - a*a;
	//so we have two solutions ... which do we use?
	//from gnuplot it looks like the two of these are the same ...
	real const r = sqrt((RSq_minus_aSq + sqrt(RSq_minus_aSq * RSq_minus_aSq + 4.*a*a*x.z*x.z)) / 2.);	//use the positive root

	//should I use the Kerr-Schild 'r' coordinate?
	//well, if 'm' is the mass enclosed within the coordinate
	// and that determines 'a', the angular momentum per mass within the coordinate (should it?)
	// then we would have a circular definition
	//real R = real3_len(x);
	real const matterRadius = min(r, radius);
	real const volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
	real const m = density * volumeOfMatterRadius;	// m^3

	real const Q = 0;	//charge
	real const H = (r*m - Q*Q/2.)/(r*r + a*a*x.z*x.z/(r*r));

	//3.4.33 through 3.4.35 of Alcubierre "Introduction to 3+1 Numerical Relativity"

	/*TODO fix this for the metric within the star
	 in other news, this is an unsolved problem!
	https://arxiv.org/pdf/1503.02172.pdf section 3.11
	https://arxiv.org/pdf/1410.2130.pdf section 4.2 last paragraph
	*/
	//metricPrims.alpha = 1./sqrt(1. + 2*H);
	gPrim.alpha = sqrt(1. - 2*H/(1+2*H) );

	real3 const l = _real3( (r*x.x + a*x.y)/(r*r + a*a), (r*x.y - a*x.x)/(r*r + a*a), x.z/r );
<?
for i=0,sDim-1 do
?>	gPrim.betaU.s<?=i?> = 2. * H * l.s<?=i?> / (1. + 2. * H);
<?
	for j=i,sDim-1 do
?>	gPrim.gammaLL.s<?=i..j?> = <?if i==j then?>1. + <?end?>2. * H * l.s<?=i?> * l.s<?=j?>;
<?	end
end
?>
	return gPrim;
}

real4s4 calc_gLL_from_gPrim(
	gPrim_t const gPrim
) {
	real const alpha = gPrim.alpha;
	real3 const betaU = gPrim.betaU;
	real3s3 const gammaLL = gPrim.gammaLL;

	real const alphaSq = alpha * alpha;
	real3 const betaL = real3s3_real3_mul(gammaLL, betaU);
	real const betaSq = real3_dot(betaL, betaU);

	real4s4 gLL;
	gLL.s00 = -alphaSq + betaSq;
	for (int i = 0; i < sDim; ++i) {
		gLL.s[sym4[0][i+1]] = betaL.s[i];
		for (int j = i; j < sDim; ++j) {
			gLL.s[sym4[i+1][j+1]] = gammaLL.s[sym3[i][j]];
		}
	}

	return gLL;
}

real4s4 calc_gUU_from_gPrim(
	gPrim_t const gPrim
) {
	real const alpha = gPrim.alpha;
	real const invAlphaSq = 1. / (alpha * alpha);
	real3 const betaU = gPrim.betaU;

	real3s3 const gammaLL = gPrim.gammaLL;
	real const det_gammaLL = real3s3_det(gammaLL);
	real3s3 const gammaUU = real3s3_inv(det_gammaLL, gammaLL);

	real4s4 gUU;
	gUU.s00 = -invAlphaSq;
	for (int i = 0; i < sDim; ++i) {
		gUU.s[sym4[0][i+1]] = betaU.s[i] * invAlphaSq;
		for (int j = i; j < sDim; ++j) {
			gUU.s[sym4[i+1][j+1]] = gammaUU.s[sym3[i][j]] - betaU.s[i] * betaU.s[j] * invAlphaSq;
		}
	}
	return gUU;
}

// init

kernel void init_gPrims(
	global gPrim_t * const gPrims
) {
	initKernel();
	real3 const x = getX(i);
	gPrims[index] = <?=solver.initConds[solver.initCond].code?>(x);
}

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
	gLLs[index] = calc_gLL_from_gPrim(gPrim);
	gUUs[index] = calc_gUU_from_gPrim(gPrim);
}

kernel void calc_GammaULLs(
	global real4x4s4 * const GammaULLs,
	global real4s4 const * const gLLs,
	global real4s4 const * const gUUs
) {
	initKernel();

	//partial_xU_of_gLL.c.ab := ∂/∂x^c(g_ab) = g_ab,c
<?=solver:finiteDifference{
	bufferName = "gLLs",
	srcType = "4s4",
	resultName = "partial_xU_of_gLL",
	--getBoundary = function(args) return "new_real4s4_Minkowski()" end,
	getBoundary = function(args) return "real4s4_Minkowski" end,
}?>

	//Γ_abc := GammaLLL.a.bc
	//Γ_abc = 1/2 (g_ab,c + g_ac,b - g_bc,a)
	real4x4s4 GammaLLL;
	for (int a = 0; a < stDim; ++a) {
		for (int b = 0; b < stDim; ++b) {
			for (int c = b; c < stDim; ++c) {
				int const bc = sym4[b][c];
				GammaLLL.s[a].s[bc] = .5 * (
					  partial_xU_of_gLL.s[c].s[sym4[a][b]]
					+ partial_xU_of_gLL.s[b].s[sym4[a][c]]
					- partial_xU_of_gLL.s[a].s[bc]
				);
			}
		}
	}

	//Γ^a_bc = GammaULL.a.bc
	//Γ^a_bc = g^ad Γ_dbc
	real4s4 const gUU = gUUs[index];
	GammaULLs[index] = real4s4_real4x4s4_mul(gUU, GammaLLL);
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

<?
--[[ TODO this is 2nd order, and the middle is missing, because it's an inverse to a discrete Laplacian solved with Jacobi iteration
=solver:finiteDifference{
	bufferName = "TPrims",
	getValue = function(args) return "TPrims["..args.index.."].JU" end,
	valueType = "real4",
}
--]]
?>

	real4 skewSum = real4_zero;
	for (int i = 0; i < sDim; ++i) {

		<?=TPrim_t?> TPrim_prev;
		if (i.s[i] > 0) {
			int4 iL = i;
			--iL.s[i];
			int const indexL = indexForInt4(iL);
			TPrim_prev = TPrims[indexL];
		} else {
			// boundary condition
			TPrim_prev = TPrims[index];
		}

		<?=TPrim_t?> TPrim_next;
		if (i.s[i] < size.s[i] - 1) {
			int4 iR = i;
			++iR.s[i];
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
				inv_dx.s[i] * inv_dx.s[i])
			);
	}

	real const diag = -2. * (0.
		+ 1. / (dx.s0 * dx.s0)
		+ 1. / (dx.s1 * dx.s1)
		+ 1. / (dx.s2 * dx.s2)
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
	EinsteinLLs[index] = calc_EinsteinLL(i, gLLs, gUUs, GammaULLs);
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
	real4s4 const EinsteinLL = calc_EinsteinLL(i, gLLs, gUUs, GammaULLs);
	real4s4 const _8piTLL = calc_8piTLL(gLL, TPrim);
	// EFEs(x) = G_ab(x) - 8 π T_ab(x)
//debugging
	//EFEs[index] = real4s4_zero;	//no nans
	//EFEs[index] = _8piTLL;	//no nans
	//EFEs[index] = EinsteinLL;	//no nans
	EFEs[index] = real4s4_sub(EinsteinLL, _8piTLL);
	/* no nans
	for (int a = 0; a < stDim; ++a) {
		for (int b = a; b < stDim; ++b) {
			EFEs[index].s[sym4[a][b]] = EinsteinLL.s[sym4[a][b]] - _8piTLL.s[sym4[a][b]];
		}
	}
	*/
	/* inlined real4s4_sub ... has nans
	EFEs[index] = (real4s4){
	.s00 = real_sub(EinsteinLL.s00, _8piTLL.s00),
	.s01 = real_sub(EinsteinLL.s01, _8piTLL.s01),
	.s02 = real_sub(EinsteinLL.s02, _8piTLL.s02),
	.s03 = real_sub(EinsteinLL.s03, _8piTLL.s03),
	.s11 = real_sub(EinsteinLL.s11, _8piTLL.s11),
	.s12 = real_sub(EinsteinLL.s12, _8piTLL.s12),
	.s13 = real_sub(EinsteinLL.s13, _8piTLL.s13),
	.s22 = real_sub(EinsteinLL.s22, _8piTLL.s22),
	.s23 = real_sub(EinsteinLL.s23, _8piTLL.s23),
	.s33 = real_sub(EinsteinLL.s33, _8piTLL.s33),
	};
	*/
	/* inlined'd real4s4_sub v2 ... has nans
	EFEs[index].s00 = real_sub(EinsteinLL.s00, _8piTLL.s00);
	EFEs[index].s01 = real_sub(EinsteinLL.s01, _8piTLL.s01);
	EFEs[index].s02 = real_sub(EinsteinLL.s02, _8piTLL.s02);
	EFEs[index].s03 = real_sub(EinsteinLL.s03, _8piTLL.s03);
	EFEs[index].s11 = real_sub(EinsteinLL.s11, _8piTLL.s11);
	EFEs[index].s12 = real_sub(EinsteinLL.s12, _8piTLL.s12);
	EFEs[index].s13 = real_sub(EinsteinLL.s13, _8piTLL.s13);
	EFEs[index].s22 = real_sub(EinsteinLL.s22, _8piTLL.s22);
	EFEs[index].s23 = real_sub(EinsteinLL.s23, _8piTLL.s23);
	EFEs[index].s33 = real_sub(EinsteinLL.s33, _8piTLL.s33);
	*/
	/* another try ... works
	for (int ab = 0; ab < 10; ++ab) {
		EFEs[index].s[ab] = real_sub(EinsteinLL.s[ab], _8piTLL.s[ab]);
	}
	*/
}

// gradient descent , but I'm moving everything into this file so the display shader can see it...

real4s4x4s4 calc_partial_gLL_of_8piTLL(
	<?=TPrim_t?> const TPrim,
	real4s4 const gLL,
	real4s4 const gUU
) {
	real4s4x4s4 partial_gLL_of_8piTLL = real4s4x4s4_zero;
<?
if solver.body.useEM then ?>
	real4 const EU = real3_to_real4(TPrim.E);
	real4 const EL = real4s4_real4_mul(gLL, EU);
	real const ESq = real4_dot(EL, EU);

	real4 const BU = real3_to_real4(TPrim.B);
	real4 const BL = real4s4_real4_mul(gLL, BU);
	real const BSq = real4_dot(BL, BU);

	real const sqrt_det_g = sqrt(fabs(real4s4_det(gLL)));
	real3 const SL = real3_real_mul(real3_cross(TPrim.E, TPrim.B), sqrt_det_g);

	for (int e = 0; e < stDim; ++e) {
		for (int f = e; f < stDim; ++f) {
			int const ef = sym4[e][f];
			if (e > 0 && f > 0) {
				partial_gLL_of_8piTLL.s[ef].s[sym4[0][0]] += TPrim.E.s[e-1] * TPrim.E.s[f-1] + TPrim.B.s[e-1] * TPrim.B.s[f-1];
			}
			for (int i = 0; i < sDim; ++i) {
				partial_gLL_of_8piTLL.s[ef].s[sym4[0][i+1]] -= SL.s[i] * gUU.s[ef];
				for (int j = i; j < sDim; ++j) {

					real sum = 0;
					if (e == i+1 && f == j+1) sum += ESq + BSq;

					if (e > 0 && f > 0) {
						sum += gLL.s[sym4[i+1][j+1]] * (TPrim.E.s[e-1] * TPrim.E.s[f-1] + TPrim.B.s[e-1] * TPrim.B.s[f-1]);
					}
					if (e == i+1) {
						sum -= 2. * (EU.s[f] * EL.s[j+1] + BU.s[f] * BL.s[j+1]);
					}
					if (e == j+1) {
						sum -= 2. * (EU.s[f] * EL.s[i+1] + BU.s[f] * BL.s[i+1]);
					}

					partial_gLL_of_8piTLL.s[ef].s[sym4[i+1][j+1]] += sum;
				}
			}
		}
	}
<?
end

if solver.body.useMatter then
	if solver.body.useVel then -- if we're using velocity ...
?>
	//set vU.t = 0 so we only lower by the spatial component of the metric.  right?
	real4 const vU = real3_to_real4(TPrim.v);
	real4 const vL = real4s4_real4_mul(gLL, vU);
	real const vLenSq = real4_dot(vL, vU);	//vU.t = 0 so we'll neglect the vL.t component
	real const W = 1. / sqrt(1. - sqrt(vLenSq));
	//real4 const uU = (real4){.s={W, W * vU.s1, W * vU.s2, W * vU.s3}};
	real4 const uU = (real4){.s0 = W, .s1 = W * vU.s1, .s2 = W * vU.s2, .s3 = W * vU.s3}};
	real4 const uL = real4s4_real4_mul(gLL, uU);
<?
	else -- otherwise uL = gLL.s0
?>
	//real4 const uL = (real4){.s={gLL.s00, gLL.s01, gLL.s02, gLL.s03}};
	//real4 const uL = (real4){.s={gLL.s[sym4[0][0]], gLL.s[sym4[0][1]], gLL.s[sym4[0][2]], gLL.s[sym4[0][3]]}};
	real4 const uL = (real4){.s0 = gLL.s00, .s1 = gLL.s01, .s2 = gLL.s02, .s3 = gLL.s03};
<?
	end
?>
	for (int e = 0; e < stDim; ++e) {
		for (int f = e; f < stDim; ++f) {
			int const ef = sym4[e][f];
			for (int a = 0; a < stDim; ++a) {
				for (int b = a; b < stDim; ++b) {
					int const ab = sym4[a][b];
					real sum = 0;
					if (e == a) sum += uL.s[b];
					if (e == b) sum += uL.s[a];
					sum *= uL.s[f] * (TPrim.rho * (1. + TPrim.eInt) + TPrim.P);
					if (ef == ab) sum += TPrim.P;
					partial_gLL_of_8piTLL.s[ef].s[ab] += sum;
				}
			}
		}
	}
<?
end
?>
	return real4s4x4s4_real_mul(partial_gLL_of_8piTLL, 8. * M_PI);
}

static constant int4 const int4_dirs[3] = {
	(int4)(1, 0, 0, 0),
	(int4)(0, 1, 0, 0),
	(int4)(0, 0, 1, 0),
};
int4 int4_dir(int dim, int offset) {
	return int4_dirs[dim] * offset;
}

real4s4 EFE_LL_minus_half_trace_at(
	int4 const i,
	global real4s4 const * const gLLs,
	global real4s4 const * const gUUs,
	global real4s4 const * const EFEs
) {
	if (i.x <= 0 || i.y <= 0 || i.z <= 0 ||
		i.x >= size.x || i.y >= size.y || i.z >= size.z
	) {
		return real4s4_zero;	//TODO ... consider boundary conditions
	}

	int const index = indexForInt4ForSize(i, size.x, size.y, size.z);
	real4s4 const gLL = gLLs[index];	// g_ab
	real4s4 const gUU = gUUs[index];	// g^ab
	real4s4 const EFE = EFEs[index];	// G_ab - 8 π T_ab

	//lower times lower, but used for minimizing frobenius norm of EFE_ab
	// Sum_ab (G_ab - 8 π T_ab) g_ab
	real const EFE_LL_dot_gLL = real4s4_dot(EFE, gLL);

	// common term in the gradient descent:
	// (G_uv - 8 π T_uv) - 1/2 (G_ab - 8 π T_ab) g_ab g^uv
	return real4s4_mul_add(EFE, gUU, -.5 * EFE_LL_dot_gLL);
}

// this is hardcoded to g_ab = η_ab at boundaries
//TODO ... consider boundary conditions
// or TODO ... put all g^ab boundary conditions here, and use this function in the finite-difference calculations
real4s4 gUU_at(
	int4 const i,
	global real4s4 const * const gUUs
) {
	if (i.x <= 0 || i.y <= 0 || i.z <= 0 ||
		i.x >= size.x || i.y >= size.y || i.z >= size.z
	) {
		return real4s4_Minkowski;
	}
	int const index = indexForInt4ForSize(i, size.x, size.y, size.z);
	return gUUs[index];
}

//GammaULL.a.b.c := Γ^a_bc
real4x4s4 GammaULL_at(
	int4 const i,
	global real4x4s4 const * const GammaULLs
) {
	if (i.x <= 0 || i.y <= 0 || i.z <= 0 ||
		i.x >= size.x || i.y >= size.y || i.z >= size.z
	) {
		return real4x4s4_zero;
	}

	int const index = indexForInt4ForSize(i, size.x, size.y, size.z);
	return GammaULLs[index];
}

//GammaUUL.a.b.c := Γ^ab_c = Γ^a_dc g^db
real4x4x4 GammaUUL_at(
	int4 const i,
	global real4s4 const * const gUUs,
	global real4x4s4 const * const GammaULLs
) {
	if (i.x <= 0 || i.y <= 0 || i.z <= 0 ||
		i.x >= size.x || i.y >= size.y || i.z >= size.z
	) {
		return real4x4x4_zero;
	}

	int const index = indexForInt4ForSize(i, size.x, size.y, size.z);
	return real4x4s4_real4s4_mul21(GammaULLs[index], gUUs[index]);
}

real4s4 calc_partial_gLL_of_Phi(
	int4 const i,
	global <?=TPrim_t?> const * const TPrims,
	global real4s4 const * const gLLs,
	global real4s4 const * const gUUs,
	global real4x4s4 const * const GammaULLs,
	global real4s4 const * const EFEs
) {
	int const index = indexForInt4ForSize(i, size.x, size.y, size.z);

	real4s4 const gLL = gLLs[index];
	real4s4 const gUU = gUUs[index];
	real4x4s4 const GammaULL = GammaULLs[index];

	//partial_xU2_of_gLL.cd.ab := g_ab,cd
<?= solver:finiteDifference2{
	bufferName = "gLLs",
	srcType = "4s4",
	resultName = "partial_xU2_of_gLL",
	getBoundary = function(args) return "real4s4_Minkowski" end,
} ?>

	//TODO or store this?
	//GammaLLL.a.bc := Γ_abc = g_au Γ^u_bc
	real4x4s4 const GammaLLL = real4s4_real4x4s4_mul(gLL, GammaULL);

	//partial_xU2_of_gLL_asym.a.b.c.d := g_ad,bc - g_bd,ac - g_ac,bd + g_bc,ad
	// TODO antisymmetric storage
	//both are T_abcd = -T_bacd = -T_abdc and T_abcd = T_cdab
	real4x4x4x4 partial_xU2_of_gLL_asym;
	for (int a = 0; a < stDim; ++a) {
		for (int b = 0; b < stDim; ++b) {
			for (int c = 0; c < stDim; ++c) {
				for (int d = 0; d < stDim; ++d) {
					partial_xU2_of_gLL_asym.s[a].s[b].s[c].s[d] =
						partial_xU2_of_gLL.s[sym4[a][d]].s[sym4[b][c]]
						+ partial_xU2_of_gLL.s[sym4[b][c]].s[sym4[a][d]]
						- partial_xU2_of_gLL.s[sym4[b][d]].s[sym4[a][c]]
						- partial_xU2_of_gLL.s[sym4[a][c]].s[sym4[b][d]];
				}
			}
		}
	}

	//GammaSq_asym_LLLL.a.b.c.d := Γ^e_ad Γ_ebc - Γ^e_ac Γ_ebd
	// not the same as the popular 2 Γ^a_e[c Γ^e_d]b used for Riemann with 2 ∂/dx^[c Γ^a_d]b
	// instead this is the one used with 2 g_[a[b,c]d]
	// TODO antisymmetric storage
	real4x4x4x4 GammaSq_asym_LLLL;
	for (int a = 0; a < stDim; ++a) {
		for (int b = 0; b < stDim; ++b) {
			for (int c = 0; c < stDim; ++c) {
				for (int d = 0; d < stDim; ++d) {
					real sum = 0.;
					for (int e = 0; e < stDim; ++e) {
						sum += GammaULL.s[e].s[sym4[a][d]] * GammaLLL.s[e].s[sym4[b][c]]
							- GammaULL.s[e].s[sym4[a][c]] * GammaLLL.s[e].s[sym4[b][d]];
					}
					GammaSq_asym_LLLL.s[a].s[b].s[c].s[d] = sum;
				}
			}
		}
	}

	// linear transform, not necessarily tensoral (since neither is anyways)
	real4x4x4x4 const gUU_times_partial_xU2_gLL_asym = real4x4x4x4_real4s4_mul_1_1(partial_xU2_of_gLL_asym, gUU);
	real4x4x4x4 const GammaSq_asym_ULLL = real4x4x4x4_real4s4_mul_1_1(GammaSq_asym_LLLL, gUU);

	//RiemannLLLL.a.b.c.d := R_abcd
	//= 1/2 (g_ad,bc - g_bd,ac - g_ac,bd + g_bc,ad) + g^fg (Γ_fad Γ_gbc - Γ_fac Γ_gbd)
	//= 1/2 (g_ad,bc - g_bd,ac - g_ac,bd + g_bc,ad) + Γ^e_ad Γ_ebc - Γ^e_ac Γ_ebd
	// TODO antisymmetric storage
	// but this requires inserting -1's for reading/writing ...

	//RiemannULLL.a.b.c.d := R^a_bcd = 1/2 g^ae ((g_ed,cb - g_bd,ce - g_ec,bd + g_bc,de) + g^fg (Γ_fed Γ_gbc - Γ_fec Γ_gbd))
	//TODO antisymmetric storage
	real4x4x4x4 const RiemannULLL = real4x4x4x4_mul_add(GammaSq_asym_ULLL, gUU_times_partial_xU2_gLL_asym, .5);

	//RiemannULUL.a.b.c.d = R^a_b^c_d = R^a_bud g^uc
	real4x4x4x4 const RiemannULUL = real4x4x4x4_real4s4_mul_3_1(RiemannULLL, gUU);

	//RicciLL.ab := R_ab = R^c_acb
	real4s4 const RicciLL = real4x4x4x4_tr13_to_real4s4(RiemannULLL);

	//RicciUL.a.b := R^a_b = g^ac R_cb
	real4x4 const RicciUL = real4s4_real4s4_mul(gUU, RicciLL);

	//RicciUU.ab := R^ab = R^a_c g^cb
	real4s4 const RicciUU = real4x4_real4s4_to_real4s4_mul(RicciUL, gUU);

	//Gaussian := R = R^a_a
	real const Gaussian = real4x4_tr(RicciUL);

	real4s4x4s4 const partial_gLL_of_8piTLL = calc_partial_gLL_of_8piTLL(TPrims[index], gLL, gUU);

	real4s4 const EFE = EFEs[index];	// G_ab - 8 π T_ab

	//lower times lower, but used for minimizing frobenius norm of EFE_ab
	// Sum_ab (G_ab - 8 π T_ab) g_ab
	real const EFE_LL_dot_gLL = real4s4_dot(EFE, gLL);

	// common term in the gradient descent:
	// (G_uv - 8 π T_uv) - 1/2 (G_ab - 8 π T_ab) g_ab g^uv
	real4s4 EFE_LL_minus_half_trace = real4s4_mul_add(EFE, gUU, -.5 * EFE_LL_dot_gLL);

	//GammaUUL.a.b.c := Γ^ab_c = Γ^a_dc g^db
	real4x4x4 const GammaUUL = real4x4s4_real4s4_mul21(GammaULL, gUU);

	//Gamma23U.a := Γ^a = Γ^au_u
	real4 const Gamma23U = real4x4x4_tr23(GammaUUL);

	// GammaSq_tr_2_2.pq.uv := Γ^pc_v Γ^q_cu - Γ^pc_c Γ^q_uv
	// symmetries?
	real4x4x4x4 GammaSq_tr_2_2;
	for (int p = 0; p < stDim; ++p) {
		for (int q = 0; q < stDim; ++q) {
			for (int u = 0; u < stDim; ++u) {
				for (int v = 0; v < stDim; ++v) {
					real sum = - Gamma23U.s[p] * GammaULL.s[q].s[sym4[u][v]];
					for (int c = 0; c < stDim; ++c) {
						sum += GammaUUL.s[p].s[c].s[v] * GammaULL.s[q].s[sym4[c][u]];
					}
					GammaSq_tr_2_2.s[p].s[q].s[u].s[v] = sum;
				}
			}
		}
	}

	//partial_gLL_of_Phi.pq := ∂Φ/∂g_pq
	real4s4 partial_gLL_of_Phi;
	// first calculate non-partial non ∂/∂g_pq terms:
	for (int p = 0; p < stDim; ++p) {
		for (int q = p; q < stDim; ++q) {
			int const pq = sym4[p][q];
			real sum = 0;
#if 1
			// -(G_pq - 8 π T_pq) R
			sum -= EFE.s[pq] * Gaussian;
#endif
#if 1
			// + Sum_ab EFE_ab 1/2 g_ab R^pq
			sum += EFE_LL_dot_gLL * .5 * RicciUU.s[pq];
#endif
#if 1
			for (int a = 0; a < stDim; ++a) {
				for (int b = 0; b < stDim; ++b) {
					int const ab = sym4[a][b];
					// - Sum_ab EFE_ab 8 π dT_ab/dg_pq
					sum -= EFE.s[ab] * partial_gLL_of_8piTLL.s[pq].s[ab];
				}
			}
#endif
#if 1
			for (int u = 0; u < stDim; ++u) {
				for (int v = 0; v < stDim; ++v) {
					int const uv = sym4[u][v];
					sum -= EFE_LL_minus_half_trace.s[uv] * (
					// - (EFE_uv - 1/2 EFE_ab g_ab g^uv) R^p_u^q_v
						RiemannULUL.s[p].s[u].s[q].s[v]
					// - (EFE_uv - 1/2 EFE_ab g_ab g^uv) (Γ^pc_v Γ^q_cu - Γ^pc_c Γ^q_uv)
						+ GammaSq_tr_2_2.s[p].s[q].s[u].s[v]
					);
				}
			}
#endif
#if 0
			// next calculate first-derivatives of ∂/∂g_pq terms:
			for (int dim = 0; dim < sDim; ++dim) {
				for (int offset = -<?=solver.diffOrder/2?>; offset <= <?=solver.diffOrder/2?>; ++offset) {
					if (offset == 0) continue;
					int4 const iofs = i + int4_dir(dim, offset);
					// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) g^cd (δ^f_c δ^g_v Γ^e_ud + δ^f_u δ^g_d Γ^e_cv - δ^f_c δ^g_d Γ^e_uv - δ6f_u δ^g_v Γ^e_cd) 1/2 (∂/∂g_pq(x') (D_g[g_ef] + D_f[g_eg] - D_e[g_fg]) = δ^p_e δ^q_f d1coeff_g + δ^p_e δ^q_g d1coeff_f - δ^p_f δ^q_g d1coeff_e)
					// ... expanding all 12 terms with symmath and simplifying...
					/*
					asking symmath to simplify this for me:
					let H^uv = (G_ab - 8 π T_ab) (δ^u_a δ^v_b - 1/2 g_ab g^uv) = EFE_uv - 1/2 EFE_ab g_ab g^uv
						> delta = Tensor.deltaSymbol()
						> g = Tensor.metricSymbol()
						> (H'^uv' * g'^cd' * (delta'^f_c' * delta'^g_v' * Gamma'^e_ud' + delta'^f_u' * delta'^g_d' * Gamma'^e_cv' - delta'^f_c' * delta'^g_d' * Gamma'^e_uv' - delta'^f_u' * delta'^g_v' * Gamma'^e_cd') * (delta'^p_e' * delta'^q_f' * phi'_g' + delta'^p_e' * delta'^q_g' * phi'_f' - delta'^p_f' * delta'^q_g' * phi'_e'))
							:simplifyMetrics()
							:symmetrizeIndexes(H, {1,2})()
							:symmetrizeIndexes(Gamma, {2,3})()
							:symmetrizeIndexes(g, {1,2})()
							:tidyIndexes()()
							:symmetrizeIndexes(H, {1,2})()
							:favorTensorVariance(phi'_a')()
					...gives...
							  φ_a * H^pq * Γ^ab_b
						-     φ_a * H^pb * Γ^aq_b
						-     φ_a * H^bq * Γ^ap_b
						+     φ_a * H^bc * Γ^a_bc * g^pq
						- 2 * φ_a * H^aq * Γ^pb_b
						+ 2 * φ_a * H^ab * Γ^pq_b
						+ 2 * φ_a * H^bq * Γ^pa_b
						- 2 * φ_a * H^bc * Γ^p_cb * g^aq
					*/
					int const a = dim+1;
					for (int b = 0; b < stDim; ++b) {
						sum += .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[pq] * GammaUUL_at(iofs, gUUs, GammaULLs).s[a].s[b].s[b] * d1coeff_for_offset(offset) * inv_dx.s[dim];
						sum -= .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[p][b]] * GammaUUL_at(iofs, gUUs, GammaULLs).s[a].s[q].s[b] * d1coeff_for_offset(offset) * inv_dx.s[dim];
						sum -= .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[b][q]] * GammaUUL_at(iofs, gUUs, GammaULLs).s[a].s[p].s[b] * d1coeff_for_offset(offset) * inv_dx.s[dim];
						for (int c = 0; c < stDim; ++c) {
							sum += .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[b][c]] * GammaULL_at(iofs, GammaULLs).s[a].s[sym4[b][c]] * gUU_at(iofs, gUUs).s[pq] * d1coeff_for_offset(offset) * inv_dx.s[dim];
						}
						sum -= EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[a][q]] * GammaUUL_at(iofs, gUUs, GammaULLs).s[p].s[b].s[b] * d1coeff_for_offset(offset) * inv_dx.s[dim];
						sum += EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[a][b]] * GammaUUL_at(iofs, gUUs, GammaULLs).s[p].s[q].s[b] * d1coeff_for_offset(offset) * inv_dx.s[dim];
						sum += EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[b][q]] * GammaUUL_at(iofs, gUUs, GammaULLs).s[p].s[a].s[b] * d1coeff_for_offset(offset) * inv_dx.s[dim];
						for (int c = 0; c < stDim; ++c) {
							sum -= EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[b][c]] * GammaULL_at(iofs, GammaULLs).s[p].s[sym4[c][b]] * gUU_at(iofs, gUUs).s[sym4[a][q]] * d1coeff_for_offset(offset) * inv_dx.s[dim];
						}
					}
				}
			}
#endif
#if 0
			// next calculate second-derivatives of ∂/∂g_pq terms:
			for (int dim1 = 0; dim1 < sDim; ++dim1) {
				for (int dim2 = 0; dim2 < sDim; ++dim2) {
					if (dim1 == dim2) {
						for (int offset = -<?=solver.diffOrder/2?>; offset <= <?=solver.diffOrder/2?>; ++offset) {
							int4 const iofs = i + int4_dir(dim1, offset);
							{
								// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^cd (∂/∂g_pq(x') D^2_ud[g_cv] = δ^p_c δ^q_v D^2 d2coeff_ud)
								// = + (EFE_uq - 1/2 EFE_ab g_ab g^uq) 1/2 g^pd (d2coeff_ud) |x'=x + dx^u=dx^d
								int const u = dim1+1;
								int const d = u;
								sum += .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[u][q]]
									* gUU_at(iofs, gUUs).s[sym4[p][d]]
									* (d2coeffs[abs(offset)] * inv_dx.s[dim1]);
							}{
								// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^cd (∂/∂g_pq(x') D^2_cv[g_ud] = δ^p_u δ^q_d D^2 d2coeff_cv)
								// = + (EFE_pv - 1/2 EFE_ab g_ab g^pv) 1/2 g^cq (d2coeff_cv) |x'=x + dx^v=dx^c
								int const v = dim1+1;
								int const c = v;
								sum += .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[p][v]]
									* gUU_at(iofs, gUUs).s[sym4[c][q]]
									* (d2coeffs[abs(offset)] * inv_dx.s[dim1]);
							}{
								// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^cd (∂/∂g_pq(x') D^2_cd[g_uv] = δ^p_u δ^q_v D^2 d2coeff_cd)
								// = + (EFE_pq - 1/2 EFE_ab g_ab g^pq) 1/2 g^cd (d2coeff_cd) |x'=x + dx^c=dx^d
								int const c = dim1+1;
								int const d = c;
								sum -= .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[pq]
									* gUU_at(iofs, gUUs).s[sym4[c][d]]
									* (d2coeffs[abs(offset)] * inv_dx.s[dim1]);
							}{
								// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^cd (∂/∂g_pq(x') D^2_uv[g_cd] = δ^p_c δ^q_d D^2 d2coeff_uv)
								// = + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^pq (d2coeff_uv) |x'=x + dx^u=dx^v
								int const u = dim1+1;
								int const v = u;
								sum -= .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[u][v]]
									* gUU_at(iofs, gUUs).s[pq]
									* (d2coeffs[abs(offset)] * inv_dx.s[dim1]);
							}
						}
					} else {
						for (int offset1 = -<?=solver.diffOrder/2?>; offset1 <= <?=solver.diffOrder/2?>; ++offset1) {
							if (offset1 == 0) continue;
							for (int offset2 = -<?=solver.diffOrder/2?>; offset2 <= <?=solver.diffOrder/2?>; ++offset2) {
								if (offset2 == 0) continue;
								int4 const iofs = i + int4_dir(dim1, offset1) + int4_dir(dim2, offset2);
								{
									// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^cd (∂/∂g_pq(x') D^2_ud[g_cv] = δ^p_c δ^q_v D^2 d1coeff_u d1coeff_d)
									// = + (EFE_uq - 1/2 EFE_ab g_ab g^uq) 1/2 g^pd (d1coeff_u d1coeff_d) |x'=x + dx^u + dx^d
									int const u = dim1+1;
									int const d = dim2+1;
									int const offset_u = offset1;
									int const offset_d = offset2;
									sum += .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[u][q]]
										* gUU_at(iofs, gUUs).s[sym4[p][d]]
										* (d1coeff_for_offset(offset_u) * d1coeff_for_offset(offset_d) * inv_dx.s[dim1] * inv_dx.s[dim2]);
								}{
									// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^cd (∂/∂g_pq(x') D^2_cv[g_ud] = δ^p_u δ^q_d D^2 d1coeff_c d1coeff_v)
									// = + (EFE_pv - 1/2 EFE_ab g_ab g^pv) 1/2 g^cq (d1coeff_c d1coeff_v) |x'=x + dx^v + dx^c
									int const v = dim1+1;
									int const c = dim2+1;
									int const offset_v = offset1;
									int const offset_c = offset2;
									sum += .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[p][v]]
										* gUU_at(iofs, gUUs).s[sym4[c][q]]
										* (d1coeff_for_offset(offset_v) * d1coeff_for_offset(offset_c) * inv_dx.s[dim1] * inv_dx.s[dim2]);
								}{
									// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^cd (∂/∂g_pq(x') D^2_cd[g_uv] = δ^p_u δ^q_v D^2 d1coeff_c d1coeff_d)
									// = + (EFE_pq - 1/2 EFE_ab g_ab g^pq) 1/2 g^cd (d1coeff_c d1coeff_d) |x'=x + dx^c + dx^d
									int const c = dim1+1;
									int const d = dim2+1;
									int const offset_c = offset1;
									int const offset_d = offset2;
									sum -= .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[pq]
										* gUU_at(iofs, gUUs).s[sym4[c][d]]
										* (d1coeff_for_offset(offset_c) * d1coeff_for_offset(offset_d) * inv_dx.s[dim1] * inv_dx.s[dim2]);
								}{
									// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^cd (∂/∂g_pq(x') D^2_uv[g_cd] = δ^p_c δ^q_d D^2 d1coeff_u d1coeff_v)
									// = + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^pq (d1coeff_u d1coeff_v) |x'=x + dx^u + dx^v
									int const u = dim1+1;
									int const v = dim2+1;
									int const offset_u = offset1;
									int const offset_v = offset2;
									sum -= .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[u][v]]
										* gUU_at(iofs, gUUs).s[pq]
										* (d1coeff_for_offset(offset_u) * d1coeff_for_offset(offset_v) * inv_dx.s[dim1] * inv_dx.s[dim2]);
								}
							}
						}
					}
				}
			}
#endif
			partial_gLL_of_Phi.s[pq] = sum;
		}
	}

	return partial_gLL_of_Phi;
}

_Static_assert(sizeof(real4s4) == sizeof(real) * 10, "here");
_Static_assert(sizeof(gPrim_t) == sizeof(real) * 10, "here");

gPrim_t calc_partial_gPrim_of_Phi(
	int4 const i,
	global <?=TPrim_t?> const * const TPrims,
	global gPrim_t const * const gPrims,
	global real4s4 const * const gLLs,
	global real4s4 const * const gUUs,
	global real4x4s4 const * const GammaULLs,
	global real4s4 const * const EFEs
) {
	int const index = indexForInt4ForSize(i, size.x, size.y, size.z);

	real4s4 const partial_gLL_of_Phi = calc_partial_gLL_of_Phi(i, TPrims, gLLs, gUUs, GammaULLs, EFEs);

	gPrim_t const gPrim = gPrims[index];
	real3s3 const gammaLL = gPrim.gammaLL;
	real3 const betaU = gPrim.betaU;
	real3 const betaL = real3s3_real3_mul(gammaLL, betaU);

	gPrim_t partial_gPrim_of_Phi;

	partial_gPrim_of_Phi.alpha = -2. * gPrim.alpha * partial_gLL_of_Phi.s00;
	for (int m = 0; m < sDim; ++m) {
		partial_gPrim_of_Phi.betaU.s[m] = 2. * partial_gLL_of_Phi.s00 * betaL.s[m];
		for (int n = 0; n < sDim; ++n) {
			partial_gPrim_of_Phi.betaU.s[m] += partial_gLL_of_Phi.s[sym4[0][n+1]] * gammaLL.s[sym3[n][m]];
		}
	}
	for (int m = 0; m < sDim; ++m) {
		for (int n = m; n < sDim; ++n) {
			int const mn = sym3[m][n];
			partial_gPrim_of_Phi.gammaLL.s[mn] =
				partial_gLL_of_Phi.s[sym4[m+1][n+1]]
				+ betaU.s[m] * (
					  partial_gLL_of_Phi.s[sym4[0][0]] * betaU.s[n]
					+ 2. * partial_gLL_of_Phi.s[sym4[0][n+1]]
				);
		}
	}

	return partial_gPrim_of_Phi;
}
