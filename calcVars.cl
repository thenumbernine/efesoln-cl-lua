// different initial conditions / boundary conditions of gPrims
// this is the only code that uses solver.body

gPrim_t calc_gPrim_flat(real3 const x) {
	return (gPrim_t){
		.alpha = 1,
		.betaU = real3_zero,
		.gammaLL = real3s3_ident,
	};
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
<?
for i=0,sDim-1 do
?>	gPrim.betaU.s<?=i?> = 0;
<?	for j=i,sDim-1 do
?>	gPrim.gammaLL.s<?=i..j?> = <?if i==j then?>1. + <?end?>x.s<?=i?>/r * x.s<?=j?>/r * 2*m/(r - 2*m);
<?	end
end
?>
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
<?
for i=0,sDim-1 do 
?>	gLL.s0<?=i+1?> = betaL.s<?=i?>;
<?	for j=i,sDim-1 do
?>	gLL.s<?=i+1?><?=j+1?> = gammaLL.s<?=i?><?=j?>;
<?	end
end
?>
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
<?
for i=0,sDim-1 do
?>	gUU.s0<?=i+1?> = betaU.s<?=i?> * invAlphaSq;
<?	for j=i,sDim-1 do
?>	gUU.s<?=i+1?><?=j+1?> = gammaUU.s<?=i?><?=j?> - betaU.s<?=i?> * betaU.s<?=j?> * invAlphaSq;
<?	end
end
?>
	return gUU;
}

real4s4 calc_gLL_flat() {
#if 1	
	return (real4s4){
<?
for a=0,3 do
	for b=a,3 do
?>		.s<?=a..b?> = <?=a==b and (a==0 and -1 or 1) or 0?>,
<?	end
end
?>	};
#else //doesn't work so well
	int4 i = globalInt4();
	real3 const x = getX(i);
	return calc_gLL_from_gPrim(calc_gPrim_flat(x));
#endif
}

// for _zero, the constant access is much faster than making a new struct in-place
constant real4s4 const gLL_flat = (real4s4){
<?
for a=0,3 do
	for b=a,3 do
?>		.s<?=a..b?> = <?=a==b and (a==0 and -1 or 1) or 0?>,
<?	end
end
?>	};

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

	//g_ab,c := dgLLL.c.ab
<?= solver:finiteDifference{
	bufferName = "gLLs",
	srcType = "4s4",
	resultName = "dgLLL",
	--boundaryCode = "calc_gLL_flat()",
	boundaryCode = "gLL_flat",
} ?>

	//Γ_abc := GammaLLL.a.bc
	//Γ_abc = 1/2 (g_ab,c + g_ac,b - g_bc,a)
	real4x4s4 const GammaLLL = (real4x4s4){
<? 
for a=0,stDim-1 do
	for b=0,stDim-1 do
		for c=b,stDim-1 do
?>		.s<?=a?>.s<?=b?><?=c?> = .5 * (dgLLL.s<?=c?>.s<?=sym(a,b)?> + dgLLL.s<?=b?>.s<?=sym(a,c)?> - dgLLL.s<?=a?>.s<?=sym(b,c)?>),
<?		end
	end
end
?>	};

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
= solver:finiteDifference{
	bufferName = "TPrims",
	getValue = function(index) return "TPrims["..index.."].JU" end,
	valueType = "real4",
} 
--]]
?>

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
