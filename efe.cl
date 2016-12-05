#define _real3(a,b,c) (real3){.s={a,b,c}}

constant real c = 299792458;			// m/s 
constant real G = 6.67384e-11;		// m^3 / (kg s^2)

static inline real real3_dot(real3 a, real3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline real3 real3_cross(real3 a, real3 b) {
	return _real3(
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x);
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

constant sym3 sym3_ident = {
<? for i=0,2 do
	for j=i,2 do ?>
		.s<?=i?><?=j?> = <?= i==j and 1 or 0 ?>,
<?	end
end ?>
};

static inline real sym3_det(sym3 m) {
	return m.s00 * m.s11 * m.s22
		+ m.s01 * m.s12 * m.s02
		+ m.s02 * m.s01 * m.s12
		- m.s02 * m.s11 * m.s02
		- m.s01 * m.s01 * m.s22
		- m.s00 * m.s12 * m.s12;
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
	return (real3){
		.x = m.xx * v.x + m.xy * v.y + m.xz * v.z,
		.y = m.xy * v.y + m.yy * v.y + m.yz * v.z,
		.z = m.xz * v.z + m.yz * v.y + m.zz * v.z,
	};
}

static inline sym4 sym4_scale(sym4 a, real s) {
	return (sym4){
<?
for a=0,3 do
	for b=a,3 do ?>
		.s<?=a?><?=b?> = a.s<?=a?><?=b?> * s,
<?	end
end ?>
	};
}

static inline sym4 sym4_sub(sym4 a, sym4 b) {
	return (sym4){
<?
for a=0,3 do
	for b=a,3 do ?>
		.s<?=a?><?=b?> = a.s<?=a?><?=b?> - b.s<?=a?><?=b?>,
<?	end
end ?>		
	};
}

static inline real4 sym4_real4_mul(sym4 m, real4 v) {
	return (real4){
<? for i=0,3 do ?>
		0.
	<? for j=0,3 do ?>
			+ m.s<?=sym(i,j)?> * v.s<?=j?>
	<? end ?>,
<? end ?>
	};
}

static inline real sym4_det(sym4 m) {
	return
		m.s03 * m.s12 * m.s12 * m.s03 - m.s02 * m.s13 * m.s12 * m.s03 -
		m.s03 * m.s11 * m.s22 * m.s03 + m.s01 * m.s13 * m.s22 * m.s03 +
		m.s02 * m.s11 * m.s23 * m.s03 - m.s01 * m.s12 * m.s23 * m.s03 -
		m.s03 * m.s12 * m.s02 * m.s13 + m.s02 * m.s13 * m.s02 * m.s13 +
		m.s03 * m.s01 * m.s22 * m.s13 - m.s00 * m.s13 * m.s22 * m.s13 -
		m.s02 * m.s01 * m.s23 * m.s13 + m.s00 * m.s12 * m.s23 * m.s13 +
		m.s03 * m.s11 * m.s02 * m.s23 - m.s01 * m.s13 * m.s02 * m.s23 -
		m.s03 * m.s01 * m.s12 * m.s23 + m.s00 * m.s13 * m.s12 * m.s23 +
		m.s01 * m.s01 * m.s23 * m.s23 - m.s00 * m.s11 * m.s23 * m.s23 -
		m.s02 * m.s11 * m.s02 * m.s33 + m.s01 * m.s12 * m.s02 * m.s33 +
		m.s02 * m.s01 * m.s12 * m.s33 - m.s00 * m.s12 * m.s12 * m.s33 -
		m.s01 * m.s01 * m.s22 * m.s33 + m.s00 * m.s11 * m.s22 * m.s33;
}

constant tensor_4sym4 tensor_4sym4_zero = (tensor_4sym4){
<? for a=0,3 do ?>
	.s<?=a?> = (sym4){
<?	for b=0,3 do
		for c=b,3 do ?>
		.s<?=b?><?=c?> = 0.,
<?		end
	end ?>
	},
<? end ?>
};

static inline tensor_4sym4 tensor_4sym4_scale(tensor_4sym4 a, real s) {
	return (tensor_4sym4){
<? for a=0,3 do ?>
		.s<?=a?> = sym4_scale(a.s<?=a?>, s),
<? end ?>	
	};
}

static inline tensor_4sym4 tensor_4sym4_sub(tensor_4sym4 a, tensor_4sym4 b) {
	return (tensor_4sym4){
<? for a=0,3 do ?>
		.s<?=a?> = sym4_sub(a.s<?=a?>, b.s<?=a?>),
<? end ?>
	};
}

constant const int dim = <?=dim?>;	
constant const int subDim = <?=subDim?>;
constant const int gridDim = <?=gridDim?>;

constant const int4 size = (int4)(<?=clnumber(size.x)?>, <?=clnumber(size.y)?>, <?=clnumber(size.z)?>, 0);
constant const int4 stepsize = (int4)(1, <?=size.x?>, <?=size.x * size.y?>, <?=size.x * size.y * size.z?>);

#define globalInt4()	(int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0)

#define indexForInt4(i) (i.x + size.x * (i.y + size.y * i.z))

#define INIT_KERNEL() \
	int4 i = globalInt4(); \
	if (i.x >= size.x || i.y >= size.y || i.z >= size.z) return; \
	int index = indexForInt4(i);

constant real3 xmin = _real3(<?=xmin.x?>, <?=xmin.y?>, <?=xmin.z?>);
constant real3 xmax = _real3(<?=xmax.x?>, <?=xmax.y?>, <?=xmax.z?>);
constant real3 dx = _real3(
	<?=tonumber(xmax.x - xmin.x) / tonumber(size.x)?>,
	<?=tonumber(xmax.x - xmin.x) / tonumber(size.x)?>,
	<?=tonumber(xmax.x - xmin.x) / tonumber(size.x)?>);

constant real3 inv_dx = _real3(
	<?=tonumber(size.x) / tonumber(xmax.x - xmin.x)?>,
	<?=tonumber(size.x) / tonumber(xmax.x - xmin.x)?>,
	<?=tonumber(size.x) / tonumber(xmax.x - xmin.x)?>);


#define getX(i) _real3( \
	xmin.x + ((real)i.x + .5)/(real)size.x * (xmax.x - xmin.x),	\
	xmin.y + ((real)i.y + .5)/(real)size.y * (xmax.y - xmin.y),	\
	xmin.z + ((real)i.z + .5)/(real)size.z * (xmax.z - xmin.z));

void calc_EinsteinLL(
	sym4* result,
	global const sym4* gLLs,
	global const sym4* gUUs,
	global const tensor_4sym4* GammaULLs
) {
	int4 i = globalInt4();
	int index = indexForInt4(i);
	
	//here's where the finite difference stuff comes in ...
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

	global const tensor_4sym4* GammaULL = GammaULLs + index;
	real4 Gamma12L = (real4){
<? for a=0,dim-1 do ?>
		0. 
	<? for b=0,dim-1 do ?>
			+ GammaULL->s<?=b?>.s<?=sym(b,a)?>
	<? end ?>
	<? if a < dim-1 then ?>,<? end ?>
<? end ?>};

	sym4 RicciLL = {
<? for a=0,dim-1 do ?>
	<? for b=a,dim-1 do ?>
		.s<?=a?><?=b?> = 0.
		<? for c=0,dim-1 do ?>
			+ dGammaLULL.s<?=c?>.s<?=c?>.s<?=a..b?> 
			- dGammaLULL.s<?=b?>.s<?=c?>.s<?=sym(c,a)?> 
			+ Gamma12L.s<?=c?> * GammaULL->s<?=c?>.s<?=a..b?>
			<? for d=0,dim-1 do ?>
			- GammaULL->s<?=c?>.s<?=sym(d,b)?> * GammaULL->s<?=d?>.s<?=sym(c,a)?> 
			<? end ?>
		<? end ?>,
	<? end ?>
<? end ?>};

	global const sym4* gUU = gUUs + index;
	real Gaussian = 0.
<? for a=0,dim-1 do ?>
	<? for b=0,dim-1 do ?>
		+ gUU->s<?=sym(a,b)?> * RicciLL.s<?=sym(a,b)?>
	<? end ?>
<? end ?>;

	global const sym4* gLL = gLLs + index;
	
<? for a=0,dim-1 do ?>
	<? for b=a,dim-1 do ?>
	result->s<?=a..b?> = RicciLL.s<?=a..b?> - .5 * gLL->s<?=a..b?> * Gaussian;
	<? end ?>
<? end ?>
}

void calc_8piTLL(
	sym4* result,
	global const sym4* gLL,
	global const TPrim_t* TPrim
) {
	
	//if we're using EM ...
	
	/*
	assume the E and B fields are upper 3-vectors
	T_ab = F_au F_b^u - 1/4 g_ab F_uv F^uv
	*/

	real4 EU = (real4)(0 <?for i=0,2 do ?>, TPrim->E.s<?=i?> <? end ?>); 
	real4 EL = sym4_real4_mul(*gLL, EU);
	real ESq = dot(EL, EU);
	
	real4 BU = (real4)(0 <?for i=0,2 do ?>, TPrim->B.s<?=i?> <? end ?>); 
	real4 BL = sym4_real4_mul(*gLL, BU);
	real BSq = dot(BL, BU);

	real sqrt_det_g = sqrt(fabs(sym4_det(*gLL)));
	real3 SL = real3_scale(real3_cross(TPrim->E, TPrim->B), sqrt_det_g);
	
	sym4 _8piT_EM_LL = {
		.s00 = ESq + BSq,
<? for i=0,subDim-1 do ?>
		.s0<?=i+1?> = -2. * SL.s<?=i?>,
	<? for j=i,subDim-1 do ?>
		.s<?=i+1?><?=j+1?> = gLL->s<?=i?><?=j?> * (ESq + BSq)
			- 2. * (
				EL.s<?=i+1?> * EL.s<?=j+1?>
				+ BL.s<?=i+1?> * BL.s<?=j+1?>
			),
	<? end ?>
<? end ?>};

	//if we're using matter ...

		//if we're using velocity ...

		//otherwise uL = gLL->s0
	real4 uL = (real4)(gLL->s00, gLL->s01, gLL->s02, gLL->s03);

	sym4 _8piT_matter_LL = (sym4){
<? for a=0,dim-1 do ?>
	<? for b=a,dim-1 do ?>
		.s<?=a..b?> = 8. * M_PI * (
			uL.s<?=a?> * uL.s<?=b?> * (TPrim->rho * (1. + TPrim->eInt) + TPrim->P)
			+ gLL->s<?=a..b?> * TPrim->P),
	<? end ?>
<? end ?>
	};

	<? for a=0,dim-1 do ?>
		<? for b=a,dim-1 do ?>
	result->s<?=a..b?> = _8piT_EM_LL.s<?=a..b?> + _8piT_matter_LL.s<?=a..b?>;
		<? end ?>
	<? end ?>
}

// init

kernel void init_gPrims(
	global gPrim_t* gPrims
) {
	INIT_KERNEL();
	
	real3 x = getX(i);
	real r = real3_len(x);

	global gPrim_t* gPrim = gPrims + index;

	//init to flat by default
	*gPrim = (gPrim_t){
		.alpha = 1,
		.betaU = _real3(0,0,0),
		.gammaLL = sym3_ident,
	};

	<?=initCond.code?>
}

kernel void init_TPrims(
	global TPrim_t* TPrims
) {
	INIT_KERNEL();

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

	<?=body.init?>
}

// compute buffers to compute EFE

kernel void calc_gLLs_and_gUUs(
	global sym4* gLLs,
	global sym4* gUUs,
	const global gPrim_t* gPrims
) {
	INIT_KERNEL();

	global const gPrim_t* gPrim = gPrims + index;

	real alphaSq = gPrim->alpha * gPrim->alpha;
	real3 betaL = sym3_real3_mul(gPrim->gammaLL, gPrim->betaU);
	real betaSq = real3_dot(betaL, gPrim->betaU);

	global sym4* gLL = gLLs + index;
	gLL->s00 = -alphaSq + betaSq;
	<? for i=0,subDim-1 do ?>
	gLL->s0<?=i+1?> = betaL.s<?=i?>;
		<? for j=i,subDim-1 do ?>
	gLL->s<?=i+1?><?=j+1?> = gPrim->gammaLL.s<?=i?><?=j?>;
		<? end ?>
	<? end ?>

	real det_gammaLL = sym3_det(gPrim->gammaLL);
	sym3 gammaUU = sym3_inv(det_gammaLL, gPrim->gammaLL);

	global sym4* gUU = gUUs + index;
	gUU->s00 = -1./alphaSq;
	<? for i=0,subDim-1 do ?>
	gUU->s0<?=i+1?> = gPrim->betaU.s<?=i?> / alphaSq;
		<? for j=i,subDim-1 do ?>
	gUU->s<?=i+1?><?=j+1?> = gammaUU.s<?=i?><?=j?> - gPrim->betaU.s<?=i?> * gPrim->betaU.s<?=j?> / alphaSq;
		<? end ?>
	<? end ?>
}

kernel void calc_GammaULLs(
	global tensor_4sym4* GammaULLs,
	const global sym4* gLLs,
	const global sym4* gUUs
) {
	INIT_KERNEL();

	//here's where the finite difference stuff comes in ...
	tensor_4sym4 dgLLL;
	dgLLL.s0 = (sym4){
<? for a=0,dim-1 do
	for b=a,dim-1 do ?>
		.s<?=a..b?> = 0.,
<?	end
end ?>
	};
	<? for i=0,gridDim-1 do ?>{
		int4 iL = i;
		iL.s<?=i?> = max(i.s<?=i?> - 1, 0);
		int indexL = indexForInt4(iL);
		global const sym4* gLL_prev = gLLs + indexL;
		
		int4 iR = i;
		iR.s<?=i?> = min(i.s<?=i?> + 1, size.s<?=i?> - 1);
		int indexR = indexForInt4(iR);
		global const sym4* gLL_next = gLLs + indexR;
		
		dgLLL.s<?=i+1?> = (sym4){
		<? for a=0,dim-1 do ?>
			<? for b=a,dim-1 do ?>
			.s<?=a..b?> = (gLL_next->s<?=a..b?> - gLL_prev->s<?=a..b?>) * .5 * inv_dx.s<?=i?>,
			<? end ?>
		<? end ?>
		};
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
}

kernel void calc_EFEs(
	global sym4* EFEs,
	global const gPrim_t* gPrims,
	global const TPrim_t* TPrims,
	global const sym4* gLLs,
	global const sym4* gUUs,
	global const tensor_4sym4* GammaULLs
) {
	INIT_KERNEL();
	
	sym4 EinsteinLL;
	calc_EinsteinLL(&EinsteinLL, gLLs, gUUs, GammaULLs);
	
	sym4 _8piTLL;
	calc_8piTLL(&_8piTLL, gLLs+index, TPrims+index);

	global sym4* EFE = EFEs + index;
	<? for a=0,dim-1 do ?>
		<? for b=a,dim-1 do ?>
	EFE->s<?=a..b?> = EinsteinLL.s<?=a..b?> - _8piTLL.s<?=a..b?>;
		<? end ?>
	<? end ?>
}

// compute aux buffers for output

kernel void calc_detGammas(
	global real* detGammas,
	global const gPrim_t* gPrims
) {
	INIT_KERNEL();
	detGammas[index] = sym3_det(gPrims[index].gammaLL);
}

kernel void calc_numericalGravity(
	global real* numericalGravity,
	global const tensor_4sym4* GammaULLs
) {
	INIT_KERNEL();
	real3 x = getX(i);
	real r = real3_len(x);
	global const tensor_4sym4* GammaULL = GammaULLs + index;
	numericalGravity[index] = (0.
		+ GammaULL->s1.s00 * x.s0 / r
		+ GammaULL->s2.s00 * x.s1 / r
		+ GammaULL->s3.s00 * x.s2 / r) * c * c;
}

kernel void calc_analyticalGravity(
	global real* analyticalGravity
) {
	INIT_KERNEL();
	real3 x = getX(i);
	real r = real3_len(x);
	real matterRadius = min(r, (real)<?=body.radius?>);
	real volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
	real m = <?=body.density?> * volumeOfMatterRadius;	// m^3
	real dm_dr = 0;
	analyticalGravity[index] = (2*m * (r - 2*m) + 2 * dm_dr * r * (2*m - r)) / (2 * r * r * r)
		* c * c;	//+9 at earth surface, without matter derivatives
}

kernel void calc_norm_EFE_tts(
	global real* norm_EFE_tts,
	global const sym4* EFEs
) {
	INIT_KERNEL();
	norm_EFE_tts[index] = EFEs[index].s00 / (8. * M_PI) * c * c / G / 1000.;
}

kernel void calc_norm_EFE_tis(
	global real* norm_EFE_tis,
	global const sym4* EFEs
) {
	INIT_KERNEL();
	global const sym4* EFE = EFEs + index;	
	norm_EFE_tis[index] = sqrt(0.
<? for i=0,subDim-1 do ?>
		+ EFE->s0<?=i+1?> * EFE->s0<?=i+1?>
<? end ?>) * c;
}

kernel void calc_norm_EFE_ijs(
	global real* norm_EFE_ijs,
	global const sym4* EFEs
) {
	INIT_KERNEL();
	global const sym4* EFE = EFEs + index;
	norm_EFE_ijs[index] = sqrt(0.
<? for i=0,subDim-1 do
	for j=0,subDim-1 do ?>
		+ EFE->s<?=sym(i+1,j+1)?> * EFE->s<?=sym(i+1,j+1)?>
<?	end
end ?>);
}

kernel void calc_norm_EinsteinLLs(
	global real* norm_EinsteinLLs,
	global const sym4* gLLs,
	global const sym4* gUUs,
	global const tensor_4sym4* GammaULLs
) {
	INIT_KERNEL();
	sym4 EinsteinLL;
	calc_EinsteinLL(&EinsteinLL, gLLs, gUUs, GammaULLs);
	norm_EinsteinLLs[index] = sqrt(0.
<? for a=0,dim-1 do
	for b=0,dim-1 do ?>
		+ EinsteinLL.s<?=sym(a,b)?> * EinsteinLL.s<?=sym(a,b)?>
<?	end
end ?>);
}

// and here is the code used for updating weights based on the EFE

kernel void calc_dPhi_dgLLs(
	global sym4* dPhi_dgLLs,
	global const TPrim_t* TPrims,
	global const sym4* gLLs,
	global const sym4* gUUs,
	global const tensor_4sym4* GammaULLs,
	global const sym4* EFEs
) {
	INIT_KERNEL();
	
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

	sym4 RicciLL = (sym4){
<? for a=0,dim-1 do ?>
	<? for b=a,dim-1 do ?>
		.s<?=a..b?> = 0.
		<? for c=0,dim-1 do ?>
			+ RiemannULLL.s<?=c?>.s<?=a?>.s<?=sym(c,b)?>
		<? end ?>,
	<? end ?>
<? end ?>};

	real RicciUL[4][4] = {
<? for a=0,dim-1 do ?>
		{
	<? for b=0,dim-1 do ?>
			0.
		<? for c=0,dim-1 do ?>
				+ gUU->s<?=sym(a,c)?> * RicciLL.s<?=sym(c,b)?>
		<? end ?>,
	<? end ?>
		},
<? end ?>
	};

	sym4 RicciUU = (sym4){
<? for a=0,dim-1 do
	for b=a,dim-1 do ?>
		.s<?=a..b?> = 0.
	<?	for c=0,dim-1 do ?>
		 + RicciUL[<?=a?>][<?=c?>] * gUU->s<?=sym(c,b)?>
	<?	end ?>,
<?	end
end ?>
	};

	real Gaussian = 0.
<? for a=0,dim-1 do ?>
		+ RicciUL[<?=a?>][<?=a?>]
<? end ?>;



	real GammaUUL[4][4][4] = {
<? for a=0,dim-1 do ?>
	{<? for b=0,dim-1 do ?>
		{<? for c=0,dim-1 do ?>
			0.
			<? for d=0,dim-1 do ?>
			+ GammaULL->s<?=a?>.s<?=sym(d,c)?> * gUU->s<?=sym(d,b)?>
			<? end ?>,
		<? end ?>},
	<? end ?>},
<? end ?>
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
					+ GammaUUL[<?=e?>][<?=c?>][<?=c?>] * GammaULL->s<?=f?>.s<?=sym(b,a)?>
					- GammaUUL[<?=e?>][<?=c?>][<?=b?>] * GammaULL->s<?=f?>.s<?=sym(c,a)?>
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
					-RicciUU.s<?=e?><?=f?>
				<? for c=0,dim-1 do ?>
					<? for d=0,dim-1 do ?>
					+ gUU->s<?=sym(c,d)?> * dRLL_dgLL.s<?=e?><?=f?>.s<?=sym(c,d)?>
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

	real sqrt_det_g = sqrt(fabs(sym4_det(*gLL)));
	real3 SL = real3_scale(real3_cross(TPrim->E, TPrim->B), sqrt_det_g);

	tensor_sym4sym4 d_8piTLL_dgLL = (tensor_sym4sym4){
<? for e=0,dim-1 do ?>
	<? for f=e,dim-1 do ?>
		.s<?=e?><?=f?> = (sym4){
			//dTEM_00/dg_<?=e?><?=f?>
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
}

kernel void update_dgLLs(
	global sym4* gLLs,
	global const sym4* dPhi_dgLLs
) {
	INIT_KERNEL();
	global sym4* gLL = gLLs + index;
	global const sym4* dPhi_dgLL = dPhi_dgLLs + index;
<? for a=0,dim-1 do ?>
	<? for b=a,dim-1 do ?>
	gLL.s<?=a..b?> -= dPhi_dgLL.s<?=a..b?> * <?=updateAlpha?>;
	<? end ?>
<? end ?>
}
