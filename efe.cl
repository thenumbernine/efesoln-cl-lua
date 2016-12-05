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

constant sym4 sym4_zero = {
<?
for a=0,3 do
	for b=a,3 do ?>
	.s<?=a?><?=b?> = 0.,
<?	end
end ?>
};

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

static inline real sym4_dot(sym4 a, sym4 b) {
	return 0
<? for a=0,3 do ?>
	<? for b=a,3 do ?>
	+ a.s<?=a?><?=b?> * b.s<?=a?><?=b?>
	<? end ?>
<? end ?>;
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

#define getX(i) _real3( \
	xmin.x + ((real)i.x + .5)/(real)size.x * (xmax.x - xmin.x),	\
	xmin.y + ((real)i.y + .5)/(real)size.y * (xmax.y - xmin.y),	\
	xmin.z + ((real)i.z + .5)/(real)size.z * (xmax.z - xmin.z));

sym4 calc_EinsteinLL(
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
		iL.s<?=i?> = min(i.s<?=i?> + 1, size.s<?=i?> - 1);
		int indexL = indexForInt4(iL);
		global const tensor_4sym4* GammaULL_prev = GammaULLs + indexL;
		
		int4 iR = i;
		iR.s<?=i?> = max(i.s<?=i?> - 1, 0);
		int indexR = indexForInt4(iR);
		global const tensor_4sym4* GammaULL_next = GammaULLs + indexR;
		
		dGammaLULL.s<?=i+1?> = tensor_4sym4_scale(
			tensor_4sym4_sub(*GammaULL_next, *GammaULL_prev),
			1. / (2. * dx.s<?=i?> ));
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
			+ dGammaLULL.s<?=c?>.s<?=c?>.s<?=a..b?> - dGammaLULL.s<?=b?>.s<?=c?>.s<?=sym(a,c)?> + GammaULL->s<?=c?>.s<?=a..b?> * Gamma12L.s<?=c?>
			<? for d=0,dim-1 do ?>
				- GammaULL->s<?=d?>.s<?=sym(a,c)?> * GammaULL->s<?=c?>.s<?=sym(b,d)?>
			<? end ?>
		<? end ?>,
	<? end ?>
<? end ?>};

	global const sym4* gUU = gUUs + index;
	real Gaussian = sym4_dot(*gUU, RicciLL);

	global const sym4* gLL = gLLs + index;
	return sym4_sub(RicciLL, sym4_scale(*gLL, -.5 * Gaussian));
}

sym4 calc_8piTLL(
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

	real sqrt_det_g = sqrt(sym4_det(*gLL));
	real3 SL = real3_scale(real3_cross(TPrim->E, TPrim->B), sqrt_det_g);
	
	sym4 T_EM_LL = {
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

	sym4 T_matter_LL = (sym4){
<? for a=0,dim-1 do ?>
	<? for b=a,dim-1 do ?>
		.s<?=a..b?> = uL.s<?=a?> * uL.s<?=b?> * (TPrim->rho * (1. + TPrim->eInt) + TPrim->P)
			+ gLL->s<?=a..b?> * TPrim->P,
	<? end ?>
<? end ?>
	};

	return T_EM_LL;
}

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
}

