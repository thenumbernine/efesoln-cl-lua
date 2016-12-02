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
	return (real3){
		.x = m.xx * v.x + m.xy * v.y + m.xz * v.z,
		.y = m.xy * v.y + m.yy * v.y + m.yz * v.z,
		.z = m.xz * v.z + m.yz * v.y + m.zz * v.z,
	};
}

constant sym4 sym4_zero = (sym4){
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
<? for i=0,9 do ?>
	+ a.s[<?=i?>] * b.s[<?=i?>]
<? end ?>;
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
	gPrim_t gPrim,
	sym4 gLL,
	TPrim_t TPrim
) {
	
	/*
	assume the E and B fields are upper 3-vectors
	T_ab = F_au F_b^u - 1/4 g_ab F_uv F^uv
	*/
	real ESq = 0
<? 
for i=0,subDim-1 do 
	for j=i,subDim-1 do ?>
		+ gLL.s<?=sym(i+1,j+1)?> * TPrim.E.s<?=i?> * TPrim.E.s<?=j?>
<?	end
end ?>;
	real BSq = 0
<? 
for i=0,subDim-1 do 
	for j=i,subDim-1 do ?>
		+ gLL.s<?=sym(i+1,j+1)?> * TPrim.B.s<?=i?> * TPrim.B.s<?=j?>
<?	end
end ?>;

	real3 S = real3_cross(TPrim.E, TPrim.B);

	return (sym4){
		.s00 = ESq + BSq,
<? for i=0,subDim-1 do ?>
		.s0<?=i+1?> = 2. * S.s<?=i?>,
	<? for j=i,subDim-1 do ?>
		.s<?=i+1?><?=j+1?> = 
			<? if i == j then ?> ESq + BSq <? else ?> 0 <? end ?>
			- 2. * (
				TPrim.E.s<?=i?> * TPrim.E.s<?=j?>
				+ TPrim.B.s<?=i?> * TPrim.B.s<?=j?>
			),
	<? end ?>
<? end ?>};
}
