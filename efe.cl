#define _real3(a,b,c) (real3){.s={a,b,c}}

constant real c = 299792458;			// m/s 
constant real G = 6.67384e-11;		// m^3 / (kg s^2)

#define real3_zero ((real3){ \
	.x = 0, \
	.y = 0, \
	.z = 0, \
})

static inline real real3_dot(real3 const a, real3 const b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline real3 real3_cross(real3 const a, real3 const b) {
	return _real3(
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x);
}

static inline real real3_lenSq(real3 const a) {
	return real3_dot(a,a);
}

static inline real real3_len(real3 const a) {
	return sqrt(real3_lenSq(a));
}

static inline real3 real3_real_mul(real3 const a, real const s) {
	return _real3(a.x * s, a.y * s, a.z * s);
}

static inline real3 real3_add(real3 const a, real3 const b) {
	return _real3(a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline real3 real3_sub(real3 const a, real3 const b) {
	return _real3(a.x - b.x, a.y - b.y, a.z - b.z);
}

constant sym3 sym3_ident = {<? 
for i=0,2 do
	for j=i,2 do 
?>	.s<?=i?><?=j?> = <?= i==j and 1 or 0 ?>,<?
	end
end 
?>};

static inline real sym3_det(sym3 const m) {
	return m.s00 * m.s11 * m.s22
		+ m.s01 * m.s12 * m.s02
		+ m.s02 * m.s01 * m.s12
		- m.s02 * m.s11 * m.s02
		- m.s01 * m.s01 * m.s22
		- m.s00 * m.s12 * m.s12;
}

static inline sym3 sym3_inv(real const d, sym3 const m) {
	return (sym3){
		.xx = (m.yy * m.zz - m.yz * m.yz) / d,
		.xy = (m.xz * m.yz - m.xy * m.zz) / d,
		.xz = (m.xy * m.yz - m.xz * m.yy) / d,
		.yy = (m.xx * m.zz - m.xz * m.xz) / d,
		.yz = (m.xz * m.xy - m.xx * m.yz) / d,
		.zz = (m.xx * m.yy - m.xy * m.xy) / d,
	};
}

static inline real3 sym3_real3_mul(sym3 const m, real3 const v) {
	return (real3){
		.x = m.xx * v.x + m.xy * v.y + m.xz * v.z,
		.y = m.xy * v.y + m.yy * v.y + m.yz * v.z,
		.z = m.xz * v.z + m.yz * v.y + m.zz * v.z,
	};
}

#define sym4_zero ((sym4){ \
	.tt = 0, .tx = 0, .ty = 0, .tz = 0, \
	.xx = 0, .xy = 0, .xz = 0, \
	.yy = 0, .yz = 0, \
	.zz = 0, \
})

static inline real sym4_dot(sym4 const a, sym4 const b) {
	return 0.<?
for a=0,stDim-1 do
	for b=0,stDim-1 do 
?> + a.s<?=sym(a,b)?> * b.s<?=sym(a,b)?><?
	end 
end ?>;
}

static inline sym4 sym4_outer(real4 const v) {
	return (sym4){
<? for a=0,stDim-1 do
	for b=a,stDim-1 do
?>		.s<?=a..b?> = v.s<?=a?> * v.s<?=b?>,
<? 	end
end
?>	};
}

static inline sym4 sym4_real_mul(sym4 const a, real const s) {
	return (sym4){
<?
for a=0,3 do
	for b=a,3 do 
?>		.s<?=a?><?=b?> = a.s<?=a?><?=b?> * s,
<?	end
end
?>	};
}

static inline sym4 sym4_add(sym4 const a, sym4 const b) {
	return (sym4){
<?
for a=0,3 do
	for b=a,3 do 
?>		.s<?=a?><?=b?> = a.s<?=a?><?=b?> + b.s<?=a?><?=b?>,
<?	end
end
?>	};
}

static inline sym4 sym4_sub(sym4 const a, sym4 const b) {
	return (sym4){
<?
for a=0,3 do
	for b=a,3 do 
?>		.s<?=a?><?=b?> = a.s<?=a?><?=b?> - b.s<?=a?><?=b?>,
<?	end
end
?>	};
}

static inline sym4 sym4_mul_add(sym4 const a, sym4 const b, real const c) {
	return sym4_add(a, sym4_real_mul(b, c));
}

static inline real4 sym4_real4_mul(sym4 const m, real4 const v) {
	return (real4){
<? for i=0,3 do
?>		0. <? 
	for j=0,3 do 
?> + m.s<?=sym(i,j)?> * v.s<?=j?><? 
	end ?>,
<? end 
?>	};
}

static inline real sym4_det(sym4 const m) {
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

#define real4x4s4_zero ((real4x4s4){ \
	.t = sym4_zero, \
	.x = sym4_zero, \
	.y = sym4_zero, \
	.z = sym4_zero, \
})

static inline real4x4s4 tensor_4sym4_real_mul(real4x4s4 const a, real const s) {
	return (real4x4s4){
<? for a=0,3 do
?>		.s<?=a?> = sym4_real_mul(a.s<?=a?>, s),
<? end 
?>	};
}

static inline real4x4s4 tensor_4sym4_sub(real4x4s4 const a, real4x4s4 const b) {
	return (real4x4s4){
<? for a=0,3 do 
?>		.s<?=a?> = sym4_sub(a.s<?=a?>, b.s<?=a?>),
<? end 
?>	};
}

static inline real4s4x4s4 tensor_sym4sym4_add(
	real4s4x4s4 const a,
	real4s4x4s4 const b
) {
	return (real4s4x4s4){
<? 
for a=0,3 do
	for b=a,3 do
?>		.s<?=a..b?> = sym4_add(a.s<?=a..b?>, b.s<?=a..b?>),
<? 	end
end 
?>	};
}

constant int const stDim = <?=stDim?>;	
constant int const sDim = <?=sDim?>;

constant real3 const xmin = _real3(<?=xmin.x?>, <?=xmin.y?>, <?=xmin.z?>);
constant real3 const xmax = _real3(<?=xmax.x?>, <?=xmax.y?>, <?=xmax.z?>);
constant real3 const dx = _real3(<?=
	tonumber(xmax.x - xmin.x) / tonumber(size.x)?>,<?=
	tonumber(xmax.x - xmin.x) / tonumber(size.x)?>,<?=
	tonumber(xmax.x - xmin.x) / tonumber(size.x)?>);

constant real3 const inv_dx = _real3(<?=
	tonumber(size.x) / tonumber(xmax.x - xmin.x)?>,<?=
	tonumber(size.x) / tonumber(xmax.x - xmin.x)?>,<?=
	tonumber(size.x) / tonumber(xmax.x - xmin.x)?>);

#define getX(i) _real3( \
	xmin.x + ((real)i.x + .5)/(real)size.x * (xmax.x - xmin.x),	\
	xmin.y + ((real)i.y + .5)/(real)size.y * (xmax.y - xmin.y),	\
	xmin.z + ((real)i.z + .5)/(real)size.z * (xmax.z - xmin.z));

#define real4x4x4s4_zero ((real4x4x4s4){ \
	.t = real4x4s4_zero, \
	.x = real4x4s4_zero, \
	.y = real4x4s4_zero, \
	.z = real4x4s4_zero, \
})

real4x4x4s4 calc_dGammaLULL(
	global real4x4s4 const * const GammaULLs
) {
	//initKernel();
	int4 i = globalInt4();
	if (i.x >= size.x || i.y >= size.y || i.z >= size.z) {
		return real4x4x4s4_zero;
	}
	int index = indexForInt4ForSize(i, size.x, size.y, size.z);

	real4x4x4s4 dGammaLULL;
	dGammaLULL.s0 = real4x4s4_zero;
	<? for i=0,sDim-1 do ?>{
		real4x4s4 GammaULL_prev;
		if (i.s<?=i?> > 0) {
			int4 iL = i;
			--iL.s<?=i?>;
			int const indexL = indexForInt4(iL);
			GammaULL_prev = GammaULLs[indexL];
		} else {
			// boundary condition
			GammaULL_prev = real4x4s4_zero;
		}
		
		real4x4s4 GammaULL_next;
		if (i.s<?=i?> < size.s<?=i?> - 1) {
			int4 iR = i;
			++iR.s<?=i?>;
			int const indexR = indexForInt4(iR);
			GammaULL_next = GammaULLs[indexR];
		} else {
			// boundary condition
			GammaULL_next = real4x4s4_zero;
		}
	
		dGammaLULL.s<?=i+1?> = tensor_4sym4_real_mul(
			tensor_4sym4_sub(GammaULL_next, GammaULL_prev),
			.5 * inv_dx.s<?=i?> );
	}<? end ?>

	return dGammaLULL;
}

sym4 calc_EinsteinLL(
	global sym4 const * const gLLs,
	global sym4 const * const gUUs,
	global real4x4s4 const * const GammaULLs
) {
	int4 const i = globalInt4();
	int const index = indexForInt4(i);
	
	real4x4x4s4 const dGammaLULL = calc_dGammaLULL(GammaULLs);

	//this Ricci calculation differs from the one in calc_dPhi_dgLLs because
	// that one can extract RiemannULLL, which can be used for RicciLL calcs
	// but this one doesn't need RiemannULLL, so we can contract one of the terms in RicciLL's calcs

	global real4x4s4 const * const GammaULL = GammaULLs + index;
	real4 const Gamma12L = (real4){
<? 
for a=0,stDim-1 do 
?>		0. <? 
	for b=0,stDim-1 do 
?> + GammaULL->s<?=b?>.s<?=sym(b,a)?><? 
	end
	if a < stDim-1 then ?>,
<? 	end
end ?>
};

	sym4 const RicciLL = {
<? 
for a=0,stDim-1 do
	for b=a,stDim-1 do
?>		.s<?=a?><?=b?> = 0.<?
		for c=0,stDim-1 do ?>
			+ dGammaLULL.s<?=c?>.s<?=c?>.s<?=a..b?> 
			- dGammaLULL.s<?=b?>.s<?=c?>.s<?=sym(c,a)?> 
			+ Gamma12L.s<?=c?> * GammaULL->s<?=c?>.s<?=a..b?><?
			for d=0,stDim-1 do ?>
			- GammaULL->s<?=c?>.s<?=sym(d,b)?> * GammaULL->s<?=d?>.s<?=sym(c,a)?><?
			end
		end 
?>,
<?	end
end ?>
	};

	real const Gaussian = sym4_dot(RicciLL, gUUs[index]);
	return sym4_mul_add(
		RicciLL,
		gLLs[index], -.5 * Gaussian
	);
}

sym4 calc_8piTLL(
	sym4 const gLL,
	<?=TPrim_t?> const TPrim
) {
	sym4 _8piTLL = sym4_zero;

<? 
if solver.body.useEM then 
?>
	
	/*
	assume the E and B fields are upper 3-vectors
	T_ab = F_au F_b^u - 1/4 g_ab F_uv F^uv
	*/

	real4 const EU = (real4)(0 <?for i=0,2 do ?>, TPrim.E.s<?=i?> <? end ?>); 
	real4 const EL = sym4_real4_mul(gLL, EU);
	real const ESq = dot(EL, EU);
	
	real4 const BU = (real4)(0 <?for i=0,2 do ?>, TPrim.B.s<?=i?> <? end ?>); 
	real4 const BL = sym4_real4_mul(gLL, BU);
	real const BSq = dot(BL, BU);

	real const sqrt_det_g = sqrt(fabs(sym4_det(gLL)));
	real3 const SL = real3_real_mul(
		real3_cross(TPrim.E, TPrim.B),
		sqrt_det_g
	);
	
	_8piTLL.s00 += ESq + BSq;
<? 
	for i=0,sDim-1 do 
?>	_8piTLL.s0<?=i+1?> += -2. * SL.s<?=i?>;
<? 		for j=i,sDim-1 do
?>	_8piTLL.s<?=i+1?><?=j+1?> += gLL.s<?=i?><?=j?> * (ESq + BSq) <?
			?>- 2. * (<?
				?>EL.s<?=i+1?> * EL.s<?=j+1?> <?
				?>+ BL.s<?=i+1?> * BL.s<?=j+1?>);
<? 		end
	end 
end
if solver.body.useMatter then
	if solver.body.useVel then 
?>	//if we're using velocity ...
	//set vU.t = 0 so we only lower by the spatial component of the metric.  right?
	real4 const vU = (real4)(0, TPrim.v.x, TPrim.v.y, TPrim.v.z);
	real4 const vL = sym4_real4_mul(gLL, vU);
	real const vLenSq = dot(vL, vU);	//vU.t = 0 so we'll neglect the vL.t component
	real const W = 1. / sqrt(1. - sqrt(vLenSq));
	real4 const uU = (real4)(W, W * vU.s1, W * vU.s2, W * vU.s3);
	real4 const uL = sym4_real4_mul(gLL, uU);
	<? else ?>//otherwise uL = gLL.s0
	real4 const uL = (real4)(gLL.s00, gLL.s01, gLL.s02, gLL.s03);
<?
	end 
?>

	//8 π T_matter_ab = 8 π (u_a u_b (ρ (1 + eInt) + P) + g_ab P)
	sym4 const _8piT_matter_LL = sym4_real_mul(
		sym4_add(
			sym4_real_mul(
				sym4_outer(uL),
				TPrim.rho * (1. + TPrim.eInt) + TPrim.P
			),
			sym4_real_mul(
				gLL,
				TPrim.P
			)
		), 8. * M_PI);
	_8piTLL = sym4_add(_8piTLL, _8piT_matter_LL);
<? end ?>

	return _8piTLL;
}
