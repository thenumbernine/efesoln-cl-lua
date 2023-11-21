// I reported this bug to intel like 3 years ago.  it's fixed, right?
//is buggy with doubles on intel opencl ubuntu compiler
//#define _real3(a,b,c)          ((real3){.s={a,b,c}})
//so we do this instead and are safe:
#define _real3(a,b,c)            ((real3){.x=a, .y=b, .z=c})

constant real const c = 299792458;			// m/s
constant real const G = 6.67384e-11;		// m^3 / (kg s^2)

constant int const sym3[3][3] = {
	{0, 1, 2},
	{1, 3, 4},
	{2, 4, 5},
};


constant int const sym4[4][4] = {
	{0, 1, 2, 3},
	{1, 4, 5, 6},
	{2, 5, 7, 8},
	{3, 6, 8, 9},
};

<?
local range = require 'ext.range'
local function makeAdds(ctype)
	local n = solver.diffOrder
	assert(n % 2 == 0)
	for i=4,n,2 do
		local xs = range(0,i-1):mapi(function(j) return "x"..j end)
?>
#define <?=ctype?>_add<?=i?>(<?=xs:concat', '?>)\
	<?=ctype?>_add2(\
		<?=ctype?>_add2(x0, x1),\
		<?=ctype?>_add<?=i-2?>(<?=xs:sub(3):concat', '?>)\
	)
<?	end
end
?>

#define new_real_zero() 0.
#define real_zero new_real_zero()

#define real_add(a,b)	((a) + (b))
#define real_add2		real_add
#define real_sub(a,b)	((a) - (b))
#define real_mul(a,b)	((a) * (b))
#define real_real_mul	real_mul
#define real_dot		real_mul
#define real_div(a,b)	((a) / (b))
#define real_real_div	real_div
<?makeAdds"real"?>

<?
local range = require 'ext.range'
local table = require 'ext.table'
?>

<?
local function makeops(ctype, fieldtype, fields)
?>
#define new_<?=ctype?>_zero() ((<?=ctype?>){\
<? for _,field in ipairs(fields) do --\
?>	.<?=field?> = new_<?=fieldtype?>_zero(),\
<? end --\
?>})
constant <?=ctype?> const <?=ctype?>_zero = new_<?=ctype?>_zero();
<?
	for _,op in ipairs{"add", "sub"} do
?>
static inline <?=ctype?> <?=ctype?>_<?=op?>(<?=ctype?> const a, <?=ctype?> const b) {
	return (<?=ctype?>){
<? 		for _,field in ipairs(fields) do
?>	.<?=field?> = <?=fieldtype?>_<?=op?>(a.<?=field?>, b.<?=field?>),
<? 		end
?>	};
}
<?	end
	for _,op in ipairs{"mul", "div"} do
?>
static inline <?=ctype?> <?=ctype?>_real_<?=op?>(<?=ctype?> const a, real const b) {
	return (<?=ctype?>){
<? 		for _,field in ipairs(fields) do
?>	.<?=field?> = <?=fieldtype?>_real_<?=op?>(a.<?=field?>, b),
<? 		end
?>	};
}
<?	end
?>

#define <?=ctype?>_add2 <?=ctype?>_add
<?makeAdds(ctype)?>

static inline <?=ctype?> <?=ctype?>_mul_add(<?=ctype?> const a, <?=ctype?> const b, real const c) {
	return <?=ctype?>_add(a, <?=ctype?>_real_mul(b, c));
}

static inline real <?=ctype?>_dot(<?=ctype?> const a, <?=ctype?> const b) {
	return real_zero
<?	for _,field in ipairs(fields) do
?>		+ <?=fieldtype?>_dot(a.<?=field?>, b.<?=field?>)
<?	end
?>	;
}

// "length" name for vectors
static inline real <?=ctype?>_lenSq(<?=ctype?> const a) {
	return <?=ctype?>_dot(a, a);
}
static inline real <?=ctype?>_len(<?=ctype?> const a) {
	return sqrt(<?=ctype?>_lenSq(a));
}

// "norm" name for vectors and tensors
#define <?=ctype?>_normSq <?=ctype?>_lenSq
#define <?=ctype?>_norm <?=ctype?>_len

<?
end
?>

<?makeops("real3", "real", {"x", "y", "z"})?>

//ε_ijk a_j b_k
static inline real3 real3_cross(real3 const a, real3 const b) {
	return _real3(
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x);
}

// ok weird convention ...
// for spacetime vars, real3 is xyz, real4 is used such that .s0123 == txyz (time dimension first)
// but for opencl vars (i, etc), real4 is named xyzw (time dimension last)
// this is going to assume the former (spacetime vars)
static inline real3 real4_to_real3(real4 const a) {
	return _real3(a.s1, a.s2, a.s3);
}

static inline real4 real3_to_real4(real3 const a) {
	return (real4){.s0 = 0, .s1 = a.s0, .s2 = a.s1, .s3 = a.s2};
}

<?makeops("real4", "real", {"s0", "s1", "s2", "s3"})?>

#define real3s3_ident ((real3s3){\
	.xx = 1, .xy = 0, .xz = 0,\
	.yy = 1, .yz = 0,\
	.zz = 1,\
})

static inline real real3s3_det(real3s3 const m) {
	return m.s00 * m.s11 * m.s22
		+ m.s01 * m.s12 * m.s02
		+ m.s02 * m.s01 * m.s12
		- m.s02 * m.s11 * m.s02
		- m.s01 * m.s01 * m.s22
		- m.s00 * m.s12 * m.s12;
}

<?makeops("real3s3", "real", {"s00", "s01", "s02", "s11", "s12", "s22"})?>

static inline real3s3 real3s3_inv(real const d, real3s3 const m) {
	return (real3s3){
		.xx = (m.yy * m.zz - m.yz * m.yz) / d,
		.xy = (m.xz * m.yz - m.xy * m.zz) / d,
		.xz = (m.xy * m.yz - m.xz * m.yy) / d,
		.yy = (m.xx * m.zz - m.xz * m.xz) / d,
		.yz = (m.xz * m.xy - m.xx * m.yz) / d,
		.zz = (m.xx * m.yy - m.xy * m.xy) / d,
	};
}

static inline real3 real3s3_real3_mul(real3s3 const m, real3 const v) {
	return (real3){
		.x = m.xx * v.x + m.xy * v.y + m.xz * v.z,
		.y = m.xy * v.y + m.yy * v.y + m.yz * v.z,
		.z = m.xz * v.z + m.yz * v.y + m.zz * v.z,
	};
}

<?makeops("real4x4", "real4", {"s0", "s1", "s2", "s3"})?>

<?makeops("real4s4", "real", {"s00", "s01", "s02", "s03", "s11", "s12", "s13", "s22", "s23", "s33"})?>

static inline real4s4 real4s4_outer(real4 const v) {
	return (real4s4){
<? for a=0,3 do
	for b=a,3 do
?>		.s<?=a..b?> = v.s<?=a?> * v.s<?=b?>,
<? 	end
end
?>	};
}

//a_i = b_ij c_j
static inline real4 real4s4_real4_mul(
	real4s4 const m,
	real4 const v
) {
	return (real4){
<? for i=0,3 do
?>		.s<?=i?> = 0. <?
	for j=0,3 do
?> + m.s<?=sym(i,j)?> * v.s<?=j?><?
	end ?>,
<? end
?>	};
}

//a_ij = b_ik c_kj
static inline real4x4 real4s4_real4s4_mul(
	real4s4 const a,
	real4s4 const b
) {
	return (real4x4){
<? for a=0,3 do
?>		.s<?=a?> = (real4){
<?
	for b=0,3 do
?>			.s<?=b?> = 0.<?
		for c=0,3 do
?> + a.s<?=sym(a,c)?> * b.s<?=sym(c,b)?><?
		end ?>,
<?
	end
?>		},
<? end
?>	};
}

static real3 real4s4_i0(
	real4s4 const a
) {
	return _real3(a.s01, a.s02, a.s03);
}

static real3s3 real4s4_ij(
	real4s4 const a
) {
	return (real3s3){
<?
for i=0,2 do
	for j=i,2 do
?>		.s<?=i..j?> = a.s<?=i+1?><?=j+1?>,
<?	end
end
?>	};
}

static inline real real4s4_det(real4s4 const m) {
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

real4s4 new_real4s4_Minkowski() {
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
constant real4s4 const real4s4_Minkowski = (real4s4){
<?
for a=0,3 do
	for b=a,3 do
?>		.s<?=a..b?> = <?=a==b and (a==0 and -1 or 1) or 0?>,
<?	end
end
?>	};



//a_ij = b_ik c_kj
static inline real4s4 real4x4_real4s4_to_real4s4_mul(
	real4x4 const a,
	real4s4 const b
) {
	return (real4s4){
<? for a=0,3 do
	for b=a,3 do
?>		.s<?=a..b?> = 0.<?
		for c=0,3 do
		?> + a.s<?=a?>.s<?=c?> * b.s<?=sym(c,b)?><?
	end ?>,
<?	end
end
?>	};
}

//a = b_ii
static inline real real4x4_tr(real4x4 const a) {
	return 0.<?
for a=0,3 do
?> + a.s<?=a?>.s<?=a?><?
end ?>;
}

<?makeops("real4x4x4", "real4x4", {"s0", "s1", "s2", "s3"})?>

<?makeops("real4x4s4", "real4s4", {"s0", "s1", "s2", "s3"})?>

static inline real3 real4x4s4_i00(real4x4s4 const a) {
	return _real3(a.s1.s00, a.s2.s00, a.s3.s00);
}

//b_i = a^j_ji
static inline real4 real4x4s4_tr12(real4x4s4 const a) {
	return (real4){
<? for i=0,3 do
?>		.s<?=i?> = 0.<?
	for j=0,3 do
?> + a.s<?=j?>.s<?=sym(j,i)?><?
	end
?>,
<? end
?>	};
}

//b^i = a^ij_j
static inline real4 real4x4x4_tr23(real4x4x4 const a) {
	return (real4){
<? for i=0,3 do
?>		.s<?=i?> = real4x4_tr(a.s<?=i?>),
<? end
?>	};
}

//c^i_jk = a^il b_ljk
static inline real4x4s4 real4s4_real4x4s4_mul(
	real4s4 const a,
	real4x4s4 const b
) {
	return (real4x4s4){
<? for i=0,3 do
?>		.s<?=i?> = (real4s4){
<?	for j=0,3 do
		for k=j,3 do
?>			.s<?=sym(j,k)?> = 0.<?
			for l=0,3 do
?> + a.s<?=sym(i,l)?> * b.s<?=l?>.s<?=j..k?><?
			end
?>,
<?		end
	end
?>		},
<? end
?>	};
}

//c^ij_k = a^i_mk b^mj
static inline real4x4x4 real4x4s4_real4s4_mul21(
	real4x4s4 const a,
	real4s4 const b
) {
	return (real4x4x4){
<? for i=0,3 do
?>		.s<?=i?> = real4s4_real4s4_mul(b, a.s<?=i?>),
<? end
?>	};
}

<?makeops("real4s4x4s4", "real4s4", {"s00", "s01", "s02", "s03", "s11", "s12", "s13", "s22", "s23", "s33"})?>

<?makeops("real4x4x4x4", "real4x4x4", {"s0", "s1", "s2", "s3"})?>

static inline real4x4x4x4 real4s4_real4x4x4x4_mul(
	real4s4 const a,
	real4x4x4x4 const b
) {
	return (real4x4x4x4){
<? for a=0,3 do
?>		.s<?=a?> = (real4x4x4){
<?	for b=0,3 do
?>			.s<?=b?> = (real4x4){
<?		for c=0,3 do
?>				.s<?=c?> = (real4){
<?			for d=0,3 do
?>					.s<?=d?> = 0.<?
				for e=0,3 do
?> + a.s<?=sym(a,e)?> * b.s<?=e?>.s<?=b?>.s<?=c?>.s<?=d?><?
				end
?>					,
<?			end
?>				},
<?		end
?>			},
<? 	end
?>		},
<? end
?>	};
}

//b_ij = a^k_ikj
// assuming b_ij = b_ji i.e. a_ijkl = a_klij
static inline real4s4 real4x4x4x4_tr13_to_real4s4(real4x4x4x4 const a) {
	return (real4s4){
<? for a=0,3 do
	for b=a,3 do
?>		.s<?=a..b?> = 0. <?
		for c=0,3 do
?> + a.s<?=c?>.s<?=a?>.s<?=c?>.s<?=b?><?
		end ?>,
<? 	end
end
?>	};
}

<?makeops("real4x4x4s4", "real4x4s4", {"s0", "s1", "s2", "s3"})?>

<?makeops("real4s4x4x4s4", "real4x4s4", {"s00", "s01", "s02", "s03", "s11", "s12", "s13", "s22", "s23", "s33"})?>


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

static inline real3 getX(int4 i) {
	return _real3(
		xmin.x + ((real)i.x + .5)/(real)size.x * (xmax.x - xmin.x),
		xmin.y + ((real)i.y + .5)/(real)size.y * (xmax.y - xmin.y),
		xmin.z + ((real)i.z + .5)/(real)size.z * (xmax.z - xmin.z));
}

//[1/m^2]
real4x4x4s4 calc_partial_xU_of_GammaULL(
	global real4x4s4 const * const GammaULLs
) {
	//initKernel();
	int4 i = globalInt4();
	if (i.x >= size.x || i.y >= size.y || i.z >= size.z) {
		return real4x4x4s4_zero;
	}
	int index = indexForInt4ForSize(i, size.x, size.y, size.z);

	//partial_xU_of_GammaULL.a.b.cd := ∂/∂x^a(Γ^b_cd)
<?=solver:finiteDifference{
	bufferName = "GammaULLs",
	srcType = "4x4s4",
	resultName = "partial_xU_of_GammaULL",
}?>

	return partial_xU_of_GammaULL;
}

//[1/m^2]
real4s4 calc_EinsteinLL(
	global real4s4 const * const gLLs,
	global real4s4 const * const gUUs,
	global real4x4s4 const * const GammaULLs
) {
	int4 const i = globalInt4();
	if (i.x >= size.x || i.y >= size.y || i.z >= size.z) {
		return real4s4_zero;
	}
	int const index = indexForInt4(i);

	real4x4x4s4 const partial_xU_of_GammaULL = calc_partial_xU_of_GammaULL(GammaULLs);

	//this Ricci calculation differs from the one in calc_partial_gLL_of_Phis because
	// that one can extract RiemannULLL, which can be used for RicciLL calcs
	// but this one doesn't need RiemannULLL, so we can contract one of the terms in RicciLL's calcs

	// TODO use 2nd derivs + 2nd-deriv-finite-difference stencil

	real4x4s4 const GammaULL = GammaULLs[index];

<? if true then -- used 2nd-deriv finite-difference stencil
?>
	real4s4 const gLL = gLLs[index];
	real4s4 const gUU = gUUs[index];

	//GammaUUL.a.b.c := Γ^ab_c = Γ^a_dc g^db
	real4x4x4 const GammaUUL = real4x4s4_real4s4_mul21(GammaULL, gUU);

	//GammaLLL.a.bc := Γ_abc = g_ad Γ^d_bc
	real4x4s4 const GammaLLL = real4s4_real4x4s4_mul(gLL, GammaULL);

	//Gamma23U.a := Γ^a = Γ^au_u
	real4 const Gamma23U = real4x4x4_tr23(GammaUUL);

	//partial_xU2_of_gLL.ab.cd := ∂_a ∂_b (g_cd)
	// = ∂^2/(∂x^a ∂x^b) (g_cd)
	// = g_cd,ab
<?= solver:finiteDifference2{
	bufferName = "gLLs",
	srcType = "4s4",
	resultName = "partial_xU2_of_gLL",
	getBoundary = function(args) return "real4s4_Minkowski" end,
} ?>

<? if false then -- testing to make sure contraction of Riemann == Ricci
-- ... looks the same as the simplified version below ?>

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
	// TODO antisymmetric storage
	real4x4x4x4 GammaSq_asym_LLLL;
	for (int a = 0; a < stDim; ++a) {
		for (int b = 0; b < stDim; ++b) {
			for (int c = 0; c < stDim; ++c) {
				for (int d = 0; d < stDim; ++d) {
					real sum = 0;
					for (int e = 0; e < stDim; ++e) {
						sum += 
							  GammaULL.s[e].s[sym4[a][d]] * GammaLLL.s[e].s[sym4[b][c]]
							- GammaULL.s[e].s[sym4[a][c]] * GammaLLL.s[e].s[sym4[b][d]];
					}
					GammaSq_asym_LLLL.s[a].s[b].s[c].s[d] = sum;
				}
			}
		}
	}

	real4x4x4x4 const gUU_times_partial_xU2_gLL_asym = real4s4_real4x4x4x4_mul(gUU, partial_xU2_of_gLL_asym);
	real4x4x4x4 const GammaSq_asym_ULLL = real4s4_real4x4x4x4_mul(gUU, GammaSq_asym_LLLL);

	//RiemannULLL.a.b.cd := R^a_bcd = 1/2 g^ae ((g_ed,cb - g_bd,ce - g_ec,bd + g_bc,de) + g^fg (Γ_fed Γ_gbc - Γ_fec Γ_gbd))
	//TODO antisymmetric storage
	real4x4x4x4 const RiemannULLL = real4x4x4x4_mul_add(GammaSq_asym_ULLL, gUU_times_partial_xU2_gLL_asym, .5);

	//RicciLL.ab := R_ab = R^c_acb
	real4s4 const RicciLL = real4x4x4x4_tr13_to_real4s4(RiemannULLL);

<? else  -- testing to make sure contraction of Riemann == Ricci ?>

	/*
	R_abcd = 1/2 (g_ad,cb - g_bd,ca - g_ac,bd + g_bc,da) + g^fg (Γ_fad Γ_gbc - Γ_fac Γ_gbd)
	R^a_bcd = g^ae (1/2 (g_ed,cb - g_bd,ce - g_ec,bd + g_bc,de) + g^fg (Γ_fed Γ_gbc - Γ_fec Γ_gbd))
	RicciLL.ab := R_ab = 1/2 g^uv (g_au,bv + g_bv,au - g_ab,uv - g_uv,ab) + Γ^uv_a Γ_uvb - Γ^uv_v Γ_uab
	*/
	real4s4 RicciLL;
	for (int a = 0; a < stDim; ++a) {
		for (int b = a; b < stDim; ++b) {
			int const ab = sym4[a][b];
			real sum = 0;
			for (int u = 0; u < stDim; ++u) {
				int const au = sym4[a][u];
				for (int v = 0; v < stDim; ++v) {
					int const uv = sym4[u][v];
					int const bv = sym4[b][v];
					sum += .5 * gUU.s[uv] * (
							  partial_xU2_of_gLL.s[au].s[bv]
							+ partial_xU2_of_gLL.s[bv].s[au]
							- partial_xU2_of_gLL.s[ab].s[uv]
							- partial_xU2_of_gLL.s[uv].s[ab]
						)
						+ GammaUUL.s[u].s[v].s[a] * GammaLLL.s[u].s[bv];
				}
				sum -= Gamma23U.s[u] * GammaLLL.s[u].s[ab];
			}
			RicciLL.s[ab] = sum;
		}
	}

<? end  -- testing to make sure contraction of Riemann == Ricci ?>

<? else -- only use 1st-deriv finite-difference stencils (and difference-of-difference for 2nd-deriv)
?>

	real4 const Gamma12L = real4x4s4_tr12(GammaULL);
	real4s4 RicciLL;
	for (int a = 0; a < stDim; ++a) {
		for (int b = a; b < stDim; ++b) {
			int const ab = sym4[a][b];
			real sum = 0;
			for (int c = 0; c < stDim; ++c) {
				sum += 
					+ partial_xU_of_GammaULL.s[c].s[c].s[ab]
					- partial_xU_of_GammaULL.s[b].s[c].s[sym4[c][a]]
					+ Gamma12L.s[c] * GammaULL.s[c].s[ab];
				for (int d = 0; d < stDim; ++d) {
					sum -= GammaULL.s[c].s[sym4[d][b]] * GammaULL.s[d].s[sym4[c][a]];
				}
			}
			RicciLL.s[ab] = sum;
		}
	}
<? end
?>

	real const Gaussian = real4s4_dot(RicciLL, gUUs[index]);

	return real4s4_mul_add(
		RicciLL,
		gLLs[index],
		-.5 * Gaussian
	);
}

//[1/m^2]
real4s4 calc_8piTLL(
	real4s4 const gLL,
	<?=TPrim_t?> const TPrim
) {
	real4s4 _8piTLL = real4s4_zero;

<?
if solver.body.useEM then
?>

	/*
	assume the E and B fields are upper 3-vectors
	T_ab = F_au F_b^u - 1/4 g_ab F_uv F^uv
	*/

	real4 const EU = real3_to_real4(TPrim.E);
	real4 const EL = real4s4_real4_mul(gLL, EU);
	real const ESq = real4_dot(EL, EU);

	real4 const BU = real3_to_real4(TPrim.B);
	real4 const BL = real4s4_real4_mul(gLL, BU);
	real const BSq = real4_dot(BL, BU);

	real const sqrt_det_g = sqrt(fabs(real4s4_det(gLL)));
	real3 const SL = real3_real_mul(
		real3_cross(TPrim.E, TPrim.B),
		sqrt_det_g
	);

	_8piTLL.s00 += ESq + BSq;
	for (int i = 0; i < sDim; ++i) {
		_8piTLL.s[sym4[0][i+1]] += -2. * SL.s[i];
		for (int j = i; j < sDim; ++j) {
			_8piTLL.s[sym4[i+1][j+1]] += gLL.s[sym4[i+1][j+1]] * (ESq + BSq) 
				- 2. * (EL.s[i+1] * EL.s[j+1] + BL.s[i+1] * BL.s[j+1]);
		}
	}
<? 
end
if solver.body.useMatter then
	if solver.body.useVel then
?>	//if we're using velocity ...
	//set vU.t = 0 so we only lower by the spatial component of the metric.  right?
	real4 const vU = real3_to_real4(TPrim.v);
	real4 const vL = real4s4_real4_mul(gLL, vU);
	real const vLenSq = real4_dot(vL, vU);	//vU.t = 0 so we'll neglect the vL.t component
	real const W = 1. / sqrt(1. - sqrt(vLenSq));
	real4 const uU = (real4){.s={W, W * vU.s1, W * vU.s2, W * vU.s3}};
	real4 const uL = real4s4_real4_mul(gLL, uU);
	<? else ?>//otherwise uL = gLL.s0
	real4 const uL = (real4){.s={gLL.s00, gLL.s01, gLL.s02, gLL.s03}};
<?
	end
?>

	//8 π T_matter_ab = 8 π (u_a u_b (ρ (1 + eInt) + P) + g_ab P)
	//[1/m^2]
	real4s4 const _8piT_matter_LL = real4s4_real_mul(
		real4s4_add(
			real4s4_real_mul(
				real4s4_outer(uL),
				TPrim.rho * (1. + TPrim.eInt) + TPrim.P
			),
			real4s4_real_mul(
				gLL,
				TPrim.P
			)
		), 8. * M_PI);
	_8piTLL = real4s4_add(_8piTLL, _8piT_matter_LL);
<? end ?>

	return _8piTLL;
}
