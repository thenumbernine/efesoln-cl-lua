#pragma once

#include "autogen.h"	//vec3sz_t
#include "math.hpp"	//real3

struct env_t {
	vec3sz_t size = {};
	vec3sz_t stepsize = {};
	real3 xmin = {};
	real3 xmax = {};
	int dim = {};
};

struct gPrim_t {
	union {
		real s[10];
		struct {
			real __attribute__((packed)) alpha = 1;
			real3 __attribute__((packed)) betaU = real3{0,0,0};
			real3s3 __attribute__((packed)) gammaLL = real3s3{1,0,0,1,0,1};
		};
	};

	gPrim_t() {}

	gPrim_t(
		real const alpha_,
		real3 const & betaU_,
		real3s3 const & gammaLL_
	) : alpha(alpha_),
		betaU(betaU_),
		gammaLL(gammaLL_)
	{}
};

static_assert(sizeof(gPrim_t) == sizeof(real) * 10, "here");
static_assert(sizeof(real4s4) == sizeof(real) * 10, "here");

// putting the type header in one place for ffi as well
// and the function / macro headers here:

extern constant real const c;
extern constant real const G;
extern constant int const sym3[3][3];
extern constant int const sym4[4][4];

<?
local d1coeffs = derivCoeffs[1][solver.diffOrder]
?>
extern constant real const d1coeffs[<?=#d1coeffs+1?>];
real d1coeff_for_offset(int offest);

<?
local d2coeffs = derivCoeffs[2][solver.diffOrder]
?>
extern constant real const d2coeffs[<?=#d2coeffs+1?>];

#define real_dot(a, b) ((a) * (b))

<?
local function makeOpsHeader(ctype, fieldtype, fields)
	for _,op in ipairs{"+", "-"} do
?>
inline <?=ctype?> operator<?=op?>(<?=ctype?> const & a, <?=ctype?> const & b) {
	return <?=ctype?>{
<?		for _,field in ipairs(fields) do
?>		a.<?=field?> <?=op?> b.<?=field?>,
<?		end
?>	};
}
inline <?=ctype?> & operator<?=op?>=(<?=ctype?> & a, <?=ctype?> const & b) {
	return a = a <?=op?> b;
}
<?
	end
	for _,op in ipairs{"*", "/"} do
?>
inline <?=ctype?> operator<?=op?>(<?=ctype?> const & a, real const b) {
	return <?=ctype?>{
<?		for _,field in ipairs(fields) do
?>		a.<?=field?> <?=op?> b,
<?		end
?>	};
}
inline <?=ctype?> & operator<?=op?>=(<?=ctype?> & a, real const b) {
	return a = a <?=op?> b;
}
<?
	end
?>

inline <?=ctype?> operator*(real const a, <?=ctype?> const & b) {
	return <?=ctype?>{
<?		for _,field in ipairs(fields) do
?>		a * b.<?=field?>,
<?		end
?>	};
}

real <?=ctype?>_dot(<?=ctype?> const & a, <?=ctype?> const & b);
real <?=ctype?>_lenSq(<?=ctype?> const & a);
real <?=ctype?>_len(<?=ctype?> const & a);

// "norm" name for vectors and tensors
#define <?=ctype?>_normSq <?=ctype?>_lenSq
#define <?=ctype?>_norm <?=ctype?>_len

<?
end
?>

<?makeOpsHeader("real3", "real", {"x", "y", "z"})?>
real3 real3_cross(real3 const & a, real3 const & b);

real3 real4_to_real3(real4 const & a);
real4 real3_to_real4(real3 const & a);
<?makeOpsHeader("real4", "real", {"s0", "s1", "s2", "s3"})?>

#define real3s3_ident (real3s3{1,0,0,1,0,1})

real real3s3_det(real3s3 const & m);
<?makeOpsHeader("real3s3", "real", {"s00", "s01", "s02", "s11", "s12", "s22"})?>
real3s3 real3s3_inv(real3s3 const & m, real const d);

inline real3 operator*(
	real3s3 const & m,
	real3 const & v
) {
	return real3{
		m.s00 * v.s0 + m.s01 * v.s1 + m.s02 * v.s2,
		m.s01 * v.s1 + m.s11 * v.s1 + m.s12 * v.s2,
		m.s02 * v.s2 + m.s12 * v.s1 + m.s22 * v.s2,
	};
}



<?makeOpsHeader("real4x4", "real4", {"s0", "s1", "s2", "s3"})?>
<?makeOpsHeader("real4s4", "real", {"s00", "s01", "s02", "s03", "s11", "s12", "s13", "s22", "s23", "s33"})?>
real4s4 real4s4_outer(real4 const & v);

//a_i = b_ij c_j
inline real4 operator*(
	real4s4 const & m,
	real4 const & v
) {
	real4 result = {};
	for (int i = 0; i < 4; ++i) {
		real sum = {};
		for (int j = 0; j < 4; ++j) {
			sum += m.s[sym4[i][j]] * v.s[j];
		}
		result.s[i] = sum;
	}
	return result;
}

//a_ij = b_ik c_kj
inline real4x4 operator*(
	real4s4 const & a,
	real4s4 const & b
) {
	real4x4 result = {};
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			real sum = {};
			for (int k = 0; k < 4; ++k) {
				sum += a.s[sym4[i][k]] * b.s[sym4[k][j]];
			}
			result.s[i].s[j] = sum;
		}
	}
	return result;
}

real3 real4s4_i0(real4s4 const & a);
real3s3 real4s4_ij(real4s4 const & a);
real real4s4_det(real4s4 const & m);
real4s4 new_real4s4_Minkowski();
extern constant real4s4 const real4s4_Minkowski;

//a_ij = b_ik c_kj
inline real4s4 real4x4_real4s4_to_real4s4_mul(
	real4x4 const & a,
	real4s4 const & b
) {
	real4s4 result = {};
	for (int i = 0; i < 4; ++i) {
		for (int j = i; j < 4; ++j) {
			real sum = {};
			for (int k = 0; k < 4; ++k) {
				sum += a.s[i].s[k] * b.s[sym4[k][j]];
			}
			result.s[sym4[i][j]] = sum;
		}
	}
	return result;
}


real real4x4_tr(real4x4 const & a);

<?makeOpsHeader("real4x4x4", "real4x4", {"s0", "s1", "s2", "s3"})?>

<?makeOpsHeader("real4x4s4", "real4s4", {"s0", "s1", "s2", "s3"})?>

real3 real4x4s4_i00(real4x4s4 const & a);
real4 real4x4s4_tr12(real4x4s4 const & a);
real4 real4x4x4_tr23(real4x4x4 const & a);
real4x4s4 real4s4_real4x4s4_mul(real4s4 const & a, real4x4s4 const & b);
real4x4x4 real4x4s4_real4s4_mul21(real4x4s4 const & a, real4s4 const & b);

<?makeOpsHeader("real4s4x4s4", "real4s4", {"s00", "s01", "s02", "s03", "s11", "s12", "s13", "s22", "s23", "s33"})?>

<?makeOpsHeader("real4x4x4x4", "real4x4x4", {"s0", "s1", "s2", "s3"})?>

real4x4x4x4 real4x4x4x4_real4s4_mul_1_1(real4x4x4x4 const & a, real4s4 const & b);
real4x4x4x4 real4x4x4x4_real4s4_mul_3_1(real4x4x4x4 const & a, real4s4 const & b);
real4s4 real4x4x4x4_tr13_to_real4s4(real4x4x4x4 const & a);

<?makeOpsHeader("real4x4x4s4", "real4x4s4", {"s0", "s1", "s2", "s3"})?>

//I don't trust constant int const ...
#define stDim		<?=stDim?>
#define sDim		<?=sDim?>
extern constant real3 const xmin;
extern constant real3 const xmax;
extern constant real3 const dx;
extern constant real3 const inv_dx;

real3 getX(
	constant env_t const * const env,
	int4 const i
);
extern constant int4 const int4_dirs[3];
int4 int4_dir(int const dim, int const offset);

real4s4 gLL_from_gPrims_at(
	constant env_t const * const env,
	global gPrim_t const * const gPrims,
	int4 const i
);

real4s4 gUU_from_gPrims_at(
	constant env_t const * const env,
	global gPrim_t const * const gPrims,
	int4 const i
);

real4s4 calc_RicciLL(
	constant env_t const * const env,
	int4 const i,
	global gPrim_t const * const gPrims,
	global real4x4s4 const * const GammaULLs
);

real4s4 calc_EinsteinLL(
	constant env_t const * const env,
	int4 const i,
	global gPrim_t const * const gPrims,
	global real4x4s4 const * const GammaULLs
);

real4s4 calc_8piTLL(
	real4s4 const gLL,
	<?=TPrim_t?> const TPrim
);

//need to be here o calcVars.clcpp can see them
// it calls them based on populated code from efe.initConds
gPrim_t calc_gPrim_flat(real3 const x);
gPrim_t calc_gPrim_stellar_Schwarzschild(real3 const x);
gPrim_t calc_gPrim_stellar_Kerr_Newman(real3 const x);

real4s4 calc_gLL_from_gPrim(gPrim_t const & gPrim);
real4s4 calc_gUU_from_gPrim(gPrim_t const & gPrim);

real4s4 calc_partial_gLL_of_Phi(
	constant env_t const * const env,
	global <?=TPrim_t?> const * const TPrims,
	global gPrim_t const * const gPrims,
	global real4x4s4 const * const GammaULLs,
	global real4s4 const * const EFEs,
	int4 const i
);

gPrim_t calc_partial_gPrim_of_Phi(
	constant env_t const * const env,
	global <?=TPrim_t?> const * const TPrims,
	global gPrim_t const * const gPrims,
	global real4x4s4 const * const GammaULLs,
	global real4s4 const * const EFEs,
	int4 const i
);
