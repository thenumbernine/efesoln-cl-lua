#pragma once

#include "autogen.h"	//vec3sz_t
#include "math.hpp"	//real3

#if defined(CLCPU_ENABLED)
#define constant
#define global
#define local

#include <iostream>
#endif

<?=solver.env_mt.cppcode?>
<?=solver.gPrim_mt.cppcode?>
static_assert(sizeof(gPrim_t) == sizeof(real) * 10);

inline std::ostream& operator<<(std::ostream & o, gPrim_t const & x) {
	return o << "("
		<< "alpha=" << x.alpha << ", "
		<< "betaU=" << x.betaU << ", "
		<< "gammaLL=" << real3x3(x.gammaLL)
		<< ")";
}

<?=solver.TPrim_mt.cppcode?>

static_assert(sizeof(real4s4) == sizeof(real) * 10);

inline std::ostream& operator<<(std::ostream & o, <?=TPrim_t?> const & x) {
	return o << "("
<?	if solver.body.useMatter then ?>
		<< "rho=" << x.rho << ", "
		<< "P=" << x.P << ", "
		<< "eInt=" << x.eInt
<?		if solver.body.useVel then ?>
		<< "v=" << x.v
<?		end ?>
<?	end ?>
<?	if solver.body.useEM then ?>
<?		if self.useFourPotential then ?>
		<< "JL=" << x.JL
		<< "AL=" << x.AL
<?		else ?>
		<< "E=" << x.E
		<< "B=" << x.B
<?		end ?>
<?	end ?>
		<< ")";
}


// putting the type header in one place for ffi as well
// and the function / macro headers here:

extern constant real const c;
extern constant real const G;

<?
local d1coeffs = derivCoeffs[1][solver.diffOrder]
?>
extern constant real const d1coeffs[<?=#d1coeffs+1?>];
real d1coeff_for_offset(int offest);

<?
local d2coeffs = derivCoeffs[2][solver.diffOrder]
?>
extern constant real const d2coeffs[<?=#d2coeffs+1?>];

real3 real4_to_real3(real4 const & a);
real4 real3_to_real4(real3 const & a);

real3 real4s4_i0(real4s4 const & a);
real3s3 real4s4_ij(real4s4 const & a);
extern constant real4s4 const real4s4_Minkowski;

real3 real4x4s4_i00(real4x4s4 const & a);
real4 real4x4s4_tr12(real4x4s4 const & a);
real4 real4x4x4_tr23(real4x4x4 const & a);
real4x4x4 real4x4s4_real4s4_mul21(real4x4s4 const & a, real4s4 const & b);

real4x4x4x4 real4x4x4x4_real4s4_mul_1_1(real4x4x4x4 const & a, real4s4 const & b);
real4x4x4x4 real4x4x4x4_real4s4_mul_3_1(real4x4x4x4 const & a, real4s4 const & b);
real4s4 real4x4x4x4_tr13_to_real4s4(real4x4x4x4 const & a);

//I don't trust constant int const ...
#define stDim		<?=stDim?>
#define sDim		<?=sDim?>

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

real4s4 RicciLL_at(
	constant env_t const * const env,
	global gPrim_t const * const gPrims,
	global real4x4s4 const * const GammaULLs,
	int4 const i
);

real4s4 EinsteinLL_at(
	constant env_t const * const env,
	global gPrim_t const * const gPrims,
	global real4x4s4 const * const GammaULLs,
	int4 const i
);

real4s4 calc_8piTLL(
	real4s4 const gLL,
	<?=TPrim_t?> const TPrim
);

//need to be here o calcVars.clcpp can see them
// it calls them based on populated code from efe.initConds
gPrim_t calc_gPrim_flat();
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
	global real4s4 const * const EFE_LL_minus_half_traces,
	int4 const i
);

gPrim_t calc_partial_gPrim_of_Phi(
	constant env_t const * const env,
	global <?=TPrim_t?> const * const TPrims,
	global gPrim_t const * const gPrims,
	global real4x4s4 const * const GammaULLs,
	global real4s4 const * const EFEs,
	global real4s4 const * const EFE_LL_minus_half_traces,
	int4 const i
);

#if defined(CLCPU_ENABLED)
#undef constant
#undef global
#undef local
#endif
