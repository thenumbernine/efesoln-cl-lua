#pragma once

#include "autogen.h"	//vec3sz_t
#include "math.hpp"	//real3

#if defined(CLCPU_ENABLED)
#define constant
#define global
#define local

#include <iostream>
#endif

//NOTICE this has to match in luajit
// maybe I should autogen it?
// but POD vs class...
struct env_t {
	vec3sz_t size = {};
	vec3sz_t stepsize = {};
	real3 xmin = {};
	real3 xmax = {};
	real3 dx = {};
	real3 invdx = {};
	int dim = {};
};

struct gPrim_t {
// this works in Tensor/Vector.h, but why not here?
//	union {
//		struct {
			real alpha;
			real3 betaU;
			real3s3 gammaLL;
//		};
//		std::array<real, 10> s = {};
//	};

	gPrim_t() {
		alpha = 1;
		betaU = real3(0,0,0);
		gammaLL = real3s3(
			1,
			0,1,
			0,0,1
		);
	}

	gPrim_t(
		real const alpha_,
		real3 const & betaU_,
		real3s3 const & gammaLL_
	) {
		alpha = alpha_;
		betaU = betaU_;
		gammaLL = gammaLL_;
	}

#if 0
	using This = gPrim_t;
	static constexpr auto fields = std::make_tuple(
		std::make_pair("alpha", &This::alpha),
		std::make_pair("betaU", &This::betaU),
		std::make_pair("gammaLL", &This::gammaLL)
	);
#endif
};

static_assert(sizeof(gPrim_t) == sizeof(real) * 10);

inline std::ostream& operator<<(std::ostream & o, gPrim_t const & x) {
	return o << "("
		<< "alpha=" << x.alpha << ", " 
		<< "betaU=" << x.betaU << ", " 
		<< "gammaLL=" << real3x3(x.gammaLL)
		<< ")";
}

//TODO I could give struct-lua an option for generating cpp classes instead of C typedef structs ...
struct <?=TPrim_t?> {
<?	if solver.body.useMatter then ?>
	real rho;
	real P;
	real eInt;
<?		if solver.body.useVel then ?>
	real3 v;
<?		end ?>
<?	end ?>
<?	if solver.body.useEM then ?>
<?		if self.useFourPotential then ?>
	real4 JL;
	real4 AL;
<?		else ?>
	real3 E;
	real3 B;
<?		end ?>
<?	end ?>
};

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

#if defined(CLCPU_ENABLED)
#undef constant
#undef global
#undef local
#endif
