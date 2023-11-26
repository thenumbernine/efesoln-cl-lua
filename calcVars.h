#pragma once

gPrim_t calc_gPrim_stellar_Schwarzschild(real3 const x);

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
