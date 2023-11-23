// here is the code used for updating weights based on the EFE

kernel void calc_partial_gPrim_of_Phis_kernel(
	global gPrim_t * const partial_gPrim_of_Phis,
	global <?=TPrim_t?> const * const TPrims,
	global gPrim_t const * const gPrims,
	global real4s4 const * const gLLs,
	global real4s4 const * const gUUs,
	global real4x4s4 const * const GammaULLs,
	global real4s4 const * const EFEs
) {
	initKernel();

	partial_gPrim_of_Phis[index] = calc_partial_gPrim_of_Phi(
		i,
		TPrims,
		gPrims,
		gLLs,
		gUUs,
		GammaULLs,
		EFEs
	);	

#if 0
	//scale up our gradient?
	//scale by c^4 / G ~ 1e+44
	// which is the units of conversion
	//c^4/G * G_ab = 8 π T_ab
	partial_gPrim_of_Phis[index].alpha *= c*c*c*c/G;
	partial_gPrim_of_Phis[index].betaU = real3_real_mul(partial_gPrim_of_Phis[index].betaU, c*c*c*c/G);
	partial_gPrim_of_Phis[index].gammaLL = real3s3_real_mul(partial_gPrim_of_Phis[index].gammaLL, c*c*c*c/G);
#endif
#if 0 //debugging
	partial_gPrim_of_Phis[index].alpha = 0;
	partial_gPrim_of_Phis[index].betaU = real3_zero;
	partial_gPrim_of_Phis[index].gammaLL = real3s3_zero;
#endif

}

kernel void update_gPrims(
	global gPrim_t * const gPrims,
	global gPrim_t const * const partial_gPrim_of_Phis,
	real const updateLambda
) {
	initKernel();
	global gPrim_t * const gPrim = gPrims + index;
	gPrim_t const partial_gPrim_of_Phi = partial_gPrim_of_Phis[index];

#if 0 //debugging
	gPrim->alpha = 1;
	gPrim->betaU = real3_zero;
	gPrim->gammaLL = real3s3_ident;
	return;
#endif

<? if solver.convergeAlpha then
?>	gPrim->alpha -= updateLambda * partial_gPrim_of_Phi.alpha;
<? end
if solver.convergeBeta then
?>	gPrim->betaU = real3_mul_add(gPrim->betaU, partial_gPrim_of_Phi.betaU, -updateLambda);
<? end
if solver.convergeGamma then
?>	gPrim->gammaLL = real3s3_mul_add(gPrim->gammaLL, partial_gPrim_of_Phi.gammaLL, -updateLambda);
<? end
?>
}
