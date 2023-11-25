// here is the code used for updating weights based on the EFE

_Static_assert(sizeof(real4s4) == sizeof(real) * 10, "here");
_Static_assert(sizeof(gPrim_t) == sizeof(real) * 10, "here");

kernel void calc_partial_gPrim_of_Phis_kernel(
	constant env_t const * const env,
	global gPrim_t * const partial_gPrim_of_Phis,
	global <?=TPrim_t?> const * const TPrims,
	global gPrim_t const * const gPrims,
	global real4x4s4 const * const GammaULLs,
	global real4s4 const * const EFEs
) {
	initKernel();

	partial_gPrim_of_Phis[index] = calc_partial_gPrim_of_Phi(
		env,
		i,
		TPrims,
		gPrims,
		GammaULLs,
		EFEs
	);

#if 0
	//scale up our gradient?
	//scale by c^4 / G ~ 1e+44
	// which is the units of conversion
	//c^4/G * G_ab = 8 π T_ab
	// but why not just do this in the lambda?
	partial_gPrim_of_Phis[index].alpha *= c*c*c*c/G;
	partial_gPrim_of_Phis[index].betaU = real3_real_mul(partial_gPrim_of_Phis[index].betaU, c*c*c*c/G);
	partial_gPrim_of_Phis[index].gammaLL = real3s3_real_mul(partial_gPrim_of_Phis[index].gammaLL, c*c*c*c/G);
#endif
}

kernel void update_gPrims(
	constant env_t const * const env,
	global gPrim_t * const gPrims,
	global gPrim_t const * const partial_gPrim_of_Phis,
	real const updateLambda
) {
	initKernel();
	global gPrim_t * const gPrim = gPrims + index;
	gPrim_t const partial_gPrim_of_Phi = partial_gPrim_of_Phis[index];

<? if solver.convergeAlpha then
?>	gPrim->alpha -= partial_gPrim_of_Phi.alpha * updateLambda;
<? end
if solver.convergeBeta then
?>	gPrim->betaU.x -= partial_gPrim_of_Phi.betaU.x * updateLambda;
	gPrim->betaU.y -= partial_gPrim_of_Phi.betaU.y * updateLambda;
	gPrim->betaU.z -= partial_gPrim_of_Phi.betaU.z * updateLambda;
<? end
if solver.convergeGamma then
?>	gPrim->gammaLL.xx -= partial_gPrim_of_Phi.gammaLL.xx * updateLambda;
	gPrim->gammaLL.xy -= partial_gPrim_of_Phi.gammaLL.xy * updateLambda;
	gPrim->gammaLL.xz -= partial_gPrim_of_Phi.gammaLL.xz * updateLambda;
	gPrim->gammaLL.yy -= partial_gPrim_of_Phi.gammaLL.yy * updateLambda;
	gPrim->gammaLL.yz -= partial_gPrim_of_Phi.gammaLL.yz * updateLambda;
	gPrim->gammaLL.zz -= partial_gPrim_of_Phi.gammaLL.zz * updateLambda;
<? end
?>
}
