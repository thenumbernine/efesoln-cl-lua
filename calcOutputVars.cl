// compute aux buffers for output

kernel void calc_detGammas(
	global real* detGammas,
	global const gPrim_t* gPrims
) {
	INIT_KERNEL();
	detGammas[index] = sym3_det(gPrims[index].gammaLL);
}

kernel void calc_numericalGravity(
	global real* numericalGravity,
	global const tensor_4sym4* GammaULLs
) {
	INIT_KERNEL();
	real3 x = getX(i);
	real r = real3_len(x);
	global const tensor_4sym4* GammaULL = GammaULLs + index;
	numericalGravity[index] = (0.
		+ GammaULL->s1.s00 * x.s0 / r
		+ GammaULL->s2.s00 * x.s1 / r
		+ GammaULL->s3.s00 * x.s2 / r) * c * c;
}

kernel void calc_analyticalGravity(
	global real* analyticalGravity
) {
	INIT_KERNEL();
	real3 x = getX(i);
	real r = real3_len(x);
	real matterRadius = min(r, (real)<?=solver.body.radius?>);
	real volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
	real m = <?=solver.body.density?> * volumeOfMatterRadius;	// m^3
	real dm_dr = 0;
	analyticalGravity[index] = (2*m * (r - 2*m) + 2 * dm_dr * r * (2*m - r)) / (2 * r * r * r)
		* c * c;	//+9 at earth surface, without matter derivatives
}

kernel void calc_norm_EFE_tts(
	global real* norm_EFE_tts,
	global const sym4* EFEs
) {
	INIT_KERNEL();
	norm_EFE_tts[index] = EFEs[index].s00 / (8. * M_PI) * c * c / G / 1000.;
}

kernel void calc_norm_EFE_tis(
	global real* norm_EFE_tis,
	global const sym4* EFEs
) {
	INIT_KERNEL();
	global const sym4* EFE = EFEs + index;	
	norm_EFE_tis[index] = sqrt(0.
<? for i=0,subDim-1 do ?>
		+ EFE->s0<?=i+1?> * EFE->s0<?=i+1?>
<? end ?>) * c;
}

kernel void calc_norm_EFE_ijs(
	global real* norm_EFE_ijs,
	global const sym4* EFEs
) {
	INIT_KERNEL();
	global const sym4* EFE = EFEs + index;
	norm_EFE_ijs[index] = sqrt(0.
<? for i=0,subDim-1 do
	for j=0,subDim-1 do ?>
		+ EFE->s<?=sym(i+1,j+1)?> * EFE->s<?=sym(i+1,j+1)?>
<?	end
end ?>);
}

kernel void calc_norm_EinsteinLLs(
	global real* norm_EinsteinLLs,
	global const sym4* gLLs,
	global const sym4* gUUs,
	global const tensor_4sym4* GammaULLs
) {
	INIT_KERNEL();
	sym4 EinsteinLL = calc_EinsteinLL(gLLs, gUUs, GammaULLs);
	norm_EinsteinLLs[index] = sqrt(sym4_dot(EinsteinLL, EinsteinLL));
}

//converts solver buffers to float[]
kernel void updateDisplayVar(
	global float* result,

	// alpha-1
	//global const gPrim_t* gPrims

	// numericalGravity
	global const tensor_4sym4* GammaULLs
) {
	INIT_KERNEL();
	
	// alpha-1
	//result[index] = gPrims[index].alpha - 1.;

	// numericalGravity 
	real3 x = getX(i);
	real r = real3_len(x);
	global const tensor_4sym4* GammaULL = GammaULLs + index;
	result[index] = (0.
		+ GammaULL->s1.s00 * x.s0 / r
		+ GammaULL->s2.s00 * x.s1 / r
		+ GammaULL->s3.s00 * x.s2 / r) * c * c;
}
