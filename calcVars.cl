// init

kernel void init_gPrims(
	global gPrim_t* gPrims
) {
	INIT_KERNEL();
	
	real3 x = getX(i);
	real r = real3_len(x);

	global gPrim_t* gPrim = gPrims + index;

	//init to flat by default
	*gPrim = (gPrim_t){
		.alpha = 1,
		.betaU = real3_zero,
		.gammaLL = sym3_ident,
	};

	<?=initCond.code?>
}

kernel void init_TPrims(
	global TPrim_t* TPrims
) {
	INIT_KERNEL();

	real3 x = getX(i);
	real r = real3_len(x);
	
	global TPrim_t* TPrim = TPrims + index;
	
	*TPrim = (TPrim_t){
<? if body.useMatter then ?>
		.rho = 0,
		.eInt = 0,
		.P = 0,
<? 	if body.useVel then ?>
		.v = real3_zero,
<? 	end
end 
if body.useEM then ?>		
		.E = real3_zero,
		.B = real3_zero,
<? end ?>
	};

	<?=body.init?>
}

// compute buffers to compute EFE

kernel void calc_gLLs_and_gUUs(
	global sym4* gLLs,
	global sym4* gUUs,
	const global gPrim_t* gPrims
) {
	INIT_KERNEL();

	global const gPrim_t* gPrim = gPrims + index;

	real alphaSq = gPrim->alpha * gPrim->alpha;
	real3 betaL = sym3_real3_mul(gPrim->gammaLL, gPrim->betaU);
	real betaSq = real3_dot(betaL, gPrim->betaU);

	global sym4* gLL = gLLs + index;
	gLL->s00 = -alphaSq + betaSq;
	<? for i=0,subDim-1 do ?>
	gLL->s0<?=i+1?> = betaL.s<?=i?>;
		<? for j=i,subDim-1 do ?>
	gLL->s<?=i+1?><?=j+1?> = gPrim->gammaLL.s<?=i?><?=j?>;
		<? end ?>
	<? end ?>

	real det_gammaLL = sym3_det(gPrim->gammaLL);
	sym3 gammaUU = sym3_inv(det_gammaLL, gPrim->gammaLL);

	global sym4* gUU = gUUs + index;
	gUU->s00 = -1./alphaSq;
	<? for i=0,subDim-1 do ?>
	gUU->s0<?=i+1?> = gPrim->betaU.s<?=i?> / alphaSq;
		<? for j=i,subDim-1 do ?>
	gUU->s<?=i+1?><?=j+1?> = gammaUU.s<?=i?><?=j?> - gPrim->betaU.s<?=i?> * gPrim->betaU.s<?=j?> / alphaSq;
		<? end ?>
	<? end ?>
}

kernel void calc_GammaULLs(
	global tensor_4sym4* GammaULLs,
	const global sym4* gLLs,
	const global sym4* gUUs
) {
	INIT_KERNEL();

	//here's where the finite difference stuff comes in ...
	tensor_4sym4 dgLLL;
	dgLLL.s0 = sym4_zero;
	<? for i=0,gridDim-1 do ?>{
		int4 iL = i;
		iL.s<?=i?> = max(i.s<?=i?> - 1, 0);
		int indexL = indexForInt4(iL);
		global const sym4* gLL_prev = gLLs + indexL;
		
		int4 iR = i;
		iR.s<?=i?> = min(i.s<?=i?> + 1, size.s<?=i?> - 1);
		int indexR = indexForInt4(iR);
		global const sym4* gLL_next = gLLs + indexR;
		
		dgLLL.s<?=i+1?> = (sym4){
		<? for a=0,dim-1 do ?>
			<? for b=a,dim-1 do ?>
			.s<?=a..b?> = (gLL_next->s<?=a..b?> - gLL_prev->s<?=a..b?>) * .5 * inv_dx.s<?=i?>,
			<? end ?>
		<? end ?>
		};
	}<? end ?>

	tensor_4sym4 GammaLLL = (tensor_4sym4){
	<? for a=0,dim-1 do ?>
		<? for b=0,dim-1 do ?>
			<? for c=b,dim-1 do ?>
		.s<?=a?>.s<?=b?><?=c?> = .5 * (
			dgLLL.s<?=c?>.s<?=sym(a,b)?>
			+ dgLLL.s<?=b?>.s<?=sym(a,c)?>
			- dgLLL.s<?=a?>.s<?=sym(b,c)?>),
			<? end ?>
		<? end ?>
	<? end ?>};

	global const sym4* gUU = gUUs + index;
	global tensor_4sym4* GammaULL = GammaULLs + index;
	<? for a=0,dim-1 do ?>
		<? for b=0,dim-1 do ?>
			<? for c=b,dim-1 do ?>
	GammaULL->s<?=a?>.s<?=b?><?=c?> = 0.
				<? for d=0,dim-1 do ?>
		+ gUU->s<?=sym(a,d)?> * GammaLLL.s<?=d?>.s<?=b?><?=c?>
				<? end ?>;
			<? end ?>
		<? end ?>
	<? end ?>
}

kernel void calc_EFEs(
	global sym4* EFEs,
	global const gPrim_t* gPrims,
	global const TPrim_t* TPrims,
	global const sym4* gLLs,
	global const sym4* gUUs,
	global const tensor_4sym4* GammaULLs
) {
	INIT_KERNEL();
	sym4 EinsteinLL = calc_EinsteinLL(gLLs, gUUs, GammaULLs);
	sym4 _8piTLL = calc_8piTLL(gLLs+index, TPrims+index);
	EFEs[index] = sym4_sub(EinsteinLL, _8piTLL);
}
