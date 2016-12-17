// here is the code used for updating weights based on the EFE

kernel void calc_dPhi_dgPrims(
	global gPrim_t* dPhi_dgPrims,
	global const <?=TPrim_t?>* TPrims,
	global const gPrim_t* gPrims,
	global const sym4* gLLs,
	global const sym4* gUUs,
	global const tensor_4sym4* GammaULLs,
	global const sym4* EFEs
) {
	initKernel();
	
	global const sym4* gLL = gLLs + index;
	global const sym4* gUU = gUUs + index;
	global const tensor_4sym4* GammaULL = GammaULLs + index;

	//this is also in the Ricci computation, but should I be storing it?  is it too big?
	tensor_44sym4 dGammaLULL;
	calc_dGammaLULL(&dGammaLULL, GammaULLs);

	tensor_44sym4 RiemannULLL = (tensor_44sym4){
<? for a=0,stDim-1 do
?>		.s<?=a?> = (tensor_4sym4){
<? 	for b=0,stDim-1 do
?>			.s<?=b?> = (sym4){
<? 		for c=0,stDim-1 do 
			for d=c,stDim-1 do 
?>				.s<?=c..d?> = dGammaLULL.s<?=c?>.s<?=a?>.s<?=sym(b,d)?>
					- dGammaLULL.s<?=d?>.s<?=a?>.s<?=sym(b,c)?><?
				for e=0,stDim-1 do ?>
					+ GammaULL->s<?=a?>.s<?=sym(e,c)?> * GammaULL->s<?=e?>.s<?=sym(b,d)?>
					- GammaULL->s<?=a?>.s<?=sym(e,d)?> * GammaULL->s<?=e?>.s<?=sym(b,c)?><? 
				end ?>,
<? 			end 
		end
?>			},
<? end
?>		},
<? end 
?>	};

	sym4 RicciLL = (sym4){
<? for a=0,stDim-1 do 
	for b=a,stDim-1 do
?>		.s<?=a..b?> = 0. <?
		for c=0,stDim-1 do 
?> + RiemannULLL.s<?=c?>.s<?=a?>.s<?=sym(c,b)?><? 
		end ?>,
<? 	end
end 
?>	};

	real RicciUL[4][4] = {
<? for a=0,stDim-1 do 
?>		{
<?
	for b=0,stDim-1 do 
?>			0.<?
		for c=0,stDim-1 do 
?> + gUU->s<?=sym(a,c)?> * RicciLL.s<?=sym(c,b)?><?
		end ?>,
<?
	end 
?>		},
<? end 
?>	};

	sym4 RicciUU = (sym4){
<? for a=0,stDim-1 do
	for b=a,stDim-1 do 
?>		.s<?=a..b?> = 0.<?
		for c=0,stDim-1 do 
		?> + RicciUL[<?=a?>][<?=c?>] * gUU->s<?=sym(c,b)?><?
	end ?>,
<?	end
end 
?>	};

	real Gaussian = 0.<? 
for a=0,stDim-1 do 
?> + RicciUL[<?=a?>][<?=a?>]<?
end ?>;

	real GammaUUL[4][4][4] = {
<? for a=0,stDim-1 do 
?>		{
<?	for b=0,stDim-1 do 
?>			{
<? 		for c=0,stDim-1 do 
?>				0.<? 
			for d=0,stDim-1 do 
?> + GammaULL->s<?=a?>.s<?=sym(d,c)?> * gUU->s<?=sym(d,b)?><?
			end ?>,
<? 		end 
?>			},
<?	 end 
?>		},
<? end 
?>	};

	//dR_ab/dg_pq
	tensor_sym4sym4 dRicciLL_dgLL = (tensor_sym4sym4){
<? for e=0,stDim-1 do 
	for f=e,stDim-1 do
?>		.s<?=e..f?> = (sym4){
<? 		for a=0,stDim-1 do
			for b=a,stDim-1 do
?>			.s<?=a..b?> = 0.<?
				for c=0,stDim-1 do ?> 
				+ GammaUUL[<?=e?>][<?=c?>][<?=c?>] * GammaULL->s<?=f?>.s<?=sym(b,a)?>
				- GammaUUL[<?=e?>][<?=c?>][<?=b?>] * GammaULL->s<?=f?>.s<?=sym(c,a)?>
				- gUU->s<?=sym(c,e)?> * RiemannULLL.s<?=f?>.s<?=a?>.s<?=sym(c,b)?><? 
				end ?>,
<? 			end
		end 
?>		},
<? 	end
end 
?>	};

	//g^cd dR_cd/dg_ef
	sym4 gUU_dRicciLL_dgLL = (sym4){
<? for e=0,stDim-1 do 
	for f=e,stDim-1 do 
?>		.s<?=e..f?> = sym4_dot(*gUU, dRicciLL_dgLL.s<?=e..f?>),
<?	end
end ?>
	};

	//dG_ab/dg_pq = tensor_sym4sym4.s[pq].s[ab]
	tensor_sym4sym4 dEinsteinLL_dgLL = (tensor_sym4sym4){
<? for e=0,stDim-1 do
	for f=e,stDim-1 do 
?>		.s<?=e..f?> = (sym4){
<? 		for a=0,stDim-1 do
			for b=a,stDim-1 do 
?>			.s<?=a..b?> = dRicciLL_dgLL.s<?=e..f?>.s<?=a..b?><?
				?> - .5 * (<?
				if e==a and f==b then ?>Gaussian<? else ?>0.<? end 
				?> + gLL->s<?=a..b?> * (<?
					?>gUU_dRicciLL_dgLL.s<?=e..f?><?
					?> - RicciUU.s<?=e..f?><?
				?>)),
<? 			end
		end
?>		},
	<? end ?>
<? end 
?>	};

	global const <?=TPrim_t?>* TPrim = TPrims + index;

	tensor_sym4sym4 d_8piTLL_dgLL = (tensor_sym4sym4){
<? for a=0,stDim-1 do
	for b=a,stDim-1 do
?>		.s<?=a..b?> = sym4_zero,
<?	end
end 
?>	};

<?
if solver.body.useEM then ?>	
	real4 EU = (real4)(0 <? for i=0,sDim-1 do ?>, TPrim->E.s<?=i?> <? end ?>); 
	real4 EL = sym4_real4_mul(*gLL, EU);
	real ESq = dot(EL, EU);
	
	real4 BU = (real4)(0 <? for i=0,sDim-1 do ?>, TPrim->B.s<?=i?> <? end ?>); 
	real4 BL = sym4_real4_mul(*gLL, BU);
	real BSq = dot(BL, BU);

	real sqrt_det_g = sqrt(fabs(sym4_det(*gLL)));
	real3 SL = real3_scale(real3_cross(TPrim->E, TPrim->B), sqrt_det_g);

<?
	for e=0,stDim-1 do
		for f=e,stDim-1 do
?>	d_8piTLL_dgLL.s<?=e..f?>.s00 += <?
			if e==0 or f==0 then 
				?>0<? 
			else 
				?>TPrim->E.s<?=e-1?> * TPrim->E.s<?=f-1?> + TPrim->B.s<?=e-1?> * TPrim->B.s<?=f-1?><? 
			end
			?>;
<? 			for i=0,sDim-1 do 
?>	d_8piTLL_dgLL.s<?=e..f?>.s0<?=i+1?> = -SL.s<?=i?> * gUU->s<?=e..f?>;
<? 				for j=i,sDim-1 do 
?>	d_8piTLL_dgLL.s<?=e..f?>.s<?=i+1?><?=j+1?> = 0.<? 
					if e==i+1 and f==j+1 then ?> + ESq + BSq <? end 
					?> + gLL->s<?=i..j?> * (<?
					if e==0 or f==0 then
						?>0<?
					else
						?>TPrim->E.s<?=e-1?> * TPrim->E.s<?=f-1?> + TPrim->B.s<?=e-1?> * TPrim->B.s<?=f-1?><?
					end
					?>) - 2. * (0.<?
					if e==i+1 then ?>
		+ EU.s<?=f?> * EL.s<?=j+1?> + BU.s<?=f?> * BL.s<?=j+1?><?
					end
					if e==j+1 then ?>
		+ EU.s<?=f?> * EL.s<?=i+1?> + BU.s<?=f?> * BL.s<?=i+1?><?
					end 
					?>);
<?
				end
			end 
		end
	end 
end
?>

<?
if solver.body.useMatter then
	if solver.body.useVel then ?>//if we're using velocity ...
	//set vU.t = 0 so we only lower by the spatial component of the metric.  right?
	real4 vU = (real4)(0, TPrim->v.x, TPrim->v.y, TPrim->v.z);
	real4 vL = sym4_real4_mul(*gLL, vU);
	real vLenSq = dot(vL, vU);	//vU.t = 0 so we'll neglect the vL.t component
	real W = 1. / sqrt(1. - sqrt(vLenSq));
	real4 uU = (real4)(W, W * vU.s1, W * vU.s2, W * vU.s3);
	real4 uL = sym4_real4_mul(*gLL, uU);
	<? else ?>//otherwise uL = gLL->s0
	real4 uL = (real4)(gLL->s00, gLL->s01, gLL->s02, gLL->s03);
	<? 
	end
	
	for e=0,stDim-1 do
		for f=e,stDim-1 do
			for a=0,stDim-1 do
				for b=a,stDim-1 do
?>	d_8piTLL_dgLL.s<?=e..f?>.s<?=a..b?> += (0. <? 
					if e==a then ?>+ uL.s<?=b?><? end 
					if e==b then ?>+ uL.s<?=a?><? end 
					?>) * uL.s<?=f?> * (TPrim->rho * (1. + TPrim->eInt) + TPrim->P)<?
					if e==a and f==b then ?> + TPrim->P<? end ?>;
<?				end
			end
		end
	end 
end
?>

	global const sym4* EFE = EFEs + index;	// G_ab - 8 pi T_ab
	sym4 dPhi_dgLL;
<? for p=0,stDim-1 do
	for q=p,stDim-1 do 
?>	dPhi_dgLL.s<?=p..q?> = sym4_dot(<?
		?>*EFE, <?
		?>sym4_sub(<?
			?>dEinsteinLL_dgLL.s<?=p..q?>, <?
			?>d_8piTLL_dgLL.s<?=p..q?><?
		?>)<?
	?>);
<? 	end
end ?>

	global gPrim_t* dPhi_dgPrim = dPhi_dgPrims + index;
	global const gPrim_t* gPrim = gPrims + index;

	sym3 gammaLL = gPrim->gammaLL;
	real3 betaU = gPrim->betaU;
	real3 betaL = sym3_real3_mul(gammaLL, betaU);

<?
if solver.convergeAlpha[0] then
?>	dPhi_dgPrim->alpha = -2. * gPrim->alpha * dPhi_dgLL.s00;
<?
end
if solver.convergeBeta[0] then
	for m=0,sDim-1 do
?>	dPhi_dgPrim->betaU.s<?=m?> = 2. * (dPhi_dgLL.s00 * betaL.s<?=m?>
<?		for n=0,sDim-1 do ?>
		+ dPhi_dgLL.s0<?=n+1?> * gammaLL.s<?=sym(n,m)?>
<? 		end ?>);
<?	end
end
if solver.convergeGamma[0] then
	for m=0,sDim-1 do
		for n=m,sDim-1 do
?>	dPhi_dgPrim->gammaLL.s<?=m..n?> = 
		betaU.s<?=m?> * (dPhi_dgLL.s00 * betaU.s<?=n?>
		+ 2. * dPhi_dgLL.s0<?=n+1?>)
		+ dPhi_dgLL.s<?=m+1?><?=n+1?>;
<?		end
	end
end ?>

	//scale up our gradient?
	//scale by c^4 / G ~ 1e+44
	// which is the units of conversion 
	//c^4/G * G_ab = 8 pi T_ab
	dPhi_dgPrim->alpha *= c*c*c*c/G;
<? for i=0,sDim-1 do
?>	dPhi_dgPrim->betaU.s<?=i?> *= c*c*c*c/G;
<? end
for i=0,sDim-1 do
	for j=i,sDim-1 do
?>	dPhi_dgPrim->gammaLL.s<?=i..j?> *= c*c*c*c/G;
<?	end
end
?>
}

kernel void update_gPrims(
	global gPrim_t* gPrims,
	global const gPrim_t* dPhi_dgPrims,
	real updateAlpha
) {
	initKernel();
	global gPrim_t* gPrim = gPrims + index;
	global const gPrim_t* dPhi_dgPrim = dPhi_dgPrims + index;

<? 
if solver.convergeAlpha[0] then
?>	gPrim->alpha -= updateAlpha * dPhi_dgPrim->alpha;
<?
end
if solver.convergeBeta[0] then
	for m=0,sDim-1 do
?>	gPrim->betaU.s<?=m?> -= updateAlpha * dPhi_dgPrim->betaU.s<?=m?>;
<? 	end 
end
if solver.convergeGamma[0] then
	for m=0,sDim-1 do
		for n=m,sDim-1 do
?>	gPrim->gammaLL.s<?=m..n?> -= updateAlpha * dPhi_dgPrim->gammaLL.s<?=m..n?>;
<?		end
	end
end
?>
}
