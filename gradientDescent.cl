// here is the code used for updating weights based on the EFE

kernel void calc_dPhi_dgLLs(
	global sym4* dPhi_dgLLs,
	global const TPrim_t* TPrims,
	global const sym4* gLLs,
	global const sym4* gUUs,
	global const tensor_4sym4* GammaULLs,
	global const sym4* EFEs
) {
	INIT_KERNEL();
	
	global const sym4* gLL = gLLs + index;
	global const sym4* gUU = gUUs + index;
	global const tensor_4sym4* GammaULL = GammaULLs + index;

	//this is also in the Ricci computation, but should I be storing it?  is it too big?
	tensor_44sym4 dGammaLULL;
	calc_dGammaLULL(&dGammaLULL, GammaULLs);

	tensor_44sym4 RiemannULLL = (tensor_44sym4){
<? for a=0,dim-1 do
?>		.s<?=a?> = (tensor_4sym4){
<? 	for b=0,dim-1 do
?>			.s<?=b?> = (sym4){
<? 		for c=0,dim-1 do 
			for d=c,dim-1 do 
?>				.s<?=c..d?> = dGammaLULL.s<?=c?>.s<?=a?>.s<?=sym(b,d)?>
					- dGammaLULL.s<?=d?>.s<?=a?>.s<?=sym(b,c)?><?
				for e=0,dim-1 do ?>
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
<? for a=0,dim-1 do 
	for b=a,dim-1 do
?>		.s<?=a..b?> = 0. <?
		for c=0,dim-1 do 
?> + RiemannULLL.s<?=c?>.s<?=a?>.s<?=sym(c,b)?><? 
		end ?>,
<? 	end
end 
?>	};

	real RicciUL[4][4] = {
<? for a=0,dim-1 do 
?>		{
<?
	for b=0,dim-1 do 
?>			0.<?
		for c=0,dim-1 do 
?> + gUU->s<?=sym(a,c)?> * RicciLL.s<?=sym(c,b)?><?
		end ?>,
<?
	end 
?>		},
<? end 
?>	};

	sym4 RicciUU = (sym4){
<? for a=0,dim-1 do
	for b=a,dim-1 do 
?>		.s<?=a..b?> = 0.<?
		for c=0,dim-1 do 
		?> + RicciUL[<?=a?>][<?=c?>] * gUU->s<?=sym(c,b)?><?
	end ?>,
<?	end
end 
?>	};

	real Gaussian = 0.<? 
for a=0,dim-1 do 
?> + RicciUL[<?=a?>][<?=a?>]<?
end ?>;

	real GammaUUL[4][4][4] = {
<? for a=0,dim-1 do 
?>		{
<?	for b=0,dim-1 do 
?>			{
<? 		for c=0,dim-1 do 
?>				0.<? 
			for d=0,dim-1 do 
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
<? for e=0,dim-1 do 
	for f=e,dim-1 do
?>		.s<?=e..f?> = (sym4){
<? 		for a=0,dim-1 do
			for b=a,dim-1 do
?>			.s<?=a..b?> = 0.<?
				for c=0,dim-1 do ?> 
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
<? for e=0,dim-1 do 
	for f=e,dim-1 do 
?>		.s<?=e..f?> = sym4_dot(*gUU, dRicciLL_dgLL.s<?=e..f?>),
<?	end
end ?>
	};

	//dG_ab/dg_pq = tensor_sym4sym4.s[pq].s[ab]
	tensor_sym4sym4 dEinsteinLL_dgLL = (tensor_sym4sym4){
<? for e=0,dim-1 do
	for f=e,dim-1 do 
?>		.s<?=e..f?> = (sym4){
<? 		for a=0,dim-1 do
			for b=a,dim-1 do 
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

	global const TPrim_t* TPrim = TPrims + index;

	real4 EU = (real4)(0 <?for i=0,2 do ?>, TPrim->E.s<?=i?> <? end ?>); 
	real4 EL = sym4_real4_mul(*gLL, EU);
	real ESq = dot(EL, EU);
	
	real4 BU = (real4)(0 <?for i=0,2 do ?>, TPrim->B.s<?=i?> <? end ?>); 
	real4 BL = sym4_real4_mul(*gLL, BU);
	real BSq = dot(BL, BU);

	real sqrt_det_g = sqrt(fabs(sym4_det(*gLL)));
	real3 SL = real3_scale(real3_cross(TPrim->E, TPrim->B), sqrt_det_g);

	tensor_sym4sym4 d_8piT_EM_LL_dgLL = (tensor_sym4sym4){
<? for e=0,dim-1 do
	for f=e,dim-1 do
?>		.s<?=e..f?> = (sym4){
			//dTEM_00/dg_<?=e..f?>
			.s00 = <?
		if e==0 or f==0 then 
			?>0<? 
		else 
			?>TPrim->E.s<?=e-1?> * TPrim->E.s<?=f-1?> + TPrim->B.s<?=e-1?> * TPrim->B.s<?=f-1?><? 
		end
			?>,
<? 		for i=0,subDim-1 do 
?>			.s0<?=i+1?> = -SL.s<?=i?> * gUU->s<?=e..f?>,
<? 			for j=i,subDim-1 do 
?>			.s<?=i+1?><?=j+1?> = 0.<? 
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
				?>),
<? 			end
		end 
?>		},
<? 	end
end 
?>	};

	//if we're using matter ...

		//if we're using velocity ...

		//otherwise uL = gLL->s0
	real4 uL = (real4)(gLL->s00, gLL->s01, gLL->s02, gLL->s03);

	tensor_sym4sym4 d_8piT_matter_LL_dgLL = (tensor_sym4sym4){
<? for e=0,dim-1 do
	for f=e,dim-1 do
?>
		.s<?=e..f?> = (sym4){
<?		for a=0,dim-1 do
			for b=a,dim-1 do
?>			.s<?=a..b?> = (0. <? 
			if e==a then ?>+ uL.s<?=b?><? end 
			if e==b then ?>+ uL.s<?=a?><? end 
			?>) * uL.s<?=f?> * (TPrim->rho * (1. + TPrim->eInt) + TPrim->P)<?
			if e==a and f==b then ?> + TPrim->P<? end ?>,
<?			end
		end
?>		},
<?	end
end 
?>	};

	tensor_sym4sym4 d_8piTLL_dgLL = tensor_sym4sym4_add(d_8piT_EM_LL_dgLL, d_8piT_matter_LL_dgLL);
	
	global const sym4* EFE = EFEs + index;	// G_ab - 8 pi T_ab
<? for p=0,dim-1 do
	for q=p,dim-1 do 
?>	dPhi_dgLLs[index].s<?=p..q?> = sym4_dot(<?
		?>*EFE, <?
		?>sym4_sub(<?
			?>dEinsteinLL_dgLL.s<?=p..q?>, <?
			?>d_8piTLL_dgLL.s<?=p..q?><?
		?>)<?
	?>);
<? 	end
end ?>
}

kernel void update_dgLLs(
	global sym4* gLLs,
	global const sym4* dPhi_dgLLs
) {
	INIT_KERNEL();
	global sym4* gLL = gLLs + index;
	global const sym4* dPhi_dgLL = dPhi_dgLLs + index;
<? for a=0,dim-1 do
	for b=a,dim-1 do
?>	gLL->s<?=a..b?> -= dPhi_dgLL->s<?=a..b?><? 
if updateAlpha~= 1 then ?>* <?=updateAlpha?><? end 
?>;
<?	end
end 
?>}
