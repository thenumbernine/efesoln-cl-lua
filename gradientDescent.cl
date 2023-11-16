// here is the code used for updating weights based on the EFE

kernel void calc_dPhi_dgPrims(
	global gPrim_t * const dPhi_dgPrims,
	global <?=TPrim_t?> const * const TPrims,
	global gPrim_t const * const gPrims,
	global real4s4 const * const gLLs,
	global real4s4 const * const gUUs,
	global real4x4s4 const * const GammaULLs,
	global real4s4 const * const EFEs
) {
	initKernel();
	
	real4s4 const gLL = gLLs[index];
	real4s4 const gUU = gUUs[index];
	real4x4s4 const GammaULL = GammaULLs[index];

	//this is also in the Ricci computation, but should I be storing it?  is it too big?
	real4x4x4s4 const dGammaLULL = calc_dGammaLULL(GammaULLs);

	//RiemannULLL.a.b.cd := R^a_bcd = Γ^a_bd,c - Γ^a_bc,d + Γ^a_ec Γ^e_bd - Γ^a_ed Γ^e_bc
	real4x4x4s4 const RiemannULLL = (real4x4x4s4){
<? for a=0,stDim-1 do
?>		.s<?=a?> = (real4x4s4){
<? 	for b=0,stDim-1 do
?>			.s<?=b?> = (real4s4){
<? 		for c=0,stDim-1 do 
			for d=c,stDim-1 do 
?>				.s<?=c..d?> = 
					  dGammaLULL.s<?=c?>.s<?=a?>.s<?=sym(b,d)?>
					- dGammaLULL.s<?=d?>.s<?=a?>.s<?=sym(b,c)?><?
				for e=0,stDim-1 do ?>
					+ GammaULL.s<?=a?>.s<?=sym(e,c)?> * GammaULL.s<?=e?>.s<?=sym(b,d)?>
					- GammaULL.s<?=a?>.s<?=sym(e,d)?> * GammaULL.s<?=e?>.s<?=sym(b,c)?><? 
				end ?>,
<? 			end 
		end
?>			},
<? end
?>		},
<? end 
?>	};

	//RicciLL.ab := R_ab = R^c_acb
	real4s4 const RicciLL = real4x4x4s4_tr13(RiemannULLL);

	//RicciUL.a.b := R^a_b = g^ac R_cb
	real4x4 const RicciUL = real4s4_real4s4_mul(gUU, RicciLL);

	//RicciUU.ab := R^ab = R^a_c g^cb
	real4s4 const RicciUU = real4x4_real4s4_to_real4s4_mul(RicciUL, gUU);

	//Gaussian := R = R^a_a
	real const Gaussian = real4x4_tr(RicciUL);

	//GammaUUL.a.b.c := Γ^ab_c = Γ^a_dc g^db
	real4x4x4 const GammaUUL = real4x4s4_real4s4_mul21(GammaULL, gUU);

	/*
	∂g_ab/∂g_pq = δ_a^p δ_b^q

	∂/∂g_pq δ^a_b = 0
	∂/∂g_pq (g^ac g_cb) = 0
	∂/∂g_pq g^ac g_cb + g^ac ∂/∂g_pq g_cb = 0
	∂/∂g_pq g^ac g_cb g^bd = -g^ac ∂/∂g_pq g_cb g^bd
	∂/∂g_pq g^ad = -g^ac ∂/∂g_pq g_cb g^bd
	∂/∂g_pq g^ad = -g^ac ∂_c^p ∂_b^q g^bd
	∂/∂g_pq g^ad = -g^ap g^qd

	∂/∂g_pq Γ^a_bc 
	= 1/2 ∂/∂g_pq g^au (g_ub,c + g_uc,b - g_bc,u)
	= 1/2 (∂/∂g_pq g^au) (g_ub,c + g_uc,b - g_bc,u)
		+ 1/2 g^au (∂/∂g_pq g_ub,c + ∂/∂g_pq g_uc,b - ∂/∂g_pq g_bc,u)
	= 1/2 (∂/∂g_pq g^au) (g_ub,c + g_uc,b - g_bc,u)
		+ 1/2 g^au ((δ_u^p δ_b^q),c + (δ_u^p δ_c^q),b - (δ_b^p δ_c^q),u)
	= -1/2 (g^ae ∂/∂g_pq g_ef g^fu) (g_ub,c + g_uc,b - g_bc,u)
	= -1/2 g^ap g^qu (g_ub,c + g_uc,b - g_bc,u)
	= -g^ap Γ^q_bc

	dRicciLL_dgLL.pq.ab := ∂R_ab/∂g_pq
	= ∂/∂g_pq R^c_acb
	= ∂/∂g_pq (Γ^c_ab,c - Γ^c_ac,b + Γ^c_ec Γ^e_ab - Γ^c_eb Γ^e_ac)
	= ∂/∂g_pq Γ^c_ab,c - ∂/∂g_pq Γ^c_ac,b 
		+ ∂/∂g_pq Γ^c_ec Γ^e_ab 
		+ Γ^c_ec ∂/∂g_pq Γ^e_ab 
		- ∂/∂g_pq Γ^c_eb Γ^e_ac
		- Γ^c_eb ∂/∂g_pq Γ^e_ac
	= -(g^cp Γ^q_ab),c + (g^cp Γ^q_ac),b 
		- g^cp Γ^q_ec Γ^e_ab 
		- g^ep Γ^q_ab Γ^c_ec
		+ g^cp Γ^q_eb Γ^e_ac
		+ g^ep Γ^q_ac Γ^c_eb
	= -g^cp_,c Γ^q_ab 
		- g^cp Γ^q_ab,c
		+ g^cp_,b Γ^q_ac 
		+ g^cp Γ^q_ac,b
		- Γ^e_ab Γ^qp_e
		- Γ^q_ab Γ^cp_c
		+ Γ^q_eb Γ^ep_a
		+ Γ^q_ac Γ^cp_b
	= -g^cp_,c Γ^q_ab 
		+ g^cp_,b Γ^q_ac 
		- Γ^q_ab Γ^cp_c
		+ Γ^q_ac Γ^cp_b
		
		- g^cp Γ^q_ab,c
		+ g^cp Γ^q_ac,b
		- g^cp Γ^q_ec Γ^e_ab 
		+ g^cp Γ^q_eb Γ^e_ac
	= g^cu g_uv,c g^vp Γ^q_ab 
		- g^cu g_uv,b g^vp Γ^q_ac 
		- Γ^q_ab Γ^cp_c
		+ Γ^q_ac Γ^cp_b
		
		- g^cp R^q_acb
	
	... using g_ab,c = Γ_acb + Γ_bca

	= g^cu (Γ_ucv + Γ_vcu) g^vp Γ^q_ab 
		- g^cu (Γ_ubv + Γ_vbu) g^vp Γ^q_ac 
		- Γ^q_ab Γ^cp_c
		+ Γ^q_ac Γ^cp_b
		
		- g^cp R^q_acb

	= Γ^cp_c Γ^q_ab 
		+ Γ^pc_c Γ^q_ab
		- Γ^cp_b Γ^q_ac
		- Γ^pc_b Γ^q_ac
		- Γ^q_ab Γ^cp_c
		+ Γ^q_ac Γ^cp_b
		
		- g^cp R^q_acb
	
	= Γ^pc_c Γ^q_ba - Γ^pc_b Γ^q_ca - g^cp R^q_acb
	*/
	real4s4x4s4 const dRicciLL_dgLL = (real4s4x4s4){
<? for p=0,stDim-1 do 
	for q=p,stDim-1 do
?>		.s<?=p..q?> = (real4s4){
<? 		for a=0,stDim-1 do
			for b=a,stDim-1 do
?>			.s<?=a..b?> = 0.<?
				for c=0,stDim-1 do ?> 
				+ GammaUUL.s<?=p?>.s<?=c?>.s<?=c?> * GammaULL.s<?=q?>.s<?=sym(b,a)?>
				- GammaUUL.s<?=p?>.s<?=c?>.s<?=b?> * GammaULL.s<?=q?>.s<?=sym(c,a)?>
				- gUU.s<?=sym(c,p)?> * RiemannULLL.s<?=q?>.s<?=a?>.s<?=sym(c,b)?><? 
				end ?>,
<? 			end
		end 
?>		},
<? 	end
end 
?>	};

	//g^ab ∂R_ab/∂g_pq
	real4s4 const gUU_dRicciLL_dgLL = (real4s4){
<? for p=0,stDim-1 do 
	for q=p,stDim-1 do 
?>		.s<?=p..q?> = real4s4_dot(gUU, dRicciLL_dgLL.s<?=p..q?>),
<?	end
end ?>
	};

	/*
	∂G_ab/∂g_pq = dEinsteinLL_dgLL.pq.ab
	= ∂/∂g_pq (R_ab - 1/2 R g_ab)
	= ∂/∂g_pq R_ab - 1/2 ∂/∂g_pq R g_ab - 1/2 R ∂/∂g_pq g_ab
	= ∂/∂g_pq R_ab - 1/2 g_ab ∂/∂g_pq (R_uv g^uv) - 1/2 R δ_a^p δ_b^q
	= ∂/∂g_pq R_ab - 1/2 g_ab (∂/∂g_pq R_uv g^uv + R_uv ∂/∂g_pq g^uv) - 1/2 R δ_a^p δ_b^q
	= ∂/∂g_pq R_ab - 1/2 g_ab (g^uv ∂/∂g_pq R_uv - R_uv g^pu g^qv) - 1/2 R δ_a^p δ_b^q
	= ∂R_ab/∂g_pq - 1/2 (R δ_a^p δ_b^q + g_ab (g^uv ∂R_uv/∂g_pq - R^pq))
	*/
	real4s4x4s4 const dEinsteinLL_dgLL = (real4s4x4s4){
<? for p=0,stDim-1 do
	for q=p,stDim-1 do 
?>		.s<?=p..q?> = (real4s4){
<? 		for a=0,stDim-1 do
			for b=a,stDim-1 do 
?>			.s<?=a..b?> = dRicciLL_dgLL.s<?=p..q?>.s<?=a..b?><?
				?> - .5 * (<?
				if p==a and q==b then ?>Gaussian<? else ?>0.<? end 
				?> + gLL.s<?=a..b?> * (<?
					?>gUU_dRicciLL_dgLL.s<?=p..q?><?
					?> - RicciUU.s<?=p..q?><?
				?>)),
<? 			end
		end
?>		},
	<? end ?>
<? end 
?>	};

	<?=TPrim_t?> const TPrim = TPrims[index];

	real4s4x4s4 d_8piTLL_dgLL = real4s4x4s4_zero;
<?
if solver.body.useEM then ?>	
	real4 const EU = (real4)(0 <? for i=0,sDim-1 do ?>, TPrim.E.s<?=i?> <? end ?>); 
	real4 const EL = real4s4_real4_mul(gLL, EU);
	real const ESq = dot(EL, EU);
	
	real4 const BU = (real4)(0 <? for i=0,sDim-1 do ?>, TPrim.B.s<?=i?> <? end ?>); 
	real4 const BL = real4s4_real4_mul(gLL, BU);
	real const BSq = dot(BL, BU);

	real const sqrt_det_g = sqrt(fabs(real4s4_det(gLL)));
	real3 const SL = real3_real_mul(real3_cross(TPrim.E, TPrim.B), sqrt_det_g);

<?
	for e=0,stDim-1 do
		for f=e,stDim-1 do
?>	d_8piTLL_dgLL.s<?=e..f?>.s00 += <?
			if e==0 or f==0 then 
				?>0<? 
			else 
				?>TPrim.E.s<?=e-1?> * TPrim.E.s<?=f-1?> + TPrim.B.s<?=e-1?> * TPrim.B.s<?=f-1?><? 
			end
			?>;
<? 			for i=0,sDim-1 do 
?>	d_8piTLL_dgLL.s<?=e..f?>.s0<?=i+1?> = -SL.s<?=i?> * gUU.s<?=e..f?>;
<? 				for j=i,sDim-1 do 
?>	d_8piTLL_dgLL.s<?=e..f?>.s<?=i+1?><?=j+1?> = 0.<? 
					if e==i+1 and f==j+1 then ?> + ESq + BSq <? end 
					?> + gLL.s<?=i..j?> * (<?
					if e==0 or f==0 then
						?>0<?
					else
						?>TPrim.E.s<?=e-1?> * TPrim.E.s<?=f-1?> + TPrim.B.s<?=e-1?> * TPrim.B.s<?=f-1?><?
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
	real4 const vU = (real4)(0, TPrim.v.x, TPrim.v.y, TPrim.v.z);
	real4 const vL = real4s4_real4_mul(gLL, vU);
	real const vLenSq = dot(vL, vU);	//vU.t = 0 so we'll neglect the vL.t component
	real const W = 1. / sqrt(1. - sqrt(vLenSq));
	real4 const uU = (real4)(W, W * vU.s1, W * vU.s2, W * vU.s3);
	real4 const uL = real4s4_real4_mul(gLL, uU);
	<? else ?>//otherwise uL = gLL.s0
	real4 const uL = (real4)(gLL.s00, gLL.s01, gLL.s02, gLL.s03);
	<? 
	end
	
	for e=0,stDim-1 do
		for f=e,stDim-1 do
			for a=0,stDim-1 do
				for b=a,stDim-1 do
?>	d_8piTLL_dgLL.s<?=e..f?>.s<?=a..b?> += (0. <? 
					if e==a then ?>+ uL.s<?=b?><? end 
					if e==b then ?>+ uL.s<?=a?><? end 
					?>) * uL.s<?=f?> * (TPrim.rho * (1. + TPrim.eInt) + TPrim.P)<?
					if e==a and f==b then ?> + TPrim.P<? end ?>;
<?				end
			end
		end
	end 
end
?>

	real4s4 const EFE = EFEs[index];	// G_ab - 8 π T_ab
	
	//∂Φ/∂g_pq = (G_ab - 8 π T_ab) (∂G_ab/∂g_pq - 8 π ∂T_ab/∂g_pq)
	real4s4 dPhi_dgLL = (real4s4){
<? for p=0,stDim-1 do
	for q=p,stDim-1 do 
?>	.s<?=p..q?> = real4s4_dot(
			EFE, 
			real4s4_sub(
				dEinsteinLL_dgLL.s<?=p..q?>, 
				d_8piTLL_dgLL.s<?=p..q?>
			)
		),
<? 	end
end
?>	};

	global gPrim_t * const dPhi_dgPrim = dPhi_dgPrims + index;
	gPrim_t const gPrim = gPrims[index];

	real3s3 const gammaLL = gPrim.gammaLL;
	real3 const betaU = gPrim.betaU;
	real3 const betaL = real3s3_real3_mul(gammaLL, betaU);

<?
if solver.convergeAlpha then
?>	dPhi_dgPrim->alpha = -2. * gPrim.alpha * dPhi_dgLL.s00;
<?
end
if solver.convergeBeta then
	for m=0,sDim-1 do
?>	dPhi_dgPrim->betaU.s<?=m?> = 2. * (dPhi_dgLL.s00 * betaL.s<?=m?>
<?		for n=0,sDim-1 do ?>
		+ dPhi_dgLL.s0<?=n+1?> * gammaLL.s<?=sym(n,m)?>
<? 		end ?>);
<?	end
end
if solver.convergeGamma then
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
	//c^4/G * G_ab = 8 π T_ab
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
	global gPrim_t * const gPrims,
	global gPrim_t const * const dPhi_dgPrims,
	real const updateLambda
) {
	initKernel();
	global gPrim_t * const gPrim = gPrims + index;
	gPrim_t const dPhi_dgPrim = dPhi_dgPrims[index];

<? 
if solver.convergeAlpha then
?>	gPrim->alpha -= updateLambda * dPhi_dgPrim.alpha;
<?
end
if solver.convergeBeta then
	for m=0,sDim-1 do
?>	gPrim->betaU.s<?=m?> -= updateLambda * dPhi_dgPrim.betaU.s<?=m?>;
<? 	end 
end
if solver.convergeGamma then
	for m=0,sDim-1 do
		for n=m,sDim-1 do
?>	gPrim->gammaLL.s<?=m..n?> -= updateLambda * dPhi_dgPrim.gammaLL.s<?=m..n?>;
<?		end
	end
end
?>
}
