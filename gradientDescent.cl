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

<? if true then -- used 2nd-deriv finite-difference stencil
?>
	//g_ab,cd := dgLLL.cd.ab
<?= solver:finiteDifference2{
	bufferName = "gLLs",
	srcType = "4s4",
	resultName = "d2gLLLL",
	getBoundary = function(args) return "real4s4_Minkowski" end,
} ?>

	//TODO or store this?
	//GammaLLL.a.bc := Γ_abc = g_au Γ^u_bc
	real4x4s4 const GammaLLL = real4s4_real4x4s4_mul(gLL, GammaULL);

	//both are [ab][cd] and (ac)(bd)
	//dg_asym_LLLL.a.b.c.d := g_ad,bc - g_bd,ac - g_ac,bd + g_bc,ad
	real4x4x4x4 const dg_asym_LLLL = (real4x4x4x4){
<? for a=0,stDim-1 do
?>		.s<?=a?> = (real4x4x4){
<?	for b=0,stDim-1 do
?>			.s<?=b?> = (real4x4){
<?		for c=0,stDim-1 do
?>				.s<?=c?> = (real4)(
<?			for d=0,stDim-1 do
?>					  d2gLLLL.s<?=sym(a,d)?>.s<?=sym(b,c)?>
					+ d2gLLLL.s<?=sym(b,c)?>.s<?=sym(a,d)?>
					- d2gLLLL.s<?=sym(b,d)?>.s<?=sym(a,c)?>
					- d2gLLLL.s<?=sym(a,c)?>.s<?=sym(b,d)?>
					<?=d < 3 and "," or ""?>
<?			end
?>				),
<?		end
?>			},
<?	end
?>		},
<? end
?>	};

	//GammaSq_asym_LLLL.a.b.c.d := Γ^e_ad Γ_ebc - Γ^e_ac Γ_ebd
	real4x4x4x4 const GammaSq_asym_LLLL = (real4x4x4x4){
<? for a=0,stDim-1 do
?>		.s<?=a?> = (real4x4x4){
<?	for b=0,stDim-1 do
?>			.s<?=b?> = (real4x4){
<?		for c=0,stDim-1 do
?>				.s<?=c?> = (real4)(
<?			for d=0,stDim-1 do
?>					0.
<?				for e=0,stDim-1 do
?>					+ GammaULL.s<?=e?>.s<?=sym(a,d)?> * GammaLLL.s<?=e?>.s<?=sym(b,c)?>
					- GammaULL.s<?=e?>.s<?=sym(a,c)?> * GammaLLL.s<?=e?>.s<?=sym(b,d)?>
<?				end
?>					<?=d < 3 and "," or ""?>
<?			end
?>				),
<?		end
?>			},
<?	end
?>		},
<? end
?>	};

	//RiemannLLLL.ab.cd := R_abcd 
	//= 1/2 (g_ad,bc - g_bd,ac - g_ac,bd + g_bc,ad) + g^fg (Γ_fad Γ_gbc - Γ_fac Γ_gbd)
	//= 1/2 (g_ad,bc - g_bd,ac - g_ac,bd + g_bc,ad) + Γ^e_ad Γ_ebc - Γ^e_ac Γ_ebd
	// TODO antisymmetric storage
	// but this requires inserting -1's for reading/writing ...
	real4x4x4x4 const RiemannLLLL = real4x4x4x4_mul_add(GammaSq_asym_LLLL, dg_asym_LLLL, .5);

	//RiemannULLL.a.b.cd := R^a_bcd = 1/2 g^ae ((g_ed,cb - g_bd,ce - g_ec,bd + g_bc,de) + g^fg (Γ_fed Γ_gbc - Γ_fec Γ_gbd))
	//TODO antisymmetric storage
	real4x4x4x4 const RiemannULLL = real4s4_real4x4x4x4_mul(gUU, RiemannLLLL);

	//RicciLL.ab := R_ab = R^c_acb
	real4s4 const RicciLL = real4x4x4x4_tr13_to_real4s4(RiemannULLL);

	//RicciUL.a.b := R^a_b = g^ac R_cb
	real4x4 const RicciUL = real4s4_real4s4_mul(gUU, RicciLL);

	//RicciUU.ab := R^ab = R^a_c g^cb
	real4s4 const RicciUU = real4x4_real4s4_to_real4s4_mul(RicciUL, gUU);

	//Gaussian := R = R^a_a
	real const Gaussian = real4x4_tr(RicciUL);

	//GammaUUL.a.b.c := Γ^ab_c = Γ^a_dc g^db
	real4x4x4 const GammaUUL = real4x4s4_real4s4_mul21(GammaULL, gUU);

	/*
	... if we use g_ab,cd and don't limit h->0 our ∂/∂g_pq(x) D_c(g_ab(x')) 's ...
	dRicciLL_dgLL.pq.ab := ∂R_ab/∂g_pq
	= ∂/∂g_pq ( g^uv (g_au,bv + g_bv,au - g_ab,uv - g_uv,ab) + Γ^uv_a Γ_uvb - Γ^uv_v Γ_uab )
	= ∂/∂g_pq ( g^uv (g_au,bv + g_bv,au - g_ab,uv - g_uv,ab) + g^cd g^ef (Γ_cea Γ_dbf - Γ_cef Γ_dba))
	=	  ∂/∂g_pq g^uv (g_au,bv + g_bv,au - g_ab,uv - g_uv,ab)
		+ g^uv (∂/∂g_pq g_au,bv + ∂/∂g_pq g_bv,au - ∂/∂g_pq g_ab,uv - ∂/∂g_pq g_uv,ab)
		+ (
			  ∂/∂g_pq g^cd g^ef
			+ g^cd ∂/∂g_pq g^ef
		) (Γ_cea Γ_dbf - Γ_cef Γ_dba)
		+ g^cd g^ef (
			  ∂/∂g_pq Γ_cea Γ_dbf + Γ_cea ∂/∂g_pq Γ_dbf
			- ∂/∂g_pq Γ_cef Γ_dba - Γ_cef ∂/∂g_pq Γ_dba
		)
	=	- g^pu g^qv (g_au,bv + g_bv,au - g_ab,uv - g_uv,ab)
		+ g^uv (
			  ∂/∂g_pq D2[g_au]_bv
			+ ∂/∂g_pq D2[g_bv]_au
			- ∂/∂g_pq D2[g_ab]_uv
			- ∂/∂g_pq D2[g_uv]_ab
		)
		- (
			  g^cp g^dq g^ef
			+ g^ep g^fq g^cd
		) (Γ_cea Γ_dbf - Γ_cef Γ_dba)
		+ g^cd g^ef (
			  Γ_dbf ∂/∂g_pq Γ_cea + Γ_cea ∂/∂g_pq Γ_dbf
			- Γ_dba ∂/∂g_pq Γ_cef - Γ_cef ∂/∂g_pq Γ_dba
		)

	=	- g^cp R^q_acb
		- g^cp g^qe g^fg (Γ_efb Γ_cag - Γ_cfg Γ_eab)		<- same as the Γ^2 term in the g_ab,cd formulation of R^q_acb ...)
		+ g^uv (
			  ∂/∂g_pq D2[g_au]_bv
			+ ∂/∂g_pq D2[g_bv]_au
			- ∂/∂g_pq D2[g_ab]_uv
			- ∂/∂g_pq D2[g_uv]_ab
		)
		+ g^cd g^ef (
			  Γ_dbf ∂/∂g_pq Γ_cea + Γ_cea ∂/∂g_pq Γ_dbf
			- Γ_dba ∂/∂g_pq Γ_cef - Γ_cef ∂/∂g_pq Γ_dba
		)
	*/
	real4s4x4s4 const dRicciLL_dgLL = (real4s4x4s4){
<? for p=0,stDim-1 do
	for q=p,stDim-1 do
?>		.s<?=p..q?> = (real4s4){
<?		for a=0,stDim-1 do
			for b=a,stDim-1 do
?>			.s<?=a..b?> = 0.
<?				for c=0,stDim-1 do
?>				- gUU.s<?=sym(c,p)?> * (
					RiemannULLL.s<?=q?>.s<?=a?>.s<?=c?>.s<?=b?>

<?					for e=0,stDim-1 do
?>					+ gUU.s<?=sym(q,e)?> * GammaSq_asym_LLLL.s<?=e?>.s<?=a?>.s<?=c?>.s<?=b?>
<?					end
?>				)
				+ GammaUUL.s<?=p?>.s<?=c?>.s<?=c?> * GammaULL.s<?=q?>.s<?=sym(b,a)?>
				- GammaUUL.s<?=p?>.s<?=c?>.s<?=b?> * GammaULL.s<?=q?>.s<?=sym(c,a)?>
<?				end ?>,
<?			end
		end
?>		},
<?	end
end
?>	};


<? else -- only use 1st-deriv finite-difference stencils (and difference-of-difference for 2nd-deriv)
?>
	
	//this is also in the Ricci computation, but should I be storing it?  is it too big?
	real4x4x4s4 const dGammaLULL = calc_dGammaLULL(GammaULLs);

	//RiemannULLL.a.b.cd := R^a_bcd = Γ^a_bd,c - Γ^a_bc,d + Γ^a_ec Γ^e_bd - Γ^a_ed Γ^e_bc
	real4x4x4x4 const RiemannULLL = (real4x4x4x4){
<? for a=0,stDim-1 do
?>		.s<?=a?> = (real4x4x4){
<?	for b=0,stDim-1 do
?>			.s<?=b?> = (real4x4){
<?		for c=0,stDim-1 do
?>				.s<?=c?> = (real4)(
<?			for d=0,stDim-1 do
?>					  dGammaLULL.s<?=c?>.s<?=a?>.s<?=sym(b,d)?>
					- dGammaLULL.s<?=d?>.s<?=a?>.s<?=sym(b,c)?><?
				for e=0,stDim-1 do
?>					+ GammaULL.s<?=a?>.s<?=sym(e,c)?> * GammaULL.s<?=e?>.s<?=sym(b,d)?>
					- GammaULL.s<?=a?>.s<?=sym(e,d)?> * GammaULL.s<?=e?>.s<?=sym(b,c)?>
<?				end
?>					<?=d < 3 and "," or ""?>
<?			end
?>				),
<?		end
?>			},
<? 	end
?>		},
<? end
?>	};

	//RicciLL.ab := R_ab = R^c_acb
	real4s4 const RicciLL = real4x4x4x4_tr13_to_real4s4(RiemannULLL);

	//RicciUL.a.b := R^a_b = g^ac R_cb
	real4x4 const RicciUL = real4s4_real4s4_mul(gUU, RicciLL);

	//RicciUU.ab := R^ab = R^a_c g^cb
	real4s4 const RicciUU = real4x4_real4s4_to_real4s4_mul(RicciUL, gUU);

	//Gaussian := R = R^a_a
	real const Gaussian = real4x4_tr(RicciUL);

	//GammaUUL.a.b.c := Γ^ab_c = Γ^a_dc g^db
	real4x4x4 const GammaUUL = real4x4s4_real4s4_mul21(GammaULL, gUU);

	/*
	∂/∂g_pq δ^a_b = 0
	∂/∂g_pq(x) g_ab(x') = δ_a^p δ_b^q δ(x - x')
	∂/∂g_pq(x) g^ad(x') = -g^ap(x) g^qd(x) δ(x - x')

	D_c(g_ab(x)) = C_m g_ab(x + m dx^c) ≈ g_ab,c
	∂/∂g_pq(x) D_c(g_ab(x'))
		= ∂/∂g_pq(g_ab) C_(x'^c - x^c)
		= δ_a^p δ_b^q C_(x'^c - x^c)

	... where C_(x'^c - x^c) is the coefficient of the derivative based on (x'^c - x^c) / h, or 0 if any of the other coefficients don't match

	∂/∂g_pq(x) Γ_abc(x')
		= 1/2 ∂/∂g_pq(x) (D_c(g_ab(x')) + D_b(g_ac(x')) - D_a(g_bc(x')))
		= 1/2 ∂/∂g_pq(x) (
			  δ_a^p δ_b^q C_(x'^c - x^c)
			+ δ_a^p δ_c^q C_(x'^b - x^b)
			- δ_b^p δ_c^q C_(x'^a - x^a)
		)

	∂/∂g_pq Γ^a_bc
	= -g^ap Γ^q_bc based on the limit, however for finite-difference it has another value ...
	∂/∂g_pq(x) Γ^a_bc(x')
	= -g^ap(x') Γ^q_bc(x') δ(x - x') + g^ad(x') ∂/∂g_pq(x) Γ_dbc(x')


	dRicciLL_dgLL.pq.ab := ∂R_ab/∂g_pq
	= Γ^pc_c Γ^q_ba - Γ^pc_b Γ^q_ca - g^cp R^q_acb
	(work done in efe.html)
	*/
	real4s4x4s4 const dRicciLL_dgLL = (real4s4x4s4){
<?
for p=0,stDim-1 do
	for q=p,stDim-1 do
?>		.s<?=p..q?> = (real4s4){
<?		for a=0,stDim-1 do
			for b=a,stDim-1 do
?>			.s<?=a..b?> = 0.
<?				for c=0,stDim-1 do ?>
				+ GammaUUL.s<?=p?>.s<?=c?>.s<?=c?> * GammaULL.s<?=q?>.s<?=sym(b,a)?>
				- GammaUUL.s<?=p?>.s<?=c?>.s<?=b?> * GammaULL.s<?=q?>.s<?=sym(c,a)?>
				- gUU.s<?=sym(c,p)?> * RiemannULLL.s<?=q?>.s<?=a?>.s<?=c?>.s<?=b?><?
				end ?>,
<?			end
		end
?>		},
<?	end
end
?>	};

<? end ?>

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
<?		for a=0,stDim-1 do
			for b=a,stDim-1 do
?>			.s<?=a..b?> = dRicciLL_dgLL.s<?=p..q?>.s<?=a..b?><?
				?> - .5 * (<?
				if p==a and q==b then ?>Gaussian<? else ?>0.<? end
				?> + gLL.s<?=a..b?> * (<?
					?>gUU_dRicciLL_dgLL.s<?=p..q?><?
					?> - RicciUU.s<?=p..q?><?
				?>)),
<?			end
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
<?			for i=0,sDim-1 do
?>	d_8piTLL_dgLL.s<?=e..f?>.s0<?=i+1?> = -SL.s<?=i?> * gUU.s<?=e..f?>;
<?				for j=i,sDim-1 do
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
<?	end
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
<?		end ?>);
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
<?	end
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
