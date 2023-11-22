// here is the code used for updating weights based on the EFE

real4s4x4s4 calc_partial_gLL_of_8piTLL(
	<?=TPrim_t?> const TPrim,
	real4s4 const gLL,
	real4s4 const gUU
) {

	real4s4x4s4 partial_gLL_of_8piTLL = real4s4x4s4_zero;
<?
if solver.body.useEM then ?>
	real4 const EU = real3_to_real4(TPrim.E);
	real4 const EL = real4s4_real4_mul(gLL, EU);
	real const ESq = real4_dot(EL, EU);

	real4 const BU = real3_to_real4(TPrim.B);
	real4 const BL = real4s4_real4_mul(gLL, BU);
	real const BSq = real4_dot(BL, BU);

	real const sqrt_det_g = sqrt(fabs(real4s4_det(gLL)));
	real3 const SL = real3_real_mul(real3_cross(TPrim.E, TPrim.B), sqrt_det_g);

	for (int e = 0; e < stDim; ++e) {
		for (int f = e; f < stDim; ++f) {
			int const ef = sym4[e][f];
			if (e > 0 && f > 0) {
				partial_gLL_of_8piTLL.s[ef].s00 += TPrim.E.s[e-1] * TPrim.E.s[f-1] + TPrim.B.s[e-1] * TPrim.B.s[f-1];
			}
			for (int i = 0; i < sDim; ++i) {
				partial_gLL_of_8piTLL.s[ef].s[sym4[0][i+1]] -= SL.s[i] * gUU.s[ef];
				for (int j = i; j < sDim; ++j) {

					real sum = 0;
					if (e == i+1 && f == j+1) sum += ESq + BSq;

					if (e > 0 && f > 0) {
						sum += gLL.s[sym4[i+1][j+1]] * (TPrim.E.s[e-1] * TPrim.E.s[f-1] + TPrim.B.s[e-1] * TPrim.B.s[f-1]);
					}
					if (e == i+1) {
						sum -= 2. * (EU.s[f] * EL.s[j+1] + BU.s[f] * BL.s[j+1]);
					}
					if (e == j+1) {
						sum -= 2. * (EU.s[f] * EL.s[i+1] + BU.s[f] * BL.s[i+1]);
					}

					partial_gLL_of_8piTLL.s[ef].s[sym4[i+1][j+1]] += sum;
				}
			}
		}
	}
<?
end

if solver.body.useMatter then
	if solver.body.useVel then -- if we're using velocity ...
?>
	//set vU.t = 0 so we only lower by the spatial component of the metric.  right?
	real4 const vU = real3_to_real4(TPrim.v);
	real4 const vL = real4s4_real4_mul(gLL, vU);
	real const vLenSq = real4_dot(vL, vU);	//vU.t = 0 so we'll neglect the vL.t component
	real const W = 1. / sqrt(1. - sqrt(vLenSq));
	real4 const uU = (real4){.s={W, W * vU.s1, W * vU.s2, W * vU.s3}};
	real4 const uL = real4s4_real4_mul(gLL, uU);
<?
	else -- otherwise uL = gLL.s0
?>
	real4 const uL = (real4){.s={gLL.s00, gLL.s01, gLL.s02, gLL.s03}};
<?
	end
?>
	for (int e = 0; e < stDim; ++e) {
		for (int f = e; f < stDim; ++f) {
			int const ef = sym4[e][f];
			for (int a = 0; a < stDim; ++a) {
				for (int b = a; b < stDim; ++b) {
					int const ab = sym4[a][b];
					real sum = 0;
					if (e == a) sum += uL.s[b];
					if (e == b) sum += uL.s[a];
					sum *= uL.s[f] * (TPrim.rho * (1. + TPrim.eInt) + TPrim.P);
					if (ef == ab) sum += TPrim.P;
					partial_gLL_of_8piTLL.s[ef].s[ab] += sum;
				}
			}
		}
	}
<?
end
?>
	return partial_gLL_of_8piTLL;
}

kernel void calc_partial_gPrim_of_Phis(
	global gPrim_t * const partial_gPrim_of_Phis,
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

<? if false then -- used 2nd-deriv finite-difference stencil
?>
	//partial_xU2_of_gLL.cd.ab := g_ab,cd
<?= solver:finiteDifference2{
	bufferName = "gLLs",
	srcType = "4s4",
	resultName = "partial_xU2_of_gLL",
	getBoundary = function(args) return "real4s4_Minkowski" end,
} ?>

	//TODO or store this?
	//GammaLLL.a.bc := Γ_abc = g_au Γ^u_bc
	real4x4s4 const GammaLLL = real4s4_real4x4s4_mul(gLL, GammaULL);

	//partial_xU2_of_gLL_asym.a.b.c.d := g_ad,bc - g_bd,ac - g_ac,bd + g_bc,ad
	// TODO antisymmetric storage
	//both are T_abcd = -T_bacd = -T_abdc and T_abcd = T_cdab
	real4x4x4x4 partial_xU2_of_gLL_asym;
	for (int a = 0; a < stDim; ++a) {
		for (int b = 0; b < stDim; ++b) {
			for (int c = 0; c < stDim; ++c) {
				for (int d = 0; d < stDim; ++d) {
					partial_xU2_of_gLL_asym.s[a].s[b].s[c].s[d] =
						partial_xU2_of_gLL.s[sym4[a][d]].s[sym4[b][c]]
						+ partial_xU2_of_gLL.s[sym4[b][c]].s[sym4[a][d]]
						- partial_xU2_of_gLL.s[sym4[b][d]].s[sym4[a][c]]
						- partial_xU2_of_gLL.s[sym4[a][c]].s[sym4[b][d]];
				}
			}
		}
	}

	//GammaSq_asym_LLLL.a.b.c.d := Γ^e_ad Γ_ebc - Γ^e_ac Γ_ebd
	// not the same as the popular 2 Γ^a_e[c Γ^e_d]b used for Riemann with 2 ∂/dx^[c Γ^a_d]b
	// instead this is the one used with 2 g_[a[b,c]d]
	// TODO antisymmetric storage
	real4x4x4x4 GammaSq_asym_LLLL;
	for (int a = 0; a < stDim; ++a) {
		for (int b = 0; b < stDim; ++b) {
			for (int c = 0; c < stDim; ++c) {
				for (int d = 0; d < stDim; ++d) {
					real sum = 0.;
					for (int e = 0; e < stDim; ++e) {
						sum += GammaULL.s[e].s[sym4[a][d]] * GammaLLL.s[e].s[sym4[b][c]]
							- GammaULL.s[e].s[sym4[a][c]] * GammaLLL.s[e].s[sym4[b][d]];
					}
					GammaSq_asym_LLLL.s[a].s[b].s[c].s[d] = sum;
				}
			}
		}
	}

	// linear transform, not necessarily tensoral (since neither is anyways)
	real4x4x4x4 const gUU_times_partial_xU2_gLL_asym = real4x4x4x4_real4s4_mul_1_1(partial_xU2_of_gLL_asym, gUU);
	real4x4x4x4 const GammaSq_asym_ULLL = real4x4x4x4_real4s4_mul_1_1(GammaSq_asym_LLLL, gUU);

	//RiemannLLLL.a.b.c.d := R_abcd
	//= 1/2 (g_ad,bc - g_bd,ac - g_ac,bd + g_bc,ad) + g^fg (Γ_fad Γ_gbc - Γ_fac Γ_gbd)
	//= 1/2 (g_ad,bc - g_bd,ac - g_ac,bd + g_bc,ad) + Γ^e_ad Γ_ebc - Γ^e_ac Γ_ebd
	// TODO antisymmetric storage
	// but this requires inserting -1's for reading/writing ...

	//RiemannULLL.a.b.c.d := R^a_bcd = 1/2 g^ae ((g_ed,cb - g_bd,ce - g_ec,bd + g_bc,de) + g^fg (Γ_fed Γ_gbc - Γ_fec Γ_gbd))
	//TODO antisymmetric storage
	real4x4x4x4 const RiemannULLL = real4x4x4x4_mul_add(GammaSq_asym_ULLL, gUU_times_partial_xU2_gLL_asym, .5);

	//RiemannULUL.a.b.c.d = R^a_b^c_d = R^a_bud g^uc
	real4x4x4x4 const RiemannULUL = real4x4x4x4_real4s4_mul_3_1(RiemannULLL, gUU);
	
	real4x4x4x4 const GammaSq_asym_ULUL = real4x4x4x4_real4s4_mul_3_1(GammaSq_asym_ULLL, gUU);

	//RicciLL.ab := R_ab = R^c_acb
	real4s4 const RicciLL = real4x4x4x4_tr13_to_real4s4(RiemannULLL);

	//RicciUL.a.b := R^a_b = g^ac R_cb
	real4x4 const RicciUL = real4s4_real4s4_mul(gUU, RicciLL);

	//RicciUU.ab := R^ab = R^a_c g^cb
	real4s4 const RicciUU = real4x4_real4s4_to_real4s4_mul(RicciUL, gUU);

	//Gaussian := R = R^a_a
	real const Gaussian = real4x4_tr(RicciUL);

<?
local math = require 'ext.math'
local order = solver.diffOrder
local d1coeffs = assert(derivCoeffs[1][order])
local d2coeffs = assert(derivCoeffs[2][order], "couldn't find d2 coeffs for order "..order)
?>

	// TODO
	real4s4x4x4s4 partial_gLL_of_DgLLL;
#if 0	
	for (int p = 0; p < stDim; ++p) {
		for (int q = p; q < stDim; ++q) {
			int const pq = sym4[p][q];
			for (int a = 0; a < stDim; ++a) {
				for (int b = 0; b < stDim; ++b) {
					for (int c = b; c < stDim; ++c) {
						int const bc = sym4[b][c];
						partial_gLL_of_DgLLL.s[pq].s[a].s[bc] = 
					}
				}
			}
		}
	}
#endif	
	real4s4x4s4x4s4 partial_gLL_of_D2gLLLL;

	//partial_gLL_of_GammaLLL.pq.a.bc =: ∂/∂g_pq(x) Γ_abc
	real4s4x4x4s4 partial_gLL_of_GammaLLL = real4s4x4x4s4_zero;
	for (int p = 0; p < stDim; ++p) {
		for (int q = p; q < stDim; ++q) {
			int const pq = sym4[p][q];
			for (int a = 0; a < stDim; ++a) {
				for (int b = 0; b < stDim; ++b) {
					for (int c = b; c < stDim; ++c) {
						int const bc = sym4[b][c];
						partial_gLL_of_GammaLLL.s[pq].s[a].s[bc] = .5 * (
							  partial_gLL_of_DgLLL.s[pq].s[c].s[sym4[a][b]]
							+ partial_gLL_of_DgLLL.s[pq].s[b].s[sym4[a][c]]
							- partial_gLL_of_DgLLL.s[pq].s[a].s[bc]
						);
					}
				}
			}
		}
	}

	/*
	... if we use g_ab,cd and don't limit h->0 our ∂/∂g_pq(x) D_c(g_ab(x')) 's ...

	∂/∂g_pq R_abcd
	= ∂/∂g_pq (1/2 (g_ad,bc - g_bd,ac - g_ac,bd + g_bc,ad) + g^fg (Γ_fad Γ_gbc - Γ_fac Γ_gbd))
	= 1/2 ∂/∂g_pq (g_ad,bc - g_bd,ac - g_ac,bd + g_bc,ad)
		+ (∂/∂g_pq g^fg) (Γ_fad Γ_gbc - Γ_fac Γ_gbd)
		+ g^fg ∂/∂g_pq (Γ_fad Γ_gbc - Γ_fac Γ_gbd)
	= 1/2 ∂/∂g_pq (g_ad,bc - g_bd,ac - g_ac,bd + g_bc,ad)
		- g^pf g^qg (Γ_fad Γ_gbc - Γ_fac Γ_gbd)
		+ g^fg ∂/∂g_pq (Γ_fad Γ_gbc - Γ_fac Γ_gbd)
	= 1/2 ∂/∂g_pq (g_ad,bc - g_bd,ac - g_ac,bd + g_bc,ad)
		- Γ^p_ad Γ^q_bc + Γ^p_ac Γ^q_bd
		+ g^ef ∂/∂g_pq (Γ_ead Γ_fbc - Γ_eac Γ_fbd)
	*/
	real4s4x4x4x4x4 partial_gLL_of_RiemannLLLL;
	for (int p = 0; p < stDim; ++p) {
		for (int q = p; q < stDim; ++q) {
			int const pq = sym4[p][q];
			for (int a = 0; a < stDim; ++a) {
				for (int b = 0; b < stDim; ++b) {
					for (int c = 0; c < stDim; ++c) {
						for (int d = 0; d < stDim; ++d) {
							real sum = .5 * (
								  partial_gLL_of_D2gLLLL.s[pq].s[sym4[a][d]].s[sym4[b][c]]
								+ partial_gLL_of_D2gLLLL.s[pq].s[sym4[b][c]].s[sym4[a][d]]
								- partial_gLL_of_D2gLLLL.s[pq].s[sym4[b][d]].s[sym4[a][c]]
								- partial_gLL_of_D2gLLLL.s[pq].s[sym4[a][c]].s[sym4[b][d]]
							)
							- GammaULL.s[p].s[sym4[a][d]] * GammaULL.s[q].s[sym4[b][c]]
							+ GammaULL.s[p].s[sym4[a][c]] * GammaULL.s[q].s[sym4[b][d]];
							for (int e = 0; e < stDim; ++e) {
								for (int f = 0; f < stDim; ++f) {
									sum += gUU.s[sym4[e][f]] * (
										  partial_gLL_of_GammaLLL.s[pq].s[e].s[sym4[a][d]] * GammaLLL.s[f].s[sym4[b][c]]
										+ GammaLLL.s[e].s[sym4[a][d]] * partial_gLL_of_GammaLLL.s[pq].s[f].s[sym4[b][c]]
										- partial_gLL_of_GammaLLL.s[pq].s[e].s[sym4[a][c]] * GammaLLL.s[f].s[sym4[b][d]]
										+ GammaLLL.s[e].s[sym4[a][c]] * partial_gLL_of_GammaLLL.s[pq].s[f].s[sym4[b][d]]
									);
								}
							}
							partial_gLL_of_RiemannLLLL.s[pq].s[a].s[b].s[c].s[d] = sum;
						}
					}
				}
			}
		}
	}

	/*
	g^uv ∂/∂g_pq R_uavb
	= 1/2 g^uv ∂/∂g_pq (g_ub,av - g_ab,uv - g_uv,ab + g_av,ub)
		- g^uv Γ^p_ub Γ^q_av + Γ^p_uv Γ^q_ab
		+ g^uv g^fg ∂/∂g_pq (Γ_fub Γ_gav - Γ_fuv Γ_gab)
	= 1/2 g^uv ∂/∂g_pq (g_ub,av - g_ab,uv - g_uv,ab + g_av,ub)
		- Γ^pu_b Γ^q_ua + Γ^pu_u Γ^q_ab							<- looks like the 2nd term of the Ricci eqn based on g_ab,cd
		+ g^uv g^fg (
			  (∂/∂g_pq Γ_fub) Γ_gav
			+ Γ_fub (∂/∂g_pq Γ_gav)
			- (∂/∂g_pq Γ_fuv) Γ_gab
			- Γ_fuv (∂/∂g_pq Γ_gab)
		)

	partial_gLL_of_RicciLL.pq.ab := ∂R_ab/∂g_pq
	= ∂/∂g_pq (g^uv R_uavb)
	= (∂/∂g_pq g^uv) R_uavb + g^uv (∂/∂g_pq R_uavb)
	= -g^up g^vq R_uavb + g^uv (∂/∂g_pq R_uavb)
	= -R^p_a^q_b + g^uv (∂/∂g_pq R_uavb)

	*/
	real4s4x4s4 partial_gLL_of_RicciLL;
	for (int p = 0; p < stDim; ++p) {
		for (int q = p; q < stDim; ++q) {
			int const pq = sym4[p][q];
			for (int a = 0; a < stDim; ++a) {
				for (int b = a; b < stDim; ++b) {
					int const ab = sym4[a][b];
					real sum = 0;
					sum -= RiemannULUL.s[p].s[a].s[q].s[b];
					for (int u = 0; u < stDim; ++u) {
						for (int v = 0; v < stDim; ++v) {
							int const uv = sym4[u][v];
							sum += gUU.s[uv] * partial_gLL_of_RiemannLLLL.s[pq].s[u].s[a].s[v].s[b];
						}
					}
					partial_gLL_of_RicciLL.s[pq].s[ab] = sum;
				}
			}
		}
	}


<? else -- only use 1st-deriv finite-difference stencils (and difference-of-difference for 2nd-deriv)
?>

	//this is also in the Ricci computation, but should I be storing it?  is it too big?
	real4x4x4s4 const partial_xU_of_GammaULL = calc_partial_xU_of_GammaULL(GammaULLs);

	//RiemannULLL.a.b.cd := R^a_bcd = Γ^a_bd,c - Γ^a_bc,d + Γ^a_ec Γ^e_bd - Γ^a_ed Γ^e_bc
	real4x4x4x4 RiemannULLL;
	for (int a = 0; a < stDim; ++a) {
		for (int b = 0; b < stDim; ++b) {
			for (int c = 0; c < stDim; ++c) {
				for (int d = 0; d < stDim; ++d) {
					real sum =
						  partial_xU_of_GammaULL.s[c].s[a].s[sym4[b][d]]
						- partial_xU_of_GammaULL.s[d].s[a].s[sym4[b][c]];
					for (int e = 0; e < stDim; ++e) {
						sum +=
							  GammaULL.s[a].s[sym4[e][c]] * GammaULL.s[e].s[sym4[b][d]]
							- GammaULL.s[a].s[sym4[e][d]] * GammaULL.s[e].s[sym4[b][c]];
					}
					RiemannULLL.s[a].s[b].s[c].s[d] = sum;
				}
			}
		}
	}

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


	partial_gLL_of_RicciLL.pq.ab := ∂R_ab/∂g_pq
	= Γ^pc_c Γ^q_ba - Γ^pc_b Γ^q_ca - g^cp R^q_acb
	(work done in efe.html)
	*/
	real4s4x4s4 partial_gLL_of_RicciLL;
	for (int p = 0; p < stDim; ++p) {
		for (int q = p; q < stDim; ++q) {
			int const pq = sym4[p][q];
			for (int a = 0; a < stDim; ++a) {
				for (int b = a; b < stDim; ++b) {
					int const ab = sym4[a][b];
					real sum = 0;
					for (int c = 0; c < stDim; ++c) {
						sum +=
							  GammaUUL.s[p].s[c].s[c] * GammaULL.s[q].s[sym4[b][a]]
							- GammaUUL.s[p].s[c].s[b] * GammaULL.s[q].s[sym4[c][a]]
							- gUU.s[sym4[c][p]] * RiemannULLL.s[q].s[a].s[c].s[b];
					}
					partial_gLL_of_RicciLL.s[pq].s[ab] = sum;
				}
			}
		}
	}

<? end ?>

	//g^ab ∂R_ab/∂g_pq
	real4s4 gUU_times_partial_gLL_of_RicciLL;
	for (int p = 0; p < stDim; ++p) {
		for (int q = p; q < stDim; ++q) {
			int const pq = sym4[p][q];
			gUU_times_partial_gLL_of_RicciLL.s[pq] = real4s4_dot(gUU, partial_gLL_of_RicciLL.s[pq]);
		}
	}

	/*
	partial_gLL_of_EinsteinLL.pq.ab := ∂G_ab(x)/∂g_pq(x')
	= 1/2 δ(x - x') (g_ab(x) R^pq(x) - δ_a^p δ_b^q R(x)) + (δ_a^u δ_b^v - 1/2 g_ab(x) g^uv(x)) ∂R_uv(x)/∂g_pq(x')
	*/
	real4s4x4s4 partial_gLL_of_EinsteinLL;
	for (int pq = 0; pq < 10; ++pq) {
		for (int ab = 0; ab < 10; ++ab) {
			real sum = partial_gLL_of_RicciLL.s[pq].s[ab];
			if (pq == ab) {
				sum -= .5 * Gaussian;
			}
			sum -= .5 * gUU.s[ab] * (
				gUU_times_partial_gLL_of_RicciLL.s[pq]
				- RicciUU.s[pq]
			);
			partial_gLL_of_EinsteinLL.s[pq].s[ab] = sum;
		}
	}

	real4s4x4s4 const partial_gLL_of_8piTLL = calc_partial_gLL_of_8piTLL(TPrims[index], gLL, gUU);
	
	real4s4 const EFE = EFEs[index];	// G_ab - 8 π T_ab

	//lower times lower, but used for minimizing frobenius norm of EFE_ab
	real const EFE_LL_dot_gLL = real4s4_dot(EFE, gLL);

	// GammaSq_tr_23_32.ab := Γ^pv_u Γ^qu_v
	real4s4 GammaSq_tr_23_32;
	for (int p = 0; p < stDim; ++p) {
		for (int q = p; q < stDim; ++q) {
			int const pq = sym4[p][q];
			real sum = 0;
			for (int u = 0; u < stDim; ++u) {
				for (int v = 0; v < stDim; ++v) {
					sum += GammaUUL.s[p].s[v].s[u] * GammaUUL.s[q].s[u].s[v];
				}
			}
			GammaSq_tr_23_32.s[pq] = sum;
		}
	}

	//partial_gLL_of_Phi.pq := ∂Φ/∂g_pq
	real4s4 partial_gLL_of_Phi;
	// first calculate non-partial non ∂/∂g_pq terms:
	for (int p = 0; p < stDim; ++p) {
		for (int q = p; q < stDim; ++q) {
			int const pq = sym4[p][q];
			real sum = 
				// -(G_pq - 8 pi T_pq)) R
				-EFE.s[pq] * Gaussian				
				// + Sum_ab EFE_ab g_ab R^pq
				+ RicciUU.s[pq] * EFE_LL_dot_gLL;	
			;
			for (int a = 0; a < stDim; ++a) {
				for (int b = 0; b < stDim; ++b) {
					int const ab = sym4[a][b];
					// - Sum_ab EFE_ab dT_ab/dg_pq
					sum -= EFE.s[ab] * 8. * M_PI * partial_gLL_of_8piTLL.s[pq].s[ab];	
					// - Sum_ab EFE_ab R^p_a^q_b
					sum -= EFE.s[ab] * RiemannULUL.s[p].s[a].s[q].s[b];			
					// - EFE_ab Γ^pc_b Γ^q_ca
					for (int c = 0; c < stDim; ++c) {
						sum -= EFE.s[ab] * GammaUUL.s[p].s[c].s[b] * GammaULL.s[q].s[sym4[c][a]];
					}
					// + EFE_ab Γ^pu_u Γ^q_ab
					sum += EFE.s[ab] * Gamma23U.s[p] * GammaULL.s[q].s[ab]
					// + 1/2 EFE_ab g_ab Γ^pv_u Γ^qu_v
					sum += .5 * EFE.s[ab] * gLL.s[ab] * GammaSq_tr_23_32.s[pq];
					// - 1/2 EFE_ab g_ab Γ^pu_u Γ^qv_v
					sum -= .5 * EFE.s[ab] * gLL.s[ab] * Gamma23U.s[p] * Gamma23U.s[q];
				
					//TODO 1st derivs and stuff here
				}
			}
			partial_gLL_of_Phi.s[pq] = sum;
		}
	}
	// next calculate first-derivatives of ∂/∂g_pq terms:
	for (int offset = -<?=solver.diffOrder/2?>; offset <= <?=solver.diffOrder/2?>; ++offset) {
		if (offset == 0) continue;
		
	}

	global gPrim_t * const partial_gPrim_of_Phi = partial_gPrim_of_Phis + index;
	gPrim_t const gPrim = gPrims[index];

	real3s3 const gammaLL = gPrim.gammaLL;
	real3 const betaU = gPrim.betaU;
	real3 const betaL = real3s3_real3_mul(gammaLL, betaU);

<?
if solver.convergeAlpha then
?>	partial_gPrim_of_Phi->alpha = -2. * gPrim.alpha * partial_gLL_of_Phi.s00;
<?
end
if solver.convergeBeta then
?>
	for (int m = 0; m < sDim; ++m) {
		partial_gPrim_of_Phi->betaU.s[m] = 2. * partial_gLL_of_Phi.s00 * betaL.s[m];
		for (int n = 0; n < sDim; ++n) {
			partial_gPrim_of_Phi->betaU.s[m] += partial_gLL_of_Phi.s[sym4[0][n+1]] * gammaLL.s[sym3[n][m]];
		}
	}
<?
end
if solver.convergeGamma then
?>
	for (int m = 0; m < sDim; ++m) {
		for (int n = m; n < sDim; ++n) {
			int const mn = sym3[m][n];
			partial_gPrim_of_Phi->gammaLL.s[mn] =
				partial_gLL_of_Phi.s[sym4[m+1][n+1]]
				+ betaU.s[m] * (
					  partial_gLL_of_Phi.s00 * betaU.s[n]
					+ 2. * partial_gLL_of_Phi.s[sym4[0][n+1]]
				);
		}
	}
<?
end
?>

	//scale up our gradient?
	//scale by c^4 / G ~ 1e+44
	// which is the units of conversion
	//c^4/G * G_ab = 8 π T_ab
	partial_gPrim_of_Phi->alpha *= c*c*c*c/G;
	partial_gPrim_of_Phi->betaU = real3_real_mul(partial_gPrim_of_Phi->betaU, c*c*c*c/G);
	partial_gPrim_of_Phi->gammaLL = real3s3_real_mul(partial_gPrim_of_Phi->gammaLL, c*c*c*c/G);
}

kernel void update_gPrims(
	global gPrim_t * const gPrims,
	global gPrim_t const * const partial_gPrim_of_Phis,
	real const updateLambda
) {
	initKernel();
	global gPrim_t * const gPrim = gPrims + index;
	gPrim_t const partial_gPrim_of_Phi = partial_gPrim_of_Phis[index];

<?
if solver.convergeAlpha then
?>
	gPrim->alpha -= updateLambda * partial_gPrim_of_Phi.alpha;
<?
end
if solver.convergeBeta then
?>	
	gPrim->betaU = real3_mul_add(
		gPrim->betaU,
		partial_gPrim_of_Phi.betaU,
		-updateLambda
	);
<?
end
if solver.convergeGamma then
?>	
	gPrim->gammaLL = real3s3_mul_add(
		gPrim->gammaLL,
		partial_gPrim_of_Phi.gammaLL,
		-updateLambda
	);
<?
end
?>
}
