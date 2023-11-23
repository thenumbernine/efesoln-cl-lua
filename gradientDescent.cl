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
				partial_gLL_of_8piTLL.s[ef].s[sym4[0][0]] += TPrim.E.s[e-1] * TPrim.E.s[f-1] + TPrim.B.s[e-1] * TPrim.B.s[f-1];
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
	//real4 const uU = (real4){.s={W, W * vU.s1, W * vU.s2, W * vU.s3}};
	real4 const uU = (real4){.s0 = W, .s1 = W * vU.s1, .s2 = W * vU.s2, .s3 = W * vU.s3}};
	real4 const uL = real4s4_real4_mul(gLL, uU);
<?
	else -- otherwise uL = gLL.s0
?>
	//real4 const uL = (real4){.s={gLL.s00, gLL.s01, gLL.s02, gLL.s03}};
	//real4 const uL = (real4){.s={gLL.s[sym4[0][0]], gLL.s[sym4[0][1]], gLL.s[sym4[0][2]], gLL.s[sym4[0][3]]}};
	real4 const uL = (real4){.s0 = gLL.s00, .s1 = gLL.s01, .s2 = gLL.s02, .s3 = gLL.s03};
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
	return real4s4x4s4_real_mul(partial_gLL_of_8piTLL, 8. * M_PI);
}

static constant int4 const int4_dirs[3] = {
	(int4)(1, 0, 0, 0),
	(int4)(0, 1, 0, 0),
	(int4)(0, 0, 1, 0),
};
static inline int4 int4_dir(int dim, int offset) {
	return int4_dirs[dim] * offset;
}

static inline real4s4 EFE_LL_minus_half_trace_at(
	int4 const i,
	global real4s4 const * const gLLs,
	global real4s4 const * const gUUs,
	global real4s4 const * const EFEs
) {
	if (i.x <= 0 || i.y <= 0 || i.z <= 0 ||
		i.x >= size.x || i.y >= size.y || i.z >= size.z
	) {
		return real4s4_zero;	//TODO ... consider boundary conditions
	}
	
	int const index = indexForInt4ForSize(i, size.x, size.y, size.z);
	real4s4 const gLL = gLLs[index];	// g_ab
	real4s4 const gUU = gUUs[index];	// g^ab
	real4s4 const EFE = EFEs[index];	// G_ab - 8 π T_ab
	
	//lower times lower, but used for minimizing frobenius norm of EFE_ab
	// Sum_ab (G_ab - 8 π T_ab) g_ab
	real const EFE_LL_dot_gLL = real4s4_dot(EFE, gLL);
	
	// common term in the gradient descent: 
	// (G_uv - 8 π T_uv) - 1/2 (G_ab - 8 π T_ab) g_ab g^uv
	return real4s4_mul_add(EFE, gUU, -.5 * EFE_LL_dot_gLL);
}

// this is hardcoded to g_ab = η_ab at boundaries
//TODO ... consider boundary conditions
// or TODO ... put all g^ab boundary conditions here, and use this function in the finite-difference calculations
static inline real4s4 gUU_at(
	int4 const i,
	global real4s4 const * const gUUs
) {
	if (i.x <= 0 || i.y <= 0 || i.z <= 0 ||
		i.x >= size.x || i.y >= size.y || i.z >= size.z
	) {
		return real4s4_Minkowski;
	}
	int const index = indexForInt4ForSize(i, size.x, size.y, size.z);
	return gUUs[index];
}

//GammaULL.a.b.c := Γ^a_bc
static inline real4x4s4 GammaULL_at(
	int4 const i,
	global real4x4s4 const * const GammaULLs
) {
	if (i.x <= 0 || i.y <= 0 || i.z <= 0 ||
		i.x >= size.x || i.y >= size.y || i.z >= size.z
	) {
		return real4x4s4_zero;
	}
	
	int const index = indexForInt4ForSize(i, size.x, size.y, size.z);
	return GammaULLs[index];
}

//GammaUUL.a.b.c := Γ^ab_c = Γ^a_dc g^db
static inline real4x4x4 GammaUUL_at(
	int4 const i,
	global real4s4 const * const gUUs,
	global real4x4s4 const * const GammaULLs
) {
	if (i.x <= 0 || i.y <= 0 || i.z <= 0 ||
		i.x >= size.x || i.y >= size.y || i.z >= size.z
	) {
		return real4x4x4_zero;
	}
	
	int const index = indexForInt4ForSize(i, size.x, size.y, size.z);
	return real4x4s4_real4s4_mul21(GammaULLs[index], gUUs[index]);
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
	
	//RicciLL.ab := R_ab = R^c_acb
	real4s4 const RicciLL = real4x4x4x4_tr13_to_real4s4(RiemannULLL);

	//RicciUL.a.b := R^a_b = g^ac R_cb
	real4x4 const RicciUL = real4s4_real4s4_mul(gUU, RicciLL);

	//RicciUU.ab := R^ab = R^a_c g^cb
	real4s4 const RicciUU = real4x4_real4s4_to_real4s4_mul(RicciUL, gUU);

	//Gaussian := R = R^a_a
	real const Gaussian = real4x4_tr(RicciUL);

	real4s4x4s4 const partial_gLL_of_8piTLL = calc_partial_gLL_of_8piTLL(TPrims[index], gLL, gUU);
	
	real4s4 const EFE = EFEs[index];	// G_ab - 8 π T_ab

	//lower times lower, but used for minimizing frobenius norm of EFE_ab
	// Sum_ab (G_ab - 8 π T_ab) g_ab
	real const EFE_LL_dot_gLL = real4s4_dot(EFE, gLL);
	
	// common term in the gradient descent: 
	// (G_uv - 8 π T_uv) - 1/2 (G_ab - 8 π T_ab) g_ab g^uv
	real4s4 EFE_LL_minus_half_trace = real4s4_mul_add(EFE, gUU, -.5 * EFE_LL_dot_gLL);
	
	//GammaUUL.a.b.c := Γ^ab_c = Γ^a_dc g^db
	real4x4x4 const GammaUUL = real4x4s4_real4s4_mul21(GammaULL, gUU);
	
	//Gamma23U.a := Γ^a = Γ^au_u
	real4 const Gamma23U = real4x4x4_tr23(GammaUUL);

	// GammaSq_tr_2_2.pq.uv := Γ^pc_v Γ^q_cu - Γ^pc_c Γ^q_uv
	// symmetries?
	real4x4x4x4 GammaSq_tr_2_2;
	for (int p = 0; p < stDim; ++p) {
		for (int q = 0; q < stDim; ++q) {
			for (int u = 0; u < stDim; ++u) {
				for (int v = 0; v < stDim; ++v) {
					real sum = - Gamma23U.s[p] * GammaULL.s[q].s[sym4[u][v]];
					for (int c = 0; c < stDim; ++c) {
						sum += GammaUUL.s[p].s[c].s[v] * GammaULL.s[q].s[sym4[c][u]];
					}
					GammaSq_tr_2_2.s[p].s[q].s[u].s[v] = sum;
				}
			}
		}
	}

	//partial_gLL_of_Phi.pq := ∂Φ/∂g_pq
	real4s4 partial_gLL_of_Phi;
	// first calculate non-partial non ∂/∂g_pq terms:
	for (int p = 0; p < stDim; ++p) {
		for (int q = p; q < stDim; ++q) {
			int const pq = sym4[p][q];
			real sum = 0;
#if 1
			// -(G_pq - 8 π T_pq) R
			sum -= EFE.s[pq] * Gaussian;
#endif 
#if 0
			// + Sum_ab EFE_ab 1/2 g_ab R^pq
			sum += EFE_LL_dot_gLL * .5 * RicciUU.s[pq];
#endif 
#if 0
			for (int a = 0; a < stDim; ++a) {
				for (int b = 0; b < stDim; ++b) {
					int const ab = sym4[a][b];
					// - Sum_ab EFE_ab 8 π dT_ab/dg_pq
					sum -= EFE.s[ab] * partial_gLL_of_8piTLL.s[pq].s[ab];
				}
			}		
#endif 
#if 0
			for (int u = 0; u < stDim; ++u) {
				for (int v = 0; v < stDim; ++v) {
					int const uv = sym4[u][v];
					sum -= EFE_LL_minus_half_trace.s[uv] * (
					// - (EFE_uv - 1/2 EFE_ab g_ab g^uv) R^p_u^q_v
						RiemannULUL.s[p].s[u].s[q].s[v]
					// - (EFE_uv - 1/2 EFE_ab g_ab g^uv) (Γ^pc_v Γ^q_cu - Γ^pc_c Γ^q_uv)
						+ GammaSq_tr_2_2.s[p].s[q].s[u].s[v]
					);
				}
			}
#endif
#if 0
			// next calculate first-derivatives of ∂/∂g_pq terms:
			for (int dim = 0; dim < sDim; ++dim) {
				for (int offset = -<?=solver.diffOrder/2?>; offset <= <?=solver.diffOrder/2?>; ++offset) {
					if (offset == 0) continue;
					int4 const iofs = i + int4_dir(dim, offset);
					// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) g^cd (δ^f_c δ^g_v Γ^e_ud + δ^f_u δ^g_d Γ^e_cv - δ^f_c δ^g_d Γ^e_uv - δ6f_u δ^g_v Γ^e_cd) 1/2 (∂/∂g_pq(x') (D_g[g_ef] + D_f[g_eg] - D_e[g_fg]) = δ^p_e δ^q_f d1coeff_g + δ^p_e δ^q_g d1coeff_f - δ^p_f δ^q_g d1coeff_e)
					// ... expanding all 12 terms with symmath and simplifying...
					/*
					asking symmath to simplify this for me:
					let H^uv = (G_ab - 8 π T_ab) (δ^u_a δ^v_b - 1/2 g_ab g^uv) = EFE_uv - 1/2 EFE_ab g_ab g^uv
						> delta = Tensor.deltaSymbol()
						> g = Tensor.metricSymbol()
						> (H'^uv' * g'^cd' * (delta'^f_c' * delta'^g_v' * Gamma'^e_ud' + delta'^f_u' * delta'^g_d' * Gamma'^e_cv' - delta'^f_c' * delta'^g_d' * Gamma'^e_uv' - delta'^f_u' * delta'^g_v' * Gamma'^e_cd') * (delta'^p_e' * delta'^q_f' * phi'_g' + delta'^p_e' * delta'^q_g' * phi'_f' - delta'^p_f' * delta'^q_g' * phi'_e'))
							:simplifyMetrics()
							:symmetrizeIndexes(H, {1,2})()
							:symmetrizeIndexes(Gamma, {2,3})()
							:symmetrizeIndexes(g, {1,2})()
							:tidyIndexes()()
							:symmetrizeIndexes(H, {1,2})()
							:favorTensorVariance(phi'_a')()
					...gives...
							  φ_a * H^pq * Γ^ab_b 
						-     φ_a * H^pb * Γ^aq_b 
						-     φ_a * H^bq * Γ^ap_b 
						+     φ_a * H^bc * Γ^a_bc * g^pq
						- 2 * φ_a * H^aq * Γ^pb_b 
						+ 2 * φ_a * H^ab * Γ^pq_b 
						+ 2 * φ_a * H^bq * Γ^pa_b 
						- 2 * φ_a * H^bc * Γ^p_cb * g^aq
					*/				
					int const a = dim+1;
					for (int b = 0; b < stDim; ++b) {
						sum += .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[pq] * GammaUUL_at(iofs, gUUs, GammaULLs).s[a].s[b].s[b] * d1coeff_for_offset(offset) * inv_dx.s[dim];
						sum -= .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[p][b]] * GammaUUL_at(iofs, gUUs, GammaULLs).s[a].s[q].s[b] * d1coeff_for_offset(offset) * inv_dx.s[dim];
						sum -= .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[b][q]] * GammaUUL_at(iofs, gUUs, GammaULLs).s[a].s[p].s[b] * d1coeff_for_offset(offset) * inv_dx.s[dim];
						for (int c = 0; c < stDim; ++c) {
							sum += .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[b][c]] * GammaULL_at(iofs, GammaULLs).s[a].s[sym4[b][c]] * gUU_at(iofs, gUUs).s[pq] * d1coeff_for_offset(offset) * inv_dx.s[dim];
						}
						sum -= EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[a][q]] * GammaUUL_at(iofs, gUUs, GammaULLs).s[p].s[b].s[b] * d1coeff_for_offset(offset) * inv_dx.s[dim];
						sum += EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[a][b]] * GammaUUL_at(iofs, gUUs, GammaULLs).s[p].s[q].s[b] * d1coeff_for_offset(offset) * inv_dx.s[dim];
						sum += EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[b][q]] * GammaUUL_at(iofs, gUUs, GammaULLs).s[p].s[a].s[b] * d1coeff_for_offset(offset) * inv_dx.s[dim];
						for (int c = 0; c < stDim; ++c) {
							sum -= EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[b][c]] * GammaULL_at(iofs, GammaULLs).s[p].s[sym4[c][b]] * gUU_at(iofs, gUUs).s[sym4[a][q]] * d1coeff_for_offset(offset) * inv_dx.s[dim];
						}
					}
				}
			}
#endif
#if 0
			// next calculate second-derivatives of ∂/∂g_pq terms:
			for (int dim1 = 0; dim1 < sDim; ++dim1) {
				for (int dim2 = 0; dim2 < sDim; ++dim2) {
					if (dim1 == dim2) {
						for (int offset = -<?=solver.diffOrder/2?>; offset <= <?=solver.diffOrder/2?>; ++offset) {
							int4 const iofs = i + int4_dir(dim1, offset);
							{
								// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^cd (∂/∂g_pq(x') D^2_ud[g_cv] = δ^p_c δ^q_v D^2 d2coeff_ud)
								// = + (EFE_uq - 1/2 EFE_ab g_ab g^uq) 1/2 g^pd (d2coeff_ud) |x'=x + dx^u=dx^d
								int const u = dim1+1;
								int const d = u;
								sum += .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[u][q]] 
									* gUU_at(iofs, gUUs).s[sym4[p][d]] 
									* (d2coeffs[abs(offset)] * inv_dx.s[dim1]);
							}{
								// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^cd (∂/∂g_pq(x') D^2_cv[g_ud] = δ^p_u δ^q_d D^2 d2coeff_cv)
								// = + (EFE_pv - 1/2 EFE_ab g_ab g^pv) 1/2 g^cq (d2coeff_cv) |x'=x + dx^v=dx^c
								int const v = dim1+1;
								int const c = v;
								sum += .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[p][v]] 
									* gUU_at(iofs, gUUs).s[sym4[c][q]] 
									* (d2coeffs[abs(offset)] * inv_dx.s[dim1]);
							}{
								// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^cd (∂/∂g_pq(x') D^2_cd[g_uv] = δ^p_u δ^q_v D^2 d2coeff_cd)
								// = + (EFE_pq - 1/2 EFE_ab g_ab g^pq) 1/2 g^cd (d2coeff_cd) |x'=x + dx^c=dx^d
								int const c = dim1+1;
								int const d = c;
								sum -= .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[pq] 
									* gUU_at(iofs, gUUs).s[sym4[c][d]] 
									* (d2coeffs[abs(offset)] * inv_dx.s[dim1]);
							}{
								// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^cd (∂/∂g_pq(x') D^2_uv[g_cd] = δ^p_c δ^q_d D^2 d2coeff_uv)
								// = + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^pq (d2coeff_uv) |x'=x + dx^u=dx^v
								int const u = dim1+1;
								int const v = u;
								sum -= .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[u][v]]
									* gUU_at(iofs, gUUs).s[pq]
									* (d2coeffs[abs(offset)] * inv_dx.s[dim1]);
							}
						}
					} else {
						for (int offset1 = -<?=solver.diffOrder/2?>; offset1 <= <?=solver.diffOrder/2?>; ++offset1) {
							if (offset1 == 0) continue;
							for (int offset2 = -<?=solver.diffOrder/2?>; offset2 <= <?=solver.diffOrder/2?>; ++offset2) {
								if (offset2 == 0) continue;
								int4 const iofs = i + int4_dir(dim1, offset1) + int4_dir(dim2, offset2);
								{
									// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^cd (∂/∂g_pq(x') D^2_ud[g_cv] = δ^p_c δ^q_v D^2 d1coeff_u d1coeff_d)
									// = + (EFE_uq - 1/2 EFE_ab g_ab g^uq) 1/2 g^pd (d1coeff_u d1coeff_d) |x'=x + dx^u + dx^d
									int const u = dim1+1;
									int const d = dim2+1;
									int const offset_u = offset1;
									int const offset_d = offset2;
									sum += .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[u][q]]
										* gUU_at(iofs, gUUs).s[sym4[p][d]] 
										* (d1coeff_for_offset(offset_u) * d1coeff_for_offset(offset_d) * inv_dx.s[dim1] * inv_dx.s[dim2]);
								}{
									// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^cd (∂/∂g_pq(x') D^2_cv[g_ud] = δ^p_u δ^q_d D^2 d1coeff_c d1coeff_v)
									// = + (EFE_pv - 1/2 EFE_ab g_ab g^pv) 1/2 g^cq (d1coeff_c d1coeff_v) |x'=x + dx^v + dx^c
									int const v = dim1+1;
									int const c = dim2+1;
									int const offset_v = offset1;
									int const offset_c = offset2;
									sum += .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[p][v]] 
										* gUU_at(iofs, gUUs).s[sym4[c][q]] 
										* (d1coeff_for_offset(offset_v) * d1coeff_for_offset(offset_c) * inv_dx.s[dim1] * inv_dx.s[dim2]);
								}{
									// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^cd (∂/∂g_pq(x') D^2_cd[g_uv] = δ^p_u δ^q_v D^2 d1coeff_c d1coeff_d)
									// = + (EFE_pq - 1/2 EFE_ab g_ab g^pq) 1/2 g^cd (d1coeff_c d1coeff_d) |x'=x + dx^c + dx^d
									int const c = dim1+1;
									int const d = dim2+1;
									int const offset_c = offset1;
									int const offset_d = offset2;
									sum -= .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[pq] 
										* gUU_at(iofs, gUUs).s[sym4[c][d]] 
										* (d1coeff_for_offset(offset_c) * d1coeff_for_offset(offset_d) * inv_dx.s[dim1] * inv_dx.s[dim2]);
								}{
									// + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^cd (∂/∂g_pq(x') D^2_uv[g_cd] = δ^p_c δ^q_d D^2 d1coeff_u d1coeff_v)
									// = + (EFE_uv - 1/2 EFE_ab g_ab g^uv) 1/2 g^pq (d1coeff_u d1coeff_v) |x'=x + dx^u + dx^v
									int const u = dim1+1;
									int const v = dim2+1;
									int const offset_u = offset1;
									int const offset_v = offset2;
									sum -= .5 * EFE_LL_minus_half_trace_at(iofs, gLLs, gUUs, EFEs).s[sym4[u][v]]
										* gUU_at(iofs, gUUs).s[pq]
										* (d1coeff_for_offset(offset_u) * d1coeff_for_offset(offset_v) * inv_dx.s[dim1] * inv_dx.s[dim2]);
								}
							}
						}
					}
				}
			}
#endif			
			partial_gLL_of_Phi.s[pq] = sum;
		}
	}

	global gPrim_t * const partial_gPrim_of_Phi = partial_gPrim_of_Phis + index;
	partial_gPrim_of_Phi->alpha = 0;
	partial_gPrim_of_Phi->betaU = real3_zero;
	partial_gPrim_of_Phi->gammaLL = real3s3_zero;

	gPrim_t const gPrim = gPrims[index];
	real3s3 const gammaLL = gPrim.gammaLL;
	real3 const betaU = gPrim.betaU;
	real3 const betaL = real3s3_real3_mul(gammaLL, betaU);

<?
if solver.convergeAlpha then
?>	partial_gPrim_of_Phi->alpha = -2. * gPrim.alpha * partial_gLL_of_Phi.s[sym4[0][0]];
<?
end
if solver.convergeBeta then
?>
	for (int m = 0; m < sDim; ++m) {
		partial_gPrim_of_Phi->betaU.s[m] = 2. * partial_gLL_of_Phi.s[sym4[0][0]] * betaL.s[m];
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
					  partial_gLL_of_Phi.s[sym4[0][0]] * betaU.s[n]
					+ 2. * partial_gLL_of_Phi.s[sym4[0][n+1]]
				);
		}
	}
<?
end
?>

#if 0
	//scale up our gradient?
	//scale by c^4 / G ~ 1e+44
	// which is the units of conversion
	//c^4/G * G_ab = 8 π T_ab
	partial_gPrim_of_Phi->alpha *= c*c*c*c/G;
	partial_gPrim_of_Phi->betaU = real3_real_mul(partial_gPrim_of_Phi->betaU, c*c*c*c/G);
	partial_gPrim_of_Phi->gammaLL = real3s3_real_mul(partial_gPrim_of_Phi->gammaLL, c*c*c*c/G);
#endif
#if 0 //debugging
	partial_gPrim_of_Phi->alpha = 0;
	partial_gPrim_of_Phi->betaU = real3_zero;
	partial_gPrim_of_Phi->gammaLL = real3s3_zero;
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
