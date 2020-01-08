#!/usr/bin/env luajit
require 'ext'
require 'symmath'.setup{implicitVars=true}
require 'symmath.tostring.MathJax'.setup{
	title = 'generating my gradient descent code',
	url = 'file:///home/chris/Projects/christopheremoore.net/MathJax/MathJax.js?config=TeX-MML-AM_CHTML',
}

local EFE_def = G'_ab':eq(8 * pi * T'_ab')
printbr(EFE_def)

local Einstein_def = G'_ab':eq(R'_ab' - frac(1,2) * g'_ab' * R'_pq' * g'^pq')
printbr(Einstein_def)

EFE_def = EFE_def:subst(Einstein_def)
printbr(EFE_def)

local Ricci_def = R'_ab':eq(R'^c_acb')
printbr(Ricci_def)

EFE_def = EFE_def:subst(Ricci_def, Ricci_def:reindex{abc='pqr'})
printbr(EFE_def)

local Riemann_def = R'^a_bcd':eq(Gamma'^a_bd,c' - Gamma'^a_bc,d' + Gamma'^a_uc' * Gamma'^u_bd' - Gamma'^a_ud' * Gamma'^u_bc')
printbr(Riemann_def)

Riemann_def = Riemann_def
	:replace(Gamma'^a_uc' * Gamma'^u_bd',
		Gamma'^a_tc' * Gamma'^t_bd' + Gamma'^a_kc' * Gamma'^k_bd')
	:replace(Gamma'^a_ud' * Gamma'^u_bc',
		Gamma'^a_td' * Gamma'^t_bc' + Gamma'^a_kd' * Gamma'^k_bc')
	:simplify()
printbr(Riemann_def)

EFE_def = EFE_def:subst(
	Riemann_def:reindex{abcd='cacb'},
	Riemann_def:reindex{abcd='rprq'}
)()
printbr(EFE_def)

-- Binary{*}[
--		Constant[2], 
--		TensorRef{
--			Variable[\Gamma], 
--			TensorIndex{^c}, 
--			TensorIndex{_k}, 
--			TensorIndex{_b}
--		},
--		TensorRef{
-- 			Variable[\Gamma], 
--			TensorIndex{^k}, 
--			TensorIndex{_a}, 
--			TensorIndex{_c}
--		}
--	]

-- TODO you can't just substitute t=c and k=c, you must only do this on terms that include a 'c'
EFE_def = EFE_def
	:replace(Gamma'^c_ab,c', Gamma'^t_ab,t' + Gamma'^k_ab,k')()
	:replace(Gamma'^c_ac,b', Gamma'^t_at,b' + Gamma'^k_ak,b')()
	:replace( (2 * Gamma'^c_kb' * Gamma'^k_ac')(),
		(2 * (Gamma'^t_kb' * Gamma'^k_at' + Gamma'^l_kb' * Gamma'^k_al'))())()
	:replace( Gamma'^c_kc', Gamma'^t_kt' + Gamma'^l_kl' )()
	:replace( (2 * Gamma'^c_tb' * Gamma'^t_ac')(),
		(2 * (Gamma'^t_tb' * Gamma'^t_at' + Gamma'^k_tb' * Gamma'^t_ak'))())()
	:replace( Gamma'^c_tc', Gamma'^t_tt' + Gamma'^k_tk' )()
	:replace( (g'_ab' * Gamma'^r_kq' * Gamma'^k_pr' * g'^pq')(),
		(g'_ab' * ((Gamma'^t_kt' * Gamma'^k_tt' + Gamma'^l_kt' * Gamma'^k_tl') * g'^tt'
				+ (Gamma'^t_kt' * Gamma'^k_mt' + Gamma'^l_kt' * Gamma'^k_ml') * g'^mt'
				+ (Gamma'^t_km' * Gamma'^k_tt' + Gamma'^l_km' * Gamma'^k_tl') * g'^tm'
				+ (Gamma'^t_kn' * Gamma'^k_mt' + Gamma'^l_kn' * Gamma'^k_ml') * g'^mn'
				))() )()
	:replace(Gamma'^r_kr', Gamma'^t_kt' + Gamma'^l_kl')()
	:replace( (g'_tt' * Gamma'^l_kl' * Gamma'^k_pq' * g'^pq')(),
		(g'_tt' * Gamma'^l_kl' * (Gamma'^k_tt' * g'^tt'
								+ Gamma'^k_mt' * g'^mt'
								+ Gamma'^k_tm' * g'^tm'
								+ Gamma'^k_mn' * g'^mn'))())()
	:replace( (g'_ab' * Gamma'^l_kl' * Gamma'^k_pq' * g'^pq')(),
		(g'_ab' * Gamma'^l_kl' * (Gamma'^k_tt' * g'^tt'
								+ Gamma'^k_mt' * g'^mt'
								+ Gamma'^k_tm' * g'^tm'
								+ Gamma'^k_mn' * g'^mn'))())()
	:replace( (g'_ab' * Gamma'^r_pq,r' * g'^pq')(),
		(g'_ab' * ((Gamma'^t_tt,t' + Gamma'^k_tt,k') * g'^tt'
				+ (Gamma'^t_tl,t' + Gamma'^k_tl,k') * g'^tl'
				+ (Gamma'^t_lt,t' + Gamma'^k_lt,k') * g'^lt'
				+ (Gamma'^t_lm,t' + Gamma'^k_lm,k') * g'^lm'))())()
	:replace( (g'_ab' * Gamma'^r_pr,q' * g'^pq')(),
		(g'_ab' * ((Gamma'^t_tt,t' + Gamma'^k_tk,t') * g'^tt'
				+ (Gamma'^t_lt,t' + Gamma'^k_lk,t') * g'^lt'
				+ (Gamma'^t_tt,l' + Gamma'^k_tk,l') * g'^tl'
				+ (Gamma'^t_lt,m' + Gamma'^k_lk,m') * g'^lm'))())()
	:replace( (g'_ab' * Gamma'^r_tq' * Gamma'^t_pr' * g'^pq')(),
		(g'_ab' * ((Gamma'^t_tt' * Gamma'^t_tt' + Gamma'^k_tt' * Gamma'^t_tk') * g'^tt'
				+ (Gamma'^t_tl' * Gamma'^t_tt' + Gamma'^k_tl' * Gamma'^t_tk') * g'^tl'
				+ (Gamma'^t_tt' * Gamma'^t_lt' + Gamma'^k_tt' * Gamma'^t_lk') * g'^lt'
				+ (Gamma'^t_tm' * Gamma'^t_lt' + Gamma'^k_tm' * Gamma'^t_lk') * g'^lm'))())()
	:replace( (g'_ab' * Gamma'^r_tr' * Gamma'^t_pq' * g'^pq')(),
		(g'_ab' * (Gamma'^t_tt' + Gamma'^k_tk') * 
			(Gamma'^t_tt' * g'^tt' 
			+ Gamma'^t_tl' * g'^tl' 
			+ Gamma'^t_lt' * g'^lt' 
			+ Gamma'^t_lm' * g'^lm'))())()
	:replace( (g'_ab' * Gamma'^t_kt' * Gamma'^k_pq' * g'^pq')(),
		(g'_ab' * Gamma'^t_kt' * (Gamma'^k_tt' * g'^tt'
								+ Gamma'^k_tl' * g'^tl'
								+ Gamma'^k_lt' * g'^lt'
								+ Gamma'^k_lm' * g'^lm'))())()
printbr(EFE_def)

local TensorRef = require 'symmath.tensor.TensorRef'
local TensorIndex = require 'symmath.tensor.TensorIndex'
-- replace Gamma^a_bc with g^ad Gamma_dbc
EFE_def = EFE_def:map(function(expr)
	if TensorRef.is(expr)
	and expr[1] == Gamma
	then
		local a,b,c = table.unpack(expr, 2)	-- upper, lower, lower 
		assert(not a.lower)
		assert(b.lower)
		assert(c.lower)
		if #expr >= 4 then
			-- the extra indexes should be derivatives
			local result = TensorRef(g, a, TensorIndex{symbol='t', lower=false})
				* TensorRef(Gamma, TensorIndex{symbol='t', lower=true}, b, c)
				+ TensorRef(g, a, TensorIndex{symbol='o', lower=false})
				* TensorRef(Gamma, TensorIndex{symbol='o', lower=true}, b, c)
			if #expr > 4 then
				result = TensorRef(result, table.unpack(expr, 5))
			end
			return result
		end
	end
end)()
-- replace g^ab_,c with g^ae g_ef,c g^fb
EFE_def = EFE_def:map(function(expr)
	if TensorRef.is(expr)
	and expr[1] == g
	and #expr == 4
	then
		local a,b,c = table.unpack(expr, 2)
		if not a.lower
		and not b.lower
		and c.lower
		and c.derivative == 'partial'
		then
			return -TensorRef(g, a, TensorIndex{symbol='p', lower=false})
					* TensorRef(g, TensorIndex{symbol='p', lower=true}, TensorIndex{symbol='q', lower=true}, c)
					* TensorRef(g, TensorIndex{symbol='q', lower=false}, b)
		end
	end
end)()
printbr(EFE_def)
-- replace Gamma_abc with 1/2 (g_ab,c + g_ac,b - g_bc,a)
-- ... this is getting not enough memory ...
EFE_def = EFE_def:map(function(expr)
	if TensorRef.is(expr)
	and expr[1] == Gamma
	then
		local a,b,c = table.unpack(expr, 2)
		assert(a.lower and not a.derivative)
		assert(b.lower and not b.derivative)
		assert(c.lower and not c.derivative)
		local result = frac(1,2) * (
			TensorRef(g, a, b, TensorIndex{symbol=c.symbol, lower=true, derivative='partial'})
			+ TensorRef(g, a, c, TensorIndex{symbol=b.symbol, lower=true, derivative='partial'})
			- TensorRef(g, b, c, TensorIndex{symbol=a.symbol, lower=true, derivative='partial'})
		)
		if #expr > 4 then
			result = TensorRef(result, table.unpack(expr, 5))
		end
		return result
	end
end)()

printbr'(note now opq are spatial)'
printbr'split into time+space:'

local EFE_tt_def = EFE_def:reindex{ab='tt'}()
printbr(EFE_tt_def)

local EFE_ti_def = EFE_def:reindex{ab='ti'}()
printbr(EFE_ti_def)

local EFE_ij_def = EFE_def:reindex{ab='ij'}()
printbr(EFE_ij_def)
