#!/usr/bin/env luajit
local table = require 'ext.table'
local range = require 'ext.range'
local tolua = require 'ext.tolua'
require 'symmath'.setup{
	implicitVars=true}
require 'symmath.export.MathJax'.setup{
	title = 'generating my gradient descent code',
}

g = var'g'
G = var'G'
R = var'R'
Gamma = var'\\Gamma'

local n = 1	-- nbhd size

local offsets = Array(2*n+1, 2*n+1, 2*n+1)

local xs = Array:lambda(offsets, function(i,j,k)
	return var('x_{'..table{i-n-1,j-n-1,k-n-1}:concat', '..'}')
end)
--printbr(xs)

local EFEs = Array:lambda(offsets, function(i,j,k)
	return var('EFE_{'..table{i-n-1, j-n-1, k-n-1}:concat', '..'}')
end)
--printbr(EFEs)

local Phi = var'\\Phi_{0,0,0}'
local sum = 0
for i=-n,n do
	for j=-n,n do
		for k=-n,n do
			sum = sum + EFEs[i+n+1][j+n+1][k+n+1]'_ab'^2
		end
	end
end
sum = (frac(1,2) * sum)()
local Phi_def = Phi:eq(sum)
printbr(Phi_def)

local function index(t, ...)
	if t == nil then return nil end
	if select('#', ...) == 0 then return t end
	local i = ...
	return index(t[i], select(2, ...))
end

local gs = Array:lambda(offsets, function(i,j,k)
	return var('g_{'..table{i-n-1,j-n-1,k-n-1}:concat', '..'}')
end)
Phi:setDependentVars(gs[n+1][n+1][n+1]'_ab')
local Gs = Array:lambda(offsets, function(i,j,k)
	local G = var('G_{'..table{i-n-1,j-n-1,k-n-1}:concat', '..'}')
	local is = {i-n-1, j-n-1, k-n-1}
	-- [[ 1st-derivative only = + pattern
	local stencil_gLLs = table()
	for offset=1,1 do
		for dir=1,3 do
			local dxi = {0,0,0}
			dxi[dir] = offset
			local iR = range(3):mapi(function(q) return is[q] + dxi[q] end)
			--local gR = gs[{iR[1]+n+1, iR[2]+n+1, iR[3]+n+1}]
			local gR = index(gs, iR[1]+n+1, iR[2]+n+1, iR[3]+n+1)
			if gR then 
				assert(Variable:isa(gR))
				stencil_gLLs:insert(gR'_ab')
			end
			local iL = range(3):mapi(function(q) return is[q] - dxi[q] end)
			-- what to do for oob defereferencing ...
			--local gL = gs[{iL[1]+n+1, iL[2]+n+1, iL[3]+n+1}]
			local gL = index(gs, iL[1]+n+1, iL[2]+n+1, iL[3]+n+1)
			if gL then
				assert(Variable:isa(gL))
				stencil_gLLs:insert(gL'_ab')
			end
		end
	end
	--]]
	-- [[ TODO 2nd derivative = Linf-norm
	--]]
	G'_ab':setDependentVars(stencil_gLLs:unpack())
	return G
end)
local Ts = Array:lambda(offsets, function(i,j,k)
	local T = var('T_{'..table{i-n-1,j-n-1,k-n-1}:concat', '..'}')
	local g = gs[i][j][k]
	T'_ab':setDependentVars(g'_ab')
	return T
end)
for i,G in Gs:iter() do
	Phi_def = Phi_def:replace(EFEs[i]'_ab', Gs[i]'_ab' - 8 * pi * Ts[i]'_ab')
end
printbr(Phi_def)

-- [[
local Rs = Array:lambda(offsets, function(i,j,k)
	local R = var('R_{'..table{i-n-1,j-n-1,k-n-1}:concat', '..'}')
	R'_ab':setDependentVars(gs[i][j][k]'_ab')
	--R'_ab':setDependentVars(all_gLLs:unpack())
	return R
end)
for i,G in Gs:iter() do
	local function makeGamma(a,b,c)
		return frac(1,2) * (
			finiteDifference(g(a)(b), i)(c)
			+ finiteDifference(g(a)(c), i)(b)
			- finiteDifference(g(b)(c), i)(a)
		)
	end
	--Phi_def = Phi_def:replace(Gs[i]'_ab', Rs[i]'_ab' - frac(1,2) * gs[i]'_ab' * Rs[i])
	Phi_def = Phi_def:replace(
		Gs[i]'_ab',
		(gs[i]'^uv' * delta'^c_a' * delta'^d_b' - frac(1,2) * gs[i]'_ab' * gs[i]'^uv' * gs[i]'^cd') * Rs[i]'_ucvd'
	)
	-- [=[
	:replace(
		Rs[i]'_ucvd',
		frac(1,2) * (
			  finiteDifference2(g'_ad', i)'_bc'
			+ finiteDifference2(g'_bc', i)'_ad'
			- finiteDifference2(g'_ac', i)'_bd'
			- finiteDifference2(g'_bd', i)'_ac'
		)
		+ gs[i]'^ef' * (
			makeGamma('_e', '_a', '_d') * makeGamma('_f', '_b', '_c')
			- makeGamma('_e', '_a', '_c') * makeGamma('_f', '_b', '_d')
		)
	)
	--]=]
end
printbr(Phi_def)
--]]

local grad_Phi_def = Phi_def:diff(gs[n+1][n+1][n+1]'_pq')()
printbr(grad_Phi_def)

local G_x = var'G(x)'
G_x'_ab':setDependentVars(g_xnot'_ab')
local T_x = var'T(x)'
T_x'_ab':setDependentVars(T_x'_ab')

local EFE_def = EFE_x'_ab':eq(G_x'_ab' - 8 * pi * T_x'_ab')
printbr(EFE_def)

local grad_EFE_def = EFE_def:diff(g_xnot'_pq')()
printbr(grad_EFE_def)

local g_x = var'g(x)'
local R_x = var'R(x)'
local Gamma_x = var'\\Gamma(x)'

do return end

local Einstein_def = G'_ab':eq(R'_ab' - frac(1,2) * g'_ab' * R'_pq' * g'^pq')
printbr(Einstein_def)

EFE_def = EFE_def:subst(Einstein_def)
printbr(EFE_def)

local Ricci_def = R'_ab':eq(R'^c_acb')
printbr(Ricci_def)

EFE_def = EFE_def:subst(Ricci_def, Ricci_def:reindex{abc='pqr'})
printbr(EFE_def)


do return end

EFE_def = EFE_def:subst(
	Riemann_def:reindex{abcd='cacb'},
	Riemann_def:reindex{abcd='rprq'}
)()
printbr(EFE_def)



local dgUUL_def = g'^uv_,c':eq(-g'^ue' * g'_ef,c' * g'^fv')
printbr(dgUUL_def)

local GammaLLL_def = Gamma'_abc':eq(frac(1,2) * (g'_ab,c' + g'_ac,b' - g'_bc,a'))
printbr(GammaLLL_def)

local dgLLL_wrt_GammaLLL_def = (GammaLLL_def + GammaLLL_def:reindex{abc='bac'}):symmetrizeIndexes(g, {1,2})():switch()
printbr(dgLLL_wrt_GammaLLL_def)

local Gamma_def = Gamma'^a_bc':eq(g'^au' * Gamma'_ubc')
printbr(Gamma_def)

local Riemann_def = R'^a_bcd':eq(Gamma'^a_bd,c' - Gamma'^a_bc,d' + Gamma'^a_uc' * Gamma'^u_bd' - Gamma'^a_ud' * Gamma'^u_bc')
printbr(Riemann_def)

Riemann_def = Riemann_def:subst(Gamma_def'_,d'(), Gamma_def'_,d'():reindex{cd='dc'})()
printbr(Riemann_def)

Riemann_def = Riemann_def:subst(
	GammaLLL_def'_,d'():reindex{abcd='ubcd'},
	GammaLLL_def'_,d'():reindex{abcd='ubdc'}
)()
	:symmetrizeIndexes(g, {1,2})
	:symmetrizeIndexes(g, {3,4})
	:simplify()
printbr(Riemann_def)

Riemann_def = Riemann_def:subst(
	dgUUL_def:reindex{uvc='auc'},
	dgUUL_def:reindex{uvc='aud'}
)()
printbr(Riemann_def)

Riemann_def = (Riemann_def * g'_va')
	:simplifyMetrics():reindex{v='a'}
printbr(Riemann_def)

Riemann_def = Riemann_def:subst(
	dgLLL_wrt_GammaLLL_def:reindex{abc='afc'},
	dgLLL_wrt_GammaLLL_def:reindex{abc='afd'}
):tidyIndexes()()
printbr(Riemann_def)

Riemann_def = (g'^va' * Riemann_def)()
	:simplifyMetrics()
	:reindex{va='ae'}
printbr(Riemann_def)

do return end

local TensorRef = require 'symmath.tensor.Ref'
local TensorIndex = require 'symmath.tensor.Index'
-- replace Gamma^a_bc with g^ad Gamma_dbc
EFE_def = EFE_def:map(function(expr)
	if TensorRef:isa(expr)
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
	if TensorRef:isa(expr)
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
	if TensorRef:isa(expr)
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
