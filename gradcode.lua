#!/usr/bin/env luajit
local table = require 'ext.table'
local range = require 'ext.range'
local tolua = require 'ext.tolua'
require 'symmath'.setup()
require 'symmath.export.MathJax'.setup{
	title = 'generating my gradient descent code',
}

local t,x,y,z = vars('t', 'x', 'y', 'z')
local C = Tensor.Chart{coords={t,x,y,z}}
local CS = Tensor.Chart{coords={x,y,z}, symbols='ijklmn'}
local Ct = Tensor.Chart{coords={t}, symbols='t'}
local Cx = Tensor.Chart{coords={x}, symbols='x'}
local Cy = Tensor.Chart{coords={y}, symbols='y'}
local Cz = Tensor.Chart{coords={z}, symbols='z'}

local g = var'g'
local G = var'G'
local R = var'R'
local Gamma = var'\\Gamma'

local n = 1	-- nbhd size

local offsets = Array(2*n+1, 2*n+1, 2*n+1)

local function index(t, ...)
	if t == nil then return nil end
	if select('#', ...) == 0 then return t end
	local i = ...
	return index(t[i], select(2, ...))
end

--[[ don't think I'm using this at the moment ...
local xs = Array:lambda(offsets, function(i,j,k)
	--[=[
	return Tensor('^i', function(a)
		return var('x^'..a..'['..table{i-n-1,j-n-1,k-n-1}:concat','..']')
	end)
	--]=]
	-- [=[ notice :makeDense uses the correct # of dims for index 'i' (being 3) but uses the letters from txyz instead of xyz ...
	return var('x['..table{i-n-1,j-n-1,k-n-1}:concat','..']')'^i':makeDense()
	--]=]
end)
printbr(xs)
do return end
--]]

local function removeSymmetricDuplicates(T)
	-- remove symmetric duplicate variables
	for a=1,3 do
		for b=a+1,4 do
			T[b][a] = T[a][b]
		end
	end
end

local Phi = var'\\Phi'
local all_gLLs = table()
local gs = Array:lambda(offsets, function(i,j,k)
	local g_ijk = var('g['..table{i-n-1,j-n-1,k-n-1}:concat', '..']')'_ab':makeDense()
	removeSymmetricDuplicates(g_ijk)
	for ab,g_ijk_ab in g_ijk:iter() do
		g_ijk_ab.isMetric = true	-- hack for differentiation
		all_gLLs:insert(g_ijk_ab)
	end
	return g_ijk
end)

-- I'm going to make a separate copy of gs in raised-form
-- for the sake of using dense-tensors in raised and lowered
-- however this might get ugly when differentiating wrt the lower specific-index metric ...
local all_gUU = table()
local gUs = Array:lambda(offsets, function(i,j,k)
	local gU_ijk = var('g['..table{i-n-1,j-n-1,k-n-1}:concat', '..']')'^ab':makeDense()
	removeSymmetricDuplicates(gU_ijk)
	for ab,g_ijk_ab in gU_ijk:iter() do
		g_ijk_ab.isMetric = true	-- hack for differentiation
		all_gUU:insert(g_ijk_ab)
	end
	return gU_ijk
end)

do
	local deps = table()
	for a=1,4 do
		for b=1,4 do
			deps:insert(gs[n+1][n+1][n+1][a][b])
		end
	end
	-- ... is phi only dependent on the center g? nah ...
	--Phi:setDependentVars(deps:unpack())
	Phi:setDependentVars(all_gLLs:unpack())
end
local Phis = Array:lambda(offsets, function(i,j,k)
	local Phi_ijk = var('\\Phi['..table{i-n-1,j-n-1,k-n-1}:concat', '..']')
	Phi_ijk:setDependentVars(all_gLLs:unpack())
	return Phi_ijk
end)
-- Phis[i][j][k]:eq(EFEs[i][j][k]'_ab'^2)

local sum = 0
for i=-n,n do
	for j=-n,n do
		for k=-n,n do
			--sum = sum + frac(1,2) * EFEs[i+n+1][j+n+1][k+n+1]'_ab'^2
			sum = sum + Phis[i+n+1][j+n+1][k+n+1]
		end
	end
end
sum = sum()
local Phi_def = Phi:eq(sum)
printbr(Phi_def)

-- hmm .... have to do this per a,b of g_ab ...
local p = 1
local q = 1
printbr("considering the gradient for", gs[n+1][n+1][n+1][p][q])
local grad_Phi_def = Phi_def:diff(gs[n+1][n+1][n+1][p][q])()
printbr(grad_Phi_def)

local EFEs = Array:lambda(offsets, function(i,j,k)
	local EFE_ijk = var('EFE['..table{i-n-1, j-n-1, k-n-1}:concat', '..']')'_ab':makeDense()
	removeSymmetricDuplicates(EFE_ijk)
	for ab,EFE_ijk_ab in EFE_ijk:iter() do
		EFE_ijk_ab:setDependentVars(all_gLLs:unpack())
	end
	return EFE_ijk
end)
-- this is our energy function to minimize, so I don't think it matters if you double up the symmetries by summing all a,b, or if you just use the upper/lower triangular
local function make_Phi_def(i,j,k)
	local sum = 0
	for a=1,4 do
		for b=a,4 do
			sum = sum + EFEs[i][j][k][a][b]^2
		end
	end
	return Phis[i][j][k]:eq(frac(1,2) * sum)
end
printbr(make_Phi_def(n+1, n+1, n+1))
--printbr(EFEs)

for i=1,2*n+1 do
	for j=1,2*n+1 do
		for k=1,2*n+1 do
			grad_Phi_def = grad_Phi_def:subst(make_Phi_def(i,j,k))
		end
	end
end

printbr(grad_Phi_def)
grad_Phi_def = grad_Phi_def()
printbr(grad_Phi_def)

local function setDependents_D2gLL(T, i, j, k)
	-- [[ 2nd derivative = Linf-norm
	local stencil_gLLs = table()
	for di=-n,n do
		for dj=-n,n do
			for dk=-n,n do
				local g_ijk = index(gs, i+di, j+dj, k+dk)
				if g_ijk then
					for a=1,4 do
						for b=1,4 do
							stencil_gLLs:insert(g_ijk[a][b])
						end
					end
				end
			end
		end
	end
	T:setDependentVars(stencil_gLLs:unpack())
	--]]
end

local Gs = Array:lambda(offsets, function(i,j,k)
	local G_ijk = var('G['..table{i-n-1,j-n-1,k-n-1}:concat', '..']')'_ab':makeDense()
	removeSymmetricDuplicates(G_ijk)
	for ab,G_ijk_ab in G_ijk:iter() do
		setDependents_D2gLL(G_ijk_ab, i, j, k)
	end
	return G_ijk
end)
local Ts = Array:lambda(offsets, function(i,j,k)
	local T_ijk = var('T['..table{i-n-1,j-n-1,k-n-1}:concat', '..']')'_ab':makeDense()
	removeSymmetricDuplicates(T_ijk)
	local deps = table()
	for u=1,4 do
		for v=1,4 do
			deps:insert(gs[i][j][k][u][v])
		end
	end
	for ab,T_ijk_ab in T_ijk:iter() do
		T_ijk_ab:setDependentVars(deps:unpack())
	end
	return T_ijk
end)
for ijkab in Gs:iter() do
	local expr = (EFEs[ijkab]:eq(Gs[ijkab] - 8 * pi * Ts[ijkab])):diff(gs[n+1][n+1][n+1][p][q])()
	grad_Phi_def = grad_Phi_def:subst(expr)
end
printbr(grad_Phi_def)
--symmath.op.mul:pushRule'Expand/apply'
--grad_Phi_def = grad_Phi_def()
--symmath.op.mul:popRule'Expand/apply'
--printbr(grad_Phi_def)

local Rs = Array:lambda(offsets, function(i,j,k)
	local R_ijk = var('R['..table{i-n-1,j-n-1,k-n-1}:concat', '..']')'_ab':makeDense()
	removeSymmetricDuplicates(R_ijk)
	for ab,R_ijk_ab in R_ijk:iter() do
		setDependents_D2gLL(R_ijk_ab, i, j, k)
		--R_ijk[ab]:setDependentVars(gs[i][j][k][ab])
		--R_ijk[ab]:setDependentVars(all_gLLs:unpack())
	end
	return R_ijk
end)
for i=1,2*n+1 do
	for j=1,2*n+1 do
		for k=1,2*n+1 do
			local ijk = {i,j,k}
			-- if you use this again, then consider caching it 
			local Gaussian_ijk = (gUs[i][j][k]'^uv' * Rs[i][j][k]'_uv')()
			for a=1,4 do
				for b=a,4 do	-- do I have to also subst symmetric or no?
					--printbr(Gs[ijkab]:eq(Rs[ijkab]  - frac(1,2) * gs[ijkab] * gUs[ijk]'^uv' * Rs[ijk]'_uv'))
					grad_Phi_def = grad_Phi_def:subst(
						Gs[i][j][k][a][b]:eq(Rs[i][j][k][a][b]  - frac(1,2) * gs[i][j][k][a][b] * Gaussian_ijk)
					)
				end
			end
		end
	end
end
printbr(grad_Phi_def)
-- TODO i don't know that d/dg_ab g^ab works anymore ...
-- now this is going too slow ...
grad_Phi_def = grad_Phi_def:simplifyAddMulDiv()
printbr(grad_Phi_def)
do return end

-- ok we might need to be expanding all indexes ...
-- because here I have to sum the particular change-in-grid-index with its associated product's tensor-index
-- i.e. f^a g_,a means Sum_j f^j[i] g[i+dx^j] - f^j[i] g[i-dx^j] = f^1[i] g[i1+1,i2] - f^1[i] g[i1-1,i2] + f^2[i] g[i1,i2+1] - f^2[i] g[i1,i2-1]
local function finDiff2(T, ij)
end

-- [[
local delta = Tensor:deltaSymbol()
for ijk,R_ijk in Rs:iter() do
	local g_ijk = gs[ijk]
	grad_Phi_def = grad_Phi_def:subst(
		R_ijk'_ab':eq(
			g_ijk'^uv' * (
				frac(1,2) * (
					  finDiff2(gs, '_ub')'_av'
					+ finDiff2(gs, '_av')'_ub'
					- finDiff2(gs, '_uv')'_ab'
					- finDiff2(gs, '_ab')'_uv'
				)
				--[=[
				+ g_ijk'^ef' * (
-- store Gamma? 
--					makeGamma('_e', '_u', '_b') * makeGamma('_f', '_a', '_v')
--					- makeGamma('_e', '_u', '_v') * makeGamma('_f', '_a', '_b')
-- or calculate it on the fly?
				)
				--]=]
			)
		)
	)

	local function makeGamma(a,b,c)
		return frac(1,2) * (
			finDiff(g(a)(b), i)(c)
			+ finDiff(g(a)(c), i)(b)
			- finDiff(g(b)(c), i)(a)
		)
	end
	grad_Phi_def = grad_Phi_def:replace(
		Gs[i]'_ab',
		(
			gs[i]'^uv' * delta'^c_a' * delta'^d_b' 
			- frac(1,2) * gs[i]'_ab' * gs[i]'^uv' * gs[i]'^cd'
		) * Rs[i]'_ucvd'
	)
	--[=[
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
printbr(grad_Phi_def)
--]]

do return end

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
