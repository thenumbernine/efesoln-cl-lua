-- generate the constraint error functions
--	this is going slow - requires some symmath optimizations (described in symmath/diffgeom.lua)

local symmath = require 'symmath'
local Tensor = symmath.Tensor
local var = symmath.var

local t,x,y,z = symmath.vars('t', 'x', 'y', 'z')
local coords = table{t,x,y,z}
local spatialCoords = table{x,y,z}
Tensor.coords{
	{variables=coords},
	{variables=spatialCoords, symbols='ijklmn'}, 
}

local alpha = var'alpha'
print('alpha',alpha)
local betaU = Tensor('^i', function(i) return var('betaU'..i) end)
print('beta^i',betaU)
local gammaLL = Tensor('_ij', function(i,j) return var('gamma_'..sym(i,j)) end)
print('gamma_ij',gammaLL)
local det_gammaLL = symmath.Matrix.determinant(gammaLL)
print('gamma', det_gammaLL)
local gammaUU = symmath.Matrix.inverse(gammaLL)
print('gamma^ij',gammaUU)
os.exit()

local betaL = (betaU'^j' * gammaLL'_ij')()
local betaSq = (betaU'^i' * betaL'_i')()

local gLL = Tensor('_ab', function(a,b)
	if a==1 and b==1 then 
		return -alpha^2 + betaSq
	elseif a==1 then
		return betaL[b-1]
	elseif b==1 then
		return betaL[a-1]
	else
		return gammaLL[a-1][b-1]
	end
end)
local gUU = Tensor('^ab', function(a,b)
	if a==1 and b==1 then
		return -1/alpha^2
	elseif a==1 then
		return betaU[b-1] / alpha^2
	elseif b==1 then
		return betaU[a-1] / alpha^2
	else
		return gammaUU[a-1][b-1] - betaU[a-1] * betaU[b-1] / alpha^2
	end
end)

symmath.tostring = require 'symmath.tostring.MultiLine'
local diffgeom = require 'symmath.diffgeom'(gLL, gUU)

local MathJax = require 'symmath.tostring.MathJax'
symmath.tostring = MathJax

local output = table()
diffgeom:print(function(s)
	output:insert(tostring(s)..'<br>\n')
end)

file['efe_gen.html'] = 
	MathJax.header
	..output:concat'\n'
	..MathJax.footer

