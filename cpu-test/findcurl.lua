#!/usr/bin/env luajit
-- find the curl of a given vector field
-- solves curl w = f for w
require 'ext'
local matrix = require 'matrix'
local n = matrix{3,3,3}
local min = matrix{-1,-1,-1}
local max = matrix{1,1,1}
local dx = (max-min):ediv(n)
local x = n:lambda(function(...) return (matrix{...}-.5):emul(dx)+min end)
local f = n:lambda(function() return matrix{1,0,0} end)
local function curl(w)
--print('w',w)	
	local Aw = n:lambda(function(...)
		local i = matrix{...}
		local dw = matrix{3}:lambda(function(j)
			local ip = matrix(i)
			local im = matrix(i)
			-- [[ fixed boundaries
			ip[j] = math.min(ip[j]+1,n[j])
			im[j] = math.max(im[j]-1,1)
			--]]
			--[[ repeating boundaries
			ip[j] = ip[j]%n[j]+1
			im[j] = (im[j]-2)%n[j]+1
			--]]
			local wp = w[ip]
			local wm = w[im]
			return .5 * dx[j] * (wp - wm)
		end)
		return matrix{
			dw[2][3] - dw[3][2],
			dw[3][1] - dw[1][3],
			dw[1][2] - dw[2][1],
		}
	end)
--print('Aw',Aw)
	return Aw
end
local function normalize(v)
	return v / v:norm()
end
--[=[ using gmres ... not working so well
local w = require 'solver.gmres'{
	--x = f,	-- f
	--x = f:size():zeros(),	-- 0
	x = f:size():lambda(function() return math.random()*2-1 end),	-- random seems to work best
	--x = xs,	-- coord
	--x = n:lambda(function(i,j,k) return matrix{ -(j-.5)*dx[2]+min[2], (i-.5)*dx[1]+min[1], (k-.5)*dx[3]+min[3] } end),	-- cheating?
	b = f,
	A = curl:o(normalize),
	errorCallback = function(err,iter)
		io.stderr:write(iter,'\t',err,'\n')
		assert(math.isfinite(err)) 
	end,
	clone = matrix,
	dot = matrix.dot,
	maxiter = 100 * n:prod(),
	restart = 10,
}
--]=]
--[=[  using 1/4pi int(r' in R3 of curl B(r') / |r-r'|)
-- the problem here is, what if B is constant?  then this gives zero, right?
-- not if the curl is of (B(r') / |r-r'|) ... does it?
for i in x:iter() do
	local xi = x[i] 
	for j in x:iter() do
		local xj = x[j]
	end
end
--]=]
-- [=[
-- note for f1,f2,f3 div-free vector field, one solution is 
--	g1 = int(t=0,1 of t z f2(tx,ty,tz) - t y f3(tx,ty,tz)
--	g2 = int(t=0,1 of t x f3(tx,ty,tz) - t z f1(tx,ty,tz)
--	g3 = int(t=0,1 of t y f1(tx,ty,tz) - t x f2(tx,ty,tz)
-- curl [g1,g2,g3] = [f1,f2,f3]
for i in x:iter() do
	local xi = x[i]
	-- now the integral from the origin to here of the vector field does have the curl we're looking for
	local samples = 1000
	local center = .5*(max - min)
	local ic = matrix{3}:zeros()
	for j=1,samples do
		local t = (j-.5)/samples
		local x2 = ((t * (xi - center) + .5)):map(function(y,k) return math.clamp(math.floor(y),1,n[k]) end)
	end
end
--]=]
assert(w)
local f2 = curl(w)
for i in x:iter() do
	print(i,j,k,f[i],w[i],f2[i])
	return 0
end
