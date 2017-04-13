#!/usr/bin/env luajit
-- solve an underconstrained system
-- x - y = c
-- x - y - c = 0
require 'ext'
local matrix = require 'matrix'
local symmath = require 'symmath'
local x,y = symmath.vars('x','y')
--[[
for multiple constraints ... a_ij x_j = y_i
i could use a pseudoinverse but then I'd have to to write the linear function as a transpose ...
what if I wrote it as a newton function?
		phi_i = (sum j (a_ij x_j) - y_i)^2
	d/dx_j phi_i = 2 a_ij (sum_k a_ik x_k - y_k) 
		= 2 a_ij (A(x)-y)_i	<- this is just A (Ax-b)^t, so practically a pseudo inverse
	so dx_i/dt = -sum_j dphi_j/dx_i ... sum across all potentials, but what if some cancel and produce a local minima?
		= -sum_j (2 a_ji (sum_k a_jk x_k - y_k))
		= -2 A^t (A(x)-b)
	once again, for implementation, a transpose is needed ...
--]]
local c = 2
local phi = (x-y-c)^2
-- but this introduces a fixed souce between |x-y|=+c and -c, at |x-y|=0
local dphi_dx = phi:diff(x)()
local dphi_dy = phi:diff(y)()
print('phi',phi)
print('dphi_dx',dphi_dx)
print('dphi_dy',dphi_dy)
local dphi_dx_f = dphi_dx:compile{x,y}
local dphi_dy_f = dphi_dy:compile{x,y}
local dphi_dp = function(p)
	return matrix{dphi_dx_f(p[1],p[2]), dphi_dy_f(p[1],p[2])}
end
--local p = matrix{2,1}	--> 2.5, .5
--local p = matrix{1.5,1.5}	--> 2.5, .5
--local p = matrix{3,0}	--> 2.5, .5
print(p)
local lambda = .1
while true do
	local dp = -dphi_dp(p)
	print('dp',dp)
	local norm = dp:norm()
	if not math.isfinite(norm) or norm < 1e-10 then break end
	p = p + lambda * dp
	print(p)
end 
