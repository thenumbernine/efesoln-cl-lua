#!/usr/bin/env luajit
require 'ext'
local matrix = require 'matrix'
matrix.__tostring = tolua
local ns = matrix{8,8,8}
print('ns',ns)
local max = matrix{1,1,1}
local min = -max
local dxs = (max - min):ediv(ns)
print('min',min)
print('max',max)
print('dxs',dxs)
local xs = ns:lambda(function(...)
	return (matrix{...} - .5):emul(dxs) + min
end)
--print('xs',xs)
local eps3 = matrix{3,3,3}:lambda(function(i,j,k)
	return (i%3+1==j and j%3+1==k) and 1 or
		((k%3+1==j and j%3+1==i) and -1 or 0)
end)
-- x cross y = eps3 * y * x

local function cross(a,b) 
	return eps3 * b * a
end

--[[
units
once again
using Planck units

space vs time:
speed of light: c = 1 = ~c~ m/s, for ~c~ = 2.99792458e+8
so ~c~ m = 1 s, 1/~c~ s = 1 m, m/s = 1/~c~
--]]
local c_const = 2.99792458e+8
print('c_const',c_const)  
local s_in_m = 1 / c_const
print('s_in_m',s_in_m) 

-- why not just change the lenght of a meter to 1/3e+8 ?
-- here is how much it would shrink by ...
--print((1/c_const - 1/3e+8)*c_const*1000)	-- difference in 1 m and 1 adj m, in terms of mm
-- 0.69180666666663 mm 

--[[
space + time vs mass:
G = 1 = ~G~ m^3 / (kg s^2), for ~G~ = 6.6740831e-11
so 1 = ~G~ m^3 / (kg s^2)
1 kg = ~G~ (m^2 / s^2) m
1 kg = ~G~/~c~^2 m
--]]
local G_const = 6.6740831e-11
print('G_const',G_const)  
local kg_in_m = G_const / c_const^2
print('kg_in_m',kg_in_m)  

--[[
space + time + mass vs charge:
ke = 1 = ~ke~ kg m^3 / (s^2 C^2), for ~ke~ = 8.9875517873681764e+9
so 1 = ~ke~ (~G~/~c~^2 m) m (1/~c~)^2 / C^2
C^2 = ~ke~ ~G~ / ~c~^4 m^2
1 C = sqrt(~ke~ ~G~) / ~c~^2 m
--]]
local ke_const = 8.9875517873681764e+9
print('ke_const',ke_const)  
local C_in_m = math.sqrt(ke_const * G_const) / c_const^2	-- m
print('C_in_m',C_in_m)  
local N_in_m = kg_in_m / s_in_m^2	-- m^0
print('N_in_m',N_in_m)  
local V_in_m = N_in_m / C_in_m	-- m^0
print('V_in_m',V_in_m)  
local Ohm_in_m = kg_in_m / (s_in_m * C_in_m^2)	-- m^0
print('Ohm_in_m',Ohm_in_m)  

--[[
ke = 1 = 1 / (4 pi eps0) <=> eps0 = 1 / (4 pi)
c^2 = 1 / (mu0 eps0) <=> mu0 = 1/eps0 = 4 pi
--]]
local mu0 = 4 * math.pi	-- m^0
print('mu0',mu0)  
local eps0 = 1 / mu0	-- m^0
print('eps0',eps0)  

--[[
1 C = ~e~ e, for e = unit of charge for an electron, and ~e~ = 6.2415093414e+18
so 1 e = 1/~e~ C = sqrt(~ke~ ~G~) / (~c~^2 ~e~) m
notice that the classical electron radius is re = ke e^2 / (me c^2) = 2.817940322719e-15 m
--]]
local e_const = 6.2415093414e+18
print('e_const',e_const)  
local e_in_m = C_in_m / e_const
print('e_in_m',e_in_m)

local h_in_m = 1.61622938e-35
print('h_in_m',h_in_m)  
-- note the ratio of the elementary charge-converted-to-meters and Planck length 
-- is the same as the ratio of the elementary charge and Planck charge
local alpha_const = (e_in_m/h_in_m)^2
print('alpha_const',alpha_const)

--[[
local A = ns:lambda(function(...)
	local i = matrix{...}
	local x,y,z = table.unpack(xs[i])	
	return matrix{
	}
end)
--]]

--[[ flat space and constant fields
local E = ns:lambda(function(...) return matrix{1,0,0} end)
local B = ns:lambda(function(...) return matrix{0,1,0} end)
--]]
-- [[ field around a wire
-- small AA battery has 1.5 Volts
-- current in amperes = voltage / resistance 

-- source: http://hyperphysics.phy-astr.gsu.edu/hbase/electric/resis.html
local wire_resistivities = table{	-- at 20' Celsius, in ohm m
	aluminum = 2.65e-8,
	copper = 1.724e-8,
	iron = 9.71e-8,
	nichrome = 1e-6,
	gold = 2.24e-8,
	silver = 1.59e-8,
	platinum = 1.06e-7,
	tungsten = 5.65e-8,
}:map(function(v) return v * Ohm_in_m end)	-- ohm m => m^0
local in_in_m = .0254 
local wire_diameters = table{	-- starts in inches
	electrical_range = .1019,
	household_circuit = .0808,
	switch_leads = .0640,
}:map(function(v) return v * in_in_m end)	-- in => m 
local wire_radius = .5 * wire_diameters.electrical_range	-- m
print('wire_radius',wire_radius)
local wire_cross_section_area = math.pi * wire_radius^2	-- m^2
print('wire_cross_section_area',wire_cross_section_area)  
local wire_length = 12 * in_in_m	-- m
print('wire_length',wire_length)  
local wire_resistivity = wire_resistivities.gold
print('wire_resistivity',wire_resistivity)  
local wire_resistance = wire_resistivity * wire_length / wire_cross_section_area	-- m^0
print('wire_resistance',wire_resistance)  
local battery_voltage = 1.5 * V_in_m	-- m^0
print('battery_voltage',battery_voltage)  
local wire_current = battery_voltage / wire_resistance	-- amps = C / s = m / m = m^0, likewise volts = m^0, ohms = m^0, so amps = volts / ohms = m^0
print('wire_current',wire_current)  
local current_velocity = 1	-- doesn't matter if lambda = 0.  units of m/s = m^0
-- so ... inside the wire we know q=0 by Kirchoff's law
-- what about J = sigma E? doesn't that mean E = J rho, for rho = resistivity?
local wire_charge_density = 0	-- C / m^3 = m^-2
print('wire_charge_density',wire_charge_density)  
local wire_charge_density_per_length = wire_charge_density * wire_cross_section_area	-- m^-2 * m^2 = m^0
print('wire_charge_density_per_length',wire_charge_density_per_length)  
local wire_surface_charge_density = 0
print('wire_surface_charge_density',wire_surface_charge_density)
local E = ns:lambda(function(...)
	local x,y,z = table.unpack(xs(...))	-- m
	local r2 = math.sqrt(x*x+y*y)	-- m
	-- hyperphysics:
	--local lambda = wire_charge_density_per_length / current_velocity	-- m^0 / m^0 = m^0
	--local Er = lambda / (2 * math.pi) / (eps0 * r2)	-- m^0 / m = m^-1
	-- https://physics.stackexchange.com/questions/291779/electric-field-outside-wire-with-stationary-current?rq=1 
	-- J = sigma E
	local Ez_int = wire_current * wire_resistivity	-- m^0 * m^0 = m^0
	local Er_int = 0
	local Er = wire_surface_charge_density * wire_radius / (eps0 * r2)	-- m^0 * m / m = m^0
	local Ez = wire_current * wire_resistivity	-- m^0
	return matrix{x/r2*Er,y/r2*Er,Ez} * math.sqrt(eps0)
	-- times sqrt(eps0) for convenience of representing T_ab
end)
print('E',E)
local B = ns:lambda(function(...)
	local x,y,z = table.unpack(xs(...))	-- m
	local r2 = math.sqrt(x*x+y*y)	-- m
	-- http://www.ifi.unicamp.br/~assis/Found-Phys-V29-p729-753(1999).pdf
	local Bt_int = mu0 * wire_current * r2 / (2 * math.pi * wire_radius^2) -- m^0 m / (m^2) = m^-1
	local Bt = mu0 * wire_current / (2 * math.pi * r2)	-- m^-1
	--print('Bt',Bt)
	return matrix{-y/r2 * Bt, x/r2 * Bt, 0} / math.sqrt(mu0)
	-- divide sqrt(mu0) for convenience ... same deal
end)
print('B',B)
--]]

local R = ns:lambda(function(...)
	local i = matrix{...}
	local Ei = E[i]
	local Bi = B[i]
	local Si = cross(Ei,Bi)
	return matrix{4,4,4,4}:lambda(function(a,b,c,d)
		local s = 1
		if a>b then a,b,s = b,a,-s end
		if c>d then c,d,s = d,c,-s end
		if a==1 and b>1 and c==1 and d>1 then
			-- R^t_itj = -E_i E_j - B_i B_j
			return -s * (Ei[b-1] * Ei[d-1] + Bi[b-1] * Bi[d-1])
		elseif a==1 and b>1 and c>1 and d>1 then
			-- R^t_ijk = gamma_ij S_k - gamma_ik S_j
			return s * (((b==c) and Si[d-1] or 0) - ((b==d) and Si[c-1] or 0))
		elseif a>1 and b>1 and c>1 and d>1 then
			-- R^i_jkl = eps^i_jm (E^m E_n + B^m B_n) eps^n_jk
			local sum = 0
			for m=1,3 do
				for n=1,3 do
					sum = sum + eps3[a-1][b-1][m] * (Ei[m] * Ei[n] + Bi[m] * Bi[n]) * eps3[n][c-1][d-1]
				end
			end
			return sum
		end
		return 0
	end)
end)
print('R',R)
--[[
R^a_bcd = Conn^a_bd,c - Conn^a_bc,d + Conn^a_ec Conn^e_bd - Conn^a_ed Conn^e_bc
Conn^a_bd,c - Conn^a_bc,d = R^a_bcd - Conn^a_ec Conn^e_bd + Conn^a_ed Conn^e_bc
c==d: Conn^a_ec Conn^e_bd - Conn^a_ed Conn^e_bc = 0
c==t: Conn^a_bi,t - Conn^a_bt,i = R^a_bti - Conn^a_et Conn^e_bi + Conn^a_ei Conn^e_bt
d/dt vec Conn^a_b - grad phi Conn^a_b = vec R^a_bt + Conn^2 terms
d/dt vec Conn^i_t - grad phi Conn^i_t = vec R^i_tt + Conn^2 terms 
	= vec R^t_it + Conn^2 terms
	= -E_i vec E - B_i vec B + Conn^2 terms
	let's say Conn^i_t = vec A
		then d/dt vec Conn^i_t - grad phi Conn^i_t = -vec E
	vec E = d/dt vec A - grad A_t
	E_i vec E = E_i d/dt vec A - E_i grad A_t
	vec B = curl vec A
	B_i vec B = B_i curl vec A
	E_i vec E + B_i vec B = E_i d/dt vec A - E_i grad A_t + B_i curl vec A
	E_i E_j + B_i B_j = E_i A_j,t - E_i A_t,j + B_i eps^jkl A_l,k
		
				-E_i A_t,j
	E_i A_j,t	B_i eps^jkl A_l,k

		= E_i (A^j_,t + A^t_,j) + B_i eps^jkl A_l,k
		using A^u_,u = 0
		= E_i (A^j_,t + A^t_,j + A^t_,t + A^j_,j) + B_i eps^jkl A_l,k


c!=d: Conn^a_bj,i - Conn^a_bi,j = R^a_bij - Conn^a_ei Conn^e_bj + Conn^a_ej Conn^e_bi
curl vec Conn^a_b = dual vec R^a_b + Conn^2 terms

I'm looking for u'^i = -Conn^i_tt (u^t)^2

Conn^t_ij,t - Conn^t_it,j = R^t_itj + Conn^2
Conn^i_tj,t - Conn^i_tt,j = R^i_ttj + Conn^2

Conn^t_ik,j - Conn^t_ij,k = R^t_ijk + Conn^2
Conn^i_tk,j - Conn^i_tj,k = R^i_tjk + Conn^2
Conn^i_jl,k - Conn^i_jk,l = R^i_jkl + Conn^2

dtConn[{x,y,z, 1,i,j}] - diff(j,Conn[{x,y,z, 1,i,1}]) = R[{x,y,z, 1,i,1,j}] + Conn^2
dtConn[{x,y,z, i,1,j}] - diff(j,Conn[{x,y,z, i,1,1}]) = R[{x,y,z, i,1,1,j}] + Conn^2

diff(j,Conn[{x,y,z, 1,i,k}]) - diff(k,Conn[{x,y,z, 1,i,j}]) = R[{x,y,z, 1,i,j,k}] + Conn^2
diff(j,Conn[{x,y,z, i,1,k}]) - diff(k,Conn[{x,y,z, i,1,j}]) = R[{x,y,z, i,1,j,k}] + Conn^2
diff(k,Conn[{x,y,z, i,j,l}]) - diff(l,Conn[{x,y,z, i,j,k}]) = R[{x,y,z, i,j,k,l}] + Conn^2

when neglecting Conn^2, each of these can be solved separate:
*) Conn^i_jk
*) Conn^t_ij and Conn^i_tj
*) Conn^i_tt and Conn^t_it
--]]
--[=[ R3 separately
local R3 = ns:lambda(function(x,y,z)
	local Ri = R[x][y][z]
	--[[
	-- R3[x,y,z,i,j] = epsilon^ikl epsilon^jmn R_klmn
	-- so R3[xyz] will hold xx = R_yzyz, xy = R_yzzx, etc
	-- so we have as many variables as degrees of freedom
	return matrix{3,3}:lambda(function(i,j)
		local i2 = i%3+1
		local i3 = i2%3+1
		local j2 = j%3+1
		local j3 = j2%3+1
		return Ri[i2+1][i3+1][j2+1][j3+1]
	end)
	--]]
	-- [[
	-- just return the whole antisymmetric thing
	-- so we have as many variables as we want to solve for
	return matrix{3,3,3,3}:lambda(function(i,j,k,l)
		return Ri[i+1][j+1][k+1][l+1]
	end)
	--]]
end)
print('R3',R3)
-- now we can solve for Conn[x,y,z]^i_jl,k - Conn[x,y,z]^i_jk,l = R3[x,y,z
-- or can we?  R3 has 3^4 variables, Conn3 has 3^3 ...
-- if we just look at free parameters, R3 has 3^2 and Conn3 has n*n*(n+1)/2 = 18
-- hmm, maybe I'll have to solve more things at once?
-- full equation: 
-- R has n^2 (n^2 - 1) / 12 = 16*15/12 = 20 free parameters
-- Conn^a_bc has n*n*(n+1)/2 = 40 free parameters
local function AConn3(Conn3)
	local result = ns:lambda(function(...)
		local xs = matrix{...}
		local Conn3i = Conn3[xs]
		assert(type(Conn3i[1][1][1]) == 'number')
		local dxConn3i = matrix{3}:lambda(function(i)
			local xp = matrix(xs)
			xp[i] = math.min(xp[i] + 1, ns[i])
			local xm = matrix(xs)
			xm[i] = math.max(xm[i] - 1, 1)
			return .5 / dxs[i] * (Conn3[xp[1]][xp[2]][xp[3]] 
								- Conn3[xm[1]][xm[2]][xm[3]])
		end)
		-- dxConn3[l][i][j][k] = Conn3^i_jk,l
		return matrix{3,3,3,3}:lambda(function(i,j,k,l)
			local Conn3Sq = 0
			--[[ notice the square term is neglecting the partial_t component 
			for m=1,3 do
				Conn3Sq = Conn3Sq + Conn3i[i][m][k] * Conn3i[m][j][l]
								- Conn3i[i][m][l] * Conn3i[m][j][k]
			end
			--]]
			return dxConn3i[k][i][j][l]
				- dxConn3i[l][i][j][k]
				+ Conn3Sq
		end)
	end)
	local diff = (result - Conn3):normSq()
	-- hmm how does this function always work out to be identity?
	print('calculating R3 of Conn3, got |R3-Conn3| =',diff)
	return result
end
local Conn3 = require 'solver.gmres'{
	b = R3,
	x = R3:size():zeros(),
	A = AConn3,
	clone = matrix,
	dot = matrix.dot,
	errorCallback = function(err, iter)
		io.stderr:write('err ',err,' iter ',iter,'\n')
		io.stderr:flush()
	end,
	epsilon = 1e-10,
	maxiter = ns:prod(),
	restart = 10,
}
print('|Conn3-R3|',(Conn3-R3):normSq())
print('|A(Conn3)-R3|', (AConn3(Conn3) - R3):normSq())
--]=]
--[=[ R completely
-- R has 4^4 = 256 components, Conn has 4^3 = 64 components
-- in addition, another 
local function AConn(Conn)
	local result = ns:lambda(function(...)
		local xs = matrix{...}
		local Conni = Conn[xs]
		assert(type(Conni[1][1][1]) == 'number')
		local dxConni = matrix{4}:lambda(function(i)
			if i==1 then return Conn[xs]:size():zeros() end	-- return dtConn[xs] end
			local xp = matrix(xs)
			xp[i-1] = math.min(xp[i-1] + 1, ns[i-1])
			local xm = matrix(xs)
			xm[i-1] = math.max(xm[i-1] - 1, 1)
			return .5 / dxs[i-1] * (Conn[xp[1]][xp[2]][xp[3]] - Conn[xm[1]][xm[2]][xm[3]])
		end)
		-- dxConn[l][i][j][k] = Conn^i_jk,l
		return matrix{4,4,4,4}:lambda(function(i,j,k,l)
			local ConnSq = 0
			-- [[ notice the square term is neglecting the partial_t component 
			for m=1,4 do
				ConnSq = ConnSq + Conni[i][m][k] * Conni[m][j][l]
								- Conni[i][m][l] * Conni[m][j][k]
			end
			--]]
			return dxConni[k][i][j][l] - dxConni[l][i][j][k] + ConnSq
		end)
	end)
	local diff = (result - Conn):normSq()
	-- hmm how does this function always work out to be identity?
	print('calculating R of Conn, got |R-Conn| =',diff)
	return result
end
local Conn = 
	require 'solver.gmres'{
--	require 'solver.conjres'{
	b = R,
	x = R:size():ones(),
	A = AConn,
	clone = matrix,
	dot = matrix.dot,
	errorCallback = function(err, iter)
		io.stderr:write('err ',err,' iter ',iter,'\n')
		io.stderr:flush()
	end,
	epsilon = 1e-10,
	maxiter = R:size():prod(),
	restart = 10,
}
print('Conn',Conn)
print('|Conn-R|',(Conn-R):normSq())
print('|A(Conn)-R|', (AConn(Conn) - R):normSq())
--]=]
--[[
but in both these cases, Conn has a different dimension than R
so you need to solve for a pseudoinverse
and transposing A is tough when A is a sparse linear function, not even a sparse matrix
so how about another approach? 
minimize Conn wrt |A(Conn)-R|^2 ?
pseudoinverse?
how do you find the pseudoinverse of a linear function of finite differences?
--]]
--[=[ 
local function AConn(Conn)
	local result = ns:lambda(function(...)
		local xs = matrix{...}
		local Conni = Conn[xs]
		assert(type(Conni[1][1][1]) == 'number')
		local dxConni = matrix{4}:lambda(function(i)
			-- Conn^a_bc,x <=> A[index(x,abc)][index(x-1,abc)] = -h/2, A[index(x,abc)][index(x+1,abc)] = h/2
			-- Conn^a_bc,i <=> A[index(x_i,abc)][index(x_i-1,abc)] = -h/2, A[index(x_i,abc)][index(x_i+1,abc)] = h/2
			-- Conn^a_bi,j - Conn^abj,i = R^a_bji <=> 	A[index(x,abji)][index(x-dxj,abi)] = -h/2, A[index(x,abji)][index(x+dxj,abi)] = h/2
			--											A[index(x,abji)][index(x-dxs,abi)] = -h/2, A[index(x,abji)][index(x+dxs,abi)] = h/2
			-- now find AA^t
			diff(j,Conn[{x,y,z, 1,i,k}]) - diff(k,Conn[{x,y,z, 1,i,j}]) = R[{x,y,z, 1,i,j,k}] + Conn^2
			diff(j,Conn[{x,y,z, i,1,k}]) - diff(k,Conn[{x,y,z, i,1,j}]) = R[{x,y,z, i,1,j,k}] + Conn^2
			diff(k,Conn[{x,y,z, i,j,l}]) - diff(l,Conn[{x,y,z, i,j,k}]) = R[{x,y,z, i,j,k,l}] + Conn^2
		
			if i==1 then return Conn[xs]:size():zeros() end	-- return dtConn[xs] end
			local xp = matrix(xs)
			xp[i-1] = math.min(xp[i-1] + 1, ns[i-1])
			local xm = matrix(xs)
			xm[i-1] = math.max(xm[i-1] - 1, 1)
			return .5 / dxs[i-1] * (Conn[xp[1]][xp[2]][xp[3]] - Conn[xm[1]][xm[2]][xm[3]])
		end)
		-- dxConn[l][i][j][k] = Conn^i_jk,l
		return matrix{4,4,4,4}:lambda(function(i,j,k,l)
			local ConnSq = 0
			-- [[ notice the square term is neglecting the partial_t component 
			for m=1,4 do
				ConnSq = ConnSq + Conni[i][m][k] * Conni[m][j][l]
								- Conni[i][m][l] * Conni[m][j][k]
			end
			--]]
			return dxConni[k][i][j][l] - dxConni[l][i][j][k] + ConnSq
		end)
	end)
	local diff = (result - Conn):normSq()
	-- hmm how does this function always work out to be identity?
	print('calculating R of Conn, got |R-Conn| =',diff)
	return result
end

local Conn = matrix.zeros(ns[1],ns[2],ns[3], 4,4,4)
print('Conn',Conn)
local R0 = AConn(Conn)
assert(R0:size() == R:size())
local phi = (R0 - R):map(function(x) return x*x end)	-- guess
print('|A(Conn)-R|', (R0 - R):normSq())
print('phi',phi)
print('R0 volume',R0:size():prod())
print('Conn volume',Conn:size():prod())
--]=]

--[[
here's a finite-difference solution, neglecting square terms
Conn^j_tt,i = Conn^j_ti,t - R^j_tti + Conn^j_at Conn^a_ti - Conn^j_ai Conn^a_tt
Conn^j_tt,ii = (Conn^j_ti,t - R^j_tti + Conn^j_at Conn^a_ti - Conn^j_ai Conn^a_tt),i
Conn^j_tt = lap^-1 grad_i (Conn^j_ti,t - R^j_tti + Conn^j_at Conn^a_ti - Conn^j_ai Conn^a_tt)

R^a_bti = Conn^a_bi,t - Conn^a_bt,i + Conn^a_et Conn^e_bi - Conn^a_ei Conn^e_bt
Conn^a_bt,i = Conn^a_bi,t - R^a_bti + Conn^a_et Conn^e_bi - Conn^a_ei Conn^e_bt
Conn^a_bt,ii = (Conn^a_bi,t - R^a_bti + Conn^a_et Conn^e_bi - Conn^a_ei Conn^e_bt),i
this can be used to solve ...
Conn^t_tt based on R^t_tti = g^tk R_ktti = -g^tk (E_k E_i + B_k B_i)
Conn^j_tt based on R^j_tti = g^jk R_ktti = -g^jk (E_k E_i + B_k B_i)
Conn^t_jt based on R^t_jti = g^tt R_tjti + g^tk R_kjti = g^tt (E_j E_i + B_j B_i) + g^tk (gamma_ik S_j - gamma_ij S_k)
Conn^j_kt based on R^j_kti = g^jt R_tkti + g^jl R_lkti = g^jk (E_k E_i + B_k B_i) + g^jl (gamma_il S_k - gamma_ik S_l)

by symmetry we can calculating the following, assuming negligible commutations / torsion (which is probably not true):
Conn^a_bc - Conn^a_cb = c_cb^a <=> Conn^a_bc = Conn^a_cb + c_cb^a
Conn^t_tj = Conn^t_jt + c_jt^t	<- which can be derived from above
Conn^k_tj = Conn^k_jt + c_jt^k	<- "
Conn^t_kj = Conn^t_jk + c_jk^t	<- which has to be solved separately
Conn^k_lj = Conn^k_jl + c_jl^k
--]]
-- [=[ this is getting finite results!  hooray!
local dt_Conn = matrix{ns[1],ns[2],ns[3],4,4,4}:zeros()
local Conn = matrix{ns[1],ns[2],ns[3],4,4,4}:zeros()
local last_Conn_norm 
for iter=1,20 do
	-- for each a,b calculate Conn^a_bt(x,y,z), stored as [a,b,x,y,z]
	local ConnUa_bt = matrix{4,4}:lambda(function(a,b)	
		-- returns Conn^a_bt,i = Conn^a_bi,t - R^a_bti + Conn^a_et Conn^e_bi - Conn^a_ei Conn^e_bt
		local function eval_ConnUa_bt_i(x,y,z,i)
			local ConnSq = 0
			-- [[ ConnSq terms
			local Conni = Conn[x][y][z]
			for c=1,4 do
				ConnSq = ConnSq + Conni[a][c][1] * Conni[c][b][i+1]
								- Conni[a][c][i+1] * Conni[c][b][1]
			end
			--]]
			local dt_Connabi = dt_Conn[x][y][z][a][b][i+1]
			-- technically the E^j is raised by g, and g determined from Conn
			-- but I'm going to assume g ~ I for now
			local Rabti = R[x][y][z][a][b][1][i+1]
			return dt_Connabi - Rabti + ConnSq
		end

		local lapConnUa_bt = ns:lambda(function(x,y,z)
			if x == 1 or x == ns[1]
			or y == 1 or y == ns[2]
			or z == 1 or z == ns[3]
			then
				return 0
			end
			return (.5 / dxs[1]) * (eval_ConnUa_bt_i(x+1,y,z,1) - eval_ConnUa_bt_i(x-1,y,z,1))
				+ (.5 / dxs[2]) * (eval_ConnUa_bt_i(x,y+1,z,2) - eval_ConnUa_bt_i(x,y-1,z,2))
				+ (.5 / dxs[2]) * (eval_ConnUa_bt_i(x,y,z+1,3) - eval_ConnUa_bt_i(x,y,z-1,3))
		end)
		-- now solve the inverse of the divergence
		-- and the gradient of that will be the inverse of the divergence

		return require 'solver.gmres'{
			x = -lapConnUa_bt,
			b = lapConnUa_bt,
			A = function(y)
				return ns:lambda(function(i,j,k)
					if i==1 or i==ns[1] or j==1 or j==ns[2] or k==1 or k==ns[3] then
						return 0
					else
						local h2 = dxs[1]^2	-- TODO this for a non-square grid
						return (y[i+1][j][k]
								+ y[i-1][j][k]
								+ y[i][j+1][k]
								+ y[i][j-1][k]
								+ y[i][j][k+1]
								+ y[i][j][k-1]
								- 6 * y[i][j][k]) / h2
					end
				end)	
			end,
			clone = matrix,
			dot = matrix.dot,
			errorCallback = function(err, iter, x, rSq, bSq)
				-- err varies from algorithm to algorithm ... hmm ... maybe it shouldn't ... 
				print('...gmres of Conn',a,b,'iter',iter, 'err',err)
			end,
			epsilon = 1e-30,
			maxiter = 4 * ns:prod(),
			restart = 50,	--gmres-only
		}
	end)
	-- now swizzle to update the Conn^a_bt(x,y,z) components
	for i in ConnUa_bt:iter() do
		local a,b,x,y,z = table.unpack(i)
		local ConnUa_bt_xi = ConnUa_bt[i]
		Conn[x][y][z][a][b][1] = ConnUa_bt_xi 
		Conn[x][y][z][a][1][b] = ConnUa_bt_xi
	end
	local Conn_norm = Conn:norm()
	local delta_Conn_norm = last_Conn_norm and (Conn_norm - last_Conn_norm) or nil
	print('iter',iter,'|Conn|',Conn_norm,'delta|Conn|',delta_Conn_norm)
	if delta_Conn_norm and math.abs(delta_Conn_norm) < 1e-50 then break end
	last_Conn_norm = Conn_norm
end
--]=] 
--[[
now to solve those last two
R^t_ijk = Conn^t_ik,j - Conn^t_ij,k + Conn^t_aj Conn^a_ik - Conn^t_ak Conn^a_ij
1/2 eps^mjk R^t_ijk = eps^mjk Conn^t_ik,j + eps^mjk Conn^t_aj Conn^a_ik
eps^mjk Conn^t_ik,j = 1/2 eps^mjk R^t_ijk - eps^mjk Conn^t_aj Conn^a_ik
curl Conn^t_ik, using the k index as the vector = 
	the dual of the jk of R^t_ijk - eps^mjk Conn^t_aj Conn^a_ik
Conn^t_ij based on R^t_ijk = g^tt R_tijk + g^tl R_lijk = g^tt (gamma_ij S_k - gamma_ik S_j) + g^tl e_lim (E^m E^n + B^m B^n) e_mjk
Conn^i_jk based on R^i_jkl = g^it R_tjkl + g^im R_mjkl = g^it (gamma_jk S_l - gamma_jl S_k) + g^im e_mjp (E^p E^q + B^p B^q) e_pkl


... so how do you perform a discrete inverse of div and curl?
https://groups.google.com/forum/#!topic/comp.soft-sys.matlab/jv_2gsSF-pE
curl B = R
div B = D
(A,phi) = Hemholtz decomposition of B
B = curl A + grad phi
first the curl ...
div A = 0 <= gauge
curl B = R <=> curl(curl A + grad phi) = R <=> curl curl A = R <=> grad div A - lap A = R <=> -lap A = R
=> A = -veclap^-1 R ... lap is applied per-component
next the div ...
div B = D <=> div (curl A + grad phi) = D <=> div grad phi = lap phi = D
=> phi = lap^-1 D
B = curl A + grad phi = grad lap^-1 D - curl veclap^-1 R

Conn^t_tt is solved by div only
Conn^k_tt is solved by div only
Conn^t_tj is solved by div & curl
Conn^k_tj is solved by div & curl
Conn^t_ij is solved by curl only
Conn^k_ij is solved by curl only

but what if we don't separate out the commutation/torsion just yet ...
Conn^t_tt is solved by div
Conn^k_tt is solved by div
Conn^t_jt is solved by div
Conn^k_jt is solved by div
Conn^t_tj is solved by div & curl
Conn^k_tj is solved by div & curl
Conn^t_ij is solved by div & curl
Conn^k_ij is solved by div & curl

hmm, running too slow as well, only getting a 8x8x8 to perform decently
time to switch over to GPU ...
... which is already implemented in efesoln ...
... why not implement this algo there?
--]]
local f = assert(io.open('out.txt','w'))
f:write'#ix\tiy\tiz\tgrav\n'
for i=1,ns[1] do
	for j=1,ns[2] do
		for k=1,ns[3] do
			local x,y,z = table.unpack(xs[i][j][k])
			local Conni = Conn[i][j][k]
			local gx, gy, gz = Conni[2][1][1], Conni[3][1][1], Conni[4][1][1]
			local grav = math.sqrt(gx*gx + gy*gy + gz*gz)
			f:write(x,'\t',y,'\t',z,'\t',grav,'\n')
		end
	end
end
f:close()
