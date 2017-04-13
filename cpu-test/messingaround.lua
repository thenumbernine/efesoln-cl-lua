#!/usr/bin/env lua
-- lua 5.3 for utf8
--[[
representing a number takes two things:
1) power spectrum - what base the digits are 10: 1,10,100,1000...
2) coefficients as to create linear combinations of that base: 0,1,2,3,4,5,6,7,8,9, ... 
what would happen if you removed the number 5 from the base-10 counting system?
can you still form linear combinations of {0..9}/5 and 10^k to represent numbers with 5?
how about 5 itself?
hmm not without introducting alternatively a digit for 1/2 (applied to next power of base) or a digit for 50 (apply to previous power of base)
--]]

require 'ext'
local numbermeta = debug.getmetatable(0)
numbermeta.__tostring = numbermeta.tostring

numbermeta.chomp = function(t,k)
	local base = numbermeta.base
	local p = base^k
	local a = math.floor(t / (base*p))
	local b = math.floor(t / p) - base * a
	local c = t % p 
	return c,b,a
end
numbermeta.build = function(c,b,a,k)
	local base = numbermeta.base
	return c + base^k * (b + base * a)
end
-- parenthesis mean implicit multiply =D
numbermeta.__call = function(a,b)	
	return a * b
end
numbermeta.__index = function(t,k)
	assert(type(t) == 'number')
	local r
	if type(k) == 'table' then
		k,r = table.unpack(k)
	end
	if type(k) == 'string' then return numbermeta[k] end
	assert(type(k) == 'number')
	assert(r == nil or type(r) == 'number')
	-- return the k'th digit ...
	local c,b,a = numbermeta.chomp(t,k)
	if r then
		return numbermeta.build(c,r,a,k)
	else
		return b
	end
end

print((.5):__tostring(10.5))
print((10.5):__tostring(11))
--[[
print(math.pi:__tostring(math.exp(1)))
print(math.exp(1):__tostring(math.pi))
print((math.pi*math.pi):__tostring(math.pi))
print((10):__tostring(math.pi))
print((2187):__tostring(3))
print((1827):__tostring(3))
--]]
--[[
print((-100):__tostring(10))
print((-100.1):__tostring(10))
print((-100.01):__tostring(10))
print((-100.001):__tostring(10))
print((-0):__tostring(10))
print((-1):__tostring(10))
print((-0.1):__tostring(10))
print((-0.01):__tostring(10))
print((-0.001):__tostring(10))
os.exit()
--]]
--[[
local n = 123454321
for i=2,50 do
	print('base '..i..': '..n:__tostring(i))
end
os.exit()
--]]
--[[
for i=0,8 do 
	print('digit '..i..' is '..n[i])
end
print(n)
print(n[{3,9}])	-- replace operator... nasty
--]]
--[[
number abuse
overload the index/newindex of numbers
what should it mean?
what does (1)[2] represent?
the 2nd digit of 1?
the 2nd modulo of 1?
the 2nd prime factor?
(1)[2] = 1 returns 101 ...
--]]
