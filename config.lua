return {
	--size = 64,
	size = 32,
	--size = 8,
	
	-- line trace amount
	updateLambda = 1e+9,

	--body = 'vacuum',
	body = 'Earth',
	--body = 'Sun',
	--body = 'EM ring',
	--body = 'EM line',
	--body = 'EM constant',
	
	bodyRadii = 2,
	
	initCond = 'flat',
	--initCond = 'stellar Schwarzschild',
	--initCond = 'stellar Kerr-Newman',

	--boundaryCond = 'g_ab = eta_ab',	-- has bad discontinuities in the gravity calc at the edge
	boundaryCond = 'g_ab,c = 0',		-- looks just as nice as providing the exact stellar-Schwarzschild, however ... Newton gradient descent might run away ...
	--boundaryCond = 'stellar Schwarzschild',
	
	outputFilename = 'out.txt',

	solver = 'Newton',	-- newton applied to coninuous equations (not as accurate as newton applied to discrete equations, which I'm working on...)
	--solver = 'ConjRes',
	--solver = 'GMRes',
	--solver = 'JFNK',

	--useLineSearch = true,

	-- finite-difference order
	diffOrder = 2,
}
