return {
	--size = 64,
	size = 8,
	
	-- line trace amount
	updateAlpha = 1,

	body = 'vacuum',
	--body = 'Earth',
	--body = 'Sun',
	--body = 'EM ring',
	--body = 'EM line',
	
	bodyRadii = 2,
	
	initCond = 'flat',
	--initCond = 'stellar Schwarzschild',
	--initCond = 'stellar Kerr-Newman',
	
	outputFilename = 'out.txt',

	--solver = 'Newton',	-- newton applied to coninuous equations (not as accurate as newton applied to discrete equations, which I'm working on...)
	--solver = 'ConjRes',
	--solver = 'GMRes',
	solver = 'JFNK',
}
