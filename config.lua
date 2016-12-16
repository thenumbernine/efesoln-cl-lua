return {
	size = 64,
	
	updateAlpha = 1e+30,
	
	body = 'Earth',
	--body = 'Sun',
	--body = 'EM Field',
	
	bodyRadii = 2,
	
	initCond = 'flat',
	--initCond = 'stellar Schwarzschild',
	--initCond = 'stellar Kerr-Newman',
	
	outputFilename = 'out.txt',

	--solver = 'Newton',
	solver = 'ConjRes',
}
