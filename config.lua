return {
	size = 64,
	
	-- line trace amount
	updateAlpha = 1,
	
	body = 'Earth',
	--body = 'Sun',
	--body = 'EM Field',
	
	bodyRadii = 2,
	
	initCond = 'flat',
	--initCond = 'stellar Schwarzschild',
	--initCond = 'stellar Kerr-Newman',
	
	outputFilename = 'out.txt',

	solver = 'Newton',
	--solver = 'ConjRes',
}
