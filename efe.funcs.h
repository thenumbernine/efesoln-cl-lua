// putting the type header in one place for ffi as well
// and the function / macro headers here:

extern constant real const c;
extern constant real const G;
extern constant int const sym3[3][3];
extern constant int const sym4[4][4];

<?
local d1coeffs = derivCoeffs[1][solver.diffOrder]
?>
extern constant real const d1coeffs[<?=#d1coeffs+1?>];
real d1coeff_for_offset(int offest);

<?
local d2coeffs = derivCoeffs[2][solver.diffOrder]
?>
extern constant real const d2coeffs[<?=#d2coeffs+1?>];

#define new_real_zero() 0.
#define real_zero new_real_zero()

#define real_add(a,b)	((a) + (b))
#define real_add2		real_add
#define real_sub(a,b)	((a) - (b))
#define real_mul(a,b)	((a) * (b))
#define real_real_mul	real_mul
#define real_dot		real_mul
#define real_div(a,b)	((a) / (b))
#define real_real_div	real_div

<?
local range = require 'ext.range'
local function makeAdds(ctype)
	local n = solver.diffOrder
	assert(n % 2 == 0)
	for i=4,n,2 do
		local xs = range(0,i-1):mapi(function(j) return "x"..j end)
?>
#define <?=ctype?>_add<?=i?>(<?=xs:concat', '?>)\
	<?=ctype?>_add2(\
		<?=ctype?>_add2(x0, x1),\
		<?=ctype?>_add<?=i-2?>(<?=xs:sub(3):concat', '?>)\
	)
<?	end
end
?>

<?makeAdds"real"?>

<?
local function makeOpsHeader(ctype, fieldtype, fields)
?>
#define new_<?=ctype?>_zero() ((<?=ctype?>){\
<? for _,field in ipairs(fields) do --\
?>	.<?=field?> = new_<?=fieldtype?>_zero(),\
<? end --\
?>})
extern constant <?=ctype?> const <?=ctype?>_zero;
<?
	for _,op in ipairs{"add", "sub"} do
?>
<?=ctype?> <?=ctype?>_<?=op?>(<?=ctype?> const a, <?=ctype?> const b);
<?
	end
	for _,op in ipairs{"mul", "div"} do
?>
<?=ctype?> <?=ctype?>_real_<?=op?>(<?=ctype?> const a, real const b);
<?
	end
?>
#define <?=ctype?>_add2 <?=ctype?>_add
<?=ctype?> <?=ctype?>_mul_add(<?=ctype?> const a, <?=ctype?> const b, real const c);

real <?=ctype?>_dot(<?=ctype?> const a, <?=ctype?> const b);
real <?=ctype?>_lenSq(<?=ctype?> const a);
real <?=ctype?>_len(<?=ctype?> const a);

// "norm" name for vectors and tensors
#define <?=ctype?>_normSq <?=ctype?>_lenSq
#define <?=ctype?>_norm <?=ctype?>_len

<?
	makeAdds(ctype)
end
?>

// I reported this bug to intel like 3 years ago.  it's fixed, right?
//is buggy with doubles on intel opencl ubuntu compiler
//#define _real3(a,b,c)          ((real3){.s={a,b,c}})
//so we do this instead and are safe:
#define _real3(a,b,c)            ((real3){.x=a, .y=b, .z=c})

<?makeOpsHeader("real3", "real", {"x", "y", "z"})?>
real3 real3_cross(real3 const a, real3 const b);

real3 real4_to_real3(real4 const a);
real4 real3_to_real4(real3 const a);
<?makeOpsHeader("real4", "real", {"s0", "s1", "s2", "s3"})?>

#define real3s3_ident ((real3s3){\
	.s00 = 1, .s01 = 0, .s02 = 0,\
	.s11 = 1, .s12 = 0,\
	.s22 = 1,\
})

real real3s3_det(real3s3 const m);
<?makeOpsHeader("real3s3", "real", {"s00", "s01", "s02", "s11", "s12", "s22"})?>
real3s3 real3s3_inv(real const d, real3s3 const m);
real3 real3s3_real3_mul(real3s3 const m, real3 const v);

<?makeOpsHeader("real4x4", "real4", {"s0", "s1", "s2", "s3"})?>
<?makeOpsHeader("real4s4", "real", {"s00", "s01", "s02", "s03", "s11", "s12", "s13", "s22", "s23", "s33"})?>
real4s4 real4s4_outer(real4 const v);
real4 real4s4_real4_mul(real4s4 const m, real4 const v);
real4x4 real4s4_real4s4_mul(real4s4 const a, real4s4 const b);
real3 real4s4_i0(real4s4 const a);
real3s3 real4s4_ij(real4s4 const a);
real real4s4_det(real4s4 const m);
real4s4 new_real4s4_Minkowski();
extern constant real4s4 const real4s4_Minkowski;


<?makeOpsHeader("real4x4x4", "real4x4", {"s0", "s1", "s2", "s3"})?>

<?makeOpsHeader("real4x4s4", "real4s4", {"s0", "s1", "s2", "s3"})?>

<?makeOpsHeader("real4s4x4s4", "real4s4", {"s00", "s01", "s02", "s03", "s11", "s12", "s13", "s22", "s23", "s33"})?>

<?makeOpsHeader("real4x4x4x4", "real4x4x4", {"s0", "s1", "s2", "s3"})?>

<?makeOpsHeader("real4x4x4s4", "real4x4s4", {"s0", "s1", "s2", "s3"})?>
