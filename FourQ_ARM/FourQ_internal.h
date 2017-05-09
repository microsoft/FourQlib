/***********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: internal header file
*
* This code is based on the paper "FourQ: four-dimensional decompositions on a 
* Q-curve over the Mersenne prime" by Craig Costello and Patrick Longa, in Advances 
* in Cryptology - ASIACRYPT, 2015.
* Preprint available at http://eprint.iacr.org/2015/565.
************************************************************************************/  

#ifndef __FOURQ_INTERNAL_H__
#define __FOURQ_INTERNAL_H__


// For C++
#ifdef __cplusplus
extern "C" {
#endif

    
#include "FourQ_api.h"


// Extended datatype definition
                     
typedef uint64_t uint128_t[2];


// Basic parameters for variable-base scalar multiplication (without using endomorphisms)
#define NPOINTS_VARBASE       (1 << (W_VARBASE-2)) 
#define t_VARBASE             ((NBITS_ORDER_PLUS_ONE+W_VARBASE-2)/(W_VARBASE-1))


// Basic parameters for fixed-base scalar multiplication
#define E_FIXEDBASE       (NBITS_ORDER_PLUS_ONE + W_FIXEDBASE*V_FIXEDBASE - 1)/(W_FIXEDBASE*V_FIXEDBASE)
#define D_FIXEDBASE       E_FIXEDBASE*V_FIXEDBASE
#define L_FIXEDBASE       D_FIXEDBASE*W_FIXEDBASE  
#define NPOINTS_FIXEDBASE V_FIXEDBASE*(1 << (W_FIXEDBASE-1))  
#define VPOINTS_FIXEDBASE (1 << (W_FIXEDBASE-1)) 
#if (NBITS_ORDER_PLUS_ONE-L_FIXEDBASE == 0)  // This parameter selection is not supported  
    #error -- "Unsupported parameter selection for fixed-base scalar multiplication"
#endif 


// Basic parameters for double scalar multiplication
#define NPOINTS_DOUBLEMUL_WP   (1 << (WP_DOUBLEBASE-2)) 
#define NPOINTS_DOUBLEMUL_WQ   (1 << (WQ_DOUBLEBASE-2)) 
   

// FourQ's point representations        

typedef struct { f2elm_t x; f2elm_t y; f2elm_t z; f2elm_t ta; f2elm_t tb; } point_extproj;  // Point representation in extended coordinates.
typedef point_extproj point_extproj_t[1];                                                              
typedef struct { f2elm_t xy; f2elm_t yx; f2elm_t z2; f2elm_t t2; } point_extproj_precomp;   // Point representation in extended coordinates (for precomputed points).
typedef point_extproj_precomp point_extproj_precomp_t[1];  
typedef struct { f2elm_t xy; f2elm_t yx; f2elm_t t2; } point_precomp;                       // Point representation in extended affine coordinates (for precomputed points).
typedef point_precomp point_precomp_t[1];


/********************** Constant-time unsigned comparisons ***********************/

// The following functions return 1 (TRUE) if condition is true, 0 (FALSE) otherwise

static __inline unsigned int is_digit_nonzero_ct(digit_t x)
{ // Is x != 0?
    return (unsigned int)((x | (0-x)) >> (RADIX-1));
}

static __inline unsigned int is_digit_zero_ct(digit_t x)
{ // Is x = 0?
    return (unsigned int)(1 ^ is_digit_nonzero_ct(x));
}

static __inline unsigned int is_digit_lessthan_ct(digit_t x, digit_t y)
{ // Is x < y?
    return (unsigned int)((x ^ ((x ^ y) | ((x - y) ^ y))) >> (RADIX-1)); 
}


/********************** Macros for digit operations **********************/

// Digit multiplication
#define MUL(multiplier, multiplicand, hi, lo)                                                     \
    digit_x_digit((multiplier), (multiplicand), &(lo));
    
// Digit addition with carry
#define ADDC(carryIn, addend1, addend2, carryOut, sumOut)                                         \
    { digit_t tempReg = (addend1) + (digit_t)(carryIn);                                           \
    (sumOut) = (addend2) + tempReg;                                                               \
    (carryOut) = (is_digit_lessthan_ct(tempReg, (digit_t)(carryIn)) | is_digit_lessthan_ct((sumOut), tempReg)); }

// Digit subtraction with borrow
#define SUBC(borrowIn, minuend, subtrahend, borrowOut, differenceOut)                             \
    { digit_t tempReg = (minuend) - (subtrahend);                                                 \
    unsigned int borrowReg = (is_digit_lessthan_ct((minuend), (subtrahend)) | ((borrowIn) & is_digit_zero_ct(tempReg)));  \
    (differenceOut) = tempReg - (digit_t)(borrowIn);                                              \
    (borrowOut) = borrowReg; }
    
// Shift right with flexible datatype
#define SHIFTR(highIn, lowIn, shift, shiftOut, DigitSize)                                         \
    (shiftOut) = ((lowIn) >> (shift)) ^ ((highIn) << (DigitSize - (shift)));

// 64x64-bit multiplication
#define MUL128(multiplier, multiplicand, product)                                                 \
    mp_mul((digit_t*)&(multiplier), (digit_t*)&(multiplicand), (digit_t*)&(product), NWORDS_FIELD/2);

// 128-bit addition, inputs < 2^127
#define ADD128(addend1, addend2, addition)                                                        \
    mp_add((digit_t*)(addend1), (digit_t*)(addend2), (digit_t*)(addition), NWORDS_FIELD);

// 128-bit addition with output carry
#define ADC128(addend1, addend2, carry, addition)                                                 \
    (carry) = mp_add((digit_t*)(addend1), (digit_t*)(addend2), (digit_t*)(addition), NWORDS_FIELD);


/**************** Function prototypes ****************/

/************* Arithmetic functions modulo the curve order **************/

// Converting to Montgomery representation
void to_Montgomery(const digit_t* ma, digit_t* c);

// Converting from Montgomery to standard representation
void from_Montgomery(const digit_t* a, digit_t* mc);

// 256-bit Montgomery multiplication modulo the curve order
void Montgomery_multiply_mod_order(const digit_t* ma, const digit_t* mb, digit_t* mc);

// Addition modulo the curve order, c = a+b mod order
void add_mod_order(const digit_t* a, const digit_t* b, digit_t* c);

// Subtraction modulo the curve order, c = a-b mod order
void subtract_mod_order(const digit_t* a, const digit_t* b, digit_t* c);

// Reduction modulo the order using Montgomery arithmetic internally
void modulo_order(digit_t* a, digit_t* c);

/************* Multiprecision functions **************/

// Check if multiprecision element is zero
bool is_zero_ct(digit_t* a, unsigned int nwords);

// Multiprecision addition, c = a+b. Returns the carry bit
unsigned int mp_add(digit_t* a, digit_t* b, digit_t* c, unsigned int nwords);

// Schoolbook multiprecision multiply, c = a*b
void mp_mul(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords);

// Multiprecision subtraction, c = a-b. Returns the borrow bit
unsigned int subtract(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords);

// Clear "nwords" integer-size digits from memory
extern void clear_words(void* mem, unsigned int nwords);

/************ Field arithmetic functions *************/

// Copy of a field element, c = a
void fpcopy1271(felm_t a, felm_t c);

// Field negation, a = -a mod p
void fpneg1271(felm_t a);

// Integer multiplication, c = a*b
void mul1271_a(felm_t a, felm_t b, felm_t c);

// Modular correction, a = a mod p
void mod1271(felm_t a);
void mod1271_a(felm_t a, felm_t c);

// Field addition, c = a+b mod p
void fpadd1271(felm_t a, felm_t b, felm_t c);
void fpadd1271_a(felm_t a, felm_t b, felm_t c);

// Field subtraction, c = a-b mod p
void fpsub1271(felm_t a, felm_t b, felm_t c);
void fpsub1271_a(felm_t a, felm_t b, felm_t c);

// Field division by two, c = a/2 mod p
void fpdiv1271(felm_t a);

// Field multiplication, c = a*b mod p
void fpmul1271(felm_t a, felm_t b, felm_t c);

// Field squaring, c = a^2 mod p
void fpsqr1271(felm_t a, felm_t c);

// Field inversion, af = a^-1 = a^(p-2) mod p
void fpinv1271(felm_t a);

// Exponentiation over GF(p), af = a^(125-1)
void fpexp1251(felm_t a, felm_t af);

/************ Quadratic extension field arithmetic functions *************/

// Zeroing a quadratic extension field element, a=0 
void fp2zero1271(f2elm_t a);

// Copy quadratic extension field element, c = a
void fp2copy1271(f2elm_t a, f2elm_t c);

// Quadratic extension field negation, a = -a in GF((2^127-1)^2)
void fp2neg1271(f2elm_t a);

// Quadratic extension field addition, c = a+b in GF((2^127-1)^2)
void fp2add1271(f2elm_t a, f2elm_t b, f2elm_t c);
void fp2add1271_a(f2elm_t a, f2elm_t b, f2elm_t c);

// Quadratic extension field subtraction, c = a-b in GF((2^127-1)^2)
void fp2sub1271(f2elm_t a, f2elm_t b, f2elm_t c);
void fp2sub1271_a(f2elm_t a, f2elm_t b, f2elm_t c);

// Quadratic extension field multiplication, c = a*b in GF((2^127-1)^2)
void fp2mul1271(f2elm_t a, f2elm_t b, f2elm_t c);
void fp2mul1271_a(f2elm_t a, f2elm_t b, f2elm_t c);

// Quadratic extension field squaring, c = a^2 in GF((2^127-1)^2)
void fp2sqr1271(f2elm_t a, f2elm_t c);
void fp2sqr1271_a(f2elm_t a, f2elm_t c);

// Quadratic extension field inversion, af = a^-1 = a^(p-2) in GF((2^127-1)^2)
void fp2inv1271(f2elm_t a);

/************ Curve and recoding functions *************/

// Normalize projective twisted Edwards point Q = (X,Y,Z) -> P = (x,y)
void eccnorm(point_extproj_t P, point_t Q);

// Conversion from representation (X,Y,Z,Ta,Tb) to (X+Y,Y-X,2Z,2dT), where T = Ta*Tb
void R1_to_R2(point_extproj_t P, point_extproj_precomp_t Q);

// Conversion from representation (X,Y,Z,Ta,Tb) to (X+Y,Y-X,Z,T), where T = Ta*Tb 
void R1_to_R3(point_extproj_t P, point_extproj_precomp_t Q);

// Conversion from representation (X+Y,Y-X,2Z,2dT) to (2X,2Y,2Z,2dT)
void R2_to_R4(point_extproj_precomp_t P, point_extproj_t Q);

// Point doubling 2P
void eccdouble(point_extproj_t P);

// Complete point addition P = P+Q or P = P+P
void eccadd(point_extproj_precomp_t Q, point_extproj_t P);
void eccadd_core(point_extproj_precomp_t P, point_extproj_precomp_t Q, point_extproj_t R);

// Psi mapping of a point, P = psi(P)
void ecc_psi(point_extproj_t P);

// Phi mapping of a point, P = phi(P)
void ecc_phi(point_extproj_t P);

// Scalar decomposition
void decompose(uint64_t* k, uint64_t* scalars);

// Recoding sub-scalars for use in the variable-base scalar multiplication
void recode(uint64_t* scalars, unsigned int* digits, unsigned int* sign_masks);

// Computes the fixed window representation of scalar
void fixed_window_recode(uint64_t* scalar, unsigned int* digits, unsigned int* sign_masks);

// Convert scalar to odd if even using the prime subgroup order r
void conversion_to_odd(digit_t* k, digit_t* k_odd);

// Co-factor clearing
void cofactor_clearing(point_extproj_t P);

// Reduction modulo the order using Montgomery arithmetic
void modulo_order(digit_t* a, digit_t* c);

// Precomputation function
void ecc_precomp(point_extproj_t P, point_extproj_precomp_t *T);

// Constant-time table lookup to extract an extended twisted Edwards point (X+Y:Y-X:2Z:2T) from the precomputed table
void table_lookup_1x8(point_extproj_precomp_t* table, point_extproj_precomp_t P, unsigned int digit, unsigned int sign_mask);
void table_lookup_1x8_a(point_extproj_precomp_t* table, point_extproj_precomp_t P, unsigned int* digit, unsigned int* sign_mask);

// Modular correction of input coordinates and conversion to representation (X,Y,Z,Ta,Tb) 
void point_setup(point_t P, point_extproj_t Q);

// Point validation: check if point lies on the curve     
bool ecc_point_validate(point_extproj_t P);

// Output error/success message for a given ECCRYPTO_STATUS
const char* FourQ_get_error_message(ECCRYPTO_STATUS Status);

// Mixed point addition P = P+Q or P = P+P
void eccmadd_ni(point_precomp_t Q, point_extproj_t P);

// Constant-time table lookup to extract a point represented as (x+y,y-x,2t)
void table_lookup_fixed_base(point_precomp_t* table, point_precomp_t P, unsigned int digit, unsigned int sign);

//  Computes the modified LSB-set representation of scalar
void mLSB_set_recode(uint64_t* scalar, unsigned int *digits);

// Generation of the precomputation table used internally by the double scalar multiplication function ecc_mul_double()
void ecc_precomp_double(point_extproj_t P, point_extproj_precomp_t* Table, unsigned int npoints);

// Computes wNAF recoding of a scalar
void wNAF_recode(uint64_t scalar, unsigned int w, int* digits);

// Encode point P
void encode(point_t P, unsigned char* Pencoded);

// Decode point P
ECCRYPTO_STATUS decode(const unsigned char* Pencoded, point_t P);


/************ Functions based on macros *************/

// Copy extended projective point Q = (X:Y:Z:Ta:Tb) to P
#define ecccopy(Q, P); fp2copy1271((Q)->x,  (P)->x);  \
                       fp2copy1271((Q)->y,  (P)->y);  \
                       fp2copy1271((Q)->z,  (P)->z);  \
                       fp2copy1271((Q)->ta, (P)->ta); \
                       fp2copy1271((Q)->tb, (P)->tb);

// Copy extended projective point Q = (X+Y,Y-X,2Z,2dT) to P
#define ecccopy_precomp(Q, P); fp2copy1271((Q)->xy, (P)->xy); \
                               fp2copy1271((Q)->yx, (P)->yx); \
                               fp2copy1271((Q)->z2, (P)->z2); \
                               fp2copy1271((Q)->t2, (P)->t2); 

// Copy extended affine point Q = (x+y,y-x,2dt) to P
#define ecccopy_precomp_fixed_base(Q, P); fp2copy1271((Q)->xy, (P)->xy); \
                                          fp2copy1271((Q)->yx, (P)->yx); \
                                          fp2copy1271((Q)->t2, (P)->t2);


#ifdef __cplusplus
}
#endif


#endif
