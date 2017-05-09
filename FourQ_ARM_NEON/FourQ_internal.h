/*****************************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: internal header file
*
* This code is based on the papers:
* [1] "FourQ: four-dimensional decompositions on a Q-curve over the Mersenne prime" 
*     by Craig Costello and Patrick Longa, ASIACRYPT2015 (http://eprint.iacr.org/2015/565).
* [2] "FourQNEON: Faster Elliptic Curve Scalar Multiplications on ARM Processors" 
*     by Patrick Longa, SAC2016 (http://eprint.iacr.org/2016/645).
******************************************************************************************/  

#ifndef __FOURQ_INTERNAL_H__
#define __FOURQ_INTERNAL_H__


// For C++
#ifdef __cplusplus
extern "C" {
#endif

    
#include "FourQ_api.h"


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


// FourQ's vectorized point representations 

typedef struct { v2elm_t x; v2elm_t y; v2elm_t z; v2elm_t ta; v2elm_t tb; } vpoint_extproj;  // Point representation in extended coordinates.
typedef vpoint_extproj vpoint_extproj_t[1];                                                                 
typedef struct { v2elm_t xy; v2elm_t yx; v2elm_t z2; v2elm_t t2; } vpoint_extproj_precomp;   // Point representation in extended coordinates (for precomputed points).
typedef vpoint_extproj_precomp vpoint_extproj_precomp_t[1];  
typedef struct { v2elm_t xy; v2elm_t yx; v2elm_t t2; } vpoint_precomp;                       // Point representation in extended affine coordinates (for precomputed points).
typedef vpoint_precomp vpoint_precomp_t[1];


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

// Multiprecision subtraction, c = a-b. Returns the borrow bit
unsigned int subtract(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords);

// Clear "nwords" integer-size digits from memory
extern void clear_words(void* mem, unsigned int nwords);

/************ Field arithmetic functions *************/
   
// Field negation, a = -a mod p
void vneg1271(velm_t a);

// Field negation over GF(2^127-1) of second element of a GF(p^2) element
void v2neg1271_felm(v2elm_t a);

// Modular correction, a = a mod p
void vmod1271(velm_t a, velm_t c);

// Field addition, c = a+b mod p
void vadd1271(velm_t a, velm_t b, velm_t c);

// Field subtraction, c = a-b mod p
void vsub1271(velm_t a, velm_t b, velm_t c);

// Field division by two, c = a/2 mod p
void vdiv1271(uint32_t* a); 

// Field squaring, c = a^2 mod p
void vsqr1271(velm_t a, velm_t c);
void vsqr1271_a(uint32_t* a, uint32_t* c);

// Field multiplication, c = a*b mod p
void vmul1271(velm_t a, velm_t b, velm_t c);
void vmul1271_a(uint32_t* a, uint32_t* b, uint32_t* c);

// Field inversion, af = a^-1 = a^(p-2) mod p
void vinv1271(velm_t a);

// Exponentiation over GF(p), af = a^(125-1)
void vexp1251(velm_t a, velm_t af);

// Conversion functions
void from_std_to_ext(f2elm_t a, v2elm_t c);
void from_ext_to_std(v2elm_t a, f2elm_t c);

/************ Quadratic extension field arithmetic functions *************/

// Convert GF(p^2) element in 23/23/26/26/26/26/26/26/26/26-bit vector representation to 
// two field elements in 23/26/26/26/26-bit vector representation
void from_v2_to_v(v2elm_t a, velm_t c0, velm_t c1);

// Convert two field elements in 23/26/26/26/26-bit vector representation to GF(p^2) element in  
// 23/23/26/26/26/26/26/26/26/26-bit vector representation
void from_v_to_v2(velm_t a0, velm_t a1, v2elm_t c);

// Zeroing a quadratic extension field element, a=0 
void v2zero1271(v2elm_t a);

// Copy quadratic extension field element, c = a
void fp2copy1271(f2elm_t a, f2elm_t c); 
void v2copy1271(v2elm_t a, v2elm_t c); 

// Quadratic extension field negation, a = -a in GF((2^127-1)^2)
void v2neg1271(v2elm_t a);

// Quadratic extension field addition, c = a+b in GF((2^127-1)^2)
void v2add1271(v2elm_t a, v2elm_t b, v2elm_t c);
void v2add1271_a(uint32_t* a, uint32_t* b, uint32_t* c);

// Quadratic extension field subtraction, c = a-b in GF((2^127-1)^2)
void v2sub1271(v2elm_t a, v2elm_t b, v2elm_t c);
void v2sub1271_a(uint32_t* a, uint32_t* b, uint32_t* c);

// Quadratic extension field addition followed by subtraction over GF(2^127-1)
void v2dblsub1271(v2elm_t a, v2elm_t b, v2elm_t c);
void v2dblsub1271_a(uint32_t* a, uint32_t* b, uint32_t* c); 

// Quadratic extension field addition and subtraction over GF(2^127-1)
void v2addsub1271(v2elm_t a, v2elm_t b, v2elm_t c, v2elm_t d);
void v2addsub1271_a(uint32_t* a, uint32_t* b, uint32_t* c, uint32_t* d); 

// Modular correction over GF(p^2)
void v2mod1271(v2elm_t a, v2elm_t c);

// Quadratic extension field multiplication, c = a*b in GF((2^127-1)^2)
void v2mul1271(v2elm_t a, v2elm_t b, v2elm_t c);
void v2mul1271_a(uint32_t* a, uint32_t* b, uint32_t* c);

// Quadratic extension field squaring, c = a^2 in GF((2^127-1)^2)
void v2sqr1271(v2elm_t a, v2elm_t c);
void v2sqr1271_a(uint32_t* a, uint32_t* c);

// Vectorized GF(p^2) squaring/addition in GF((2^127-1)^2) 
void v2sqradd1271(v2elm_t a, v2elm_t c, v2elm_t d, v2elm_t e, v2elm_t f);   
void v2sqradd1271_a(uint32_t* a, uint32_t* c, uint32_t* d, uint32_t* e, uint32_t* f); 

// Vectorized GF(p^2) squaring/addition/subtraction in GF((2^127-1)^2) 
void v2sqraddsub1271(v2elm_t a, v2elm_t c, v2elm_t d, v2elm_t e, v2elm_t f, v2elm_t g);   
void v2sqraddsub1271_a(uint32_t* a, uint32_t* c, uint32_t* d, uint32_t* e, uint32_t* f, uint32_t* g); 

// Vectorized GF(p^2) multiplication/addition in GF((2^127-1)^2) 
void v2muladd1271(v2elm_t a, v2elm_t b, v2elm_t c, v2elm_t d, v2elm_t e, v2elm_t f);   
void v2muladd1271_a(uint32_t* a, uint32_t* b, uint32_t* c, uint32_t* d, uint32_t* e, uint32_t* f);    

// Vectorized GF(p^2) multiplication/subtraction in GF((2^127-1)^2) 
void v2mulsub1271(v2elm_t a, v2elm_t b, v2elm_t c, v2elm_t d, v2elm_t e, v2elm_t f);   
void v2mulsub1271_a(uint32_t* a, uint32_t* b, uint32_t* c, uint32_t* d, uint32_t* e, uint32_t* f);    

// Vectorized GF(p^2) multiplication/addition/subtraction in GF((2^127-1)^2) 
void v2muladdsub1271(v2elm_t a, v2elm_t b, v2elm_t c, v2elm_t d, v2elm_t e, v2elm_t f, v2elm_t g);   
void v2muladdsub1271_a(uint32_t* a, uint32_t* b, uint32_t* c, uint32_t* d, uint32_t* e, uint32_t* f, uint32_t* g);       

// Vectorized GF(p^2) multiplication/addition/subtraction in GF((2^127-1)^2) 
void v2muldlbsub1271(v2elm_t a, v2elm_t b, v2elm_t c, v2elm_t d, v2elm_t e, v2elm_t f);   
void v2muldblsub1271_a(uint32_t* a, uint32_t* b, uint32_t* c, uint32_t* d, uint32_t* e, uint32_t* f);           

// Vectorized GF(p^2) inversion, af = a^-1 = a^(p-2) in GF((2^127-1)^2)
void v2inv1271(v2elm_t a);

/************ Curve and recoding functions *************/

// Normalize projective twisted Edwards point Q = (X,Y,Z) -> P = (x,y)
void eccnorm(vpoint_extproj_t P, vpoint_t Q);

// Conversion from representation (X,Y,Z,Ta,Tb) to (X+Y,Y-X,2Z,2dT), where T = Ta*Tb
void R1_to_R2(vpoint_extproj_t P, vpoint_extproj_precomp_t Q);

// Conversion from representation (X,Y,Z,Ta,Tb) to (X+Y,Y-X,Z,T), where T = Ta*Tb 
void R1_to_R3(vpoint_extproj_t P, vpoint_extproj_precomp_t Q);  
 
// Conversion from representation (X+Y,Y-X,2Z,2dT) to (2X,2Y,2Z,2dT)
void R2_to_R4(vpoint_extproj_precomp_t P, vpoint_extproj_t Q);        

// Point doubling 2P
void eccdouble(vpoint_extproj_t P);

// Complete point addition P = P+Q or P = P+P
void eccadd(vpoint_extproj_precomp_t Q, vpoint_extproj_t P);
void eccadd_core(vpoint_extproj_precomp_t P, vpoint_extproj_precomp_t Q, vpoint_extproj_t R); 

// Psi mapping of a point, P = psi(P)
void ecc_psi(vpoint_extproj_t P); 

// Phi mapping of a point, P = phi(P)
void ecc_phi(vpoint_extproj_t P);

// Scalar decomposition
void decompose(uint64_t* k, uint64_t* scalars);

// 256-bit multiplication with truncation for the scalar decomposition
// Outputs 64-bit value c = (uint64_t)((a * b) >> 256)
void mul_truncate_a(uint32_t* a, uint32_t* b, uint32_t* c);

// Recoding sub-scalars for use in the variable-base scalar multiplication
void recode(uint64_t* scalars, unsigned int* digits, unsigned int* sign_masks);

// Convert scalar to odd if even using the prime subgroup order r
void conversion_to_odd(digit_t* k, digit_t* k_odd);

// Co-factor clearing
void cofactor_clearing(vpoint_extproj_t P);

// Reduction modulo the order using Montgomery arithmetic
void modulo_order(digit_t* a, digit_t* c);

// Precomputation function
void ecc_precomp(vpoint_extproj_t P, vpoint_extproj_precomp_t *T);

// Constant-time table lookup to extract an extended twisted Edwards point (X+Y:Y-X:2Z:2T) from the precomputed table
void table_lookup_1x8(vpoint_extproj_precomp_t* table, vpoint_extproj_precomp_t P, unsigned int digit, unsigned int sign_mask);

// Modular correction of input coordinates and conversion to representation (X,Y,Z,Ta,Tb) 
void point_setup(point_t P, vpoint_extproj_t Q);
    
// Point validation: check if point lies on the curve     
bool ecc_point_validate(vpoint_extproj_t P);

// Output error/success message for a given ECCRYPTO_STATUS
const char* FourQ_get_error_message(ECCRYPTO_STATUS Status);

// Constant-time table lookup to extract a point represented as (x+y,y-x,2t)
void table_lookup_fixed_base(vpoint_precomp_t* table, vpoint_precomp_t P, unsigned int digit, unsigned int sign);

//  Computes the modified LSB-set representation of scalar
void mLSB_set_recode(uint64_t* scalar, unsigned int *digits);

// Generation of the precomputation table used internally by the double scalar multiplication function ecc_mul_double()
void ecc_precomp_double(vpoint_extproj_t P, vpoint_extproj_precomp_t* Table, unsigned int npoints);

// Computes wNAF recoding of a scalar
void wNAF_recode(uint64_t scalar, unsigned int w, int* digits);

// Encode point P
void encode(point_t P, unsigned char* Pencoded);

// Decode point P
ECCRYPTO_STATUS decode(const unsigned char* Pencoded, point_t P);


/************ Functions based on macros *************/

// Copy extended projective point Q = (X:Y:Z:Ta:Tb) to P
#define ecccopy(Q, P); v2copy1271((Q)->x,  (P)->x);  \
                       v2copy1271((Q)->y,  (P)->y);  \
                       v2copy1271((Q)->z,  (P)->z);  \
                       v2copy1271((Q)->ta, (P)->ta); \
                       v2copy1271((Q)->tb, (P)->tb);

// Copy extended affine point Q = (x+y,y-x,2dt) to P
#define ecccopy_precomp_fixed_base(Q, P); v2copy1271((Q)->xy, (P)->xy); \
                                          v2copy1271((Q)->yx, (P)->yx); \
                                          v2copy1271((Q)->t2, (P)->t2);

// Vectorize extended projective point Q = (X:Y:Z:Ta:Tb)
#define point_from_std_to_ext(Q, P); from_std_to_ext(Q->x,  P->x);  \
                                     from_std_to_ext(Q->y,  P->y);  \
                                     from_std_to_ext(Q->z,  P->z);  \
                                     from_std_to_ext(Q->ta, P->ta); \
                                     from_std_to_ext(Q->tb, P->tb);

#ifdef __cplusplus
}
#endif


#endif
