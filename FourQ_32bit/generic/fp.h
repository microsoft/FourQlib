/***********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: modular arithmetic and other low-level operations for 32-bit platforms
************************************************************************************/

#ifndef __FP_H__
#define __FP_H__


// For C++
#ifdef __cplusplus
extern "C" {
#endif


#include "../table_lookup.h"
#include "../FourQ_params.h"
#if (TARGET == TARGET_x86) && (COMPILER == COMPILER_VC)
    #include "intrin.h"
#endif

#define mask_26        (((uint32_t)1 << 26) - 1)
#define mask_23        (((uint32_t)1 << 23) - 1)


void digit_x_digit(digit_t a, digit_t b, digit_t* c)
{ // Digit multiplication, digit * digit -> 2-digit result    
    register digit_t al, ah, bl, bh, temp;
    digit_t albl, albh, ahbl, ahbh, res1, res2, res3, carry;
    digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t)*4), mask_high = (digit_t)(-1) << (sizeof(digit_t)*4);

    al = a & mask_low;                        // Low part
    ah = a >> (sizeof(digit_t) * 4);          // High part
    bl = b & mask_low;
    bh = b >> (sizeof(digit_t) * 4);

    albl = al*bl;
    albh = al*bh;
    ahbl = ah*bl;
    ahbh = ah*bh;
    c[0] = albl & mask_low;                   // C00

    res1 = albl >> (sizeof(digit_t) * 4);
    res2 = ahbl & mask_low;
    res3 = albh & mask_low;  
    temp = res1 + res2 + res3;
    carry = temp >> (sizeof(digit_t) * 4);
    c[0] ^= temp << (sizeof(digit_t) * 4);    // C01   

    res1 = ahbl >> (sizeof(digit_t) * 4);
    res2 = albh >> (sizeof(digit_t) * 4);
    res3 = ahbh & mask_low;
    temp = res1 + res2 + res3 + carry;
    c[1] = temp & mask_low;                   // C10 
    carry = temp & mask_high; 
    c[1] ^= (ahbh & mask_high) + carry;       // C11
}


static __inline void fpcopy1271(felm_t a, felm_t c)
{ // Copy of a field element, c = a
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++)
        c[i] = a[i];
}


void vadd1271(velm_t a, velm_t b, velm_t c) 
{ // Field addition over GF(2^127-1)
  // Redundant representation: 23/26/26/26/26-bit
    
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
    c[3] = a[3] + b[3];
    c[4] = a[4] + b[4];
}


void vsub1271(velm_t a, velm_t b, velm_t c) 
{ // Field subtraction over GF(2^127-1)
  // Redundant representation: 23/26/26/26/26-bit
    
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
    c[3] = a[3] - b[3];
    c[4] = a[4] - b[4];
}


void vneg1271(velm_t a) 
{ // Field negation over GF(2^127-1)
  // Redundant representation: 23/26/26/26/26-bit
    
    a[0] = mask_26 - a[0];
    a[1] = mask_26 - a[1];
    a[2] = mask_26 - a[2];
    a[3] = mask_26 - a[3];
    a[4] = mask_23 - a[4];
}


void vmul1271(velm_t a, velm_t b, velm_t c)
{ // Field multiplication, c = a*b mod p
   int32_t r0, r1, a0, a1, a2, a3, a4, b0, b1, b2, b3, b4; 
   int64_t c0, c1, c2, c3, c4;

#if (TARGET == TARGET_x86) && (COMPILER == COMPILER_VC)
   a0 = a[0]; a1 = a[1]; a2 = a[2]; a3 = a[3]; a4 = a[4]; 
   b0 = b[0]; b1 = b[1]; b2 = b[2]; b3 = b[3]; b4 = b[4];

   c0 = __emul((int)a0, (int)b0) + (__emul((int)a1, (int)b4) << 3) + (__emul((int)a4, (int)b1) << 3) + (__emul((int)a2, (int)b3) << 3) + (__emul((int)a3, (int)b2) << 3);
   c1 = __emul((int)a0, (int)b1) + __emul((int)a1, (int)b0) + (__emul((int)a2, (int)b4) << 3) + (__emul((int)a4, (int)b2) << 3) + (__emul((int)a3, (int)b3) << 3);
   c2 = __emul((int)a0, (int)b2) + __emul((int)a2, (int)b0) + __emul((int)a1, (int)b1) + (__emul((int)a3, (int)b4) << 3) + (__emul((int)a4, (int)b3) << 3);
   c3 = __emul((int)a0, (int)b3) + __emul((int)a3, (int)b0) + __emul((int)a1, (int)b2) + __emul((int)a2, (int)b1) + (__emul((int)a4, (int)b4) << 3);
   c4 = __emul((int)a0, (int)b4) + __emul((int)a4, (int)b0) + __emul((int)a1, (int)b3) + __emul((int)a3, (int)b1) + __emul((int)a2, (int)b2);
#else
   int64_t t1, t2, t3, t4;

   a0 = a[0]; a1 = a[1]; a2 = a[2]; a3 = a[3]; a4 = a[4]; 
   b0 = b[0]; b1 = b[1]; b2 = b[2]; b3 = b[3]; b4 = b[4];

   t1 = (int64_t)a1 << 3;
   t2 = (int64_t)a2 << 3;
   t3 = (int64_t)a3 << 3;
   t4 = (int64_t)a4 << 3;

   c0 = (int64_t)a0*b0 + (int64_t)t1*b4 + (int64_t)t4*b1 + (int64_t)t2*b3 + (int64_t)t3*b2;
   c1 = (int64_t)a0*b1 + (int64_t)a1*b0 + (int64_t)t2*b4 + (int64_t)t4*b2 + (int64_t)t3*b3;
   c2 = (int64_t)a0*b2 + (int64_t)a2*b0 + (int64_t)a1*b1 + (int64_t)t3*b4 + (int64_t)t4*b3;
   c3 = (int64_t)a0*b3 + (int64_t)a3*b0 + (int64_t)a1*b2 + (int64_t)a2*b1 + (int64_t)t4*b4;
   c4 = (int64_t)a0*b4 + (int64_t)a4*b0 + (int64_t)a1*b3 + (int64_t)a3*b1 + (int64_t)a2*b2;
#endif
   
                    r0   = c0 & mask_26; 
    c1 += c0 >> 26; r1   = c1 & mask_26;
    c2 += c1 >> 26; c[2] = c2 & mask_26; 
    c3 += c2 >> 26; c[3] = c3 & mask_26; 
    c4 += c3 >> 26; c[4] = c4 & mask_23; 
//    c4 += c3 >> 26; c[4] = c4 & mask_26; 
    
    c0   = r0 + (c4 >> 23);
//    c0   = r0 + ((c4 >> 26) << 3);                  
    c[0] = (int32_t)c0 & mask_26;
    c[1] = r1 + (int32_t)(c0 >> 26);
}


void vsqr1271(velm_t a, velm_t c)
{ // Field squaring, c = a*b mod p
   int32_t r0, r1, a0, a1, a2, a3, a4; 
   int64_t c0, c1, c2, c3, c4; 

#if (TARGET == TARGET_x86) && (COMPILER == COMPILER_VC)
   a0 = a[0]; a1 = a[1]; a2 = a[2]; a3 = a[3]; a4 = a[4];

   c0 = __emul((int)a0, (int)a0) + (__emul((int)a4, (int)a1) << 4) + (__emul((int)a2, (int)a3) << 4);
   c1 = (__emul((int)a0, (int)a1) << 1) + (__emul((int)a3, (int)a3) << 3) + (__emul((int)a4, (int)a2) << 4);
   c2 = (__emul((int)a0, (int)a2) << 1) + __emul((int)a1, (int)a1) + (__emul((int)a4, (int)a3) << 4);
   c3 = (__emul((int)a0, (int)a3) << 1) + (__emul((int)a1, (int)a2) << 1) + (__emul((int)a4, (int)a4) << 3);
   c4 = (__emul((int)a0, (int)a4) << 1) + (__emul((int)a1, (int)a3) << 1) + __emul((int)a2, (int)a2);

#else
   int64_t t0, t1, t2, t3, t4;

   a0 = a[0]; a1 = a[1]; a2 = a[2]; a3 = a[3]; a4 = a[4];

   t0 = (int64_t)a0 << 1;
   t1 = (int64_t)a1 << 1;
   t2 = (int64_t)a2 << 4;
   t3 = (int64_t)a3 << 3;
   t4 = (int64_t)a4 << 4;

   c0 = (int64_t)a0*a0 + (int64_t)t4*a1 + (int64_t)t2*a3;
   c1 = (int64_t)t0*a1 + (int64_t)t3*a3 + (int64_t)t4*a2;
   c2 = (int64_t)t0*a2 + (int64_t)a1*a1 + (int64_t)t4*a3;
   c3 = (int64_t)t0*a3 + (int64_t)t1*a2 + ((int64_t)a4 << 3)*a4;
   c4 = (int64_t)t0*a4 + (int64_t)t1*a3 + (int64_t)a2*a2;
#endif
   
                    r0   = c0 & mask_26; 
    c1 += c0 >> 26; r1   = c1 & mask_26;
    c2 += c1 >> 26; c[2] = c2 & mask_26; 
    c3 += c2 >> 26; c[3] = c3 & mask_26; 
    c4 += c3 >> 26; c[4] = c4 & mask_23; 
//    c4 += c3 >> 26; c[4] = c4 & mask_26; 
    
    c0   = r0 + (c4 >> 23);
//    c0   = r0 + ((c4 >> 26) << 3);                  
    c[0] = (int32_t)c0 & mask_26;
    c[1] = r1 + (int32_t)(c0 >> 26);
}


void vmod1271_incomplete(velm_t a, velm_t c)
{ // Reduce field element modulo 2^127-1
  // Redundant representation: 23/26/26/26/26-bit
  // Output is in the range [0, 2^127-1]
    int32_t t0, t1, t2, t3, t4;
    uint32_t rem;

    t0 = a[0]; t1 = a[1]; t2 = a[2]; t3 = a[3]; t4 = a[4];   

    // Carry propagation
    t1 += (t0 >> 26); t0 &= mask_26; 
    t2 += (t1 >> 26); t1 &= mask_26; 
    t3 += (t2 >> 26); t2 &= mask_26; 
    t4 += (t3 >> 26); t3 &= mask_26;
    rem = (t4 >> 23); t4 &= mask_23;
    
    // Correction
    t0 += rem; 
    t1 += (t0 >> 26); t0 &= mask_26; 
    t2 += (t1 >> 26); t1 &= mask_26; 
    t3 += (t2 >> 26); t2 &= mask_26; 
    t4 += (t3 >> 26); t3 &= mask_26;
    rem = (t4 >> 23); t4 &= mask_23;
    t0 += rem; 

    c[0] = t0; c[1] = t1; c[2] = t2; c[3] = t3; c[4] = t4;
}


void vmod1271(velm_t a, velm_t c)
{ // Reduce field element modulo 2^127-1
  // Output is in the range [0, 2^127-2]
  // Redundant representation: 23/26/26/26/26-bit
    int32_t t0, t1, t2, t3, t4;
    uint32_t mask, rem;

    t0 = a[0]; t1 = a[1]; t2 = a[2]; t3 = a[3]; t4 = a[4];   

    // First carry propagation
    t1 += (t0 >> 26); t0 &= mask_26; 
    t2 += (t1 >> 26); t1 &= mask_26; 
    t3 += (t2 >> 26); t2 &= mask_26; 
    t4 += (t3 >> 26); t3 &= mask_26;
    rem = (t4 >> 23); t4 &= mask_23;
    
    // First correction adding rem+1
    t0 += rem + 1; 
    t1 += (t0 >> 26); t0 &= mask_26; 
    t2 += (t1 >> 26); t1 &= mask_26; 
    t3 += (t2 >> 26); t2 &= mask_26; 
    t4 += (t3 >> 26); t3 &= mask_26;
    rem = (t4 >> 23); t4 &= mask_23;

    // If final carry = 0 then subtract 1
    mask = rem - 1;
    t0 -= (mask & 1);
    t1 += (t0 >> 26); t0 &= mask_26; 
    t2 += (t1 >> 26); t1 &= mask_26; 
    t3 += (t2 >> 26); t2 &= mask_26; 
    t4 += (t3 >> 26); t3 &= mask_26;
    
    c[0] = t0; c[1] = t1; c[2] = t2; c[3] = t3; c[4] = t4;
}


__inline void vexp1251(felm_t a, felm_t af)
{ // Exponentiation over GF(p), af = a^(125-1)
	int i;
	velm_t t1, t2, t3, t4, t5;

	vsqr1271(a, t2);
	vmul1271(a, t2, t2);
	vsqr1271(t2, t3);
	vsqr1271(t3, t3);
	vmul1271(t2, t3, t3);
	vsqr1271(t3, t4);
	vsqr1271(t4, t4);
	vsqr1271(t4, t4);
	vsqr1271(t4, t4);
	vmul1271(t3, t4, t4);
	vsqr1271(t4, t5);
	for (i = 0; i<7; i++) vsqr1271(t5, t5);
	vmul1271(t4, t5, t5);
	vsqr1271(t5, t2);
	for (i = 0; i<15; i++) vsqr1271(t2, t2);
	vmul1271(t5, t2, t2);
	vsqr1271(t2, t1);
	for (i = 0; i<31; i++) vsqr1271(t1, t1);
	vmul1271(t2, t1, t1);
	for (i = 0; i<32; i++) vsqr1271(t1, t1);
	vmul1271(t1, t2, t1);
	for (i = 0; i<16; i++) vsqr1271(t1, t1);
	vmul1271(t5, t1, t1);
	for (i = 0; i<8; i++) vsqr1271(t1, t1);
	vmul1271(t4, t1, t1);
	for (i = 0; i<4; i++) vsqr1271(t1, t1);
	vmul1271(t3, t1, t1);
	vsqr1271(t1, t1);
	vmul1271(a, t1, af);
}


void vinv1271(felm_t a)
{ // Field inversion, af = a^-1 = a^(p-2) mod p
  // Hardcoded for p = 2^127-1
	velm_t t;

	vexp1251(a, t);
	vsqr1271(t, t);
	vsqr1271(t, t);
	vmul1271(a, t, a);
}


void from_std_to_ext(f2elm_t a, v2elm_t c)
{ // Expand GF(p^2) element represented with two 4 32-bit digits to 23/26/26/26/26/23/26/26/26/26-bit vector representation
  // Assumes fully reduced input in [0, 2^127-1]  
    const uint32_t mask_8  = ((uint32_t)1 <<  8) - 1;   
    const uint32_t mask_14 = ((uint32_t)1 << 14) - 1;  
    const uint32_t mask_20 = ((uint32_t)1 << 20) - 1; 

    c[0] = a[0][0] & mask_26;
    c[1] = (a[0][0] >> 26) | ((a[0][1] & mask_20) <<  6);
    c[2] = (a[0][1] >> 20) | ((a[0][2] & mask_14) << 12);
    c[3] = (a[0][2] >> 14) | ((a[0][3] & mask_8 ) << 18);
    c[4] = (a[0][3] >>  8) & mask_23;

    c[5] = a[1][0] & mask_26;
    c[6] = (a[1][0] >> 26) | ((a[1][1] & mask_20) <<  6);
    c[7] = (a[1][1] >> 20) | ((a[1][2] & mask_14) << 12);
    c[8] = (a[1][2] >> 14) | ((a[1][3] & mask_8 ) << 18);
    c[9] = (a[1][3] >>  8) & mask_23;
}


void from_ext_to_std(v2elm_t a, f2elm_t c)
{ // Contract GF(p^2) element in 23/26/26/26/26/23/26/26/26/26-bit vector representation to two 4 32-bit digits
  // Assumes fully reduced input in [0, 2^127-1]
        
    c[0][0]  = (a[1] << 26) |  a[0];
    c[0][1]  = (a[2] << 20) | (a[1] >>  6);
    c[0][2]  = (a[3] << 14) | (a[2] >> 12);
    c[0][3]  = (a[4] <<  8) | (a[3] >> 18);
        
    c[1][0]  = (a[6] << 26) |  a[5];
    c[1][1]  = (a[7] << 20) | (a[6] >>  6);
    c[1][2]  = (a[8] << 14) | (a[7] >> 12);
    c[1][3]  = (a[9] <<  8) | (a[8] >> 18);
}


void mp_mul(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{ // Schoolbook multiprecision multiply, c = a*b   
    unsigned int i, j;
    digit_t u, v, UV[2];
    unsigned int carry = 0;

    for (i = 0; i < (2*nwords); i++) c[i] = 0;

    for (i = 0; i < nwords; i++) {
        u = 0;
        for (j = 0; j < nwords; j++) {
            MUL(a[i], b[j], UV+1, UV[0]); 
            ADDC(0, UV[0], u, carry, v); 
            u = UV[1] + carry;
            ADDC(0, c[i+j], v, carry, v); 
            u = u + carry;
            c[i+j] = v;
        }
        c[nwords+i] = u;
    }
}


unsigned int mp_add(digit_t* a, digit_t* b, digit_t* c, unsigned int nwords)
{ // Multiprecision addition, c = a+b, where lng(a) = lng(b) = nwords. Returns the carry bit 
    unsigned int i, carry = 0;

    for (i = 0; i < nwords; i++) {
        ADDC(carry, a[i], b[i], carry, c[i]);
    }
    
    return carry;
}


static __inline void multiply(const digit_t* a, const digit_t* b, digit_t* c)
{ // Schoolbook multiprecision multiply, c = a*b 

    mp_mul(a, b, c, NWORDS_ORDER);
}


static __inline unsigned int add(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{ // Multiprecision addition, c = a+b, where lng(a) = lng(b) = nwords. Returns the carry bit 
    unsigned int i, carry = 0;

    for (i = 0; i < nwords; i++) {
        ADDC(carry, a[i], b[i], carry, c[i]);
    }
    
    return carry;
}


unsigned int subtract(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{ // Multiprecision subtraction, c = a-b, where lng(a) = lng(b) = nwords. Returns the borrow bit
    unsigned int i, borrow = 0;

    for (i = 0; i < nwords; i++) {
        SUBC(borrow, a[i], b[i], borrow, c[i]);
    }

    return borrow;
}


void subtract_mod_order(const digit_t* a, const digit_t* b, digit_t* c)
{ // Subtraction modulo the curve order, c = a-b mod order
	digit_t mask, carry = 0;
	digit_t* order = (digit_t*)curve_order;
	unsigned int i, bout;

	bout = subtract(a, b, c, NWORDS_ORDER);            // (bout, c) = a - b
	mask = 0 - (digit_t)bout;                          // if bout = 0 then mask = 0x00..0, else if bout = 1 then mask = 0xFF..F

	for (i = 0; i < NWORDS_ORDER; i++) {               // c = c + (mask & order)
		ADDC(carry, c[i], mask & order[i], carry, c[i]);
	}
}


void add_mod_order(const digit_t* a, const digit_t* b, digit_t* c)
{ // Addition modulo the curve order, c = a+b mod order

	add(a, b, c, NWORDS_ORDER);                        // c = a + b
	subtract_mod_order(c, (digit_t*)&curve_order, c);  // if c >= order then c = c - order
}
 

void Montgomery_multiply_mod_order(const digit_t* ma, const digit_t* mb, digit_t* mc)
{ // 256-bit Montgomery multiplication modulo the curve order, mc = ma*mb*r' mod order, where ma,mb,mc in [0, order-1]
  // ma, mb and mc are assumed to be in Montgomery representation
  // The Montgomery constant r' = -r^(-1) mod 2^(log_2(r)) is the global value "Montgomery_rprime", where r is the order   
    unsigned int i;
    digit_t mask, P[2*NWORDS_ORDER], Q[2*NWORDS_ORDER], temp[2*NWORDS_ORDER];
	digit_t* order = (digit_t*)curve_order;
    unsigned int cout = 0, bout = 0;           

    multiply(ma, mb, P);                               // P = ma * mb
    multiply(P, (digit_t*)&Montgomery_rprime, Q);      // Q = P * r' mod 2^(log_2(r))
    multiply(Q, (digit_t*)&curve_order, temp);         // temp = Q * r
    cout = add(P, temp, temp, 2*NWORDS_ORDER);         // (cout, temp) = P + Q * r     

    for (i = 0; i < NWORDS_ORDER; i++) {               // (cout, mc) = (P + Q * r)/2^(log_2(r))
        mc[i] = temp[NWORDS_ORDER + i];
    }

    // Final, constant-time subtraction     
    bout = subtract(mc, (digit_t*)&curve_order, mc, NWORDS_ORDER);    // (cout, mc) = (cout, mc) - r
    mask = (digit_t)cout - (digit_t)bout;              // if (cout, mc) >= 0 then mask = 0x00..0, else if (cout, mc) < 0 then mask = 0xFF..F

    for (i = 0; i < NWORDS_ORDER; i++) {               // temp = mask & r
        temp[i] = (order[i] & mask);
    }
    add(mc, temp, mc, NWORDS_ORDER);                   //  mc = mc + (mask & r)

    return;
}


void modulo_order(digit_t* a, digit_t* c)
{ // Reduction modulo the order using Montgomery arithmetic
  // ma = a*Montgomery_Rprime mod r, where a,ma in [0, r-1], a,ma,r < 2^256
  // c = ma*1*Montgomery_Rprime^(-1) mod r, where ma,c in [0, r-1], ma,c,r < 2^256
	digit_t ma[NWORDS_ORDER], one[NWORDS_ORDER] = {0};
    
    one[0] = 1;
    Montgomery_multiply_mod_order(a, (digit_t*)&Montgomery_Rprime, ma);
    Montgomery_multiply_mod_order(ma, one, c);
}


void conversion_to_odd(digit_t* k, digit_t* k_odd)
{ // Convert scalar to odd if even using the prime subgroup order r
    digit_t mask;
    digit_t* order = (digit_t*)curve_order;
    unsigned int i, carry = 0;

    mask = ~(0 - (k[0] & 1));     

    for (i = 0; i < NWORDS_ORDER; i++) {   // If (k is odd) then k_odd = k else k_odd = k + r
        ADDC(carry, order[i] & mask, k[i], carry, k_odd[i]);
    }
}


__inline void vdiv1271(uint32_t* a) 
{ // GF(p) division by two, c = a/2 mod p
  // Redundant representation: 23/26/26/26/26-bit
    digit_t mask;

    mask = (0 - (a[0] & 1)) >> 6;  // if a[0] is odd then mask = 2^26-1, else mask = 0
    
    a[0] += mask;
    a[1] += mask;
    a[2] += mask;
    a[3] += mask;
    a[4] += (mask >> 3);

    a[0] = ((sdigit_t)a[0] >> 1) + ((a[1] & 1) << 25);
    a[1] = ((sdigit_t)a[1] >> 1) + ((a[2] & 1) << 25);
    a[2] = ((sdigit_t)a[2] >> 1) + ((a[3] & 1) << 25);
    a[3] = ((sdigit_t)a[3] >> 1) + ((a[4] & 1) << 25);
    a[4] = ((sdigit_t)a[4] >> 1);
}


#ifdef __cplusplus
}
#endif


#endif
