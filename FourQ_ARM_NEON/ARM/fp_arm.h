/***********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: modular arithmetic and other low-level operations for ARM processors
************************************************************************************/

#ifndef __FP_ARM_H__
#define __FP_ARM_H__


// For C++
#ifdef __cplusplus
extern "C" {
#endif


#include "../FourQ_params.h"

#define mask_26        (((uint32_t)1 << 26) - 1)
#define mask_23        (((uint32_t)1 << 23) - 1)


static __inline void fpcopy1271(felm_t a, felm_t c)
{ // Copy of a field element, c = a
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++)
        c[i] = a[i];
}


void vadd1271(velm_t a, velm_t b, velm_t c) 
{ // Vectorized field addition over GF(2^127-1)
  // Representation: 23/26/26/26/26-bit
    
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
    c[3] = a[3] + b[3];
    c[4] = a[4] + b[4];
}


void vsub1271(velm_t a, velm_t b, velm_t c) 
{ // Vectorized field subtraction over GF(2^127-1)
  // Representation: 23/26/26/26/26-bit
    
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
    c[3] = a[3] - b[3];
    c[4] = a[4] - b[4];
}


void vneg1271(velm_t a) 
{ // Vectorized field negation over GF(2^127-1)
  // Representation: 23/26/26/26/26-bit
    
    a[0] = mask_26 - a[0];
    a[1] = mask_26 - a[1];
    a[2] = mask_26 - a[2];
    a[3] = mask_26 - a[3];
    a[4] = mask_23 - a[4];
}


void v2neg1271_felm(v2elm_t a) 
{ // Vectorized field negation over GF(2^127-1) of second element of a GF(p^2) element
  // Representation: 23/26/26/26/26-bit
    
    a[1] = mask_26 - a[1];
    a[3] = mask_26 - a[3];
    a[5] = mask_26 - a[5];
    a[7] = mask_26 - a[7];
    a[9] = mask_23 - a[9];
}


void vmul1271(velm_t a, velm_t b, velm_t c)
{ // Vectorized field multiplication, c = a*b mod p
    vmul1271_a(a, b, c);
}


void vsqr1271(velm_t a, velm_t c)
{ // Vectorized field squaring, c = a*b mod p
    vsqr1271_a(a, c);
}


void vmod1271_incomplete(velm_t a, velm_t c)
{ // Reduce vectorized field element modulo 2^127-1
  // Representation: 23/26/26/26/26-bit
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
{ // Reduce vectorized field element modulo 2^127-1
  // Output is in the range [0, 2^127-2]
  // Representation: 23/26/26/26/26-bit
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


void v2mod1271_incomplete(v2elm_t a, v2elm_t c)
{ // Reduce vectorized GF(p^2) element modulo 2^127-1
  // Representation: 23/23/26/26/26/26/26/26/26/26-bit
  // Output is in the range [0, 2^127-1]
    int32_t t0, t1, t2, t3, t4;
    uint32_t rem;

    t0 = a[0]; t1 = a[2]; t2 = a[4]; t3 = a[6]; t4 = a[8];   

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

    c[0] = t0; c[2] = t1; c[4] = t2; c[6] = t3; c[8] = t4;
    t0 = a[1]; t1 = a[3]; t2 = a[5]; t3 = a[7]; t4 = a[9];   

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

    c[1] = t0; c[3] = t1; c[5] = t2; c[7] = t3; c[9] = t4;
}


void v2mod1271(velm_t a, velm_t c)
{ // Reduce vectorized GF(p^2) element modulo 2^127-1
  // Output is in the range [0, 2^127-2]
  // Representation: 23/23/26/26/26/26/26/26/26/26-bit
    int32_t t0, t1, t2, t3, t4;
    uint32_t mask, rem;
    
    t0 = a[0]; t1 = a[2]; t2 = a[4]; t3 = a[6]; t4 = a[8];   

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
    
    c[0] = t0; c[2] = t1; c[4] = t2; c[6] = t3; c[8] = t4;
    t0 = a[1]; t1 = a[3]; t2 = a[5]; t3 = a[7]; t4 = a[9];    

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

    c[1] = t0; c[3] = t1; c[5] = t2; c[7] = t3; c[9] = t4;
}


__inline void vexp1251(felm_t a, felm_t af)
{ // Vectorized exponentiation over GF(p), af = a^(125-1)
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
{ // Vectorized field inversion, af = a^-1 = a^(p-2) mod p
  // Hardcoded for p = 2^127-1
	velm_t t;

	vexp1251(a, t);
	vsqr1271(t, t);
	vsqr1271(t, t);
	vmul1271(a, t, a);
}


void from_std_to_ext(f2elm_t a, v2elm_t c)
{ // Expand GF(p^2) element represented with two 4 32-bit digits to 23/23/26/26/26/26/26/26/26/26-bit vector representation
  // Assumes fully reduced input in [0, 2^127-1]  
    const uint32_t mask_8  = ((uint32_t)1 <<  8) - 1;   
    const uint32_t mask_14 = ((uint32_t)1 << 14) - 1;  
    const uint32_t mask_20 = ((uint32_t)1 << 20) - 1; 

    c[0] = a[0][0] & mask_26;
    c[2] = (a[0][0] >> 26) | ((a[0][1] & mask_20) <<  6);
    c[4] = (a[0][1] >> 20) | ((a[0][2] & mask_14) << 12);
    c[6] = (a[0][2] >> 14) | ((a[0][3] & mask_8 ) << 18);
    c[8] = (a[0][3] >>  8) & mask_23;

    c[1] = a[1][0] & mask_26;
    c[3] = (a[1][0] >> 26) | ((a[1][1] & mask_20) <<  6);
    c[5] = (a[1][1] >> 20) | ((a[1][2] & mask_14) << 12);
    c[7] = (a[1][2] >> 14) | ((a[1][3] & mask_8 ) << 18);
    c[9] = (a[1][3] >>  8) & mask_23;
}


void from_ext_to_std(v2elm_t a, f2elm_t c)
{ // Contract GF(p^2) element in 23/23/26/26/26/26/26/26/26/26-bit vector representation to two 4 32-bit digits
  // Assumes fully reduced input in [0, 2^127-1]
        
    c[0][0]  = (a[2] << 26) |  a[0];
    c[0][1]  = (a[4] << 20) | (a[2] >>  6);
    c[0][2]  = (a[6] << 14) | (a[4] >> 12);
    c[0][3]  = (a[8] <<  8) | (a[6] >> 18);
        
    c[1][0]  = (a[3] << 26) |  a[1];
    c[1][1]  = (a[5] << 20) | (a[3] >>  6);
    c[1][2]  = (a[7] << 14) | (a[5] >> 12);
    c[1][3]  = (a[9] <<  8) | (a[7] >> 18);
}


void from_v2_to_v(v2elm_t a, velm_t c0, velm_t c1)
{ // Convert GF(p^2) element in 23/23/26/26/26/26/26/26/26/26-bit vector representation to 
  // two field elements in 23/26/26/26/26-bit vector representation
        
    c0[0]  = a[0]; 
    c0[1]  = a[2]; 
    c0[2]  = a[4]; 
    c0[3]  = a[6]; 
    c0[4]  = a[8]; 

    c1[0]  = a[1];
    c1[1]  = a[3]; 
    c1[2]  = a[5]; 
    c1[3]  = a[7]; 
    c1[4]  = a[9]; 
}


void from_v_to_v2(velm_t a0, velm_t a1, v2elm_t c)
{ // Convert two field elements in 23/26/26/26/26-bit vector representation to GF(p^2) element in  
  // 23/23/26/26/26/26/26/26/26/26-bit vector representation
        
    c[0]  = a0[0]; 
    c[1]  = a1[0]; 
    c[2]  = a0[1]; 
    c[3]  = a1[1]; 
    c[4]  = a0[2]; 
    c[5]  = a1[2];
    c[6]  = a0[3]; 
    c[7]  = a1[3]; 
    c[8]  = a0[4]; 
    c[9]  = a1[4]; 
}


static __inline void MULADDADD(digit_t a, digit_t b, digit_t* u, digit_t* c)
{ // Multiplication-addition-addition
     
asm volatile(
    "push   {r4-r7}           \n\t"
    "umull  r4, r5, %0, %1    \n\t"
    "ldr    r6, [%2]          \n\t"     
    "ldr    r7, [%3]          \n\t" 
    "adds   r4, r6            \n\t"     
    "adcs   r5, #0            \n\t" 
    "adds   r4, r7            \n\t"     
    "adcs   r5, #0            \n\t" 
    "str    r4, [%3]          \n\t"     
    "str    r5, [%2]          \n\t"
    "pop    {r4-r7}           \n\t"
    :
    :"r"(a), "r"(b), "r"(&u[0]), "r"(&c[0])
    :"memory","r4","r5","r6","r7"
	);    
    return;
}


static __inline void multiply(const digit_t* a, const digit_t* b, digit_t* c)
{ // Schoolbook multiprecision multiply, c = a*b  
    unsigned int i, j;
    digit_t u, v, UV[2];
    unsigned int carry = 0;

     for (i = 0; i < (2*NWORDS_ORDER); i++) c[i] = 0;

     for (i = 0; i < NWORDS_ORDER; i++) {
          u = 0;
          for (j = 0; j < NWORDS_ORDER; j++) {
              MULADDADD(a[i], b[j], &u, &c[i+j]);
          }
          c[NWORDS_ORDER+i] = u;
     }
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


void vdiv1271(uint32_t* a) 
{ // GF(p) division by two c = a/2 mod p
  // Representation: 23/26/26/26/26-bit
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


void v2div1271(uint32_t* a) 
{ // GF(p^2) division by two c = a/2 mod p
  // Representation: 23/23/26/26/26/26/26/26/26/26-bit
    digit_t mask;

    mask = (0 - (a[0] & 1)) >> 6;  // if a[0] is odd then mask = 2^26-1, else mask = 0
    
    a[0] += mask;
    a[2] += mask;
    a[4] += mask;
    a[6] += mask;
    a[8] += (mask >> 3);

    a[0] = ((sdigit_t)a[0] >> 1) + ((a[2] & 1) << 25);
    a[2] = ((sdigit_t)a[2] >> 1) + ((a[4] & 1) << 25);
    a[4] = ((sdigit_t)a[4] >> 1) + ((a[6] & 1) << 25);
    a[6] = ((sdigit_t)a[6] >> 1) + ((a[8] & 1) << 25);
    a[8] = ((sdigit_t)a[8] >> 1);

    mask = (0 - (a[1] & 1)) >> 6;  // if a[0] is odd then mask = 2^26-1, else mask = 0
    
    a[1] += mask;
    a[3] += mask;
    a[5] += mask;
    a[7] += mask;
    a[9] += (mask >> 3);

    a[1] = ((sdigit_t)a[1] >> 1) + ((a[3] & 1) << 25);
    a[3] = ((sdigit_t)a[3] >> 1) + ((a[5] & 1) << 25);
    a[5] = ((sdigit_t)a[5] >> 1) + ((a[7] & 1) << 25);
    a[7] = ((sdigit_t)a[7] >> 1) + ((a[9] & 1) << 25);
    a[9] = ((sdigit_t)a[9] >> 1);
}


#ifdef __cplusplus
}
#endif


#endif
