/****************************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: ECC operations over GF(p^2) exploiting endomorphisms
*
* This code is based on the papers:
* [1] "FourQ: four-dimensional decompositions on a Q-curve over the Mersenne prime" 
*     by Craig Costello and Patrick Longa, ASIACRYPT2015 (http://eprint.iacr.org/2015/565).
* [2] "FourQNEON: Faster Elliptic Curve Scalar Multiplications on ARM Processors" 
*     by Patrick Longa, SAC2016 (http://eprint.iacr.org/2016/645).
******************************************************************************************/ 

#include "FourQ_internal.h"
#include "FourQ_params.h"
#include "FourQ_tables.h"
#include "ARM/fp_arm.h"
#include <arm_neon.h>


/***********************************************/
/************* GF(p^2) FUNCTIONS ***************/

void fp2copy1271(f2elm_t a, f2elm_t c)
{ // Copy of a GF(p^2) element, c = a
    fpcopy1271(a[0], c[0]);
    fpcopy1271(a[1], c[1]);
}


void v2copy1271(v2elm_t a, v2elm_t c) 
{ // Copy vectorized GF(p^2) element, c <- a

    c[0] = a[0]; c[1] = a[1]; c[2] = a[2]; c[3] = a[3]; c[4] = a[4]; 
    c[5] = a[5]; c[6] = a[6]; c[7] = a[7]; c[8] = a[8]; c[9] = a[9];
}


void v2zero1271(v2elm_t a) 
{ // Zeroing vectorized GF(p^2) element, a = 0
    
    a[0] = 0; a[1] = 0; a[2] = 0; a[3] = 0; a[4] = 0; 
    a[5] = 0; a[6] = 0; a[7] = 0; a[8] = 0; a[9] = 0;
}


__inline void v2add1271(v2elm_t a, v2elm_t b, v2elm_t c)
{ // Vectorized GF(p^2) addition, c = a+b in GF((2^127-1)^2)       
    v2add1271_a((uint32_t*)a, (uint32_t*)b, (uint32_t*)c);
}


__inline void v2sub1271(v2elm_t a, v2elm_t b, v2elm_t c)
{ // Vectorized GF(p^2) subtraction, c = a-b in GF((2^127-1)^2)      
    v2sub1271_a((uint32_t*)a, (uint32_t*)b, (uint32_t*)c);
}


void v2dblsub1271(v2elm_t a, v2elm_t b, v2elm_t c)
{ // Vectorized GF(p^2) addition followed by subtraction, c = 2a-b in GF((2^127-1)^2)       
    v2dblsub1271_a((uint32_t*)a, (uint32_t*)b, (uint32_t*)c);
}


void v2addsub1271(v2elm_t a, v2elm_t b, v2elm_t c, v2elm_t d)
{ // Vectorized GF(p^2) addition and subtraction, c = a+b, d = a-b in GF((2^127-1)^2)       
    v2addsub1271_a((uint32_t*)a, (uint32_t*)b, (uint32_t*)c, (uint32_t*)d);
}


void v2neg1271(v2elm_t a) 
{ // Vectorized GF(p^2) negation
  // Representation: 23/23/26/26/26/26/26/26/26/26-bit
    
    a[0] = mask_26 - a[0];
    a[1] = mask_26 - a[1];
    a[2] = mask_26 - a[2];
    a[3] = mask_26 - a[3];
    a[4] = mask_26 - a[4];
    a[5] = mask_26 - a[5];
    a[6] = mask_26 - a[6];
    a[7] = mask_26 - a[7];
    a[8] = mask_23 - a[8];
    a[9] = mask_23 - a[9];
}


void v2mul1271(v2elm_t a, v2elm_t b, v2elm_t c)
{ // Vectorized GF(p^2) multiplication, c = a*b in GF((2^127-1)^2)     
    v2mul1271_a((uint32_t*)a, (uint32_t*)b, (uint32_t*)c);
}


void v2sqr1271(v2elm_t a, v2elm_t c)
{ // Vectorized GF(p^2) squaring, c = a^2 in GF((2^127-1)^2)       
    v2sqr1271_a((uint32_t*)a, (uint32_t*)c);
}


#if defined(MIX_ARM_NEON)

void v2muladd1271(v2elm_t a, v2elm_t b, v2elm_t c, v2elm_t d, v2elm_t e, v2elm_t f)
{ // Vectorized GF(p^2) multiplication/addition in GF((2^127-1)^2)      
    v2muladd1271_a((uint32_t*)a, (uint32_t*)b, (uint32_t*)c, (uint32_t*)d, (uint32_t*)e, (uint32_t*)f);
}


void v2mulsub1271(v2elm_t a, v2elm_t b, v2elm_t c, v2elm_t d, v2elm_t e, v2elm_t f)
{ // Vectorized GF(p^2) multiplication/subtraction in GF((2^127-1)^2)      
    v2mulsub1271_a((uint32_t*)a, (uint32_t*)b, (uint32_t*)c, (uint32_t*)d, (uint32_t*)e, (uint32_t*)f);
}

void v2muladdsub1271(v2elm_t a, v2elm_t b, v2elm_t c, v2elm_t d, v2elm_t e, v2elm_t f, v2elm_t g)
{ // Vectorized GF(p^2) multiplication/addition/subtraction in GF((2^127-1)^2)      
    v2muladdsub1271_a((uint32_t*)a, (uint32_t*)b, (uint32_t*)c, (uint32_t*)d, (uint32_t*)e, (uint32_t*)f, (uint32_t*)g);
}


void v2muldblsub1271(v2elm_t a, v2elm_t b, v2elm_t c, v2elm_t d, v2elm_t e, v2elm_t f)
{ // Vectorized GF(p^2) multiplication/addition/subtraction in GF((2^127-1)^2)      
    v2muldblsub1271_a((uint32_t*)a, (uint32_t*)b, (uint32_t*)c, (uint32_t*)d, (uint32_t*)e, (uint32_t*)f);
}


void v2sqradd1271(v2elm_t a, v2elm_t c, v2elm_t d, v2elm_t e, v2elm_t f)
{ // Vectorized GF(p^2) squaring/addition in GF((2^127-1)^2)      
    v2sqradd1271_a((uint32_t*)a, (uint32_t*)c, (uint32_t*)d, (uint32_t*)e, (uint32_t*)f);
}


void v2sqraddsub1271(v2elm_t a, v2elm_t c, v2elm_t d, v2elm_t e, v2elm_t f, v2elm_t g)
{ // Vectorized GF(p^2) squaring/addition/subtraction in GF((2^127-1)^2)      
    v2sqraddsub1271_a((uint32_t*)a, (uint32_t*)c, (uint32_t*)d, (uint32_t*)e, (uint32_t*)f, (uint32_t*)g);
}

#endif


void v2inv1271(v2elm_t a)
{ // Vectorized GF(p^2) inversion, a = (a0-i*a1)/(a0^2+a1^2)
    velm_t a0, a1, t0, t1; 

    from_v2_to_v(a, a0, a1);
    vsqr1271(a0, t0);                // t0 = a0^2
    vsqr1271(a1, t1);                // t1 = a1^2 
    vadd1271(t0, t1, t0);            // t0 = a0^2+a1^2 
    vinv1271(t0);                    // t0 = (a0^2+a1^2)^-1 
    vneg1271(a1);                    // a = a0-i*a1 
    vmul1271(a0, t0, a0);
    vmul1271(a1, t0, a1);            // a = (a0-i*a1)*(a0^2+a1^2)^-1
    from_v_to_v2(a0, a1, a);
}


__inline void clear_words(void* mem, unsigned int nwords)
{ // Clear integer-size digits from memory. "nwords" indicates the number of integer digits to be zeroed.
  // This function uses the volatile type qualifier to inform the compiler not to optimize out the memory clearing.
  // It has been tested with MSVS 2013 and GNU GCC 4.6.3, 4.7.3, 4.8.2 and 4.8.4. Users are responsible for verifying correctness with different compilers.  
  // See "Compliant Solution (C99)" at https://www.securecoding.cert.org/confluence/display/c/MSC06-C.+Beware+of+compiler+optimizations 
	unsigned int i;
	volatile unsigned int *v = mem;

	for (i = 0; i < nwords; i++)
		v[i] = 0;
}


#if (USE_ENDO == true)

// Fixed GF(p^2) constants for the endomorphisms
static v2elm_t ctau1     = {0x3ce74c3, 0x12, 0x3355f3a, 0x0, 0x120c74d, 0xc000, 0xb0ebeb, 0x0, 0x1964de, 0x0};         
static v2elm_t ctaudual1 = {0x2cdf034, 0x11, 0x2a9b677, 0x0, 0x6529ec, 0x3ff4000, 0x3ac8c16, 0x3ffffff, 0x4aa740, 0x7fffff};
static v2elm_t cphi0 = {0x3fffff7,0x366f81a, 0x3ffffff, 0x154db3b, 0x5fff, 0x3294f6, 0x0, 0x1d6460b, 0x0, 0x2553a0};
static v2elm_t cphi1 = {0x7, 0x28296f9, 0x0, 0x3643a78, 0x5000, 0x22cf334, 0x0, 0x2831431, 0x0, 0x62c8ca};
static v2elm_t cphi2 = {0x15, 0x31df391, 0x0, 0x32dc553, 0xf000, 0x1c982c2, 0x0, 0xadb26d, 0x0, 0x78df26};
static v2elm_t cphi3 = {0x3,0x3962ea4, 0x0, 0x10115e9, 0x2000, 0x342a924, 0x0, 0x12475d8, 0x0, 0x5084c6};
static v2elm_t cphi4 = {0x3, 0x2ec6855, 0x0, 0x263248e, 0x3000, 0x2ea4a10, 0x0, 0x15e9e58, 0x0, 0x124404};
static v2elm_t cphi5 = {0xf, 0x1052df3, 0x0, 0x2c874f1, 0xa000, 0x59e669, 0x0, 0x1062863, 0x0, 0x459195};
static v2elm_t cphi6 = {0x18, 0x20a5be7, 0x0, 0x190e9e2, 0x12000, 0xb3ccd3, 0x0, 0x20c50c6, 0x0, 0xb232a};
static v2elm_t cphi7 = {0x23, 0x348781a, 0x0, 0x60c0d7, 0x18000, 0x2a1a66c, 0x0, 0x72678b, 0x0, 0x3963bc};
static v2elm_t cphi8 = {0xf0, 0x35d0ef0, 0x0, 0x94560a, 0xaa000, 0xbe544e, 0x0, 0x2180c5b, 0x0, 0x1f529f};
static v2elm_t cphi9 = {0xbef, 0x36e2505, 0x0, 0x34f9225, 0x870000, 0x375b014, 0x0, 0x273f800, 0x0, 0xfd52e};
static v2elm_t cpsi1 = {0x3e346ef, 0x13a, 0x1fd1d9, 0x0, 0xa02edf, 0xde000, 0x26a0f55, 0x0, 0x2af99e, 0x0};
static v2elm_t cpsi2 = {0x143, 0x203f372, 0x0, 0x37addc3, 0xe4000, 0x1f034c7, 0x0, 0x1ee66a0, 0x0, 0x21b8d0};
static v2elm_t cpsi3 = {0x9, 0x1e73a61, 0x0, 0x39aaf9d, 0x6000, 0x29063a6, 0x0, 0x5875f5, 0x0, 0x4cb26f};
static v2elm_t cpsi4 = {0x3fffff6, 0x218c59e, 0x3ffffff, 0x655062, 0x3ff9fff, 0x16f9c59, 0x3ffffff, 0x3a78a0a, 0x7fffff, 0x334d90};

// Fixed integer constants for the decomposition
// Close "offset" vector
static uint64_t c1  = {0x72482C5251A4559C};
static uint64_t c2  = {0x59F95B0ADD276F6C};
static uint64_t c3  = {0x7DD2D17C4625FA78};
static uint64_t c4  = {0x6BC57DEF56CE8877};
// Optimal basis vectors 
static uint64_t b11 = {0x0906FF27E0A0A196};   
static uint64_t b12 = {0x1363E862C22A2DA0};                                              
static uint64_t b13 = {0x07426031ECC8030F};                                              
static uint64_t b14 = {0x084F739986B9E651};   
static uint64_t b21 = {0x1D495BEA84FCC2D4};
static uint64_t b24 = {0x25DBC5BC8DD167D0};
static uint64_t b31 = {0x17ABAD1D231F0302};
static uint64_t b32 = {0x02C4211AE388DA51};
static uint64_t b33 = {0x2E4D21C98927C49F};
static uint64_t b34 = {0x0A9E6F44C02ECD97};
static uint64_t b41 = {0x136E340A9108C83F};
static uint64_t b42 = {0x3122DF2DC3E0FF32};
static uint64_t b43 = {0x068A49F02AA8A9B5};
static uint64_t b44 = {0x18D5087896DE0AEA};
// Precomputed integers for fast-Babai rounding
static uint64_t ell1[4] = {0x259686E09D1A7D4F, 0xF75682ACE6A6BD66, 0xFC5BB5C5EA2BE5DF, 0x07};
static uint64_t ell2[4] = {0xD1BA1D84DD627AFB, 0x2BD235580F468D8D, 0x8FD4B04CAA6C0F8A, 0x03};
static uint64_t ell3[4] = {0x9B291A33678C203C, 0xC42BD6C965DCA902, 0xD038BF8D0BFFBAF6, 0x00};
static uint64_t ell4[4] = {0x12E5666B77E7FDC0, 0x81CBDC3714983D82, 0x1B073877A22D8410, 0x03};


/***********************************************/
/**********  CURVE/SCALAR FUNCTIONS  ***********/

static __inline void ecc_tau(vpoint_extproj_t P)
{ // Apply tau mapping to a point, P = tau(P)
  // Input: P = (X1:Y1:Z1) on E in twisted Edwards coordinates
  // Output: P = (Xfinal:Yfinal:Zfinal) on Ehat in twisted Edwards coordinates
    v2elm_t t0, t1; 

    v2sqr1271(P->x, t0);                     // t0 = X1^2                   
    v2sqr1271(P->y, t1);                     // t1 = Y1^2
    v2mul1271(P->x, P->y, P->x);             // X = X1*Y1
    v2sqr1271(P->z, P->y);                   // Y = Z1^2
    v2addsub1271(t1, t0, P->z, t0);          // Z = X1^2+Y1^2, t0 = Y1^2-X1^2
    v2mul1271(P->x, t0, P->x);               // X = X1*Y1*(Y1^2-X1^2)
    v2dblsub1271(P->y, t0, P->y);            // Y = 2*Z1^2-(Y1^2-X1^2)
    v2mul1271(P->x, ctau1, P->x);            // Xfinal = X*ctau1
    v2mul1271(P->y, P->z, P->y);             // Yfinal = Y*Z
    v2mul1271(P->z, t0, P->z);               // Zfinal = t0*Z
}


static __inline void ecc_tau_dual(vpoint_extproj_t P)
{ // Apply tau_dual mapping to a point, P = tau_dual(P)
  // Input: P = (X1:Y1:Z1) on Ehat in twisted Edwards coordinates
  // Output: P = (Xfinal,Yfinal,Zfinal,Tafinal,Tbfinal) on E, where Tfinal = Tafinal*Tbfinal,
  //         corresponding to (Xfinal:Yfinal:Zfinal:Tfinal) in extended twisted Edwards coordinates
    v2elm_t t0, t1, t2;
    
    v2sqr1271(P->x, t0);                          // t0 = X1^2
    v2sqr1271(P->z, t2);                          // t2 = Z1^2
    v2sqr1271(P->y, t1);                          // t1 = Y1^2
    v2addsub1271(t1, t0, t0, P->ta);              // t0 = X1^2+Y1^2, Tafinal = Y1^2-X1^2
    v2mul1271(P->x, P->y, P->x);                  // X = X1*Y1
    v2dblsub1271(t2, P->ta, P->z);                // Z = 2*Z1^2-(Y1^2-X1^2)
    v2mul1271(P->x, ctaudual1, P->tb);            // Tbfinal = ctaudual1*X1*X1
    v2mul1271(P->z, P->ta, P->y);                 // Yfinal = Z*Tafinal
    v2mul1271(P->tb, t0, P->x);                   // Xfinal = Tbfinal*t0
    v2mul1271(P->z, t0, P->z);                    // Zfinal = Z*t0 
}


static __inline void ecc_delphidel(vpoint_extproj_t P)
{ // Apply delta_phi_delta mapping to a point, P = delta(phi_W(delta_inv(P))), 
  // where phi_W is the endomorphism on the Weierstrass form.
  // Input: P = (X1:Y1:Z1) on Ehat in twisted Edwards coordinates
  // Output: P = (Xfinal:Yfinal:Zfinal) on Ehat in twisted Edwards coordinates
    v2elm_t t0, t1, t2, t3, t4, t5, t6; 

    v2sqr1271(P->z, t4);                          // t4 = Z1^2
    v2mul1271(P->y, P->z, t3);                    // t3 = Y1*Z1
    v2mul1271(t4, cphi4, t0);                     // t0 = cphi4*t4
    v2sqr1271(P->y, t2);                          // t2 = Y1^2
    v2add1271(t0, t2, t0);                        // t0 = t0+t2
    v2mul1271(t3, cphi3, t1);                     // t1 = cphi3*t3
    v2addsub1271(t0, t1, t0, t5);                 // t0 = t0+t1, t5 = t0-t1
    v2mul1271(t0, P->z, t0);                      // t0 = t0*Z1
    v2mul1271(t3, cphi1, t1);                     // t1 = cphi1*t3
    v2mul1271(t0, t5, t0);                        // t0 = t0*t5
    v2mul1271(t4, cphi2, t5);                     // t5 = cphi2*t4
    v2add1271(t2, t5, t5);                        // t5 = t2+t5
    v2addsub1271(t1, t5, t1, t6);                 // t1 = t1+t5, t6 = t1-t5
    v2mul1271(t6, t1, t6);                        // t6 = t1*t6
    v2mul1271(t6, cphi0, t6);                     // t6 = cphi0*t6
    v2mul1271(P->x, t6, P->x);                    // X = X1*t6
    v2sqr1271(t2, t6);                            // t6 = t2^2
    v2sqr1271(t3, t2);                            // t2 = t3^2
    v2sqr1271(t4, t3);                            // t3 = t4^2
    v2mul1271(t2, cphi8, t1);                     // t1 = cphi8*t2
    v2mul1271(t3, cphi9, t5);                     // t5 = cphi9*t3
    v2add1271(t1, t6, t1);                        // t1 = t1+t6
    v2mul1271(t2, cphi6, t2);                     // t2 = cphi6*t2
    v2mul1271(t3, cphi7, t3);                     // t3 = cphi7*t3
    v2add1271(t1, t5, t1);                        // t1 = t1+t5
    v2add1271(t2, t3, t2);                        // t2 = t2+t3
    v2mul1271(t1, P->y, t1);                      // t1 = Y1*t1
    v2add1271(t6, t2, P->y);                      // Y = t6+t2
    v2mul1271(P->x, t1, P->x);                    // X = X*t1
    v2mul1271(P->y, cphi5, P->y);                 // Y = cphi5*Y
    v2neg1271_felm(P->x);                         // Xfinal = X^p
    v2mul1271(P->y, P->z, P->y);                  // Y = Y*Z1
    v2mul1271(t0, t1, P->z);                      // Z = t0*t1
    v2mul1271(P->y, t0, P->y);                    // Y = Y*t0
    v2neg1271_felm(P->z);                         // Zfinal = Z^p
    v2neg1271_felm(P->y);                         // Yfinal = Y^p 
}


static __inline void ecc_delpsidel(vpoint_extproj_t P)
{ // Apply delta_psi_delta mapping to a point, P = delta(psi_W(delta_inv(P))), 
  // where psi_W is the endomorphism on the Weierstrass form.
  // Input: P = (X1:Y1:Z1) on Ehat in twisted Edwards coordinates
  // Output: P = (Xfinal:Yfinal:Zfinal) on Ehat in twisted Edwards coordinates
    v2elm_t t0, t1, t2; 
    
    v2neg1271_felm(P->x);                         // X = X1^p
    v2neg1271_felm(P->z);                         // Z = Z1^p
    v2neg1271_felm(P->y);                         // Y = Y1^p
    v2sqr1271(P->z, t2);                          // t2 = Z1^p^2
    v2sqr1271(P->x, t0);                          // t0 = X1^p^2
    v2mul1271(P->x, t2, P->x);                    // X = X1^p*Z1^p^2
    v2mul1271(t2, cpsi2, P->z);                   // Z = cpsi2*Z1^p^2
    v2mul1271(t2, cpsi3, t1);                     // t1 = cpsi3*Z1^p^2
    v2mul1271(t2, cpsi4, t2);                     // t2 = cpsi4*Z1^p^2
    v2add1271(t0, P->z, P->z);                    // Z = X1^p^2 + cpsi2*Z1^p^2
    v2add1271(t0, t2, t2);                        // t2 = X1^p^2 + cpsi4*Z1^p^2
    v2add1271(t0, t1, t1);                        // t1 = X1^p^2 + cpsi3*Z1^p^2
    v2neg1271(t2);                                // t2 = -(X1^p^2 + cpsi4*Z1^p^2)
    v2mul1271(P->z, P->y, P->z);                  // Z = Y1^p*(X1^p^2 + cpsi2*Z1^p^2)
    v2mul1271(P->x, t2, P->x);                    // X = -X1^p*Z1^p^2*(X1^p^2 + cpsi4*Z1^p^2)
    v2mul1271(t1, P->z, P->y);                    // Yfinal = t1*Z
    v2mul1271(P->x, cpsi1, P->x);                 // Xfinal = cpsi1*X
    v2mul1271(P->z, t2, P->z);                    // Zfinal = Z*t2 
}


void ecc_psi(vpoint_extproj_t P)
{ // Apply psi mapping to a point, P = psi(P)
  // Input: P = (X1:Y1:Z1) on E in twisted Edwards coordinates
  // Output: P = (Xfinal,Yfinal,Zfinal,Tafinal,Tbfinal) on E, where Tfinal = Tafinal*Tbfinal,
  //         corresponding to (Xfinal:Yfinal:Zfinal:Tfinal) in extended twisted Edwards coordinates

    ecc_tau(P);                            
    ecc_delpsidel(P);                      		
    ecc_tau_dual(P);                        
}


void ecc_phi(vpoint_extproj_t P)
{ // Apply phi mapping to a point, P = phi(P)
  // Input: P = (X1:Y1:Z1) on E in twisted Edwards coordinates
  // Output: P = (Xfinal,Yfinal,Zfinal,Tafinal,Tbfinal) on E, where Tfinal = Tafinal*Tbfinal,
  //         corresponding to (Xfinal:Yfinal:Zfinal:Tfinal) in extended twisted Edwards coordinates

    ecc_tau(P);                            
    ecc_delphidel(P);                      		
    ecc_tau_dual(P);  
}


void ecc_precomp(vpoint_extproj_t P, vpoint_extproj_precomp_t *T)
{ // Generation of the precomputation table used by the variable-base scalar multiplication ecc_mul().
  // Input: P = (X1,Y1,Z1,Ta,Tb), where T1 = Ta*Tb, corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
  // Output: table T containing 8 points: P, P+phi(P), P+psi(P), P+phi(P)+psi(P), P+psi(phi(P)), P+phi(P)+psi(phi(P)), P+psi(P)+psi(phi(P)), P+phi(P)+psi(P)+psi(phi(P))
  // Precomputed points use the representation (X+Y,Y-X,2Z,2dT) corresponding to (X:Y:Z:T) in extended twisted Edwards coordinates
    vpoint_extproj_precomp_t Q, R, S;
    vpoint_extproj_t PP; 

    // Generating Q = phi(P) = (XQ+YQ,YQ-XQ,ZQ,TQ)
    ecccopy(P, PP);
    ecc_phi(PP);
    R1_to_R3(PP, Q);                       // Converting from (X,Y,Z,Ta,Tb) to (X+Y,Y-X,Z,T) 

    // Generating S = psi(Q) = (XS+YS,YS-XS,ZS,TS)
    ecc_psi(PP);  
    R1_to_R3(PP, S);                       // Converting from (X,Y,Z,Ta,Tb) to (X+Y,Y-X,Z,T) 

    // Generating T[0] = P = (XP+YP,YP-XP,2ZP,2dTP) 
    R1_to_R2(P, T[0]);                     // Converting from (X,Y,Z,Ta,Tb) to (X+Y,Y-X,2Z,2dT)

    // Generating R = psi(P) = (XR+YR,YR-XR,ZR,TR)
    ecc_psi(P); 
    R1_to_R3(P, R);                        // Converting from (X,Y,Z,Ta,Tb) to (X+Y,Y-X,Z,T)  

    eccadd_core(T[0], Q, PP);              // T[1] = P+Q using the representations (X,Y,Z,Ta,Tb) <- (X+Y,Y-X,2Z,2dT) + (X+Y,Y-X,Z,T)
    R1_to_R2(PP, T[1]);                    // Converting from (X,Y,Z,Ta,Tb) to (X+Y,Y-X,2Z,2dT)
    eccadd_core(T[0], R, PP);              // T[2] = P+R 
    R1_to_R2(PP, T[2]);
    eccadd_core(T[1], R, PP);              // T[3] = P+Q+R 
    R1_to_R2(PP, T[3]);
    eccadd_core(T[0], S, PP);              // T[4] = P+S 
    R1_to_R2(PP, T[4]);
    eccadd_core(T[1], S, PP);              // T[5] = P+Q+S 
    R1_to_R2(PP, T[5]);
    eccadd_core(T[2], S, PP);              // T[6] = P+R+S 
    R1_to_R2(PP, T[6]);
    eccadd_core(T[3], S, PP);              // T[7] = P+Q+R+S 
    R1_to_R2(PP, T[7]);              
}


static __inline void mul_truncate(uint32_t* s, uint32_t* C, uint32_t* out)       
{ // 256-bit multiplication with truncation for the scalar decomposition
  // Outputs 64-bit value "out" = (uint64_t)((s * C) >> 256).

   mul_truncate_a(s, C, out);
}


void decompose(uint64_t* k, uint64_t* scalars)
{ // Scalar decomposition for the variable-base scalar multiplication
  // Input: scalar in the range [0, 2^256-1].
  // Output: 4 64-bit sub-scalars. 
    uint64_t a1, a2, a3, a4, temp, mask;

    mul_truncate_a((uint32_t*)k, (uint32_t*)ell1, (uint32_t*)&a1);
    mul_truncate_a((uint32_t*)k, (uint32_t*)ell2, (uint32_t*)&a2);
    mul_truncate_a((uint32_t*)k, (uint32_t*)ell3, (uint32_t*)&a3);
    mul_truncate_a((uint32_t*)k, (uint32_t*)ell4, (uint32_t*)&a4);

    temp = k[0] - (uint64_t)a1*b11 - (uint64_t)a2*b21 - (uint64_t)a3*b31 - (uint64_t)a4*b41 + c1;
    mask = ~(0 - (temp & 1));      // If temp is even then mask = 0xFF...FF, else mask = 0
    
    scalars[0] = temp + (mask & b41);
    scalars[1] = (uint64_t)a1*b12 + (uint64_t)a2     - (uint64_t)a3*b32 - (uint64_t)a4*b42 + c2 + (mask & b42);
    scalars[2] = (uint64_t)a3*b33 - (uint64_t)a1*b13 - (uint64_t)a2     + (uint64_t)a4*b43 + c3 - (mask & b43);
    scalars[3] = (uint64_t)a1*b14 - (uint64_t)a2*b24 - (uint64_t)a3*b34 + (uint64_t)a4*b44 + c4 - (mask & b44);
}


void recode(uint64_t* scalars, unsigned int* digits, unsigned int* sign_masks)
{ // Recoding sub-scalars for use in the variable-base scalar multiplication. See Algorithm 1 in "Efficient and Secure Methods for GLV-Based Scalar 
  // Multiplication and their Implementation on GLV-GLS Curves (Extended Version)", A. Faz-Hernandez, P. Longa, and A.H. Sanchez, in Journal
  // of Cryptographic Engineering, Vol. 5(1), 2015.
  // Input: 4 64-bit sub-scalars passed through "scalars", which are obtained after calling decompose().
  // Outputs: "digits" array with 65 nonzero entries. Each entry is in the range [0, 7], corresponding to one entry in the precomputed table.
  //          "sign_masks" array with 65 entries storing the signs for their corresponding digits in "digits". 
  //          Notation: if the corresponding digit > 0 then sign_mask = 0xFF...FF, else if digit < 0 then sign_mask = 0.
    unsigned int i, bit, bit0, carry;
    sign_masks[64] = (unsigned int)-1; 

    for (i = 0; i < 64; i++)
    {
        scalars[0] >>= 1;
        bit0 = (unsigned int)scalars[0] & 1;
        sign_masks[i] = 0 - bit0;

        bit = (unsigned int)scalars[1] & 1;
        carry = (bit0 | bit) ^ bit0; 
        scalars[1] = (scalars[1] >> 1) + (uint64_t)carry; 
        digits[i] = bit;

        bit = (unsigned int)scalars[2] & 1;
        carry = (bit0 | bit) ^ bit0; 
        scalars[2] = (scalars[2] >> 1) + (uint64_t)carry; 
        digits[i] += (bit << 1);

        bit = (unsigned int)scalars[3] & 1;
        carry = (bit0 | bit) ^ bit0; 
        scalars[3] = (scalars[3] >> 1) + (uint64_t)carry; 
        digits[i] += (bit << 2);
    }
    digits[64] = (unsigned int)(scalars[1] + (scalars[2] << 1) + (scalars[3] << 2));
}


void cofactor_clearing(vpoint_extproj_t P)
{ // Co-factor clearing
  // Input: P = (X1,Y1,Z1,Ta,Tb), where T1 = Ta*Tb, corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates.
  // Output: P = 392*P = (Xfinal,Yfinal,Zfinal,Tafinal,Tbfinal), where Tfinal = Tafinal*Tbfinal,
  //         corresponding to (Xfinal:Yfinal:Zfinal:Tfinal) in extended twisted Edwards coordinates
    vpoint_extproj_precomp_t Q;
    
    R1_to_R2(P, Q);            // Converting from (X,Y,Z,Ta,Tb) to (X+Y,Y-X,2Z,2dT)
    eccdouble(P);              // P = 2*P using representations (X,Y,Z,Ta,Tb) <- 2*(X,Y,Z)
    eccadd(Q, P);              // P = P+Q using representations (X,Y,Z,Ta,Tb) <- (X,Y,Z,Ta,Tb) + (X+Y,Y-X,2Z,2dT)
    eccdouble(P);
    eccdouble(P);
    eccdouble(P);
    eccdouble(P);
    eccadd(Q, P);
    eccdouble(P);
    eccdouble(P);
    eccdouble(P);
}


bool ecc_mul(point_t P, digit_t* k, point_t Q, bool clear_cofactor)
{ // Variable-base scalar multiplication Q = k*P using a 4-dimensional decomposition
  // Inputs: scalar "k" in [0, 2^256-1],
  //         point P = (x,y) in affine coordinates,
  //         clear_cofactor = 1 (TRUE) or 0 (FALSE) whether cofactor clearing is required or not, respectively.
  // Output: Q = k*P in affine coordinates (x,y).
  // This function performs point validation and (if selected) cofactor clearing.
    vpoint_t A;
    vpoint_extproj_t R;
    vpoint_extproj_precomp_t S, Table[8];
    uint64_t scalars[NWORDS64_ORDER];
    unsigned int digits[65], sign_masks[65];
    int i;

    point_setup(P, R);                                        // Convert to vectorized representation (X,Y,1,Ta,Tb)
    
    if (ecc_point_validate(R) == false) {                     // Check if point lies on the curve
        return false;
    }
    
    decompose((uint64_t*)k, scalars);                         // Scalar decomposition
    if (clear_cofactor == true) {
        cofactor_clearing(R);
    }
    recode(scalars, digits, sign_masks);                      // Scalar recoding
    ecc_precomp(R, Table);                                    // Precomputation
    table_lookup_1x8(Table, S, digits[64], sign_masks[64]);   // Extract initial point in (X+Y,Y-X,2Z,2dT) representation
    R2_to_R4(S, R);                                           // Conversion to representation (2X,2Y,2Z)
    
    for (i = 63; i >= 0; i--)
    {
        table_lookup_1x8(Table, S, digits[i], sign_masks[i]); // Extract point S in (X+Y,Y-X,2Z,2dT) representation
        eccdouble(R);                                         // P = 2*P using representations (X,Y,Z,Ta,Tb) <- 2*(X,Y,Z)
        eccadd(S, R);                                         // P = P+S using representations (X,Y,Z,Ta,Tb) <- (X,Y,Z,Ta,Tb) + (X+Y,Y-X,2Z,2dT) 
    }
    eccnorm(R, A);                                            // Conversion to affine coordinates (x,y) and modular correction.
    from_ext_to_std(A->x, Q->x);
    from_ext_to_std(A->y, Q->y); 

    return true;
}

#endif


void eccset(point_t P)
{ // Set generator  
  // Output: P = (x,y)

	fp2copy1271((felm_t*)&GENERATOR_x, P->x);    // X1
	fp2copy1271((felm_t*)&GENERATOR_y, P->y);    // Y1
}


__inline void eccnorm(vpoint_extproj_t P, vpoint_t Q)
{ // Normalize a projective point (X1:Y1:Z1), including full reduction
  // Input: P = (X1:Y1:Z1) in twisted Edwards coordinates    
  // Output: Q = (X1/Z1,Y1/Z1), corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
    
    v2inv1271(P->z);                      // Z1 = Z1^-1
    v2mul1271(P->x, P->z, Q->x);          // X1 = X1/Z1
    v2mul1271(P->y, P->z, Q->y);          // Y1 = Y1/Z1
    v2mod1271(Q->x, Q->x); 
    v2mod1271(Q->y, Q->y); 
}


void R1_to_R2(vpoint_extproj_t P, vpoint_extproj_precomp_t Q) 
{ // Conversion from representation (X,Y,Z,Ta,Tb) to (X+Y,Y-X,2Z,2dT), where T = Ta*Tb
  // Input:  P = (X1,Y1,Z1,Ta,Tb), where T1 = Ta*Tb, corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
  // Output: Q = (X1+Y1,Y1-X1,2Z1,2dT1) corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
    
    v2add1271(P->ta, P->ta, Q->t2);                   // T = 2*Ta
    v2add1271(P->x, P->y, Q->xy);                     // QX = X+Y
    v2sub1271(P->y, P->x, Q->yx);                     // QY = Y-X 
    v2mul1271(Q->t2, P->tb, Q->t2);                   // T = 2*T
    v2add1271(P->z, P->z, Q->z2);                     // QZ = 2*Z
    v2mul1271(Q->t2, (digit_t*)&PARAMETER_d, Q->t2);  // QT = 2d*T
}


void R1_to_R3(vpoint_extproj_t P, vpoint_extproj_precomp_t Q)      
{ // Conversion from representation (X,Y,Z,Ta,Tb) to (X+Y,Y-X,Z,T), where T = Ta*Tb 
  // Input:  P = (X1,Y1,Z1,Ta,Tb), where T1 = Ta*Tb, corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
  // Output: Q = (X1+Y1,Y1-X1,Z1,T1) corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates 
    
    v2addsub1271(P->y, P->x, Q->xy, Q->yx);  // XQ = (X1+Y1), YQ = (Y1-X1)
    v2mul1271(P->ta, P->tb, Q->t2);          // TQ = T1
    v2copy1271(P->z, Q->z2);                 // ZQ = Z1 
}


void R2_to_R4(vpoint_extproj_precomp_t P, vpoint_extproj_t Q)      
{ // Conversion from representation (X+Y,Y-X,2Z,2dT) to (2X,2Y,2Z,2dT) 
  // Input:  P = (X1+Y1,Y1-X1,2Z1,2dT1) corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
  // Output: Q = (2X1,2Y1,2Z1) corresponding to (X1:Y1:Z1) in twisted Edwards coordinates 
    
    v2addsub1271(P->xy, P->yx, Q->y, Q->x);  // YQ = 2*Y1, XQ = 2*X1
    v2copy1271(P->z2, Q->z);                 // ZQ = 2*Z1
    v2mod1271_incomplete(Q->x, Q->x);   
    v2mod1271_incomplete(Q->y, Q->y);   
}


void eccdouble(vpoint_extproj_t P)
{ // Point doubling 2P
  // Input: P = (X1:Y1:Z1) in twisted Edwards coordinates
  // Output: 2P = (Xfinal,Yfinal,Zfinal,Tafinal,Tbfinal), where Tfinal = Tafinal*Tbfinal,
  //         corresponding to (Xfinal:Yfinal:Zfinal:Tfinal) in extended twisted Edwards coordinates
    v2elm_t t1, t2; 

#if !defined(MIX_ARM_NEON)
    v2sqr1271(P->x, t1);                              // t1 = X1^2 
    v2sqr1271(P->y, t2);                              // t2 = Y1^2  
    v2add1271(P->x, P->y, P->x);                      // X = X1+Y1
    v2addsub1271(t2, t1, P->tb, t1);                  // Tbfinal = X1^2+Y1^2, t1 = Y1^2-X1^2
    v2sqr1271(P->z, t2);                              // t2 = Z1^2 
    v2sqr1271(P->x, P->ta);                           // Ta = (X1+Y1)^2 
    v2dblsub1271(t2, t1, t2);                         // t2 = 2Z1^2-(Y1^2-X1^2) 
    v2sub1271(P->ta, P->tb, P->ta);                   // Tafinal = 2X1*Y1 = (X1+Y1)^2-(X1^2+Y1^2) 
    v2mul1271(t1, P->tb, P->y);                       // Yfinal = (X1^2+Y1^2)(Y1^2-X1^2) 
    v2mul1271(t2, P->ta, P->x);                       // Xfinal = 2X1*Y1*[2Z1^2-(Y1^2-X1^2)]
    v2mul1271(t1, t2, P->z);                          // Zfinal = (Y1^2-X1^2)[2Z1^2-(Y1^2-X1^2)] 
#else
    v2sqradd1271(P->x, t1, P->x, P->y, P->x);         // t1 = X1^2, X = X1+Y1
    v2sqr1271(P->y, t2);                              // t2 = Y1^2
    v2sqraddsub1271(P->z, P->z, t2, t1, P->tb, t1);   // Z = Z1^2, Tbfinal = X1^2+Y1^2, t1 = Y1^2-X1^2
    v2sqr1271(P->x, P->ta);                           // Ta = (X1+Y1)^2
    v2muldblsub1271(t1, P->tb, P->y, P->z, t1, t2);   // Yfinal = (X1^2+Y1^2)(Y1^2-X1^2), t2 = 2Z1^2-(Y1^2-X1^2)
    v2mulsub1271(t1, t2, P->z, P->ta, P->tb, P->ta);  // Zfinal = (Y1^2-X1^2)[2Z1^2-(Y1^2-X1^2)], Tafinal = 2X1*Y1 = (X1+Y1)^2-(X1^2+Y1^2) 
    v2mul1271(t2, P->ta, P->x);                       // Xfinal = 2X1*Y1*[2Z1^2-(Y1^2-X1^2)]
#endif
}


__inline void eccadd_core(vpoint_extproj_precomp_t P, vpoint_extproj_precomp_t Q, vpoint_extproj_t R)      
{ // Basic point addition R = P+Q or R = P+P
  // Inputs: P = (X1+Y1,Y1-X1,2Z1,2dT1) corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
  //         Q = (X2+Y2,Y2-X2,Z2,T2) corresponding to (X2:Y2:Z2:T2) in extended twisted Edwards coordinates    
  // Output: R = (Xfinal,Yfinal,Zfinal,Tafinal,Tbfinal), where Tfinal = Tafinal*Tbfinal,
  //         corresponding to (Xfinal:Yfinal:Zfinal:Tfinal) in extended twisted Edwards coordinates
    v2elm_t t1, t2; 
          
#if !defined(MIX_ARM_NEON)
    v2mul1271(P->t2, Q->t2, R->z);                            // Z = 2dT1*T2 
    v2mul1271(P->z2, Q->z2, t1);                              // t1 = 2Z1*Z2 
    v2mul1271(P->xy, Q->xy, R->x);                            // X = (X1+Y1)(X2+Y2)
    v2mul1271(P->yx, Q->yx, R->y);                            // Y = (Y1-X1)(Y2-X2)
    v2addsub1271(R->x, R->y, R->ta, R->tb);                   // Tafinal = omega, Tbfinal = beta
    v2addsub1271(t1, R->z, t1, t2);                           // t1 = alpha, t2 = theta
    v2mul1271(R->tb, t2, R->x);                               // Xfinal = beta*theta
    v2mul1271(t1, t2, R->z);                                  // Zfinal = theta*alpha
    v2mul1271(R->ta, t1, R->y);                               // Yfinal = alpha*omega
#else   
    v2mul1271(P->t2, Q->t2, R->z);                            // Z = 2dT1*T2 
    v2mul1271(P->z2, Q->z2, t1);                              // t1 = 2Z1*Z2 
    v2mul1271(P->xy, Q->xy, R->x);                            // X = (X1+Y1)(X2+Y2) 
    v2muladdsub1271(P->yx, Q->yx, R->y, t1, R->z, t1, t2);    // Y = (Y1-X1)(Y2-X2), t1 = alpha, t2 = theta 
    v2muladdsub1271(t1, t2, R->z, R->x, R->y, R->ta, R->tb);  // Zfinal = theta*alpha, Tafinal = omega, Tbfinal = beta
    v2mul1271(t2, R->tb, R->x);                               // Xfinal = beta*theta
    v2mul1271(t1, R->ta, R->y);                               // Yfinal = alpha*omega  
#endif
}


void eccadd(vpoint_extproj_precomp_t Q, vpoint_extproj_t P)      
{ // Complete point addition P = P+Q or P = P+P
  // Inputs: P = (X1,Y1,Z1,Ta,Tb), where T1 = Ta*Tb, corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
  //         Q = (X2+Y2,Y2-X2,2Z2,2dT2) corresponding to (X2:Y2:Z2:T2) in extended twisted Edwards coordinates   
  // Output: P = (Xfinal,Yfinal,Zfinal,Tafinal,Tbfinal), where Tfinal = Tafinal*Tbfinal, 
  //         corresponding to (Xfinal:Yfinal:Zfinal:Tfinal) in extended twisted Edwards coordinates
    
#if !defined(MIX_ARM_NEON)
    vpoint_precomp_t R;
    v2elm_t t1; 

    v2addsub1271(P->y, P->x, R->xy, R->yx);  // XR = (X1+Y1), YR = (Y1-X1)
    v2mul1271(P->ta, P->tb, R->t2);          // TR = T1 
    v2mul1271(Q->z2, P->z, t1);              // t1 = 2Z1*Z2      
    v2mul1271(Q->t2, R->t2, P->z);           // Z = 2dT1*T2 
    v2mul1271(Q->xy, R->xy, P->x);           // X = (X1+Y1)(X2+Y2)
    v2mul1271(Q->yx, R->yx, P->y);           // Y = (Y1-X1)(Y2-X2)
    v2addsub1271(P->x, P->y, P->ta, P->tb);  // Tafinal = omega, Tbfinal = beta
    v2addsub1271(t1, P->z, t1, R->t2);       // t1 = alpha, TR = theta
    v2mul1271(P->tb, R->t2, P->x);           // Xfinal = beta*theta
    v2mul1271(t1, R->t2, P->z);              // Zfinal = theta*alpha
    v2mul1271(P->ta, t1, P->y);              // Yfinal = alpha*omega
#else  
    v2elm_t t1, t2; 

    v2muladdsub1271(P->ta, P->tb, P->ta, P->y, P->x, P->tb, P->y);  // Ta = T1, Tb = (X1+Y1), Y = (Y1-X1) 
    v2mul1271(P->z, Q->z2, t1);                                     // t1 = 2Z1*Z2 
    v2mul1271(P->ta, Q->t2, P->z);                                  // Z = 2dT1*T2 
    v2muladdsub1271(P->tb, Q->xy, P->x, t1, P->z, t1, t2);          // X = (X1+Y1)(X2+Y2), t1 = alpha, t2 = theta 
    v2mul1271(P->y, Q->yx, P->y);                                   // Y = (Y1-X1)(Y2-X2) 
    v2muladdsub1271(t2, t1, P->z, P->x, P->y, P->ta, P->tb);        // Zfinal = theta*alpha, Tafinal = omega, Tbfinal = beta
    v2mul1271(t2, P->tb, P->x);                                     // Xfinal = beta*theta
    v2mul1271(t1, P->ta, P->y);                                     // Yfinal = alpha*omega  
#endif
}


void point_setup(point_t P, vpoint_extproj_t Q)
{ // Point conversion to vectorized representation (X,Y,Z,Ta,Tb) 
  // Input: P = (x,y) in affine coordinates
  // Output: P = (X,Y,1,Ta,Tb), where Ta=X, Tb=Y and T=Ta*Tb, corresponding to (X:Y:Z:T) in extended twisted Edwards coordinates

    from_std_to_ext(P->x, Q->x);
    from_std_to_ext(P->y, Q->y);
    v2copy1271(Q->x, Q->ta);              // Ta = X1
    v2copy1271(Q->y, Q->tb);              // Tb = Y1
    v2zero1271(Q->z); Q->z[0]=1;          // Z1 = 1
}


bool ecc_point_validate(vpoint_extproj_t P)
{ // Point validation: check if point lies on the curve
  // Input: P = (x,y) in affine coordinates, where x, y in [0, 2^127-1].
  // Output: TRUE (1) if point lies on the curve E: -x^2+y^2-1-dx^2*y^2 = 0, FALSE (0) otherwise.
  // SECURITY NOTE: this function does not run in constant time (input point P is assumed to be public).
    v2elm_t t1, t2, t3;
    unsigned int i;

    v2sqr1271(P->y, t1);  
    v2sqr1271(P->x, t2);
    v2sub1271(t1, t2, t3);                      // -x^2 + y^2 
    v2mul1271(t1, t2, t1);                      // x^2*y^2
    v2mul1271((digit_t*)&PARAMETER_d, t1, t2);  // dx^2*y^2
    v2zero1271(t1);  t1[0] = 1;                 // t1 = 1
    v2add1271(t2, t1, t2);                      // 1 + dx^2*y^2
    v2sub1271(t3, t2, t1);                      // -x^2 + y^2 - 1 - dx^2*y^2
    v2mod1271(t1, t1);
    
    for (i = 0; i < 2*VWORDS_FIELD-1; i++) {
        if (t1[i] != 0) return false;
    }
    return true; 
}


static __inline void R5_to_R1(vpoint_precomp_t P, vpoint_extproj_t Q)      
{ // Conversion from representation (x+y,y-x,2dt) to (X,Y,Z,Ta,Tb) 
  // Input:  P = (x1+y1,y1-x1,2dt1) corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates, where Z1=1
  // Output: Q = (x1,y1,z1,x1,y1), where z1=1, corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates 
    
    v2addsub1271(P->xy, P->yx, Q->y, Q->x);  // 2*y1, 2*x1
    v2zero1271(Q->z); Q->z[0]=1;             // ZQ = 1
    v2div1271(Q->x);                         // XQ = x1
    v2div1271(Q->y);                         // YQ = y1 
    v2copy1271(Q->x, Q->ta);                 // TaQ = x1
    v2copy1271(Q->y, Q->tb);                 // TbQ = y1
}


static __inline void eccmadd(vpoint_precomp_t Q, vpoint_extproj_t P)
{ // Mixed point addition P = P+Q or P = P+P
  // Inputs: P = (X1,Y1,Z1,Ta,Tb), where T1 = Ta*Tb, corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
  //         Q = (x2+y2,y2-x2,2dt2) corresponding to (X2:Y2:Z2:T2) in extended twisted Edwards coordinates, where Z2=1  
  // Output: P = (Xfinal,Yfinal,Zfinal,Tafinal,Tbfinal), where Tfinal = Tafinal*Tbfinal, 
  //         corresponding to (Xfinal:Yfinal:Zfinal:Tfinal) in extended twisted Edwards coordinates 
    v2elm_t t1, t2;

#if !defined(MIX_ARM_NEON)    
    v2mul1271(P->ta, P->tb, P->ta);                       // Ta = T1
    v2add1271(P->z, P->z, t1);                            // t1 = 2Z1        
    v2mul1271(P->ta, Q->t2, P->ta);                       // Ta = 2dT1*t2 
    v2addsub1271(P->y, P->x, P->z, P->tb);                // Z = (X1+Y1), Tb = (Y1-X1) 
    v2addsub1271(t1, P->ta, t1, t2);                      // t1 = alpha, t2 = theta
    v2mul1271(Q->xy, P->z, P->ta);                        // Ta = (X1+Y1)(x2+y2)
    v2mul1271(Q->yx, P->tb, P->x);                        // X = (Y1-X1)(y2-x2)
    v2addsub1271(P->ta, P->x, P->ta, P->tb);              // Tafinal = omega, Tbfinal = beta
    v2mul1271(t1, t2, P->z);                              // Zfinal = theta*alpha
    v2mul1271(P->tb, t2, P->x);                           // Xfinal = beta*theta
    v2mul1271(P->ta, t1, P->y);                           // Yfinal = alpha*omega
#else    
    v2muladd1271(P->ta, P->tb, P->ta, P->z, P->z, t1);              // Ta = T1, t1 = 2Z1
    v2muladdsub1271(P->ta, Q->t2, P->ta, P->y, P->x, P->z, P->tb);  // Ta = 2dT1*t2, Z = (X1+Y1), Tb = (Y1-X1) 
    v2muladdsub1271(Q->xy, P->z, P->y, t1, P->ta, t1, t2);          // Y = (X1+Y1)(x2+y2), t1 = alpha, t2 = theta
    v2mul1271(Q->yx, P->tb, P->x);                                  // X = (Y1-X1)(y2-x2)
    v2muladdsub1271(t1, t2, P->z, P->y, P->x, P->ta, P->tb);        // Zfinal = theta*alpha, Tafinal = omega, Tbfinal = beta
    v2mul1271(P->tb, t2, P->x);                                     // Xfinal = beta*theta
    v2mul1271(P->ta, t1, P->y);                                     // Yfinal = alpha*omega
#endif
}


bool ecc_mul_fixed(digit_t* k, point_t Q)
{ // Fixed-base scalar multiplication Q = k*G, where G is the generator. FIXED_BASE_TABLE stores v*2^(w-1) = 80 multiples of G.
  // Inputs: scalar "k" in [0, 2^256-1].
  // Output: Q = k*G in affine coordinates (x,y).
  // The function is based on the modified LSB-set comb method, which converts the scalar to an odd signed representation
  // with (bitlength(order)+w*v) digits.
    unsigned int j, w = W_FIXEDBASE, v = V_FIXEDBASE, d = D_FIXEDBASE, e = E_FIXEDBASE;
    unsigned int digit = 0, digits[NBITS_ORDER_PLUS_ONE+(W_FIXEDBASE*V_FIXEDBASE)-1] = {0}; 
	digit_t temp[NWORDS_ORDER];
    vpoint_t A;
    vpoint_extproj_t R;
    vpoint_precomp_t S;
    int i, ii;

	modulo_order(k, temp);                                      // temp = k mod (order) 
	conversion_to_odd(temp, temp);                              // Converting scalar to odd using the prime subgroup order
	mLSB_set_recode((uint64_t*)temp, digits);                   // Scalar recoding

    // Extracting initial digit 
    digit = digits[w*d-1];
    for (i = (int)((w-1)*d-1); i >= (int)(2*d-1); i = i-d)           
    {
        digit = 2*digit + digits[i];
    }
    // Initialize R = (x+y,y-x,2dt) with a point from the table
    table_lookup_fixed_base(((vpoint_precomp_t*)&FIXED_BASE_TABLE)+(v-1)*(1 << (w-1)), S, digit, digits[d-1]);
    R5_to_R1(S, R);                                             // Converting to representation (X:Y:1:Ta:Tb)

    for (j = 0; j < (v-1); j++)
    {
        digit = digits[w*d-(j+1)*e-1];
        for (i = (int)((w-1)*d-(j+1)*e-1); i >= (int)(2*d-(j+1)*e-1); i = i-d)           
        {
            digit = 2*digit + digits[i];
        }
        // Extract point in (x+y,y-x,2dt) representation
        table_lookup_fixed_base(((vpoint_precomp_t*)&FIXED_BASE_TABLE)+(v-j-2)*(1 << (w-1)), S, digit, digits[d-(j+1)*e-1]);   
        eccmadd(S, R);                                          // R = R+S using representations (X,Y,Z,Ta,Tb) <- (X,Y,Z,Ta,Tb) + (x+y,y-x,2dt) 
    }

    for (ii = (e-2); ii >= 0; ii--)
    {
        eccdouble(R);                                           // R = 2*R using representations (X,Y,Z,Ta,Tb) <- 2*(X,Y,Z)
        for (j = 0; j < v; j++)
        {
            digit = digits[w*d-j*e+ii-e];
            for (i = (int)((w-1)*d-j*e+ii-e); i >= (int)(2*d-j*e+ii-e); i = i-d)           
            {
                digit = 2*digit + digits[i];
            }
            // Extract point in (x+y,y-x,2dt) representation
            table_lookup_fixed_base(((vpoint_precomp_t*)&FIXED_BASE_TABLE)+(v-j-1)*(1 << (w-1)), S, digit, digits[d-j*e+ii-e]); 
            eccmadd(S, R);                                      // R = R+S using representations (X,Y,Z,Ta,Tb) <- (X,Y,Z,Ta,Tb) + (x+y,y-x,2dt)
        }        
    }     
    eccnorm(R, A);                                              // Conversion to affine coordinates (x,y) and modular correction. 
    from_ext_to_std(A->x, Q->x);
    from_ext_to_std(A->y, Q->y); 
    
    return true;
}


void mLSB_set_recode(uint64_t* scalar, unsigned int *digits)
{ // Computes the modified LSB-set representation of a scalar
  // Inputs: scalar in [0, order-1], where the order of FourQ's subgroup is 246 bits.
  // Output: digits, where the first "d" values (from index 0 to (d-1)) store the signs for the recoded values using the convention: -1 (negative), 0 (positive), and
  //         the remaining values (from index d to (l-1)) store the recoded values in mLSB-set representation, excluding their sign, 
  //         where l = d*w and d = ceil(bitlength(order)/(w*v))*v. The values v and w are fixed and must be in the range [1, 10] (see FourQ.h); they determine the size 
  //         of the precomputed table "FIXED_BASE_TABLE" used by ecc_mul_fixed(). 
    unsigned int i, j, d = D_FIXEDBASE, l = L_FIXEDBASE;
    uint64_t temp, carry;
    
    digits[d-1] = 0;

    // Shift scalar to the right by 1   
    for (j = 0; j < (NWORDS64_ORDER-1); j++) {
        SHIFTR(scalar[j+1], scalar[j], 1, scalar[j], RADIX64);
    }
    scalar[NWORDS64_ORDER-1] >>= 1;

    for (i = 0; i < (d-1); i++)
    {
        digits[i] = (unsigned int)((scalar[0] & 1) - 1);  // Convention for the "sign" row: 
                                                          // if scalar_(i+1) = 0 then digit_i = -1 (negative), else if scalar_(i+1) = 1 then digit_i = 0 (positive)
        // Shift scalar to the right by 1   
        for (j = 0; j < (NWORDS64_ORDER-1); j++) {
            SHIFTR(scalar[j+1], scalar[j], 1, scalar[j], RADIX64);
        }
        scalar[NWORDS64_ORDER-1] >>= 1;
    } 

    for (i = d; i < l; i++)
    {
        digits[i] = (unsigned int)(scalar[0] & 1);        // digits_i = k mod 2. Sign is determined by the "sign" row

        // Shift scalar to the right by 1  
        for (j = 0; j < (NWORDS64_ORDER-1); j++) {
            SHIFTR(scalar[j+1], scalar[j], 1, scalar[j], RADIX64);
        }
        scalar[NWORDS64_ORDER-1] >>= 1;

        temp = (0 - digits[i-(i/d)*d]) & digits[i];       // if (digits_i=0 \/ 1) then temp = 0, else if (digits_i=-1) then temp = 1 
            
        // floor(scalar/2) + temp
        scalar[0] = scalar[0] + temp;
        carry = (temp & (uint64_t)is_digit_zero_ct((digit_t)scalar[0]));       // carry = (scalar[0] < temp);
        for (j = 1; j < NWORDS64_ORDER; j++)
        {
            scalar[j] = scalar[j] + carry; 
            carry = (carry & (uint64_t)is_digit_zero_ct((digit_t)scalar[j]));  // carry = (scalar[j] < temp);
        }
    } 
    return;              
}


static __inline void eccneg_extproj_precomp(vpoint_extproj_precomp_t P, vpoint_extproj_precomp_t Q)
{ // Point negation
  // Input : point P in coordinates (X+Y,Y-X,2Z,2dT)
  // Output: point Q = -P = (Y-X,X+Y,2Z,-2dT)
    v2copy1271(P->t2, Q->t2);
    v2copy1271(P->xy, Q->yx);
    v2copy1271(P->yx, Q->xy);
    v2copy1271(P->z2, Q->z2);
    v2neg1271(Q->t2);
}


static __inline void eccneg_precomp(vpoint_precomp_t P, vpoint_precomp_t Q)
{ // Point negation
  // Input : point P in coordinates (x+y,y-x,2dt)
  // Output: point Q = -P = (y-x,x+y,-2dt)
    v2copy1271(P->t2, Q->t2);
    v2copy1271(P->xy, Q->yx);
    v2copy1271(P->yx, Q->xy);
    v2neg1271(Q->t2);
}


bool ecc_mul_double(digit_t* k, point_t Q, digit_t* l, point_t R)
{ // Double scalar multiplication R = k*G + l*Q, where the G is the generator. Uses DOUBLE_SCALAR_TABLE, which contains multiples of G, Phi(G), Psi(G) and Phi(Psi(G)).
  // Inputs: point Q in affine coordinates,
  //         Scalars "k" and "l" in [0, 2^256-1].
  // Output: R = k*G + l*Q in affine coordinates (x,y).
  // The function uses wNAF with interleaving.
    vpoint_t A;

  // SECURITY NOTE: this function is intended for a non-constant-time operation such as signature verification. 

#if (USE_ENDO == true)
    unsigned int position;
    int i, digits_k1[65] = {0}, digits_k2[65] = {0}, digits_k3[65] = {0}, digits_k4[65] = {0};
    int digits_l1[65] = {0}, digits_l2[65] = {0}, digits_l3[65] = {0}, digits_l4[65] = {0};
    vpoint_precomp_t V;
    vpoint_extproj_t Q1, Q2, Q3, Q4, T; 
    vpoint_extproj_precomp_t U, Q_table1[NPOINTS_DOUBLEMUL_WQ], Q_table2[NPOINTS_DOUBLEMUL_WQ], Q_table3[NPOINTS_DOUBLEMUL_WQ], Q_table4[NPOINTS_DOUBLEMUL_WQ];
    uint64_t k_scalars[4], l_scalars[4];
    
    point_setup(Q, Q1);                                        // Convert to representation (X,Y,1,Ta,Tb)
    
    if (ecc_point_validate(Q1) == false) {                     // Check if point lies on the curve
        return false;
    }
    
    // Computing endomorphisms over point Q
    ecccopy(Q1, Q2);
    ecc_phi(Q2);
    ecccopy(Q1, Q3);    
    ecc_psi(Q3); 
    ecccopy(Q2, Q4); 
    ecc_psi(Q4);  
    
    decompose((uint64_t*)k, k_scalars);                        // Scalar decomposition
    decompose((uint64_t*)l, l_scalars);  
    wNAF_recode(k_scalars[0], WP_DOUBLEBASE, digits_k1);       // Scalar recoding
    wNAF_recode(k_scalars[1], WP_DOUBLEBASE, digits_k2);
    wNAF_recode(k_scalars[2], WP_DOUBLEBASE, digits_k3);
    wNAF_recode(k_scalars[3], WP_DOUBLEBASE, digits_k4);
    wNAF_recode(l_scalars[0], WQ_DOUBLEBASE, digits_l1);      
    wNAF_recode(l_scalars[1], WQ_DOUBLEBASE, digits_l2);
    wNAF_recode(l_scalars[2], WQ_DOUBLEBASE, digits_l3);
    wNAF_recode(l_scalars[3], WQ_DOUBLEBASE, digits_l4);
    ecc_precomp_double(Q1, Q_table1, NPOINTS_DOUBLEMUL_WQ);    // Precomputation
    ecc_precomp_double(Q2, Q_table2, NPOINTS_DOUBLEMUL_WQ); 
    ecc_precomp_double(Q3, Q_table3, NPOINTS_DOUBLEMUL_WQ); 
    ecc_precomp_double(Q4, Q_table4, NPOINTS_DOUBLEMUL_WQ); 

    v2zero1271(T->x);                                          // Initialize T as the neutral point (0:1:1)
    v2zero1271(T->y); T->y[0] = 1; 
    v2zero1271(T->z); T->z[0] = 1;     

    for (i = 64; i >= 0; i--)
    {   
        eccdouble(T);                                          // Double (X_T,Y_T,Z_T,Ta_T,Tb_T) = 2(X_T,Y_T,Z_T,Ta_T,Tb_T)
        if (digits_l1[i] < 0) {
            position = (-digits_l1[i])/2;                      
            eccneg_extproj_precomp(Q_table1[position], U);     // Load and negate U = (X_U,Y_U,Z_U,Td_U) <- -(X+Y,Y-X,2Z,2dT) from a point in the precomputed table 
            eccadd(U, T);                                      // T = T+U = (X_T,Y_T,Z_T,Ta_T,Tb_T) = (X_T,Y_T,Z_T,Ta_T,Tb_T) + (X_U,Y_U,Z_U,Td_U) 
        } else if (digits_l1[i] > 0) {            
            position = (digits_l1[i])/2;                       // Take U = (X_U,Y_U,Z_U,Td_U) <- (X+Y,Y-X,2Z,2dT) from a point in the precomputed table
            eccadd(Q_table1[position], T);                     // T = T+U = (X_T,Y_T,Z_T,Ta_T,Tb_T) = (X_T,Y_T,Z_T,Ta_T,Tb_T) + (X_U,Y_U,Z_U,Td_U) 
        }                                          
        if (digits_l2[i] < 0) {
            position = (-digits_l2[i])/2;                      
            eccneg_extproj_precomp(Q_table2[position], U);      
            eccadd(U, T);                                
        } else if (digits_l2[i] > 0) {            
            position = (digits_l2[i])/2;                       
            eccadd(Q_table2[position], T);               
        }                                        
        if (digits_l3[i] < 0) {
            position = (-digits_l3[i])/2;                      
            eccneg_extproj_precomp(Q_table3[position], U);      
            eccadd(U, T);                                
        } else if (digits_l3[i] > 0) {            
            position = (digits_l3[i])/2;                       
            eccadd(Q_table3[position], T);               
        }                                        
        if (digits_l4[i] < 0) {
            position = (-digits_l4[i])/2;                      
            eccneg_extproj_precomp(Q_table4[position], U);      
            eccadd(U, T);                                
        } else if (digits_l4[i] > 0) {            
            position = (digits_l4[i])/2;                       
            eccadd(Q_table4[position], T);               
        }

        if (digits_k1[i] < 0) {
            position = (-digits_k1[i])/2;                      
            eccneg_precomp(((vpoint_precomp_t*)&DOUBLE_SCALAR_TABLE)[position], V);    // Load and negate V = (X_V,Y_V,Z_V,Td_V) <- -(x+y,y-x,2dt) from a point in the precomputed table 
            eccmadd(V, T);                                                             // T = T+V = (X_T,Y_T,Z_T,Ta_T,Tb_T) = (X_T,Y_T,Z_T,Ta_T,Tb_T) + (X_V,Y_V,Z_V,Td_V) 
        } else if (digits_k1[i] > 0) {            
            position = (digits_k1[i])/2;                                               // Take V = (X_V,Y_V,Z_V,Td_V) <- (x+y,y-x,2dt) from a point in the precomputed table
            eccmadd(((vpoint_precomp_t*)&DOUBLE_SCALAR_TABLE)[position], T);           // T = T+V = (X_T,Y_T,Z_T,Ta_T,Tb_T) = (X_T,Y_T,Z_T,Ta_T,Tb_T) + (X_V,Y_V,Z_V,Td_V) 
        }
        if (digits_k2[i] < 0) {
            position = (-digits_k2[i])/2;                      
            eccneg_precomp(((vpoint_precomp_t*)&DOUBLE_SCALAR_TABLE)[NPOINTS_DOUBLEMUL_WP+position], V);              
            eccmadd(V, T);                              
        } else if (digits_k2[i] > 0) {            
            position = (digits_k2[i])/2;                       
            eccmadd(((vpoint_precomp_t*)&DOUBLE_SCALAR_TABLE)[NPOINTS_DOUBLEMUL_WP+position], T);               
        }
        if (digits_k3[i] < 0) {
            position = (-digits_k3[i])/2;                      
            eccneg_precomp(((vpoint_precomp_t*)&DOUBLE_SCALAR_TABLE)[2*NPOINTS_DOUBLEMUL_WP+position], V);              
            eccmadd(V, T);                              
        } else if (digits_k3[i] > 0) {            
            position = (digits_k3[i])/2;                       
            eccmadd(((vpoint_precomp_t*)&DOUBLE_SCALAR_TABLE)[2*NPOINTS_DOUBLEMUL_WP+position], T);               
        }
        if (digits_k4[i] < 0) {
            position = (-digits_k4[i])/2;                      
            eccneg_precomp(((vpoint_precomp_t*)&DOUBLE_SCALAR_TABLE)[3*NPOINTS_DOUBLEMUL_WP+position], V);              
            eccmadd(V, T);                              
        } else if (digits_k4[i] > 0) {            
            position = (digits_k4[i])/2;                       
            eccmadd(((vpoint_precomp_t*)&DOUBLE_SCALAR_TABLE)[3*NPOINTS_DOUBLEMUL_WP+position], T);               
        }
    }

#else
	point_t B;
	vpoint_extproj_t T;
	vpoint_extproj_precomp_t S;

	if (ecc_mul(Q, l, B, false) == false) {
		return false;
	}
	point_setup(B, T);
	R1_to_R2(T, S);

	ecc_mul_fixed(k, B);
	point_setup(B, T);
	eccadd(S, T);
#endif
    eccnorm(T, A);                                             // Conversion to affine coordinates (x,y) and modular correction. 
    from_ext_to_std(A->x, R->x);
    from_ext_to_std(A->y, R->y); 
    
    return true;
}


void ecc_precomp_double(vpoint_extproj_t P, vpoint_extproj_precomp_t* Table, unsigned int npoints)
{ // Generation of the precomputation table used internally by the double scalar multiplication function ecc_mul_double().  
  // Inputs: point P in representation (X,Y,Z,Ta,Tb),
  //         Table with storage for npoints, 
  //         number of points "npoints".
  // Output: Table containing multiples of the base point P using representation (X+Y,Y-X,2Z,2dT).
    vpoint_extproj_t Q;
    vpoint_extproj_precomp_t PP;
    unsigned int i; 
           
    R1_to_R2(P, Table[0]);                     // Precomputed point Table[0] = P in coordinates (X+Y,Y-X,2Z,2dT)
    eccdouble(P);                              // A = 2*P in (X,Y,Z,Ta,Tb)
    R1_to_R3(P, PP);                           // Converting from (X,Y,Z,Ta,Tb) to (X+Y,Y-X,Z,T) 
    
    for (i = 1; i < npoints; i++) {
        eccadd_core(Table[i-1], PP, Q);        // Table[i] = Table[i-1]+2P using the representations (X,Y,Z,Ta,Tb) <- (X+Y,Y-X,2Z,2dT) + (X+Y,Y-X,Z,T)
        R1_to_R2(Q, Table[i]);                 // Converting from (X,Y,Z,Ta,Tb) to (X+Y,Y-X,2Z,2dT)
    }
    
    return;
}


void wNAF_recode(uint64_t scalar, unsigned int w, int* digits)
{ // Computes wNAF recoding of a scalar, where digits are in set {0,+-1,+-3,...,+-(2^(w-1)-1)}
    unsigned int i;
    int digit, index = 0; 
    int val1 = (int)(1 << (w-1)) - 1;                  // 2^(w-1) - 1
    int val2 = (int)(1 << w);                          // 2^w;
    uint64_t k = scalar, mask = (uint64_t)val2 - 1;    // 2^w - 1 

    while (k != 0)
    {
        digit = (int)(k & 1); 

        if (digit == 0) {                         
            k >>= 1;                 // Shift scalar to the right by 1
            digits[index] = 0;
        } else {
            digit = (int)(k & mask); 
            k >>= w;                 // Shift scalar to the right by w            

            if (digit > val1) {
                digit -= val2; 
            }
            if (digit < 0) {         // scalar + 1
                k += 1;
            }
            digits[index] = digit; 
                       
            if (k != 0) {            // Check if scalar != 0
                for (i = 0; i < (w-1); i++) 
                {     
                    index++; 
                    digits[index] = 0;
                }
            }
        }
        index++;
    } 
    return;
}

