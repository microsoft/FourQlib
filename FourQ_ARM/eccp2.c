/***********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: ECC operations over GF(p^2) exploiting endomorphisms
*
* This code is based on the paper "FourQ: four-dimensional decompositions on a 
* Q-curve over the Mersenne prime" by Craig Costello and Patrick Longa, in Advances 
* in Cryptology - ASIACRYPT, 2015.
* Preprint available at http://eprint.iacr.org/2015/565.
************************************************************************************/

#include "FourQ_internal.h"
#include "FourQ_params.h"
#include "FourQ_tables.h"
#include "ARM/fp_arm.h"


/***********************************************/
/************* GF(p^2) FUNCTIONS ***************/

void fp2copy1271(f2elm_t a, f2elm_t c)
{// Copy of a GF(p^2) element, c = a
    fpcopy1271(a[0], c[0]);
    fpcopy1271(a[1], c[1]);
}


void fp2zero1271(f2elm_t a)
{// Zeroing a GF(p^2) element, a = 0
    fpzero1271(a[0]);
    fpzero1271(a[1]);
}


void fp2neg1271(f2elm_t a)
{// GF(p^2) negation, a = -a in GF((2^127-1)^2)
    fpneg1271(a[0]);
    fpneg1271(a[1]);
}


void fp2sqr1271(f2elm_t a, f2elm_t c)
{// GF(p^2) squaring, c = a^2 in GF((2^127-1)^2)

    fp2sqr1271_a(a, c);
}


void fp2mul1271(f2elm_t a, f2elm_t b, f2elm_t c)
{// GF(p^2) multiplication, c = a*b in GF((2^127-1)^2)
       
    fp2mul1271_a(a, b, c);
}


void fp2add1271(f2elm_t a, f2elm_t b, f2elm_t c)
{// GF(p^2) addition, c = a+b in GF((2^127-1)^2)

    fp2add1271_a(a, b, c);
}


void fp2sub1271(f2elm_t a, f2elm_t b, f2elm_t c)
{// GF(p^2) subtraction, c = a-b in GF((2^127-1)^2) 

    fp2sub1271_a(a, b, c); 
}


void fp2inv1271(f2elm_t a)
{// GF(p^2) inversion, a = (a0-i*a1)/(a0^2+a1^2)
    f2elm_t t1;

    fpsqr1271(a[0], t1[0]);             // t10 = a0^2
    fpsqr1271(a[1], t1[1]);             // t11 = a1^2
    fpadd1271(t1[0], t1[1], t1[0]);     // t10 = a0^2+a1^2
    fpinv1271(t1[0]);                   // t10 = (a0^2+a1^2)^-1
    fpneg1271(a[1]);                    // a = a0-i*a1
    fpmul1271(a[0], t1[0], a[0]);
    fpmul1271(a[1], t1[0], a[1]);       // a = (a0-i*a1)*(a0^2+a1^2)^-1
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
static uint64_t ctau1[4]     = {0x74DCD57CEBCE74C3, 0x1964DE2C3AFAD20C, 0x12, 0x0C};         
static uint64_t ctaudual1[4] = {0x9ECAA6D9DECDF034, 0x4AA740EB23058652, 0x11, 0x7FFFFFFFFFFFFFF4};
static uint64_t cphi0[4] = {0xFFFFFFFFFFFFFFF7, 0x05, 0x4F65536CEF66F81A, 0x2553A0759182C329};
static uint64_t cphi1[4] = {0x07, 0x05, 0x334D90E9E28296F9, 0x62C8CAA0C50C62CF};
static uint64_t cphi2[4] = {0x15, 0x0F, 0x2C2CB7154F1DF391, 0x78DF262B6C9B5C98};
static uint64_t cphi3[4] = {0x03, 0x02, 0x92440457A7962EA4, 0x5084C6491D76342A};
static uint64_t cphi4[4] = {0x03, 0x03, 0xA1098C923AEC6855, 0x12440457A7962EA4};
static uint64_t cphi5[4] = {0x0F, 0x0A, 0x669B21D3C5052DF3, 0x459195418A18C59E};
static uint64_t cphi6[4] = {0x18, 0x12, 0xCD3643A78A0A5BE7, 0x0B232A8314318B3C};
static uint64_t cphi7[4] = {0x23, 0x18, 0x66C183035F48781A, 0x3963BC1C99E2EA1A};
static uint64_t cphi8[4] = {0xF0, 0xAA, 0x44E251582B5D0EF0, 0x1F529F860316CBE5};
static uint64_t cphi9[4] = {0xBEF, 0x870, 0x14D3E48976E2505, 0xFD52E9CFE00375B};
static uint64_t cpsi1[4] = {0xEDF07F4767E346EF, 0x2AF99E9A83D54A02, 0x13A, 0xDE};
static uint64_t cpsi2[4] = {0x143, 0xE4, 0x4C7DEB770E03F372, 0x21B8D07B99A81F03};
static uint64_t cpsi3[4] = {0x09, 0x06, 0x3A6E6ABE75E73A61, 0x4CB26F161D7D6906};
static uint64_t cpsi4[4] = {0xFFFFFFFFFFFFFFF6, 0x7FFFFFFFFFFFFFF9, 0xC59195418A18C59E, 0x334D90E9E28296F9};

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

static __inline void ecc_tau(point_extproj_t P)
{ // Apply tau mapping to a point, P = tau(P)
  // Input: P = (X1:Y1:Z1) on E in twisted Edwards coordinates
  // Output: P = (Xfinal:Yfinal:Zfinal) on Ehat in twisted Edwards coordinates
    f2elm_t t0, t1; 

    fp2sqr1271(P->x, t0);                     // t0 = X1^2
    fp2sqr1271(P->y, t1);                     // t1 = Y1^2
    fp2mul1271(P->x, P->y, P->x);             // X = X1*Y1
    fp2sqr1271(P->z, P->y);                   // Y = Z1^2
    fp2add1271(t0, t1, P->z);                 // Z = X1^2+Y1^2
    fp2add1271(P->y, P->y, P->y);             // Y = 2*Z1^2
    fp2sub1271(t1, t0, t0);                   // t0 = Y1^2-X1^2
    fp2mul1271(P->x, t0, P->x);               // X = X1*Y1*(Y1^2-X1^2)
    fp2sub1271(P->y, t0, P->y);               // Y = 2*Z1^2-(Y1^2-X1^2)
    fp2mul1271(P->x, (felm_t*)&ctau1, P->x);  // Xfinal = X*ctau1
    fp2mul1271(P->y, P->z, P->y);             // Yfinal = Y*Z
    fp2mul1271(P->z, t0, P->z);               // Zfinal = t0*Z
}


static __inline void ecc_tau_dual(point_extproj_t P)
{ // Apply tau_dual mapping to a point, P = tau_dual(P)
  // Input: P = (X1:Y1:Z1) on Ehat in twisted Edwards coordinates
  // Output: P = (Xfinal,Yfinal,Zfinal,Tafinal,Tbfinal) on E, where Tfinal = Tafinal*Tbfinal,
  //         corresponding to (Xfinal:Yfinal:Zfinal:Tfinal) in extended twisted Edwards coordinates
    f2elm_t t0, t1;

    fp2sqr1271(P->x, t0);                          // t0 = X1^2
    fp2sqr1271(P->z, P->ta);                       // Ta = Z1^2
    fp2sqr1271(P->y, t1);                          // t1 = Y1^2
    fp2add1271(P->ta, P->ta, P->z);                // Z = 2*Z1^2
    fp2sub1271(t1, t0, P->ta);                     // Tafinal = Y1^2-X1^2
    fp2add1271(t0, t1, t0);                        // t0 = X1^2+Y1^2
    fp2mul1271(P->x, P->y, P->x);                  // X = X1*Y1
    fp2sub1271(P->z, P->ta, P->z);                 // Z = 2*Z1^2-(Y1^2-X1^2)
    fp2mul1271(P->x, (felm_t*)&ctaudual1, P->tb);  // Tbfinal = ctaudual1*X1*X1
    fp2mul1271(P->z, P->ta, P->y);                 // Yfinal = Z*Tafinal
    fp2mul1271(P->tb, t0, P->x);                   // Xfinal = Tbfinal*t0
    fp2mul1271(P->z, t0, P->z);                    // Zfinal = Z*t0
}


static __inline void ecc_delphidel(point_extproj_t P)
{ // Apply delta_phi_delta mapping to a point, P = delta(phi_W(delta_inv(P))), 
  // where phi_W is the endomorphism on the Weierstrass form.
  // Input: P = (X1:Y1:Z1) on Ehat in twisted Edwards coordinates
  // Output: P = (Xfinal:Yfinal:Zfinal) on Ehat in twisted Edwards coordinates
    f2elm_t t0, t1, t2, t3, t4, t5, t6; 

    fp2sqr1271(P->z, t4);                          // t4 = Z1^2
    fp2mul1271(P->y, P->z, t3);                    // t3 = Y1*Z1
    fp2mul1271(t4, (felm_t*)&cphi4, t0);           // t0 = cphi4*t4
    fp2sqr1271(P->y, t2);                          // t2 = Y1^2
    fp2add1271(t0, t2, t0);                        // t0 = t0+t2
    fp2mul1271(t3, (felm_t*)&cphi3, t1);           // t1 = cphi3*t3
    fp2sub1271(t0, t1, t5);                        // t5 = t0-t1
    fp2add1271(t0, t1, t0);                        // t0 = t0+t1
    fp2mul1271(t0, P->z, t0);                      // t0 = t0*Z1
    fp2mul1271(t3, (felm_t*)&cphi1, t1);           // t1 = cphi1*t3
    fp2mul1271(t0, t5, t0);                        // t0 = t0*t5
    fp2mul1271(t4, (felm_t*)&cphi2, t5);           // t5 = cphi2*t4
    fp2add1271(t2, t5, t5);                        // t5 = t2+t5
    fp2sub1271(t1, t5, t6);                        // t6 = t1-t5
    fp2add1271(t1, t5, t1);                        // t1 = t1+t5
    fp2mul1271(t6, t1, t6);                        // t6 = t1*t6
    fp2mul1271(t6, (felm_t*)&cphi0, t6);           // t6 = cphi0*t6
    fp2mul1271(P->x, t6, P->x);                    // X = X1*t6
    fp2sqr1271(t2, t6);                            // t6 = t2^2
    fp2sqr1271(t3, t2);                            // t2 = t3^2
    fp2sqr1271(t4, t3);                            // t3 = t4^2
    fp2mul1271(t2, (felm_t*)&cphi8, t1);           // t1 = cphi8*t2
    fp2mul1271(t3, (felm_t*)&cphi9, t5);           // t5 = cphi9*t3
    fp2add1271(t1, t6, t1);                        // t1 = t1+t6
    fp2mul1271(t2, (felm_t*)&cphi6, t2);           // t2 = cphi6*t2
    fp2mul1271(t3, (felm_t*)&cphi7, t3);           // t3 = cphi7*t3
    fp2add1271(t1, t5, t1);                        // t1 = t1+t5
    fp2add1271(t2, t3, t2);                        // t2 = t2+t3
    fp2mul1271(t1, P->y, t1);                      // t1 = Y1*t1
    fp2add1271(t6, t2, P->y);                      // Y = t6+t2
    fp2mul1271(P->x, t1, P->x);                    // X = X*t1
    fp2mul1271(P->y, (felm_t*)&cphi5, P->y);       // Y = cphi5*Y
    fpneg1271(P->x[1]);                            // Xfinal = X^p
    fp2mul1271(P->y, P->z, P->y);                  // Y = Y*Z1
    fp2mul1271(t0, t1, P->z);                      // Z = t0*t1
    fp2mul1271(P->y, t0, P->y);                    // Y = Y*t0
    fpneg1271(P->z[1]);                            // Zfinal = Z^p
    fpneg1271(P->y[1]);                            // Yfinal = Y^p
}


static __inline void ecc_delpsidel(point_extproj_t P)
{ // Apply delta_psi_delta mapping to a point, P = delta(psi_W(delta_inv(P))), 
  // where psi_W is the endomorphism on the Weierstrass form.
  // Input: P = (X1:Y1:Z1) on Ehat in twisted Edwards coordinates
  // Output: P = (Xfinal:Yfinal:Zfinal) on Ehat in twisted Edwards coordinates
    f2elm_t t0, t1, t2; 

    fpneg1271(P->x[1]);                            // X = X1^p
    fpneg1271(P->z[1]);                            // Z = Z1^p
    fpneg1271(P->y[1]);                            // Y = Y1^p
    fp2sqr1271(P->z, t2);                          // t2 = Z1^p^2
    fp2sqr1271(P->x, t0);                          // t0 = X1^p^2
    fp2mul1271(P->x, t2, P->x);                    // X = X1^p*Z1^p^2
    fp2mul1271(t2, (felm_t*)&cpsi2, P->z);         // Z = cpsi2*Z1^p^2
    fp2mul1271(t2, (felm_t*)&cpsi3, t1);           // t1 = cpsi3*Z1^p^2
    fp2mul1271(t2, (felm_t*)&cpsi4, t2);           // t2 = cpsi4*Z1^p^2
    fp2add1271(t0, P->z, P->z);                    // Z = X1^p^2 + cpsi2*Z1^p^2
    fp2add1271(t0, t2, t2);                        // t2 = X1^p^2 + cpsi4*Z1^p^2
    fp2add1271(t0, t1, t1);                        // t1 = X1^p^2 + cpsi3*Z1^p^2
    fp2neg1271(t2);                                // t2 = -(X1^p^2 + cpsi4*Z1^p^2)
    fp2mul1271(P->z, P->y, P->z);                  // Z = Y1^p*(X1^p^2 + cpsi2*Z1^p^2)
    fp2mul1271(P->x, t2, P->x);                    // X = -X1^p*Z1^p^2*(X1^p^2 + cpsi4*Z1^p^2)
    fp2mul1271(t1, P->z, P->y);                    // Yfinal = t1*Z
    fp2mul1271(P->x, (felm_t*)&cpsi1, P->x);       // Xfinal = cpsi1*X
    fp2mul1271(P->z, t2, P->z);                    // Zfinal = Z*t2
}


void ecc_psi(point_extproj_t P)
{ // Apply psi mapping to a point, P = psi(P)
  // Input: P = (X1:Y1:Z1) on E in twisted Edwards coordinates
  // Output: P = (Xfinal,Yfinal,Zfinal,Tafinal,Tbfinal) on E, where Tfinal = Tafinal*Tbfinal,
  //         corresponding to (Xfinal:Yfinal:Zfinal:Tfinal) in extended twisted Edwards coordinates

    ecc_tau(P);                            
    ecc_delpsidel(P);                              
    ecc_tau_dual(P);                        
}


void ecc_phi(point_extproj_t P)
{ // Apply phi mapping to a point, P = phi(P)
  // Input: P = (X1:Y1:Z1) on E in twisted Edwards coordinates
  // Output: P = (Xfinal,Yfinal,Zfinal,Tafinal,Tbfinal) on E, where Tfinal = Tafinal*Tbfinal,
  //         corresponding to (Xfinal:Yfinal:Zfinal:Tfinal) in extended twisted Edwards coordinates

    ecc_tau(P);                            
    ecc_delphidel(P);                              
    ecc_tau_dual(P);  
}


static __inline void mul_truncate(uint64_t* s, uint64_t* C, uint64_t* out)       
{ // 256-bit multiplication with truncation for the scalar decomposition
  // Outputs 64-bit value "out" = (uint64_t)((s*C) >> 256).
    uint128_t tt1, tt2;
    unsigned int carry1, carry2;
    uint64_t temp;

    MUL128(s[0], C[0], tt2);   
    tt2[0] = tt2[1];
    tt2[1] = 0;
    MUL128(s[1], C[0], tt1); 
    ADD128(tt1, tt2, tt1);
    MUL128(s[0], C[1], tt2); 
    ADC128(tt1, tt2, carry1, tt1);
    tt1[0] = tt1[1];
    tt1[1] = (uint64_t)(carry1);
    MUL128(s[2], C[0], tt2); 
    ADD128(tt1, tt2, tt1);
    MUL128(s[0], C[2], tt2); 
    ADC128(tt1, tt2, carry1, tt1);
    MUL128(s[1], C[1], tt2); 
    ADC128(tt1, tt2, carry2, tt1);
    tt1[0] = tt1[1];
    tt1[1] = (uint64_t)carry1 + (uint64_t)carry2;
    MUL128(s[0], C[3], tt2); 
    ADD128(tt1, tt2, tt1);
    MUL128(s[3], C[0], tt2); 
    ADC128(tt1, tt2, carry1, tt1);
    MUL128(s[1], C[2], tt2); 
    ADC128(tt1, tt2, carry2, tt1);
    temp = (uint64_t)carry1 + (uint64_t)carry2;
    MUL128(s[2], C[1], tt2); 
    ADC128(tt1, tt2, carry2, tt1);
    tt1[0] = tt1[1];
    tt1[1] = temp + (uint64_t)carry2;
    MUL128(s[1], C[3], tt2); 
    ADD128(tt1, tt2, tt1);
    MUL128(s[3], C[1], tt2); 
    ADD128(tt1, tt2, tt1);
    MUL128(s[2], C[2], tt2); 
    ADD128(tt1, tt2, tt1);
    *out = tt1[0];
}


void ecc_precomp(point_extproj_t P, point_extproj_precomp_t *T)
{ // Generation of the precomputation table used by the variable-base scalar multiplication ecc_mul().
  // Input: P = (X1,Y1,Z1,Ta,Tb), where T1 = Ta*Tb, corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
  // Output: table T containing 8 points: P, P+phi(P), P+psi(P), P+phi(P)+psi(P), P+psi(phi(P)), P+phi(P)+psi(phi(P)), P+psi(P)+psi(phi(P)), P+phi(P)+psi(P)+psi(phi(P))
  // Precomputed points use the representation (X+Y,Y-X,2Z,2dT) corresponding to (X:Y:Z:T) in extended twisted Edwards coordinates
    point_extproj_precomp_t Q, R, S;
    point_extproj_t PP;

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


void decompose(uint64_t* k, uint64_t* scalars)
{ // Scalar decomposition for the variable-base scalar multiplication
  // Input: scalar in the range [0, 2^256-1].
  // Output: 4 64-bit sub-scalars. 
    uint64_t a1, a2, a3, a4, temp, mask;

    mul_truncate(k, ell1, &a1);
    mul_truncate(k, ell2, &a2);
    mul_truncate(k, ell3, &a3);
    mul_truncate(k, ell4, &a4);

    temp = (uint64_t)k[0] - (uint64_t)a1*b11 - (uint64_t)a2*b21 - (uint64_t)a3*b31 - (uint64_t)a4*b41 + c1;
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


void cofactor_clearing(point_extproj_t P)
{ // Co-factor clearing
  // Input: P = (X1,Y1,Z1,Ta,Tb), where T1 = Ta*Tb, corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
  // Output: P = 392*P = (Xfinal,Yfinal,Zfinal,Tafinal,Tbfinal), where Tfinal = Tafinal*Tbfinal,
  //         corresponding to (Xfinal:Yfinal:Zfinal:Tfinal) in extended twisted Edwards coordinates
    point_extproj_precomp_t Q;

    R1_to_R2(P, Q);                      // Converting from (X,Y,Z,Ta,Tb) to (X+Y,Y-X,2Z,2dT)
    eccdouble(P);                        // P = 2*P using representations (X,Y,Z,Ta,Tb) <- 2*(X,Y,Z)
    eccadd(Q, P);                        // P = P+Q using representations (X,Y,Z,Ta,Tb) <- (X,Y,Z,Ta,Tb) + (X+Y,Y-X,2Z,2dT)
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
    point_extproj_t R;
    point_extproj_precomp_t S, Table[8];
    uint64_t scalars[NWORDS64_ORDER];
    unsigned int digits[65], sign_masks[65];
    int i;

    DISABLE_CACHE;
    point_setup(P, R);                                        // Convert to representation (X,Y,1,Ta,Tb)
    decompose((uint64_t*)k, scalars);                         // Scalar decomposition

    if (ecc_point_validate(R) == false) {                     // Check if point lies on the curve
        return false;
    }

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
    eccnorm(R, Q);                                            // Conversion to affine coordinates (x,y) and modular correction. 

    return true;
}

#endif


void eccset(point_t P)
{ // Set generator  
  // Output: P = (x,y)

    fp2copy1271((felm_t*)&GENERATOR_x, P->x);    // X1
    fp2copy1271((felm_t*)&GENERATOR_y, P->y);    // Y1
}


void eccnorm(point_extproj_t P, point_t Q)
{ // Normalize a projective point (X1:Y1:Z1), including full reduction
  // Input: P = (X1:Y1:Z1) in twisted Edwards coordinates    
  // Output: Q = (X1/Z1,Y1/Z1), corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
    
    fp2inv1271(P->z);                      // Z1 = Z1^-1
    fp2mul1271(P->x, P->z, Q->x);          // X1 = X1/Z1
    fp2mul1271(P->y, P->z, Q->y);          // Y1 = Y1/Z1
    mod1271(Q->x[0]); mod1271(Q->x[1]); 
    mod1271(Q->y[0]); mod1271(Q->y[1]); 
}


void R1_to_R2(point_extproj_t P, point_extproj_precomp_t Q)
{ // Conversion from representation (X,Y,Z,Ta,Tb) to (X+Y,Y-X,2Z,2dT), where T = Ta*Tb
  // Input:  P = (X1,Y1,Z1,Ta,Tb), where T1 = Ta*Tb, corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
  // Output: Q = (X1+Y1,Y1-X1,2Z1,2dT1) corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates

    fp2add1271(P->ta, P->ta, Q->t2);                  // T = 2*Ta
    fp2add1271(P->x, P->y, Q->xy);                    // QX = X+Y
    fp2sub1271(P->y, P->x, Q->yx);                    // QY = Y-X 
    fp2mul1271(Q->t2, P->tb, Q->t2);                  // T = 2*T
    fp2add1271(P->z, P->z, Q->z2);                    // QZ = 2*Z
    fp2mul1271(Q->t2, (felm_t*)&PARAMETER_d, Q->t2);  // QT = 2d*T
}


__inline void R1_to_R3(point_extproj_t P, point_extproj_precomp_t Q)      
{ // Conversion from representation (X,Y,Z,Ta,Tb) to (X+Y,Y-X,Z,T), where T = Ta*Tb 
  // Input:  P = (X1,Y1,Z1,Ta,Tb), where T1 = Ta*Tb, corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
  // Output: Q = (X1+Y1,Y1-X1,Z1,T1) corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates 
    
    fp2add1271(P->x, P->y, Q->xy);         // XQ = (X1+Y1) 
    fp2sub1271(P->y, P->x, Q->yx);         // YQ = (Y1-X1) 
    fp2mul1271(P->ta, P->tb, Q->t2);       // TQ = T1
    fp2copy1271(P->z, Q->z2);              // ZQ = Z1 
}


void R2_to_R4(point_extproj_precomp_t P, point_extproj_t Q)      
{ // Conversion from representation (X+Y,Y-X,2Z,2dT) to (2X,2Y,2Z,2dT) 
  // Input:  P = (X1+Y1,Y1-X1,2Z1,2dT1) corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
  // Output: Q = (2X1,2Y1,2Z1) corresponding to (X1:Y1:Z1) in twisted Edwards coordinates 
    
    fp2sub1271(P->xy, P->yx, Q->x);        // XQ = 2*X1
    fp2add1271(P->xy, P->yx, Q->y);        // YQ = 2*Y1
    fp2copy1271(P->z2, Q->z);              // ZQ = 2*Z1
}


void eccdouble(point_extproj_t P)
{ // Point doubling 2P
  // Input: P = (X1:Y1:Z1) in twisted Edwards coordinates
  // Output: 2P = (Xfinal,Yfinal,Zfinal,Tafinal,Tbfinal), where Tfinal = Tafinal*Tbfinal,
  //         corresponding to (Xfinal:Yfinal:Zfinal:Tfinal) in extended twisted Edwards coordinates
    f2elm_t t1, t2;  

    fp2sqr1271(P->x, t1);                  // t1 = X1^2
    fp2sqr1271(P->y, t2);                  // t2 = Y1^2
    fp2add1271(P->x, P->y, P->x);          // t3 = X1+Y1
    fp2add1271(t1, t2, P->tb);             // Tbfinal = X1^2+Y1^2      
    fp2sub1271(t2, t1, t1);                // t1 = Y1^2-X1^2  
    fp2sqr1271(P->z, t2);                  // t2 = Z1^2       
    fp2sqr1271(P->x, P->ta);               // Ta = (X1+Y1)^2 
    fp2add1271(t2, t2, t2);                // t2 = 2Z1^2 
    fp2sub1271(P->ta, P->tb, P->ta);       // Tafinal = 2X1*Y1 = (X1+Y1)^2-(X1^2+Y1^2) 
    fp2sub1271(t2, t1, t2);                // t2 = 2Z1^2-(Y1^2-X1^2) 
    fp2mul1271(t1, P->tb, P->y);           // Yfinal = (X1^2+Y1^2)(Y1^2-X1^2)  
    fp2mul1271(t2, P->ta, P->x);           // Xfinal = 2X1*Y1*[2Z1^2-(Y1^2-X1^2)]
    fp2mul1271(t1, t2, P->z);              // Zfinal = (Y1^2-X1^2)[2Z1^2-(Y1^2-X1^2)]
}


__inline void eccadd_core(point_extproj_precomp_t P, point_extproj_precomp_t Q, point_extproj_t R)      
{ // Basic point addition R = P+Q or R = P+P
  // Inputs: P = (X1+Y1,Y1-X1,2Z1,2dT1) corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
  //         Q = (X2+Y2,Y2-X2,Z2,T2) corresponding to (X2:Y2:Z2:T2) in extended twisted Edwards coordinates    
  // Output: R = (Xfinal,Yfinal,Zfinal,Tafinal,Tbfinal), where Tfinal = Tafinal*Tbfinal,
  //         corresponding to (Xfinal:Yfinal:Zfinal:Tfinal) in extended twisted Edwards coordinates
    f2elm_t t1, t2; 
          
    fp2mul1271(P->t2, Q->t2, R->z);        // Z = 2dT1*T2 
    fp2mul1271(P->z2, Q->z2, t1);          // t1 = 2Z1*Z2  
    fp2mul1271(P->xy, Q->xy, R->x);        // X = (X1+Y1)(X2+Y2) 
    fp2mul1271(P->yx, Q->yx, R->y);        // Y = (Y1-X1)(Y2-X2) 
    fp2sub1271(t1, R->z, t2);              // t2 = theta
    fp2add1271(t1, R->z, t1);              // t1 = alpha
    fp2sub1271(R->x, R->y, R->tb);         // Tbfinal = beta
    fp2add1271(R->x, R->y, R->ta);         // Tafinal = omega
    fp2mul1271(R->tb, t2, R->x);           // Xfinal = beta*theta
    fp2mul1271(t1, t2, R->z);              // Zfinal = theta*alpha
    fp2mul1271(R->ta, t1, R->y);           // Yfinal = alpha*omega
}


void eccadd(point_extproj_precomp_t Q, point_extproj_t P)      
{ // Complete point addition P = P+Q or P = P+P
  // Inputs: P = (X1,Y1,Z1,Ta,Tb), where T1 = Ta*Tb, corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
  //         Q = (X2+Y2,Y2-X2,2Z2,2dT2) corresponding to (X2:Y2:Z2:T2) in extended twisted Edwards coordinates   
  // Output: P = (Xfinal,Yfinal,Zfinal,Tafinal,Tbfinal), where Tfinal = Tafinal*Tbfinal, 
  //         corresponding to (Xfinal:Yfinal:Zfinal:Tfinal) in extended twisted Edwards coordinates
    point_extproj_precomp_t R;
    
    R1_to_R3(P, R);                        // R = (X1+Y1,Y1-Z1,Z1,T1)
    eccadd_core(Q, R, P);                  // P = (X2+Y2,Y2-X2,2Z2,2dT2) + (X1+Y1,Y1-Z1,Z1,T1)
}


void point_setup(point_t P, point_extproj_t Q)
{ // Point conversion to representation (X,Y,Z,Ta,Tb) 
  // Input: P = (x,y) in affine coordinates
  // Output: P = (X,Y,1,Ta,Tb), where Ta=X, Tb=Y and T=Ta*Tb, corresponding to (X:Y:Z:T) in extended twisted Edwards coordinates

    fp2copy1271(P->x, Q->x);
    fp2copy1271(P->y, Q->y);
    fp2copy1271(Q->x, Q->ta);              // Ta = X1
    fp2copy1271(Q->y, Q->tb);              // Tb = Y1
    fp2zero1271(Q->z); Q->z[0][0]=1;       // Z1 = 1
}


bool ecc_point_validate(point_extproj_t P)
{ // Point validation: check if point lies on the curve
  // Input: P = (x,y) in affine coordinates, where x, y in [0, 2^127-1]. 
  // Output: TRUE (1) if point lies on the curve E: -x^2+y^2-1-dx^2*y^2 = 0, FALSE (0) otherwise. 
  // SECURITY NOTE: this function does not run in constant time (input point P is assumed to be public).
    f2elm_t t1, t2, t3;

    fp2sqr1271(P->y, t1);
    fp2sqr1271(P->x, t2);
    fp2sub1271(t1, t2, t3);                     // -x^2 + y^2 
    fp2mul1271(t1, t2, t1);                     // x^2*y^2
    fp2mul1271((felm_t*)&PARAMETER_d, t1, t2);  // dx^2*y^2
    fp2zero1271(t1);  t1[0][0] = 1;             // t1 = 1
    fp2add1271(t2, t1, t2);                     // 1 + dx^2*y^2
    fp2sub1271(t3, t2, t1);                     // -x^2 + y^2 - 1 - dx^2*y^2 
    
    return ((is_digit_zero_ct(t1[0][0] | t1[0][1] | t1[0][2] | t1[0][3]) || is_digit_zero_ct((t1[0][0]+1) | (t1[0][1]+1) | (t1[0][2]+1) | (t1[0][3]+1))) &&
            (is_digit_zero_ct(t1[1][0] | t1[1][1] | t1[1][2] | t1[1][3]) || is_digit_zero_ct((t1[1][0]+1) | (t1[1][1]+1) | (t1[1][2]+1) | (t1[1][3]+1))));
}


__inline void R5_to_R1(point_precomp_t P, point_extproj_t Q)
{ // Conversion from representation (x+y,y-x,2dt) to (X,Y,Z,Ta,Tb) 
  // Input:  P = (x1+y1,y1-x1,2dt1) corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates, where Z1=1
  // Output: Q = (x1,y1,z1,x1,y1), where z1=1, corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates 

    fp2sub1271(P->xy, P->yx, Q->x);        // 2*x1
    fp2add1271(P->xy, P->yx, Q->y);        // 2*y1
    fp2div1271(Q->x);                      // XQ = x1
    fp2div1271(Q->y);                      // YQ = y1 
    fp2zero1271(Q->z); Q->z[0][0] = 1;     // ZQ = 1
    fp2copy1271(Q->x, Q->ta);              // TaQ = x1
    fp2copy1271(Q->y, Q->tb);              // TbQ = y1
}

static __inline void eccmadd(point_precomp_t Q, point_extproj_t P)
{ // Mixed point addition P = P+Q or P = P+P
  // Inputs: P = (X1,Y1,Z1,Ta,Tb), where T1 = Ta*Tb, corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
  //         Q = (x2+y2,y2-x2,2dt2) corresponding to (X2:Y2:Z2:T2) in extended twisted Edwards coordinates, where Z2=1  
  // Output: P = (Xfinal,Yfinal,Zfinal,Tafinal,Tbfinal), where Tfinal = Tafinal*Tbfinal, 
  //         corresponding to (Xfinal:Yfinal:Zfinal:Tfinal) in extended twisted Edwards coordinates
    f2elm_t t1, t2;
    
    fp2mul1271(P->ta, P->tb, P->ta);        // Ta = T1
    fp2add1271(P->z, P->z, t1);             // t1 = 2Z1        
    fp2mul1271(P->ta, Q->t2, P->ta);        // Ta = 2dT1*t2 
    fp2add1271(P->x, P->y, P->z);           // Z = (X1+Y1) 
    fp2sub1271(P->y, P->x, P->tb);          // Tb = (Y1-X1)
    fp2sub1271(t1, P->ta, t2);              // t2 = theta
    fp2add1271(t1, P->ta, t1);              // t1 = alpha
    fp2mul1271(Q->xy, P->z, P->ta);         // Ta = (X1+Y1)(x2+y2)
    fp2mul1271(Q->yx, P->tb, P->x);         // X = (Y1-X1)(y2-x2)
    fp2mul1271(t1, t2, P->z);               // Zfinal = theta*alpha
    fp2sub1271(P->ta, P->x, P->tb);         // Tbfinal = beta
    fp2add1271(P->ta, P->x, P->ta);         // Tafinal = omega
    fp2mul1271(P->tb, t2, P->x);            // Xfinal = beta*theta
    fp2mul1271(P->ta, t1, P->y);            // Yfinal = alpha*omega
}


bool ecc_mul_fixed(digit_t* k, point_t Q)
{ // Fixed-base scalar multiplication Q = k*G, where G is the generator. FIXED_BASE_TABLE stores v*2^(w-1) = 80 multiples of G.
  // Inputs: scalar "k" in [0, 2^256-1].
  // Output: Q = k*G in affine coordinates (x,y).
  // The function is based on the modified LSB-set comb method, which converts the scalar to an odd signed representation
  // with (bitlength(order)+w*v) digits.
    unsigned int j, w = W_FIXEDBASE, v = V_FIXEDBASE, d = D_FIXEDBASE, e = E_FIXEDBASE;
    unsigned int digit = 0, digits[NBITS_ORDER_PLUS_ONE + (W_FIXEDBASE*V_FIXEDBASE) - 1] = {0};
    digit_t temp[NWORDS_ORDER];
    point_extproj_t R;
    point_precomp_t S;
    int i, ii;
    
    DISABLE_CACHE;
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
    table_lookup_fixed_base(((point_precomp_t*)&FIXED_BASE_TABLE)+(v-1)*(1 << (w-1)), S, digit, digits[d-1]);
    R5_to_R1(S, R);                                             // Converting to representation (X:Y:1:Ta:Tb)

    for (j = 0; j < (v-1); j++)
    {
        digit = digits[w*d-(j+1)*e-1];
        for (i = (int)((w-1)*d-(j+1)*e-1); i >= (int)(2*d-(j+1)*e-1); i = i-d)
        {
            digit = 2*digit + digits[i];
        }
        // Extract point in (x+y,y-x,2dt) representation
        table_lookup_fixed_base(((point_precomp_t*)&FIXED_BASE_TABLE)+(v-j-2)*(1 << (w-1)), S, digit, digits[d-(j+1)*e-1]);
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
            table_lookup_fixed_base(((point_precomp_t*)&FIXED_BASE_TABLE)+(v-j-1)*(1 << (w-1)), S, digit, digits[d-j*e+ii-e]);
            eccmadd(S, R);                                      // R = R+S using representations (X,Y,Z,Ta,Tb) <- (X,Y,Z,Ta,Tb) + (x+y,y-x,2dt)
        }
    }
    eccnorm(R, Q);                                              // Conversion to affine coordinates (x,y) and modular correction.

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


static __inline void eccneg_extproj_precomp(point_extproj_precomp_t P, point_extproj_precomp_t Q)
{ // Point negation
  // Input : point P in coordinates (X+Y,Y-X,2Z,2dT)
  // Output: point Q = -P = (Y-X,X+Y,2Z,-2dT)
    fp2copy1271(P->t2, Q->t2);
    fp2copy1271(P->xy, Q->yx);
    fp2copy1271(P->yx, Q->xy);
    fp2copy1271(P->z2, Q->z2);
    fp2neg1271(Q->t2);
}


static __inline void eccneg_precomp(point_precomp_t P, point_precomp_t Q)
{ // Point negation
  // Input : point P in coordinates (x+y,y-x,2dt)
  // Output: point Q = -P = (y-x,x+y,-2dt)
    fp2copy1271(P->t2, Q->t2);
    fp2copy1271(P->xy, Q->yx);
    fp2copy1271(P->yx, Q->xy);
    fp2neg1271(Q->t2);
}


bool ecc_mul_double(digit_t* k, point_t Q, digit_t* l, point_t R)
{ // Double scalar multiplication R = k*G + l*Q, where the G is the generator. Uses DOUBLE_SCALAR_TABLE, which contains multiples of G, Phi(G), Psi(G) and Phi(Psi(G)).
  // Inputs: point Q in affine coordinates,
  //         scalars "k" and "l" in [0, 2^256-1].
  // Output: R = k*G + l*Q in affine coordinates (x,y).
  // The function uses wNAF with interleaving.

  // SECURITY NOTE: this function is intended for a non-constant-time operation such as signature verification.

#if (USE_ENDO == true)
    unsigned int position;
    int i, digits_k1[65] = {0}, digits_k2[65] = {0}, digits_k3[65] = {0}, digits_k4[65] = {0};
    int digits_l1[65] = {0}, digits_l2[65] = {0}, digits_l3[65] = {0}, digits_l4[65] = {0};
    point_precomp_t V;
    point_extproj_t Q1, Q2, Q3, Q4, T;
    point_extproj_precomp_t U, Q_table1[NPOINTS_DOUBLEMUL_WQ], Q_table2[NPOINTS_DOUBLEMUL_WQ], Q_table3[NPOINTS_DOUBLEMUL_WQ], Q_table4[NPOINTS_DOUBLEMUL_WQ];
    uint64_t k_scalars[4], l_scalars[4];
     
    ENABLE_CACHE;
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

    fp2zero1271(T->x);                                         // Initialize T as the neutral point (0:1:1)
    fp2zero1271(T->y); T->y[0][0] = 1;
    fp2zero1271(T->z); T->z[0][0] = 1;

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
            eccneg_precomp(((point_precomp_t*)&DOUBLE_SCALAR_TABLE)[position], V);    // Load and negate V = (X_V,Y_V,Z_V,Td_V) <- -(x+y,y-x,2dt) from a point in the precomputed table 
            eccmadd(V, T);                                                            // T = T+V = (X_T,Y_T,Z_T,Ta_T,Tb_T) = (X_T,Y_T,Z_T,Ta_T,Tb_T) + (X_V,Y_V,Z_V,Td_V) 
        } else if (digits_k1[i] > 0) {
            position = (digits_k1[i])/2;                                              // Take V = (X_V,Y_V,Z_V,Td_V) <- (x+y,y-x,2dt) from a point in the precomputed table
            eccmadd(((point_precomp_t*)&DOUBLE_SCALAR_TABLE)[position], T);           // T = T+V = (X_T,Y_T,Z_T,Ta_T,Tb_T) = (X_T,Y_T,Z_T,Ta_T,Tb_T) + (X_V,Y_V,Z_V,Td_V) 
        }
        if (digits_k2[i] < 0) {
            position = (-digits_k2[i])/2;
            eccneg_precomp(((point_precomp_t*)&DOUBLE_SCALAR_TABLE)[NPOINTS_DOUBLEMUL_WP+position], V);
            eccmadd(V, T);
        } else if (digits_k2[i] > 0) {
            position = (digits_k2[i])/2;
            eccmadd(((point_precomp_t*)&DOUBLE_SCALAR_TABLE)[NPOINTS_DOUBLEMUL_WP+position], T);
        }
        if (digits_k3[i] < 0) {
            position = (-digits_k3[i])/2;
            eccneg_precomp(((point_precomp_t*)&DOUBLE_SCALAR_TABLE)[2*NPOINTS_DOUBLEMUL_WP+position], V);
            eccmadd(V, T);
        } else if (digits_k3[i] > 0) {
            position = (digits_k3[i])/2;
            eccmadd(((point_precomp_t*)&DOUBLE_SCALAR_TABLE)[2*NPOINTS_DOUBLEMUL_WP+position], T);
        }
        if (digits_k4[i] < 0) {
            position = (-digits_k4[i])/2;
            eccneg_precomp(((point_precomp_t*)&DOUBLE_SCALAR_TABLE)[3*NPOINTS_DOUBLEMUL_WP+position], V);
            eccmadd(V, T);
        } else if (digits_k4[i] > 0) {
            position = (digits_k4[i])/2;
            eccmadd(((point_precomp_t*)&DOUBLE_SCALAR_TABLE)[3*NPOINTS_DOUBLEMUL_WP+position], T);
        }
    }

#else
    point_t A;
    point_extproj_t T;
    point_extproj_precomp_t S;
     
    ENABLE_CACHE;
    if (ecc_mul(Q, l, A, false) == false) {
        return false;
    }
    point_setup(A, T);
    R1_to_R2(T, S);

    ecc_mul_fixed(k, A);
    point_setup(A, T);
    eccadd(S, T);
#endif
    eccnorm(T, R);                                             // Output R = (x,y)
    DISABLE_CACHE;

    return true;
}


void ecc_precomp_double(point_extproj_t P, point_extproj_precomp_t* Table, unsigned int npoints)
{ // Generation of the precomputation table used internally by the double scalar multiplication function ecc_mul_double().  
  // Inputs: point P in representation (X,Y,Z,Ta,Tb),
  //         Table with storage for npoints, 
  //         number of points "npoints".
  // Output: Table containing multiples of the base point P using representation (X+Y,Y-X,2Z,2dT).
    point_extproj_t Q;
    point_extproj_precomp_t PP;
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
