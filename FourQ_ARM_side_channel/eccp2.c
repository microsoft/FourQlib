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
#include "../random/random.h"
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


void felmt_randomize_point(point_extedwards_t P, felm_t random)
{ // Randomization of point coordinates using a random field element
  // Input:  P = (X,Y,Z,T) in extended twisted Edwards coordinates.
  // Output: P = (X,Y,Z,T) in extended twisted Edwards coordinates, 
  //         where X = random*(Xa+Xb*i), Y = random*(Ya+Yb*i), Z = random*(Za+Zb*i), T = random*(Ta+Tb*i). 

    fpmul1271(P->x[0], random, P->x[0]);
    fpmul1271(P->x[1], random, P->x[1]);
    fpmul1271(P->y[0], random, P->y[0]);
    fpmul1271(P->y[1], random, P->y[1]);
    fpmul1271(P->z[0], random, P->z[0]);
    fpmul1271(P->z[1], random, P->z[1]);
    fpmul1271(P->t[0], random, P->t[0]);
    fpmul1271(P->t[1], random, P->t[1]);
}


void randomize_table(point_extedwards_t* Table, felm_t random)
{ // Randomization of all the point coordinates in the precomputed table using a random field element    ///// NOTE: SHOULD WE RANDOMIZE POINTS WITH AN ELEMENT IN GF(P^2)
    unsigned int i;

    for (i = 0; i < 16; i++) {
        felmt_randomize_point(Table[i], random);
    }
}


void random_felmt(felm_t random)
{ // Generate random field element in [0, 2^127-1]
    unsigned char* rand_byte = (unsigned char*)&random[0];

    RandomBytesFunction(rand_byte, 16);
    rand_byte[15] &= 0x7F;
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

// Fixed integer constants for scalar randomization 
static uint64_t p11 = {0x190BE2D7F2E68811};   
static uint64_t p12 = {0x2E5EBE12E05824E0};                                              
static uint64_t p13 = {0x27C2D7D95E7F1AEB};                                              
static uint64_t p14 = {0x02684DFF36C48F4F};   
static uint64_t p21 = {0x36553EC277E34AE5};
static uint64_t p22 = {0x2E5EBE12E05824DF};
static uint64_t p23 = {0x27C2D7D95E7F1AEC};
static uint64_t p24 = {0x284413BBC495F71F};
static uint64_t p31 = {0x39BE8F1CF6A62CA9};
static uint64_t p32 = {0x1DBEF6CB01B6D191};
static uint64_t p33 = {0x00B81641C21F595B};
static uint64_t p34 = {0x04B749AA70397695};
static uint64_t p41 = {0x3F5C3DEA5883EC7B};
static uint64_t p42 = {0x1AFAD5B01E2DF73F};
static uint64_t p43 = {0x2F05380B4B471DFB};
static uint64_t p44 = {0x1FF4A0223DDC10CE};


/***********************************************/
/**********  CURVE/SCALAR FUNCTIONS  ***********/

static __inline void ecc_tau(point_extedwards_t P)
{ // Apply tau mapping to a point, P = tau(P)
  // Input: P = (X1:Y1:Z1) on E in twisted Edwards coordinates
  // Output: P = (Xfinal:Yfinal:Zfinal) on Ehat in twisted Edwards coordinates
    f2elm_t t0; 

    fp2sqr1271(P->x, t0);                     // t0 = X1^2
    fp2sqr1271(P->y, P->t);                   // T = Y1^2
    fp2mul1271(P->x, P->y, P->x);             // X = X1*Y1
    fp2sqr1271(P->z, P->y);                   // Y = Z1^2
    fp2add1271(t0, P->t, P->z);               // Z = X1^2+Y1^2
    fp2add1271(P->y, P->y, P->y);             // Y = 2*Z1^2
    fp2sub1271(t0, P->t, t0);                 // t0 = X1^2-Y1^2
    fp2neg1271(P->y);                         // Y = -2*Z1^2
    fp2mul1271(P->x, t0, P->x);               // X = X1*Y1*(X1^2-Y1^2)
    fp2sub1271(P->y, t0, P->y);               // Y = -2*Z1^2-(X1^2-Y1^2)
    fp2mul1271(P->x, (felm_t*)&ctau1, P->x);  // Xfinal = X*ctau1
    fp2mul1271(P->y, P->z, P->y);             // Yfinal = Y*Z
    fp2mul1271(P->z, t0, P->z);               // Zfinal = t0*Z
}


static __inline void ecc_tau_dual(point_extedwards_t P)
{ // Apply tau_dual mapping to a point, P = tau_dual(P)
  // Input: P = (X1:Y1:Z1:T1) on Ehat in extended twisted Edwards coordinates
  // Output: P = (Xfinal,Yfinal,Zfinal,Tfinal) on E in extended twisted Edwards coordinates
    f2elm_t t0, t1;

    fp2sqr1271(P->x, t0);                          // t0 = X1^2
    fp2sqr1271(P->z, P->t);                        // T = Z1^2
    fp2sqr1271(P->y, t1);                          // t1 = Y1^2
    fp2add1271(P->t, P->t, P->z);                  // Z = 2*Z1^2
    fp2sub1271(t1, t0, P->t);                      // T = Y1^2-X1^2
    fp2add1271(t0, t1, t0);                        // t0 = X1^2+Y1^2
    fp2mul1271(P->x, P->y, P->x);                  // X = X1*Y1
    fp2sub1271(P->z, P->t, P->z);                  // Z = 2*Z1^2-(Y1^2-X1^2)
    fp2mul1271(P->x, (felm_t*)&ctaudual1, t1);     // t1 = ctaudual1*X1
    fp2mul1271(P->z, P->t, P->y);                  // Yfinal = Z*T
    fp2mul1271(t0, t1, P->x);                      // Xfinal = t0*t1
    fp2mul1271(P->z, t0, P->z);                    // Zfinal = Z*t0
    fp2mul1271(P->t, t1, P->t);                    // Tfinal = T*t1
}


static __inline void ecc_delphidel(point_extedwards_t P)
{ // Apply delta_phi_delta mapping to a point, P = delta(phi_W(delta_inv(P))), 
  // where phi_W is the endomorphism on the Weierstrass form.
  // Input: P = (X1:Y1:Z1) on Ehat in twisted Edwards coordinates
  // Output: P = (Xfinal:Yfinal:Zfinal) on Ehat in twisted Edwards coordinates
    f2elm_t t0, t1, t2, t3, t4, t5; 
    
    fp2sqr1271(P->y, t2);                          // t2 = Y1^2
    fp2sqr1271(P->z, t4);                          // t4 = Z1^2
    fp2mul1271(t4, (felm_t*)&cphi4, t0);           // t0 = cphi4*t4
    fp2mul1271(P->y, P->z, t3);                    // t3 = Y1*Z1
    fp2add1271(t0, t2, t0);                        // t0 = t0+t2
    fp2mul1271(t3, (felm_t*)&cphi3, t1);           // t1 = cphi3*t3
    fp2sub1271(t0, t1, t5);                        // t5 = t0-t1
    fp2add1271(t0, t1, t0);                        // t0 = t0+t1
    fp2mul1271(t0, P->z, t0);                      // t0 = t0*Z1
    fp2mul1271(t3, (felm_t*)&cphi1, t1);           // t1 = cphi1*t3
    fp2mul1271(t0, t5, t0);                        // t0 = t0*t5
    fp2mul1271(t4, (felm_t*)&cphi2, t5);           // t5 = cphi2*t4
    fp2add1271(t2, t5, t5);                        // t5 = t2+t5
    fp2sub1271(t1, t5, P->t);                      // T = t1-t5
    fp2add1271(t1, t5, t1);                        // t1 = t1+t5
    fp2mul1271(P->t, t1, P->t);                    // T = t1*P->t
    fp2mul1271(P->t, (felm_t*)&cphi0, P->t);       // T = cphi0*T
    fp2mul1271(P->x, P->t, P->x);                  // X = X1*T
    fp2sqr1271(t2, P->t);                          // T = t2^2
    fp2sqr1271(t3, t2);                            // t2 = t3^2
    fp2sqr1271(t4, t3);                            // t3 = t4^2
    fp2mul1271(t2, (felm_t*)&cphi8, t1);           // t1 = cphi8*t2
    fp2mul1271(t3, (felm_t*)&cphi9, t5);           // t5 = cphi9*t3
    fp2add1271(t1, P->t, t1);                      // t1 = t1+T
    fp2mul1271(t2, (felm_t*)&cphi6, t2);           // t2 = cphi6*t2
    fp2mul1271(t3, (felm_t*)&cphi7, t3);           // t3 = cphi7*t3
    fp2add1271(t1, t5, t1);                        // t1 = t1+t5
    fp2add1271(t2, t3, t2);                        // t2 = t2+t3
    fp2mul1271(t1, P->y, t1);                      // t1 = Y1*t1
    fp2add1271(P->t, t2, P->y);                    // Y = T+t2
    fp2mul1271(P->x, t1, P->x);                    // X = X*t1
    fp2mul1271(P->y, (felm_t*)&cphi5, P->y);       // Y = cphi5*Y
    fpneg1271(P->x[1]);                            // Xfinal = X^p
    fp2mul1271(P->y, P->z, P->y);                  // Y = Y*Z1
    fp2mul1271(t0, t1, P->z);                      // Z = t0*t1
    fp2mul1271(P->y, t0, P->y);                    // Y = Y*t0
    fpneg1271(P->z[1]);                            // Zfinal = Z^p
    fpneg1271(P->y[1]);                            // Yfinal = Y^p
}


static __inline void ecc_delpsidel(point_extedwards_t P)
{ // Apply delta_psi_delta mapping to a point, P = delta(psi_W(delta_inv(P))), 
  // where psi_W is the endomorphism on the Weierstrass form.
  // Input: P = (X1:Y1:Z1) on Ehat in twisted Edwards coordinates
  // Output: P = (Xfinal:Yfinal:Zfinal) on Ehat in twisted Edwards coordinates
    f2elm_t t0, t1; 

    fpneg1271(P->x[1]);                            // X = X1^p
    fpneg1271(P->z[1]);                            // Z = Z1^p
    fpneg1271(P->y[1]);                            // Y = Y1^p
    fp2sqr1271(P->z, P->t);                        // T = Z1^p^2
    fp2sqr1271(P->x, t0);                          // t0 = X1^p^2
    fp2mul1271(P->x, P->t, P->x);                  // X = X1^p*Z1^p^2
    fp2mul1271(P->t, (felm_t*)&cpsi2, P->z);       // Z = cpsi2*Z1^p^2
    fp2mul1271(P->t, (felm_t*)&cpsi3, t1);         // t1 = cpsi3*Z1^p^2
    fp2mul1271(P->t, (felm_t*)&cpsi4, P->t);       // T = cpsi4*Z1^p^2
    fp2add1271(t0, P->z, P->z);                    // Z = X1^p^2 + cpsi2*Z1^p^2
    fp2add1271(t0, P->t, P->t);                    // T = X1^p^2 + cpsi4*Z1^p^2
    fp2add1271(t0, t1, t1);                        // t1 = X1^p^2 + cpsi3*Z1^p^2
    fp2neg1271(P->t);                              // T = -(X1^p^2 + cpsi4*Z1^p^2)
    fp2mul1271(P->z, P->y, P->z);                  // Z = Y1^p*(X1^p^2 + cpsi2*Z1^p^2)
    fp2mul1271(P->x, P->t, P->x);                  // X = -X1^p*Z1^p^2*(X1^p^2 + cpsi4*Z1^p^2)
    fp2mul1271(t1, P->z, P->y);                    // Yfinal = t1*Z
    fp2mul1271(P->x, (felm_t*)&cpsi1, P->x);       // Xfinal = cpsi1*X
    fp2mul1271(P->z, P->t, P->z);                  // Zfinal = Z*T
}


void ecc_psi(point_extedwards_t P)
{ // Apply psi mapping to a point, P = psi(P)
  // Input: P = (X1:Y1:Z1) on E in twisted Edwards coordinates
  // Output: P = (Xfinal,Yfinal,Zfinal,Tfinal) on E in extended twisted Edwards coordinates

    ecc_tau(P);                            
    ecc_delpsidel(P);                      		
    ecc_tau_dual(P);                        
}


void ecc_phi(point_extedwards_t P)
{ // Apply phi mapping to a point, P = phi(P)
  // Input: P = (X1:Y1:Z1) on E in twisted Edwards coordinates
  // Output: P = (Xfinal,Yfinal,Zfinal,Tfinal) on E in extended twisted Edwards coordinates

    ecc_tau(P);                            
    ecc_delphidel(P);                      		
    ecc_tau_dual(P);  
}


static __inline void mul_truncate(uint64_t* s, uint64_t* C, uint64_t* out)       
{ // 256-bit multiplication with truncation for the scalar decomposition
  // Outputs 64-bit value "out" = (uint64_t)((s * C) >> 256).
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


void ecc_precomp(point_extedwards_t P, point_extedwards_t R, point_extedwards_t *Table)
{ // Generation of the precomputation table used by the variable-base scalar multiplication ecc_mul_SCA_secure().
  // Input: P = (XP,YP,ZP,TP) in extended twisted Edwards coordinates
  //        R = (XR,YR,ZR,TR) in extended twisted Edwards coordinates. R is used for point blinding
  // Output: Table containing 16 points: P, P+phi(P), P+psi(P), P+phi(P)+psi(P), P+psi(phi(P)), P+phi(P)+psi(phi(P)), P+psi(P)+psi(phi(P)), P+phi(P)+psi(P)+psi(phi(P))
  // Precomputed points use the representation (X+Y,Y-X,2Z,2dT) in extended twisted Edwards coordinates
    point_extedwards_t S, T, U;
    unsigned int i;

    // Negating the random point R, Table[0] = -R
    ecccopy(R, Table[0]);
    fp2neg1271(Table[0]->x);
    fp2neg1271(Table[0]->t);

    // Generating S = phi(P) = (XS,YS,ZS,TS)
    ecccopy(P, S);
    ecc_phi(S);

    // Generating T = psi(T) = (XT,YT,ZT,TT)
    ecccopy(P, T);
    ecc_psi(T); 

    // Generating U = psi(phi(P)) = (XU,YU,ZU,TU)
    ecccopy(S, U);
    ecc_psi(U); 

    eccadd(Table[0], P, Table[1]);         // Table[1] = -R+P
    eccadd(Table[0], S, Table[2]);         // Table[2] = -R+S
    eccadd(Table[1], S, Table[3]);         // Table[3] = -R+P+S
    eccadd(Table[0], T, Table[4]);         // Table[4] = -R+T
    eccadd(Table[4], P, Table[5]);         // Table[5] = -R+P+T
    eccadd(Table[4], S, Table[6]);         // Table[6] = -R+S+T
    eccadd(Table[6], P, Table[7]);         // Table[7] = -R+P+S+T
    eccadd(Table[0], U, Table[8]);         // Table[8] = -R+U
    eccadd(Table[1], U, Table[9]);         // Table[9] = -R+P+U
    eccadd(Table[2], U, Table[10]);        // Table[10] = -R+S+U
    eccadd(Table[3], U, Table[11]);        // Table[11] = -R+P+S+U
    eccadd(Table[4], U, Table[12]);        // Table[12] = -R+T+U
    eccadd(Table[5], U, Table[13]);        // Table[13] = -R+P+T+U
    eccadd(Table[6], U, Table[14]);        // Table[14] = -R+S+T+U
    eccadd(Table[7], U, Table[15]);        // Table[15] = -R+P+S+T+U
    
    for (i = 0; i < 16; i++) {
        fp2add1271(Table[i]->x, Table[i]->y, S->x);                 // Converting to coordinates (X+Y,Y-X,2*Z,2*d*T)
        fp2sub1271(Table[i]->x, Table[i]->y, Table[i]->y);
        fp2copy1271(S->x, Table[i]->x);
        fp2add1271(Table[i]->z, Table[i]->z, Table[i]->z);
        fp2mul1271(Table[i]->t, (felm_t*)&PARAMETER_d, Table[i]->t);    
        fp2add1271(Table[i]->t, Table[i]->t, Table[i]->t);
    }
}


void decompose(uint64_t* k, uint64_t* scalars)
{ // Scalar decomposition for the variable-base scalar multiplication
  // Input: scalar in the range [0, 2^256-1].
  // Output: 4 64-bit sub-scalars. 
    uint64_t a1, a2, a3, a4;

    mul_truncate(k, ell1, &a1);
    mul_truncate(k, ell2, &a2);
    mul_truncate(k, ell3, &a3);
    mul_truncate(k, ell4, &a4);
    
    scalars[0] = (uint64_t)k[0] - (uint64_t)a1*b11 - (uint64_t)a2*b21 - (uint64_t)a3*b31 - (uint64_t)a4*b41 + c1;
    scalars[1] =                  (uint64_t)a1*b12 + (uint64_t)a2     - (uint64_t)a3*b32 - (uint64_t)a4*b42 + c2;
    scalars[2] =                  (uint64_t)a3*b33 - (uint64_t)a1*b13 - (uint64_t)a2     + (uint64_t)a4*b43 + c3;
    scalars[3] =                  (uint64_t)a1*b14 - (uint64_t)a2*b24 - (uint64_t)a3*b34 + (uint64_t)a4*b44 + c4;
}


void randomize(uint64_t* scalars, unsigned char* r, uint128_t* random_scalars)
{ // Scalar randomization for the variable-base scalar multiplication
  // Input: 4 64-bit sub-scalars, and 4 16-bit random values.
  // Output: 4 80-bit sub-scalars. 
    uint128_t tt0;
    uint64_t r0 = (uint64_t)r[0], r1 = (uint64_t)r[1];
    uint64_t r2 = (uint64_t)r[2], r3 = (uint64_t)r[3];

    random_scalars[0][1] = 0; random_scalars[0][0] = scalars[0];
    random_scalars[1][1] = 0; random_scalars[1][0] = scalars[1];
    random_scalars[2][1] = 0; random_scalars[2][0] = scalars[2];
    random_scalars[3][1] = 0; random_scalars[3][0] = scalars[3];

    MUL128(r0, p11, tt0);
    ADD128((digit_t*)&random_scalars[0], tt0, (digit_t*)&random_scalars[0]);
    MUL128(r1, p21, tt0);
    ADD128((digit_t*)&random_scalars[0], tt0, (digit_t*)&random_scalars[0]);
    MUL128(r2, p31, tt0);
    ADD128((digit_t*)&random_scalars[0], tt0, (digit_t*)&random_scalars[0]);
    MUL128(r3, p41, tt0);
    ADD128((digit_t*)&random_scalars[0], tt0, (digit_t*)&random_scalars[0]);

    MUL128(r0, p12, tt0);
    ADD128((digit_t*)&random_scalars[1], tt0, (digit_t*)&random_scalars[1]);
    MUL128(r1, p22, tt0);
    ADD128((digit_t*)&random_scalars[1], tt0, (digit_t*)&random_scalars[1]);
    MUL128(r2, p32, tt0);
    ADD128((digit_t*)&random_scalars[1], tt0, (digit_t*)&random_scalars[1]);
    MUL128(r3, p42, tt0);
    ADD128((digit_t*)&random_scalars[1], tt0, (digit_t*)&random_scalars[1]);

    MUL128(r0, p13, tt0);
    ADD128((digit_t*)&random_scalars[2], tt0, (digit_t*)&random_scalars[2]);
    MUL128(r1, p23, tt0);
    ADD128((digit_t*)&random_scalars[2], tt0, (digit_t*)&random_scalars[2]);
    MUL128(r2, p33, tt0);
    ADD128((digit_t*)&random_scalars[2], tt0, (digit_t*)&random_scalars[2]);
    MUL128(r3, p43, tt0);
    ADD128((digit_t*)&random_scalars[2], tt0, (digit_t*)&random_scalars[2]);

    MUL128(r0, p14, tt0);
    ADD128((digit_t*)&random_scalars[3], tt0, (digit_t*)&random_scalars[3]);
    MUL128(r1, p24, tt0);
    ADD128((digit_t*)&random_scalars[3], tt0, (digit_t*)&random_scalars[3]);
    MUL128(r2, p34, tt0);
    ADD128((digit_t*)&random_scalars[3], tt0, (digit_t*)&random_scalars[3]);
    MUL128(r3, p44, tt0);
    ADD128((digit_t*)&random_scalars[3], tt0, (digit_t*)&random_scalars[3]);
}


void recode(uint128_t* scalars, unsigned int* digits)
{ // Recoding sub-scalars for use in the variable-base scalar multiplication. 
  // Input: 4 80-bit sub-scalars passed through "scalars", which are obtained after calling randomize().
  // Outputs: "digits" array with 80 nonzero entries. Each entry is in the range [0, 15], corresponding to one entry in the precomputed table.
    unsigned int i, bit;

    for (i = 0; i < 64; i++)
    {
        bit = (unsigned int)(scalars[0][0]) & 1;
        scalars[0][0] >>= 1;
        digits[i] = bit;
        
        bit = (unsigned int)(scalars[1][0]) & 1;
        scalars[1][0] >>= 1;
        digits[i] += (bit << 1);
        
        bit = (unsigned int)(scalars[2][0]) & 1;
        scalars[2][0] >>= 1;
        digits[i] += (bit << 2);
        
        bit = (unsigned int)(scalars[3][0]) & 1;
        scalars[3][0] >>= 1;
        digits[i] += (bit << 3);
    }

    for (i = 64; i < 80; i++)
    {
        bit = (unsigned int)scalars[0][1] & 1;
        scalars[0][1] >>= 1;
        digits[i] = bit;
        
        bit = (unsigned int)scalars[1][1] & 1;
        scalars[1][1] >>= 1;
        digits[i] += (bit << 1);
        
        bit = (unsigned int)scalars[2][1] & 1;
        scalars[2][1] >>= 1;
        digits[i] += (bit << 2);
        
        bit = (unsigned int)scalars[3][1] & 1;
        scalars[3][1] >>= 1;
        digits[i] += (bit << 3);
    }
}


void cofactor_clearing(point_extedwards_t P)
{ // Co-factor clearing
  // Input: P = (X1,Y1,Z1,T1) in extended twisted Edwards coordinates
  // Output: P = 392*P = (Xfinal,Yfinal,Zfinal,Tfinal) in extended twisted Edwards coordinates
    point_extedwards_t Q;
     
    ecccopy(P, Q);                
    eccdouble(P);                        // P = 2*P using representations (X,Y,Z,T) <- 2*(X,Y,Z,T)
    eccadd(P, Q, P);                     // P = P+Q using representations (X,Y,Z,T) <- (X,Y,Z,T) + (X,Y,Z,T)
    eccdouble(P);
    eccdouble(P);
    eccdouble(P);
    eccdouble(P);
    eccadd(P, Q, P);
    eccdouble(P);
    eccdouble(P);
    eccdouble(P);
}


void select_f2elm(f2elm_t a, f2elm_t b, digit_t bit, f2elm_t c)      
{ // Select c <- a if bit = 0, c <- b if bit = 1
    unsigned int j;
    digit_t mask, temp, value = 0xAAAAAAAA;

    // If digit=0 mask = 0x55...5 else mask = 0xAA...A
    mask = bit - 1;
    mask = (mask & ~value) | (~mask & value);

    for (j = 0; j < NWORDS_FIELD; j++) {
        temp = a[0][j] ^ b[0][j];
        c[0][j] = ((mask & temp) ^ a[0][j]) ^ (value & temp); 
        temp = a[1][j] ^ b[1][j];
        c[1][j] = ((mask & temp) ^ a[1][j]) ^ (value & temp);
    }        
}


bool ecc_mul_SCA_secure(point_t P, point_t R, digit_t* k, point_t Q, bool clear_cofactor)
{ // Variable-base scalar multiplication Q = k*P using a 4-dimensional decomposition and protected against side-channel attacks
  // The computation is executed as Q = (k*P + R) - R with a random point R used for point blinding
  // Inputs: scalar "k" in [0, 2^256-1],
  //         point P = (x,y) in affine coordinates,
  //         random point R = (x,y) in affine coordinates,
  //         clear_cofactor = 1 (TRUE) or 0 (FALSE) whether cofactor clearing is required or not, respectively.
  // Output: Q = k*P in affine coordinates (x,y).
  //         R is updated to 2*R.
  // This function performs point validation and (if selected) cofactor clearing.
    point_extedwards_t PP, RR, S, Table[16];
    uint64_t scalars[NWORDS64_ORDER];
    uint128_t rand_scalars[NWORDS64_ORDER];
    unsigned int digits[80];
    unsigned char rand_bytes[8];
    digit_t bit;
    f2elm_t Ta, Tb;
    felm_t rand_felmt[82]; 
    int i;

    point_setup(P, PP);                                       // Convert to representation (X,Y,1,T)    
    if (ecc_point_validate(PP) == false) {                    // Check if point lies on the curve
        return false;
    }
    if (clear_cofactor == true) {
        cofactor_clearing(PP);
    }

    point_setup(R, RR);                                       // Convert to representation (X,Y,1,T)    
    if (ecc_point_validate(RR) == false) {                    // Check if blinding point lies on the curve
        return false;
    }

    RandomBytesFunction((unsigned char*)rand_felmt[0], 82*16);
    bit = rand_felmt[81][NWORDS_FIELD-1] >> (RADIX-1);
    for (i = 0; i < 82; i++) {    
        rand_felmt[i][NWORDS_FIELD-1] &= (digit_t)(-1) >> 1;
    }
    
    felmt_randomize_point(RR, rand_felmt[81]);                // Randomization of R's coordinates  
    ecccopy(RR, S);                                           // Update R = ((-1)^b*3)*R  
    eccdouble(RR);                                            
    eccadd(S, RR, S);    
    ecccopy(S, RR);                                           
    fp2neg1271(S->y);                                          
    fp2neg1271(S->t); 
    select_f2elm(RR->y, S->y, bit, RR->y); 
    select_f2elm(RR->t, S->t, bit, RR->t);                                     
    felmt_randomize_point(PP, rand_felmt[80]);                // Randomization of P's coordinates                
     
    decompose((uint64_t*) k, scalars);                        // Scalar decomposition
    RandomBytesFunction(&rand_bytes[0], 8);
    randomize(scalars, rand_bytes, rand_scalars);             // Scalar randomization
    recode(rand_scalars, digits);                             // Scalar recoding
    ecc_precomp(PP, RR, Table);                               // Precomputation
    ecccopy(RR, PP);                                          
    
    for (i = 79; i >= 0; i--)
    {
        eccdouble(PP);                                        // P = 2*P using representations (X,Y,Z,T) <- 2*(X,Y,Z)
#ifdef FULL_TABLE_RANDOMIZATION
        randomize_table(Table, rand_felmt[i]);                // Randomization of the full table
        table_lookup_1x16(Table, S, digits[i]);               // Extract point S in (X+Y,Y-X,2Z,2dT) representation
#else   
        table_lookup_1x16(Table, S, digits[i]);               // Extract point S in (X+Y,Y-X,2Z,2dT) representation
        felmt_randomize_point(S, rand_felmt[i]);              // Randomization of the extracted point
#endif
        eccadd_core(PP, S, PP, Ta, Tb);                       // P = P+S using representations (X,Y,Z,Ta,Tb) <- (X,Y,Z,T) + (X+Y,Y-X,2Z,2dT)
    }
                                                  
    fp2mul1271(Ta, Tb, PP->t);
    eccadd_core(PP, Table[0], PP, Ta, Tb);                    // Final correction: (k*P + R) - R
    eccnorm2(PP, Q, RR, R);                                   // Converting output and blinding point to affine coordinates (x,y). 

    return true;
}


void eccset(point_t P)
{ // Set generator  
  // Output: P = (x,y)

	fp2copy1271((felm_t*)&GENERATOR_x, P->x);    // X1
	fp2copy1271((felm_t*)&GENERATOR_y, P->y);    // Y1
}


void eccnorm(point_extedwards_t P, point_t Q)
{ // Normalize a projective point (X1:Y1:Z1), including full reduction
  // Input: P = (X1:Y1:Z1) in twisted Edwards coordinates    
  // Output: Q = (X1/Z1,Y1/Z1), corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
    
    fp2inv1271(P->z);                      // Z1 = Z1^-1
    fp2mul1271(P->x, P->z, Q->x);          // X1 = X1/Z1
    fp2mul1271(P->y, P->z, Q->y);          // Y1 = Y1/Z1
    mod1271(Q->x[0]); mod1271(Q->x[1]); 
    mod1271(Q->y[0]); mod1271(Q->y[1]); 
}


void eccnorm2(point_extedwards_t P, point_t Q, point_extedwards_t R, point_t S)
{ // Normalize two projective points, including full reduction
  // Inputs:  P = (X1:Y1:Z1) and R = (X2:Y2:Z2) in twisted Edwards coordinates    
  // Outputs: Q = (X1/Z1,Y1/Z1) and S = (X2/Z2,Y2/Z2) in affine coordinates
    f2elm_t t1;
    
    fp2mul1271(P->z, R->z, t1);            // t1 = Z1*Z2
    fp2inv1271(t1);                        // t1 = (Z1*Z2)^-1
    fp2mul1271(R->z, t1, Q->y);            // yQ = Z1^-1
    fp2mul1271(Q->y, P->x, Q->x);          // xQ = X1/Z1
    fp2mul1271(Q->y, P->y, Q->y);          // yQ = Y1/Z1
    mod1271(Q->x[0]); mod1271(Q->x[1]); 
    mod1271(Q->y[0]); mod1271(Q->y[1]); 
    fp2mul1271(t1, P->z, t1);              // t1 = Z2^-1
    fp2mul1271(R->x, t1, S->x);            // xS = X2/Z2
    fp2mul1271(R->y, t1, S->y);            // yS = Y2/Z2
    mod1271(S->x[0]); mod1271(S->x[1]); 
    mod1271(S->y[0]); mod1271(S->y[1]); 
}


void eccdouble(point_extedwards_t P)
{ // Point doubling 2P
  // Input: P = (X1:Y1:Z1:T1) in twisted Edwards coordinates
  // Output: 2P = (Xfinal,Yfinal,Zfinal,Tfinal) in extended twisted Edwards coordinates
    f2elm_t t1, t2, t3;  

    fp2sqr1271(P->x, t1);                  // t1 = X1^2
    fp2sqr1271(P->y, P->t);                // T = Y1^2
    fp2add1271(P->x, P->y, P->x);          // X1 = X1+Y1
    fp2add1271(t1, P->t, t2);              // t2 = X1^2+Y1^2      
    fp2sub1271(P->t, t1, t1);              // t1 = Y1^2-X1^2  
    fp2sqr1271(P->z, P->t);                // T = Z1^2       
    fp2sqr1271(P->x, t3);                  // t3 = (X1+Y1)^2  
    fp2add1271(P->t, P->t, P->t);          // T = 2Z1^2   
    fp2sub1271(t3, t2, t3);                // t3 = 2X1*Y1 = (X1+Y1)^2-(X1^2+Y1^2)
    fp2sub1271(P->t, t1, P->t);            // T = 2Z1^2-(Y1^2-X1^2) 
    fp2mul1271(t1, t2, P->y);              // Yfinal = (X1^2+Y1^2)(Y1^2-X1^2)  
    fp2mul1271(P->t, t3, P->x);            // Xfinal = 2X1Y1*[2Z1^2-(Y1^2-X1^2)]
    fp2mul1271(t1, P->t, P->z);            // Zfinal = (Y1^2-X1^2)[2Z1^2-(Y1^2-X1^2)]
    fp2mul1271(t2, t3, P->t);              // Tfinal = 2X1Y1*(X1^2+Y1^2)
}


void eccadd_core(point_extedwards_t P, point_extedwards_t Q, point_extedwards_t R, f2elm_t Ta, f2elm_t Tb)      
{ // Basic point addition R = P+Q or R = P+P
  // Inputs: P = (X1,Y1,Z1,T1) in extended twisted Edwards coordinates
  //         Q = (X2+Y2,Y2-X2,2Z2,2dT2) in extended twisted Edwards coordinates    
  // Output: R = (Xfinal,Yfinal,Zfinal,Tfinal) in extended twisted Edwards coordinates. Tfinal is only used for temporary storage
  //         Tafinal and Tbfinal such that Tfinal = Tafinal*Tbfinal  
    f2elm_t t1; 
          
    fp2mul1271(P->z, Q->z, R->z);                   // Z = 2Z1*Z2  
    fp2mul1271(P->t, Q->t, R->t);                   // T = 2dT1*T2 
    fp2add1271(P->x, P->y, Ta);                     // Ta = X1+Y1
    fp2mul1271(Ta, Q->x, t1);                       // t1 = (X1+Y1)(X2+Y2) 
    fp2sub1271(P->x, P->y, Ta);                     // Ta = X1-Y1
    fp2add1271(R->z, R->t, R->y);                   // Y = alpha
    fp2mul1271(Ta, Q->y, Ta);                       // Ta = (X1-Y1)(X2-Y2) 
    fp2sub1271(R->z, R->t, R->z);                   // Z = theta
    fp2sub1271(t1, Ta, Tb);                         // Tbfinal = beta = (X1+Y1)(X2+Y2)-(X1-Y1)(X2-Y2)
    fp2add1271(t1, Ta, Ta);                         // Tafinal = omega = (X1+Y1)(X2+Y2)+(X1-Y1)(X2-Y2)
    fp2mul1271(Tb, R->z, R->x);                     // Xfinal = beta*theta
    fp2mul1271(R->z, R->y, R->z);                   // Zfinal = theta*alpha
    fp2mul1271(R->y, Ta, R->y);                     // Yfinal = alpha*omega
}


void eccadd(point_extedwards_t P, point_extedwards_t Q, point_extedwards_t R)      
{ // Complete point addition P = P+Q or P = P+P
  // Inputs: P = (X1,Y1,Z1,T1) in extended twisted Edwards coordinates
  //         Q = (X2,Y2,Z2,T2) in extended twisted Edwards coordinates   
  // Output: P = (Xfinal,Yfinal,Zfinal,Tfinal) in extended twisted Edwards coordinates
    f2elm_t Ta, Tb, t1; 
          
    fp2mul1271(P->z, Q->z, R->z);                   // Z = Z1*Z2  
    fp2mul1271(P->t, Q->t, R->t);                   // T = T1*T2 
    fp2add1271(R->z, R->z, R->z);                   // Z = 2Z1*Z2
    fp2add1271(R->t, R->t, R->t);                   // T = 2T1*T2
    fp2add1271(P->x, P->y, Ta);                     // Ta = X1+Y1
    fp2add1271(Q->x, Q->y, Tb);                     // Tb = X2+Y2 
    fp2mul1271(R->t, (felm_t*)&PARAMETER_d, R->t);  // T = 2d*T1*T2 
    fp2mul1271(Ta, Tb, t1);                         // t1 = (X1+Y1)(X2+Y2) 
    fp2sub1271(P->x, P->y, Ta);                     // Ta = X1-Y1
    fp2sub1271(Q->x, Q->y, Tb);                     // Tb = X2-Y2 
    fp2add1271(R->z, R->t, R->y);                   // Y = alpha
    fp2mul1271(Ta, Tb, Ta);                         // Ta = (X1-Y1)(X2-Y2) 
    fp2sub1271(R->z, R->t, R->z);                   // Z = theta
    fp2sub1271(t1, Ta, Tb);                         // Tb = beta = (X1+Y1)(X2+Y2)-(X1-Y1)(X2-Y2)
    fp2add1271(t1, Ta, Ta);                         // Ta = omega = (X1+Y1)(X2+Y2)+(X1-Y1)(X2-Y2)
    fp2mul1271(Tb, R->z, R->x);                     // Xfinal = beta*theta
    fp2mul1271(R->z, R->y, R->z);                   // Zfinal = theta*alpha
    fp2mul1271(R->y, Ta, R->y);                     // Yfinal = alpha*omega
    fp2mul1271(Ta, Tb, R->t);                       // Tfinal = Ta*Tb
}


void point_setup(point_t P, point_extedwards_t Q)
{ // Point conversion to representation (X,Y,Z,T) 
  // Input: P = (x,y) in affine coordinates
  // Output: P = (X,Y,1,T) corresponding to (X:Y:Z:T) in extended twisted Edwards coordinates

    fp2copy1271(P->x, Q->x);
    fp2copy1271(P->y, Q->y);
    fp2mul1271(P->x, P->y, Q->t);          // T = X*Y
    fp2zero1271(Q->z); Q->z[0][0]=1;       // Z = 1
}


bool ecc_point_validate(point_extedwards_t P)
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


static __inline void eccmadd(point_precomp_t Q, point_extedwards_t P)
{ // Mixed point addition P = P+Q or P = P+P
  // Inputs: P = (X1,Y1,Z1,T1) in extended twisted Edwards coordinates
  //         Q = (x2+y2,y2-x2,2dt2) corresponding to (X2:Y2:Z2:T2) in extended twisted Edwards coordinates, where Z2=1  
  // Output: P = (Xfinal,Yfinal,Zfinal,Tfinal) in extended twisted Edwards coordinates
    f2elm_t t1, t2, t3;
    
    fp2add1271(P->z, P->z, t1);            // t1 = 2Z1        
    fp2mul1271(P->t, Q->t2, P->t);         // T = 2dT1*t2 
    fp2add1271(P->x, P->y, P->z);          // Z = (X1+Y1) 
    fp2sub1271(P->y, P->x, t3);            // t3 = (Y1-X1)
    fp2sub1271(t1, P->t, t2);              // t2 = theta
    fp2add1271(t1, P->t, t1);              // t1 = alpha
    fp2mul1271(Q->xy, P->z, P->t);         // T = (X1+Y1)(x2+y2)
    fp2mul1271(Q->yx, t3, P->x);           // X = (Y1-X1)(y2-x2)
    fp2mul1271(t1, t2, P->z);              // Zfinal = theta*alpha
    fp2sub1271(P->t, P->x, t3);            // t3 = beta
    fp2add1271(P->t, P->x, P->t);          // T = omega
    fp2mul1271(t3, t2, P->x);              // Xfinal = beta*theta
    fp2mul1271(P->t, t1, P->y);            // Yfinal = alpha*omega
    fp2mul1271(P->t, t3, P->t);            // Tfinal = beta*omega
}


static __inline void eccneg_extproj_precomp(point_extedwards_t P, point_extedwards_t Q)
{ // Point negation
  // Input : point P in coordinates (X+Y,X-Y,2Z,2dT)
  // Output: point Q = -P = (Y-X,X+Y,2Z,-2dT)
    fp2copy1271(P->t, Q->t);
    fp2copy1271(P->x, Q->y);
    fp2copy1271(P->y, Q->x);
    fp2copy1271(P->z, Q->z);
    fp2neg1271(Q->x);
    fp2neg1271(Q->y);
    fp2neg1271(Q->t);
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

	unsigned int position;
	int i, digits_k1[65] = {0}, digits_k2[65] = {0}, digits_k3[65] = {0}, digits_k4[65] = {0};
	int digits_l1[65] = {0}, digits_l2[65] = {0}, digits_l3[65] = {0}, digits_l4[65] = {0};
	point_precomp_t V;
	point_extedwards_t Q1, Q2, Q3, Q4, T;
	point_extedwards_t U, Q_table1[NPOINTS_DOUBLEMUL_WQ], Q_table2[NPOINTS_DOUBLEMUL_WQ], Q_table3[NPOINTS_DOUBLEMUL_WQ], Q_table4[NPOINTS_DOUBLEMUL_WQ];
    f2elm_t t0, t1;
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

	fp2zero1271(T->x);                                         // Initialize T as the neutral point (0:1:1)
	fp2zero1271(T->y); T->y[0][0] = 1;
	fp2zero1271(T->z); T->z[0][0] = 1;

	for (i = 64; i >= 0; i--)
	{
		eccdouble(T);                                          // Double (X_T,Y_T,Z_T,T_T) = 2(X_T,Y_T,Z_T,T_T)
		if (digits_l1[i] < 0) {
			position = (-digits_l1[i])/2;
			eccneg_extproj_precomp(Q_table1[position], U);     // Load and negate U = (X_U,Y_U,Z_U,Td_U) <- -(X+Y,X-Y,2Z,2dT) from a point in the precomputed table 
			eccadd_core(T, U, T, t0, t1);                      // T = T+U = (X_T,Y_T,Z_T,T_T) = (X_T,Y_T,Z_T,T_T) + (X_U,Y_U,Z_U,Td_U) 
            fp2mul1271(t0, t1, T->t);
		} else if (digits_l1[i] > 0) {
			position = (digits_l1[i])/2;                       // Take U = (X_U,Y_U,Z_U,Td_U) <- (X+Y,X-Y,2Z,2dT) from a point in the precomputed table
			eccadd_core(T, Q_table1[position], T, t0, t1);     // T = T+U = (X_T,Y_T,Z_T,T_T) = (X_T,Y_T,Z_T,T_T) + (X_U,Y_U,Z_U,Td_U) 
            fp2mul1271(t0, t1, T->t);
		}
		if (digits_l2[i] < 0) {
			position = (-digits_l2[i])/2;
			eccneg_extproj_precomp(Q_table2[position], U);
			eccadd_core(T, U, T, t0, t1); 
            fp2mul1271(t0, t1, T->t);
		} else if (digits_l2[i] > 0) {
			position = (digits_l2[i])/2;
			eccadd_core(T, Q_table2[position], T, t0, t1);
            fp2mul1271(t0, t1, T->t);
		}
		if (digits_l3[i] < 0) {
			position = (-digits_l3[i])/2;
			eccneg_extproj_precomp(Q_table3[position], U);
			eccadd_core(T, U, T, t0, t1); 
            fp2mul1271(t0, t1, T->t);
		} else if (digits_l3[i] > 0) {
			position = (digits_l3[i])/2;
			eccadd_core(T, Q_table3[position], T, t0, t1);
            fp2mul1271(t0, t1, T->t);
		}
		if (digits_l4[i] < 0) {
			position = (-digits_l4[i])/2;
			eccneg_extproj_precomp(Q_table4[position], U);
			eccadd_core(T, U, T, t0, t1); 
            fp2mul1271(t0, t1, T->t);
		} else if (digits_l4[i] > 0) {
			position = (digits_l4[i])/2;
			eccadd_core(T, Q_table4[position], T, t0, t1);
            fp2mul1271(t0, t1, T->t);
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
	eccnorm(T, R);                                             // Output R = (x,y)

	return true;
}


static __inline void R_to_R2(point_extedwards_t P, point_extedwards_t Q) 
{ // Conversion from representation (X,Y,Z,T) to (X+Y,X-Y,2Z,2dT), where T = Ta*Tb
  // Input:  P = (X1,Y1,Z1,T1) in extended twisted Edwards coordinates
  // Output: Q = (X1+Y1,Y1-X1,2Z1,2dT1) corresponding to (X1:Y1:Z1:T1) in extended twisted Edwards coordinates
    
    fp2add1271(P->t, P->t, Q->t);                    // T = 2*T
    fp2add1271(P->x, P->y, Q->x);                    // QX = X+Y
    fp2sub1271(P->x, P->y, Q->y);                    // QY = X-Y
    fp2add1271(P->z, P->z, Q->z);                    // QZ = 2*Z
	fp2mul1271(Q->t, (felm_t*)&PARAMETER_d, Q->t);   // QT = 2d*T
}


void ecc_precomp_double(point_extedwards_t P, point_extedwards_t* Table, unsigned int npoints)
{ // Generation of the precomputation table used internally by the double scalar multiplication function ecc_mul_double().  
  // Inputs: point P in representation (X,Y,Z,T),
  //         Table with storage for npoints, 
  //         number of points "npoints".
  // Output: Table containing multiples of the base point P using representation (X+Y,X-Y,2Z,2dT).
	point_extedwards_t Q;
	unsigned int i;

	R_to_R2(P, Table[0]);                     // Precomputed point Table[0] = P in coordinates (X+Y,X-Y,2Z,2dT)
    ecccopy(P, Q);
	eccdouble(P);                             // 2*P in (X,Y,Z,T)

	for (i = 1; i < npoints; i++) {
		eccadd(Q, P, Q);                      // Table[i] = Table[i-1]+2P
		R_to_R2(Q, Table[i]);                 // Converting from (X,Y,Z,T) to (X+Y,X-Y,2Z,2dT)
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