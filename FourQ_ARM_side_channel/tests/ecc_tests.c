/***********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: testing code for FourQ's curve arithmetic 
************************************************************************************/   

#include "../FourQ_api.h"
#include "../FourQ_params.h"
#include "../FourQ_tables.h"
#include "test_extras.h"
#include <stdio.h>


// Benchmark and test parameters 
#define BENCH_LOOPS       100       // Number of iterations per bench
#define SHORT_BENCH_LOOPS 10        // Number of iterations per bench (for expensive operations)
#define TEST_LOOPS        1000      // Number of iterations per test


bool ecc_test()
{
    bool clear_cofactor, OK = true;
    unsigned int n;
    int passed;
    point_t A, R;
    point_extedwards_t P, Q;
    f2elm_t t1;
    uint64_t scalar[4], res_x[4], res_y[4];

    
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Testing FourQ's curve arithmetic: \n\n"); 

    // Point doubling
    passed = 1;
    eccset(A); 
    point_setup(A, P);

    for (n=0; n<TEST_LOOPS; n++)
    {
        eccdouble(P);                      // 2*P
    }
    eccnorm(P, A);
    mod1271(A->x[0]); mod1271(A->x[1]);    // Fully reduced P
    mod1271(A->y[0]); mod1271(A->y[1]);  

    // Result
    res_x[0] = 0xC9099C54855859D6; res_x[1] = 0x2C3FD8822C82270F; res_x[2] = 0xA7B3F6E2043E8E68; res_x[3] = 0x4DA5B9E83AA7A1B2;
    res_y[0] = 0x3EE089F0EB49AA14; res_y[1] = 0x2001EB3A57688396; res_y[2] = 0x1FEE5617A7E954CD; res_y[3] = 0x0FFDB0D761421F50;

    if (fp2compare64((uint64_t*)A->x, res_x)!=0 || fp2compare64((uint64_t*)A->y, res_y)!=0) passed=0;
    if (passed==1) printf("  Point doubling tests .................................................................... PASSED");
    else { printf("  Point doubling tests ... FAILED"); printf("\n"); return false; }
    printf("\n");
   
    // Point addition
    eccset(A); 
    point_setup(A, P);
    
    for (n=0; n<TEST_LOOPS; n++)
    {
        ecccopy(P, Q);
        eccadd(P, Q, P);         // 2*P
    }
    eccnorm(P, A);
    mod1271(A->x[0]); mod1271(A->x[1]);    // Fully reduced P
    mod1271(A->y[0]); mod1271(A->y[1]);  

    // Result
    res_x[0] = 0xC9099C54855859D6; res_x[1] = 0x2C3FD8822C82270F; res_x[2] = 0xA7B3F6E2043E8E68; res_x[3] = 0x4DA5B9E83AA7A1B2;
    res_y[0] = 0x3EE089F0EB49AA14; res_y[1] = 0x2001EB3A57688396; res_y[2] = 0x1FEE5617A7E954CD; res_y[3] = 0x0FFDB0D761421F50;

    if (fp2compare64((uint64_t*)A->x, res_x)!=0 || fp2compare64((uint64_t*)A->y, res_y)!=0) passed=0;

    eccset(A); 
    point_setup(A, P);
    ecccopy(P, Q);
    eccdouble(P);                          // P = 2P 

    for (n=0; n<TEST_LOOPS; n++)
    {
        eccadd(P, Q, P);         // P = P+Q
    }    
    eccnorm(P, A);
    mod1271(A->x[0]); mod1271(A->x[1]);    // Fully reduced P
    mod1271(A->y[0]); mod1271(A->y[1]);    

    // Result
    res_x[0] = 0x6480B1EF0A151DB0; res_x[1] = 0x3E243958590C4D90; res_x[2] = 0xAA270F644A65D473; res_x[3] = 0x5327AF7D84238CD0;
    res_y[0] = 0x5E06003D73C43EB1; res_y[1] = 0x3EF69A49CB7E0237; res_y[2] = 0x4E752648AC2EF0AB; res_y[3] = 0x293EB1E26DD23B4E;

    if (fp2compare64((uint64_t*)A->x, res_x)!=0 || fp2compare64((uint64_t*)A->y, res_y)!=0) passed=0;

    if (passed==1) printf("  Point addition tests .................................................................... PASSED");
    else { printf("  Point addition tests ... FAILED"); printf("\n"); return false; }
    printf("\n");
   
    // Psi endomorphism
    eccset(A); 
    point_setup(A, P);

    for (n=0; n<TEST_LOOPS; n++)
    {
        ecc_psi(P);                        // P = Psi(P)
    }    
    eccnorm(P, A);
    mod1271(A->x[0]); mod1271(A->x[1]);    // Fully reduced P
    mod1271(A->y[0]); mod1271(A->y[1]);   

    // Result
    res_x[0] = 0xD8F3C8C24A2BC7E2; res_x[1] = 0x75AF54EDB41A2B93; res_x[2] = 0x4DE2466701F009A9; res_x[3] = 0x065249F9EDE0C798;
    res_y[0] = 0x1C6E119ADD608104; res_y[1] = 0x06DBB85BFFB7C21E; res_y[2] = 0xFD234D6C4CFA3EC1; res_y[3] = 0x060A30903424BF13;

    if (fp2compare64((uint64_t*)A->x, res_x)!=0 || fp2compare64((uint64_t*)A->y, res_y)!=0) passed=0;

    if (passed==1) printf("  Psi endomorphism tests .................................................................. PASSED");
    else { printf("  Psi endomorphism tests ... FAILED"); printf("\n"); return false; }
    printf("\n");
   
    // Phi endomorphism
    eccset(A); 
    point_setup(A, P);

    for (n=0; n<TEST_LOOPS; n++)
    {
        ecc_phi(P);                        // P = Phi(P)
        eccnorm(P, A);
        point_setup(A, P);
    }    
    mod1271(A->x[0]); mod1271(A->x[1]);    // Fully reduced P
    mod1271(A->y[0]); mod1271(A->y[1]); 

    // Result
    res_x[0] = 0xD5B5A3061287DB16; res_x[1] = 0x5550AAB9E7A620EE; res_x[2] = 0xEC321E6CF33610FC; res_x[3] = 0x3E61EBB9A1CB0210;
    res_y[0] = 0x7E2851D5A8E83FB9; res_y[1] = 0x5474BF8EC55603AE; res_y[2] = 0xA5077613491788D5; res_y[3] = 0x5476093DBF8BF6BF;

    if (fp2compare64((uint64_t*)A->x, res_x)!=0 || fp2compare64((uint64_t*)A->y, res_y)!=0) passed=0;
    if (passed==1) printf("  Phi endomorphism tests .................................................................. PASSED");
    else { printf("  Phi endomorphism tests ... FAILED"); printf("\n"); return false; }
    printf("\n");
    
    // Scalar multiplication
    eccset(A); 
    eccset(R); 
    clear_cofactor = false;
    scalar[0] = 0x3AD457AB55456230; scalar[1] = 0x3A8B3C2C6FD86E0C; scalar[2] = 0x7E38F7C9CFBB9166; scalar[3] = 0x0028FD6CBDA458F0;
    
    for (n=0; n<TEST_LOOPS; n++)
    {
        scalar[1] = scalar[2];
        scalar[2] += scalar[0];
        
        ecc_mul_SCA_secure(A, R, (digit_t*)scalar, A, clear_cofactor);
    }

    res_x[0] = 0xDFD2B477BD494BEF; res_x[1] = 0x257C122BBFC94A1B; res_x[2] = 0x769593547237C459; res_x[3] = 0x469BF80CB5B11F01;
    res_y[0] = 0x281C5067996F3344; res_y[1] = 0x0901B3817C0E936C; res_y[2] = 0x4FE8C429915F1245; res_y[3] = 0x570B948EACACE210;
    if (fp2compare64((uint64_t*)A->x, res_x)!=0 || fp2compare64((uint64_t*)A->y, res_y)!=0) passed=0;

        
    eccset(A); 
    eccset(R); 
    clear_cofactor = true;
    scalar[0] = 0x3AD457AB55456230; scalar[1] = 0x3A8B3C2C6FD86E0C; scalar[2] = 0x7E38F7C9CFBB9166; scalar[3] = 0x0028FD6CBDA458F0;
    
    for (n=0; n<TEST_LOOPS; n++)
    {
        scalar[1] = scalar[2];
        scalar[2] += scalar[0];
        
        ecc_mul_SCA_secure(A, R, (digit_t*)scalar, A, clear_cofactor);
    }

    res_x[0] = 0x85CF54A3BEE3FD23; res_x[1] = 0x7A7EC43976FAAD92; res_x[2] = 0x7697567B785E2327; res_x[3] = 0x4CBDAB448B1539F2;
    res_y[0] = 0xE9193B41CDDF94D0; res_y[1] = 0x5AA6C859ECC810D5; res_y[2] = 0xAA876E760AA8B331; res_y[3] = 0x320C53F02230094A;
    if (fp2compare64((uint64_t*)A->x, res_x)!=0 || fp2compare64((uint64_t*)A->y, res_y)!=0) passed=0;
        
    if (passed==1) printf("  Scalar multiplication tests ............................................................. PASSED");
    else { printf("  Scalar multiplication tests ... FAILED"); printf("\n"); return false; }
    printf("\n");

    {    
    point_t PP, QQ, RR, UU, TT; 
    point_extedwards_t AA, BB;
    uint64_t k[4], l[4], kk[4];

    // Double scalar multiplication
    eccset(QQ); 
    eccset(PP);
    eccset(R);
    
    for (n=0; n<TEST_LOOPS; n++)
    {
        random_scalar_test(kk); 
        ecc_mul_SCA_secure(QQ, R, (digit_t*)kk, QQ, false);
        random_scalar_test(k); 
        random_scalar_test(l);
        ecc_mul_double((digit_t*)k, QQ, (digit_t*)l, RR);
        ecc_mul_SCA_secure(PP, R, (digit_t*)k, UU, false);
        ecc_mul_SCA_secure(QQ, R, (digit_t*)l, TT, false);    
        point_setup(UU, AA);
        point_setup(TT, BB);
        eccadd(BB, AA, BB);
        eccnorm(BB, UU);
        
        if (fp2compare64((uint64_t*)UU->x,(uint64_t*)RR->x)!=0 || fp2compare64((uint64_t*)UU->y,(uint64_t*)RR->y)!=0) { passed=0; break; }
    }
    
    if (passed==1) printf("  Double scalar multiplication tests ...................................................... PASSED");
    else { printf("  Double scalar multiplication tests ... FAILED"); printf("\n"); return false; }
    printf("\n");
    }
    
    return OK;
}


bool ecc_run()
{
    bool OK = true;
    unsigned int n, i, digit=1;
    unsigned long long nsec, nsec1, nsec2;
    point_t A, B, R;
    point_extedwards_t P, Q, RR, Table[16];
    f2elm_t t1;    
    uint64_t scalar[4];
        
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Benchmarking FourQ's curve arithmetic \n\n"); 

    // Point doubling (twisted Edwards a=-1)
    eccset(A);
    point_setup(A, P);

    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        nsec1 = cpu_nseconds();
        for (i = 0; i < 1000; i++) {
            eccdouble(P);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    printf("  Point doubling runs in ...                                       %8lld nsec", nsec/(BENCH_LOOPS*1000));
    printf("\n");

    // Point addition (twisted Edwards a=-1)
    eccset(A);
    point_setup(A, P);
    ecccopy(P, Q);
    eccdouble(P);                       // P = 2P 

    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        nsec1 = cpu_nseconds();
        for (i = 0; i < 1000; i++) {
            eccadd(P, Q, P);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    printf("  Point addition runs in ...                                       %8lld nsec", nsec/(BENCH_LOOPS*1000));
    printf("\n");

    // Psi endomorphism
    eccset(A);
    point_setup(A, P);

    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        nsec1 = cpu_nseconds();
        for (i = 0; i < 1000; i++) {
            ecc_psi(P);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    printf("  Psi mapping runs in ...                                          %8lld nsec", nsec/(BENCH_LOOPS*1000));
    printf("\n");
   
    // Phi endomorphism
    eccset(A);
    point_setup(A, P);

    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        nsec1 = cpu_nseconds();
        for (i = 0; i < 1000; i++) {
            ecc_phi(P);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    printf("  Phi mapping runs in ...                                          %8lld nsec", nsec/(BENCH_LOOPS*1000));
    printf("\n");
   
    // Scalar decomposition
    {
    uint64_t scalars[4];
    random_scalar_test(scalar); 
    
    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        nsec1 = cpu_nseconds();
        for (i = 0; i < 1000; i++) {
            decompose(scalar, scalars);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    printf("  Scalar decomposition runs in ...                                 %8lld nsec", nsec/(BENCH_LOOPS*1000));
    printf("\n");
    }

    // Scalar recoding
    {
    unsigned int digits[72];
    uint128_t rand_scalars[4];
    unsigned char rand_bytes[4] = {1};

    random_scalar_test(scalar);
    randomize(scalar, rand_bytes, rand_scalars);

    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        nsec1 = cpu_nseconds();
        nsec1 = cpu_nseconds();
        for (i = 0; i < 1000; i++) {
            recode(rand_scalars, digits);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    printf("  Scalar recoding runs in ...                                      %8lld nsec", nsec/(BENCH_LOOPS*1000));
    printf("\n");  
    }
    
    // Precomputation
    eccset(A);
    point_setup(A, P);
    ecccopy(P, RR);

    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        nsec1 = cpu_nseconds();
        for (i = 0; i < 100; i++) {
            ecc_precomp(P, RR, Table);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }    
    printf("  Precomputation runs in ...                                       %8lld nsec", nsec/(BENCH_LOOPS*100));
    printf("\n");  

    // Table lookup
    eccset(A);
    point_setup(A, P);
    ecccopy(P, RR);
    ecc_precomp(P, RR, Table);

    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        nsec1 = cpu_nseconds();
        for (i = 0; i < 1000; i++) {
            table_lookup_1x16(Table, Q, digit);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    } 
    
    printf("  Table lookup runs in ...                                         %8lld nsec", nsec/(BENCH_LOOPS*1000));
    printf("\n");  

    // Scalar multiplication
    random_scalar_test(scalar); 

    nsec = 0;
    for (n=0; n<SHORT_BENCH_LOOPS; n++)
    {
        eccset(A);
        eccset(R);
        nsec1 = cpu_nseconds();
        for (i = 0; i < 100; i++) {
            ecc_mul_SCA_secure(A, R, (digit_t*)scalar, B, false);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }    
    printf("  Scalar multiplication (without clearing cofactor) runs in ...... %8lld nsec", nsec/(SHORT_BENCH_LOOPS*100));
    printf("\n"); 
    
    random_scalar_test(scalar); 

    nsec = 0;
    for (n=0; n<SHORT_BENCH_LOOPS; n++)
    {
        eccset(A);
        eccset(R);
        nsec1 = cpu_nseconds();
        for (i = 0; i < 100; i++) {
            ecc_mul_SCA_secure(A, R, (digit_t*)scalar, B, true);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }    
    printf("  Scalar multiplication (including clearing cofactor) runs in .... %8lld nsec", nsec/(SHORT_BENCH_LOOPS*100));
    printf("\n"); 

    {    
    point_t PP, QQ, RR;
    uint64_t k[4], l[4], kk[4];

    // Double scalar multiplication
    eccset(QQ); 
    eccset(PP);
    random_scalar_test(kk); 
    ecc_mul_SCA_secure(QQ, PP, (digit_t*)kk, QQ, false);
    
    nsec = 0;
    for (n=0; n<SHORT_BENCH_LOOPS; n++)
    {        
        random_scalar_test(k); 
        random_scalar_test(l);  
        nsec1 = cpu_nseconds();
        ecc_mul_double((digit_t*)k, QQ, (digit_t*)l, RR);
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    
    printf("  Double scalar mul runs in ...                                    %8lld nsec with wP=%d and wQ=%d", nsec/SHORT_BENCH_LOOPS, WP_DOUBLEBASE, WQ_DOUBLEBASE);
    printf("\n"); 
    }  
    
    return OK;
} 


int main()
{
    bool OK = true;

    OK = OK && ecc_test();         // Test FourQ's curve functions
    OK = OK && ecc_run();          // Benchmark FourQ's curve functions
    
    return OK;
}
