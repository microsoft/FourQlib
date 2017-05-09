/***********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: testing code for FourQ's curve arithmetic 
************************************************************************************/   

#include "../FourQ_internal.h"
#include "../FourQ_params.h"
#include "../FourQ_tables.h"
#include "test_extras.h"
#include "../stm32f4_wrapper.h" 
#include <stdio.h>

static unsigned int *DWT_CYCCNT = (unsigned int*)0xE0001004;
static unsigned int *DWT_CTRL = (unsigned int*)0xE0001000;
static unsigned int *SCB_DEMCR = (unsigned int*)0xE000EDFC;

// Benchmark and test parameters 
#define BENCH_LOOPS       10       // Number of iterations per bench
#define TEST_LOOPS        10       // Number of iterations per test

#define cpucycles() (*DWT_CYCCNT);


static void print_test(const char *text)
{
    unsigned char output[100];

    sprintf((char*)output, "%s", text);
    send_USART_str(output);
}


static void print_bench(const char *s, unsigned int cycles)
{
  unsigned char output[100];

  sprintf((char*)output, "%s %8u cycles", s, cycles);
  send_USART_str(output);
}


bool ecc_test()
{
    bool clear_cofactor, OK = true;
    unsigned int n;
    int passed;
    point_t A, R;
    point_extedwards_t P, Q;
    f2elm_t t1;
    uint64_t scalar[4], res_x[4], res_y[4];
        
    print_test("\n--------------------------------------------------------------------------------------------------------\n"); 
    print_test("Testing FourQ's curve arithmetic: \n"); 

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
    res_x[0] = 0xFA4FAD9EC7732700; res_x[1] = 0x619F5D1FD93BC4F5; res_x[2] = 0x814B78DADF6A9024; res_x[3] = 0x72EC1D429F026578;
    res_y[0] = 0x7FF28C92C8CEF9DE; res_y[1] = 0x799208A76EAD2BA3; res_y[2] = 0x9B1AE60FFFCB520A; res_y[3] = 0x051698145D42F3E2;

    if (fp2compare64((uint64_t*)A->x, res_x)!=0 || fp2compare64((uint64_t*)A->y, res_y)!=0) passed=0;
    if (passed==1) print_test("  Point doubling tests .................................................................... PASSED");
    else { print_test("  Point doubling tests ... FAILED"); print_test("\n"); return false; }
   
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
    res_x[0] = 0xFA4FAD9EC7732700; res_x[1] = 0x619F5D1FD93BC4F5; res_x[2] = 0x814B78DADF6A9024; res_x[3] = 0x72EC1D429F026578;
    res_y[0] = 0x7FF28C92C8CEF9DE; res_y[1] = 0x799208A76EAD2BA3; res_y[2] = 0x9B1AE60FFFCB520A; res_y[3] = 0x051698145D42F3E2;

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
    res_x[0] = 0xB92B573D2C4B06FF; res_x[1] = 0x6B62D585800A9F6A; res_x[2] = 0xECB6DFB3FA1ACB7C; res_x[3] = 0xD9D9F54A8335E2B;
    res_y[0] = 0xDF3BD744D9BB783D; res_y[1] = 0x2B827EEDA23988A6; res_y[2] = 0x947C187247366CDD; res_y[3] = 0x3B7E00BA2F9525B3;

    if (fp2compare64((uint64_t*)A->x, res_x)!=0 || fp2compare64((uint64_t*)A->y, res_y)!=0) passed=0;

    if (passed==1) print_test("  Point addition tests .................................................................... PASSED");
    else { print_test("  Point addition tests ... FAILED"); print_test("\n"); return false; }
   
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
    res_x[0] = 0xABC340A7DDC08580; res_x[1] = 0x6B74D34E155D2119; res_x[2] = 0x1B6E0A6DC6A5BC70; res_x[3] = 0x5CAE354597C9106A;
    res_y[0] = 0xE276B58944E2D60B; res_y[1] = 0x1812145CDE0E8DCB; res_y[2] = 0xF4D6895A6375AA22; res_y[3] = 0x1A593C1711EEBCDE;

    if (fp2compare64((uint64_t*)A->x, res_x)!=0 || fp2compare64((uint64_t*)A->y, res_y)!=0) passed=0;

    if (passed==1) print_test("  Psi endomorphism tests .................................................................. PASSED");
    else { print_test("  Psi endomorphism tests ... FAILED"); print_test("\n"); return false; }
   
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
    res_x[0] = 0x1365D931AFEBC83E; res_x[1] = 0x1873BB71FF4FFF87; res_x[2] = 0x7BF9ACB5C770F61F; res_x[3] = 0x773EA05D9B4B0D62;
    res_y[0] = 0xCFFDD1A374E18F42; res_y[1] = 0x369B19C1F39C1A97; res_y[2] = 0x38B8E623E4E0049A; res_y[3] = 0x12435E356960429A;

    if (fp2compare64((uint64_t*)A->x, res_x)!=0 || fp2compare64((uint64_t*)A->y, res_y)!=0) passed=0;
    if (passed==1) print_test("  Phi endomorphism tests .................................................................. PASSED");
    else { print_test("  Phi endomorphism tests ... FAILED"); print_test("\n"); return false; }
    
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

    res_x[0] = 0x8F7033298B9CD5A4; res_x[1] = 0x6A60DF430E52E299; res_x[2] = 0x51D6EAFEEA829A8B; res_x[3] = 0x56F40C1CE3C3CD34;
    res_y[0] = 0x5B611ABE0387F840; res_y[1] = 0x59C6A5C83477F57C; res_y[2] = 0xF33C879AB74E2490; res_y[3] = 0x12C18E67FB2A3A9D;
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

    res_x[0] = 0x10EA7CF4F502CF5E; res_x[1] = 0x4FC1A5862ACAF69B; res_x[2] = 0x886D85328FB1E1A9; res_x[3] = 0x6F134E7E5129772A;
    res_y[0] = 0x35FFAD6E8F0681DC; res_y[1] = 0x681067510F99389E; res_y[2] = 0xA4BE7A70A1820895; res_y[3] = 0x34C0A821F434D672;
    if (fp2compare64((uint64_t*)A->x, res_x)!=0 || fp2compare64((uint64_t*)A->y, res_y)!=0) passed=0;
        
    if (passed==1) print_test("  Scalar multiplication tests ............................................................. PASSED");
    else { print_test("  Scalar multiplication tests ... FAILED"); print_test("\n"); return false; }
     
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
    
    if (passed==1) print_test("  Double scalar multiplication tests ...................................................... PASSED");
    else { print_test("  Double scalar multiplication tests ... FAILED"); print_test("\n"); return false; }
    }
    
    return OK;
}


bool ecc_run()
{
    bool OK = true;
    unsigned int n, i, digit=1;
    unsigned int cycles, cycles1, cycles2;
    point_t A, B, R;
    point_extedwards_t P, Q, RR, Table[16];
    f2elm_t t1;    
    uint64_t scalar[4];
        
    print_test("\n--------------------------------------------------------------------------------------------------------\n"); 
    print_test("Benchmarking FourQ's curve arithmetic \n"); 

    // Point doubling (twisted Edwards a=-1)
    eccset(A);
    point_setup(A, P);

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        for (i = 0; i < 10; i++) {
            eccdouble(P);
        }
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    print_bench("  Point doubling runs in ...                                       ", cycles/(BENCH_LOOPS*10));

    // Point addition (twisted Edwards a=-1)
    eccset(A);
    point_setup(A, P);
    ecccopy(P, Q);
    eccdouble(P);                       // P = 2P 

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        for (i = 0; i < 10; i++) {
            eccadd(P, Q, P);
        }
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    print_bench("  Point addition runs in ...                                       ", cycles/(BENCH_LOOPS*10));

    // Psi endomorphism
    eccset(A);
    point_setup(A, P);

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        for (i = 0; i < 10; i++) {
            ecc_psi(P);
        }
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    print_bench("  Psi mapping runs in ...                                          ", cycles/(BENCH_LOOPS*10));
   
    // Phi endomorphism
    eccset(A);
    point_setup(A, P);

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        for (i = 0; i < 10; i++) {
            ecc_phi(P);
        }
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    print_bench("  Phi mapping runs in ...                                          ", cycles/(BENCH_LOOPS*10));
   
    // Scalar decomposition
    {
    uint64_t scalars[4];
    random_scalar_test(scalar); 
    
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        for (i = 0; i < 10; i++) {
            decompose(scalar, scalars);
        }
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }

    print_bench("  Scalar decomposition runs in ...                                 ", cycles/(BENCH_LOOPS*10));
    }

    // Scalar recoding
    {
    unsigned int digits[72];
    uint128_t rand_scalars[4];
    unsigned char rand_bytes[4] = {1};

    random_scalar_test(scalar);
    randomize(scalar, rand_bytes, rand_scalars);

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        for (i = 0; i < 10; i++) {
            recode(rand_scalars, digits);
        }
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }

    print_bench("  Scalar recoding runs in ...                                      ", cycles/(BENCH_LOOPS*10));
    }
    
    // Precomputation
    eccset(A);
    point_setup(A, P);
    ecccopy(P, RR);

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        ecc_precomp(P, RR, Table);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }    

    print_bench("  Precomputation runs in ...                                       ", cycles/BENCH_LOOPS);

    // Table lookup
    eccset(A);
    point_setup(A, P);
    ecccopy(P, RR);
    ecc_precomp(P, RR, Table);

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        for (i = 0; i < 10; i++) {
            table_lookup_1x16(Table, Q, digit);
        }
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    } 
    
    print_bench("  Table lookup runs in ...                                         ", cycles/(BENCH_LOOPS*10));

    // Scalar multiplication
    random_scalar_test(scalar); 

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        eccset(A);
        eccset(R);
        cycles1 = cpucycles();
        ecc_mul_SCA_secure(A, R, (digit_t*)scalar, B, false);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }   
     
    print_bench("  Scalar multiplication (without clearing cofactor) runs in ...    ", cycles/BENCH_LOOPS);
    
    random_scalar_test(scalar); 

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        eccset(A);
        eccset(R);
        cycles1 = cpucycles();
        ecc_mul_SCA_secure(A, R, (digit_t*)scalar, B, true);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }   
     
    print_bench("  Scalar multiplication (including clearing cofactor) runs in ...  ", cycles/BENCH_LOOPS);

    {    
    point_t PP, QQ, RR;
    uint64_t k[4], l[4], kk[4];

    // Double scalar multiplication
    eccset(QQ); 
    eccset(PP);
    random_scalar_test(kk); 
    ecc_mul_SCA_secure(QQ, PP, (digit_t*)kk, QQ, false);
    
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {        
        random_scalar_test(k); 
        random_scalar_test(l);  
        cycles1 = cpucycles();
        ecc_mul_double((digit_t*)k, QQ, (digit_t*)l, RR);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    
    print_bench("  Double scalar mul runs in ...                                    ", cycles/BENCH_LOOPS);
    } 
    
    return OK;
} 


int main()
{
    clock_setup();
    gpio_setup();
    usart_setup(115200);
    rng_setup();

    *SCB_DEMCR = *SCB_DEMCR | 0x01000000;
    *DWT_CYCCNT = 0;               // reset the counter
    *DWT_CTRL = *DWT_CTRL | 1 ;    // enable the counter
    bool OK = true;

    OK = OK && ecc_test();         // Test FourQ's curve functions
    OK = OK && ecc_run();          // Benchmark FourQ's curve functions
    signal_host();
    
    return OK;
}
