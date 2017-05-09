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
    point_t A;
    point_extproj_t P;
    point_extproj_precomp_t Q;
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
        fp2copy1271((felm_t*)&PARAMETER_d, t1);
        fp2mul1271(t1, P->ta, t1);         // d*ta
        fp2add1271(t1, t1, t1);            // 2*d*ta
        fp2mul1271(t1, P->tb, Q->t2);      // 2*d*t
        fp2add1271(P->x, P->y, Q->xy);     // x+y    
        fp2sub1271(P->y, P->x, Q->yx);     // y-x
        fp2copy1271(P->z, Q->z2); 
        fp2add1271(Q->z2, Q->z2, Q->z2);   // 2*z
        eccadd(Q, P);                      // 2*P
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
    fp2copy1271((felm_t*)&PARAMETER_d, t1);
    fp2mul1271(t1, P->x, t1);              // d*x
    fp2add1271(t1, t1, t1);                // 2*d*x
    fp2mul1271(t1, P->y, Q->t2);           // 2*d*t
    fp2add1271(P->x, P->y, Q->xy);         // x+y    
    fp2sub1271(P->y, P->x, Q->yx);         // y-x
    fp2zero1271(Q->z2); *Q->z2[0] = 2;     // 2*z
    eccdouble(P);                          // P = 2P 

    for (n=0; n<TEST_LOOPS; n++)
    {
        eccadd(Q, P);                      // P = P+Q
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
   
#if (USE_ENDO == true)
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
    {        
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
    
    // Scalar decomposition and recoding
    {        
    uint64_t acc1, acc2, acc3, acc4, scalars[4];
    unsigned int digits[65], sign_masks[65];
    uint64_t k[4];
    int i;

    for (n=0; n<TEST_LOOPS*10; n++)
    {
        random_scalar_test(k);
        decompose(k, scalars);  
        fp2copy1271((felm_t*)scalars, (felm_t*)scalar);
        recode(scalars, digits, sign_masks); 

        acc1 = acc2 = acc3 = acc4 = 0; 

        for (i = 64; i >= 0; i--)
        {
            acc1 = 2*acc1; acc2 = 2*acc2; acc3 = 2*acc3; acc4 = 2*acc4; 
            if (sign_masks[i] == (unsigned int)-1) {
                acc1 += 1;
                acc2 += (digits[i] & 1);
                acc3 += ((digits[i] >> 1) & 1);
                acc4 += ((digits[i] >> 2) & 1);
            } else if (sign_masks[i] == 0) {
                acc1 -= 1;
                acc2 -= (digits[i] & 1);
                acc3 -= ((digits[i] >> 1) & 1);
                acc4 -= ((digits[i] >> 2) & 1);
            }
        }   
        if (scalar[0] != acc1 || scalar[1] != acc2  || scalar[2] != acc3 || scalar[3] != acc4) { passed=0; break; }
    }
    
    if (passed==1) print_test("  Recoding and decomposition tests ........................................................ PASSED");
    else { print_test("  Recoding and decomposition tests ... FAILED"); print_test("\n"); return false; }
    }
    }
#endif

    // Scalar multiplication
    eccset(A); 
    clear_cofactor = false;
    scalar[0] = 0x3AD457AB55456230; scalar[1] = 0x3A8B3C2C6FD86E0C; scalar[2] = 0x7E38F7C9CFBB9166; scalar[3] = 0x0028FD6CBDA458F0;
    
    for (n=0; n<TEST_LOOPS; n++)
    {
        scalar[1] = scalar[2];
        scalar[2] += scalar[0];

        ecc_mul(A, (digit_t*)scalar, A, clear_cofactor);
    }

    res_x[0] = 0x8F7033298B9CD5A4; res_x[1] = 0x6A60DF430E52E299; res_x[2] = 0x51D6EAFEEA829A8B; res_x[3] = 0x56F40C1CE3C3CD34;
    res_y[0] = 0x5B611ABE0387F840; res_y[1] = 0x59C6A5C83477F57C; res_y[2] = 0xF33C879AB74E2490; res_y[3] = 0x12C18E67FB2A3A9D;
    if (fp2compare64((uint64_t*)A->x, res_x)!=0 || fp2compare64((uint64_t*)A->y, res_y)!=0) passed=0;

        
    eccset(A); 
    clear_cofactor = true;
    scalar[0] = 0x3AD457AB55456230; scalar[1] = 0x3A8B3C2C6FD86E0C; scalar[2] = 0x7E38F7C9CFBB9166; scalar[3] = 0x0028FD6CBDA458F0;
    
    for (n=0; n<TEST_LOOPS; n++)
    {
        scalar[1] = scalar[2];
        scalar[2] += scalar[0];

        ecc_mul(A, (digit_t*)scalar, A, clear_cofactor);
    }

    res_x[0] = 0x10EA7CF4F502CF5E; res_x[1] = 0x4FC1A5862ACAF69B; res_x[2] = 0x886D85328FB1E1A9; res_x[3] = 0x6F134E7E5129772A;
    res_y[0] = 0x35FFAD6E8F0681DC; res_y[1] = 0x681067510F99389E; res_y[2] = 0xA4BE7A70A1820895; res_y[3] = 0x34C0A821F434D672;
    if (fp2compare64((uint64_t*)A->x, res_x)!=0 || fp2compare64((uint64_t*)A->y, res_y)!=0) passed=0;
        
    if (passed==1) print_test("  Scalar multiplication tests ............................................................. PASSED");
    else { print_test("  Scalar multiplication tests ... FAILED"); print_test("\n"); return false; }
     
    {    
    point_t AA, B, C; 
    unsigned int j, w, v, e, d;
    uint64_t k[4];
    unsigned int digits_fixed[NBITS_ORDER_PLUS_ONE+(W_FIXEDBASE*V_FIXEDBASE)-1] = {0};
        
    // Scalar recoding using the mLSB-set representation    
    w = W_FIXEDBASE;
    v = V_FIXEDBASE;
    e = E_FIXEDBASE;     
    d = D_FIXEDBASE;      
          
    for (n=0; n<TEST_LOOPS; n++) 
    {
        random_scalar_test(scalar);
        modulo_order((digit_t*)scalar, (digit_t*)scalar);    // k = k mod (order) 
        conversion_to_odd((digit_t*)scalar, (digit_t*)k);                       
        for (j = 0; j < NWORDS64_ORDER; j++) scalar[j] = k[j];
        mLSB_set_recode(k, digits_fixed);
            
        if (verify_mLSB_recoding(scalar, (int*)digits_fixed)==false) { passed=0; break; }
    }
    
    if (passed==1) print_test("  mLSB-set recoding tests ................................................................. PASSED");
    else { print_test("  mLSB-set recoding tests ... FAILED"); print_test("\n"); return false; }

    // Fixed-base scalar multiplication
    eccset(AA); 
    
    for (n=0; n<TEST_LOOPS; n++)
    {
        random_scalar_test(scalar); 
        ecc_mul_fixed((digit_t*)scalar, B);
        ecc_mul(AA, (digit_t*)scalar, C, false);
        
        if (fp2compare64((uint64_t*)B->x,(uint64_t*)C->x)!=0 || fp2compare64((uint64_t*)B->y,(uint64_t*)C->y)!=0) { passed=0; break; }
    }

    if (passed==1) print_test("  Fixed-base scalar multiplication tests .................................................. PASSED");
    else { print_test("  Fixed-base scalar multiplication tests ... FAILED"); print_test("\n"); return false; }
    }
          
    {    
    point_t PP, QQ, RR, UU, TT; 
    point_extproj_precomp_t AA;
    point_extproj_t BB;
    uint64_t k[4], l[4], kk[4];

    // Double scalar multiplication
    eccset(QQ); 
    eccset(PP);
    
    for (n=0; n<TEST_LOOPS; n++)
    {
        random_scalar_test(kk); 
        ecc_mul(QQ, (digit_t*)kk, QQ, false);
        random_scalar_test(k); 
        random_scalar_test(l); 
        ecc_mul_double((digit_t*)k, QQ, (digit_t*)l, RR);
        ecc_mul(PP, (digit_t*)k, UU, false);
        ecc_mul(QQ, (digit_t*)l, TT, false);
    
        fp2add1271(UU->x, UU->y, AA->xy); 
        fp2sub1271(UU->y, UU->x, AA->yx); 
        fp2mul1271(UU->x, UU->y, AA->t2);    
        fp2add1271(AA->t2, AA->t2, AA->t2); 
        fp2mul1271(AA->t2, (felm_t*)&PARAMETER_d, AA->t2); 
        fp2zero1271(AA->z2); AA->z2[0][0] = 2;
        point_setup(TT, BB);             

        eccadd(AA, BB);
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
    unsigned int n;
    unsigned int cycles, cycles1, cycles2;
    point_t A, B;
    point_extproj_t P;
    point_extproj_precomp_t Q, Table[8];
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
        eccdouble(P);
        eccdouble(P);
        eccdouble(P);
        eccdouble(P);
        eccdouble(P);
        eccdouble(P);
        eccdouble(P);
        eccdouble(P);
        eccdouble(P);
        eccdouble(P);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    print_bench("  Point doubling runs in ...                                       ", cycles/(BENCH_LOOPS*10));

    // Point addition (twisted Edwards a=-1)
    eccset(A);
    point_setup(A, P);
    fp2copy1271((felm_t*)&PARAMETER_d, t1);
    fp2mul1271(t1, P->x, t1);              // d*x
    fp2add1271(t1, t1, t1);                // 2*d*x
    fp2mul1271(t1, P->y, Q->t2);           // 2*d*t
    fp2add1271(P->x, P->y, Q->xy);         // x+y    
    fp2sub1271(P->y, P->x, Q->yx);         // y-x
    fp2zero1271(Q->z2); *Q->z2[0] = 2;     // 2*z
    eccdouble(P);                          // P = 2P 

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        eccadd(Q, P);
        eccadd(Q, P);
        eccadd(Q, P);
        eccadd(Q, P);
        eccadd(Q, P);
        eccadd(Q, P);
        eccadd(Q, P);
        eccadd(Q, P);
        eccadd(Q, P);
        eccadd(Q, P);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    print_bench("  Point addition runs in ...                                       ", cycles/(BENCH_LOOPS*10));
   
#if (USE_ENDO == true)
    // Psi endomorphism
    eccset(A);
    point_setup(A, P);

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        ecc_psi(P);
        ecc_psi(P);
        ecc_psi(P);
        ecc_psi(P);
        ecc_psi(P);
        ecc_psi(P);
        ecc_psi(P);
        ecc_psi(P);
        ecc_psi(P);
        ecc_psi(P);
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
        ecc_phi(P);
        ecc_phi(P);
        ecc_phi(P);
        ecc_phi(P);
        ecc_phi(P);
        ecc_phi(P);
        ecc_phi(P);
        ecc_phi(P);
        ecc_phi(P);
        ecc_phi(P);
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
        decompose(scalar, scalars);
        decompose(scalar, scalars);
        decompose(scalar, scalars);
        decompose(scalar, scalars);
        decompose(scalar, scalars);
        decompose(scalar, scalars);
        decompose(scalar, scalars);
        decompose(scalar, scalars);
        decompose(scalar, scalars);
        decompose(scalar, scalars);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    
    print_bench("  Scalar decomposition runs in ...                                 ", cycles/(BENCH_LOOPS*10));
    }

    // Scalar recoding
    {
    unsigned int digits[65], sign_masks[65];
    random_scalar_test(scalar); 

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        recode(scalar, digits, sign_masks);
        recode(scalar, digits, sign_masks);
        recode(scalar, digits, sign_masks);
        recode(scalar, digits, sign_masks);
        recode(scalar, digits, sign_masks);
        recode(scalar, digits, sign_masks);
        recode(scalar, digits, sign_masks);
        recode(scalar, digits, sign_masks);
        recode(scalar, digits, sign_masks);
        recode(scalar, digits, sign_masks);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    
    print_bench("  Scalar recoding runs in ...                                      ", cycles/(BENCH_LOOPS*10));
    }
#endif

    // Precomputation
    eccset(A);
    point_setup(A, P);

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        ecc_precomp(P, Table);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    
    print_bench("  Precomputation runs in ...                                       ", cycles/BENCH_LOOPS);

    // Table lookup
    eccset(A);
    point_setup(A, P);
    ecc_precomp(P, Table);

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        table_lookup_1x8(Table, Q, 0, 0);
        table_lookup_1x8(Table, Q, 1, (unsigned int)-1);
        table_lookup_1x8(Table, Q, 2, 0);
        table_lookup_1x8(Table, Q, 3, (unsigned int)-1);
        table_lookup_1x8(Table, Q, 4, 0);
        table_lookup_1x8(Table, Q, 5, (unsigned int)-1);
        table_lookup_1x8(Table, Q, 6, 0);
        table_lookup_1x8(Table, Q, 7, (unsigned int)-1);
        table_lookup_1x8(Table, Q, 0, 0);
        table_lookup_1x8(Table, Q, 1, (unsigned int)-1);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    } 
    
    print_bench("  Table lookup runs in ...                                         ", cycles/(BENCH_LOOPS*10));

    // Scalar multiplication
    random_scalar_test(scalar); 

    for (n=0; n<BENCH_LOOPS; n++)
    {
        eccset(A);
        ecc_mul(A, (digit_t*)scalar, B, false);
    }
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        eccset(A);
        cycles1 = cpucycles();
        ecc_mul(A, (digit_t*)scalar, B, false);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    
    print_bench("  Scalar multiplication (without clearing cofactor) runs in ...    ", cycles/BENCH_LOOPS);
    
    random_scalar_test(scalar); 

    for (n=0; n<BENCH_LOOPS; n++)
    {
        eccset(A);
        ecc_mul(A, (digit_t*)scalar, B, true);
    }
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        eccset(A);
        cycles1 = cpucycles();
        ecc_mul(A, (digit_t*)scalar, B, true);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    
    print_bench("  Scalar multiplication (including clearing cofactor) runs in ...  ", cycles/BENCH_LOOPS);
     
    {      
    point_precomp_t T;
    unsigned int digits_fixed[256+(W_FIXEDBASE*V_FIXEDBASE)-1] = {0};

    // Reduction modulo the order 
    random_scalar_test(scalar); 

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        modulo_order((digit_t*)scalar, (digit_t*)scalar);
        modulo_order((digit_t*)scalar, (digit_t*)scalar);
        modulo_order((digit_t*)scalar, (digit_t*)scalar);
        modulo_order((digit_t*)scalar, (digit_t*)scalar);
        modulo_order((digit_t*)scalar, (digit_t*)scalar);
        modulo_order((digit_t*)scalar, (digit_t*)scalar);
        modulo_order((digit_t*)scalar, (digit_t*)scalar);
        modulo_order((digit_t*)scalar, (digit_t*)scalar);
        modulo_order((digit_t*)scalar, (digit_t*)scalar);
        modulo_order((digit_t*)scalar, (digit_t*)scalar);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    
    print_bench("  Reduction modulo the order runs in ...                           ", cycles/(BENCH_LOOPS*10));

    // Scalar recoding using the mLSB-set representation 
    random_scalar_test(scalar); 

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        mLSB_set_recode(scalar, digits_fixed);
        mLSB_set_recode(scalar, digits_fixed);
        mLSB_set_recode(scalar, digits_fixed);
        mLSB_set_recode(scalar, digits_fixed);
        mLSB_set_recode(scalar, digits_fixed);
        mLSB_set_recode(scalar, digits_fixed);
        mLSB_set_recode(scalar, digits_fixed);
        mLSB_set_recode(scalar, digits_fixed);
        mLSB_set_recode(scalar, digits_fixed);
        mLSB_set_recode(scalar, digits_fixed);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    
    print_bench("  Fixed-base recoding runs in ...                                  ", cycles/(BENCH_LOOPS*10));

    // Table lookup for fixed-base scalar multiplication
    eccset(A);

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        table_lookup_fixed_base((point_precomp_t*)&FIXED_BASE_TABLE, T, 1, 0);
        table_lookup_fixed_base((point_precomp_t*)&FIXED_BASE_TABLE, T, 2, (unsigned int)-1);
        table_lookup_fixed_base((point_precomp_t*)&FIXED_BASE_TABLE, T, 1, 0);
        table_lookup_fixed_base((point_precomp_t*)&FIXED_BASE_TABLE, T, 2, (unsigned int)-1);
        table_lookup_fixed_base((point_precomp_t*)&FIXED_BASE_TABLE, T, 1, 0);
        table_lookup_fixed_base((point_precomp_t*)&FIXED_BASE_TABLE, T, 2, (unsigned int)-1);
        table_lookup_fixed_base((point_precomp_t*)&FIXED_BASE_TABLE, T, 1, 0);
        table_lookup_fixed_base((point_precomp_t*)&FIXED_BASE_TABLE, T, 2, (unsigned int)-1);
        table_lookup_fixed_base((point_precomp_t*)&FIXED_BASE_TABLE, T, 1, 0);
        table_lookup_fixed_base((point_precomp_t*)&FIXED_BASE_TABLE, T, 2, (unsigned int)-1);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    
    print_bench("  Fixed-base table lookup runs in ...                              ", cycles/(BENCH_LOOPS*10));

    // Fixed-base scalar multiplication  
    eccset(A);

    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {        
        random_scalar_test(scalar); 
        cycles1 = cpucycles();
        ecc_mul_fixed((digit_t*)scalar, B);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    } 
    
    print_bench("  Fixed-base scalar mul runs in ...                                ", cycles/BENCH_LOOPS);
    } 
              
    {    
    point_t PP, QQ, RR; 
    uint64_t k[4], l[4], kk[4];

    // Double scalar multiplication
    eccset(QQ); 
    eccset(PP);
    random_scalar_test(kk); 
    ecc_mul(QQ, (digit_t*)kk, QQ, false);
    
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
