/***********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: testing code for FourQ's field arithmetic 
************************************************************************************/    

#include "../FourQ_internal.h"
#include "test_extras.h"
#include "../stm32f4_wrapper.h"
#include <stdio.h>

static unsigned int *DWT_CYCCNT = (unsigned int*)0xE0001004;
static unsigned int *DWT_CTRL = (unsigned int*)0xE0001000;
static unsigned int *SCB_DEMCR = (unsigned int*)0xE000EDFC;

// Benchmark and test parameters  
#define BENCH_LOOPS       100       // Number of iterations per bench
#define TEST_LOOPS        100       // Number of iterations per test

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


bool fp2_test()
{ // Tests for the quadratic extension field arithmetic
    bool OK = true;
    int n, passed;
    f2elm_t a, b, c, d, e, f;

    print_test("\n--------------------------------------------------------------------------------------------------------\n"); 
    print_test("Testing quadratic extension field arithmetic over GF((2^127-1)^2): \n"); 

    // GF(p^2) multiplication using p = 2^127-1
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {    
        fp2random1271_test(a); fp2random1271_test(b); fp2random1271_test(c); 
        fp2random1271_test(d); fp2random1271_test(e); fp2random1271_test(f);
       
        fp2mul1271(a, b, d);                                             // d = a*b 
        mod1271(d[0]); mod1271(d[1]);
        fp2mul1271(b, a, e);                                             // e = b*a 
        mod1271(e[0]); mod1271(e[1]);
        if (fp2compare64((uint64_t*)d,(uint64_t*)e)!=0) { passed=0; break; }
        
        fp2mul1271(a, b, d); fp2mul1271(d, c, e);                        // e = (a*b)*c
        mod1271(e[0]); mod1271(e[1]);
        fp2mul1271(b, c, d); fp2mul1271(d, a, f);                        // f = a*(b*c)
        mod1271(f[0]); mod1271(f[1]);
        if (fp2compare64((uint64_t*)e,(uint64_t*)f)!=0) { passed=0; break; }
      
        fp2add1271(b, c, d); fp2mul1271(a, d, e);                        // e = a*(b+c)
        mod1271(e[0]); mod1271(e[1]);
        fp2mul1271(a, b, d); fp2mul1271(a, c, f); fp2add1271(d, f, f);   // f = a*b+a*c
        mod1271(f[0]); mod1271(f[1]);
        if (fp2compare64((uint64_t*)e,(uint64_t*)f)!=0) { passed=0; break; }
        
        fp2zero1271(b); b[0][0] = 1;
        fp2mul1271(a, b, d);                                             // d = a*1 
        mod1271(d[0]); mod1271(d[1]);                                                                     
        if (fp2compare64((uint64_t*)a,(uint64_t*)d)!=0) { passed=0; break; }
        
        fp2zero1271(b);
        fp2mul1271(a, b, d);                                             // d = a*0  
        mod1271(d[0]); mod1271(d[1]); 
        if (fp2compare64((uint64_t*)b,(uint64_t*)d)!=0) { passed=0; break; }
    }
    if (passed==1) print_test("  GF(p^2) multiplication tests .................................................................... PASSED");
    else { print_test("  GF(p^2) multiplication tests... FAILED"); print_test("\n"); return false; }
    
    // GF(p^2) squaring using p = 2^127-1
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fp2random1271_test(a); fp2random1271_test(b); fp2random1271_test(c);

        fp2sqr1271(a, b);                                            // b = a^2
        fp2mul1271(a, a, c);                                         // c = a*a 
        if (fp2compare64((uint64_t*)b,(uint64_t*)c)!=0) { passed=0; break; }
        
        fp2zero1271(a);
        fp2sqr1271(a, d);                                            // d = 0^2 
        if (fp2compare64((uint64_t*)a,(uint64_t*)d)!=0) { passed=0; break; }
    }
    if (passed==1) print_test("  GF(p^2) squaring tests........................................................................... PASSED");
    else { print_test("  GF(p^2) squaring tests... FAILED"); print_test("\n"); return false; }

    // GF(p^2) inversion using p = 2^127-1
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fp2random1271_test(a);  
        
        fp2zero1271(d); d[0][0] = 1;           
        fp2copy1271(a, b);                        
        fp2inv1271(a);                                
        fp2mul1271(a, b, c);                                        // c = a*a^-1 = 1
        mod1271(c[0]); mod1271(c[1]);
        if (fp2compare64((uint64_t*)c,(uint64_t*)d)!=0) { passed=0; break; }
    }
    if (passed==1) print_test("  GF(p^2) inversion tests.......................................................................... PASSED");
    else { print_test("  GF(p^2) inversion tests... FAILED"); print_test("\n"); return false; }
    
    return OK;
}


bool fp2_run()
{
    bool OK = true;
    int n, i;
    unsigned long long cycles, cycles1, cycles2;
    f2elm_t a, b, c, d, e, f;
        
    print_test("\n--------------------------------------------------------------------------------------------------------\n"); 
    print_test("Benchmarking quadratic extension field arithmetic over GF((2^127-1)^2): \n");

    // GF(p^2) addition using p = 2^127-1
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        fp2random1271_test(a); fp2random1271_test(b);

        cycles1 = cpucycles();
        for (i = 0; i < 1000; i++) {
            fp2add1271(a, b, c);
        }
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    print_bench("  GF(p^2) addition runs in ............... ", cycles/(BENCH_LOOPS*1000));

    // GF(p^2) subtraction using p = 2^127-1
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        fp2random1271_test(a); fp2random1271_test(b);  
         
        cycles1 = cpucycles();
        for (i = 0; i < 1000; i++) {
            fp2sub1271(a, b, c);
        }
        cycles2 = cpucycles();

        cycles = cycles+(cycles2-cycles1);
    }
    print_bench("  GF(p^2) subtraction runs in ............ ", cycles/(BENCH_LOOPS*1000));

    // GF(p^2) squaring using p = 2^127-1
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        fp2random1271_test(a); 
        
        cycles1 = cpucycles();
        for (i = 0; i < 1000; i++) {
            fp2sqr1271(a, b);
        }
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    print_bench("  GF(p^2) squaring runs in ............... ", cycles/(BENCH_LOOPS*1000));

    // GF(p^2) multiplication using p = 2^127-1
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        fp2random1271_test(a); fp2random1271_test(b); fp2random1271_test(c);

        cycles1 = cpucycles();
        for (i = 0; i < 1000; i++) {
            fp2mul1271(a, b, c);
        }
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    print_bench("  GF(p^2) multiplication runs in ......... ", cycles/(BENCH_LOOPS*1000));

    // GF(p^2) inversion using p = 2^127-1
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        fp2random1271_test(a); fp2random1271_test(b); fp2random1271_test(c);  

        cycles1 = cpucycles();
        for (i = 0; i < 100; i++) {
            fp2inv1271(a);
        }
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    print_bench("  GF(p^2) inversion runs in .............. ", cycles/(BENCH_LOOPS*100)); 

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

    OK = OK && fp2_test();     // Test quadratic extension field operations using p = 2^127-1
    OK = OK && fp2_run();      // Benchmark quadratic extension field operations using p = 2^127-1
    signal_host();
    
    return OK;
}
