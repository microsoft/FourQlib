/***********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: testing code for FourQ's field arithmetic 
************************************************************************************/   

#include "../FourQ_internal.h"
#include "test_extras.h"
#include <stdio.h>


// Benchmark and test parameters 
#define BENCH_LOOPS       10000      // Number of iterations per bench
#define SHORT_BENCH_LOOPS 1000       // Number of iterations per bench (for expensive operations)
#define TEST_LOOPS        1000       // Number of iterations per test


bool fp_test()
{ // Tests for the quadratic extension field arithmetic
    bool OK = true;
    int n, i, passed;
    f2elm_t b, c, d, e, f;
    velm_t va, vb, vc, vd, ve, vf, vz = {0}, vo = {0};
    v2elm_t vv;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Testing field arithmetic over GF(2^127-1): \n\n"); 
    
    // Field multiplication using p = 2^127-1
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {    
        vrandom1271_test(va); vrandom1271_test(vb); vrandom1271_test(vc); 
        vrandom1271_test(vd); vrandom1271_test(ve); vrandom1271_test(vf);
       
        vmul1271(va, vb, vd);                                                        // d = a*b 
        vmod1271(vd, vd); from_v_to_v2(vd, vz, vv); from_ext_to_std(vv, d);
        vmul1271(vb, va, ve);                                                        // e = b*a 
        vmod1271(ve, ve); from_v_to_v2(ve, vz, vv); from_ext_to_std(vv, e);
        if (fpcompare64((uint64_t*)d,(uint64_t*)e)!=0) { passed=0; break; }

        vmul1271(va, vb, vd); vmul1271(vd, vc, ve);                                  // e = (a*b)*c
        vmod1271(ve, ve); from_v_to_v2(ve, vz, vv); from_ext_to_std(vv, e);
        vmul1271(vb, vc, vd); vmul1271(vd, va, vf);                                  // f = a*(b*c)
        vmod1271(vf, vf); from_v_to_v2(vf, vz, vv); from_ext_to_std(vv, f);
        if (fpcompare64((uint64_t*)e,(uint64_t*)f)!=0) { passed=0; break; }
      
        vadd1271(vb, vc, vd); vmul1271(va, vd, ve);                                  // e = a*(b+c)
        vmod1271(ve, ve); from_v_to_v2(ve, vz, vv); from_ext_to_std(vv, e);
        vmul1271(va, vb, vd); vmul1271(va, vc, vf); vadd1271(vd, vf, vf);            // f = a*b+a*c
        vmod1271(vf, vf); from_v_to_v2(vf, vz, vv); from_ext_to_std(vv, f);
        if (fpcompare64((uint64_t*)e,(uint64_t*)f)!=0) { passed=0; break; }
        
        vo[0] = 1;
        vmul1271(va, vo, vd);                                                        // d = a*1 
        if (fpcompare64((uint64_t*)va,(uint64_t*)vd)!=0) { passed=0; break; }
        
        vmul1271(va, vz, vd); from_v_to_v2(vd, vz, vv); from_ext_to_std(vv, e);      // d = a*0 
        if (fpcompare64((uint64_t*)e,(uint64_t*)vz)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  Field multiplication tests ...................................................................... PASSED");
    else { printf("  Field multiplication tests... FAILED"); printf("\n"); return false; }
    printf("\n");
    
    // Field squaring using p = 2^127-1
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        vrandom1271_test(va); vrandom1271_test(vb); vrandom1271_test(vc);
        
        vsqr1271(va, vb);                                                            // b = a^2
        vmod1271(vb, vb); from_v_to_v2(vb, vz, vv); from_ext_to_std(vv, b);
        vmul1271(va, va, vc);                                                        // c = a*a 
        vmod1271(vc, vc); from_v_to_v2(vc, vz, vv); from_ext_to_std(vv, c);
        if (fpcompare64((uint64_t*)b,(uint64_t*)c)!=0) { passed=0; break; }
        
        vsqr1271(vz, vd); from_v_to_v2(vd, vz, vv); from_ext_to_std(vv, e);          // d = 0^2 
        if (fpcompare64((uint64_t*)e,(uint64_t*)vz)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  Field squaring tests............................................................................. PASSED");
    else { printf("  Field squaring tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // Field inversion using p = 2^127-1
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        vrandom1271_test(va);  
        
        vo[0] = 1;           
        for (i = 0; i < VWORDS_FIELD; i++) { vb[i] = va[i]; }                        
        vinv1271(va);                                
        vmul1271(va, vb, vc);                                                        // c = a*a^-1 = 1
        vmod1271(vc, vc); from_v_to_v2(vc, vz, vv); from_ext_to_std(vv, c);
        if (fpcompare64((uint64_t*)c,(uint64_t*)vo)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  Field inversion tests............................................................................ PASSED");
    else { printf("  Field inversion tests... FAILED"); printf("\n"); return false; }
    printf("\n");
    
    return OK;
}


bool fp2_test()
{ // Tests for the quadratic extension field arithmetic
    bool OK = true;
    int n, passed;
    f2elm_t a, b, c, d, e, f;
    v2elm_t va, vb, vc, vd, ve, vf;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Testing quadratic extension field arithmetic over GF((2^127-1)^2): \n\n"); 
    
    // GF(p^2) multiplication using p = 2^127-1
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {    
        v2random1271_test(va); v2random1271_test(vb); v2random1271_test(vc); 
        v2random1271_test(vd); v2random1271_test(ve); v2random1271_test(vf);
       
        v2mul1271(va, vb, vd);                                                // d = a*b 
        v2mod1271(vd, vd); from_ext_to_std(vd, d);
        v2mul1271(vb, va, ve);                                                // e = b*a 
        v2mod1271(ve, ve); from_ext_to_std(ve, e);
        if (fp2compare64((uint64_t*)d,(uint64_t*)e)!=0) { passed=0; break; }
            
        v2mul1271(va, vb, vd); v2mul1271(vd, vc, ve);                         // e = (a*b)*c
        v2mod1271(ve, ve); from_ext_to_std(ve, e);
        v2mul1271(vb, vc, vd); v2mul1271(vd, va, vf);                         // f = a*(b*c)
        v2mod1271(vf, vf); from_ext_to_std(vf, f);
        if (fp2compare64((uint64_t*)e,(uint64_t*)f)!=0) { passed=0; break; }
      
        v2add1271(vb, vc, vd); v2mul1271(va, vd, ve);                         // e = a*(b+c)
        v2mod1271(ve, ve); from_ext_to_std(ve, e);
        v2mul1271(va, vb, vd); v2mul1271(va, vc, vf); v2add1271(vd, vf, vf);  // f = a*b+a*c
        v2mod1271(vf, vf); from_ext_to_std(vf, f);
        if (fp2compare64((uint64_t*)e,(uint64_t*)f)!=0) { passed=0; break; }
        
        v2zero1271(vb); vb[0] = 1;
        v2mul1271(va, vb, vd);                                                // d = a*1  
        v2mod1271(va, va); from_ext_to_std(va, a); 
        v2mod1271(vd, vd); from_ext_to_std(vd, d);                                                                    
        if (fp2compare64((uint64_t*)a,(uint64_t*)d)!=0) { passed=0; break; }
        
        v2zero1271(vb);
        from_ext_to_std(vb, b);
        v2mul1271(va, vb, vd);                                                // d = a*0  
        v2mod1271(vd, vd); from_ext_to_std(vd, d); 
        if (fp2compare64((uint64_t*)b,(uint64_t*)d)!=0) { passed=0; break; } 
    }
    if (passed==1) printf("  GF(p^2) multiplication tests .................................................................... PASSED");
    else { printf("  GF(p^2) multiplication tests... FAILED"); printf("\n"); return false; }
    printf("\n");
    
    // GF(p^2) squaring using p = 2^127-1
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        v2random1271_test(va); v2random1271_test(vb); v2random1271_test(vc);
            
        v2sqr1271(va, vb);                                            // b = a^2
        v2mod1271(vb, vb); from_ext_to_std(vb, b);
        v2mul1271(va, va, vc);                                        // c = a*a 
        v2mod1271(vc, vc); from_ext_to_std(vc, c);
        if (fp2compare64((uint64_t*)b,(uint64_t*)c)!=0) { passed=0; break; }
        
        v2zero1271(va);
        from_ext_to_std(va, a);
        v2sqr1271(va, vd);                                            // d = 0^2 
        v2mod1271(vd, vd); from_ext_to_std(vd, d);
        if (fp2compare64((uint64_t*)a,(uint64_t*)d)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p^2) squaring tests........................................................................... PASSED");
    else { printf("  GF(p^2) squaring tests... FAILED"); printf("\n"); return false; }
    printf("\n");
    
    // GF(p^2) inversion using p = 2^127-1
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        v2random1271_test(va);         
        
        v2zero1271(vd); vd[0] = 1;  
        from_ext_to_std(vd, d);         
        v2copy1271(va, vb);                        
        v2inv1271(va);                                
        v2mul1271(va, vb, vc);                                        // c = a*a^-1 = 1
        v2mod1271(vc, vc); from_ext_to_std(vc, c);
        if (fp2compare64((uint64_t*)c,(uint64_t*)d)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p^2) inversion tests.......................................................................... PASSED");
    else { printf("  GF(p^2) inversion tests... FAILED"); printf("\n"); return false; }
    printf("\n");
    
    return OK;
}


bool fp_run()
{
    bool OK = true;
    int n, i;
    unsigned long long nsec, nsec1, nsec2;
    velm_t va, vb, vc;
        
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Benchmarking field arithmetic over GF(2^127-1): \n\n"); 

    // Field squaring over p = 2^127-1
    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        vrandom1271_test(va);

        nsec1 = cpu_nseconds();
        for (i = 0; i < 1000; i++) {
            vsqr1271(va, vc);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    printf("  Field squaring runs in ................. %8lld nsec", nsec/(BENCH_LOOPS*1000));
    printf("\n"); 

    // Field multiplication over p = 2^127-1
    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        vrandom1271_test(va); vrandom1271_test(vb); vrandom1271_test(vc); 

        nsec1 = cpu_nseconds();
        for (i = 0; i < 1000; i++) {
            vmul1271(va, vb, vc);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    printf("  Field multiplication runs in ........... %8lld nsec", nsec/(BENCH_LOOPS*1000));
    printf("\n");

    // Field inversion using p = 2^127-1
    nsec = 0;
    for (n=0; n<SHORT_BENCH_LOOPS; n++)
    {
        vrandom1271_test(va);

        nsec1 = cpu_nseconds();
        for (i = 0; i < 100; i++) {
            vinv1271(va);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    printf("  Field inversion runs in ................ %8lld nsec", nsec/(SHORT_BENCH_LOOPS*100));
    printf("\n");
    
    return OK;
}


bool fp2_run()
{
    bool OK = true;
    int n, i;
    unsigned long long nsec, nsec1, nsec2;
    v2elm_t va, vb, vc, vd, ve, vf, vg;
        
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Benchmarking quadratic extension field arithmetic over GF((2^127-1)^2): \n\n"); 

    // GF(p^2) addition using p = 2^127-1
    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        v2random1271_test(va); v2random1271_test(vb); v2random1271_test(vc); 

        nsec1 = cpu_nseconds();
        for (i = 0; i < 1000; i++) {
            v2add1271(va, vb, vc);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    printf("  GF(p^2) addition runs in ............... %8lld nsec", nsec/(BENCH_LOOPS*1000));
    printf("\n"); 

    // GF(p^2) subtraction using p = 2^127-1
    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        v2random1271_test(va); v2random1271_test(vb); v2random1271_test(vc); 

        nsec1 = cpu_nseconds();
        for (i = 0; i < 1000; i++) {
            v2sub1271(va, vb, vc);
        }
        nsec2 = cpu_nseconds();

        nsec = nsec+(nsec2-nsec1);
    }
    printf("  GF(p^2) subtraction runs in ............ %8lld nsec", nsec/(BENCH_LOOPS*1000));
    printf("\n"); 

    // GF(p^2) squaring using p = 2^127-1
    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        v2random1271_test(va); v2random1271_test(vb); 

        nsec1 = cpu_nseconds();
        for (i = 0; i < 1000; i++) {
            v2sqr1271(va, vb);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    printf("  GF(p^2) squaring runs in ............... %8lld nsec", nsec/(BENCH_LOOPS*1000));
    printf("\n");

    // GF(p^2) multiplication using p = 2^127-1
    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        v2random1271_test(va); v2random1271_test(vb); v2random1271_test(vc); 

        nsec1 = cpu_nseconds();
        for (i = 0; i < 1000; i++) {
            v2mul1271(va, vb, vc);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    printf("  GF(p^2) multiplication runs in ......... %8lld nsec", nsec/(BENCH_LOOPS*1000));
    printf("\n");
    
#if defined(MIX_ARM_NEON)

    // GF(p^2) squaring/addition/subtraction using p = 2^127-1
    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        v2random1271_test(va); v2random1271_test(vc); v2random1271_test(vd); 

        nsec1 = cpu_nseconds();
        for (i = 0; i < 1000; i++) {
            v2sqraddsub1271(va, vb, vc, vd, ve, vf);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    printf("  GF(p^2) squaring/add/sub runs in ....... %8lld nsec", nsec/(BENCH_LOOPS*1000));
    printf("\n");

    // GF(p^2) multiplication/subtraction using p = 2^127-1
    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        v2random1271_test(va); v2random1271_test(vb); v2random1271_test(vd); v2random1271_test(ve); 

        nsec1 = cpu_nseconds();
        for (i = 0; i < 1000; i++) {
            v2mulsub1271(va, vb, vc, vd, ve, vf);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    printf("  GF(p^2) mul/subtraction runs in ........ %8lld nsec", nsec/(BENCH_LOOPS*1000));
    printf("\n");

    // GF(p^2) multiplication/addition/subtraction using p = 2^127-1
    nsec = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        v2random1271_test(va); v2random1271_test(vb); v2random1271_test(vc); v2random1271_test(vd); v2random1271_test(ve);

        nsec1 = cpu_nseconds();
        for (i = 0; i < 1000; i++) {
            v2muladdsub1271(va, vb, vc, vd, ve, vf, vg);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    printf("  GF(p^2) mul/add/sub runs in ............ %8lld nsec", nsec/(BENCH_LOOPS*1000));
    printf("\n");

#endif

    // GF(p^2) inversion using p = 2^127-1
    nsec = 0;
    for (n=0; n<SHORT_BENCH_LOOPS; n++)
    {
        v2random1271_test(va); v2random1271_test(vb); v2random1271_test(vc); 
        v2random1271_test(vd); v2random1271_test(ve); v2random1271_test(vf); 

        nsec1 = cpu_nseconds();
        for (i = 0; i < 100; i++) {
            v2inv1271(va);
        }
        nsec2 = cpu_nseconds();
        nsec = nsec+(nsec2-nsec1);
    }
    printf("  GF(p^2) inversion runs in .............. %8lld nsec", nsec/(SHORT_BENCH_LOOPS*100));
    printf("\n");
    
    return OK;
}


int main()
{
    bool OK = true;
    
    OK = OK && fp_test();      // Test field operations using p = 2^127-1
    OK = OK && fp_run();       // Benchmark field operations using p = 2^127-1
    OK = OK && fp2_test();     // Test quadratic extension field operations using p = 2^127-1
    OK = OK && fp2_run();      // Benchmark quadratic extension field operations using p = 2^127-1
    
    return OK;
}
