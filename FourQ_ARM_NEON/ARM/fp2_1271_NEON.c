/***********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: arithmetic over GF(p) and GF(p^2) using ARM/NEON inlined assembly
************************************************************************************/ 

#include "../FourQ_internal.h"


void v2mul1271_a(uint32_t* a, uint32_t* b, uint32_t* c)  
{ // Multiplication over GF((2^127-1)^2) using NEON instructions
  // Operation: c [r3] = a [r1] * b [r2]
  // Representation: 23/23/26/26/26/26/26/26/26/26-bit

asm volatile(
    // q0-q2[0] <- A0|A1, q3-q5[0] <- B0|B1
    "vld1.8       {d0,d1}, [%0]!            \n\t"
    "vld1.8       {d2,d3}, [%0]!            \n\t"
    "vld1.8       {d4}, [%0]!               \n\t"
#if defined(INTERLEAVE)
    "vshl.i32     q6, q0, #3                \n\t"    // q6-q7 = 8*(a0->a8)  
    "vpush        {q4-q7}                   \n\t"  
    "vld1.8       {d6,d7}, [%1]!          	\n\t"
    "vshl.i32     d5, d4, #3                \n\t"    // q2[1] = 8*(a4|a9)   
    "vld1.8       {d8,d9}, [%1]!          	\n\t"
    "vshl.i32     q7, q1, #3                \n\t" 
    "vld1.8       {d10}, [%1]!              \n\t"
#else
    "vshl.i32     q6, q0, #3                \n\t"    // q6-q7 = 8*(a0->a8)  
    "vld1.8       {d6,d7}, [%1]!          	\n\t"
    "vpush        {q4-q7}                   \n\t"  
    "vld1.8       {d8,d9}, [%1]!          	\n\t"
    "vld1.8       {d10}, [%1]!              \n\t"
    "vshl.i32     d5, d4, #3                \n\t"    // q2[1] = 8*(a4|a9)   
    "vshl.i32     q7, q1, #3                \n\t" 
#endif

    // q11-q15 <- A0*B0|A0*B1, q8-q10 <- A1*B0|A1*B1 
    "vmull.s32    q8, d6, d0[1]             \n\t"    // q8 = a5*b0 | a5*b5
    "vmlal.s32    q8, d10, d13[1]           \n\t"    // q8 = a5*b0 + 8*(a6*b4) | a5*b5 + 8*(a6*b9)
    "vmlal.s32    q8, d7, d5[1]             \n\t"    // q8 = a5*b0 + 8*(a6*b4 + a9*b1) | a5*b5 + 8*(a6*b9 + a9*b6)
    "vmlal.s32    q8, d9, d14[1]            \n\t"    // q8 = a5*b0 + 8*(a6*b4 + a9*b1 + a7*b3) | a5*b5 + 8*(a6*b9 + a9*b6 + a7*b8)
    "vmlal.s32    q8, d8, d15[1]            \n\t"    // q8 = a5*b0 + 8*(a6*b4 + a9*b1 + a7*b3 + a8*b2) | a5*b5 + 8*(a6*b9 + a9*b6 + a7*b8 + a8*b7)
    "vmull.s32    q9, d7, d0[1]             \n\t"    // q9 = a5*b1 | a5*b6
    "vmlal.s32    q9, d6, d1[1]             \n\t"    // q9 = a5*b1 + a6*b0 | a5*b6 + a6*b5
    "vmlal.s32    q9, d10, d14[1]           \n\t"    // q9 = a5*b1 + a6*b0 + 8*(a7*b4) | a5*b6 + a6*b5 + 8*(a7*b9)
    "vmlal.s32    q9, d8, d5[1]             \n\t"    // q9 = a5*b1 + a6*b0 + 8*(a7*b4 + a9*b2) | a5*b6 + a6*b5 + 8*(a7*b9 + a9*b7)
    "vmlal.s32    q9, d9, d15[1]            \n\t"    // q9 = a5*b1 + a6*b0 + 8*(a7*b4 + a9*b2 + a8*b3) | a5*b6 + a6*b5 + 8*(a7*b9 + a9*b7 + a8*b8)
    "vmull.s32    q10, d8, d0[1]            \n\t"    // q10 = a5*b2 | a5*b7
    "vmlal.s32    q10, d6, d2[1]            \n\t"    // q10 = a5*b2 + a7*b0 | a5*b7 + a7*b5
    "vmlal.s32    q10, d7, d1[1]            \n\t"    // q10 = a5*b2 + a7*b0 + a6*b1 | a5*b7 + a7*b5 + a6*b6
    "vmlal.s32    q10, d10, d15[1]          \n\t"    // q10 = a5*b2 + a7*b0 + a6*b1 + 8*(a8*b4) | a5*b7 + a7*b5 + a6*b6 + 8*(a8*b9)
    "vmlal.s32    q10, d9, d5[1]            \n\t"    // q10 = a5*b2 + a7*b0 + a6*b1 + 8*(a8*b4 + a9*b3) | a5*b7 + a7*b5 + a6*b6 + 8*(a8*b9 + a9*b8)
    
    "vmull.s32    q11, d6, d0[0]             \n\t"    // q11 = a0*b0 | a0*b5
    "vmlal.s32    q11, d10, d13[0]           \n\t"    // q11 = a0*b0 + 8*(a1*b4) | a0*b5 + 8*(a1*b9)
    "vmlal.s32    q11, d7, d5[0]             \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1) | a0*b5 + 8*(a1*b9 + a4*b6)
    "vmlal.s32    q11, d9, d14[0]            \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3) | a0*b5 + 8*(a1*b9 + a4*b6 + a2*b8)
    "vmlal.s32    q11, d8, d15[0]            \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3 + a3*b2) | a0*b5 + 8*(a1*b9 + a4*b6 + a2*b8 + a3*b7)
    "vmull.s32    q12, d7, d0[0]             \n\t"    // q12 = a0*b1 | a0*b6
    "vmlal.s32    q12, d6, d1[0]             \n\t"    // q12 = a0*b1 + a1*b0 | a0*b6 + a1*b5
    "vmlal.s32    q12, d10, d14[0]           \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4) | a0*b6 + a1*b5 + 8*(a2*b9)
    "vmlal.s32    q12, d8, d5[0]             \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2) | a0*b6 + a1*b5 + 8*(a2*b9 + a4*b7)
    "vmlal.s32    q12, d9, d15[0]            \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2 + a3*b3) | a0*b6 + a1*b5 + 8*(a2*b9 + a4*b7 + a3*b8)
    "vmull.s32    q13, d8, d0[0]             \n\t"    // q13 = a0*b2 | a0*b7
    "vmlal.s32    q13, d6, d2[0]             \n\t"    // q13 = a0*b2 + a2*b0 | a0*b7 + a2*b5
    "vmlal.s32    q13, d7, d1[0]             \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 | a0*b7 + a2*b5 + a1*b6
    "vmlal.s32    q13, d10, d15[0]           \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4) | a0*b7 + a2*b5 + a1*b6 + 8*(a3*b9)
    "vmlal.s32    q13, d9, d5[0]             \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4 + a4*b3) | a0*b7 + a2*b5 + a1*b6 + 8*(a3*b9 + a4*b8)
    "vmull.s32    q14, d9, d0[0]             \n\t"    // q14 = a0*b3 | a0*b8
    "vmlal.s32    q14, d6, d3[0]             \n\t"    // q14 = a0*b3 + a3*b0 | a0*b8 + a3*b5
    "vmlal.s32    q14, d8, d1[0]             \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 | a0*b8 + a3*b5 + a1*b7
    "vmlal.s32    q14, d7, d2[0]             \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 + a2*b1 | a0*b8 + a3*b5 + a1*b7 + a2*b6
    "vmlal.s32    q14, d10, d5[0]            \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 + a2*b1 + 8*(a4*b4) | a0*b8 + a3*b5 + a1*b7 + a2*b6 + 8*(a4*b9)
    "vmull.s32    q15, d10, d0[0]            \n\t"    // q15 = a0*b4 | a0*b9
    "vmlal.s32    q15, d6, d4[0]             \n\t"    // q15 = a0*b4 + a4*b0 | a0*b9 + a4*b5
    "vmlal.s32    q15, d9, d1[0]             \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 | a0*b9 + a4*b5 + a1*b8
    "vmlal.s32    q15, d7, d3[0]             \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 + a3*b1 | a0*b9 + a4*b5 + a1*b8 + a3*b6
    "vmlal.s32    q15, d8, d2[0]             \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 + a3*b1 + a2*b2 | a0*b9 + a4*b5 + a1*b8 + a3*b6 + a2*b7
    
    // q11[0]->q13[0] = A0*B0 - A1*B1, q11[1]->q13[1] = A0*B0 + A1*B1 
    "vsub.s64     d22, d22, d17              \n\t"
    "vsub.s64     d24, d24, d19              \n\t"
    "vsub.s64     d26, d26, d21              \n\t"
    "vadd.s64     d23, d23, d16              \n\t"
    "vadd.s64     d25, d25, d18              \n\t"
    "vadd.s64     d27, d27, d20              \n\t"
    "vshr.s64     q10, q11, #26              \n\t"    ///
        
    // Complete computation A1*B0|A1*B1, q8-q9 <- A1*B0|A1*B1 
    "vmull.s32    q8, d9, d0[1]             \n\t"    // q8 = a5*b3 | a5*b8
    "vmlal.s32    q8, d6, d3[1]             \n\t"    // q8 = a5*b3 + a8*b0 | a5*b8 + a8*b5
    "vmlal.s32    q8, d8, d1[1]             \n\t"    // q8 = a5*b3 + a8*b0 + a6*b2 | a5*b8 + a8*b5 + a6*b7
    "vmlal.s32    q8, d7, d2[1]             \n\t"    // q8 = a5*b3 + a8*b0 + a6*b2 + a7*b1 | a5*b8 + a8*b5 + a6*b7 + a7*b6
    "vmlal.s32    q8, d10, d5[1]            \n\t"    // q8 = a5*b3 + a8*b0 + a6*b2 + a7*b1 + 8*(a9*b4) | a5*b8 + a8*b5 + a6*b7 + a7*b6 + 8*(a9*b9)
    "vmull.s32    q9, d10, d0[1]            \n\t"    // q9 = a5*b4 | a5*b9
    "vmlal.s32    q9, d6, d4[1]             \n\t"    // q9 = a5*b4 + a9*b0 | a5*b9 + a9*b5
    "vmlal.s32    q9, d9, d1[1]             \n\t"    // q9 = a5*b4 + a9*b0 + a6*b3 | a5*b9 + a9*b5 + a6*b8
    "vmlal.s32    q9, d7, d3[1]             \n\t"    // q9 = a5*b4 + a9*b0 + a6*b3 + a8*b1 | a5*b9 + a9*b5 + a6*b8 + a8*b6
    "vmlal.s32    q9, d8, d2[1]             \n\t"    // q9 = a5*b4 + a9*b0 + a6*b3 + a8*b1 + a7*b2 | a5*b9 + a9*b5 + a6*b8 + a8*b6 + a7*b7
    "vadd.s64     q10, q12, q10             \n\t"    ///
    
    // Complete q14[0]->q15[0] = A0*B0 - A1*B1, q14[1]->q15[0] = A0*B1 + A1*B0 
    "vmov.i64     q7, 0xFFFFFFFF            \n\t"    // mask_26
    "vadd.s64     d29, d29, d16             \n\t"
    "vsub.s64     d28, d28, d17             \n\t"
    "vshr.s64     q8, q10, #26              \n\t"    ///
    "vadd.s64     d31, d31, d18             \n\t"
    "vadd.s64     q8, q13, q8               \n\t"    ///
    "vsub.s64     d30, d30, d19             \n\t"
    
    // Reduction  
    "vshr.s64     q9, q8, #26               \n\t"
	"vshr.u64     q7, q7, #6                \n\t"
    "vadd.s64     q9, q14, q9               \n\t"
    "vand.u64     q1, q10, q7               \n\t"
    "vand.u64     q0, q11, q7               \n\t"    
    "vshr.s64     q10, q9, #26              \n\t"
    "vand.u64     q2, q8, q7                \n\t"
    "vadd.s64     q10, q15, q10             \n\t"
	"vshr.u64     q6, q7, #3                \n\t"    // mask_23
    "vand.u64     q3, q9, q7                \n\t"    
    "vshr.s64     q8, q10, #23              \n\t"
    "vand.u64     q11, q10, q6              \n\t"
    "vand.u64     q12, q8, q7               \n\t"
    "vshr.s64     q6, q8, #26               \n\t"
    "vadd.s64     q0, q0, q12               \n\t"
    "vadd.s64     q1, q1, q6                \n\t"
    "vst2.32      {d0[0],d1[0]}, [%2]!      \n\t"
    "vst2.32      {d2[0],d3[0]}, [%2]!      \n\t"
#if defined(INTERLEAVE)
    "vpop         {q4-q7}                   \n\t"
#endif
    "vst2.32      {d4[0],d5[0]}, [%2]!      \n\t"    
    "vst2.32      {d6[0],d7[0]}, [%2]!      \n\t"
    "vst2.32      {d22[0],d23[0]}, [%2]!    \n\t"
#if !defined(INTERLEAVE)
    "vpop         {q4-q7}                   \n\t"
#endif
    :
    :"r"(&a[0]), "r"(&b[0]), "r"(&c[0])
	);
	return;
}


void v2sqr1271_a(uint32_t* a, uint32_t* c) 
{ // Squaring over GF((2^127-1)^2) using NEON instructions
  // Operation: c = a^2 in GF((2^127-1)^2)
  // Representation: 23/23/26/26/26/26/26/26/26/26-bit

asm volatile(
    // q0-q2[0] <- A0|A1
    "vld1.8       {d0,d1}, [%0]!            \n\t"
    "vld1.8       {d2,d3}, [%0]!            \n\t"
#if defined(INTERLEAVE)
    "vshr.u64     q3, q0, #32               \n\t"    // q3-q5[0] = (0|a9->0|a5) 
    "vld1.8       {d4}, [%0]!               \n\t"
    "vpush        {q4-q7}                   \n\t" 
    "vmov.i64     q8, 0xFFFFFFFF            \n\t"    // mask_26
#else
    "vld1.8       {d4}, [%0]!               \n\t"
    "vpush        {q4-q7}                   \n\t" 
    "vmov.i64     q8, 0xFFFFFFFF            \n\t"    // mask_26
    "vshr.u64     q3, q0, #32               \n\t"    // q3-q5[0] = (0|a9->0|a5) 
#endif
    "vshr.u64     q4, q1, #32               \n\t" 
    "vshr.u64     d10, d4, #32              \n\t"    
	"vsub.s32     q13, q0, q3               \n\t"    // q13-q15[0] = (a9|a4-a9->a5|a0-a5)
	"vsub.s32     q14, q1, q4               \n\t" 
	"vsub.s32     d30, d4, d10              \n\t"     
	"vadd.s32     q3, q0, q3                \n\t"    // q3-q5[0] = (a9|a4+a9->a5|a0+a5)
	"vadd.s32     q4, q1, q4                \n\t" 
	"vadd.s32     d10, d4, d10              \n\t" 
    "vshl.i64     q0, q0, #33               \n\t"    // q0-q2[0] = (2*a4|0->2*a0|0) 
    "vshl.i64     q1, q1, #33               \n\t" 
    "vshl.i64     d4, d4, #33               \n\t"      
    "vbit         q0, q13, q8               \n\t"    // q0-q2[0] = (2*a4|a4-a9->2*a0|a0-a5)
    "vbit         q1, q14, q8               \n\t"  
    "vbit         d4, d30, d16              \n\t"     
    "vshl.i32     q6, q0, #3                \n\t"    // 8*(2*a4|a4-a9->2*a0|a0-a5)   
    "vshl.i32     q7, q1, #3                \n\t" 
    "vshl.i32     d5, d4, #3                \n\t"    
        
    // q11-q15 <- 2*A0*B0|(A0+A1)*(A0-A1)   
    "vmull.s32    q11, d6, d0               \n\t"    // q11 = a0*b0 | a0*b5
    "vmlal.s32    q11, d10, d13             \n\t"    // q11 = a0*b0 + 8*(a1*b4) | a0*b5 + 8*(a1*b9)
    "vmlal.s32    q11, d7, d5               \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1) | a0*b5 + 8*(a1*b9 + a4*b6)
    "vmlal.s32    q11, d9, d14              \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3) | a0*b5 + 8*(a1*b9 + a4*b6 + a2*b8)
    "vmlal.s32    q11, d8, d15              \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3 + a3*b2) | a0*b5 + 8*(a1*b9 + a4*b6 + a2*b8 + a3*b7)
    "vmull.s32    q12, d7, d0               \n\t"    // q12 = a0*b1 | a0*b6
    "vmlal.s32    q12, d6, d1               \n\t"    // q12 = a0*b1 + a1*b0 | a0*b6 + a1*b5
    "vmlal.s32    q12, d10, d14             \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4) | a0*b6 + a1*b5 + 8*(a2*b9)
    "vmlal.s32    q12, d8, d5               \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2) | a0*b6 + a1*b5 + 8*(a2*b9 + a4*b7)
    "vmlal.s32    q12, d9, d15              \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2 + a3*b3) | a0*b6 + a1*b5 + 8*(a2*b9 + a4*b7 + a3*b8)
    "vmull.s32    q13, d8, d0               \n\t"    // q13 = a0*b2 | a0*b7
    "vmlal.s32    q13, d6, d2               \n\t"    // q13 = a0*b2 + a2*b0 | a0*b7 + a2*b5
    "vmlal.s32    q13, d7, d1               \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 | a0*b7 + a2*b5 + a1*b6
    "vmlal.s32    q13, d10, d15             \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4) | a0*b7 + a2*b5 + a1*b6 + 8*(a3*b9)
    "vmlal.s32    q13, d9, d5               \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4 + a4*b3) | a0*b7 + a2*b5 + a1*b6 + 8*(a3*b9 + a4*b8)
    "vshr.s64     q10, q11, #26             \n\t"    ///
    "vmull.s32    q14, d9, d0               \n\t"    // q14 = a0*b3 | a0*b8
    "vadd.s64     q10, q12, q10             \n\t"    ///
    "vmlal.s32    q14, d6, d3               \n\t"    // q14 = a0*b3 + a3*b0 | a0*b8 + a3*b5
    "vmlal.s32    q14, d8, d1               \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 | a0*b8 + a3*b5 + a1*b7
    "vmlal.s32    q14, d7, d2               \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 + a2*b1 | a0*b8 + a3*b5 + a1*b7 + a2*b6
    "vmlal.s32    q14, d10, d5              \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 + a2*b1 + 8*(a4*b4) | a0*b8 + a3*b5 + a1*b7 + a2*b6 + 8*(a4*b9)
    "vmull.s32    q15, d10, d0              \n\t"    // q15 = a0*b4 | a0*b9
    "vshr.s64     q5, q10, #26              \n\t"    ///
    "vmlal.s32    q15, d6, d4               \n\t"    // q15 = a0*b4 + a4*b0 | a0*b9 + a4*b5
    "vadd.s64     q5, q13, q5               \n\t"    ///
    "vmlal.s32    q15, d9, d1               \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 | a0*b9 + a4*b5 + a1*b8
    "vmlal.s32    q15, d7, d3               \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 + a3*b1 | a0*b9 + a4*b5 + a1*b8 + a3*b6
    "vmlal.s32    q15, d8, d2               \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 + a3*b1 + a2*b2 | a0*b9 + a4*b5 + a1*b8 + a3*b6 + a2*b7
       
    // Reduction  
    "vshr.s64     q9, q5, #26               \n\t"
	"vshr.u64     q8, q8, #6                \n\t"
    "vadd.s64     q9, q14, q9               \n\t"
    "vand.u64     q1, q10, q8               \n\t"
    "vand.u64     q0, q11, q8               \n\t"    
    "vshr.s64     q10, q9, #26              \n\t"
    "vand.u64     q2, q5, q8                \n\t"
    "vadd.s64     q10, q15, q10             \n\t"
	"vshr.u64     q6, q8, #3                \n\t"    // mask_23    
    "vshr.s64     q5, q10, #23              \n\t"
    "vand.u64     q3, q9, q8                \n\t"
    "vand.u64     q11, q10, q6              \n\t"
    "vand.u64     q12, q5, q8               \n\t"
    "vshr.s64     q6, q5, #26               \n\t"
    "vadd.s64     q0, q0, q12               \n\t"
    "vadd.s64     q1, q1, q6                \n\t"
    "vst2.32      {d0[0],d1[0]}, [%1]!      \n\t"
    "vst2.32      {d2[0],d3[0]}, [%1]!      \n\t"
#if defined(INTERLEAVE)
    "vpop         {q4-q7}                   \n\t"
#endif
    "vst2.32      {d4[0],d5[0]}, [%1]!      \n\t"    
    "vst2.32      {d6[0],d7[0]}, [%1]!      \n\t"
    "vst2.32      {d22[0],d23[0]}, [%1]!    \n\t"
#if !defined(INTERLEAVE)
    "vpop         {q4-q7}                   \n\t"
#endif
    :
    :"r"(&a[0]), "r"(&c[0])
	);
	return;
}


#if defined(MIX_ARM_NEON)

void v2muladdsub1271_a(uint32_t* a, uint32_t* b, uint32_t* c, uint32_t* d, uint32_t* e, uint32_t* f, uint32_t* g)  
{ // Multiplication/addition over GF((2^127-1)^2) combining ARM and NEON instructions
  // Operation: c = a * b, f = d + e, g = d - e in GF((2^127-1)^2)
  // Representation: 23/23/26/26/26/26/26/26/26/26-bit

asm volatile(
    // q0-q2[0] <- A0|A1, q3-q5[0] <- B0|B1
    "vld1.8       {d0,d1}, [%0]!            \n\t"
    "vld1.8       {d2,d3}, [%0]!            \n\t"
    "vld1.8       {d4}, [%0]!               \n\t" 
    "vshl.i32     q6, q0, #3                \n\t"    // q6-q7 = 8*(a0->a8)  
    "vpush        {q4-q7}                   \n\t" 
    "vld1.8       {d6,d7}, [%1]!          	\n\t"
    "vshl.i32     d5, d4, #3                \n\t"    // q2[1] = 8*(a4|a9)   
    "vld1.8       {d8,d9}, [%1]!          	\n\t"
    "vshl.i32     q7, q1, #3                \n\t"  
    "vld1.8       {d10}, [%1]!              \n\t"

    // q11-q15 <- A0*B0|A0*B1, q8-q10 <- A1*B0|A1*B1 
    "vmull.s32    q8, d6, d0[1]             \n\t"    // q8 = a5*b0 | a5*b5
    "vmlal.s32    q8, d10, d13[1]           \n\t"    // q8 = a5*b0 + 8*(a6*b4) | a5*b5 + 8*(a6*b9)
    "vmlal.s32    q8, d7, d5[1]             \n\t"    // q8 = a5*b0 + 8*(a6*b4 + a9*b1) | a5*b5 + 8*(a6*b9 + a9*b6)
    "vmlal.s32    q8, d9, d14[1]            \n\t"    // q8 = a5*b0 + 8*(a6*b4 + a9*b1 + a7*b3) | a5*b5 + 8*(a6*b9 + a9*b6 + a7*b8)
    "vmlal.s32    q8, d8, d15[1]            \n\t"    // q8 = a5*b0 + 8*(a6*b4 + a9*b1 + a7*b3 + a8*b2) | a5*b5 + 8*(a6*b9 + a9*b6 + a7*b8 + a8*b7)
    "vmull.s32    q9, d7, d0[1]             \n\t"    // q9 = a5*b1 | a5*b6
    "vmlal.s32    q9, d6, d1[1]             \n\t"    // q9 = a5*b1 + a6*b0 | a5*b6 + a6*b5
    "vmlal.s32    q9, d10, d14[1]           \n\t"    // q9 = a5*b1 + a6*b0 + 8*(a7*b4) | a5*b6 + a6*b5 + 8*(a7*b9)
    "vmlal.s32    q9, d8, d5[1]             \n\t"    // q9 = a5*b1 + a6*b0 + 8*(a7*b4 + a9*b2) | a5*b6 + a6*b5 + 8*(a7*b9 + a9*b7)
    "vmlal.s32    q9, d9, d15[1]            \n\t"    // q9 = a5*b1 + a6*b0 + 8*(a7*b4 + a9*b2 + a8*b3) | a5*b6 + a6*b5 + 8*(a7*b9 + a9*b7 + a8*b8)
    "vmull.s32    q10, d8, d0[1]            \n\t"    // q10 = a5*b2 | a5*b7
    "vmlal.s32    q10, d6, d2[1]            \n\t"    // q10 = a5*b2 + a7*b0 | a5*b7 + a7*b5
    "vmlal.s32    q10, d7, d1[1]            \n\t"    // q10 = a5*b2 + a7*b0 + a6*b1 | a5*b7 + a7*b5 + a6*b6
    "vmlal.s32    q10, d10, d15[1]          \n\t"    // q10 = a5*b2 + a7*b0 + a6*b1 + 8*(a8*b4) | a5*b7 + a7*b5 + a6*b6 + 8*(a8*b9)
    "vmlal.s32    q10, d9, d5[1]            \n\t"    // q10 = a5*b2 + a7*b0 + a6*b1 + 8*(a8*b4 + a9*b3) | a5*b7 + a7*b5 + a6*b6 + 8*(a8*b9 + a9*b8)
         
    "push         {r7-r12}                 	\n\t"
    "ldm          %3!, {r7,r8}              \n\t"    //***** 
    "ldm          %4!, {r9,r10}             \n\t"    //***** 
    "add          r11, r7, r9               \n\t"    //***** 
    "add          r12, r8, r10              \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "sub          r8, r10                   \n\t"    //***** 
    "stm          %5!, {r11,r12}            \n\t"    //*****
    "stm          %6!, {r7,r8}              \n\t"    //*****
    
    "vmull.s32    q11, d6, d0[0]             \n\t"    // q11 = a0*b0 | a0*b5
    "vmlal.s32    q11, d10, d13[0]           \n\t"    // q11 = a0*b0 + 8*(a1*b4) | a0*b5 + 8*(a1*b9)
    "vmlal.s32    q11, d7, d5[0]             \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1) | a0*b5 + 8*(a1*b9 + a4*b6)
    "vmlal.s32    q11, d9, d14[0]            \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3) | a0*b5 + 8*(a1*b9 + a4*b6 + a2*b8)
    "vmlal.s32    q11, d8, d15[0]            \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3 + a3*b2) | a0*b5 + 8*(a1*b9 + a4*b6 + a2*b8 + a3*b7)
    "vmull.s32    q12, d7, d0[0]             \n\t"    // q12 = a0*b1 | a0*b6
    "vmlal.s32    q12, d6, d1[0]             \n\t"    // q12 = a0*b1 + a1*b0 | a0*b6 + a1*b5
    "vmlal.s32    q12, d10, d14[0]           \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4) | a0*b6 + a1*b5 + 8*(a2*b9)
    "vmlal.s32    q12, d8, d5[0]             \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2) | a0*b6 + a1*b5 + 8*(a2*b9 + a4*b7)
    "vmlal.s32    q12, d9, d15[0]            \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2 + a3*b3) | a0*b6 + a1*b5 + 8*(a2*b9 + a4*b7 + a3*b8)
    "vmull.s32    q13, d8, d0[0]             \n\t"    // q13 = a0*b2 | a0*b7
    "vmlal.s32    q13, d6, d2[0]             \n\t"    // q13 = a0*b2 + a2*b0 | a0*b7 + a2*b5
    "vmlal.s32    q13, d7, d1[0]             \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 | a0*b7 + a2*b5 + a1*b6
    "vmlal.s32    q13, d10, d15[0]           \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4) | a0*b7 + a2*b5 + a1*b6 + 8*(a3*b9)
    "vmlal.s32    q13, d9, d5[0]             \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4 + a4*b3) | a0*b7 + a2*b5 + a1*b6 + 8*(a3*b9 + a4*b8)
    "vmull.s32    q14, d9, d0[0]             \n\t"    // q14 = a0*b3 | a0*b8
    "vmlal.s32    q14, d6, d3[0]             \n\t"    // q14 = a0*b3 + a3*b0 | a0*b8 + a3*b5
    "vmlal.s32    q14, d8, d1[0]             \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 | a0*b8 + a3*b5 + a1*b7
    "vmlal.s32    q14, d7, d2[0]             \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 + a2*b1 | a0*b8 + a3*b5 + a1*b7 + a2*b6
    "vmlal.s32    q14, d10, d5[0]            \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 + a2*b1 + 8*(a4*b4) | a0*b8 + a3*b5 + a1*b7 + a2*b6 + 8*(a4*b9)
    "vmull.s32    q15, d10, d0[0]            \n\t"    // q15 = a0*b4 | a0*b9
    "vmlal.s32    q15, d6, d4[0]             \n\t"    // q15 = a0*b4 + a4*b0 | a0*b9 + a4*b5
    "vmlal.s32    q15, d9, d1[0]             \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 | a0*b9 + a4*b5 + a1*b8
    "vmlal.s32    q15, d7, d3[0]             \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 + a3*b1 | a0*b9 + a4*b5 + a1*b8 + a3*b6
    "vmlal.s32    q15, d8, d2[0]             \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 + a3*b1 + a2*b2 | a0*b9 + a4*b5 + a1*b8 + a3*b6 + a2*b7
         
    "ldm          %3!, {r7,r8}              \n\t"    //***** 
    "ldm          %4!, {r9,r10}             \n\t"    //***** 
    "add          r11, r7, r9               \n\t"    //***** 
    "add          r12, r8, r10              \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "sub          r8, r10                   \n\t"    //***** 
    "stm          %5!, {r11,r12}            \n\t"    //*****
    "stm          %6!, {r7,r8}              \n\t"    //*****
    
    // q11[0]->q13[0] = A0*B0 - A1*B1, q11[1]->q13[1] = A0*B0 + A1*B1 
    "vsub.s64     d22, d22, d17              \n\t"
    "vsub.s64     d24, d24, d19              \n\t"
    "vsub.s64     d26, d26, d21              \n\t"
    "vadd.s64     d23, d23, d16              \n\t"
    "vadd.s64     d25, d25, d18              \n\t"
    "vadd.s64     d27, d27, d20              \n\t"
    
    "ldm          %3!, {r7,r8}              \n\t"    //***** 
    "ldm          %4!, {r9,r10}             \n\t"    //***** 
    "add          r11, r7, r9               \n\t"    //***** 
    "add          r12, r8, r10              \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "sub          r8, r10                   \n\t"    //***** 
    "stm          %5!, {r11,r12}            \n\t"    //*****
    "stm          %6!, {r7,r8}              \n\t"    //*****
        
    // Complete computation A1*B0|A1*B1, q8-q9 <- A1*B0|A1*B1 
    "vmull.s32    q8, d9, d0[1]             \n\t"    // q8 = a5*b3 | a5*b8
    "vmlal.s32    q8, d6, d3[1]             \n\t"    // q8 = a5*b3 + a8*b0 | a5*b8 + a8*b5
    "vmlal.s32    q8, d8, d1[1]             \n\t"    // q8 = a5*b3 + a8*b0 + a6*b2 | a5*b8 + a8*b5 + a6*b7
    "vmlal.s32    q8, d7, d2[1]             \n\t"    // q8 = a5*b3 + a8*b0 + a6*b2 + a7*b1 | a5*b8 + a8*b5 + a6*b7 + a7*b6
    "vmlal.s32    q8, d10, d5[1]            \n\t"    // q8 = a5*b3 + a8*b0 + a6*b2 + a7*b1 + 8*(a9*b4) | a5*b8 + a8*b5 + a6*b7 + a7*b6 + 8*(a9*b9)
    "vmull.s32    q9, d10, d0[1]            \n\t"    // q9 = a5*b4 | a5*b9
    "vmlal.s32    q9, d6, d4[1]             \n\t"    // q9 = a5*b4 + a9*b0 | a5*b9 + a9*b5
    "vmlal.s32    q9, d9, d1[1]             \n\t"    // q9 = a5*b4 + a9*b0 + a6*b3 | a5*b9 + a9*b5 + a6*b8
    "vmlal.s32    q9, d7, d3[1]             \n\t"    // q9 = a5*b4 + a9*b0 + a6*b3 + a8*b1 | a5*b9 + a9*b5 + a6*b8 + a8*b6
    "vmlal.s32    q9, d8, d2[1]             \n\t"    // q9 = a5*b4 + a9*b0 + a6*b3 + a8*b1 + a7*b2 | a5*b9 + a9*b5 + a6*b8 + a8*b6 + a7*b7
    
    "ldm          %3!, {r7,r8}              \n\t"    //***** 
    "ldm          %4!, {r9,r10}             \n\t"    //***** 
    "add          r11, r7, r9               \n\t"    //***** 
    "add          r12, r8, r10              \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "sub          r8, r10                   \n\t"    //***** 
    "stm          %5!, {r11,r12}            \n\t"    //*****
    "stm          %6!, {r7,r8}              \n\t"    //*****
    
    // Complete q14[0]->q15[0] = A0*B0 - A1*B1, q14[1]->q15[0] = A0*B1 + A1*B0 
    "vmov.i64     q7, 0xFFFFFFFF            \n\t"    // mask_26
    "vadd.s64     d29, d29, d16             \n\t"
    "vsub.s64     d28, d28, d17             \n\t"
	"vshr.u64     q7, q7, #6                \n\t"
    "vadd.s64     d31, d31, d18             \n\t"
    "vsub.s64     d30, d30, d19             \n\t"
    
    "ldm          %3!, {r7,r8}              \n\t"    //***** 
    "ldm          %4!, {r9,r10}             \n\t"    //***** 
    "add          r11, r7, r9               \n\t"    //***** 
    "add          r12, r8, r10              \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "sub          r8, r10                   \n\t"    //***** 
    "stm          %5!, {r11,r12}            \n\t"    //*****
    "stm          %6!, {r7,r8}              \n\t"    //*****
    "pop          {r7-r12}                 	\n\t"
    
    // Reduction   
	"vshr.u64     q6, q7, #3                \n\t"    // mask_23
    "vshr.s64     q10, q11, #26             \n\t"
    "vand.u64     q0, q11, q7               \n\t"

    "vshr.s64     q9, q14, #26              \n\t"
    "vand.u64     q3, q14, q7               \n\t"
    "vadd.s64     q10, q12, q10             \n\t"
    "vadd.s64     q9, q15, q9               \n\t"
    "vand.u64     q1, q10, q7               \n\t"
    "vand.u64     q11, q9, q6               \n\t" 

    "vshr.s64     q10, q10, #26             \n\t"
    "vshr.s64     q9, q9, #23               \n\t"
    "vadd.s64     q10, q13, q10             \n\t"
    "vand.u64     q12, q9, q7               \n\t"
    "vand.u64     q2, q10, q7               \n\t"

    "vshr.s64     q10, q10, #26             \n\t"
    "vshr.s64     q6, q9, #26               \n\t"
    "vadd.s64     q10, q3, q10              \n\t"
    "vadd.s64     q0, q0, q12               \n\t"
    "vadd.s64     q1, q1, q6                \n\t"
    
    "vst2.32      {d0[0],d1[0]}, [%2]!      \n\t"
    "vshr.s64     q9, q10, #26              \n\t"
    "vst2.32      {d2[0],d3[0]}, [%2]!      \n\t"
    "vand.u64     q3, q10, q7               \n\t"
    "vst2.32      {d4[0],d5[0]}, [%2]!      \n\t" 
    "vpop         {q4-q7}                   \n\t"
    "vadd.s64     q11, q11, q9              \n\t"   
    "vst2.32      {d6[0],d7[0]}, [%2]!      \n\t"
    "vst2.32      {d22[0],d23[0]}, [%2]!    \n\t" 
    :
    :"r"(&a[0]), "r"(&b[0]), "r"(&c[0]), "r"(&d[0]), "r"(&e[0]), "r"(&f[0]), "r"(&g[0])
    :"memory","r7","r8","r9","r10","r11","r12"
	);
	return;
}


void v2muladd1271_a(uint32_t* a, uint32_t* b, uint32_t* c, uint32_t* d, uint32_t* e, uint32_t* f)  
{ // Multiplication/addition over GF((2^127-1)^2) combining ARM and NEON instructions
  // Operation: c = a * b, f = d + e in GF((2^127-1)^2)
  // Representation: 23/23/26/26/26/26/26/26/26/26-bit

asm volatile(
    // q0-q2[0] <- A0|A1, q3-q5[0] <- B0|B1
    "vld1.8       {d0,d1}, [%0]!            \n\t"
    "vld1.8       {d2,d3}, [%0]!            \n\t"
    "vld1.8       {d4}, [%0]!               \n\t" 
    "vshl.i32     q6, q0, #3                \n\t"    // q6-q7 = 8*(a0->a8)  
    "vpush        {q4-q7}                   \n\t" 
    "vld1.8       {d6,d7}, [%1]!          	\n\t"
    "vshl.i32     d5, d4, #3                \n\t"    // q2[1] = 8*(a4|a9)   
    "vld1.8       {d8,d9}, [%1]!          	\n\t"
    "vshl.i32     q7, q1, #3                \n\t"  
    "vld1.8       {d10}, [%1]!              \n\t"

    // q11-q15 <- A0*B0|A0*B1, q8-q10 <- A1*B0|A1*B1 
    "vmull.s32    q8, d6, d0[1]             \n\t"    // q8 = a5*b0 | a5*b5
    "vmlal.s32    q8, d10, d13[1]           \n\t"    // q8 = a5*b0 + 8*(a6*b4) | a5*b5 + 8*(a6*b9)
    "vmlal.s32    q8, d7, d5[1]             \n\t"    // q8 = a5*b0 + 8*(a6*b4 + a9*b1) | a5*b5 + 8*(a6*b9 + a9*b6)
    "vmlal.s32    q8, d9, d14[1]            \n\t"    // q8 = a5*b0 + 8*(a6*b4 + a9*b1 + a7*b3) | a5*b5 + 8*(a6*b9 + a9*b6 + a7*b8)
    "vmlal.s32    q8, d8, d15[1]            \n\t"    // q8 = a5*b0 + 8*(a6*b4 + a9*b1 + a7*b3 + a8*b2) | a5*b5 + 8*(a6*b9 + a9*b6 + a7*b8 + a8*b7)
    "vmull.s32    q9, d7, d0[1]             \n\t"    // q9 = a5*b1 | a5*b6
    "vmlal.s32    q9, d6, d1[1]             \n\t"    // q9 = a5*b1 + a6*b0 | a5*b6 + a6*b5
    "vmlal.s32    q9, d10, d14[1]           \n\t"    // q9 = a5*b1 + a6*b0 + 8*(a7*b4) | a5*b6 + a6*b5 + 8*(a7*b9)
    "vmlal.s32    q9, d8, d5[1]             \n\t"    // q9 = a5*b1 + a6*b0 + 8*(a7*b4 + a9*b2) | a5*b6 + a6*b5 + 8*(a7*b9 + a9*b7)
    "vmlal.s32    q9, d9, d15[1]            \n\t"    // q9 = a5*b1 + a6*b0 + 8*(a7*b4 + a9*b2 + a8*b3) | a5*b6 + a6*b5 + 8*(a7*b9 + a9*b7 + a8*b8)
    "vmull.s32    q10, d8, d0[1]            \n\t"    // q10 = a5*b2 | a5*b7
    "vmlal.s32    q10, d6, d2[1]            \n\t"    // q10 = a5*b2 + a7*b0 | a5*b7 + a7*b5
    "vmlal.s32    q10, d7, d1[1]            \n\t"    // q10 = a5*b2 + a7*b0 + a6*b1 | a5*b7 + a7*b5 + a6*b6
    "vmlal.s32    q10, d10, d15[1]          \n\t"    // q10 = a5*b2 + a7*b0 + a6*b1 + 8*(a8*b4) | a5*b7 + a7*b5 + a6*b6 + 8*(a8*b9)
    "vmlal.s32    q10, d9, d5[1]            \n\t"    // q10 = a5*b2 + a7*b0 + a6*b1 + 8*(a8*b4 + a9*b3) | a5*b7 + a7*b5 + a6*b6 + 8*(a8*b9 + a9*b8)
         
    "push         {r6-r9}                 	\n\t"
    "ldm          %3!, {r6,r7}              \n\t"    //***** 
    "ldm          %4!, {r8,r9}              \n\t"    //***** 
    "add          r6, r8                    \n\t"    //***** 
    "add          r7, r9                    \n\t"    //***** 
    "stm          %5!, {r6,r7}              \n\t"    //*****
    
    "vmull.s32    q11, d6, d0[0]             \n\t"    // q11 = a0*b0 | a0*b5
    "vmlal.s32    q11, d10, d13[0]           \n\t"    // q11 = a0*b0 + 8*(a1*b4) | a0*b5 + 8*(a1*b9)
    "vmlal.s32    q11, d7, d5[0]             \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1) | a0*b5 + 8*(a1*b9 + a4*b6)
    "vmlal.s32    q11, d9, d14[0]            \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3) | a0*b5 + 8*(a1*b9 + a4*b6 + a2*b8)
    "vmlal.s32    q11, d8, d15[0]            \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3 + a3*b2) | a0*b5 + 8*(a1*b9 + a4*b6 + a2*b8 + a3*b7)
    "vmull.s32    q12, d7, d0[0]             \n\t"    // q12 = a0*b1 | a0*b6
    "vmlal.s32    q12, d6, d1[0]             \n\t"    // q12 = a0*b1 + a1*b0 | a0*b6 + a1*b5
    "vmlal.s32    q12, d10, d14[0]           \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4) | a0*b6 + a1*b5 + 8*(a2*b9)
    "vmlal.s32    q12, d8, d5[0]             \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2) | a0*b6 + a1*b5 + 8*(a2*b9 + a4*b7)
    "vmlal.s32    q12, d9, d15[0]            \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2 + a3*b3) | a0*b6 + a1*b5 + 8*(a2*b9 + a4*b7 + a3*b8)
    "vmull.s32    q13, d8, d0[0]             \n\t"    // q13 = a0*b2 | a0*b7
    "vmlal.s32    q13, d6, d2[0]             \n\t"    // q13 = a0*b2 + a2*b0 | a0*b7 + a2*b5
    "vmlal.s32    q13, d7, d1[0]             \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 | a0*b7 + a2*b5 + a1*b6
    "vmlal.s32    q13, d10, d15[0]           \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4) | a0*b7 + a2*b5 + a1*b6 + 8*(a3*b9)
    "vmlal.s32    q13, d9, d5[0]             \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4 + a4*b3) | a0*b7 + a2*b5 + a1*b6 + 8*(a3*b9 + a4*b8)
    "vmull.s32    q14, d9, d0[0]             \n\t"    // q14 = a0*b3 | a0*b8
    "vmlal.s32    q14, d6, d3[0]             \n\t"    // q14 = a0*b3 + a3*b0 | a0*b8 + a3*b5
    "vmlal.s32    q14, d8, d1[0]             \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 | a0*b8 + a3*b5 + a1*b7
    "vmlal.s32    q14, d7, d2[0]             \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 + a2*b1 | a0*b8 + a3*b5 + a1*b7 + a2*b6
    "vmlal.s32    q14, d10, d5[0]            \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 + a2*b1 + 8*(a4*b4) | a0*b8 + a3*b5 + a1*b7 + a2*b6 + 8*(a4*b9)
    "vmull.s32    q15, d10, d0[0]            \n\t"    // q15 = a0*b4 | a0*b9
    "vmlal.s32    q15, d6, d4[0]             \n\t"    // q15 = a0*b4 + a4*b0 | a0*b9 + a4*b5
    "vmlal.s32    q15, d9, d1[0]             \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 | a0*b9 + a4*b5 + a1*b8
    "vmlal.s32    q15, d7, d3[0]             \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 + a3*b1 | a0*b9 + a4*b5 + a1*b8 + a3*b6
    "vmlal.s32    q15, d8, d2[0]             \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 + a3*b1 + a2*b2 | a0*b9 + a4*b5 + a1*b8 + a3*b6 + a2*b7
         
    "ldm          %3!, {r6,r7}              \n\t"    //***** 
    "ldm          %4!, {r8,r9}              \n\t"    //***** 
    "add          r6, r8                    \n\t"    //***** 
    "add          r7, r9                    \n\t"    //***** 
    "stm          %5!, {r6,r7}              \n\t"    //*****
    
    // q11[0]->q13[0] = A0*B0 - A1*B1, q11[1]->q13[1] = A0*B0 + A1*B1 
    "vsub.s64     d22, d22, d17              \n\t"
    "vsub.s64     d24, d24, d19              \n\t"
    "vsub.s64     d26, d26, d21              \n\t"
    "vadd.s64     d23, d23, d16              \n\t"
    "vadd.s64     d25, d25, d18              \n\t"
    "vadd.s64     d27, d27, d20              \n\t"
    
    "ldm          %3!, {r6,r7}              \n\t"    //***** 
    "ldm          %4!, {r8,r9}              \n\t"    //***** 
    "add          r6, r8                    \n\t"    //***** 
    "add          r7, r9                    \n\t"    //***** 
    "stm          %5!, {r6,r7}              \n\t"    //*****
        
    // Complete computation A1*B0|A1*B1, q8-q9 <- A1*B0|A1*B1 
    "vmull.s32    q8, d9, d0[1]             \n\t"    // q8 = a5*b3 | a5*b8
    "vmlal.s32    q8, d6, d3[1]             \n\t"    // q8 = a5*b3 + a8*b0 | a5*b8 + a8*b5
    "vmlal.s32    q8, d8, d1[1]             \n\t"    // q8 = a5*b3 + a8*b0 + a6*b2 | a5*b8 + a8*b5 + a6*b7
    "vmlal.s32    q8, d7, d2[1]             \n\t"    // q8 = a5*b3 + a8*b0 + a6*b2 + a7*b1 | a5*b8 + a8*b5 + a6*b7 + a7*b6
    "vmlal.s32    q8, d10, d5[1]            \n\t"    // q8 = a5*b3 + a8*b0 + a6*b2 + a7*b1 + 8*(a9*b4) | a5*b8 + a8*b5 + a6*b7 + a7*b6 + 8*(a9*b9)
    "vmull.s32    q9, d10, d0[1]            \n\t"    // q9 = a5*b4 | a5*b9
    "vmlal.s32    q9, d6, d4[1]             \n\t"    // q9 = a5*b4 + a9*b0 | a5*b9 + a9*b5
    "vmlal.s32    q9, d9, d1[1]             \n\t"    // q9 = a5*b4 + a9*b0 + a6*b3 | a5*b9 + a9*b5 + a6*b8
    "vmlal.s32    q9, d7, d3[1]             \n\t"    // q9 = a5*b4 + a9*b0 + a6*b3 + a8*b1 | a5*b9 + a9*b5 + a6*b8 + a8*b6
    "vmlal.s32    q9, d8, d2[1]             \n\t"    // q9 = a5*b4 + a9*b0 + a6*b3 + a8*b1 + a7*b2 | a5*b9 + a9*b5 + a6*b8 + a8*b6 + a7*b7
    
    "ldm          %3!, {r6,r7}              \n\t"    //***** 
    "ldm          %4!, {r8,r9}              \n\t"    //***** 
    "add          r6, r8                    \n\t"    //***** 
    "add          r7, r9                    \n\t"    //***** 
    "stm          %5!, {r6,r7}              \n\t"    //*****
    
    // Complete q14[0]->q15[0] = A0*B0 - A1*B1, q14[1]->q15[0] = A0*B1 + A1*B0 
    "vmov.i64     q7, 0xFFFFFFFF            \n\t"    // mask_26
    "vadd.s64     d29, d29, d16             \n\t"
    "vsub.s64     d28, d28, d17             \n\t"
	"vshr.u64     q7, q7, #6                \n\t"
    "vadd.s64     d31, d31, d18             \n\t"
    "vsub.s64     d30, d30, d19             \n\t"
    
    "ldm          %3!, {r6,r7}              \n\t"    //***** 
    "ldm          %4!, {r8,r9}              \n\t"    //***** 
    "add          r6, r8                    \n\t"    //***** 
    "add          r7, r9                    \n\t"    //***** 
    "stm          %5!, {r6,r7}              \n\t"    //*****
    "pop          {r6-r9}                 	\n\t"

    // Reduction   
	"vshr.u64     q6, q7, #3                \n\t"    // mask_23
    "vshr.s64     q10, q11, #26             \n\t"
    "vand.u64     q0, q11, q7               \n\t"

    "vshr.s64     q9, q14, #26              \n\t"
    "vand.u64     q3, q14, q7               \n\t"
    "vadd.s64     q10, q12, q10             \n\t"
    "vadd.s64     q9, q15, q9               \n\t"
    "vand.u64     q1, q10, q7               \n\t"
    "vand.u64     q11, q9, q6               \n\t" 

    "vshr.s64     q10, q10, #26             \n\t"
    "vshr.s64     q9, q9, #23               \n\t"
    "vadd.s64     q10, q13, q10             \n\t"
    "vand.u64     q12, q9, q7               \n\t"
    "vand.u64     q2, q10, q7               \n\t"

    "vshr.s64     q10, q10, #26             \n\t"
    "vshr.s64     q6, q9, #26               \n\t"
    "vadd.s64     q10, q3, q10              \n\t"
    "vadd.s64     q0, q0, q12               \n\t"
    "vadd.s64     q1, q1, q6                \n\t"
    
    "vst2.32      {d0[0],d1[0]}, [%2]!      \n\t"
    "vshr.s64     q9, q10, #26              \n\t"
    "vst2.32      {d2[0],d3[0]}, [%2]!      \n\t"
    "vand.u64     q3, q10, q7               \n\t"
    "vst2.32      {d4[0],d5[0]}, [%2]!      \n\t" 
    "vpop         {q4-q7}                   \n\t"
    "vadd.s64     q11, q11, q9              \n\t"   
    "vst2.32      {d6[0],d7[0]}, [%2]!      \n\t"
    "vst2.32      {d22[0],d23[0]}, [%2]!    \n\t"
    :
    :"r"(&a[0]), "r"(&b[0]), "r"(&c[0]), "r"(&d[0]), "r"(&e[0]), "r"(&f[0])
    :"memory","r6","r7","r8","r9"
	);
	return;
}


void v2mulsub1271_a(uint32_t* a, uint32_t* b, uint32_t* c, uint32_t* d, uint32_t* e, uint32_t* f)  
{ // Multiplication/subtraction over GF((2^127-1)^2) combining ARM and NEON instructions 
  // Operation: c = a * b, f = d - e in GF((2^127-1)^2)
  // Representation: 23/23/26/26/26/26/26/26/26/26-bit

asm volatile(
    // q0-q2[0] <- A0|A1, q3-q5[0] <- B0|B1
    "vld1.8       {d0,d1}, [%0]!            \n\t"
    "vld1.8       {d2,d3}, [%0]!            \n\t"
    "vld1.8       {d4}, [%0]!               \n\t" 
    "vshl.i32     q6, q0, #3                \n\t"    // q6-q7 = 8*(a0->a8)  
    "vpush        {q4-q7}                   \n\t" 
    "vld1.8       {d6,d7}, [%1]!          	\n\t"
    "vshl.i32     d5, d4, #3                \n\t"    // q2[1] = 8*(a4|a9)   
    "vld1.8       {d8,d9}, [%1]!          	\n\t"
    "vshl.i32     q7, q1, #3                \n\t"  
    "vld1.8       {d10}, [%1]!              \n\t"

    // q11-q15 <- A0*B0|A0*B1, q8-q10 <- A1*B0|A1*B1 
    "vmull.s32    q8, d6, d0[1]             \n\t"    // q8 = a5*b0 | a5*b5
    "vmlal.s32    q8, d10, d13[1]           \n\t"    // q8 = a5*b0 + 8*(a6*b4) | a5*b5 + 8*(a6*b9)
    "vmlal.s32    q8, d7, d5[1]             \n\t"    // q8 = a5*b0 + 8*(a6*b4 + a9*b1) | a5*b5 + 8*(a6*b9 + a9*b6)
    "vmlal.s32    q8, d9, d14[1]            \n\t"    // q8 = a5*b0 + 8*(a6*b4 + a9*b1 + a7*b3) | a5*b5 + 8*(a6*b9 + a9*b6 + a7*b8)
    "vmlal.s32    q8, d8, d15[1]            \n\t"    // q8 = a5*b0 + 8*(a6*b4 + a9*b1 + a7*b3 + a8*b2) | a5*b5 + 8*(a6*b9 + a9*b6 + a7*b8 + a8*b7)
    "vmull.s32    q9, d7, d0[1]             \n\t"    // q9 = a5*b1 | a5*b6
    "vmlal.s32    q9, d6, d1[1]             \n\t"    // q9 = a5*b1 + a6*b0 | a5*b6 + a6*b5
    "vmlal.s32    q9, d10, d14[1]           \n\t"    // q9 = a5*b1 + a6*b0 + 8*(a7*b4) | a5*b6 + a6*b5 + 8*(a7*b9)
    "vmlal.s32    q9, d8, d5[1]             \n\t"    // q9 = a5*b1 + a6*b0 + 8*(a7*b4 + a9*b2) | a5*b6 + a6*b5 + 8*(a7*b9 + a9*b7)
    "vmlal.s32    q9, d9, d15[1]            \n\t"    // q9 = a5*b1 + a6*b0 + 8*(a7*b4 + a9*b2 + a8*b3) | a5*b6 + a6*b5 + 8*(a7*b9 + a9*b7 + a8*b8)
    "vmull.s32    q10, d8, d0[1]            \n\t"    // q10 = a5*b2 | a5*b7
    "vmlal.s32    q10, d6, d2[1]            \n\t"    // q10 = a5*b2 + a7*b0 | a5*b7 + a7*b5
    "vmlal.s32    q10, d7, d1[1]            \n\t"    // q10 = a5*b2 + a7*b0 + a6*b1 | a5*b7 + a7*b5 + a6*b6
    "vmlal.s32    q10, d10, d15[1]          \n\t"    // q10 = a5*b2 + a7*b0 + a6*b1 + 8*(a8*b4) | a5*b7 + a7*b5 + a6*b6 + 8*(a8*b9)
    "vmlal.s32    q10, d9, d5[1]            \n\t"    // q10 = a5*b2 + a7*b0 + a6*b1 + 8*(a8*b4 + a9*b3) | a5*b7 + a7*b5 + a6*b6 + 8*(a8*b9 + a9*b8)
         
    "push         {r6-r9}                 	\n\t"
    "ldm          %3!, {r6,r7}              \n\t"    //***** 
    "ldm          %4!, {r8,r9}              \n\t"    //***** 
    "sub          r6, r8                    \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "stm          %5!, {r6,r7}              \n\t"    //*****
    
    "vmull.s32    q11, d6, d0[0]             \n\t"    // q11 = a0*b0 | a0*b5
    "vmlal.s32    q11, d10, d13[0]           \n\t"    // q11 = a0*b0 + 8*(a1*b4) | a0*b5 + 8*(a1*b9)
    "vmlal.s32    q11, d7, d5[0]             \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1) | a0*b5 + 8*(a1*b9 + a4*b6)
    "vmlal.s32    q11, d9, d14[0]            \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3) | a0*b5 + 8*(a1*b9 + a4*b6 + a2*b8)
    "vmlal.s32    q11, d8, d15[0]            \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3 + a3*b2) | a0*b5 + 8*(a1*b9 + a4*b6 + a2*b8 + a3*b7)
    "vmull.s32    q12, d7, d0[0]             \n\t"    // q12 = a0*b1 | a0*b6
    "vmlal.s32    q12, d6, d1[0]             \n\t"    // q12 = a0*b1 + a1*b0 | a0*b6 + a1*b5
    "vmlal.s32    q12, d10, d14[0]           \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4) | a0*b6 + a1*b5 + 8*(a2*b9)
    "vmlal.s32    q12, d8, d5[0]             \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2) | a0*b6 + a1*b5 + 8*(a2*b9 + a4*b7)
    "vmlal.s32    q12, d9, d15[0]            \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2 + a3*b3) | a0*b6 + a1*b5 + 8*(a2*b9 + a4*b7 + a3*b8)
    "vmull.s32    q13, d8, d0[0]             \n\t"    // q13 = a0*b2 | a0*b7
    "vmlal.s32    q13, d6, d2[0]             \n\t"    // q13 = a0*b2 + a2*b0 | a0*b7 + a2*b5
    "vmlal.s32    q13, d7, d1[0]             \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 | a0*b7 + a2*b5 + a1*b6
    "vmlal.s32    q13, d10, d15[0]           \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4) | a0*b7 + a2*b5 + a1*b6 + 8*(a3*b9)
    "vmlal.s32    q13, d9, d5[0]             \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4 + a4*b3) | a0*b7 + a2*b5 + a1*b6 + 8*(a3*b9 + a4*b8)
    "vmull.s32    q14, d9, d0[0]             \n\t"    // q14 = a0*b3 | a0*b8
    "vmlal.s32    q14, d6, d3[0]             \n\t"    // q14 = a0*b3 + a3*b0 | a0*b8 + a3*b5
    "vmlal.s32    q14, d8, d1[0]             \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 | a0*b8 + a3*b5 + a1*b7
    "vmlal.s32    q14, d7, d2[0]             \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 + a2*b1 | a0*b8 + a3*b5 + a1*b7 + a2*b6
    "vmlal.s32    q14, d10, d5[0]            \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 + a2*b1 + 8*(a4*b4) | a0*b8 + a3*b5 + a1*b7 + a2*b6 + 8*(a4*b9)
    "vmull.s32    q15, d10, d0[0]            \n\t"    // q15 = a0*b4 | a0*b9
    "vmlal.s32    q15, d6, d4[0]             \n\t"    // q15 = a0*b4 + a4*b0 | a0*b9 + a4*b5
    "vmlal.s32    q15, d9, d1[0]             \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 | a0*b9 + a4*b5 + a1*b8
    "vmlal.s32    q15, d7, d3[0]             \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 + a3*b1 | a0*b9 + a4*b5 + a1*b8 + a3*b6
    "vmlal.s32    q15, d8, d2[0]             \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 + a3*b1 + a2*b2 | a0*b9 + a4*b5 + a1*b8 + a3*b6 + a2*b7
         
    "ldm          %3!, {r6,r7}              \n\t"    //***** 
    "ldm          %4!, {r8,r9}              \n\t"    //***** 
    "sub          r6, r8                    \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "stm          %5!, {r6,r7}              \n\t"    //*****
    
    // q11[0]->q13[0] = A0*B0 - A1*B1, q11[1]->q13[1] = A0*B0 + A1*B1 
    "vsub.s64     d22, d22, d17              \n\t"
    "vsub.s64     d24, d24, d19              \n\t"
    "vsub.s64     d26, d26, d21              \n\t"
    "vadd.s64     d23, d23, d16              \n\t"
    "vadd.s64     d25, d25, d18              \n\t"
    "vadd.s64     d27, d27, d20              \n\t"
    
    "ldm          %3!, {r6,r7}              \n\t"    //***** 
    "ldm          %4!, {r8,r9}              \n\t"    //***** 
    "sub          r6, r8                    \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "stm          %5!, {r6,r7}              \n\t"    //*****
        
    // Complete computation A1*B0|A1*B1, q8-q9 <- A1*B0|A1*B1 
    "vmull.s32    q8, d9, d0[1]             \n\t"    // q8 = a5*b3 | a5*b8
    "vmlal.s32    q8, d6, d3[1]             \n\t"    // q8 = a5*b3 + a8*b0 | a5*b8 + a8*b5
    "vmlal.s32    q8, d8, d1[1]             \n\t"    // q8 = a5*b3 + a8*b0 + a6*b2 | a5*b8 + a8*b5 + a6*b7
    "vmlal.s32    q8, d7, d2[1]             \n\t"    // q8 = a5*b3 + a8*b0 + a6*b2 + a7*b1 | a5*b8 + a8*b5 + a6*b7 + a7*b6
    "vmlal.s32    q8, d10, d5[1]            \n\t"    // q8 = a5*b3 + a8*b0 + a6*b2 + a7*b1 + 8*(a9*b4) | a5*b8 + a8*b5 + a6*b7 + a7*b6 + 8*(a9*b9)
    "vmull.s32    q9, d10, d0[1]            \n\t"    // q9 = a5*b4 | a5*b9
    "vmlal.s32    q9, d6, d4[1]             \n\t"    // q9 = a5*b4 + a9*b0 | a5*b9 + a9*b5
    "vmlal.s32    q9, d9, d1[1]             \n\t"    // q9 = a5*b4 + a9*b0 + a6*b3 | a5*b9 + a9*b5 + a6*b8
    "vmlal.s32    q9, d7, d3[1]             \n\t"    // q9 = a5*b4 + a9*b0 + a6*b3 + a8*b1 | a5*b9 + a9*b5 + a6*b8 + a8*b6
    "vmlal.s32    q9, d8, d2[1]             \n\t"    // q9 = a5*b4 + a9*b0 + a6*b3 + a8*b1 + a7*b2 | a5*b9 + a9*b5 + a6*b8 + a8*b6 + a7*b7
    
    "ldm          %3!, {r6,r7}              \n\t"    //***** 
    "ldm          %4!, {r8,r9}              \n\t"    //***** 
    "sub          r6, r8                    \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "stm          %5!, {r6,r7}              \n\t"    //*****
    
    // Complete q14[0]->q15[0] = A0*B0 - A1*B1, q14[1]->q15[0] = A0*B1 + A1*B0 
    "vmov.i64     q7, 0xFFFFFFFF            \n\t"    // mask_26
    "vadd.s64     d29, d29, d16             \n\t"
    "vsub.s64     d28, d28, d17             \n\t"
	"vshr.u64     q7, q7, #6                \n\t"
    "vadd.s64     d31, d31, d18             \n\t"
    "vsub.s64     d30, d30, d19             \n\t"
    
    "ldm          %3!, {r6,r7}              \n\t"    //***** 
    "ldm          %4!, {r8,r9}              \n\t"    //***** 
    "sub          r6, r8                    \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "stm          %5!, {r6,r7}              \n\t"    //*****
    "pop          {r6-r9}                 	\n\t"

    // Reduction   
	"vshr.u64     q6, q7, #3                \n\t"    // mask_23
    "vshr.s64     q10, q11, #26             \n\t"
    "vand.u64     q0, q11, q7               \n\t"

    "vshr.s64     q9, q14, #26              \n\t"
    "vand.u64     q3, q14, q7               \n\t"
    "vadd.s64     q10, q12, q10             \n\t"
    "vadd.s64     q9, q15, q9               \n\t"
    "vand.u64     q1, q10, q7               \n\t"
    "vand.u64     q11, q9, q6               \n\t" 

    "vshr.s64     q10, q10, #26             \n\t"
    "vshr.s64     q9, q9, #23               \n\t"
    "vadd.s64     q10, q13, q10             \n\t"
    "vand.u64     q12, q9, q7               \n\t"
    "vand.u64     q2, q10, q7               \n\t"

    "vshr.s64     q10, q10, #26             \n\t"
    "vshr.s64     q6, q9, #26               \n\t"
    "vadd.s64     q10, q3, q10              \n\t"
    "vadd.s64     q0, q0, q12               \n\t"
    "vadd.s64     q1, q1, q6                \n\t"
    
    "vst2.32      {d0[0],d1[0]}, [%2]!      \n\t"
    "vshr.s64     q9, q10, #26              \n\t"
    "vst2.32      {d2[0],d3[0]}, [%2]!      \n\t"
    "vand.u64     q3, q10, q7               \n\t"
    "vst2.32      {d4[0],d5[0]}, [%2]!      \n\t" 
    "vpop         {q4-q7}                   \n\t"
    "vadd.s64     q11, q11, q9              \n\t"   
    "vst2.32      {d6[0],d7[0]}, [%2]!      \n\t"
    "vst2.32      {d22[0],d23[0]}, [%2]!    \n\t"
    :
    :"r"(&a[0]), "r"(&b[0]), "r"(&c[0]), "r"(&d[0]), "r"(&e[0]), "r"(&f[0])
    :"memory","r6","r7","r8","r9"
	);
	return;
}


void v2muldblsub1271_a(uint32_t* a, uint32_t* b, uint32_t* c, uint32_t* d, uint32_t* e, uint32_t* f)  
{ // Multiplication/addition/subtraction over GF((2^127-1)^2) combining ARM and NEON instructions 
  // Operation: c = a * b, f = 2*d - e in GF((2^127-1)^2)
  // Representation: 23/23/26/26/26/26/26/26/26/26-bit

asm volatile(
    // q0-q2[0] <- A0|A1, q3-q5[0] <- B0|B1
    "vld1.8       {d0,d1}, [%0]!            \n\t"
    "vld1.8       {d2,d3}, [%0]!            \n\t"
    "vld1.8       {d4}, [%0]!               \n\t" 
    "vshl.i32     q6, q0, #3                \n\t"    // q6-q7 = 8*(a0->a8)  
    "vpush        {q4-q7}                   \n\t" 
    "vld1.8       {d6,d7}, [%1]!          	\n\t"
    "vshl.i32     d5, d4, #3                \n\t"    // q2[1] = 8*(a4|a9)   
    "vld1.8       {d8,d9}, [%1]!          	\n\t"
    "vshl.i32     q7, q1, #3                \n\t"  
    "vld1.8       {d10}, [%1]!              \n\t"

    // q11-q15 <- A0*B0|A0*B1, q8-q10 <- A1*B0|A1*B1 
    "vmull.s32    q8, d6, d0[1]             \n\t"    // q8 = a5*b0 | a5*b5
    "vmlal.s32    q8, d10, d13[1]           \n\t"    // q8 = a5*b0 + 8*(a6*b4) | a5*b5 + 8*(a6*b9)
    "vmlal.s32    q8, d7, d5[1]             \n\t"    // q8 = a5*b0 + 8*(a6*b4 + a9*b1) | a5*b5 + 8*(a6*b9 + a9*b6)
    "vmlal.s32    q8, d9, d14[1]            \n\t"    // q8 = a5*b0 + 8*(a6*b4 + a9*b1 + a7*b3) | a5*b5 + 8*(a6*b9 + a9*b6 + a7*b8)
    "vmlal.s32    q8, d8, d15[1]            \n\t"    // q8 = a5*b0 + 8*(a6*b4 + a9*b1 + a7*b3 + a8*b2) | a5*b5 + 8*(a6*b9 + a9*b6 + a7*b8 + a8*b7)
    "vmull.s32    q9, d7, d0[1]             \n\t"    // q9 = a5*b1 | a5*b6
    "vmlal.s32    q9, d6, d1[1]             \n\t"    // q9 = a5*b1 + a6*b0 | a5*b6 + a6*b5
    "vmlal.s32    q9, d10, d14[1]           \n\t"    // q9 = a5*b1 + a6*b0 + 8*(a7*b4) | a5*b6 + a6*b5 + 8*(a7*b9)
    "vmlal.s32    q9, d8, d5[1]             \n\t"    // q9 = a5*b1 + a6*b0 + 8*(a7*b4 + a9*b2) | a5*b6 + a6*b5 + 8*(a7*b9 + a9*b7)
    "vmlal.s32    q9, d9, d15[1]            \n\t"    // q9 = a5*b1 + a6*b0 + 8*(a7*b4 + a9*b2 + a8*b3) | a5*b6 + a6*b5 + 8*(a7*b9 + a9*b7 + a8*b8)
    "vmull.s32    q10, d8, d0[1]            \n\t"    // q10 = a5*b2 | a5*b7
    "vmlal.s32    q10, d6, d2[1]            \n\t"    // q10 = a5*b2 + a7*b0 | a5*b7 + a7*b5
    "vmlal.s32    q10, d7, d1[1]            \n\t"    // q10 = a5*b2 + a7*b0 + a6*b1 | a5*b7 + a7*b5 + a6*b6
    "vmlal.s32    q10, d10, d15[1]          \n\t"    // q10 = a5*b2 + a7*b0 + a6*b1 + 8*(a8*b4) | a5*b7 + a7*b5 + a6*b6 + 8*(a8*b9)
    "vmlal.s32    q10, d9, d5[1]            \n\t"    // q10 = a5*b2 + a7*b0 + a6*b1 + 8*(a8*b4 + a9*b3) | a5*b7 + a7*b5 + a6*b6 + 8*(a8*b9 + a9*b8)
         
    "push         {r6-r9}                 	\n\t"
    "ldm          %3!, {r6,r7}              \n\t"    //***** 
    "ldm          %4!, {r8,r9}              \n\t"    //***** 
    "add          r6, r6                    \n\t"    //***** 
    "add          r7, r7                    \n\t"    //***** 
    "sub          r6, r8                    \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "stm          %5!, {r6,r7}              \n\t"    //*****
    
    "vmull.s32    q11, d6, d0[0]             \n\t"    // q11 = a0*b0 | a0*b5
    "vmlal.s32    q11, d10, d13[0]           \n\t"    // q11 = a0*b0 + 8*(a1*b4) | a0*b5 + 8*(a1*b9)
    "vmlal.s32    q11, d7, d5[0]             \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1) | a0*b5 + 8*(a1*b9 + a4*b6)
    "vmlal.s32    q11, d9, d14[0]            \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3) | a0*b5 + 8*(a1*b9 + a4*b6 + a2*b8)
    "vmlal.s32    q11, d8, d15[0]            \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3 + a3*b2) | a0*b5 + 8*(a1*b9 + a4*b6 + a2*b8 + a3*b7)
    "vmull.s32    q12, d7, d0[0]             \n\t"    // q12 = a0*b1 | a0*b6
    "vmlal.s32    q12, d6, d1[0]             \n\t"    // q12 = a0*b1 + a1*b0 | a0*b6 + a1*b5
    "vmlal.s32    q12, d10, d14[0]           \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4) | a0*b6 + a1*b5 + 8*(a2*b9)
    "vmlal.s32    q12, d8, d5[0]             \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2) | a0*b6 + a1*b5 + 8*(a2*b9 + a4*b7)
    "vmlal.s32    q12, d9, d15[0]            \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2 + a3*b3) | a0*b6 + a1*b5 + 8*(a2*b9 + a4*b7 + a3*b8)
    "vmull.s32    q13, d8, d0[0]             \n\t"    // q13 = a0*b2 | a0*b7
    "vmlal.s32    q13, d6, d2[0]             \n\t"    // q13 = a0*b2 + a2*b0 | a0*b7 + a2*b5
    "vmlal.s32    q13, d7, d1[0]             \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 | a0*b7 + a2*b5 + a1*b6
    "vmlal.s32    q13, d10, d15[0]           \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4) | a0*b7 + a2*b5 + a1*b6 + 8*(a3*b9)
    "vmlal.s32    q13, d9, d5[0]             \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4 + a4*b3) | a0*b7 + a2*b5 + a1*b6 + 8*(a3*b9 + a4*b8)
    "vmull.s32    q14, d9, d0[0]             \n\t"    // q14 = a0*b3 | a0*b8
    "vmlal.s32    q14, d6, d3[0]             \n\t"    // q14 = a0*b3 + a3*b0 | a0*b8 + a3*b5
    "vmlal.s32    q14, d8, d1[0]             \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 | a0*b8 + a3*b5 + a1*b7
    "vmlal.s32    q14, d7, d2[0]             \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 + a2*b1 | a0*b8 + a3*b5 + a1*b7 + a2*b6
    "vmlal.s32    q14, d10, d5[0]            \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 + a2*b1 + 8*(a4*b4) | a0*b8 + a3*b5 + a1*b7 + a2*b6 + 8*(a4*b9)
    "vmull.s32    q15, d10, d0[0]            \n\t"    // q15 = a0*b4 | a0*b9
    "vmlal.s32    q15, d6, d4[0]             \n\t"    // q15 = a0*b4 + a4*b0 | a0*b9 + a4*b5
    "vmlal.s32    q15, d9, d1[0]             \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 | a0*b9 + a4*b5 + a1*b8
    "vmlal.s32    q15, d7, d3[0]             \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 + a3*b1 | a0*b9 + a4*b5 + a1*b8 + a3*b6
    "vmlal.s32    q15, d8, d2[0]             \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 + a3*b1 + a2*b2 | a0*b9 + a4*b5 + a1*b8 + a3*b6 + a2*b7
         
    "ldm          %3!, {r6,r7}              \n\t"    //***** 
    "ldm          %4!, {r8,r9}              \n\t"    //***** 
    "add          r6, r6                    \n\t"    //***** 
    "add          r7, r7                    \n\t"    //***** 
    "sub          r6, r8                    \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "stm          %5!, {r6,r7}              \n\t"    //*****
    
    // q11[0]->q13[0] = A0*B0 - A1*B1, q11[1]->q13[1] = A0*B0 + A1*B1 
    "vsub.s64     d22, d22, d17              \n\t"
    "vsub.s64     d24, d24, d19              \n\t"
    "vsub.s64     d26, d26, d21              \n\t"
    "vadd.s64     d23, d23, d16              \n\t"
    "vadd.s64     d25, d25, d18              \n\t"
    "vadd.s64     d27, d27, d20              \n\t"
    
    "ldm          %3!, {r6,r7}              \n\t"    //***** 
    "ldm          %4!, {r8,r9}              \n\t"    //***** 
    "add          r6, r6                    \n\t"    //***** 
    "add          r7, r7                    \n\t"    //***** 
    "sub          r6, r8                    \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "stm          %5!, {r6,r7}              \n\t"    //*****
        
    // Complete computation A1*B0|A1*B1, q8-q9 <- A1*B0|A1*B1 
    "vmull.s32    q8, d9, d0[1]             \n\t"    // q8 = a5*b3 | a5*b8
    "vmlal.s32    q8, d6, d3[1]             \n\t"    // q8 = a5*b3 + a8*b0 | a5*b8 + a8*b5
    "vmlal.s32    q8, d8, d1[1]             \n\t"    // q8 = a5*b3 + a8*b0 + a6*b2 | a5*b8 + a8*b5 + a6*b7
    "vmlal.s32    q8, d7, d2[1]             \n\t"    // q8 = a5*b3 + a8*b0 + a6*b2 + a7*b1 | a5*b8 + a8*b5 + a6*b7 + a7*b6
    "vmlal.s32    q8, d10, d5[1]            \n\t"    // q8 = a5*b3 + a8*b0 + a6*b2 + a7*b1 + 8*(a9*b4) | a5*b8 + a8*b5 + a6*b7 + a7*b6 + 8*(a9*b9)
    "vmull.s32    q9, d10, d0[1]            \n\t"    // q9 = a5*b4 | a5*b9
    "vmlal.s32    q9, d6, d4[1]             \n\t"    // q9 = a5*b4 + a9*b0 | a5*b9 + a9*b5
    "vmlal.s32    q9, d9, d1[1]             \n\t"    // q9 = a5*b4 + a9*b0 + a6*b3 | a5*b9 + a9*b5 + a6*b8
    "vmlal.s32    q9, d7, d3[1]             \n\t"    // q9 = a5*b4 + a9*b0 + a6*b3 + a8*b1 | a5*b9 + a9*b5 + a6*b8 + a8*b6
    "vmlal.s32    q9, d8, d2[1]             \n\t"    // q9 = a5*b4 + a9*b0 + a6*b3 + a8*b1 + a7*b2 | a5*b9 + a9*b5 + a6*b8 + a8*b6 + a7*b7
    
    "ldm          %3!, {r6,r7}              \n\t"    //***** 
    "ldm          %4!, {r8,r9}              \n\t"    //***** 
    "add          r6, r6                    \n\t"    //***** 
    "add          r7, r7                    \n\t"    //***** 
    "sub          r6, r8                    \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "stm          %5!, {r6,r7}              \n\t"    //*****
    
    // Complete q14[0]->q15[0] = A0*B0 - A1*B1, q14[1]->q15[0] = A0*B1 + A1*B0 
    "vmov.i64     q7, 0xFFFFFFFF            \n\t"    // mask_26
    "vadd.s64     d29, d29, d16             \n\t"
    "vsub.s64     d28, d28, d17             \n\t"
	"vshr.u64     q7, q7, #6                \n\t"
    "vadd.s64     d31, d31, d18             \n\t"
    "vsub.s64     d30, d30, d19             \n\t"
    
    "ldm          %3!, {r6,r7}              \n\t"    //***** 
    "ldm          %4!, {r8,r9}              \n\t"    //***** 
    "add          r6, r6                    \n\t"    //***** 
    "add          r7, r7                    \n\t"    //***** 
    "sub          r6, r8                    \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "stm          %5!, {r6,r7}              \n\t"    //*****
    "pop          {r6-r9}                 	\n\t"

    // Reduction   
	"vshr.u64     q6, q7, #3                \n\t"    // mask_23
    "vshr.s64     q10, q11, #26             \n\t"
    "vand.u64     q0, q11, q7               \n\t"

    "vshr.s64     q9, q14, #26              \n\t"
    "vand.u64     q3, q14, q7               \n\t"
    "vadd.s64     q10, q12, q10             \n\t"
    "vadd.s64     q9, q15, q9               \n\t"
    "vand.u64     q1, q10, q7               \n\t"
    "vand.u64     q11, q9, q6                \n\t" 

    "vshr.s64     q10, q10, #26             \n\t"
    "vshr.s64     q9, q9, #23               \n\t"
    "vadd.s64     q10, q13, q10             \n\t"
    "vand.u64     q12, q9, q7               \n\t"
    "vand.u64     q2, q10, q7               \n\t"

    "vshr.s64     q10, q10, #26             \n\t"
    "vshr.s64     q6, q9, #26               \n\t"
    "vadd.s64     q10, q3, q10              \n\t"
    "vadd.s64     q0, q0, q12               \n\t"
    "vadd.s64     q1, q1, q6                \n\t"
    
    "vst2.32      {d0[0],d1[0]}, [%2]!      \n\t"
    "vshr.s64     q9, q10, #26              \n\t"
    "vst2.32      {d2[0],d3[0]}, [%2]!      \n\t"
    "vand.u64     q3, q10, q7               \n\t"
    "vst2.32      {d4[0],d5[0]}, [%2]!      \n\t" 
    "vpop         {q4-q7}                   \n\t"
    "vadd.s64     q11, q11, q9              \n\t"   
    "vst2.32      {d6[0],d7[0]}, [%2]!      \n\t"
    "vst2.32      {d22[0],d23[0]}, [%2]!    \n\t"
    :
    :"r"(&a[0]), "r"(&b[0]), "r"(&c[0]), "r"(&d[0]), "r"(&e[0]), "r"(&f[0])
    :"memory","r6","r7","r8","r9"
	);
	return;
}


void v2sqradd1271_a(uint32_t* a, uint32_t* c, uint32_t* d, uint32_t* e, uint32_t* f) 
{ // Squaring/addition over GF((2^127-1)^2) combining ARM and NEON instructions
  // Operation: c = a^2, f = d + e in GF((2^127-1)^2)
  // Representation: 23/23/26/26/26/26/26/26/26/26-bit

asm volatile(
    // q0-q2[0] <- A0|A1
    "vld1.8       {d0,d1}, [%0]!            \n\t"
    "vld1.8       {d2,d3}, [%0]!            \n\t"
    "vshr.u64     q3, q0, #32               \n\t"    // q3-q5[0] = (0|a9->0|a5) 
    "vld1.8       {d4}, [%0]!               \n\t"
    "vpush        {q4-q7}                   \n\t" 
    "vmov.i64     q8, 0xFFFFFFFF            \n\t"    // mask_26 
    "vshr.u64     q4, q1, #32               \n\t" 
    "vshr.u64     d10, d4, #32              \n\t"    
	"vsub.s32     q13, q0, q3               \n\t"    // q13-q15[0] = (a9|a4-a9->a5|a0-a5)
	"vsub.s32     q14, q1, q4               \n\t" 
	"vsub.s32     d30, d4, d10              \n\t"     
	"vadd.s32     q3, q0, q3                \n\t"    // q3-q5[0] = (a9|a4+a9->a5|a0+a5)
	"vadd.s32     q4, q1, q4                \n\t" 
	"vadd.s32     d10, d4, d10              \n\t" 
    "vshl.i64     q0, q0, #33               \n\t"    // q0-q2[0] = (2*a4|0->2*a0|0) 
    "vshl.i64     q1, q1, #33               \n\t" 
    "vshl.i64     d4, d4, #33               \n\t"      
    "vbit         q0, q13, q8               \n\t"    // q0-q2[0] = (2*a4|a4-a9->2*a0|a0-a5)
    "vbit         q1, q14, q8               \n\t"  
    "vbit         d4, d30, d16              \n\t"     
    "vshl.i32     q6, q0, #3                \n\t"    // 8*(2*a4|a4-a9->2*a0|a0-a5)   
    "vshl.i32     q7, q1, #3                \n\t" 
    "vshl.i32     d5, d4, #3                \n\t"    
    
    "push         {r5-r8}                 	\n\t"
    "ldm          %2!, {r5,r6}              \n\t"    //***** 
    "ldm          %3!, {r7,r8}              \n\t"    //***** 
    "add          r5, r7                    \n\t"    //***** 
    "add          r6, r8                    \n\t"    //***** 
    "stm          %4!, {r5,r6}              \n\t"    //*****
        
    // q11-q15 <- 2*A0*B0|(A0+A1)*(A0-A1)  
	"vshr.u64     q8, q8, #6                \n\t"    
    "vmull.s32    q11, d6, d0               \n\t"    // q11 = a0*b0 | a0*b5
    "vmlal.s32    q11, d10, d13             \n\t"    // q11 = a0*b0 + 8*(a1*b4) | a0*b5 + 8*(a1*b9)
    "vmlal.s32    q11, d7, d5               \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1) | a0*b5 + 8*(a1*b9 + a4*b6)
    "vmlal.s32    q11, d9, d14              \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3) | a0*b5 + 8*(a1*b9 + a4*b6 + a2*b8)
    "vmlal.s32    q11, d8, d15              \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3 + a3*b2) | a0*b5 + 8*(a1*b9 + a4*b6 + a2*b8 + a3*b7)
    "vmull.s32    q12, d7, d0               \n\t"    // q12 = a0*b1 | a0*b6 
    "ldm          %2!, {r5,r6}              \n\t"    //***** 
    "ldm          %3!, {r7,r8}              \n\t"    //***** 
    "add          r5, r7                    \n\t"    //***** 
    "add          r6, r8                    \n\t"    //***** 
    "stm          %4!, {r5,r6}              \n\t"    //*****
    "vmlal.s32    q12, d6, d1               \n\t"    // q12 = a0*b1 + a1*b0 | a0*b6 + a1*b5
    "vmlal.s32    q12, d10, d14             \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4) | a0*b6 + a1*b5 + 8*(a2*b9)
    "vmlal.s32    q12, d8, d5               \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2) | a0*b6 + a1*b5 + 8*(a2*b9 + a4*b7)
    "vmlal.s32    q12, d9, d15              \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2 + a3*b3) | a0*b6 + a1*b5 + 8*(a2*b9 + a4*b7 + a3*b8)
    "vmull.s32    q13, d8, d0               \n\t"    // q13 = a0*b2 | a0*b7
    "vmlal.s32    q13, d6, d2               \n\t"    // q13 = a0*b2 + a2*b0 | a0*b7 + a2*b5
    "vmlal.s32    q13, d7, d1               \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 | a0*b7 + a2*b5 + a1*b6
    "vmlal.s32    q13, d10, d15             \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4) | a0*b7 + a2*b5 + a1*b6 + 8*(a3*b9)
    "vmlal.s32    q13, d9, d5               \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4 + a4*b3) | a0*b7 + a2*b5 + a1*b6 + 8*(a3*b9 + a4*b8)
    "vmull.s32    q14, d9, d0               \n\t"    // q14 = a0*b3 | a0*b8 
    "ldm          %2!, {r5,r6}              \n\t"    //***** 
    "ldm          %3!, {r7,r8}              \n\t"    //***** 
    "add          r5, r7                    \n\t"    //***** 
    "add          r6, r8                    \n\t"    //***** 
    "stm          %4!, {r5,r6}              \n\t"    //*****
    "vmlal.s32    q14, d6, d3               \n\t"    // q14 = a0*b3 + a3*b0 | a0*b8 + a3*b5
    "vmlal.s32    q14, d8, d1               \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 | a0*b8 + a3*b5 + a1*b7
    "vmlal.s32    q14, d7, d2               \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 + a2*b1 | a0*b8 + a3*b5 + a1*b7 + a2*b6
    "vmlal.s32    q14, d10, d5              \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 + a2*b1 + 8*(a4*b4) | a0*b8 + a3*b5 + a1*b7 + a2*b6 + 8*(a4*b9)
    "vmull.s32    q15, d10, d0              \n\t"    // q15 = a0*b4 | a0*b9
    "vmlal.s32    q15, d6, d4               \n\t"    // q15 = a0*b4 + a4*b0 | a0*b9 + a4*b5
    "vmlal.s32    q15, d9, d1               \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 | a0*b9 + a4*b5 + a1*b8
    "vmlal.s32    q15, d7, d3               \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 + a3*b1 | a0*b9 + a4*b5 + a1*b8 + a3*b6
    "vmlal.s32    q15, d8, d2               \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 + a3*b1 + a2*b2 | a0*b9 + a4*b5 + a1*b8 + a3*b6 + a2*b7
     
    "ldm          %2!, {r5,r6}              \n\t"    //***** 
    "ldm          %3!, {r7,r8}              \n\t"    //***** 
    "add          r5, r7                    \n\t"    //***** 
    "add          r6, r8                    \n\t"    //***** 
    "stm          %4!, {r5,r6}              \n\t"    //*****

    // Reduction   
	"vshr.u64     q6, q8, #3                \n\t"    // mask_23
    "vshr.s64     q10, q11, #26             \n\t"
    "vand.u64     q0, q11, q8               \n\t"

    "vshr.s64     q9, q14, #26              \n\t"
    "vand.u64     q3, q14, q8               \n\t"
    "vadd.s64     q10, q12, q10             \n\t"
    "vadd.s64     q9, q15, q9               \n\t"
    "vand.u64     q1, q10, q8               \n\t"
    "vand.u64     q11, q9, q6               \n\t" 

    "ldm          %2!, {r5,r6}              \n\t"    //***** 
    "ldm          %3!, {r7,r8}              \n\t"    //***** 
    "add          r5, r7                    \n\t"    //***** 
    "add          r6, r8                    \n\t"    //***** 
    "stm          %4!, {r5,r6}              \n\t"    //***** 
    "pop          {r5-r8}                 	\n\t"

    "vshr.s64     q10, q10, #26             \n\t"
    "vshr.s64     q9, q9, #23               \n\t"
    "vadd.s64     q10, q13, q10             \n\t"
    "vand.u64     q12, q9, q8               \n\t"
    "vand.u64     q2, q10, q8               \n\t"

    "vshr.s64     q10, q10, #26             \n\t"
    "vshr.s64     q6, q9, #26               \n\t"
    "vadd.s64     q10, q3, q10              \n\t"
    "vadd.s64     q0, q0, q12               \n\t"    
    "vadd.s64     q1, q1, q6                \n\t"
    
    "vst2.32      {d0[0],d1[0]}, [%1]!      \n\t"
    "vshr.s64     q9, q10, #26              \n\t"
    "vst2.32      {d2[0],d3[0]}, [%1]!      \n\t"
    "vand.u64     q3, q10, q8               \n\t"
    "vst2.32      {d4[0],d5[0]}, [%1]!      \n\t" 
    "vpop         {q4-q7}                   \n\t"
    "vadd.s64     q11, q11, q9              \n\t"   
    "vst2.32      {d6[0],d7[0]}, [%1]!      \n\t"
    "vst2.32      {d22[0],d23[0]}, [%1]!    \n\t"
    :
    :"r"(&a[0]), "r"(&c[0]), "r"(&d[0]), "r"(&e[0]), "r"(&f[0])
    :"memory","r5","r6","r7","r8"
	);
	return;
}


void v2sqraddsub1271_a(uint32_t* a, uint32_t* c, uint32_t* d, uint32_t* e, uint32_t* f, uint32_t* g)
{ // Squaring/addition/subtraction over GF((2^127-1)^2) combining ARM and NEON instructions
  // Operation: c = a^2, f = d + e, g = d - e in GF((2^127-1)^2)
  // Representation: 23/23/26/26/26/26/26/26/26/26-bit

asm volatile(
    // q0-q2[0] <- A0|A1
    "vld1.8       {d0,d1}, [%0]!            \n\t"
    "vld1.8       {d2,d3}, [%0]!            \n\t"
    "vshr.u64     q3, q0, #32               \n\t"    // q3-q5[0] = (0|a9->0|a5) 
    "vld1.8       {d4}, [%0]!               \n\t"
    "vpush        {q4-q7}                   \n\t" 
    "vmov.i64     q8, 0xFFFFFFFF            \n\t"    // mask_26 
    "vshr.u64     q4, q1, #32               \n\t" 
    "vshr.u64     d10, d4, #32              \n\t"    
	"vsub.s32     q13, q0, q3               \n\t"    // q13-q15[0] = (a9|a4-a9->a5|a0-a5)
	"vsub.s32     q14, q1, q4               \n\t" 
	"vsub.s32     d30, d4, d10              \n\t"     
	"vadd.s32     q3, q0, q3                \n\t"    // q3-q5[0] = (a9|a4+a9->a5|a0+a5)
	"vadd.s32     q4, q1, q4                \n\t" 
	"vadd.s32     d10, d4, d10              \n\t" 
    "vshl.i64     q0, q0, #33               \n\t"    // q0-q2[0] = (2*a4|0->2*a0|0) 
    "vshl.i64     q1, q1, #33               \n\t" 
    "vshl.i64     d4, d4, #33               \n\t"      
    "vbit         q0, q13, q8               \n\t"    // q0-q2[0] = (2*a4|a4-a9->2*a0|a0-a5)
    "vbit         q1, q14, q8               \n\t"  
    "vbit         d4, d30, d16              \n\t"     
    "vshl.i32     q6, q0, #3                \n\t"    // 8*(2*a4|a4-a9->2*a0|a0-a5)   
    "vshl.i32     q7, q1, #3                \n\t" 
    "vshl.i32     d5, d4, #3                \n\t"    
    
    "push         {r6-r11}                 	\n\t"
    "ldm          %2!, {r6,r7}              \n\t"    //***** 
    "ldm          %3!, {r8,r9}              \n\t"    //***** 
    "add          r10, r6, r8               \n\t"    //***** 
    "add          r11, r7, r9               \n\t"    //***** 
    "sub          r6, r8                    \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "stm          %4!, {r10,r11}            \n\t"    //*****
    "stm          %5!, {r6,r7}              \n\t"    //*****
        
    // q11-q15 <- 2*A0*B0|(A0+A1)*(A0-A1)  
	"vshr.u64     q8, q8, #6                \n\t"    
    "vmull.s32    q11, d6, d0               \n\t"    // q11 = a0*b0 | a0*b5
    "vmlal.s32    q11, d10, d13             \n\t"    // q11 = a0*b0 + 8*(a1*b4) | a0*b5 + 8*(a1*b9)
    "vmlal.s32    q11, d7, d5               \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1) | a0*b5 + 8*(a1*b9 + a4*b6)
    "vmlal.s32    q11, d9, d14              \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3) | a0*b5 + 8*(a1*b9 + a4*b6 + a2*b8)
    "vmlal.s32    q11, d8, d15              \n\t"    // q11 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3 + a3*b2) | a0*b5 + 8*(a1*b9 + a4*b6 + a2*b8 + a3*b7)
    "vmull.s32    q12, d7, d0               \n\t"    // q12 = a0*b1 | a0*b6 
    "ldm          %2!, {r6,r7}              \n\t"    //***** 
    "ldm          %3!, {r8,r9}              \n\t"    //***** 
    "add          r10, r6, r8               \n\t"    //***** 
    "add          r11, r7, r9               \n\t"    //***** 
    "sub          r6, r8                    \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "stm          %4!, {r10,r11}            \n\t"    //*****
    "stm          %5!, {r6,r7}              \n\t"    //*****
    "vmlal.s32    q12, d6, d1               \n\t"    // q12 = a0*b1 + a1*b0 | a0*b6 + a1*b5
    "vmlal.s32    q12, d10, d14             \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4) | a0*b6 + a1*b5 + 8*(a2*b9)
    "vmlal.s32    q12, d8, d5               \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2) | a0*b6 + a1*b5 + 8*(a2*b9 + a4*b7)
    "vmlal.s32    q12, d9, d15              \n\t"    // q12 = a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2 + a3*b3) | a0*b6 + a1*b5 + 8*(a2*b9 + a4*b7 + a3*b8)
    "vmull.s32    q13, d8, d0               \n\t"    // q13 = a0*b2 | a0*b7
    "vmlal.s32    q13, d6, d2               \n\t"    // q13 = a0*b2 + a2*b0 | a0*b7 + a2*b5
    "vmlal.s32    q13, d7, d1               \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 | a0*b7 + a2*b5 + a1*b6
    "vmlal.s32    q13, d10, d15             \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4) | a0*b7 + a2*b5 + a1*b6 + 8*(a3*b9)
    "vmlal.s32    q13, d9, d5               \n\t"    // q13 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4 + a4*b3) | a0*b7 + a2*b5 + a1*b6 + 8*(a3*b9 + a4*b8)
    "vmull.s32    q14, d9, d0               \n\t"    // q14 = a0*b3 | a0*b8 
    "ldm          %2!, {r6,r7}              \n\t"    //***** 
    "ldm          %3!, {r8,r9}              \n\t"    //***** 
    "add          r10, r6, r8               \n\t"    //***** 
    "add          r11, r7, r9               \n\t"    //***** 
    "sub          r6, r8                    \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "stm          %4!, {r10,r11}            \n\t"    //*****
    "stm          %5!, {r6,r7}              \n\t"    //*****
    "vmlal.s32    q14, d6, d3               \n\t"    // q14 = a0*b3 + a3*b0 | a0*b8 + a3*b5
    "vmlal.s32    q14, d8, d1               \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 | a0*b8 + a3*b5 + a1*b7
    "vmlal.s32    q14, d7, d2               \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 + a2*b1 | a0*b8 + a3*b5 + a1*b7 + a2*b6
    "vmlal.s32    q14, d10, d5              \n\t"    // q14 = a0*b3 + a3*b0 + a1*b2 + a2*b1 + 8*(a4*b4) | a0*b8 + a3*b5 + a1*b7 + a2*b6 + 8*(a4*b9)
    "vmull.s32    q15, d10, d0              \n\t"    // q15 = a0*b4 | a0*b9
    "vmlal.s32    q15, d6, d4               \n\t"    // q15 = a0*b4 + a4*b0 | a0*b9 + a4*b5
    "vmlal.s32    q15, d9, d1               \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 | a0*b9 + a4*b5 + a1*b8
    "vmlal.s32    q15, d7, d3               \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 + a3*b1 | a0*b9 + a4*b5 + a1*b8 + a3*b6
    "vmlal.s32    q15, d8, d2               \n\t"    // q15 = a0*b4 + a4*b0 + a1*b3 + a3*b1 + a2*b2 | a0*b9 + a4*b5 + a1*b8 + a3*b6 + a2*b7
     
    "ldm          %2!, {r6,r7}              \n\t"    //***** 
    "ldm          %3!, {r8,r9}              \n\t"    //***** 
    "add          r10, r6, r8               \n\t"    //***** 
    "add          r11, r7, r9               \n\t"    //***** 
    "sub          r6, r8                    \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "stm          %4!, {r10,r11}            \n\t"    //*****
    "stm          %5!, {r6,r7}              \n\t"    //*****

    // Reduction   
	"vshr.u64     q6, q8, #3                \n\t"    // mask_23
    "vshr.s64     q10, q11, #26             \n\t"
    "vand.u64     q0, q11, q8               \n\t"

    "vshr.s64     q9, q14, #26              \n\t"
    "vand.u64     q3, q14, q8               \n\t"
    "vadd.s64     q10, q12, q10             \n\t"
    "vadd.s64     q9, q15, q9               \n\t"
    "vand.u64     q1, q10, q8               \n\t"
    "vand.u64     q11, q9, q6               \n\t" 
    
    "ldm          %2!, {r6,r7}              \n\t"    //***** 
    "ldm          %3!, {r8,r9}              \n\t"    //***** 
    "add          r10, r6, r8               \n\t"    //***** 
    "add          r11, r7, r9               \n\t"    //***** 
    "sub          r6, r8                    \n\t"    //***** 
    "sub          r7, r9                    \n\t"    //***** 
    "stm          %4!, {r10,r11}            \n\t"    //*****
    "stm          %5!, {r6,r7}              \n\t"    //*****
    "pop          {r6-r11}                 	\n\t"

    "vshr.s64     q10, q10, #26             \n\t"
    "vshr.s64     q9, q9, #23               \n\t"
    "vadd.s64     q10, q13, q10             \n\t"
    "vand.u64     q12, q9, q8               \n\t"
    "vand.u64     q2, q10, q8               \n\t"

    "vshr.s64     q10, q10, #26             \n\t"
    "vshr.s64     q6, q9, #26               \n\t"
    "vadd.s64     q10, q3, q10              \n\t"
    "vadd.s64     q0, q0, q12               \n\t"
    "vadd.s64     q1, q1, q6                \n\t"
    
    "vst2.32      {d0[0],d1[0]}, [%1]!      \n\t"
    "vshr.s64     q9, q10, #26              \n\t"
    "vst2.32      {d2[0],d3[0]}, [%1]!      \n\t"
    "vand.u64     q3, q10, q8               \n\t"
    "vst2.32      {d4[0],d5[0]}, [%1]!      \n\t" 
    "vpop         {q4-q7}                   \n\t"
    "vadd.s64     q11, q11, q9              \n\t"   
    "vst2.32      {d6[0],d7[0]}, [%1]!      \n\t"
    "vst2.32      {d22[0],d23[0]}, [%1]!    \n\t"
    :
    :"r"(&a[0]), "r"(&c[0]), "r"(&d[0]), "r"(&e[0]), "r"(&f[0]), "r"(&g[0])
    :"memory","r6","r7","r8","r9","r10","r11"
	);
	return;
}

#endif


void v2add1271_a(uint32_t* a, uint32_t* b, uint32_t* c) 
{ // GF(p^2) addition, c = a+b in GF((2^127-1)^2)
  // Representation: 23/23/26/26/26/26/26/26/26/26-bit

asm volatile( 
    "vld1.8       {d0,d1}, [%0]!            \n\t"
    "vld1.8       {d2,d3}, [%1]!          	\n\t"    
	"vadd.s32     q0, q1                    \n\t"
    
    "vld1.8       {d2,d3}, [%0]!            \n\t"
    "vld1.8       {d4,d5}, [%1]!          	\n\t"    
	"vadd.s32     q1, q2                    \n\t"
    "vst1.64      {d0,d1}, [%2]!            \n\t"

    "vld1.8       {d0}, [%0]!               \n\t"
    "vld1.8       {d1}, [%1]!            	\n\t"    
	"vadd.s32     d0, d1                    \n\t"
    "vst1.64      {d2,d3}, [%2]!            \n\t"
    "vst1.64      {d0}, [%2]!               \n\t" 
    :
    :"r"(&a[0]), "r"(&b[0]), "r"(&c[0])
	);
	return; 
}


void v2sub1271_a(uint32_t* a, uint32_t* b, uint32_t* c) 
{ // GF(p^2) subtraction, c = a-b in GF((2^127-1)^2)
  // Representation: 23/23/26/26/26/26/26/26/26/26-bit
  
asm volatile(  
    "vld1.8       {d0,d1}, [%0]!            \n\t"
    "vld1.8       {d2,d3}, [%1]!          	\n\t"    
	"vsub.s32     q0, q1                    \n\t"
    
    "vld1.8       {d2,d3}, [%0]!            \n\t"
    "vld1.8       {d4,d5}, [%1]!          	\n\t"    
	"vsub.s32     q1, q2                    \n\t"
    "vst1.64      {d0,d1}, [%2]!            \n\t"

    "vld1.8       {d0}, [%0]!               \n\t"
    "vld1.8       {d1}, [%1]!            	\n\t"    
	"vsub.s32     d0, d1                    \n\t"
    "vst1.64      {d2,d3}, [%2]!            \n\t"
    "vst1.64      {d0}, [%2]!               \n\t" 
    :
    :"r"(&a[0]), "r"(&b[0]), "r"(&c[0])
	);
	return; 
}


void v2dblsub1271_a(uint32_t* a, uint32_t* b, uint32_t* c) 
{ // GF(p^2) addition followed by subtraction, c = 2a-c in GF((2^127-1)^2)
  // Representation: 23/23/26/26/26/26/26/26/26/26-bit
     
asm volatile( 
    "vld1.8       {d0,d1}, [%0]!            \n\t"           
    "vld1.8       {d2,d3}, [%1]!            \n\t"        	    
	"vadd.s32     q0, q0                    \n\t"     
    
    "vld1.8       {d4,d5}, [%0]!            \n\t"                       
	"vsub.s32     q0, q1                    \n\t"      
    "vld1.8       {d6,d7}, [%1]!            \n\t"         	    
	"vadd.s32     q2, q2                    \n\t"              
    "vst1.64      {d0,d1}, [%2]!            \n\t"                    
	"vsub.s32     q2, q3                    \n\t"         

    "vld1.8       {d0}, [%0]!               \n\t"              
    "vld1.8       {d1}, [%1]!               \n\t"                  
	"vadd.s32     d0, d0                    \n\t"                
    "vst1.64      {d4,d5}, [%2]!            \n\t"                         
	"vsub.s32     d0, d1                    \n\t"                    
    "vst1.64      {d0}, [%2]!               \n\t"              
    :
    :"r"(&a[0]), "r"(&b[0]), "r"(&c[0])
	);
	return; 
}


void v2addsub1271_a(uint32_t* a, uint32_t* b, uint32_t* c, uint32_t* d) 
{ // GF(p^2) addition and subtraction, c = a+b and d = a-b in GF((2^127-1)^2)
  // Representation: 23/23/26/26/26/26/26/26/26/26-bit
        
asm volatile(  
    "vld1.8       {d0,d1}, [%0]!            \n\t"            
    "vld1.8       {d2,d3}, [%1]!            \n\t"          	  
	"vsub.s32     q2, q0, q1                \n\t"                    
	"vadd.s32     q0, q1                    \n\t"     
    "vld1.8       {d2,d3}, [%0]!            \n\t"               
    "vst1.64      {d4,d5}, [%3]!            \n\t" 
    "vld1.8       {d4,d5}, [%1]!            \n\t"                   
    "vst1.64      {d0,d1}, [%2]!            \n\t"                           	  
	"vsub.s32     q3, q1, q2                \n\t"                       
	"vadd.s32     q1, q2                    \n\t"     
    "vld1.8       {d0}, [%0]!               \n\t"              
    "vst1.64      {d6,d7}, [%3]!            \n\t"  
    "vld1.8       {d1}, [%1]!               \n\t"                   
    "vst1.64      {d2,d3}, [%2]!            \n\t"                                	      
	"vsub.s32     d4, d0, d1                \n\t"                    
	"vadd.s32     d0, d1                    \n\t"                  
    "vst1.64      {d4}, [%3]!               \n\t"                   
    "vst1.64      {d0}, [%2]!               \n\t"              
    :
    :"r"(&a[0]), "r"(&b[0]), "r"(&c[0]), "r"(&d[0])
	);
	return; 
}


void table_lookup_1x8(vpoint_extproj_precomp_t* Table, vpoint_extproj_precomp_t point, unsigned int digit, unsigned int sign_mask)
{ // Constant-time table lookup to extract a point represented as (X+Y,Y-X,2Z,2dT) corresponding to extended twisted Edwards coordinates (X:Y:Z:T)
  // Inputs: sign_mask, digit, table containing 8 points
  // Output: P = sign*table[digit], where sign=1 if sign_mask=0xFF...FF and sign=-1 if sign_mask=0
    uint32_t* table = (uint32_t*)Table;
    uint32_t* P = (uint32_t*)point;
    
asm volatile(
    // q0-q11 <- table[0]
    "vpush        {q4-q7}                   \n\t" 
    "vldm         %0!, {q0-q4}              \n\t"
    "push         {r4-r6}                  	\n\t"
    "vldm         %0!, {q5-q9}              \n\t"
    "sub          %2, %2, #1                \n\t"
    "mov          r4, #6                    \n\t"

"loop0:                                     \n\t"                     
    // If digit>=0 then mask = 0x00...0 else mask = 0xFF...F
    "asr          r5, %2, #31               \n\t"
    "vdup.32      q13, r5                   \n\t"
    
    "vld1.8       {d24,d25}, [%0]!          \n\t"
    "vbif         q0, q12, q13              \n\t"
    "vld1.8       {d24,d25}, [%0]!          \n\t"
    "vbif         q1, q12, q13              \n\t"
    "vld1.8       {d24,d25}, [%0]!          \n\t"
    "vbif         q2, q12, q13              \n\t"
    "vld1.8       {d24,d25}, [%0]!          \n\t"
    "vbif         q3, q12, q13              \n\t"
    "vld1.8       {d24,d25}, [%0]!          \n\t"
    "vbif         q4, q12, q13              \n\t"
    "vld1.8       {d24,d25}, [%0]!          \n\t"
    "vbif         q5, q12, q13              \n\t"
    "vld1.8       {d24,d25}, [%0]!          \n\t"
    "vbif         q6, q12, q13              \n\t"
    "vld1.8       {d24,d25}, [%0]!          \n\t"
    "vbif         q7, q12, q13              \n\t"
    "vld1.8       {d24,d25}, [%0]!          \n\t"
    "vbif         q8, q12, q13              \n\t"
    "vld1.8       {d24,d25}, [%0]!          \n\t"
    "vbif         q9, q12, q13              \n\t"
    "sub          %2, %2, #1                \n\t"
    
    "subs         r4, r4, #1                \n\t"
    "bpl          loop0                     \n\t"
    
    // If sign_mask = 0 then choose negative of the point
    "vdup.32      q15, %3                   \n\t"
    "vmov         d28, d5                   \n\t"
    "vmov         q12, q3                   \n\t"
    "vmov         q13, q4                   \n\t"
    "vmov.i32     d20, 0xFFFFFFFF           \n\t"    // mask_26
    "vbif         d5, d0, d30               \n\t"
    "vbif         d6, d1, d30               \n\t"
    "vbif         d7, d2, d30               \n\t"
    "vbif         d8, d3, d30               \n\t"
    "vbif         d9, d4, d30               \n\t"
	"vshr.u32     d20, d20, #6              \n\t"   
    "vbif         d0, d28, d30              \n\t"
    "vbif         d1, d24, d30              \n\t"
    "vbif         d2, d25, d30              \n\t"
    "vbif         d3, d26, d30              \n\t"
    "vbif         d4, d27, d30              \n\t"
    "vstmia       %1!, {d0-d14}             \n\t"
	"vshr.u32     d21, d20, #3              \n\t"    // mask_23
	"vsub.s32     d24, d20, d15             \n\t"    // Negate t2 coordinate 
	"vsub.s32     d25, d20, d16             \n\t"   
	"vsub.s32     d26, d20, d17             \n\t"
	"vsub.s32     d27, d20, d18             \n\t" 
	"vsub.s32     d28, d21, d19             \n\t"
    "vbif         d15, d24, d30             \n\t"
    "vbif         d16, d25, d30             \n\t"
    "vbif         d17, d26, d30             \n\t"
    "vbif         d18, d27, d30             \n\t"
    "vbif         d19, d28, d30             \n\t"
    "vstmia       %1!, {d15-d19}            \n\t"
    "pop          {r4-r6}                 	\n\t"
    "vpop         {q4-q7}                   \n\t"
    :
    :"r"(&table[0]), "r"(&P[0]), "r"(digit), "r"(sign_mask)
    :"memory","r4","r5","r6"
	);
	return;
}


void vmul1271_a(uint32_t* a, uint32_t* b, uint32_t* c) 
{ // Vectorized field multiplication over GF(2^127-1)
  // Operation: c [r3] = a [r1] * b [r2]

asm volatile(
    // q0-q1[0] <- a, q3-q4[0] <- b
    "vld1.8       {d0,d1}, [%0]!            \n\t"
    "vld1.8       {d2}, [%0]!               \n\t"
    "vpush        {q4-q7}                   \n\t"  
    "vshl.i32     q6, q0, #3                \n\t"    // q6-q7[0] = 8*(a0-a4)   
    "vshl.i32     d14, d2, #3               \n\t" 
    "vld1.8       {d6,d7}, [%1]!          	\n\t"
    "vld1.8       {d8}, [%1]!            	\n\t"

    // q13-q15[0] <- a*b 
    "vmull.s32    q15, d12, d8[0]           \n\t"    // q15 = 8*(a0*b4 | a1*b4)
    "vmlal.s32    q15, d6, d14[0]           \n\t"    // q15 = 8*(a0*b4 + a4*b0 | a1*b4 + a4*b1)
    "vmlal.s32    q15, d7, d13[0]           \n\t"    // q15 = 8*(a0*b4 + a4*b0 + a2*b2 | a1*b4 + a4*b1 + a2*b3)
    "vmull.s32    q13, d13, d8[0]           \n\t"    // q13 = 8*(a2*b4 | a3*b4)
    "vmlal.s32    q13, d7, d14[0]           \n\t"    // q13 = 8*(a2*b4 + a4*b2 | a3*b4 + a4*b3)
    "vmlal.s32    q13, d6, d0[1]            \n\t"    // q13 = a1*b0 + 8*(a2*b4 + a4*b2) | a1*b1 + 8*(a3*b4 + a4*b3)
    "vmull.s32    q14, d8, d14[0]           \n\t"    // q14 = 8*a4*b4 | x 
	"vshr.s64     d30, d30, #3              \n\t"    // q15 = a0*b4 + a4*b0 + a2*b2 | 8*(a1*b4 + a4*b1 + a2*b3)
    "vmov         d29, d30                  \n\t"  
    "vmlal.s32    q14, d7, d0[1]            \n\t"    // q14 = a1*b2 + 8*(a4*b4) | a1*b3
    "vmlal.s32    q14, d1, d6[1]            \n\t"    // q14 = a1*b2 + a2*b1 + 8*(a4*b4) | a0*b4 + a4*b0 + a1*b3 + a3*b1 + a2*b2
    "vmov         d30, d29                  \n\t"  
    "vext.s32     q14, q13, q14, #2         \n\t"  
    "vext.s32     q13, q15, q13, #2         \n\t"    // result: q13-q15[0] 
    "vmlal.s32    q13, d6, d0[0]	        \n\t"    // q13 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3) | a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2)
    "vmlal.s32    q13, d7, d13[1]	       	\n\t"    // q13 = a0*b0 + 8*(a1*b4 + a4*b1 + a2*b3 + a3*b2) | a0*b1 + a1*b0 + 8*(a2*b4 + a4*b2 + a3*b3)
    "vmlal.s32    q14, d7, d0[0]	        \n\t"    // q14 = a0*b2 + a1*b1 + 8*(a3*b4 + a4*b3) | a0*b3 + a1*b2 + a2*b1 + 8*(a4*b4)
    "vmlal.s32    q14, d1, d6[0]	        \n\t"    // q14 = a0*b2 + a2*b0 + a1*b1 + 8*(a3*b4 + a4*b3) | a0*b3 + a3*b0 + a1*b2 + a2*b1 + 8*(a4*b4)
        
    // Reduction
    "vmov.i64     q9, 0xFFFFFFFF            \n\t"    // mask_26
	"vshr.u64     q9, q9, #6                \n\t"
    "vmov.i64     d31, #0                   \n\t"
    "vswp         d29, d27                  \n\t"
    "vswp         d29, d28                  \n\t"  
    "vshr.s64     q10, q13, #26             \n\t"
    "vswp         d30, d29                  \n\t"
    "vadd.s64     q10, q14, q10             \n\t"  
    "vand.u64     q0, q13, q9               \n\t"
    "vshr.s64     q8, q10, #26              \n\t"
    "vand.u64     q1, q10, q9               \n\t"
    "vadd.s64     q8, q15, q8               \n\t"
    "vand.u64     q2, q8, q9                \n\t"
	"vshl.i64     d22, d17, #3              \n\t"    
    "vshr.s64     d23, d16, #26             \n\t"
    "vand.u64     q5, q11, q9               \n\t" 
    "vshr.s64     q6, q11, #26              \n\t"
    "vadd.s64     q0, q0, q5                \n\t"
    "vadd.s64     q1, q1, q6                \n\t"

    "vst2.32      {d0[0],d2[0]}, [%2]!      \n\t"
    "vst1.32      {d4[0]}, [%2]!            \n\t"
    "vst2.32      {d1[0],d3[0]}, [%2]!      \n\t"
    "vpop         {q4-q7}                   \n\t"
    :
    :"r"(&a[0]), "r"(&b[0]), "r"(&c[0])
	);
	return;
}


void vsqr1271_a(uint32_t* a, uint32_t* c) 
{ // Vectorized field squaring over GF(2^127-1)
  // Operation: c [r2] = a^2 [r1]

asm volatile( 
    // q0-q1[0] <- a
    "vld1.8       {d0,d1}, [%0]!            \n\t"
    "vld1.8       {d2}, [%0]!               \n\t"
    "vext.s32     q2, q0, q1, #1            \n\t"    // q2: a1|a2|a3|a4 
    "vpush        {q4-q7}                   \n\t" 

    // q13-q15[0] <- a^2 
    "vqdmull.s32  q14, d4, d0[1]            \n\t"    // q14 = 2*(a1^2 | a1*a2)
    "vext.s32     d6, d5, d5, #1            \n\t"     
    "vqdmull.s32  q15, d6, d0               \n\t"    // q15 = 2*(a0*a4 | a1*a3)
    "vshl.i32     d2, d2, #1                \n\t"    // d2: 2*a4|x
	"vshr.s64     d28, d28, #1              \n\t"    // q14 = a1^2 | 2*a1*a2
    "vadd.s64     d30, d30, d31             \n\t"    // q15 = 2*(a0*a4 + a1*a3) | x
    "vtrn.32      d2, d6                    \n\t"    // d2: 2*a4|a4    
    "vqdmull.s32  q13, d0, d0[0]            \n\t"    // q13 = 2*(a0^2 | a0*a1)
    "vshl.i32     q2, q2, #3                \n\t"    // q2: 8*(a1|a2|a3|a4) 
    "vmlal.s32    q15, d1, d1               \n\t"    // q15 = 2*(a0*a4 + a1*a3) + a2^2 | x
    "vmov.i64     q9, 0xFFFFFFFF            \n\t"    // mask_26
	"vshr.s64     d26, d26, #1              \n\t"    // q13 = a0^2 | 2*a0*a1
    "vmlal.s32    q13, d4, d2[0]            \n\t"    // q13 = a0^2 + 16*a1*a4 | 2*a0*a1 +16*a2*a4
    "vshr.s32     d4, d4, #2                \n\t"    // q2: 2*(a1|a2)|8*(a3|a4) 
    "vmlal.s32    q14, d5, d2               \n\t"    // q14 = a1^2 + 16*a3*a4 | 2*a1*a2 + 8*a4^2
    "vqdmlal.s32  q14, d1, d0[0]            \n\t"    // q14 = 2*a0*a2 + a1^2 + 16*a3*a4 | 2*a0*a3 + 2*a1*a2 + 8*a4^2
    "vtrn.32      d4, d1                    \n\t"    // d1: 2*a2|a3   
    "vmlal.s32    q13, d1, d5[0]            \n\t"    // q13 = a0^2 + 16*a1*a4 + 16*a2*a3 | 2*a0*a1 + 16*a2*a4 + 8*a3^2
    
    // Reduction
	"vshr.u64     q9, q9, #6                \n\t"
    "vmov.i64     d31, #0                   \n\t"
    "vswp         d29, d27                  \n\t"
    "vswp         d29, d28                  \n\t"  
    "vshr.s64     q10, q13, #26             \n\t"
    "vswp         d30, d29                  \n\t"
    "vadd.s64     q10, q14, q10             \n\t"  
    "vand.u64     q0, q13, q9               \n\t"
    "vshr.s64     q8, q10, #26              \n\t"
    "vand.u64     q1, q10, q9               \n\t"
    "vadd.s64     q8, q15, q8               \n\t"
    "vand.u64     q2, q8, q9                \n\t"
	"vshl.i64     d22, d17, #3              \n\t"    
    "vshr.s64     d23, d16, #26             \n\t"
    "vand.u64     q5, q11, q9               \n\t" 
    "vshr.s64     q6, q11, #26              \n\t"
    "vadd.s64     q0, q0, q5                \n\t"
    "vadd.s64     q1, q1, q6                \n\t"

    "vst2.32      {d0[0],d2[0]}, [%1]!      \n\t"
    "vst1.32      {d4[0]}, [%1]!            \n\t"
    "vst2.32      {d1[0],d3[0]}, [%1]!      \n\t"
    "vpop         {q4-q7}                   \n\t"
    :
    :"r"(&a[0]), "r"(&c[0])
	);
	return;
}


void mul_truncate_a(uint32_t* a, uint32_t* b, uint32_t* c) 
{ // 256-bit multiplication with truncation for the scalar decomposition
  // Outputs 64-bit value c = (uint64_t)((a * b) >> 256).

asm volatile(   
    "push         {r4-r12}                 	\n\t"

    "ldm          %0!, {r3-r4}              \n\t"    
    "ldr          r5, [%1, #0]              \n\t"       
    "mov          r8, #0                    \n\t" 
    "sub          r13, r13, #12             \n\t"    // Allocating space in the stack
    "umull        r7, r6, r5, r3            \n\t" 
    "umlal        r6, r8, r5, r4            \n\t"
    "ldm          %0!, {r3-r4}              \n\t"      
    "mov          r9, #0                    \n\t"     
    "mov          r10, #0                   \n\t" 
    "str          r6, [r13], #4             \n\t"    // Store in stack  
    "umlal        r8, r9, r5, r3            \n\t"
    "umlal        r9, r10, r5, r4           \n\t"
    "ldm          %0!, {r3-r4}              \n\t"      
    "mov          r11, #0                   \n\t"       
    "mov          r12, #0                   \n\t"  
    "umlal        r10, r11, r5, r3          \n\t"
    "umlal        r11, r12, r5, r4          \n\t"
    "ldm          %0!, {r3-r4}              \n\t"     
    "mov          r6, #0                    \n\t"     
    "mov          r7, #0                    \n\t" 
    "sub          %0, %0, #32               \n\t"  
    "umlal        r12, r6, r5, r3           \n\t"
    "umlal        r6, r7, r5, r4            \n\t"
    "stm          r13, {r6-r7}              \n\t"    // Store in stack 
  
    "sub          r13, r13, #4              \n\t" 
    "ldm          %0!, {r3-r4}              \n\t"    
    "ldr          r5, [%1, #4]              \n\t"    
    "ldr          r6, [r13]                 \n\t"      
    "mov          r7, #0                    \n\t"   
    "umlal        r6, r7, r5, r3            \n\t" 
    "umaal        r7, r8, r5, r4            \n\t"
    "ldm          %0!, {r3-r4}              \n\t" 
    "str          r7, [r13], #4             \n\t"    // Store in stack
    "umaal        r8, r9, r5, r3            \n\t"
    "umaal        r9, r10, r5, r4           \n\t"
    "ldm          %0!, {r3-r4}              \n\t"  
    "umaal        r10, r11, r5, r3          \n\t"
    "umaal        r11, r12, r5, r4          \n\t"
    "ldm          %0!, {r3-r4}              \n\t"    
    "ldm          r13, {r6-r7}              \n\t" 
    "sub          %0, %0, #32               \n\t"
    "umaal        r12, r6, r5, r3           \n\t"
    "umaal        r6, r7, r5, r4            \n\t" 
    "stm          r13, {r6-r7}              \n\t"    // Store in stack 
   
    "sub          r13, r13, #4              \n\t" 
    "ldm          %0!, {r3-r4}              \n\t"    
    "ldr          r5, [%1, #8]              \n\t"    
    "ldr          r7, [r13]                 \n\t"        
    "mov          r6, #0                    \n\t"    
    "umlal        r7, r6, r5, r3            \n\t" 
    "umaal        r6, r8, r5, r4            \n\t"
    "ldm          %0!, {r3-r4}              \n\t" 
    "str          r6, [r13], #4             \n\t"    // Store in stack
    "umaal        r8, r9, r5, r3            \n\t"
    "umaal        r9, r10, r5, r4           \n\t"
    "ldm          %0!, {r3-r4}              \n\t"  
    "umaal        r10, r11, r5, r3          \n\t"
    "umaal        r11, r12, r5, r4          \n\t"
    "ldm          %0!, {r3-r4}              \n\t" 
    "ldm          r13, {r6-r7}              \n\t" 
    "sub          %0, %0, #32               \n\t"    
    "umaal        r12, r6, r5, r3           \n\t"
    "umaal        r6, r7, r5, r4            \n\t"
    "stm          r13, {r6-r7}              \n\t"    // Store in stack 
  
    "sub          r13, r13, #4              \n\t" 
    "ldm          %0!, {r3-r4}              \n\t"    
    "ldr          r5, [%1, #12]             \n\t"    
    "ldr          r6, [r13]                 \n\t"          
    "mov          r7, #0                    \n\t"    
    "umlal        r6, r7, r5, r3            \n\t" 
    "umaal        r7, r8, r5, r4            \n\t"
    "ldm          %0!, {r3-r4}              \n\t" 
    "str          r7, [r13], #4             \n\t"    // Store in stack
    "umaal        r8, r9, r5, r3            \n\t"
    "umaal        r9, r10, r5, r4           \n\t"
    "ldm          %0!, {r3-r4}              \n\t"  
    "umaal        r10, r11, r5, r3          \n\t"
    "umaal        r11, r12, r5, r4          \n\t"
    "ldm          %0!, {r3-r4}              \n\t"  
    "ldm          r13, {r6-r7}              \n\t"  
    "sub          %0, %0, #32               \n\t" 
    "umaal        r12, r6, r5, r3           \n\t"   

    "sub          r13, r13, #4              \n\t"    
    "ldm          %0!, {r3-r4}              \n\t"    
    "ldr          r5, [%1, #16]             \n\t"   
    "ldr          r7, [r13]                 \n\t"           
    "mov          r6, #0                    \n\t"      
    "umlal        r7, r6, r5, r3            \n\t" 
    "umaal        r6, r8, r5, r4            \n\t"
    "ldm          %0!, {r3-r4}              \n\t"
    "str          r6, [r13]                 \n\t"    // Store in stack 
    "umaal        r8, r9, r5, r3            \n\t"
    "umaal        r9, r10, r5, r4           \n\t"
    "ldm          %0!, {r3-r4}              \n\t" 
    "sub          %0, %0, #24               \n\t" 
    "umaal        r10, r11, r5, r3          \n\t"
    "umaal        r11, r12, r5, r4          \n\t"    
  
    "ldm          %0!, {r3-r4}              \n\t"    
    "ldr          r5, [%1, #20]             \n\t"    
    "ldr          r6, [r13], #12            \n\t"             
    "mov          r7, #0                    \n\t"     
    "umlal        r6, r7, r5, r3            \n\t" 
    "umaal        r7, r8, r5, r4            \n\t"
    "ldm          %0!, {r3-r4}              \n\t" 
    "umaal        r8, r9, r5, r3            \n\t"
    "umaal        r9, r10, r5, r4           \n\t"
    "ldm          %0!, {r3-r4}              \n\t" 
    "sub          %0, %0, #24               \n\t"  
    "umaal        r10, r11, r5, r3          \n\t"        
  
    "ldm          %0!, {r3-r4}              \n\t"    
    "ldr          r5, [%1, #24]             \n\t" 
    "ldm          %0!, {r11-r12}            \n\t"        
    "mov          r6, #0                    \n\t"    
    "umlal        r7, r6, r5, r3            \n\t" 
    "umaal        r6, r8, r5, r4            \n\t"
    "umaal        r8, r9, r5, r11           \n\t"
    "umaal        r9, r10, r5, r12          \n\t"     
    "stm          %2!, {r8-r9}              \n\t" 
    "pop          {r4-r12}                 	\n\t"
    :
    :"r"(&a[0]), "r"(&b[0]), "r"(&c[0])
    :"memory","r3","r4","r5","r6","r7","r8","r9","r10","r11","r12"
	);
	return; 
}


void table_lookup_fixed_base(vpoint_precomp_t* Table, vpoint_precomp_t point, unsigned int digit, unsigned int sign)
{ // Constant-time table lookup to extract a point represented as (x+y,y-x,2t) corresponding to extended twisted Edwards coordinates (X:Y:Z:T) with Z=1
  // Inputs: sign, digit, table containing VPOINTS_FIXEDBASE = 2^(W_FIXEDBASE-1) points
  // Output: if sign=0 then P = table[digit], else if (sign=-1) then P = -table[digit]
    uint32_t* table = (uint32_t*)Table;
    uint32_t* P = (uint32_t*)point;
    uint32_t npoints = VPOINTS_FIXEDBASE;
    
asm volatile(
    // q0-q8 <- table[0]
    "vpush        {q4-q7}                   \n\t" 
    "vldm         %0!, {q0-q4}              \n\t"
    "push         {r5-r6}                  	\n\t"
    "vldm         %0!, {q5-q6}              \n\t"
    "vld1.8       {d14}, [%0]!              \n\t"
    "sub          %4, %4, #2                \n\t"

"loop1:                                     \n\t"                   
    // If digit>=0 then mask = 0x00...0 else mask = 0xFF...F
    "sub          %2, %2, #1                \n\t"
    "asr          r5, %2, #31               \n\t"
    "vdup.32      q14, r5                   \n\t"
    
    "vld1.8       {d30,d31}, [%0]!          \n\t"
    "vbif         q0, q15, q14              \n\t"
    "vld1.8       {d30,d31}, [%0]!          \n\t"
    "vbif         q1, q15, q14              \n\t"
    "vld1.8       {d30,d31}, [%0]!          \n\t"
    "vbif         q2, q15, q14              \n\t"
    "vld1.8       {d30,d31}, [%0]!          \n\t"
    "vbif         q3, q15, q14              \n\t"
    "vld1.8       {d30,d31}, [%0]!          \n\t"
    "vbif         q4, q15, q14              \n\t"
    "vld1.8       {d30,d31}, [%0]!          \n\t"
    "vbif         q5, q15, q14              \n\t"
    "vld1.8       {d30,d31}, [%0]!          \n\t"
    "vbif         q6, q15, q14              \n\t"
    "vld1.8       {d30}, [%0]!              \n\t"
    "vbif         d14, d30, d28             \n\t"
    
    "subs         %4, %4, #1                \n\t"
    "bpl          loop1                     \n\t"
    
    // If sign = 0 then choose positive of the point
    "vdup.32      q15, %3                   \n\t"    
    "vmov         d28, d5                   \n\t"
    "vmov         q12, q3                   \n\t"
    "vmov         q13, q4                   \n\t"
    "vbit         d5, d0, d30               \n\t"
    "vbit         d6, d1, d30               \n\t"
    "vbit         d7, d2, d30               \n\t"
    "vbit         d8, d3, d30               \n\t"
    "vbit         d9, d4, d30               \n\t"  
    "vbit         d0, d28, d30              \n\t"
    "vbit         d1, d24, d30              \n\t"
    "vbit         d2, d25, d30              \n\t"
    "vbit         d3, d26, d30              \n\t"
    "vbit         d4, d27, d30              \n\t"
    "vstmia       %1!, {d0-d9}              \n\t"
    "vmov.i32     q0, 0xFFFFFFFF            \n\t"    // mask_26
	"vshr.u32     q0, q0, #6                \n\t"
	"vsub.s32     q12, q0, q5               \n\t"    // Negate t2 coordinate    
	"vsub.s32     q13, q0, q6               \n\t"  
	"vshr.u32     q0, q0, #3                \n\t"    // mask_23
	"vsub.s32     d28, d0, d14              \n\t"  
    "vbit         q5, q12, q15              \n\t"
    "vbit         q6, q13, q15              \n\t"
    "vbit         d14, d28, d30             \n\t"
    "vstmia       %1!, {d10-d14}            \n\t"
    "pop          {r5-r6}                 	\n\t"
    "vpop         {q4-q7}                   \n\t"
    :
    :"r"(&table[0]), "r"(&P[0]), "r"(digit), "r"(sign), "r"(npoints)
    :"memory","r5","r6"
	);
	return;
}