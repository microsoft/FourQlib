//***********************************************************************************
// FourQlib: a high-performance crypto library based on the elliptic curve FourQ
//
//   Copyright (c) Microsoft Corporation. All rights reserved.
//
// Abstract: arithmetic over GF(p^2) using ARM assembly for Linux
//***********************************************************************************


.text
.align	4  
//***********************************************************************
//  Multiplication over GF(p^2)
//  Operation: c [reg_p3] = a [reg_p1] * b [reg_p2]
//  Restriction : b != c
//*********************************************************************** 
.global fp2mul1271_a
.type	fp2mul1271_a, %function
fp2mul1271_a:
  push   {r4-r12}
  sub    r13, r13, #44        // Allocating space in the stack
  
  // T0 = a0 * b0
  ldm    r0!, {r3-r6}
  ldr    r7, [r1]    
  mov    r10, #0
  mov    r11, #0
  mov    r12, #0
  umull  r9, r8, r7, r3      
  umlal  r8, r10, r7, r4
  umlal  r10, r11, r7, r5
  umlal  r11, r12, r7, r6
  ldr    r7, [r1, #4]
  str    r9, [r13], #4        // Store in stack
  mov    r9, #0
  
  umlal  r8, r9, r7, r3      
  umaal  r9, r10, r7, r4 
  umaal  r10, r11, r7, r5 
  umaal  r11, r12, r7, r6 
  ldr    r7, [r1, #8]
  str    r8, [r13], #4        // Store in stack
  mov    r8, #0
  
  umlal  r9, r8, r7, r3      
  umaal  r8, r10, r7, r4 
  umaal  r10, r11, r7, r5 
  umaal  r11, r12, r7, r6 
  ldr    r7, [r1, #12]
  str    r9, [r13], #4        // Store in stack
  mov    r9, #0
  
  umlal  r8, r9, r7, r3      
  umaal  r9, r10, r7, r4 
  umaal  r10, r11, r7, r5 
  umaal  r11, r12, r7, r6
  
  // T1 = a1 * b1
  ldm    r0!, {r3-r6}
  ldr    r7, [r1, #16]   
  stmia  r13!, {r8-r12}       // Store in stack  
  mov    r10, #0
  mov    r11, #0
  mov    r12, #0
  umull  r9, r8, r7, r3      
  umlal  r8, r10, r7, r4
  umlal  r10, r11, r7, r5
  umlal  r11, r12, r7, r6
  ldr    r7, [r1, #20]
  str    r9, [r13], #4        // Store in stack
  mov    r9, #0
  
  umlal  r8, r9, r7, r3      
  umaal  r9, r10, r7, r4 
  umaal  r10, r11, r7, r5 
  umaal  r11, r12, r7, r6 
  ldr    r7, [r1, #24]
  str    r8, [r13], #4        // Store in stack
  mov    r8, #0
  
  umlal  r9, r8, r7, r3      
  umaal  r8, r10, r7, r4 
  umaal  r10, r11, r7, r5 
  umaal  r11, r12, r7, r6 
  ldr    r7, [r1, #28]
  str    r9, [r13], #4        // Store in stack
  mov    r9, #0
   
  sub    r13, r13, #52 
  stmia  r13!, {r0-r1}
  umlal  r8, r9, r7, r3      
  umaal  r9, r10, r7, r4 
  umaal  r10, r11, r7, r5 
  umaal  r11, r12, r7, r6 

  // c0 = T0 - T1 = a0*b0 - a1*b1
  ldmia  r13!, {r3-r6}
  ldr    r0, [r13, #16]
  ldr    r1, [r13, #20]
  ldr    r7, [r13, #24]
  subs   r3, r3, r0
  sbcs   r4, r4, r1 
  sbcs   r5, r5, r7   
  sbcs   r6, r6, r8
  ldmia  r13!, {r7-r8}
  ldmia  r13!, {r0-r1}
  sbcs   r7, r7, r9
  sbcs   r8, r8, r10
  sbcs   r9, r0, r11
  sbcs   r10, r1, r12
  sub    r13, #40 
  pop    {r0-r1}

  // Reducing and storing c0
  lsl    r10, r10, #1
  orr    r10, r10, r9, lsr #31
  lsl    r9, r9, #1
  orr    r9, r9, r8, lsr #31
  lsl    r8, r8, #1
  orr    r8, r8, r7, lsr #31
  lsl    r7, r7, #1
  orr    r7, r7, r6, lsr #31    
  and    r6, r6, #0x7FFFFFFF
  
  subs   r7, r7, r10, lsr #31        
  sbcs   r8, r8, #0                 
  sbcs   r9, r9, #0 
  sbcs   r10, r10, #0
  and    r10, r10, #0x7FFFFFFF 

  adds   r3, r3, r7
  adcs   r4, r4, r8                 
  adcs   r5, r5, r9 
  adcs   r6, r6, r10
  
  and    r7, r6, #0x7FFFFFFF 
  adds   r3, r3, r6, lsr #31        
  adcs   r4, r4, #0                 
  adcs   r5, r5, #0 
  adcs   r6, r7, #0  
  stmia  r13!, {r3-r6}
  
  // T2 = a0 * b1
  sub    r0, r0, #32
  ldm    r0!, {r3-r6}
  ldr    r7, [r1, #16]    
  mov    r10, #0
  mov    r11, #0
  mov    r12, #0
  umull  r9, r8, r7, r3      
  umlal  r8, r10, r7, r4
  umlal  r10, r11, r7, r5
  umlal  r11, r12, r7, r6
  ldr    r7, [r1, #20]
  str    r9, [r13], #4         // Store in stack
  mov    r9, #0
  
  umlal  r8, r9, r7, r3      
  umaal  r9, r10, r7, r4 
  umaal  r10, r11, r7, r5 
  umaal  r11, r12, r7, r6 
  ldr    r7, [r1, #24]
  str    r8, [r13], #4        // Store in stack
  mov    r8, #0
  
  umlal  r9, r8, r7, r3      
  umaal  r8, r10, r7, r4 
  umaal  r10, r11, r7, r5 
  umaal  r11, r12, r7, r6 
  ldr    r7, [r1, #28]
  str    r9, [r13], #4        // Store in stack
  mov    r9, #0
  
  umlal  r8, r9, r7, r3      
  umaal  r9, r10, r7, r4 
  umaal  r10, r11, r7, r5 
  umaal  r11, r12, r7, r6 
  str    r8, [r13], #4
  
  sub    r13, r13, #32
  ldm    r13, {r3-r6}
  stmia  r13!, {r9-r12}       // Store in stack 
  ldr    r7, [r1]    
  stmia  r2!, {r3-r6}         // Storing c0  
    
  // T3 = a1 * b0
  add    r13, r13, #16
  ldm    r0!, {r3-r6}
  mov    r10, #0
  mov    r11, #0
  mov    r12, #0
  umull  r9, r8, r7, r3      
  umlal  r8, r10, r7, r4
  umlal  r10, r11, r7, r5
  umlal  r11, r12, r7, r6
  ldr    r7, [r1, #4]
  str    r9, [r13], #4        // Store in stack
  mov    r9, #0
  
  umlal  r8, r9, r7, r3      
  umaal  r9, r10, r7, r4 
  umaal  r10, r11, r7, r5 
  umaal  r11, r12, r7, r6
  ldr    r7, [r1, #8] 
  str    r8, [r13], #4        // Store in stack
  mov    r8, #0
  
  umlal  r9, r8, r7, r3      
  umaal  r8, r10, r7, r4 
  umaal  r10, r11, r7, r5 
  umaal  r11, r12, r7, r6 
  ldr    r7, [r1, #12]
  str    r9, [r13], #4        // Store in stack
  mov    r9, #0
  
  umlal  r8, r9, r7, r3      
  umaal  r9, r10, r7, r4 
  umaal  r10, r11, r7, r5 
  umaal  r11, r12, r7, r6 
    
  // c1 = T2 + T3 = a0*b1 + a1*b0 
  sub    r13, r13, #28    
  ldmia  r13!, {r3-r6}
  ldmia  r13, {r0-r1,r7}
  sub    r13, r13, #32
  adds   r3, r3, r0
  adcs   r4, r4, r1
  adcs   r5, r5, r7   
  adcs   r6, r6, r8
  ldmia  r13!, {r7-r8}
  ldmia  r13!, {r0-r1}
  adcs   r7, r7, r9
  adcs   r8, r8, r10
  adcs   r9, r0, r11
  adcs   r10, r1, r12
    
  // Reducing and storing c1
  lsl    r10, r10, #1
  orr    r10, r10, r9, lsr #31
  lsl    r9, r9, #1
  orr    r9, r9, r8, lsr #31
  lsl    r8, r8, #1
  orr    r8, r8, r7, lsr #31
  lsl    r7, r7, #1
  orr    r7, r7, r6, lsr #31    
  and    r6, r6, #0x7FFFFFFF
  
  lsls   r10, r10, #1
  adcs   r3, r3, r7      
  adcs   r4, r4, r8                 
  adcs   r5, r5, r9
  adcs   r6, r6, r10, lsr #1 
  
  and    r7, r6, #0x7FFFFFFF 
  adds   r3, r3, r6, lsr #31        
  adcs   r4, r4, #0                 
  adcs   r5, r5, #0 
  adcs   r6, r7, #0    
  stmia  r2!, {r3-r6}         // Storing c1 
  
  add    r13, r13, #28        // Restoring stack  
  pop    {r4-r12}
  bx     lr


//***********************************************************************
//  Squaring over GF(p^2)
//  Operation: c [reg_p2] = a^2 [reg_p1]
//*********************************************************************** 
.global fp2sqr1271_a
.type	fp2sqr1271_a, %function
fp2sqr1271_a:
  push   {r4-r12}
  sub    r13, r13, #16        // Allocating space in the stack
  
  // t0 = a0 + a1
  ldm    r0, {r2-r9}
  adds   r2, r2, r6
  adcs   r3, r3, r7
  adcs   r4, r4, r8
  adcs   r5, r5, r9
  stm    r13, {r2-r5}         // Store in stack
  
  // t1 = a0 - a1
  ldm    r0, {r2-r5}
  subs   r2, r2, r6
  sbcs   r3, r3, r7
  sbcs   r4, r4, r8
  sbcs   r5, r5, r9 
  
  and    r6, r5, #0x7FFFFFFF 
  subs   r2, r2, r5, lsr #31        
  sbcs   r3, r3, #0                 
  sbcs   r4, r4, #0 
  sbcs   r5, r6, #0
  
  // T0 = t0 * t1 = (a0 + a1)*(a0 - a1)
  ldr    r12, [r13]    
  mov    r9, #0
  mov    r10, #0
  mov    r11, #0
  umull  r6, r7, r12, r2      
  umlal  r7, r9, r12, r3
  umlal  r9, r10, r12, r4
  umlal  r10, r11, r12, r5
  ldr    r12, [r13, #4]
  str    r6, [r13]            // Store in stack
  mov    r6, #0
  
  umlal  r7, r6, r12, r2      
  umaal  r6, r9, r12, r3 
  umaal  r9, r10, r12, r4 
  umaal  r10, r11, r12, r5 
  ldr    r12, [r13, #8]
  str    r7, [r13, #4]        // Store in stack
  mov    r7, #0
  
  umlal  r6, r7, r12, r2      
  umaal  r7, r9, r12, r3 
  umaal  r9, r10, r12, r4 
  umaal  r10, r11, r12, r5 
  ldr    r12, [r13, #12]
  mov    r8, #0
  
  umlal  r7, r8, r12, r2    
  umaal  r8, r9, r12, r3 
  umaal  r9, r10, r12, r4 
  umaal  r10, r11, r12, r5 
  
  // Reduction
  lsl    r11, r11, #1
  orr    r11, r11, r10, lsr #31
  lsl    r10, r10, #1
  orr    r10, r10, r9, lsr #31
  lsl    r9, r9, #1
  orr    r9, r9, r8, lsr #31
  lsl    r8, r8, #1
  orr    r8, r8, r7, lsr #31 
  ldm    r13, {r4,r5}  
  and    r7, r7, #0x7FFFFFFF
  
  lsls   r11, r11, #1 
  adcs   r4, r4, r8      
  adcs   r5, r5, r9                 
  adcs   r6, r6, r10
  adcs   r7, r7, r11, lsr #1 
  
  and    r8, r7, #0x7FFFFFFF 
  adds   r9, r4, r7, lsr #31        
  adcs   r10, r5, #0 
  ldm    r0, {r2-r5}                
  adcs   r11, r6, #0 
  adcs   r12, r8, #0 
  
  // t0 = 2*a0
  adds   r2, r2, r2
  adcs   r3, r3, r3
  adcs   r4, r4, r4
  adcs   r5, r5, r5
  
  stmia  r1!, {r9-r12}        // Storing c0
  
  // T1 = 2a0 * a1
  ldr    r12, [r0, #16]    
  mov    r9, #0
  mov    r10, #0
  mov    r11, #0
  umull  r6, r7, r12, r2      
  umlal  r7, r9, r12, r3
  umlal  r9, r10, r12, r4
  umlal  r10, r11, r12, r5
  ldr    r12, [r0, #20]
  str    r6, [r13]            // Store in stack
  mov    r6, #0
  
  umlal  r7, r6, r12, r2      
  umaal  r6, r9, r12, r3 
  umaal  r9, r10, r12, r4 
  umaal  r10, r11, r12, r5 
  ldr    r12, [r0, #24]
  str    r7, [r13, #4]        // Store in stack
  mov    r7, #0
  
  umlal  r6, r7, r12, r2      
  umaal  r7, r9, r12, r3 
  umaal  r9, r10, r12, r4 
  umaal  r10, r11, r12, r5 
  ldr    r12, [r0, #28]
  mov    r8, #0
  
  umlal  r7, r8, r12, r2    
  umaal  r8, r9, r12, r3 
  umaal  r9, r10, r12, r4 
  umaal  r10, r11, r12, r5 
  
  // Reduction
  lsl    r11, r11, #1
  orr    r11, r11, r10, lsr #31
  lsl    r10, r10, #1
  orr    r10, r10, r9, lsr #31
  lsl    r9, r9, #1
  orr    r9, r9, r8, lsr #31
  lsl    r8, r8, #1
  orr    r8, r8, r7, lsr #31 
  ldm    r13, {r4,r5}  
  and    r7, r7, #0x7FFFFFFF
  
  lsls   r11, r11, #1
  adcs   r4, r4, r8      
  adcs   r5, r5, r9                 
  adcs   r6, r6, r10
  adcs   r7, r7, r11, lsr #1 
  
  and    r8, r7, #0x7FFFFFFF 
  adds   r4, r4, r7, lsr #31        
  adcs   r5, r5, #0                 
  adcs   r6, r6, #0 
  adcs   r7, r8, #0 
  stmia  r1!, {r4-r7}         // Storing c1  
  
  add    r13, r13, #16        // Restoring stack  
  pop    {r4-r12}
  bx     lr


//***********************************************************************
//  Addition over GF(p^2)
//  Operation: c [reg_p3] = a [reg_p1] + b [reg_p2] mod p, p = 2^127-1
//*********************************************************************** 
.global fp2add1271_a
.type	fp2add1271_a, %function
fp2add1271_a:
  push   {r4-r10}
  ldm    r0!, {r3-r4}
  ldm    r1!, {r7-r8}
  ldm    r0!, {r5-r6}
  ldm    r1!, {r9-r10}
  adds   r3, r3, r7
  adcs   r4, r4, r8
  adcs   r5, r5, r9
  adcs   r6, r6, r10 
  
  ldm    r1!, {r7-r10}
  adds   r3, r3, r6, lsr #31        
  adcs   r4, r4, #0                 
  adcs   r5, r5, #0 
  adcs   r6, r6, #0
  and    r6, r6, #0x7FFFFFFF 

  stmia  r2!, {r3-r4} 
  ldm    r0!, {r3-r4}
  stmia  r2!, {r5-r6} 
  ldm    r0!, {r5-r6}
  adds   r3, r3, r7
  adcs   r4, r4, r8
  adcs   r5, r5, r9
  adcs   r6, r6, r10 
  
  and    r7, r6, #0x7FFFFFFF 
  adds   r3, r3, r6, lsr #31        
  adcs   r4, r4, #0                 
  adcs   r5, r5, #0 
  adcs   r6, r7, #0
  stmia  r2!, {r3-r6}
  pop    {r4-r10}
  bx     lr
  

//***********************************************************************
//  Subtraction over GF(p^2)
//  Operation: c [reg_p3] = a [reg_p1] - b [reg_p2] mod p, p = 2^127-1
//*********************************************************************** 
.global fp2sub1271_a
.type	fp2sub1271_a, %function
fp2sub1271_a:
  push   {r4-r10}
  ldm    r0!, {r3-r4}
  ldm    r1!, {r7-r8}
  ldm    r0!, {r5-r6}
  ldm    r1!, {r9-r10}
  subs   r3, r3, r7
  sbcs   r4, r4, r8
  sbcs   r5, r5, r9
  sbcs   r6, r6, r10 
  
  ldm    r1!, {r7-r10}
  subs   r3, r3, r6, lsr #31
  sbcs   r4, r4, #0 
  sbcs   r5, r5, #0 
  sbcs   r6, r6, #0
  and    r6, r6, #0x7FFFFFFF 

  stmia  r2!, {r3-r4} 
  ldm    r0!, {r3-r4}
  stmia  r2!, {r5-r6} 
  ldm    r0!, {r5-r6}  
  subs   r3, r3, r7
  sbcs   r4, r4, r8
  sbcs   r5, r5, r9
  sbcs   r6, r6, r10 
  
  and    r7, r6, #0x7FFFFFFF 
  subs   r3, r3, r6, lsr #31
  sbcs   r4, r4, #0 
  sbcs   r5, r5, #0 
  sbcs   r6, r7, #0
  stmia  r2!, {r3-r6}
  pop    {r4-r10}
  bx     lr
    

//***********************************************************************
//  Integer multiplication
//  Operation: c [reg_p3] = a [reg_p1] * b [reg_p2]
//  Restriction : b != c
//*********************************************************************** 
.global mul1271_a
.type	mul1271_a, %function
mul1271_a:
  push   {r4-r12}
  ldm    r0!, {r3-r6}
  ldr    r7, [r1]
  mov    r10, #0
  mov    r11, #0
  mov    r12, #0
  umull  r9, r8, r7, r3      
  umlal  r8, r10, r7, r4
  umlal  r10, r11, r7, r5
  umlal  r11, r12, r7, r6
  ldr    r7, [r1, #4]
  str    r9, [r2], #4
  mov    r9, #0
  
  umlal  r8, r9, r7, r3      
  umaal  r9, r10, r7, r4 
  umaal  r10, r11, r7, r5 
  umaal  r11, r12, r7, r6 
  ldr    r7, [r1, #8]
  str    r8, [r2], #4
  mov    r8, #0
  
  umlal  r9, r8, r7, r3      
  umaal  r8, r10, r7, r4 
  umaal  r10, r11, r7, r5 
  umaal  r11, r12, r7, r6 
  ldr    r7, [r1, #12] 
  str    r9, [r2], #4
  mov    r9, #0
  
  umlal  r8, r9, r7, r3      
  umaal  r9, r10, r7, r4 
  umaal  r10, r11, r7, r5 
  umaal  r11, r12, r7, r6 
  stmia  r2!, {r8-r12}
  pop    {r4-r12}
  bx     lr
  

//***********************************************************************
//  Modular reduction
//  Operation: c [reg_p3] = a [reg_p1] mod p, p = 2^127-1
//*********************************************************************** 
.global mod1271_a
.type	mod1271_a, %function
mod1271_a:
  push   {r4-r9}
  ldm    r0!, {r2-r9}
  lsl    r9, r9, #1
  orr    r9, r9, r8, lsr #31
  lsl    r8, r8, #1
  orr    r8, r8, r7, lsr #31
  lsl    r7, r7, #1
  orr    r7, r7, r6, lsr #31
  lsl    r6, r6, #1
  orr    r6, r6, r5, lsr #31
    
  and    r5, r5, #0x7FFFFFFF 
  adds   r2, r2, r6
  adcs   r3, r3, r7                 
  adcs   r4, r4, r8 
  adcs   r5, r5, r9
  
  and    r6, r5, #0x7FFFFFFF 
  adds   r2, r2, r5, lsr #31        
  adcs   r3, r3, #0                 
  adcs   r4, r4, #0 
  adcs   r5, r6, #0   
  stmia  r1!, {r2-r5}
  pop    {r4-r9}
  bx     lr


//***********************************************************************
//  Field addition
//  Operation: c [reg_p3] = a [reg_p1] + b [reg_p2] mod p, p = 2^127-1
//*********************************************************************** 
.global fpadd1271_a
.type	fpadd1271_a, %function
fpadd1271_a:
  push   {r4-r10}
  ldm    r0!, {r3-r4}
  ldm    r1!, {r7-r8}
  ldm    r0!, {r5-r6}
  ldm    r1!, {r9-r10}
  adds   r3, r3, r7
  adcs   r4, r4, r8
  adcs   r5, r5, r9
  adcs   r6, r6, r10 
  
  and    r7, r6, #0x7FFFFFFF 
  adds   r3, r3, r6, lsr #31        
  adcs   r4, r4, #0                 
  adcs   r5, r5, #0 
  adcs   r6, r7, #0
  stmia  r2!, {r3-r6}
  pop    {r4-r10}
  bx     lr
  

//***********************************************************************
//  Field subtraction
//  Operation: c [reg_p3] = a [reg_p1] - b [reg_p2] mod p, p = 2^127-1
//*********************************************************************** 
.global fpsub1271_a
.type	fpsub1271_a, %function
fpsub1271_a:
  push   {r4-r10}
  ldm    r0!, {r3-r4}
  ldm    r1!, {r7-r8}
  ldm    r0!, {r5-r6}
  ldm    r1!, {r9-r10}
  subs   r3, r3, r7
  sbcs   r4, r4, r8
  sbcs   r5, r5, r9
  sbcs   r6, r6, r10 
  
  and    r7, r6, #0x7FFFFFFF 
  subs   r3, r3, r6, lsr #31
  sbcs   r4, r4, #0 
  sbcs   r5, r5, #0 
  sbcs   r6, r7, #0
  stmia  r2!, {r3-r6}
  pop    {r4-r10}
  bx     lr