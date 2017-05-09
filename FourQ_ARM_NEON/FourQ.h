/*****************************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: main header file
*
* This code is based on the papers:
* [1] "FourQ: four-dimensional decompositions on a Q-curve over the Mersenne prime" 
*     by Craig Costello and Patrick Longa, ASIACRYPT2015 (http://eprint.iacr.org/2015/565).
* [2] "FourQNEON: Faster Elliptic Curve Scalar Multiplications on ARM Processors" 
*     by Patrick Longa, SAC2016 (http://eprint.iacr.org/2016/645).
******************************************************************************************/  

#ifndef __FOURQ_H__
#define __FOURQ_H__


// For C++
#ifdef __cplusplus
extern "C" {
#endif


#include <stdint.h>
#include <stdbool.h>


// Definition of operating system

#define OS_LINUX     1

#if defined(__LINUX__)                    // Linux OS
    #define OS_TARGET OS_LINUX 
#else
    #error -- "Unsupported OS"
#endif


// Definition of compiler

#define COMPILER_GCC     1
#define COMPILER_CLANG   2

#if defined(__GNUC__)           // GNU GCC compiler
    #define COMPILER COMPILER_GCC   
#elif defined(__clang__)        // Clang compiler
    #define COMPILER COMPILER_CLANG   
#else
    #error -- "Unsupported COMPILER"
#endif


// Definition of the targeted architecture and basic data types
    
#define TARGET_ARM          1

#if defined(_ARM_)
    #define TARGET TARGET_ARM
    #define RADIX           32
    typedef uint32_t        digit_t;      // Unsigned 32-bit digit
    typedef int32_t         sdigit_t;     // Signed 32-bit digit
    #define NWORDS_FIELD    4             
    #define NWORDS_ORDER    8 
#else
    #error -- "Unsupported ARCHITECTURE"
#endif


// Constants

#define RADIX64         64
#define NWORDS64_FIELD  2                 // Number of 64-bit words of a field element 
#define NWORDS64_ORDER  4                 // Number of 64-bit words of an element in Z_r 


// Definition of complementary cryptographic functions

#define RandomBytesFunction     random_bytes    
#define CryptoHashFunction      crypto_sha512        // Use SHA-512 by default


// Detect if additional optimizationes are enabled

#if defined(_INTERLEAVE_)
    #define INTERLEAVE                    // Interleaving of instructions 
#endif

#if defined(_MIX_ARM_NEON_)
    #define MIX_ARM_NEON                  // Mix ARM/NEON instructions 
#endif


// Basic parameters for variable-base scalar multiplication (without using endomorphisms)
#define W_VARBASE             5 
#define NBITS_ORDER_PLUS_ONE  246+1


// Basic parameters for fixed-base scalar multiplication
#define W_FIXEDBASE       5              // Memory requirement: 7.5KB (storage for 80 points).
#define V_FIXEDBASE       5              

// Basic parameters for double scalar multiplication
#define WP_DOUBLEBASE     8              // Memory requirement: 24KB (storage for 256 points).
#define WQ_DOUBLEBASE     4  
   

// FourQ's basic element definitions and point representations

typedef digit_t felm_t[NWORDS_FIELD];                      // Datatype for representing 128-bit field elements 
typedef felm_t f2elm_t[2];                                 // Datatype for representing quadratic extension field elements
        
typedef struct { f2elm_t x; f2elm_t y; } point_affine;     // Point representation in affine coordinates.
typedef point_affine point_t[1]; 


// FourQ's vectorized element definitions and point representations 

#define VWORDS_FIELD 5                                     // Number of 32-bit words of a vectorized field element

typedef uint32_t velm_t[VWORDS_FIELD];                     // Datatype for representing 128-bit vectorized field elements 
typedef uint32_t v2elm_t[2*VWORDS_FIELD];                  // Datatype for representing vectorized quadratic extension field elements 

typedef struct { v2elm_t x; v2elm_t y; } vpoint_affine;    // Point representation in affine coordinates.
typedef vpoint_affine vpoint_t[1]; 


// Definitions of the error-handling type and error codes

typedef enum {
	ECCRYPTO_ERROR,                            // 0x00
	ECCRYPTO_SUCCESS,                          // 0x01
	ECCRYPTO_ERROR_DURING_TEST,                // 0x02
	ECCRYPTO_ERROR_UNKNOWN,                    // 0x03
	ECCRYPTO_ERROR_NOT_IMPLEMENTED,            // 0x04
	ECCRYPTO_ERROR_NO_MEMORY,                  // 0x05
	ECCRYPTO_ERROR_INVALID_PARAMETER,          // 0x06
	ECCRYPTO_ERROR_SHARED_KEY,                 // 0x07
	ECCRYPTO_ERROR_SIGNATURE_VERIFICATION,     // 0x08
	ECCRYPTO_ERROR_END_OF_LIST
} ECCRYPTO_STATUS;

#define ECCRYPTO_STATUS_TYPE_SIZE (ECCRYPTO_ERROR_END_OF_LIST)


// Error message definitions

#define ECCRYPTO_MSG_ERROR                                  "ECCRYPTO_ERROR"
#define ECCRYPTO_MSG_SUCCESS                                "ECCRYPTO_SUCCESS"
#define ECCRYPTO_MSG_ERROR_DURING_TEST                      "ECCRYPTO_ERROR_DURING_TEST"
#define ECCRYPTO_MSG_ERROR_UNKNOWN                          "ECCRYPTO_ERROR_UNKNOWN"
#define ECCRYPTO_MSG_ERROR_NOT_IMPLEMENTED                  "ECCRYPTO_ERROR_NOT_IMPLEMENTED"
#define ECCRYPTO_MSG_ERROR_NO_MEMORY                        "ECCRYPTO_ERROR_NO_MEMORY"
#define ECCRYPTO_MSG_ERROR_INVALID_PARAMETER                "ECCRYPTO_ERROR_INVALID_PARAMETER"
#define ECCRYPTO_MSG_ERROR_SHARED_KEY                       "ECCRYPTO_ERROR_SHARED_KEY"
#define ECCRYPTO_MSG_ERROR_SIGNATURE_VERIFICATION           "ECCRYPTO_ERROR_SIGNATURE_VERIFICATION"


#ifdef __cplusplus
}
#endif


#endif
