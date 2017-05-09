/***********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: API header file
*
* This code is based on the paper "FourQ: four-dimensional decompositions on a 
* Q-curve over the Mersenne prime" by Craig Costello and Patrick Longa, in Advances 
* in Cryptology - ASIACRYPT, 2015.
* Preprint available at http://eprint.iacr.org/2015/565.
************************************************************************************/  

#ifndef __FOURQ_API_H__
#define __FOURQ_API_H__


// For C++
#ifdef __cplusplus
extern "C" {
#endif


#include "FourQ.h"


/**************** Public ECC API ****************/

// Set generator G = (x,y)
void eccset(point_t G);

// Variable-base scalar multiplication Q = k*P using a 4-dimensional decomposition and protected against side-channel attacks
bool ecc_mul_SCA_secure(point_t P, point_t R, digit_t* k, point_t Q, bool clear_cofactor);

// Double scalar multiplication R = k*G + l*Q, where G is the generator
bool ecc_mul_double(digit_t* k, point_t Q, digit_t* l, point_t R);


/**************** Public API for SchnorrQ ****************/

// SchnorrQ public key generation
// It produces a blinding point BlindingPoint and a public key PublicKey, which is the encoding of P = s*G, where G is the generator and
// s is the output of hashing SecretKey and taking the least significant 32 bytes of the result.
// Input:  32-byte SecretKey
// Output: 32-byte PublicKey and 64-byte BlindingPoint
ECCRYPTO_STATUS SchnorrQ_KeyGeneration_SCA_secure(const unsigned char* SecretKey, unsigned char* PublicKey, unsigned char* BlindingPoint);

// SchnorrQ keypair generation
// It produces a blinding point BlindingPoint, a private key SecretKey and computes the public key PublicKey, which is the encoding of P = s*G, 
// where G is the generator and s is the output of hashing SecretKey and taking the least significant 32 bytes of the result.
// Outputs: 32-byte SecretKey, 32-byte PublicKey and 64-byte BlindingPoint
ECCRYPTO_STATUS SchnorrQ_FullKeyGeneration_SCA_secure(unsigned char* SecretKey, unsigned char* PublicKey, unsigned char* BlindingPoint);

// SchnorrQ signature generation
// It produces the signature Signature of a message Message of size SizeMessage in bytes
// Inputs: 32-byte SecretKey, 32-byte PublicKey, Message of size SizeMessage in bytes, and 64-byte BlindingPoint
// Output: 64-byte Signature 
ECCRYPTO_STATUS SchnorrQ_Sign_SCA_secure(const unsigned char* SecretKey, const unsigned char* PublicKey, const unsigned char* Message, const unsigned int SizeMessage, unsigned char* Signature, unsigned char* BlindingPoint);

// SchnorrQ signature verification
// It verifies the signature Signature of a message Message of size SizeMessage in bytes
// Inputs: 32-byte PublicKey, 64-byte Signature, and Message of size SizeMessage in bytes
// Output: true (valid signature) or false (invalid signature)
ECCRYPTO_STATUS SchnorrQ_Verify(const unsigned char* PublicKey, const unsigned char* Message, const unsigned int SizeMessage, const unsigned char* Signature, unsigned int* valid);


/**************** Public API for co-factor ECDH key exchange with compressed, 32-byte public keys ****************/

// Compressed public key generation for key exchange
// It produces a public key PublicKey, which is the encoding of P = SecretKey*G (G is the generator), and a blinding point BlindingPoint.
// Input:  32-byte SecretKey
// Output: 32-byte PublicKey and 64-byte BlindingPoint
ECCRYPTO_STATUS CompressedPublicKeyGeneration_SCA_secure(const unsigned char* SecretKey, unsigned char* PublicKey, unsigned char* BlindingPoint);

// Keypair generation for key exchange. Public key is compressed to 32 bytes
// It produces a private key SecretKey, a public key PublicKey, which is the encoding of P = SecretKey*G (G is the generator), and a blinding point BlindingPoint.
// Outputs: 32-byte SecretKey, 32-byte PublicKey and 64-byte BlindingPoint
ECCRYPTO_STATUS CompressedKeyGeneration_SCA_secure(unsigned char* SecretKey, unsigned char* PublicKey, unsigned char* BlindingPoint);

// Secret agreement computation for key exchange using a compressed, 32-byte public key
// The output is the y-coordinate of SecretKey*A, where A is the decoding of the public key PublicKey. 
// Inputs: 32-byte SecretKey, 32-byte PublicKey and 64-byte BlindingPoint
// Output: 32-byte SharedSecret
ECCRYPTO_STATUS CompressedSecretAgreement_SCA_secure(const unsigned char* SecretKey, const unsigned char* PublicKey, unsigned char* SharedSecret, unsigned char* BlindingPoint);


/**************** Public API for co-factor ECDH key exchange with uncompressed, 64-byte public keys ****************/

// Public key generation for key exchange
// It produces the public key PublicKey = SecretKey*G, where G is the generator, and a blinding point BlindingPoint.
// Input:  32-byte SecretKey
// Output: 64-byte PublicKey and 64-byte BlindingPoint
ECCRYPTO_STATUS PublicKeyGeneration_SCA_secure(const unsigned char* SecretKey, unsigned char* PublicKey, unsigned char* BlindingPoint);

// Keypair generation for key exchange
// It produces a private key SecretKey, the public key PublicKey = SecretKey*G, where G is the generator, and a blinding point BlindingPoint.
// Outputs: 32-byte SecretKey, 64-byte PublicKey and 64-byte BlindingPoint
ECCRYPTO_STATUS KeyGeneration_SCA_secure(unsigned char* SecretKey, unsigned char* PublicKey, unsigned char* BlindingPoint);

// Secret agreement computation for key exchange
// The output is the y-coordinate of SecretKey*PublicKey. 
// Inputs: 32-byte SecretKey, 64-byte PublicKey and 64-byte BlindingPoint
// Output: 32-byte SharedSecret
ECCRYPTO_STATUS SecretAgreement_SCA_secure(const unsigned char* SecretKey, const unsigned char* PublicKey, unsigned char* SharedSecret, unsigned char* BlindingPoint);


#ifdef __cplusplus
}
#endif


#endif
