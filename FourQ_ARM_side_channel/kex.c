/********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: Diffie-Hellman key exchange based on FourQ, including countermeasures
*           against side-channel attacks
*           option 1: co-factor ecdh using compressed 32-byte public keys,
*           (see https://datatracker.ietf.org/doc/draft-ladd-cfrg-4q/).         
*           option 2: co-factor ecdh using uncompressed, 64-byte public keys.         
*********************************************************************************/

#include "FourQ_internal.h"
#include "FourQ_params.h"
#include "../random/random.h"
#include <string.h>


static __inline bool is_neutral_point(point_t P)
{ // Is P the neutral point (0,1)?
  // SECURITY NOTE: this function does not run in constant time (input point P is assumed to be public).

    if (is_zero_ct((digit_t*)P->x, 2*NWORDS_FIELD) && is_zero_ct(&((digit_t*)P->y)[1], 2*NWORDS_FIELD-1) && is_digit_zero_ct(P->y[0][0] - 1)) {  
		return true;
    }
    return false;
}


/*************** ECDH USING COMPRESSED, 32-BYTE PUBLIC KEYS ***************/

ECCRYPTO_STATUS CompressedPublicKeyGeneration_SCA_secure(const unsigned char* SecretKey, unsigned char* PublicKey, unsigned char* BlindingPoint)
{ // Compressed public key generation for key exchange
  // It produces a public key PublicKey, which is the encoding of P = SecretKey*G (G is the generator), and a blinding point BlindingPoint.
  // Input:  32-byte SecretKey
  // Output: 32-byte PublicKey and 64-byte BlindingPoint
    point_t G, R; 
    point_extedwards_t S;
    unsigned char SecretBlinding[32];
	ECCRYPTO_STATUS Status = ECCRYPTO_ERROR_UNKNOWN;

	Status = RandomBytesFunction(SecretBlinding, 32);
	if (Status != ECCRYPTO_SUCCESS) {
		goto cleanup;
	}
        
    // Set up an initial "weak" blinding point R
    fp2copy1271((felm_t*)&GENERATOR_x[0], G->x);
    fp2copy1271((felm_t*)&GENERATOR_y[0], G->y);   
    point_setup(G, S);                                                             
    eccdouble(S);       
    eccnorm(S, R);                 
    
    // Computing an initial blinding point. This computation itself is not protected with a secure point blinding
	Status = ecc_mul_SCA_secure(G, R, (digit_t*)SecretBlinding, (point_affine*)BlindingPoint, false);
	if (Status != ECCRYPTO_SUCCESS) {
		goto cleanup;
	}

	Status = ecc_mul_SCA_secure(G, (point_affine*)BlindingPoint, (digit_t*)SecretKey, R, false);  // Compute public key
	if (Status != ECCRYPTO_SUCCESS) {
		goto cleanup;
	}                             
	encode(R, PublicKey);                   // Encode public key
    
// Cleanup
	clear_words((unsigned int*)SecretBlinding, 256/(sizeof(unsigned int)*8));
    return ECCRYPTO_SUCCESS;

cleanup:
	clear_words((unsigned int*)SecretBlinding, 256/(sizeof(unsigned int)*8));
	clear_words((unsigned int*)BlindingPoint, 512/(sizeof(unsigned int)*8));
	clear_words((unsigned int*)PublicKey, 256/(sizeof(unsigned int)*8));

	return Status;
}


ECCRYPTO_STATUS CompressedKeyGeneration_SCA_secure(unsigned char* SecretKey, unsigned char* PublicKey, unsigned char* BlindingPoint)
{ // Keypair generation for key exchange. Public key is compressed to 32 bytes
  // It produces a private key SecretKey, a public key PublicKey, which is the encoding of P = SecretKey*G (G is the generator), and a blinding point BlindingPoint.
  // Outputs: 32-byte SecretKey, 32-byte PublicKey and 64-byte BlindingPoint
    ECCRYPTO_STATUS Status = ECCRYPTO_ERROR_UNKNOWN;

	Status = RandomBytesFunction(SecretKey, 32);
	if (Status != ECCRYPTO_SUCCESS) {
		goto cleanup;
	}
  
    Status = CompressedPublicKeyGeneration_SCA_secure(SecretKey, PublicKey, BlindingPoint);
    if (Status != ECCRYPTO_SUCCESS) {
        goto cleanup;
    }

    return ECCRYPTO_SUCCESS;

cleanup:
    clear_words((unsigned int*)SecretKey, 256/(sizeof(unsigned int)*8));
    clear_words((unsigned int*)PublicKey, 256/(sizeof(unsigned int)*8));
	clear_words((unsigned int*)BlindingPoint, 512/(sizeof(unsigned int)*8));

    return Status;
}


ECCRYPTO_STATUS CompressedSecretAgreement_SCA_secure(const unsigned char* SecretKey, const unsigned char* PublicKey, unsigned char* SharedSecret, unsigned char* BlindingPoint)
{ // Secret agreement computation for key exchange using a compressed, 32-byte public key
  // The output is the y-coordinate of SecretKey*A, where A is the decoding of the public key PublicKey. 
  // Inputs: 32-byte SecretKey, 32-byte PublicKey and 64-byte BlindingPoint
  // Output: 32-byte SharedSecret and updated BlindingPoint
    point_t A;
    ECCRYPTO_STATUS Status = ECCRYPTO_ERROR_UNKNOWN;

    if ((PublicKey[15] & 0x80) != 0) {  // Is bit128(PublicKey) = 0?
		Status = ECCRYPTO_ERROR_INVALID_PARAMETER;
		goto cleanup;
    }

	Status = decode(PublicKey, A);    // Also verifies that A is on the curve. If it is not, it fails
	if (Status != ECCRYPTO_SUCCESS) {
		goto cleanup;
	}
         
    Status = ecc_mul_SCA_secure(A, (point_affine*)BlindingPoint, (digit_t*)SecretKey, A, true);
	if (Status != ECCRYPTO_SUCCESS) {
		goto cleanup;
	}

    if (is_neutral_point(A)) {  // Is output = neutral point (0,1)?
		Status = ECCRYPTO_ERROR_SHARED_KEY;
		goto cleanup;
    }
  
	memmove(SharedSecret, (unsigned char*)A->y, 32);

	return ECCRYPTO_SUCCESS;
    
cleanup:
    clear_words((unsigned int*)SharedSecret, 256/(sizeof(unsigned int)*8));
    
    return Status;
}


/*************** ECDH USING UNCOMPRESSED PUBLIC KEYS ***************/

ECCRYPTO_STATUS PublicKeyGeneration_SCA_secure(const unsigned char* SecretKey, unsigned char* PublicKey, unsigned char* BlindingPoint)
{ // Public key generation for key exchange
  // It produces the public key PublicKey = SecretKey*G, where G is the generator, and a blinding point BlindingPoint.
  // Input:  32-byte SecretKey
  // Output: 64-byte PublicKey and 64-byte BlindingPoint
    point_t G, R; 
    point_extedwards_t S;
    unsigned char SecretBlinding[32];
	ECCRYPTO_STATUS Status = ECCRYPTO_ERROR_UNKNOWN;

	Status = RandomBytesFunction(SecretBlinding, 32);
	if (Status != ECCRYPTO_SUCCESS) {
		goto cleanup;
	}
        
    // Set up an initial "weak" blinding point R
    fp2copy1271((felm_t*)&GENERATOR_x[0], G->x);
    fp2copy1271((felm_t*)&GENERATOR_y[0], G->y);   
    point_setup(G, S);                                                             
    eccdouble(S);       
    eccnorm(S, R);                 
    
    // Computing an initial blinding point. This computation itself is not protected with a secure point blinding
	Status = ecc_mul_SCA_secure(G, R, (digit_t*)SecretBlinding, (point_affine*)BlindingPoint, false);
	if (Status != ECCRYPTO_SUCCESS) {
		goto cleanup;
	}

	Status = ecc_mul_SCA_secure(G, (point_affine*)BlindingPoint, (digit_t*)SecretKey, (point_affine*)PublicKey, false);  // Compute public key
	if (Status != ECCRYPTO_SUCCESS) {
		goto cleanup;
	}
    
// Cleanup
	clear_words((unsigned int*)SecretBlinding, 256/(sizeof(unsigned int)*8));
    return ECCRYPTO_SUCCESS;

cleanup:
	clear_words((unsigned int*)SecretBlinding, 256/(sizeof(unsigned int)*8));
	clear_words((unsigned int*)BlindingPoint, 512/(sizeof(unsigned int)*8));
	clear_words((unsigned int*)PublicKey, 512/(sizeof(unsigned int)*8));

	return Status;
}


ECCRYPTO_STATUS KeyGeneration_SCA_secure(unsigned char* SecretKey, unsigned char* PublicKey, unsigned char* BlindingPoint)
{ // Keypair generation for key exchange
  // It produces a private key SecretKey, the public key PublicKey = SecretKey*G, where G is the generator, and a blinding point BlindingPoint.
  // Outputs: 32-byte SecretKey, 64-byte PublicKey and 64-byte BlindingPoint
	ECCRYPTO_STATUS Status = ECCRYPTO_ERROR_UNKNOWN;

	Status = RandomBytesFunction(SecretKey, 32);
	if (Status != ECCRYPTO_SUCCESS) {
		goto cleanup;
	}

	Status = PublicKeyGeneration_SCA_secure(SecretKey, PublicKey, BlindingPoint);
	if (Status != ECCRYPTO_SUCCESS) {
		goto cleanup;
	}

	return ECCRYPTO_SUCCESS;

cleanup:
	clear_words((unsigned int*)SecretKey, 256/(sizeof(unsigned int)*8));
	clear_words((unsigned int*)PublicKey, 512/(sizeof(unsigned int)*8));
	clear_words((unsigned int*)BlindingPoint, 512/(sizeof(unsigned int)*8));

	return Status;
}


ECCRYPTO_STATUS SecretAgreement_SCA_secure(const unsigned char* SecretKey, const unsigned char* PublicKey, unsigned char* SharedSecret, unsigned char* BlindingPoint)
{ // Secret agreement computation for key exchange
  // The output is the y-coordinate of SecretKey*PublicKey.  
  // Inputs: 32-byte SecretKey, 64-byte PublicKey and 64-byte BlindingPoint
  // Output: 32-byte SharedSecret and updated BlindingPoint
	point_t A;
	ECCRYPTO_STATUS Status = ECCRYPTO_ERROR_UNKNOWN;

    if (((PublicKey[15] & 0x80) != 0) || ((PublicKey[31] & 0x80) != 0) || ((PublicKey[47] & 0x80) != 0) || ((PublicKey[63] & 0x80) != 0)) {  // Are PublicKey_x[i] and PublicKey_y[i] < 2^127?
		Status = ECCRYPTO_ERROR_INVALID_PARAMETER;
		goto cleanup;
    }

	Status = ecc_mul_SCA_secure((point_affine*)PublicKey, (point_affine*)BlindingPoint, (digit_t*)SecretKey, A, true);  // Also verifies that PublicKey and BlindingPoint are points on the curve. If not, it fails
	if (Status != ECCRYPTO_SUCCESS) {
		goto cleanup;
	}

    if (is_neutral_point(A)) {  // Is output = neutral point (0,1)?
		Status = ECCRYPTO_ERROR_SHARED_KEY;
		goto cleanup;
    }
  
	memmove(SharedSecret, (unsigned char*)A->y, 32);

	return ECCRYPTO_SUCCESS;

cleanup:
	clear_words((unsigned int*)SharedSecret, 256/(sizeof(unsigned int)*8));

	return Status;
}