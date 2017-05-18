/**********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: digital signature SchnorrQ
*
* See "SchnorrQ: Schnorr signatures on FourQ" by Craig Costello and Patrick Longa,
* MSR Technical Report, 2016. Available at: 
* https://www.microsoft.com/en-us/research/wp-content/uploads/2016/07/SchnorrQ.pdf.
***********************************************************************************/ 

#include "FourQ_internal.h"
#include "FourQ_params.h"
#include "../random/random.h"
#include "../sha512/sha512.h"
#include <malloc.h>
#include <string.h>


ECCRYPTO_STATUS SchnorrQ_KeyGeneration_SCA_secure(const unsigned char* SecretKey, unsigned char* PublicKey, unsigned char* BlindingPoint)
{ // SchnorrQ public key generation
  // It produces a blinding point BlindingPoint and a public key PublicKey, which is the encoding of P = s*G, where G is the generator and
  // s is the output of hashing SecretKey and taking the least significant 32 bytes of the result.
  // Input:  32-byte SecretKey
  // Output: 32-byte PublicKey and 64-byte BlindingPoint
    point_t G, R;
    point_extedwards_t S;
    unsigned char k[64], SecretBlinding[32];
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

    if (CryptoHashFunction(SecretKey, 32, k) != 0) {   
        Status = ECCRYPTO_ERROR;
        goto cleanup;
    }

    Status = ecc_mul_SCA_secure(G, (point_affine*)BlindingPoint, (digit_t*)k, R, false);  // Compute public key
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
    clear_words((unsigned int*)k, 512/(sizeof(unsigned int)*8));
    clear_words((unsigned int*)PublicKey, 256/(sizeof(unsigned int)*8));

    return Status;
}


ECCRYPTO_STATUS SchnorrQ_FullKeyGeneration_SCA_secure(unsigned char* SecretKey, unsigned char* PublicKey, unsigned char* BlindingPoint)
{ // SchnorrQ keypair generation
  // It produces a blinding point BlindingPoint, a private key SecretKey and computes the public key PublicKey, which is the encoding of P = s*G, 
  // where G is the generator and s is the output of hashing SecretKey and taking the least significant 32 bytes of the result.
  // Outputs: 32-byte SecretKey, 32-byte PublicKey and 64-byte BlindingPoint
    ECCRYPTO_STATUS Status = ECCRYPTO_ERROR_UNKNOWN;

    Status = RandomBytesFunction(SecretKey, 32);
    if (Status != ECCRYPTO_SUCCESS) {
        goto cleanup;
    }
  
    Status = SchnorrQ_KeyGeneration_SCA_secure(SecretKey, PublicKey, BlindingPoint);   
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


ECCRYPTO_STATUS SchnorrQ_Sign_SCA_secure(const unsigned char* SecretKey, const unsigned char* PublicKey, const unsigned char* Message, const unsigned int SizeMessage, unsigned char* Signature, unsigned char* BlindingPoint)
{ // SchnorrQ signature generation
  // It produces the signature Signature of a message Message of size SizeMessage in bytes
  // Inputs: 32-byte SecretKey, 32-byte PublicKey, Message of size SizeMessage in bytes, and 64-byte BlindingPoint
  // Output: 64-byte Signature and updated BlindingPoint 
    point_t G, R;
    unsigned char k[64], r[64], h[64], *temp = NULL;
    digit_t* H = (digit_t*)h;
    digit_t* S1 = (digit_t*)(Signature+32);
    digit_t S2[256/(sizeof(digit_t)*8)];
    ECCRYPTO_STATUS Status = ECCRYPTO_ERROR_UNKNOWN;
      
    if (CryptoHashFunction(SecretKey, 32, k) != 0) {   
        Status = ECCRYPTO_ERROR;
        goto cleanup;
    }
    
    temp = (unsigned char*)calloc(1, SizeMessage+64);
    if (temp == NULL) {
        Status = ECCRYPTO_ERROR_NO_MEMORY;
        goto cleanup;
    }
    
    memmove(temp+32, k+32, 32);
    memmove(temp+64, Message, SizeMessage);
  
    if (CryptoHashFunction(temp+32, SizeMessage+32, r) != 0) {   
        Status = ECCRYPTO_ERROR;
        goto cleanup;
    }
    
    fp2copy1271((felm_t*)&GENERATOR_x[0], G->x);
    fp2copy1271((felm_t*)&GENERATOR_y[0], G->y);   
    
    Status = ecc_mul_SCA_secure(G, (point_affine*)BlindingPoint, (digit_t*)r, R, false);  // Also verifies that BlindingPoint is a point on the curve. If not, it fails
    if (Status != ECCRYPTO_SUCCESS) {
        goto cleanup;
    }
    encode(R, Signature);                   // Encode lowest 32 bytes of signature
    memmove(temp, Signature, 32);
    memmove(temp+32, PublicKey, 32);
  
    if (CryptoHashFunction(temp, SizeMessage+64, h) != 0) {   
        Status = ECCRYPTO_ERROR;
        goto cleanup;
    }    

    Status = RandomBytesFunction((unsigned char*)S2, 32);
    if (Status != ECCRYPTO_SUCCESS) {
        goto cleanup;
    }
    modulo_order(S2, S2);
    subtract_mod_order((digit_t*)k, S2, S1);
    modulo_order((digit_t*)r, (digit_t*)r);
    modulo_order(H, H);
    to_Montgomery(S1, S1);                    // Converting to Montgomery representation
    to_Montgomery(S2, S2);
    to_Montgomery(H, H);                      // Converting to Montgomery representation
    Montgomery_multiply_mod_order(S1, H, S1);
    Montgomery_multiply_mod_order(S2, H, S2);
    from_Montgomery(S1, S1);                  // Converting back to standard representation
    from_Montgomery(S2, S2);                 
    subtract_mod_order((digit_t*)r, S1, S1);          
    subtract_mod_order(S1, S2, S1);
    Status = ECCRYPTO_SUCCESS;
    
cleanup:
    if (temp != NULL)
        free(temp);
    clear_words((unsigned int*)k, 512/(sizeof(unsigned int)*8));
    clear_words((unsigned int*)r, 512/(sizeof(unsigned int)*8));
    clear_words((unsigned int*)S2, 256/(sizeof(unsigned int)*8));
    
    return Status;
}


ECCRYPTO_STATUS SchnorrQ_Verify(const unsigned char* PublicKey, const unsigned char* Message, const unsigned int SizeMessage, const unsigned char* Signature, unsigned int* valid)
{ // SchnorrQ signature verification
  // It verifies the signature Signature of a message Message of size SizeMessage in bytes
  // Inputs: 32-byte PublicKey, 64-byte Signature, and Message of size SizeMessage in bytes
  // Output: true (valid signature) or false (invalid signature)
    point_t A;
    unsigned char *temp, h[64];
    unsigned int i;
    ECCRYPTO_STATUS Status = ECCRYPTO_ERROR_UNKNOWN;  

    *valid = false;

    temp = (unsigned char*)calloc(1, SizeMessage+64);
    if (temp == NULL) {
        Status = ECCRYPTO_ERROR_NO_MEMORY;
        goto cleanup;
    }

    if (((PublicKey[15] & 0x80) != 0) || ((Signature[15] & 0x80) != 0) || (Signature[63] != 0) || ((Signature[62] & 0xC0) != 0)) {  // Are bit128(PublicKey) = bit128(Signature) = 0 and Signature+32 < 2^246?
        Status = ECCRYPTO_ERROR_INVALID_PARAMETER;
        goto cleanup;
    }
    
    Status = decode(PublicKey, A);    // Also verifies that A is on the curve. If it is not, it fails  
    if (Status != ECCRYPTO_SUCCESS) {
        goto cleanup;                            
    }

    memmove(temp, Signature, 32);
    memmove(temp+32, PublicKey, 32);
    memmove(temp+64, Message, SizeMessage);
  
    if (CryptoHashFunction(temp, SizeMessage+64, h) != 0) {   
        Status = ECCRYPTO_ERROR;
        goto cleanup;
    }

    Status = ecc_mul_double((digit_t*)(Signature+32), A, (digit_t*)h, A);      
    if (Status != ECCRYPTO_SUCCESS) {                                                
        goto cleanup;
    }
    
    encode(A, (unsigned char*)A);

    for (i = 0; i < NWORDS_ORDER; i++) {
        if (((digit_t*)A)[i] != ((digit_t*)Signature)[i]) {
            goto cleanup;   
        }
    }
    *valid = true;

cleanup:
    if (temp != NULL)
        free(temp);
    
    return Status;
}
