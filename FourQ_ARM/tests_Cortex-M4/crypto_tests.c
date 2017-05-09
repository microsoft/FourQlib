/***********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: testing code for cryptographic functions based on FourQ 
************************************************************************************/   

#include "../FourQ_internal.h"
#include "../FourQ_params.h"
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


ECCRYPTO_STATUS SchnorrQ_test()
{ // Test the SchnorrQ digital signature scheme
    int n, passed;       
    void *msg = NULL; 
    unsigned int len, valid = false;
    unsigned char SecretKey[32], PublicKey[32], Signature[64];
    ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;

    print_test("\n--------------------------------------------------------------------------------------------------------\n"); 
    print_test("Testing the SchnorrQ signature scheme: \n"); 

    passed = 1;
    for (n = 0; n < TEST_LOOPS; n++)
    {    
        // Signature key generation
        Status = SchnorrQ_FullKeyGeneration(SecretKey, PublicKey);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }  

        // Signature computation
        msg = "a";  
        len = 1;
        Status = SchnorrQ_Sign(SecretKey, PublicKey, msg, len, Signature);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }    

        // Valid signature test
        Status = SchnorrQ_Verify(PublicKey, msg, len, Signature, &valid);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }    
        if (valid == false) {
            passed = 0;
            break;
        }

        // Invalid signature test (flipping one bit of the message)
        msg = "b";  
        Status = SchnorrQ_Verify(PublicKey, msg, len, Signature, &valid);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }    
        if (valid == true) {
            passed = 0;
            break;
        }
    } 
    if (passed==1) print_test("  Signature tests.................................................................. PASSED");
    else { print_test("  Signature tests... FAILED"); print_test("\n"); Status = ECCRYPTO_ERROR_SIGNATURE_VERIFICATION; }
    
    return Status;
}


ECCRYPTO_STATUS SchnorrQ_run()
{ // Benchmark the SchnorrQ digital signature scheme 
    int n;
    unsigned int cycles, cycles1, cycles2;   
    void *msg = NULL;
    unsigned int len = 0, valid = false;
    unsigned char SecretKey[32], PublicKey[32], Signature[64];
    ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;

    print_test("\n--------------------------------------------------------------------------------------------------------\n"); 
    print_test("Benchmarking the SchnorrQ signature scheme: \n"); 
    
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        Status = SchnorrQ_FullKeyGeneration(SecretKey, PublicKey);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }    
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    print_bench("  SchnorrQ's key generation runs in ............................................... ", cycles/BENCH_LOOPS);
    
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        Status = SchnorrQ_Sign(SecretKey, PublicKey, msg, len, Signature);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }    
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    print_bench("  SchnorrQ's signing runs in ...................................................... ", cycles/BENCH_LOOPS);
    
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        Status = SchnorrQ_Verify(PublicKey, msg, len, Signature, &valid);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }    
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    print_bench("  SchnorrQ's verification runs in ................................................. ", cycles/BENCH_LOOPS);
    
    return Status;
}


ECCRYPTO_STATUS compressedkex_test()
{ // Test ECDH key exchange based on FourQ
	int n, passed;
	unsigned int i;
	unsigned char SecretKeyA[32], PublicKeyA[32], SecretAgreementA[32];
	unsigned char SecretKeyB[32], PublicKeyB[32], SecretAgreementB[32];
	ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;

	print_test("\n--------------------------------------------------------------------------------------------------------\n");
	print_test("Testing DH key exchange using compressed, 32-byte public keys: \n");

	passed = 1;
	for (n = 0; n < TEST_LOOPS; n++)
	{
		// Alice's keypair generation
		Status = CompressedKeyGeneration(SecretKeyA, PublicKeyA);
		if (Status != ECCRYPTO_SUCCESS) {
			return Status;
		}
		// Bob's keypair generation
		Status = CompressedKeyGeneration(SecretKeyB, PublicKeyB);
		if (Status != ECCRYPTO_SUCCESS) {
			return Status;
		}

		// Alice's shared secret computation
		Status = CompressedSecretAgreement(SecretKeyA, PublicKeyB, SecretAgreementA);
		if (Status != ECCRYPTO_SUCCESS) {
			return Status;
		}
		// Bob's shared secret computation
		Status = CompressedSecretAgreement(SecretKeyB, PublicKeyA, SecretAgreementB);
		if (Status != ECCRYPTO_SUCCESS) {
			return Status;
		}

		for (i = 0; i < 32; i++) {
		    if (SecretAgreementA[i] != SecretAgreementB[i]) {
				passed = 0;
				break;
			}
		}
	}
	if (passed==1) print_test("  DH key exchange tests............................................................ PASSED");
	else { print_test("  DH key exchange tests... FAILED"); print_test("\n"); Status = ECCRYPTO_ERROR_SHARED_KEY; }

	return Status;
}


ECCRYPTO_STATUS compressedkex_run()
{ // Benchmark ECDH key exchange based on FourQ 
	int n;
	unsigned int cycles, cycles1, cycles2;
	unsigned char SecretKeyA[32], PublicKeyA[32], SecretAgreementA[32];
	unsigned char SecretKeyB[32], PublicKeyB[32];
	ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;

	print_test("\n--------------------------------------------------------------------------------------------------------\n");
	print_test("Benchmarking DH key exchange using compressed, 32-byte public keys: \n");

	cycles = 0;
	for (n = 0; n < BENCH_LOOPS; n++)
	{
		cycles1 = cpucycles();
		Status = CompressedKeyGeneration(SecretKeyA, PublicKeyA);
		if (Status != ECCRYPTO_SUCCESS) {
			return Status;
		}
		cycles2 = cpucycles();
		cycles = cycles + (cycles2 - cycles1);
	}
	print_bench("  Keypair generation runs in ...................................................... ", cycles/BENCH_LOOPS);
	
	Status = CompressedKeyGeneration(SecretKeyB, PublicKeyB);
	cycles = 0;
	for (n = 0; n < BENCH_LOOPS; n++)
	{
		cycles1 = cpucycles();
		Status = CompressedSecretAgreement(SecretKeyA, PublicKeyB, SecretAgreementA);
		if (Status != ECCRYPTO_SUCCESS) {
			return Status;
		}
		cycles2 = cpucycles();
		cycles = cycles + (cycles2 - cycles1);
	}
	print_bench("  Secret agreement runs in ........................................................ ", cycles/BENCH_LOOPS);

	return Status;
}


ECCRYPTO_STATUS kex_test()
{ // Test ECDH key exchange based on FourQ
	int n, passed;
	unsigned int i;
	unsigned char SecretKeyA[32], PublicKeyA[64], SecretAgreementA[32];
	unsigned char SecretKeyB[32], PublicKeyB[64], SecretAgreementB[32];
	ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;

	print_test("\n--------------------------------------------------------------------------------------------------------\n");
	print_test("Testing DH key exchange using uncompressed, 64-byte public keys: \n");

	passed = 1;
	for (n = 0; n < TEST_LOOPS; n++)
	{
		// Alice's keypair generation
		Status = KeyGeneration(SecretKeyA, PublicKeyA);
		if (Status != ECCRYPTO_SUCCESS) {
			return Status;
		}
		// Bob's keypair generation
		Status = KeyGeneration(SecretKeyB, PublicKeyB);
		if (Status != ECCRYPTO_SUCCESS) {
			return Status;
		}

		// Alice's shared secret computation
		Status = SecretAgreement(SecretKeyA, PublicKeyB, SecretAgreementA);
		if (Status != ECCRYPTO_SUCCESS) {
			return Status;
		}
		// Bob's shared secret computation
		Status = SecretAgreement(SecretKeyB, PublicKeyA, SecretAgreementB);
		if (Status != ECCRYPTO_SUCCESS) {
			return Status;
		}

		for (i = 0; i < 32; i++) {
			if (SecretAgreementA[i] != SecretAgreementB[i]) {
				passed = 0;
				break;
			}
		}
	}
	if (passed==1) print_test("  DH key exchange tests............................................................ PASSED");
	else { print_test("  DH key exchange tests... FAILED"); print_test("\n"); Status = ECCRYPTO_ERROR_SHARED_KEY; }

	return Status;
}


ECCRYPTO_STATUS kex_run()
{ // Benchmark ECDH key exchange based on FourQ 
	int n;
	unsigned int cycles, cycles1, cycles2;
	unsigned char SecretKeyA[32], PublicKeyA[64], SecretAgreementA[32];
	unsigned char SecretKeyB[32], PublicKeyB[64];
	ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;

	print_test("\n--------------------------------------------------------------------------------------------------------\n");
	print_test("Benchmarking DH key exchange using uncompressed, 64-byte public keys: \n");

	cycles = 0;
	for (n = 0; n < BENCH_LOOPS; n++)
	{
		cycles1 = cpucycles();
		Status = KeyGeneration(SecretKeyA, PublicKeyA);
		if (Status != ECCRYPTO_SUCCESS) {
			return Status;
		}
		cycles2 = cpucycles();
		cycles = cycles + (cycles2 - cycles1);
	}
	print_bench("  Keypair generation runs in ...................................................... ", cycles/BENCH_LOOPS);

	Status = KeyGeneration(SecretKeyB, PublicKeyB);
	cycles = 0;
	for (n = 0; n < BENCH_LOOPS; n++)
	{
		cycles1 = cpucycles();
		Status = SecretAgreement(SecretKeyA, PublicKeyB, SecretAgreementA);
		if (Status != ECCRYPTO_SUCCESS) {
			return Status;
		}
		cycles2 = cpucycles();
		cycles = cycles + (cycles2 - cycles1);
	}
	print_bench("  Secret agreement runs in ........................................................ ", cycles/BENCH_LOOPS);

	return Status;
}


int main()
{
    clock_setup();
    gpio_setup();
    usart_setup(115200);
    rng_setup();

    *SCB_DEMCR = *SCB_DEMCR | 0x01000000;
    *DWT_CYCCNT = 0;             // reset the counter
    *DWT_CTRL = *DWT_CTRL | 1 ;  // enable the counter
    ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;

    Status = SchnorrQ_test();         // Test SchnorrQ signature scheme
    if (Status != ECCRYPTO_SUCCESS) {
        print_test("\n\n   Error detected \n\n");
        return false;
    }
    Status = SchnorrQ_run();          // Benchmark SchnorrQ signature scheme
    if (Status != ECCRYPTO_SUCCESS) {
        print_test("\n\n   Error detected \n\n");
        return false;
    }

	Status = compressedkex_test();    // Test Diffie-Hellman key exchange using compressed public keys
	if (Status != ECCRYPTO_SUCCESS) {
        print_test("\n\n   Error detected \n\n");
		return false;
	}
	Status = compressedkex_run();     // Benchmark Diffie-Hellman key exchange using compressed public keys
	if (Status != ECCRYPTO_SUCCESS) {
        print_test("\n\n   Error detected \n\n");
		return false;
	}

	Status = kex_test();              // Test Diffie-Hellman key exchange using uncompressed public keys
	if (Status != ECCRYPTO_SUCCESS) {
        print_test("\n\n   Error detected \n\n");
		return false;
	}
	Status = kex_run();               // Benchmark Diffie-Hellman key exchange using uncompressed public keys
	if (Status != ECCRYPTO_SUCCESS) {
        print_test("\n\n   Error detected \n\n");
		return false;
	}
    signal_host();

    return 0;
}