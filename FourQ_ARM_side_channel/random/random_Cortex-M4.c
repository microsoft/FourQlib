/***********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: pseudorandom generation functions
************************************************************************************/ 

#include "../FourQ_internal.h"
#include "../stm32f4_wrapper.h"


ECCRYPTO_STATUS random_bytes(unsigned char* random_array, unsigned int nbytes)
{ // Generation of "nbytes" of random bytes

    if ((nbytes & 0x03) != 0) return ECCRYPTO_ERROR_INVALID_PARAMETER;

    rng_setup();
    random_int((uint32_t*)random_array, nbytes/4);
	return ECCRYPTO_SUCCESS;
}
