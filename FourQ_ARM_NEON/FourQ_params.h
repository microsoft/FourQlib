/******************************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: FourQ's curve parameters
*
* This code is based on the papers:
* [1] "FourQ: four-dimensional decompositions on a Q-curve over the Mersenne prime" 
*     by Craig Costello and Patrick Longa, ASIACRYPT2015 (http://eprint.iacr.org/2015/565).
* [2] "FourQNEON: Faster Elliptic Curve Scalar Multiplications on ARM Processors" 
*     by Patrick Longa, SAC2016 (http://eprint.iacr.org/2016/645).
*******************************************************************************************/ 

#ifndef __FOURQ_PARAMS_H__
#define __FOURQ_PARAMS_H__

#include "FourQ_internal.h"


// Encoding of field elements, elements over Z_r and elements over GF(p^2):
// -----------------------------------------------------------------------
// Elements over GF(p) and Z_r are encoded with the least significant digit located in the leftmost position (i.e., little endian format). 
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as a||b, with a in the least significant position.
// Parameter "d" is encoded as an interleaved vector using the representation b_4|a_4|...|b_0|a_0 <- 23|23|26|26|26|26|26|26|26|26-bit,
// where the 23-bit digits are the most significant digits. Digits are stored in 32-bit words.

static const uint32_t PARAMETER_d[10]      = { 0x00000142, 0x01FC0C8D, 0x00000000, 0x0085223C, 0x000E4000, 0x020FCB38, 0x00000000, 0x0211995F, 0x00000000, 0x005E472F };
static const uint64_t GENERATOR_x[4]       = { 0x286592AD7B3833AA, 0x1A3472237C2FB305, 0x96869FB360AC77F6, 0x1E1F553F2878AA9C };
static const uint64_t GENERATOR_y[4]       = { 0xB924A2462BCBB287, 0x0E3FEE9BA120785A, 0x49A7C344844C8B5C, 0x6E1C4AF8630E0242 };
static const uint64_t curve_order[4]       = { 0x2FB2540EC7768CE7, 0xDFBD004DFE0F7999, 0xF05397829CBC14E5, 0x0029CBC14E5E0A72 };
static const uint64_t Montgomery_Rprime[4] = { 0xC81DB8795FF3D621, 0x173EA5AAEA6B387D, 0x3D01B7C72136F61C, 0x0006A5F16AC8F9D3 };
static const uint64_t Montgomery_rprime[4] = { 0xE12FE5F079BC3929, 0xD75E78B8D1FCDCF3, 0xBCE409ED76B5DB21, 0xF32702FDAFC1C074 };

                         
#endif