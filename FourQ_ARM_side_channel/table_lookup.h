/***********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: table lookup functions
************************************************************************************/

#ifndef __TABLE_LOOKUP_H__
#define __TABLE_LOOKUP_H__


// For C++
#ifdef __cplusplus
extern "C" {
#endif


#include "FourQ_internal.h"


void table_lookup_1x16(point_extedwards_t* table, point_extedwards_t P, unsigned int digit)
{ // Constant-time table lookup using "interleaved" masking to extract a point represented as (X,Y,Z,T) in extended twisted Edwards coordinates  
  // Inputs: digit and table containing 16 points
  // Output: P = table[digit]
    point_extedwards_t point;
    unsigned int i, j;
    digit_t mask, temp, value = 0xAAAAAAAA;
    f2elm_t tt;
                                  
    ecccopy(table[0], point);                                                 // point = table[0]

    for (i = 1; i < 16; i++)
    {
        digit--;
        // While digit>=0 mask = 0x55...5 else mask = 0xAA...A
        mask = ((digit_t)digit >> (RADIX-1)) - 1;
        mask = (mask & ~value) | (~mask & value);

        fp2copy1271(table[i]->x, tt);                                         // tt = table[i]
        for (j = 0; j < NWORDS_FIELD; j++) {
            temp = point->x[0][j] ^ tt[0][j];
            point->x[0][j] = ((mask & temp) ^ point->x[0][j]) ^ (value & temp); 
            temp = point->x[1][j] ^ tt[1][j];
            point->x[1][j] = ((mask & temp) ^ point->x[1][j]) ^ (value & temp);
        }        
        fp2copy1271(table[i]->y, tt);    
        for (j = 0; j < NWORDS_FIELD; j++) {         
            temp = point->y[0][j] ^ tt[0][j];
            point->y[0][j] = ((mask & temp) ^ point->y[0][j]) ^ (value & temp);
            temp = point->y[1][j] ^ tt[1][j];
            point->y[1][j] = ((mask & temp) ^ point->y[1][j]) ^ (value & temp);
        }        
        fp2copy1271(table[i]->z, tt);      
        for (j = 0; j < NWORDS_FIELD; j++) {            
            temp = point->z[0][j] ^ tt[0][j];
            point->z[0][j] = ((mask & temp) ^ point->z[0][j]) ^ (value & temp);
            temp = point->z[1][j] ^ tt[1][j];
            point->z[1][j] = ((mask & temp) ^ point->z[1][j]) ^ (value & temp);
        }          
        fp2copy1271(table[i]->t, tt);   
        for (j = 0; j < NWORDS_FIELD; j++) {             
            temp = point->t[0][j] ^ tt[0][j];
            point->t[0][j] = ((mask & temp) ^ point->t[0][j]) ^ (value & temp);
            temp = point->t[1][j] ^ tt[1][j];
            point->t[1][j] = ((mask & temp) ^ point->t[1][j]) ^ (value & temp);
        }
    }                          
    ecccopy(point, P); 
}


#ifdef __cplusplus
}
#endif


#endif