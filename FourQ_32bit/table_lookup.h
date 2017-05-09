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


void table_lookup_1x8(vpoint_extproj_precomp_t* table, vpoint_extproj_precomp_t P, unsigned int digit, unsigned int sign_mask)
{ // Constant-time table lookup to extract a point represented as (X+Y,Y-X,2Z,2dT) corresponding to extended twisted Edwards coordinates (X:Y:Z:T)
  // Inputs: sign_mask, digit, table containing 8 points
  // Output: P = sign*table[digit], where sign=1 if sign_mask=0xFF...FF and sign=-1 if sign_mask=0

#if defined(NO_CACHE_MEM)
    vpoint_extproj_precomp_t temp_point[2];

    ecccopy_precomp(table[digit], temp_point[1]);                            // temp_point[1] = table[digit]
    v2copy1271(temp_point[1]->xy, temp_point[0]->yx);                        // temp_point[0] = -table[digit], i.e., coordinates (y-x,x+y,2*z,-2dt)
    v2copy1271(temp_point[1]->yx, temp_point[0]->xy);     
    v2copy1271(temp_point[1]->t2, temp_point[0]->t2);                                   
    v2copy1271(temp_point[1]->z2, temp_point[0]->z2);                                                                 
    vneg1271(&temp_point[0]->t2[0]);                                          
    vneg1271(&temp_point[0]->t2[VWORDS_FIELD]);   
    ecccopy_precomp(temp_point[sign_mask & 1], P);           

#else
    vpoint_extproj_precomp_t point, temp_point;
    unsigned int i, j;
    digit_t mask;
                                  
    ecccopy_precomp(table[0], point);                                        // point = table[0]

    for (i = 1; i < 8; i++)
    {
        digit--;
        // While digit>=0 mask = 0xFF...F else sign = 0x00...0
        mask = ((digit_t)digit >> (RADIX-1)) - 1;
        ecccopy_precomp(table[i], temp_point);                               // temp_point = table[i] 
        // If mask = 0x00...0 then point = point, else if mask = 0xFF...F then point = temp_point            
        point->xy[0] = (mask & (point->xy[0] ^ temp_point->xy[0])) ^ point->xy[0];
        point->xy[1] = (mask & (point->xy[1] ^ temp_point->xy[1])) ^ point->xy[1];
        point->xy[2] = (mask & (point->xy[2] ^ temp_point->xy[2])) ^ point->xy[2];
        point->xy[3] = (mask & (point->xy[3] ^ temp_point->xy[3])) ^ point->xy[3];
        point->xy[4] = (mask & (point->xy[4] ^ temp_point->xy[4])) ^ point->xy[4];
        point->xy[5] = (mask & (point->xy[5] ^ temp_point->xy[5])) ^ point->xy[5];
        point->xy[6] = (mask & (point->xy[6] ^ temp_point->xy[6])) ^ point->xy[6];
        point->xy[7] = (mask & (point->xy[7] ^ temp_point->xy[7])) ^ point->xy[7];
        point->xy[8] = (mask & (point->xy[8] ^ temp_point->xy[8])) ^ point->xy[8];
        point->xy[9] = (mask & (point->xy[9] ^ temp_point->xy[9])) ^ point->xy[9];
        point->yx[0] = (mask & (point->yx[0] ^ temp_point->yx[0])) ^ point->yx[0];
        point->yx[1] = (mask & (point->yx[1] ^ temp_point->yx[1])) ^ point->yx[1];
        point->yx[2] = (mask & (point->yx[2] ^ temp_point->yx[2])) ^ point->yx[2];
        point->yx[3] = (mask & (point->yx[3] ^ temp_point->yx[3])) ^ point->yx[3];
        point->yx[4] = (mask & (point->yx[4] ^ temp_point->yx[4])) ^ point->yx[4];
        point->yx[5] = (mask & (point->yx[5] ^ temp_point->yx[5])) ^ point->yx[5];
        point->yx[6] = (mask & (point->yx[6] ^ temp_point->yx[6])) ^ point->yx[6];
        point->yx[7] = (mask & (point->yx[7] ^ temp_point->yx[7])) ^ point->yx[7];
        point->yx[8] = (mask & (point->yx[8] ^ temp_point->yx[8])) ^ point->yx[8];
        point->yx[9] = (mask & (point->yx[9] ^ temp_point->yx[9])) ^ point->yx[9];
        point->z2[0] = (mask & (point->z2[0] ^ temp_point->z2[0])) ^ point->z2[0];
        point->z2[1] = (mask & (point->z2[1] ^ temp_point->z2[1])) ^ point->z2[1];
        point->z2[2] = (mask & (point->z2[2] ^ temp_point->z2[2])) ^ point->z2[2];
        point->z2[3] = (mask & (point->z2[3] ^ temp_point->z2[3])) ^ point->z2[3];
        point->z2[4] = (mask & (point->z2[4] ^ temp_point->z2[4])) ^ point->z2[4];
        point->z2[5] = (mask & (point->z2[5] ^ temp_point->z2[5])) ^ point->z2[5];
        point->z2[6] = (mask & (point->z2[6] ^ temp_point->z2[6])) ^ point->z2[6];
        point->z2[7] = (mask & (point->z2[7] ^ temp_point->z2[7])) ^ point->z2[7];
        point->z2[8] = (mask & (point->z2[8] ^ temp_point->z2[8])) ^ point->z2[8];
        point->z2[9] = (mask & (point->z2[9] ^ temp_point->z2[9])) ^ point->z2[9];
        point->t2[0] = (mask & (point->t2[0] ^ temp_point->t2[0])) ^ point->t2[0];
        point->t2[1] = (mask & (point->t2[1] ^ temp_point->t2[1])) ^ point->t2[1];
        point->t2[2] = (mask & (point->t2[2] ^ temp_point->t2[2])) ^ point->t2[2];
        point->t2[3] = (mask & (point->t2[3] ^ temp_point->t2[3])) ^ point->t2[3];
        point->t2[4] = (mask & (point->t2[4] ^ temp_point->t2[4])) ^ point->t2[4];
        point->t2[5] = (mask & (point->t2[5] ^ temp_point->t2[5])) ^ point->t2[5];
        point->t2[6] = (mask & (point->t2[6] ^ temp_point->t2[6])) ^ point->t2[6];
        point->t2[7] = (mask & (point->t2[7] ^ temp_point->t2[7])) ^ point->t2[7];
        point->t2[8] = (mask & (point->t2[8] ^ temp_point->t2[8])) ^ point->t2[8];
        point->t2[9] = (mask & (point->t2[9] ^ temp_point->t2[9])) ^ point->t2[9];
    }
    
    v2copy1271(point->t2, temp_point->t2);
    v2copy1271(point->xy, temp_point->yx);                                   // point: x+y,y-x,2dt coordinate, temp_point: y-x,x+y,-2dt coordinate
    v2copy1271(point->yx, temp_point->xy);                                   
    vneg1271(&temp_point->t2[0]);                                            // Negate 2dt coordinate
    vneg1271(&temp_point->t2[VWORDS_FIELD]);             
    for (j = 0; j < 2*VWORDS_FIELD; j++) {                                   // If sign_mask = 0 then choose negative of the point
        point->xy[j] = ((digit_t)((int)sign_mask) & (point->xy[j] ^ temp_point->xy[j])) ^ temp_point->xy[j];
        point->yx[j] = ((digit_t)((int)sign_mask) & (point->yx[j] ^ temp_point->yx[j])) ^ temp_point->yx[j];
        point->t2[j] = ((digit_t)((int)sign_mask) & (point->t2[j] ^ temp_point->t2[j])) ^ temp_point->t2[j];
    }                                
    ecccopy_precomp(point, P); 
#endif
}


void table_lookup_fixed_base(vpoint_precomp_t* table, vpoint_precomp_t P, unsigned int digit, unsigned int sign)
{ // Constant-time table lookup to extract a point represented as (x+y,y-x,2t) corresponding to extended twisted Edwards coordinates (X:Y:Z:T) with Z=1
  // Inputs: sign, digit, table containing VPOINTS_FIXEDBASE = 2^(W_FIXEDBASE-1) points
  // Output: if sign=0 then P = table[digit], else if (sign=-1) then P = -table[digit]

#if defined(NO_CACHE_MEM)
    vpoint_precomp_t temp_point[2];

    ecccopy_precomp_fixed_base(table[digit], temp_point[0]);                 // temp_point[0] = table[digit]
    v2copy1271(temp_point[0]->xy, temp_point[1]->yx);                        // temp_point[1] = -table[digit], i.e., coordinates (y-x,x+y,2*z,-2dt)
    v2copy1271(temp_point[0]->yx, temp_point[1]->xy);     
    v2copy1271(temp_point[0]->t2, temp_point[1]->t2);                                                                 
    vneg1271(&temp_point[1]->t2[0]);                                          
    vneg1271(&temp_point[1]->t2[VWORDS_FIELD]);   
    ecccopy_precomp_fixed_base(temp_point[sign & 1], P);    

#else
    vpoint_precomp_t point, temp_point;
    unsigned int i, j;
    digit_t mask;
                                   
    ecccopy_precomp_fixed_base(table[0], point);                             // point = table[0]

    for (i = 1; i < VPOINTS_FIXEDBASE; i++)
    {
        digit--;
        // While digit>=0 mask = 0xFF...F else sign = 0x00...0
        mask = ((digit_t)digit >> (RADIX-1)) - 1;
        ecccopy_precomp_fixed_base(table[i], temp_point);                    // temp_point = table[i] 
        // If mask = 0x00...0 then point = point, else if mask = 0xFF...F then point = temp_point
        for (j = 0; j < 2*VWORDS_FIELD; j++) {
            point->xy[j] = (mask & (point->xy[j] ^ temp_point->xy[j])) ^ point->xy[j];
            point->yx[j] = (mask & (point->yx[j] ^ temp_point->yx[j])) ^ point->yx[j];
            point->t2[j] = (mask & (point->t2[j] ^ temp_point->t2[j])) ^ point->t2[j];
        }
    }
    
    v2copy1271(point->t2, temp_point->t2);
    v2copy1271(point->xy, temp_point->yx);                                  // point: x+y,y-x,2dt coordinate, temp_point: y-x,x+y,-2dt coordinate
    v2copy1271(point->yx, temp_point->xy);                                   
    vneg1271(&temp_point->t2[0]);                                            // Negate 2dt coordinate
    vneg1271(&temp_point->t2[VWORDS_FIELD]);             
    for (j = 0; j < 2*VWORDS_FIELD; j++) {                                     // If sign = 0xFF...F then choose negative of the point
        point->xy[j] = ((digit_t)((int)sign) & (point->xy[j] ^ temp_point->xy[j])) ^ point->xy[j];
        point->yx[j] = ((digit_t)((int)sign) & (point->yx[j] ^ temp_point->yx[j])) ^ point->yx[j];
        point->t2[j] = ((digit_t)((int)sign) & (point->t2[j] ^ temp_point->t2[j])) ^ point->t2[j];
    }                                  
    ecccopy_precomp_fixed_base(point, P);
#endif
}


#ifdef __cplusplus
}
#endif


#endif
