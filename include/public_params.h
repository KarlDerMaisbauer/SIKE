

#ifndef PUBLIC_PARAMS_H_
#define PUBLIC_PARAMS_H_

#include "fp.h"
#include "fp2.h"

typedef struct{
    fp_t e2;
    fp_t e3;
    fp_t p;
    fp2_t xP2;
    fp2_t xQ2;
    fp2_t xR2;
    fp2_t xP3;
    fp2_t xQ3;
    fp2_t xR3;
} public_params_t;


#endif