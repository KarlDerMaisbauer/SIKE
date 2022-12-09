

#ifndef PUBLIC_KEY_H_
#define PUBLIC_KEY_H_

#include "fp2.h"

typedef struct{
    fp2_t x1;
    fp2_t x2;
    fp2_t x3;
} public_key_t;

#endif