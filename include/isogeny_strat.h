

#ifndef ISOGENY_STRAT_H_
#define ISOGENY_STRAT_H_

#include "stdint.h"

typedef struct{
   int64_t* strategy;
   int64_t strat_len;
   int64_t queue_len;
} isogeny_strat_t;

extern const isogeny_strat_t p434_e2;
extern const isogeny_strat_t p434_e3;

extern const isogeny_strat_t p503_e2;
extern const isogeny_strat_t p503_e3;

extern const isogeny_strat_t p610_e2;
extern const isogeny_strat_t p610_e3;

extern const isogeny_strat_t p751_e2;
extern const isogeny_strat_t p751_e3;

#endif