

#ifndef POINT_ARITHMETIC_HELPER_H_
#define POINT_ARITHMETIC_HELPER_H_

#include"montgomery_curve.h"

void LADDER_inner_part(proj_point_t* P0, proj_point_t* P1, proj_point_t* P2, proj_point_t* a24, uint64_t condition, fp_t* mod);

#endif