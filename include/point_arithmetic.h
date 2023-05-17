

#ifndef POINT_ARITHMETIC_H_
#define POINT_ARITHMETIC_H_

#include"fp2.h"
#include"fp.h"
#include"montgomery_curve.h"




void xDBL(proj_point_t* P, mont_curve_t* A24, fp_t* mod, proj_point_t* P2);

void xDBL_no_redc(proj_point_t* P, mont_curve_t* A24, fp_t* mod, proj_point_t* P2);

void xDBLe(proj_point_t* P, mont_curve_t* A24, int64_t e, fp_t* mod, proj_point_t* P2e);

void xDBLADD(proj_point_t* P, proj_point_t* Q, proj_point_t* PmQ, fp_t* mod, mont_curve_t* A24, proj_point_t* P2, proj_point_t* PpQ);

void xTPL(proj_point_t* P, mont_curve_t* A24, fp_t* mod, proj_point_t* P3);

void xTPL_no_redc(proj_point_t* P, mont_curve_t* A24, fp_t* mod, proj_point_t* P3);

void xTPLe(proj_point_t* P, mont_curve_t* A24, int64_t e, fp_t* mod, proj_point_t* P3e);


void Ladder3p(fp2_t* P, fp2_t* Q, fp2_t* QmP, fp_t* m, proj_point_t* PpmQ, mont_curve_t* A, fp_t* mod);


#endif