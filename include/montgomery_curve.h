

#ifndef MONTGOMERY_CURVE_H_
#define MONTGOMERY_CURVE_H_


#include"fp2.h"
#include"fp.h"
#include"isogeny_strat.h"

typedef struct{
    fp2_t X;
    fp2_t Z;
} proj_point_t;

typedef struct{
    fp2_t A;
    fp2_t C;

} mont_curve_t;

typedef struct{
    fp2_t x;
    fp2_t y;

} point_t;



void jinvariant(fp2_t* A, fp2_t* C, fp_t* mod, fp2_t* j);

void get_A(fp2_t* xp, fp2_t* xq, fp2_t* xqp, fp_t* mod, fp2_t* A);

void iso_2_curve(proj_point_t* P2, fp_t* mod, mont_curve_t* A24);

void iso_2_eval(proj_point_t* P2,proj_point_t* Q, fp_t* mod, proj_point_t* Q_strich);

void iso_4_curve(proj_point_t* P4, fp_t* mod, mont_curve_t* A24, fp2_t* K1, fp2_t* K2, fp2_t* K3);

void iso_4_eval(fp2_t* K1, fp2_t* K2, fp2_t* K3, proj_point_t* Q, fp_t* mod, proj_point_t* Q_strich);

void iso_3_curve(proj_point_t* P3, fp_t* mod, mont_curve_t* A24, fp2_t* K1, fp2_t* K2);

void iso_3_eval(fp2_t* K1, fp2_t* K2, proj_point_t* Q, fp_t* mod, proj_point_t* Q_strich);

//---------------------------------------------------------------
// computes large isogeny
// note: opt_input will be changed in computation
//---------------------------------------------------------------
void e_2_iso(mont_curve_t* A24p, proj_point_t* S, proj_point_t* opt_input, int64_t opt_input_size, fp_t* mod, isogeny_strat_t strategy ,mont_curve_t* A24p_new);

//---------------------------------------------------------------
// computes large isogeny
// note: opt_input will be changed in computation
//---------------------------------------------------------------
void e_3_iso(mont_curve_t* A24, proj_point_t* S, proj_point_t* opt_input, int64_t opt_input_size, fp_t* mod, isogeny_strat_t strategy ,mont_curve_t* A24_new);



//helper functions

void proj_pt_copy_masked(proj_point_t* a, proj_point_t* b, uint64_t mask);

void proj_pt_copy(proj_point_t* a, proj_point_t* b);

void proj_pt_zero(proj_point_t* a);

void curve_copy(mont_curve_t* a, mont_curve_t* b);

void curve_zero(mont_curve_t* a);

#endif