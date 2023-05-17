
#ifndef MONTGOMORY_REDC_H_
#define MONTGOMORY_REDC_H_

#include "fp_helper.h"
#include "fp2.h"
#include "montgomery_curve.h"


fp_t n_strich;
fp2_t n_strich2;
fp_t r;
fp_t r_minus;
fp_t r2;
fp2_t r22;
fp2_t one_mont;
fp2_t two_mont;
fp2_t four_mont;
fp2_t six_mont;
fp2_t eight_mont;
f2p_t validity_tester;


//---------------------------------------------------------
//
// calculates greatest common divisor 
//
//
// @param m:     first fp_t
// @param n:     second fp_t
// @param m_mul: multiplier of m
// @param n_mul: rmultiplier of n
//
//---------------------------------------------------------
void gcd(fp_t* m, fp_t* n, fp_t* m_mul, fp_t* n_mul);

//---------------------------------------------------------
//
// initializes values for montgomor reduction
//
//
// @param m: modulo
//
//---------------------------------------------------------
void init(fp_t* mod);

//---------------------------------------------------------
//
// montgomory reduction for fp_t
//
// not working
//
//
// @param T:   value to be reduced
// @param mod: modulo
// @param res: result of reduction
//
//---------------------------------------------------------
void REDC(fp_t* T, fp_t* mod, fp_t* res);


//---------------------------------------------------------
//
// montgomory reduction for f2p_t
//
//
// @param T:   value to be reduced
// @param mod: modulo
// @param res: result of reduction
//
//---------------------------------------------------------
void REDCL(f2p_t* T, fp_t* mod, fp_t* res);

//---------------------------------------------------------
//
// montgomory reduction for f2p2_t
//
//
// @param T:   value to be reduced
// @param mod: modulo
// @param res: result of reduction
//
//---------------------------------------------------------
void REDCL2(f2p2_t* T, fp_t* mod, fp2_t* res);


//---------------------------------------------------------
//
// montgomory reduction for fp2_t
//
//
// @param T:   value to be reduced
// @param mod: modulo
// @param res: result of reduction
//
//---------------------------------------------------------
void REDC2(fp2_t* T, fp_t* mod, fp2_t* res);


//---------------------------------------------------------
//
// performs montogomory multiplication
//
//
// @param a:   first value
// @param b:   second value
// @param mod: modulo
// @param res: result of reduction
//
//---------------------------------------------------------
void MODMUL(fp_t* a, fp_t* b, fp_t* mod, fp_t* res);

//---------------------------------------------------------
//
// performs montogomory multiplication (oter version)
//
//
// @param a:   first value
// @param b:   second value
// @param mod: modulo
// @param res: result of reduction
//
//---------------------------------------------------------
void ModMul(fp_t* a, fp_t* b, fp_t* mod, fp_t* res);

//---------------------------------------------------------
// Conversion


// To montgomory form
//---------------------------------------------------------
//
// converts fp_t into montgomory form
//
//
// @param a:        fp_t
// @param mod:      modulo
// @param a_redc:   fp_t in montgomory form
//
//---------------------------------------------------------
void fp_to_mont(fp_t* a, fp_t* mod, fp_t* a_redc);

//---------------------------------------------------------
//
// converts fp2_t into montgomory form
//
//
// @param a:        fp2_t
// @param mod:      modulo
// @param a_redc:   fp2_t in montgomory form
//
//---------------------------------------------------------
void fp2_to_mont(fp2_t* a, fp_t* mod, fp2_t* a_redc);

//---------------------------------------------------------
//
// converts proj_pt_t into montgomory form
//
//
// @param a:        proj_point_t
// @param mod:      modulo
// @param a_redc:   proj_point_t in montgomory form
//
//---------------------------------------------------------
void pt_to_mont(proj_point_t* a, fp_t* mod, proj_point_t* a_redc);

//From montogomry form
//---------------------------------------------------------
//
// converts fp_t from montgomory form
//
//
// @param a_redc:   fp_t in montgomory form
// @param mod:      modulo
// @param a:        fp_t
//
//---------------------------------------------------------
void fp_from_mont(fp_t* a_redc, fp_t* mod, fp_t* a);

//---------------------------------------------------------
//
// converts f2p_t from montgomory form
//
//
// @param a_redc:   f2p_t in montgomory form
// @param mod:      modulo
// @param a:        f2p_t
//
//---------------------------------------------------------
void f2p_from_mont(f2p_t* a_redc, fp_t* mod, fp_t* a);

//---------------------------------------------------------
//
// converts fp2_t from montgomory form
//
//
// @param a_redc:   fp2_t in montgomory form
// @param mod:      modulo
// @param a:        fp2_t
//
//---------------------------------------------------------
void fp2_from_mont(fp2_t* a_redc, fp_t* mod, fp2_t* a);

//---------------------------------------------------------
//
// converts f2p2_t from montgomory form
//
//
// @param a_redc:   f2p2_t in montgomory form
// @param mod:      modulo
// @param a:        fp2_t
//
//---------------------------------------------------------
void f2p2_from_mont(f2p2_t* a_redc, fp_t* mod, fp2_t* a);

//---------------------------------------------------------
//
// converts proj_point_t from montgomory form
//
//
// @param a_redc:   proj_point_t in montgomory form
// @param mod:      modulo
// @param a:        proj_point_t
//
//---------------------------------------------------------
void pt_from_mont(proj_point_t* a_redc, fp_t* mod, proj_point_t* a);




// helper functions
void MonPro(fp_t* a, fp_t* b, fp_t* mod, fp_t* res);
int get_shift(fp_t* a);
#endif