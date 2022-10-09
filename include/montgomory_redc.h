
#ifndef MONTGOMORY_REDC_H_
#define MONTGOMORY_REDC_H_

#include"fp_helper.h"
#include"fp2.h"


fp_t n_strich;
fp_t r;
fp_t r_minus;
fp_t r2;
fp2_t one_mont;


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


// helper functions
void MonPro(fp_t* a, fp_t* b, fp_t* mod, fp_t* res);
int get_shift(fp_t* a);
#endif