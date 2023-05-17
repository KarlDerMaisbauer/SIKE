

#ifndef FP2_H_
#define FP2_H_

#include"fp.h"


#include <stdio.h>

typedef struct{
    fp_t real;
    fp_t img;
}fp2_t;

typedef struct{
    f2p_t real;
    f2p_t img;
}f2p2_t;


//---------------------------------------------------------
//
// calculates
//
// (a + b) % mod = res
//
// @param a:   first fp2_t to be added
// @param b:   second fp2_t tp be added
// @param mod: modulo
// @param res: saves the result of the addition
//
//---------------------------------------------------------
void fp2_addm(fp2_t* a, fp2_t*b, fp_t* mod, fp2_t* res);

//---------------------------------------------------------
//
// calculates
//
// (a - b) % mod = res
//
// @param a:   first fp2_t
// @param b:   fp2_t tp be subsracted
// @param mod: modulo
// @param res: result of substraction
//
//---------------------------------------------------------
void fp2_subm(fp2_t* a, fp2_t* b, fp_t*mod, fp2_t *res);

//---------------------------------------------------------
//
// calculates
//
// (a * b) % mod = res
//
// @param a:   first fp2_t
// @param b:   second fp2_t
// @param mod: modulo
// @param res: result of multiplication
//
//---------------------------------------------------------
void fp2_mulm(fp2_t* a, fp2_t* b, fp_t*mod, fp2_t* res);

//---------------------------------------------------------
//
// calculates
//
// -a
//
// @param a:   dividend
// @param mod: modulo
// @param res: additive inverse of a
//
//---------------------------------------------------------
void fp2_add_inv(fp2_t* a, fp_t*mod, fp2_t* res);

//---------------------------------------------------------
//
// calculates
//
// a^(-1)
//
// @param a:   fp2
// @param mod: modulo
// @param res: multiplicative inverse of a
//
//---------------------------------------------------------
void fp2_mult_inv(fp2_t* a, fp_t*mod, fp2_t* res);

//---------------------------------------------------------
//
// calculates
//
// a * b^(-1) = res in fp2
//
// @param a:   dividend
// @param b:   divisor
// @param mod: modulo
// @param res: result
//
//---------------------------------------------------------
void fp2_divm(fp2_t* a, fp2_t* b, fp_t*mod, fp2_t* res);

//---------------------------------------------------------
//
// calculates
//
// (a + b) = res
//
// @param a:   first fp2_t to be added
// @param b:   second fp2_t tp be added
// @param res: saves the result of the addition
//
//---------------------------------------------------------
void fp2_add(fp2_t* a, fp2_t*b, fp2_t* res);

//---------------------------------------------------------
//
// calculates
//
// (a - b) = res
//
// @param a:   first fp2_t
// @param b:   fp2_t tp be subsracted
// @param res: result of substraction
//
//---------------------------------------------------------
void fp2_sub(fp2_t* a, fp2_t* b, fp2_t *res);

//---------------------------------------------------------
//
// performs multiplication 
//
// a * b = res
//
// @param a:   first fp2_t
// @param b:   second fp2_t
// @param res: result of multiplication
//
//---------------------------------------------------------
void fp2_mul(fp2_t* a, fp2_t* b, f2p2_t* res);

//---------------------------------------------------------
//
// performs multiplication in montgomory form
//
// a * b = res
//
// @param a:   first fp2_t in montgomery form
// @param b:   second fp2_t in montgomery form
// @param mod: modulo 
// @param res: result of multiplication
//
//---------------------------------------------------------
void fp2_mul_mont(fp2_t* a, fp2_t* b, fp_t* mod, fp2_t* res);

//---------------------------------------------------------
//
// calculates
//
// a^(-1)
//
// in montgomory form
//
// @param a:   fp2
// @param mod: modulo
// @param res: multiplicative inverse of a
//
//---------------------------------------------------------
void fp2_mult_inv_mont(fp2_t* a, fp_t*mod, fp2_t* res);

//---------------------------------------------------------
//
// performs division in montgomory form
//
// a / b = res
//
// @param a:   first fp2_t in montgomery form
// @param b:   second fp2_t in montgomery form
// @param mod: modulo 
// @param res: result of division
//
//---------------------------------------------------------
void fp2_div_mont(fp2_t* a, fp2_t* b, fp_t* mod, fp2_t* res);

//---------------------------------------------------------
//
// sets fp2_t to 0
//
//
// @param a:   fp2_t to be set to zero
//
//---------------------------------------------------------
void fp2_zero(fp2_t* a);

//---------------------------------------------------------
//
// copies values from a into b
//
// @param a:        number which is copied
// @param a_copy:   copy of a
//
//---------------------------------------------------------
void fp2_copy(fp2_t* a, fp2_t* a_copy);

//---------------------------------------------------------
//
// copies values from a into b defined by mask
//
// @param a:    number which is copied
// @param b:    masked copy of a
// @param mask: mask (if pos is 1 => bit is copied)
//
//---------------------------------------------------------
void fp2_copy_masked(fp2_t* a, fp2_t* b, uint64_t mask);

//---------------------------------------------------------
//
// checks if a >= b
// only for positive fp2_t 
//
// @param a: first number to be tested
// @param b: secont number to be tested
//
// @return 1 if a >= b, 0 else
// 
//---------------------------------------------------------
int fp2_greater_equ_pos(fp2_t* a, fp2_t* b);

int fp2_equ(fp2_t* a, fp2_t* b);

void fp2_to_f2p2(fp2_t* a, f2p2_t* res);


void f2p2_zero(f2p2_t* a);

void fp2_fprintf(FILE* file, fp2_t* p);

void f2p2_fprintf(FILE* file, f2p2_t* p);

void fp2_print(fp2_t* a);


//---------------------------------------------------------
//
// calculates
//
// (a + b) = res
//
// @param a:   first f2p2_t to be added
// @param b:   second f2p2_t tp be added
// @param res: saves the result of the addition
//
//---------------------------------------------------------
void f2p2_add(f2p2_t* a, f2p2_t*b, f2p2_t* res);

//---------------------------------------------------------
//
// calculates
//
// a^(-1)
//
// @param a:   fp2
// @param mod: modulo
// @param res: multiplicative inverse of a
//
//---------------------------------------------------------
void f2p2_mult_inv(f2p2_t* a, fp_t*mod, f2p2_t* res);


void f2p2_mul(f2p2_t* a, f2p2_t* b, f2p2_t* res);


//---------------------------------------------------------
//
// calculates
//
// a * b^(-1) = res in fp2
//
// @param a:   dividend
// @param b:   divisor
// @param mod: modulo
// @param res: result
//
//---------------------------------------------------------
void f2p2_div(f2p2_t* a, f2p2_t* b, fp_t*mod, f2p2_t* res);

void f2p2_mod(f2p2_t*a , fp_t*mod, fp2_t* res);


#endif

