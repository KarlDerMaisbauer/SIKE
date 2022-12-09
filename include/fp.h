

//gcc  -Wall -std=c11 -I ../include/ -o fp fp.c fp_helper.c
#ifndef FP_H_
#define FP_H_

#include<stdint.h>



#define  WORDS 12

typedef uint64_t fp_t[WORDS];
typedef uint64_t f2p_t[4*WORDS];


//---------------------------------------------------------
//
// initializes fp_T element with the value
// stored in val
//
// @param val: value of element as string
// @param mod: modulo
// @param a  : fp element to be initialized
//
//---------------------------------------------------------
void fp_init(const char* val, fp_t* mod, fp_t* a);

//---------------------------------------------------------
//
// calculates
//
// (a + b) % mod = res
//
// @param a:   first fp_t to be added
// @param b:   second fp_t tp be added
// @param mod: modulo
// @param res: saves the result of the addition
//
//---------------------------------------------------------
void fp_addm(fp_t* a, fp_t* b, fp_t* mod, fp_t* res);

//---------------------------------------------------------
//
// calculates
//
// (a - b) % mod = res
//
// @param a:   first fp_t
// @param b:   fp_t tp be subsracted
// @param mod: modulo
// @param res: result of substraction
//
//---------------------------------------------------------
void fp_subm(fp_t* a, fp_t* b, fp_t*mod, fp_t *res);

//---------------------------------------------------------
//
// calculates
//
// (a * b) % mod = res
//
// @param a:   first fp_t
// @param b:   second fp_t
// @param mod: modulo
// @param res: result of multiplication
//
//---------------------------------------------------------
void fp_mulm(fp_t* a, fp_t* b, fp_t*mod, fp_t* res);

//---------------------------------------------------------
//
// calculates
//
// (a / b) % mod = res
//
// @param a:   dividend
// @param b:   divisor
// @param mod: modulo
// @param res: result of integer division
//
//---------------------------------------------------------
void fp_divm(fp_t* a, fp_t* b, fp_t*mod, fp_t* res);


#endif



