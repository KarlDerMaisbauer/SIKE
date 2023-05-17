
#ifndef FP_HELPER_H_
#define FP_HELPER_H_


#include "fp.h"
#include "stdio.h"

#define LINE printf("\n");

//
//
//  Helper functions for the fp_t datatype
//
//


//---------------------------------------------------------
//
// sets given variable to zero
//
// @param a: pointer to fp_t to be set to zero
//
//---------------------------------------------------------
void fp_zero(fp_t* a);


//---------------------------------------------------------
//
// prints given variable
//
// @param a: pointer to fp_t to be printed
//
//---------------------------------------------------------
void fp_print(fp_t* p);

//---------------------------------------------------------
//
// checks if fp_t is smaller zero 
//
// @param a: number to be tested
//
// @return 1 for smaller zero, 0 else
//
//---------------------------------------------------------
int fp_smaller_zero(fp_t* a);

//---------------------------------------------------------
//
// checks if a >= b
// only for positive fp_t 
//
// @param a: first number to be tested
// @param b: secont number to be tested
//
// @return 1 if a >= b, 0 else
// 
//---------------------------------------------------------
int fp_greater_equ_pos(fp_t* a, fp_t* b);


//---------------------------------------------------------
//
// checks if a == b
//
// @param a: first number to be tested
// @param b: secont number to be tested
//
// @return 1 if a == b, 0 else
// 
//---------------------------------------------------------
int fp_equ(fp_t* a, fp_t* b);

//---------------------------------------------------------
//
// checks if a > b
// a, b must be positive
//
// @param a: first number to be tested
// @param b: secont number to be tested
//
// @return 1 if a > b, 0 else
// 
//---------------------------------------------------------
int fp_greater_pos(fp_t* a, fp_t* b);

//---------------------------------------------------------
//
// copies values from a into b
//
// @param a: number which is copied
// @param b: copy of a
//
//---------------------------------------------------------
void fp_copy(fp_t* a, fp_t* b);

//---------------------------------------------------------
//
// copies values from a into b with respect to mask
//
// @param a:    number which is copied
// @param b:    copy of a
// @param mask: defines which parts should be copied
//
//---------------------------------------------------------
void fp_copy_masked(fp_t* a, fp_t* b, uint64_t mask);

//---------------------------------------------------------
//
// shifts given number shift times to the left
//
// @param a:     pointer to fp_t to be shifted
// @param shift: how often the number has to be shifted
//
//---------------------------------------------------------
void fp_lshift(fp_t *a, uint8_t shift);

//---------------------------------------------------------
//
// shifts given number shift times to the right
//
// @param a:     pointer to fp_t to be shifted
// @param shift: how often the number has to be shifted
//
//---------------------------------------------------------
void fp_rshift(fp_t *a, uint8_t shift);


//---------------------------------------------------------
//
// gives number of bits used by a
//
// @param a:     fp_t
//
// @return       number of bits required by a
//
//---------------------------------------------------------
int64_t fp_get_len(fp_t* a);


//
//
//  basic arithmetic functions for the fp_t datatype
//
//


//---------------------------------------------------------
//
// adds two fp_t together 
//
// a + b = res
//
// @param a:   first fp_t to be added
// @param b:   second fp_t tp be added
// @param res: saves the result of the addition
//
//---------------------------------------------------------
void fp_add(fp_t* a, fp_t* b, fp_t* res);

//---------------------------------------------------------
//
// calculates additive inverse of a
//
// -a = res
//
// @param a:   number to be inverted
// @param res: saves the additive inverse of a
//
//---------------------------------------------------------
void fp_additive_inverse(fp_t* a, fp_t* res);

//---------------------------------------------------------
//
// substracts two fp_t and saves result in res
//
// a - b = res
//
// @param a:   first fp_t
// @param b:   fp_t tp be subsracted
// @param res: result of substraction
//
//---------------------------------------------------------
void fp_sub(fp_t* a, fp_t* b, fp_t *res);


//---------------------------------------------------------
//
// multiplies two fp_t 
//
// a * b = res
//
// @param a:   first fp_t
// @param b:   second fp_t
// @param res: result of multiplication
//             (has doubled size in bit)
//
//---------------------------------------------------------
void fp_mul(fp_t* a, fp_t* b, f2p_t* res);

//---------------------------------------------------------
//
// performs integer division of two fp_t
//
// a / b = res
//
// @param a:   dividend
// @param b:   divisor
// @param res: result of integer division
//
//---------------------------------------------------------
void fp_div(fp_t* a, fp_t* b, fp_t* res);

//---------------------------------------------------------
//
// performs modulo operation
//
// a % mod = res
//
// @param a:   fp_t
// @param mod: modulo
// @param res: result of modulo
//
//---------------------------------------------------------
void fp_mod(fp_t* a, fp_t* mod, fp_t* res);

//---------------------------------------------------------
//
// performs integer division of two fp_t
//
// a / b = res
//
// @param a:   dividend
// @param b:   divisor
// @param res: result of integer division
//
//-----------------------fp_ron where mod is 2^k
//
// a % mod = res
//
// @param a:   fp_t
// @param mod: modulo
// @param res: result of modulo
//
//---------------------------------------------------------
void fp_mod_2k(fp_t*a, fp_t* mod, fp_t* res);


//
// Functions for f2p_t datatype
// (same as fp_t but with doubled bit)
//


//
//
//  Helper functions for the f2p_t datatype
//
//


//---------------------------------------------------------
//
// sets given variable to zero
//
// @param a: pointer to f2p_t to be set to zero
//
//---------------------------------------------------------
void f2p_zero(f2p_t* a);

//---------------------------------------------------------
//
// prints given variable
//
// @param a: pointer to f2p_t to be printed
//
//---------------------------------------------------------
void f2p_print(f2p_t* p);

//---------------------------------------------------------
//
// shifts given number shift times to the left
//
// @param a:     pointer to f2p_t to be shifted
// @param shift: how often the number has to be shifted
//
//---------------------------------------------------------
void f2p_lshift(f2p_t *a, uint8_t shift);

//---------------------------------------------------------
//
// turns f2p_t onto fp_t
// safes lower part of f2p_t into fp_t 
//
// @param a: larger f2p_t type
// @param b: smaller fp_t type
//
//---------------------------------------------------------
void f2p_to_fp(f2p_t* a, fp_t* b);

//---------------------------------------------------------
//
// turns fp_t into f2p_t
// safes lower part of f2p_t into fp_t 
//
// @param a: smaller fp_t type
// @param b: larger f2p_t type
//
//---------------------------------------------------------
void fp_to_f2p(fp_t* a, f2p_t* b);

//---------------------------------------------------------
//
// checks if f2p_t is smaller zero 
//
// @param a: number to be tested
//
// @return 1 for smaller zero, 0 else
//
//---------------------------------------------------------
int f2p_smaller_zero(f2p_t* a);

//---------------------------------------------------------
//
// checks if a >= b
// only for positive fp_t 
//
// @param a: first number to be tested
// @param b: secont number to be tested
//
// @return 1 if a >= b, 0 else
// 
//---------------------------------------------------------
int f2p_greater_equ_pos(f2p_t* a, f2p_t* b);

//---------------------------------------------------------
//
// checks if a == b
//
// @param a: first number to be tested
// @param b: secont number to be tested
//
// @return 1 if a == b, 0 else
// 
//---------------------------------------------------------
int f2p_equ(f2p_t*a, f2p_t*b);

//---------------------------------------------------------
//
// copies values from a into b
//
// @param a: number which is copied
// @param b: copy of a
//
//---------------------------------------------------------
void f2p_copy(f2p_t* a, f2p_t* b);

//---------------------------------------------------------
//
// shifts given number shift times to the right
//
// @param a:     pointer to f2p_t to be shifted
// @param shift: how often the number has to be shifted
//
//---------------------------------------------------------
void f2p_rshift(f2p_t *a, uint8_t shift);



//
//
//  basic arithmetic functions for the fp_t datatype
//
//


//---------------------------------------------------------
//
// adds two f2p_t together 
//
// a + b = res
//
// @param a:   first f2p_t to be added
// @param b:   second f2p_t tp be added
// @param res: saves the result of the addition
//
//---------------------------------------------------------
void f2p_add(f2p_t* a, f2p_t* b, f2p_t* res);

//---------------------------------------------------------
//
// adds two f2p_t together 
//
// a + b = res
//
// @param a:   first f2p_t to be added
// @param b:   second f2p_t tp be added
// @param res: saves the result of the addition
//
//---------------------------------------------------------
void f2p_addm(f2p_t* a, f2p_t* b, fp_t* mod, fp_t* res);


//---------------------------------------------------------
//
// calculates additive inverse of a
//
// -a = res
//
// @param a:   number to be inverted
// @param res: saves the additive inverse of a
//
//---------------------------------------------------------
void f2p_additive_inverse(f2p_t* a, f2p_t* res);

//---------------------------------------------------------
//
// substracts two fp_t and saves result in res
//
// a - b = res
//
// @param a:   first f2p_t
// @param b:   f2p_t tp be subsracted
// @param res: result of substraction
//
//---------------------------------------------------------
void f2p_sub(f2p_t*a, f2p_t*b, f2p_t*res);

//---------------------------------------------------------
//
// multiplies two f2p_t 
//
// a * b = res
//
// @param a:   first f2p_t
// @param b:   second f2p_t
// @param res: result of multiplication
//             
//
//---------------------------------------------------------
void f2p_mul(f2p_t* a, f2p_t* b, f2p_t* res);

//---------------------------------------------------------
//
// multiplies two f2p_t 
//
// a * b % mod = res
//
// @param a:   first f2p_t
// @param b:   second f2p_t
// @param res: result of multiplication
//             
//
//---------------------------------------------------------
void f2p_mulm(f2p_t* a, f2p_t* b, fp_t* mod, fp_t* res);

//---------------------------------------------------------
//
// performs integer division of two f2p_t
//
// a / b = res
//
// @param a:   dividend
// @param b:   divisor
// @param res: result of integer division
//
//---------------------------------------------------------
void f2p_div(f2p_t* a, f2p_t* b, f2p_t* res);

//---------------------------------------------------------
//
// performs modulo operation
//
// a % mod = res
//
// @param a:   f2p_t
// @param mod: modulo
// @param res: result of modulo
//
//---------------------------------------------------------
void f2p_mod(f2p_t* a, fp_t* mod, fp_t* res);


//---------------------------------------------------------
//
// performs integer division of two f2p_t where b = 2^k
//
// a / b = res
//
// @param a:   dividend
// @param b:   divisor
// @param res: result of integer division
//
//---------------------------------------------------------
void f2p_div_2k(f2p_t* a, f2p_t* div, f2p_t* res);


//---------------------------------------------------------
//
// performs modulo operation where mod = 2^k
//
// a % mod = res
//
// @param a:   f2p_t
// @param mod: modulo
// @param res: result of modulo
//
//---------------------------------------------------------
void f2p_mod_2k(f2p_t*a, fp_t* mod, fp_t* res);


void fp_bitwise_and(fp_t* a, fp_t* b, fp_t* res);

void f2p_bitwise_and(f2p_t* a, f2p_t* b, f2p_t* res);


//
//
//  File stuff for testing
//
//

void fp_fprintf(FILE* file, fp_t* p);

void f2p_fprintf(FILE* file, f2p_t* p);

int fp_is_2k(fp_t* a);

void fp_mul_small(fp_t* a, uint64_t b, f2p_t* res);
#endif