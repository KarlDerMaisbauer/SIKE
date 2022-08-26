

#include<stdio.h>
#include <stdlib.h>
#include<stdint.h>


#include"fp_helper.h"


//
//
//  Helper functions for the fp_t datatype
//
//

void fp_zero(fp_t* a)
{
    for(int i = 0; i < WORDS; i++)
    {
        (*a)[i] = 0;
    }
}

void fp_print(fp_t* p)
{
    for(int i = WORDS-1; i >= 0; i--)
    {
        printf("%16lx ", (*p)[i]);
    }
}

int fp_smaller_zero(fp_t* a)
{
    return ((int64_t)((*a)[WORDS-1]) < 0) ? 1 : 0;

}

// only for positive numbers
//retrun a >= b
int fp_greater_equ_pos(fp_t* a, fp_t* b)
{
    for(int8_t i = WORDS-1; i >= 0; i--)
    {
        if((*a)[i] > (*b)[i])
            return 1;
        if((*a)[i] < (*b)[i])
            return 0;
    }
    return 1;
}

int fp_equ(fp_t* a, fp_t* b)
{
    for(uint8_t i = 0; i < WORDS; i++)
    {
        if((*a)[i] != (*b)[i])
            return 0;
    }
    return 1;
}

int fp_greater_pos(fp_t* a, fp_t* b)
{
    for(int i = WORDS-1; i >= 0; i--)
    {
        if(!(*a)[i] && !(*b)[i])
            continue;
        if((*a)[i] <= (*b)[i])
            return 0;
    }
    return 1;
}
void fp_copy(fp_t* a, fp_t* b)
{
    for(int i = 0; i < WORDS; i++)
    {
        (*b)[i] = (*a)[i];
    }
}

void fp_lshift(fp_t *a, uint8_t shift)
{
    for(uint8_t i = 0; i < shift; i++)
    {
        uint64_t carry = 0;
        for(uint8_t ctr = 0; ctr < WORDS; ctr++)
        {
            uint64_t new_number = ((*a)[ctr] << 1) + carry;
            carry = (*a)[ctr]>>63;
            (*a)[ctr] = new_number;
        }
    }
}

void fp_rshift(fp_t *a, uint8_t shift)
{
    for(uint8_t i = 0; i < shift; i++)
    {
        uint64_t carry = 0;
        for(int16_t ctr = WORDS-1; ctr >= 0; ctr--)
        {
            uint64_t new_number = ((*a)[ctr] >> 1) + carry;
            carry = (*a)[ctr]<<63;
            (*a)[ctr] = new_number;
        }
    }
}




//
//
//  basic arithmetic functions for the fp_t datatype
//
//
void fp_add(fp_t* a, fp_t* b, fp_t* res)
{
    (*res)[0] = (*a)[0] + (*b)[0];
    uint64_t carry = 0;
    if((uint64_t)(*a)[0] > (uint64_t)(*res)[0] || (uint64_t)(*b)[0] > (uint64_t)(*res)[0])
        carry = 1;
    for(int i = 1; i < WORDS; i++)
    {
        uint64_t temp = (*a)[i] + (*b)[i] + carry;
        if(((uint64_t)(*a)[i]) > ((uint64_t)temp) || ((uint64_t)(*b)[i]) > ((uint64_t)temp) ||
          (carry && (((uint64_t)(*a)[i]) >= ((uint64_t)temp) || ((uint64_t)(*b)[i]) >= ((uint64_t)temp))) )
            carry = 1;
        else
            carry = 0;
        ((*res)[i]) = temp;
    }
}

void fp_additive_inverse(fp_t* a, fp_t* res)
{
    uint64_t inverter = 0xFFFFFFFFFFFFFFFF;
    fp_t addone;
    fp_t res_temp;
    fp_zero(&addone);
    for(int i = 0; i < WORDS; i++)
    {
        res_temp[i] = (*a)[i] ^ inverter;
    }
    addone[0] = 1;
    fp_add(&res_temp, &addone, res);
}

void fp_sub(fp_t* a, fp_t* b, fp_t *res)
{
    fp_t b_inv;
    fp_additive_inverse(b, &b_inv);
    fp_add(a, &b_inv, res);
}

void fp_mul(fp_t* a, fp_t* b, f2p_t* res)
{

    int change_sign = 0;
    fp_t a_inner;
    fp_t b_inner;
    if(fp_smaller_zero(a) && fp_smaller_zero(b))
    {
        fp_additive_inverse(a, &a_inner);
        fp_additive_inverse(b, &b_inner);
    }
    else if(fp_smaller_zero(a))
    {
       fp_additive_inverse(a, &a_inner); 
       fp_copy(b, &b_inner);
       change_sign = 1;
    }
    else if(fp_smaller_zero(b))
    {
       fp_additive_inverse(b, &b_inner); 
       fp_copy(a, &a_inner);
       change_sign = 1;
    }
    else
    {
        fp_copy(a, &a_inner);
        fp_copy(b, &b_inner);
    }
    f2p_zero(res);
    f2p_t pre_res;
    f2p_zero(&pre_res);
    for(uint8_t a_iter = 0; a_iter < WORDS; a_iter++)
    {
       for(uint8_t b_iter = 0; b_iter < WORDS; b_iter++)
        {
            uint8_t shift = a_iter + b_iter;
            f2p_t adder;
            f2p_zero(&adder);
            uint64_t res_up   = ((a_inner[a_iter])>>32) * ((b_inner[b_iter])>>32);
            uint64_t res_down = a_inner[a_iter] * b_inner[b_iter];
            adder[shift] = res_down;
            adder[shift+1] = res_up;
            f2p_add(&pre_res, &adder, &pre_res);
            
        } 
    }
    f2p_copy(&pre_res, res);
    if(change_sign)
        f2p_additive_inverse(res, res);
}

void fp_div(fp_t* a, fp_t* b, fp_t* res)
{
    fp_t part_res;
    fp_t a_inner;
    fp_t b_inner;
    fp_t b_shift;
    fp_zero(&part_res);
    fp_zero(res);
    int change_sign = 0;
    if(fp_smaller_zero(a) && fp_smaller_zero(b))
    {
        fp_additive_inverse(a, &a_inner);
        fp_additive_inverse(b, &b_inner);
    }
    else if(fp_smaller_zero(a))
    {
       fp_additive_inverse(a, &a_inner); 
       fp_copy(b, &b_inner);
       change_sign = 1;
    }
    else if(fp_smaller_zero(b))
    {
       fp_additive_inverse(b, &b_inner); 
       fp_copy(a, &a_inner);
       change_sign = 1;
    }
    else
    {
        fp_copy(a, &a_inner);
        fp_copy(b, &b_inner);
    }
    fp_copy(&b_inner, &b_shift);
    fp_lshift(&b_shift, 1);
    while(fp_greater_equ_pos(&a_inner, &b_inner))
    {
        fp_t part_er;
        fp_zero(&part_er);
        part_er[0] = 1;

        while(fp_greater_equ_pos(&a_inner, &b_shift))
        {
            fp_lshift(&part_er, 1);
            fp_lshift(&b_shift, 1);
            

        }
        fp_rshift(&b_shift, 1);
        fp_t a_sub;
        fp_zero(&a_sub);
        fp_sub(&a_inner, &b_shift, &a_sub);
        fp_copy(&a_sub, &a_inner);
        fp_copy(&b_inner, &b_shift);;
        fp_lshift(&b_shift, 1);
        fp_add(res, &part_er, res);
        fp_zero(&part_er);

    }
    if(change_sign){
        fp_additive_inverse(res, res);}
}

void fp_mod(fp_t* a, fp_t* mod, fp_t* res)
{
    fp_t fraction;
    fp_t multiplicator;
    f2p_t mod_mul_l;
    fp_t mod_mul;
    fp_div(a, mod, &fraction);
    if(fp_smaller_zero(a))
    {
        fp_t one;
        fp_zero(&one);
        one[0] = 1;
        fp_add(&fraction, &one, &multiplicator);
    }
    else
    {
        fp_copy(&fraction, &multiplicator);
    }
    fp_mul(mod, &multiplicator, &mod_mul_l);
    f2p_to_fp(&mod_mul_l, &mod_mul);
    fp_sub(a, &mod_mul, res);
}

void fp_mod_2k(fp_t*a, fp_t* mod, fp_t* res)
{
    int change = 0;
    fp_t a_used;
    fp_t res_part;
    fp_zero(&res_part);
    if(fp_smaller_zero(a))
    {
        fp_additive_inverse(a, &a_used);
        change = 1;
    }
    else
    {
        fp_copy(a, &a_used);
    }
    int ctr = 0;
    while((*mod)[ctr] == 0 && ctr < WORDS)
    {
        res_part[ctr] = a_used[ctr];
        ctr++;
    }
    

    uint64_t numbertocheck = a_used[ctr];
    uint64_t checker = 1;
    uint64_t newnum = 0;
    while(checker != (*mod)[ctr])
    {
        newnum |= checker;
        checker = (checker << 1);
        
    }
    uint64_t mod_num = numbertocheck & newnum;
    res_part[ctr] = mod_num;
    ctr++;

    while(ctr < WORDS)
    {
        res_part[ctr] = 0;
        ctr++;
    }
    if(change)
    {
        fp_sub(mod, &res_part, res);
    }
    else
    {
        fp_copy(&res_part, res);
    }
}


//
// Functions for f2p_t datatype
// (same as fp_t but with doubled bit)
//


//
//
//  Helper functions for the f2p_t datatype
//
//

void f2p_zero(f2p_t* a)
{
    for(int i = 0; i < 2*WORDS; i++)
    {
        (*a)[i] = 0;
    }
}

void f2p_print(f2p_t* p)
{
    for(int i = 2*WORDS-1; i >= 0; i--)
    {
        printf("%16lx ", (*p)[i]);
        if(i % WORDS == 0)
            printf("\n");
    }
}

void f2p_lshift(f2p_t *a, uint8_t shift)
{
    for(uint8_t i = 0; i < shift; i++)
    {
        uint64_t carry = 0;
        for(uint8_t ctr = 0; ctr < 2*WORDS; ctr++)
        {
            uint64_t new_number = ((*a)[ctr] << 1) + carry;
            carry = (*a)[ctr]>>63;
            (*a)[ctr] = new_number;
        }
    }
}

void f2p_to_fp(f2p_t* a, fp_t* b)
{
    for(uint8_t i = 0; i < WORDS; i++)
    {
        (*b)[i] = (*a)[i];
    }
}

int f2p_smaller_zero(f2p_t* a)
{
    return ((int64_t)((*a)[2*WORDS-1]) < 0) ? 1 : 0;

}

int f2p_greater_equ_pos(f2p_t* a, f2p_t* b)
{
    for(int8_t i = 2*WORDS-1; i >= 0; i--)
    {
        if((*a)[i] > (*b)[i])
            return 1;
        if((*a)[i] < (*b)[i])
            return 0;
    }
    return 1;
}

int f2p_equ(f2p_t*a, f2p_t*b)
{
    for(int i = 0; i < 2*WORDS;i++)
    {
        if((*a)[i] != (*b)[i])
            return 0;
    }
    return 1;
}

void f2p_copy(f2p_t* a, f2p_t* b)
{
    for(int i = 0; i < 2*WORDS; i++)
    {
        (*b)[i] = (*a)[i];
    }
}

void f2p_rshift(f2p_t *a, uint8_t shift)
{
    for(uint8_t i = 0; i < shift; i++)
    {
        uint64_t carry = 0;
        for(int16_t ctr = 2*WORDS-1; ctr >= 0; ctr--)
        {
            uint64_t new_number = ((*a)[ctr] >> 1) + carry;
            carry = (*a)[ctr]<<63;
            (*a)[ctr] = new_number;
        }
    }
}


//
//
//  basic arithmetic functions for the fp_t datatype
//
//


void f2p_add(f2p_t* a, f2p_t* b, f2p_t* res)
{
    (*res)[0] = (*a)[0] + (*b)[0];
    uint64_t carry = 0;
    if((uint64_t)(*a)[0] > (uint64_t)(*res)[0] || (uint64_t)(*b)[0] > (uint64_t)(*res)[0])
        carry = 1;
    for(int i = 1; i < 2*WORDS; i++)
    {
        uint64_t temp = (*a)[i] + (*b)[i] + carry;
        if(((uint64_t)(*a)[i]) > ((uint64_t)temp) || ((uint64_t)(*b)[i]) > ((uint64_t)temp) ||
          (carry && (((uint64_t)(*a)[i]) >= ((uint64_t)temp) || ((uint64_t)(*b)[i]) >= ((uint64_t)temp))) )
            carry = 1;
        else
            carry = 0;
        ((*res)[i]) = temp;
    }
}

void f2p_additive_inverse(f2p_t* a, f2p_t* res)
{
    uint64_t inverter = 0xFFFFFFFFFFFFFFFF;
    f2p_t addone;
    f2p_t res_temp;
    f2p_zero(&addone);
    for(int i = 0; i < 2*WORDS; i++)
    {
        res_temp[i] = (*a)[i] ^ inverter;
    }
    addone[0] = 1;
    f2p_add(&res_temp, &addone, res);
}

void f2p_sub(f2p_t*a, f2p_t*b, f2p_t*res)
{
    f2p_t inv;
    f2p_additive_inverse(b, &inv);
    f2p_add(a, &inv, res);
}

void f2p_mul(f2p_t* a, f2p_t* b, f2p_t* res)
{

    int change_sign = 0;
    f2p_t a_inner;
    f2p_t b_inner;
    if(f2p_smaller_zero(a) && f2p_smaller_zero(b))
    {
        f2p_additive_inverse(a, &a_inner);
        f2p_additive_inverse(b, &b_inner);
    }
    else if(f2p_smaller_zero(a))
    {
       f2p_additive_inverse(a, &a_inner); 
       f2p_copy(b, &b_inner);
       change_sign = 1;
    }
    else if(f2p_smaller_zero(b))
    {
       f2p_additive_inverse(b, &b_inner); 
       f2p_copy(a, &a_inner);
       change_sign = 1;
    }
    else
    {
        f2p_copy(a, &a_inner);
        f2p_copy(b, &b_inner);
    }
    f2p_zero(res);
    for(uint8_t a_iter = 0; a_iter < 2*WORDS; a_iter++)
    {
       for(uint8_t b_iter = 0; b_iter < 2*WORDS; b_iter++)
        {
            uint8_t shift = a_iter + b_iter;
            if(shift >= 2*WORDS)
                continue;
            f2p_t adder;
            f2p_zero(&adder);
            uint64_t res_up   = ((a_inner[a_iter])>>32) * ((b_inner[b_iter])>>32);
            uint64_t res_down = a_inner[a_iter] * b_inner[b_iter];
            adder[shift] = res_down;
            if(shift+1 < 2*WORDS)
                adder[shift+1] = res_up;
            f2p_add(res, &adder, res);
            
        } 
    }
    if(change_sign)
        f2p_additive_inverse(res, res);
}

void f2p_div(f2p_t* a, f2p_t* b, f2p_t* res)
{
    f2p_t part_res;
    f2p_t a_inner;
    f2p_t b_inner;
    f2p_t b_shift;
    f2p_zero(&part_res);
    f2p_zero(res);
    int change_sign = 0;
    if(f2p_smaller_zero(a) && f2p_smaller_zero(b))
    {
        f2p_additive_inverse(a, &a_inner);
        f2p_additive_inverse(b, &b_inner);
    }
    else if(f2p_smaller_zero(a))
    {
       f2p_additive_inverse(a, &a_inner); 
       f2p_copy(b, &b_inner);
       change_sign = 1;
    }
    else if(f2p_smaller_zero(b))
    {
       f2p_additive_inverse(b, &b_inner); 
       f2p_copy(a, &a_inner);
       change_sign = 1;
    }
    else
    {
        f2p_copy(a, &a_inner);
        f2p_copy(b, &b_inner);
    }
    f2p_copy(&b_inner, &b_shift);
    f2p_lshift(&b_shift, 1);
    while(f2p_greater_equ_pos(&a_inner, &b_inner))
    {
        f2p_t part_er;
        f2p_zero(&part_er);
        part_er[0] = 1;
        while(f2p_greater_equ_pos(&a_inner, &b_shift))
        {
            f2p_lshift(&part_er, 1);
            f2p_lshift(&b_shift, 1);
        }
        f2p_rshift(&b_shift, 1);
        f2p_t a_sub;
        f2p_zero(&a_sub);
        f2p_sub(&a_inner, &b_shift, &a_sub);
        f2p_copy(&a_sub, &a_inner);
        f2p_copy(&b_inner, &b_shift);;
        f2p_lshift(&b_shift, 1);
        f2p_add(res, &part_er, res);
        f2p_zero(&part_er);
    }
    if(change_sign){
        f2p_additive_inverse(res, res);}
}

void f2p_mod(f2p_t* a, fp_t* mod, fp_t* res)
{
    f2p_t MOD;
    f2p_zero(&MOD);
    for(int i = 0; i < WORDS; i++)
    {
        MOD[i] = (*mod)[i];
    }
    f2p_t fraction;
    f2p_t multiplicator;
    f2p_t mod_mul_l;
    f2p_div(a, &MOD, &fraction);
    if(f2p_smaller_zero(&fraction))
    {
        f2p_t one;
        f2p_zero(&one);
        one[0] = 1;
        f2p_sub(&fraction, &one, &multiplicator);
    }
    else
    {
        f2p_copy(&fraction, &multiplicator);
    }
    f2p_mul(&MOD, &multiplicator, &mod_mul_l);
    


    f2p_t l_res;
    f2p_sub(a, &mod_mul_l, &l_res);
    f2p_to_fp(&l_res, res);
}


void f2p_div_2k(f2p_t* a, f2p_t* div, f2p_t* res)
{
    f2p_t a_inner;
    f2p_t div_inner;
    f2p_copy(a, &a_inner);
    f2p_copy(div, &div_inner);
    while(div_inner[0] != 1)
    {
        f2p_rshift(&a_inner, 1);
        f2p_rshift(&div_inner, 1);
    }
    f2p_copy(&a_inner, res);
}

void f2p_mod_2k(f2p_t*a, fp_t* mod, fp_t* res)
{
    int change = 0;
    f2p_t a_used;
    f2p_t res_part;
    f2p_t big_mod;
    f2p_zero(&big_mod);
    for(int i = 0; i < WORDS; i++)
    {
        big_mod[i] = (*mod)[i];
    }
    f2p_zero(&res_part);
    if(f2p_smaller_zero(a))
    {
        f2p_additive_inverse(a, &a_used);
        change = 1;
    }
    else
    {
        f2p_copy(a, &a_used);
    }
    int ctr = 0;
    while(big_mod[ctr] == 0 && ctr < 2*WORDS)
    {
        res_part[ctr] = a_used[ctr];
        ctr++;
    }


    uint64_t numbertocheck = a_used[ctr];
    uint64_t checker = 1;
    uint64_t newnum = 0;
    while(checker != big_mod[ctr])
    {
        newnum = newnum | checker;
        checker = (checker << 1);
        
    }
    uint64_t mod_num = numbertocheck & newnum;
    res_part[ctr] = mod_num;
    ctr++;

    if(change)
    {
        f2p_t sub;
        f2p_sub(&big_mod, &res_part, &sub);
        f2p_to_fp(&sub, res);
    }
    else
    {
        f2p_to_fp(&res_part, res);
    }
}


































