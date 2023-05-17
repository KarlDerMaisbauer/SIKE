

#include<stdio.h>
#include <stdlib.h>
#include<stdint.h>


#include "fp_helper.h"
#include "stdio.h"


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
        printf("%16lx", (*p)[i]);
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
    //printf("\n\n");
    for(int i = WORDS-1; i >= 0; i--)
    {
        //if(!(*a)[i] && !(*b)[i])
        //    continue;
            /*
        if((*a)[i] <= (*b)[i])
        {
            printf("%lx <= %lx\n",(*a)[i], (*b)[i]);
            return 0;
        }
        */
        if((uint64_t)(*a)[i] < (uint64_t)(*b)[i])
            return 0;
        if((uint64_t)(*a)[i] > (uint64_t)(*b)[i])
            return 1;
    }
    return 0;
}
void fp_copy(fp_t* a, fp_t* b)
{
    for(int i = 0; i < WORDS; i++)
    {
        (*b)[i] = (*a)[i];
    }
}

void fp_copy_masked(fp_t* a, fp_t* b, uint64_t mask)
{
    for(int i = 0; i < WORDS; i++)
    {
        (*b)[i] = ((*b)[i] & (~mask)) | ((*a)[i] & mask);
    }
}


void fp_lshift(fp_t *a, uint8_t shift)
{
    for(uint8_t i = 0; i < shift; i++)
    {
        uint64_t carry = 0;
        for(uint8_t ctr = 0; ctr < WORDS; ctr++)
        {
            uint64_t new_number = ((*a)[ctr] << 1) | carry;
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
            uint64_t new_number = ((*a)[ctr] >> 1) | carry;
            carry = (*a)[ctr]<<63;
            (*a)[ctr] = new_number;
        }
    }
}

int64_t fp_get_len(fp_t* a)
{
    int64_t len = 0;
    int break_point = 0;
    for(int i = WORDS - 1; i >= 0; i--)
    {
        if((*a)[i] == 0)
            continue;
        else
        {
            uint64_t num = (*a)[i];
            for(int pos = 63; pos >= 0; pos--)
            {
                uint64_t mask = (uint64_t)1 << pos;
                if(num & mask)
                {
                    len = pos + 1 + i * 64;
                    break_point = 1;
                    break;
                }
                else
                    continue;
            }
        }
        if(break_point)
            break;
    }
    return len;
}


//
//
//  basic arithmetic functions for the fp_t datatype
//
//
void fp_add(fp_t* a, fp_t* b, fp_t* res)
{
    /*
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
    }*/
    uint64_t carry;
    (*res)[0] = (*a)[0] + (*b)[0];
    carry = ((*a)[0] > ((uint64_t)-1) - (*b)[0]) ? 1 : 0;

    for(int i = 1; i < WORDS; i++)
    {
        (*res)[i] = (*a)[i] + (*b)[i] + carry;
        uint64_t test_carry = carry;

        carry = ((*a)[i] > ((uint64_t)-1) - (*b)[i]) ? 1 : 0;
        carry = ((((*a)[i]) + (*b)[i]) > (((uint64_t)-1) - test_carry)) ? 1 : carry;
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

void fp_mul_small(fp_t* a, uint64_t b, f2p_t* res)
{
    f2p_t res_pre;
    f2p_zero(&res_pre);
    uint64_t b_lower  = 0x00000000FFFFFFFF & b;
    uint64_t b_higher = b >> 32; 


    
    for(int a_iter = 0; a_iter < WORDS; a_iter++)
    {
        /*
        printf("a from small mul\n");
        fp_print(a);
        LINE 
        printf("i = %d\n", a_iter);*/
        f2p_t adder;
        f2p_zero(&adder);
        
        uint64_t a_part = (*a)[a_iter];
        //printf("val of a at %d, is %lx\n", a_iter, a_part);
        uint64_t a_lower  = 0x00000000FFFFFFFF & a_part;
        uint64_t a_higher = a_part >> 32;
        //printf("a_lower = %lx\n", a_lower);
        //printf("a_highe = %lx\n", a_higher);
        

        uint64_t res_middle1 = a_lower  * b_higher;
        uint64_t res_middle2 = a_higher * b_lower;

        uint64_t res_lower = a_lower * b_lower;
        uint64_t carry = 0;

        if((res_middle1 << 32) > (((uint64_t) -1) - (res_middle2 << 32)))
        {
            carry++;
        }


        if(((res_middle1 << 32) + (res_middle2 << 32)) > (((uint64_t) -1) - res_lower))
        {
            carry++;
        }
        res_lower += res_middle1 << 32;
        res_lower += res_middle2 << 32;
        uint64_t res_higher  = a_higher * b_higher;
        res_higher += carry;
        res_higher += res_middle1 >> 32;
        res_higher += res_middle2 >> 32;

        adder[a_iter] = res_lower;
        adder[a_iter+1] = res_higher;
        /*printf("adder:\n");
        f2p_print(&adder);
        LINE LINE LINE*/
        f2p_add(&adder, &res_pre, res);
        f2p_copy(res, &res_pre);


    }
    //LINE LINE LINE LINE LINE LINE LINE LINE LINE
}


void fp_karatsuba(fp_t* a, int a_size, fp_t* b, int b_size, f2p_t* res)
{
    //printf("a_size = %d\n", a_size);
    //printf("b_size = %d\n", b_size);
    /*
    if(a_size < 32 ||  b_size < 32)
    {
        *(res[0]) = *(a[0]) * *(b[0]);
        printf("reached condition\n");
        f2p_print(res);
        LINE LINE LINE LINE
        return; 
    }*/
    f2p_zero(res);
    if(a_size < 32)
    {
        fp_mul_small(b, (*a)[0], res);
        return;
    }
    if(b_size < 32)
    {
        fp_mul_small(a, (*b)[0], res);
        return;
    }

    int max_size = a_size >= b_size ? a_size : b_size;
    int half_size = max_size / 2;
    int shift = half_size * 2;
    //printf("maxsize %d\n", max_size);
    //printf("size %d\n", half_size);
    //printf("shift %d\n", shift);

    fp_t a_low;
    fp_t a_high;

    fp_t b_low;
    fp_t b_high;

  //  fp_zero(&a_low);
  //  fp_zero(&a_high);

  //  fp_zero(&b_low);
  //  fp_zero(&b_high);

    int size_words = half_size / 64;
    //printf("half_size is %d\n", half_size);
    //printf("size_words is %d\n", size_words);
    int remaining_bits = half_size - (size_words * 64);
    //printf("remaining_bits is %d\n", remaining_bits);
    uint64_t mask = 0;
    for(int i = 0; i < remaining_bits; i++)
    {
        mask = mask << 1;
        mask = mask | (uint64_t)1;
    }
    //printf("mask is %ld\n", mask);
    for(int i = 0; i < size_words; i++)
    {
        a_low[i] = (*a)[i];
        b_low[i] = (*b)[i];
        //printf("a_low at %d is %ld\n", i, a_low[i]);
    }
    a_low[size_words] = (*a)[size_words] & mask;
    b_low[size_words] = (*b)[size_words] & mask;

    for(int i = size_words + 1; i < WORDS; i++)
    {
        a_low[i] = 0;
        b_low[i] = 0;
    }
    /*
    printf("the duallity of a\n");
    fp_print(a);
    LINE
    fp_print(&a_low);
    LINE
    fp_print(&a_high);
    LINE*/
    fp_t temp;
    fp_copy(a, &temp);
    fp_rshift(&temp, half_size);
    fp_copy(&temp, &a_high);

    fp_copy(b, &temp);
    fp_rshift(&temp, half_size);
    fp_copy(&temp, &b_high);

    f2p_t z0;
    f2p_t z1;
    f2p_t z2;

    f2p_zero(&z0);
    f2p_zero(&z1);
    f2p_zero(&z2);

    fp_karatsuba(&a_low, half_size, &b_low, half_size, &z0);

    fp_t z1_input1;
    fp_zero(&z1_input1);
    fp_add(&a_low, &a_high, &z1_input1);

    fp_t z1_input2;
    fp_zero(&z1_input2);
    fp_add(&b_low, &b_high, &z1_input2);

    fp_karatsuba(&z1_input1, fp_get_len(&z1_input1), &z1_input2, fp_get_len(&z1_input2), &z1);

    fp_karatsuba(&a_high, fp_get_len(&a_high), &b_high, fp_get_len(&b_high), &z2);
    

    f2p_t res_temp1;
    f2p_t res_temp2;

    // 1. part of result
    f2p_copy(&z2, &res_temp1);
    f2p_lshift(&res_temp1, shift);


    // 2. part of result
    f2p_zero(&res_temp2);

    f2p_t subtractor;
    f2p_zero(&subtractor);
    f2p_add(&z2, &z0, &subtractor);
    f2p_sub(&z1, &subtractor, &res_temp2);
    f2p_lshift(&res_temp2, half_size);

    f2p_t res_temp3;
    f2p_zero(&res_temp3);
    f2p_add(&res_temp1, &res_temp2, &res_temp3);

    f2p_add(&z0, &res_temp3, res);

    /*printf("intermediate res\n");
    f2p_print(res);
    LINE LINE*/



    
    
    
    /*
    f2p_t res_temp;
    f2p_copy(&z2, res);
    f2p_lshift(res, 2*half_size);
    printf("res 1. part\n");
    f2p_print(res);
    LINE
    f2p_t mid_erg;
    f2p_t mid_erg2;
    f2p_zero(&mid_erg);
    f2p_zero(&mid_erg2);
    f2p_sub(&z1, &z2, &mid_erg);
    f2p_sub(&mid_erg, &z0, &mid_erg2);

    f2p_add(res, &mid_erg2, &res_temp);
    f2p_lshift(&res_temp, half_size);


    printf("res 2. part\n");
    f2p_print(&res_temp);
    LINE

    printf("res 3. part\n");
    f2p_add(&res_temp, &z0, res);
    f2p_print(res);
    LINE
    LINE LINE LINE*/
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
            f2p_t pre_res_2;
            f2p_zero(&pre_res_2);
            f2p_copy(&pre_res, &pre_res_2);
            f2p_zero(&adder);


            uint64_t a_lower  = 0x00000000FFFFFFFF & a_inner[a_iter];
            uint64_t a_higher = a_inner[a_iter] >> 32;
            uint64_t b_lower  = 0x00000000FFFFFFFF & b_inner[b_iter];
            uint64_t b_higher = b_inner[b_iter] >> 32;


            uint64_t res_middle1 = a_lower  * b_higher;
            uint64_t res_middle2 = a_higher * b_lower;


            
            uint64_t res_lower = a_lower * b_lower;

            uint64_t carry = 0;

            if((res_middle1 << 32) > (((uint64_t) -1) - (res_middle2 << 32)))
            {
                carry++;
            }


            if(((res_middle1 << 32) + (res_middle2 << 32)) > (((uint64_t) -1) - res_lower))
            {
                carry++;
            }

            res_lower += res_middle1 << 32;
            res_lower += res_middle2 << 32;

            uint64_t res_higher  = a_higher * b_higher;
            res_higher += carry;
            res_higher += res_middle1 >> 32;
            res_higher += res_middle2 >> 32;

            adder[shift] = res_lower;
            adder[shift+1] = res_higher;


            f2p_add(&pre_res_2, &adder, &pre_res);
            
        } 
    }

    //fp_karatsuba(&a_inner, fp_get_len(&a_inner), &b_inner, fp_get_len(&b_inner), &pre_res);

    //f2p_copy(&pre_res, res);
    if(change_sign)
        f2p_additive_inverse(&pre_res, res);
    else
       f2p_copy(&pre_res, res); 
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
    //printf("b_inner:  ");
    //fp_print(&b_inner);
    //LINE
    if(fp_is_2k(&b_inner) && b_inner[0] != 1)
    {
        //printf("\n\n 2k \n\n");
        fp_t b_shift;
        fp_zero(&b_shift);
        b_shift[0] = 1;
        while(fp_greater_equ_pos(&a_inner, &b_inner))
        {
            fp_lshift(&b_inner, 1);
            fp_lshift(&b_shift, 1);
        }
        fp_rshift(&b_inner, 1);
        fp_rshift(&b_shift, 1);
        //fp_sub(&a_inner, &b_inner, res);
        fp_copy(&b_shift, res);
    }
    else
    {
        //printf("\n\n no 2k \n\n");
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
            fp_copy(&b_inner, &b_shift);
            fp_lshift(&b_shift, 1);
            fp_add(res, &part_er, res);
            fp_zero(&part_er);

        }
    }
    fp_t res_pre;
    fp_zero(&res_pre);
    fp_copy(res, &res_pre);
    if(change_sign)
        fp_additive_inverse(&res_pre, res);
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
        fp_sub(&fraction, &one, &multiplicator);
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
    fp_t one;
    fp_zero(&one);
    one[0] = 1;
    fp_t new_mod;
    fp_sub(mod, &one, &new_mod);
    fp_bitwise_and(a, &new_mod, res);
    /*
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
    newnum |= checker;
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
        //fp_sub(mod, &res_part, res);
        fp_t sub;
        fp_additive_inverse(&res_part, &sub);
        fp_add(mod, &sub, res);
    }
    else
    {
        fp_copy(&res_part, res);
    }*/
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
    for(int i = 0; i < DWORDS; i++)
    {
        (*a)[i] = 0;
    }
}

void f2p_print(f2p_t* p)
{
    for(int i = DWORDS-1; i >= 0; i--)
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
        for(uint8_t ctr = 0; ctr < DWORDS; ctr++)
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

void fp_to_f2p(fp_t* a, f2p_t* b)
{
    //f2p_zero(b);
    for(uint8_t i = 0; i < WORDS; i++)
    {
        (*b)[i] = (*a)[i];
    }
    uint64_t mask = (uint64_t)1 << 63;
    //fp_print(a);
    //LINE
    //printf("mask is: %lx\n", mask);
    //printf("last is: %lx\n",  (*a)[WORDS-1]);
    //LINE
    if(mask & (*a)[WORDS-1])
    {
        //printf("mask is true\n");
        for(uint8_t i = WORDS; i < DWORDS; i++)
        {
            (*b)[i] = 0xffffffffffffffffULL;
        }
    }
    else
    {
        for(uint8_t i = WORDS; i < DWORDS; i++)
        {
            (*b)[i] = 0;
        }
    }
    //f2p_print(b);
    //LINE
    //LINE
    //LINE
    //LINE
}

int f2p_smaller_zero(f2p_t* a)
{
    return ((int64_t)((*a)[DWORDS-1]) < 0) ? 1 : 0;

}

int f2p_greater_equ_pos(f2p_t* a, f2p_t* b)
{
    for(int8_t i = DWORDS-1; i >= 0; i--)
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
    for(int i = 0; i < DWORDS;i++)
    {
        if((*a)[i] != (*b)[i])
            return 0;
    }
    return 1;
}

void f2p_copy(f2p_t* a, f2p_t* b)
{
    for(int i = 0; i < DWORDS; i++)
    {
        (*b)[i] = (*a)[i];
    }
}

void f2p_rshift(f2p_t *a, uint8_t shift)
{
    for(uint8_t i = 0; i < shift; i++)
    {
        uint64_t carry = 0;
        for(int16_t ctr = DWORDS-1; ctr >= 0; ctr--)
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
    f2p_zero(res);
    uint64_t carry;
    (*res)[0] = (*a)[0] + (*b)[0];
    carry = ((*a)[0] > ((uint64_t)-1) - (*b)[0]) ? 1 : 0;
    for(int i = 1; i < DWORDS; i++)
    {
        (*res)[i] = (*a)[i] + (*b)[i] + carry;
        uint64_t test_carry = carry;

        carry = ((*a)[i] > ((uint64_t)-1) - (*b)[i]) ? 1 : 0;
        carry = ((((*a)[i]) + (*b)[i]) > (((uint64_t)-1) - test_carry)) ? 1 : carry;
    }
}

void f2p_additive_inverse(f2p_t* a, f2p_t* res)
{
    uint64_t inverter = 0xFFFFFFFFFFFFFFFF;
    f2p_t addone;
    f2p_t res_temp;
    //f2p_zero(res);
    f2p_zero(&res_temp);
    f2p_zero(&addone);
    for(int i = 0; i < DWORDS; i++)
    {
        res_temp[i] = (*a)[i] ^ inverter;
    }
    addone[0] = 1;
    f2p_add(&res_temp, &addone, res);
}

void f2p_addm(f2p_t* a, f2p_t* b, fp_t* mod, fp_t* res)
{
    f2p_t temp;
    f2p_zero(&temp);
    f2p_add(a, b, &temp);
    f2p_mod(&temp, mod, res);
}


void f2p_sub(f2p_t*a, f2p_t*b, f2p_t*res)
{
    f2p_t inv;
    f2p_zero(&inv);
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
    f2p_t pre_res;
    f2p_zero(&pre_res);
    for(uint8_t a_iter = 0; a_iter < DWORDS; a_iter++)
    {
       for(uint8_t b_iter = 0; b_iter < DWORDS; b_iter++)
        {
            uint8_t shift = a_iter + b_iter;
            if(shift >= DWORDS)
                continue;
            f2p_t adder;
            f2p_t pre_res_2;
            f2p_zero(&pre_res_2);
            f2p_copy(&pre_res, &pre_res_2);
            f2p_zero(&adder);


            uint64_t a_lower  = 0x00000000FFFFFFFF & a_inner[a_iter];
            uint64_t a_higher = a_inner[a_iter] >> 32;
            uint64_t b_lower  = 0x00000000FFFFFFFF & b_inner[b_iter];
            uint64_t b_higher = b_inner[b_iter] >> 32;


            uint64_t res_middle1 = a_lower  * b_higher;
            uint64_t res_middle2 = a_higher * b_lower;


            
            uint64_t res_lower = a_lower * b_lower;

            uint64_t carry = 0;

            if((res_middle1 << 32) > (((uint64_t) -1) - (res_middle2 << 32)))
            {
                carry++;
            }


            if(((res_middle1 << 32) + (res_middle2 << 32)) > (((uint64_t) -1) - res_lower))
            {
                carry++;
            }

            res_lower += res_middle1 << 32;
            res_lower += res_middle2 << 32;

            uint64_t res_higher  = a_higher * b_higher;
            res_higher += carry;
            res_higher += res_middle1 >> 32;
            res_higher += res_middle2 >> 32;

            adder[shift] = res_lower;
            if(shift+1 < DWORDS)
                adder[shift+1] = res_higher;


            f2p_add(&pre_res_2, &adder, &pre_res);
            
        } 
    }
    if(change_sign)
        f2p_additive_inverse(&pre_res, res);
    else
        f2p_copy(&pre_res, res);
}


void f2p_mulm(f2p_t* a, f2p_t* b, fp_t* mod, fp_t* res)
{
    f2p_t temp;
    f2p_zero(&temp);
    f2p_mul(a, b, &temp);
    f2p_mod(&temp, mod, res);
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
        f2p_t res_copy;
        f2p_copy(res, &res_copy);
        f2p_add(&res_copy, &part_er, res);
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
    f2p_t zero;
    f2p_t multiplicator;
    f2p_t mod_mul_l;
    f2p_zero(&zero);
    f2p_div(a, &MOD, &fraction);
    if(f2p_smaller_zero(&fraction))
    {
        f2p_t one;
        f2p_zero(&one);
        one[0] = 1;
        f2p_sub(&fraction, &one, &multiplicator);
    }
    else if(f2p_equ(&fraction, &zero) && f2p_smaller_zero(a))
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


    int change_sign = 0;
    if(f2p_smaller_zero(a))
    {
        f2p_additive_inverse(a, &a_inner);
        change_sign = 1;
    }
    else
    {
        f2p_copy(a, &a_inner);

    }
    f2p_copy(div, &div_inner);
    while(div_inner[0] != 1)
    {
        f2p_rshift(&a_inner, 1);
        f2p_rshift(&div_inner, 1);
    }
    f2p_copy(&a_inner, res);
    if(change_sign)
        f2p_additive_inverse(res,res);
}

void f2p_mod_2k(f2p_t*a, fp_t* mod, fp_t* res)
{
    f2p_t big_mod;
    f2p_t one;
    f2p_t res_temp;
    f2p_zero(&big_mod);
    f2p_zero(&one);
    fp_to_f2p(mod, &big_mod);
    one[0] = 1;
    f2p_t new_mod;
    f2p_sub(&big_mod, &one, &new_mod);
    f2p_bitwise_and(a, &new_mod, &res_temp);
    f2p_to_fp(&res_temp, res);
/*
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
    while(big_mod[ctr] == 0 && ctr < DWORDS)
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
        f2p_t temp;
        //f2p_sub(&big_mod, &res_part, &sub);
        f2p_additive_inverse(&res_part, &sub);
        f2p_add(&big_mod, &sub, &temp);
        f2p_to_fp(&temp, res);
    }
    else
    {
        f2p_to_fp(&res_part, res);
    }*/
}

void fp_bitwise_and(fp_t* a, fp_t* b, fp_t* res)
{
    for(int i = 0; i < WORDS; i++)
    {
        (*res)[i] = (*a)[i] & (*b)[i];
    }
}


void f2p_bitwise_and(f2p_t* a, f2p_t* b, f2p_t* res)
{
    for(int i = 0; i < DWORDS; i++)
    {
        (*res)[i] = (*a)[i] & (*b)[i];
    }
}

//
//
//  File stuff for testing
//
//




void fp_fprintf(FILE* file, fp_t* p)
{
    fprintf(file, "0x");
    for(int i = WORDS-1; i >= 0; i--)
    {
        fprintf(file, "%016lx", (*p)[i]);
    }
    fprintf(file, "\n");
}


void f2p_fprintf(FILE* file, f2p_t* p)
{
    fprintf(file, "0x");
    for(int i = DWORDS-1; i >= 0; i--)
    {
        fprintf(file, "%016lx", (*p)[i]);
    }
    fprintf(file, "\n");
}


int fp_is_2k(fp_t* a)
{
    int bitcounter = 0;
    for(int i = 0; i < WORDS; i++)
    {
        uint64_t mask = 1;
        uint64_t num = (*a)[i];
        for(int pos = 0; pos < 64; pos++)
        {
            bitcounter += 1 && (mask & num);
            mask = mask << 1;
        }
    }
    //printf("num bits is %d\n", bitcounter);
    return bitcounter == 1;
}