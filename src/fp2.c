
#include"fp2.h"
#include"fp.h"
#include"fp_helper.h"
#include"montgomory_redc.h"


void fp2_addm(fp2_t* a, fp2_t*b, fp_t* mod, fp2_t* res)
{
    fp_addm(&a->real, &b->real, mod, &res->real);
    fp_addm(&a->img, &b->img, mod, &res->img);
}

void fp2_subm(fp2_t* a, fp2_t*b, fp_t* mod, fp2_t* res)
{
    fp_subm(&a->real, &b->real, mod, &res->real);
    fp_subm(&a->img, &b->img, mod, &res->img);
}


void fp2_mulm(fp2_t* a, fp2_t* b, fp_t*mod, fp2_t* res)
{
    fp_t temp_1;
    fp_t temp_2;
    fp_mulm(&a->real, &b->real, mod, &temp_1);
    fp_mulm(&a->img, &b->img, mod, &temp_2);
    fp_subm(&temp_1, &temp_2, mod, &res->real);

    fp_mulm(&a->real, &b->img, mod, &temp_1);
    fp_mulm(&a->img, &b->real, mod, &temp_2);
    fp_addm(&temp_1, &temp_2, mod, &res->real);
}

void fp2_add_inv(fp2_t* a, fp_t*mod, fp2_t* res)
{
    fp_t temp_inv;
    fp_additive_inverse(&a->real, &temp_inv);
    fp_mod(&temp_inv, mod, &res->real);

    fp_additive_inverse(&a->img, &temp_inv);
    fp_mod(&temp_inv, mod, &res->img);
}

void fp2_mult_inv(fp2_t* a, fp_t*mod, fp2_t* res)
{
    fp_t a0q;           //a0^2
    fp_t a1q;           //a1^2
    fp_t a0qpa1q;       //a0^2+a1^2
    fp_t a0qpa1q_inv;   //(a0^2+a1^2)^(-1)
    fp_t mod_inv;
    fp_t a1_inv;        //-a1
    fp_mulm(&a->real, &a->real, mod, &a0q);
    fp_mulm(&a->img, &a->img, mod, &a1q);
    fp_addm(&a0q, &a1q, mod, &a0qpa1q);
    gcd(&a0qpa1q, mod, &a0qpa1q_inv, &mod_inv);
    fp_mulm(&a0q, &a0qpa1q_inv, mod, &res->real);
    fp_additive_inverse(&a->img, &a1_inv);
    fp_mulm(&a1_inv, &a0qpa1q_inv, mod, &res->img);
}

void fp2_divm(fp2_t* a, fp2_t* b, fp_t*mod, fp2_t* res)
{
    fp2_t b_inv;
    fp2_mult_inv(b, mod, &b_inv);
    fp2_mulm(a, &b_inv, mod, res);
}

void fp2_add(fp2_t* a, fp2_t*b, fp2_t* res)
{
    fp_add(&a->real, &b->real, &res->real);
    fp_add(&a->img,  &b->img,  &res->img);
}

void fp2_sub(fp2_t* a, fp2_t* b, fp2_t *res)
{
    fp_sub(&a->real, &b->real, &res->real);
    fp_sub(&a->img,  &b->img,  &res->img);
}

void fp2_mul(fp2_t* a, fp2_t* b, f2p2_t* res)
{
    f2p_t temp_1;
    f2p_t temp_2;
    fp_mul(&a->real, &b->real, &temp_1);
    fp_mul(&a->img,  &b->img,  &temp_2);
    f2p_sub(&temp_1, &temp_2, &res->real);

    fp_mul(&a->real, &b->img, &temp_1);
    fp_mul(&a->img, &b->real, &temp_2);
    f2p_add(&temp_1, &temp_2, &res->real);
}

void fp2_mul_mont(fp2_t* a, fp2_t* b, fp_t* mod, fp2_t* res)
{
    f2p2_t temp;
    fp2_mul(a, b, &temp);
    REDCL2(&temp, mod, res);
}


void fp2_zero(fp2_t* a)
{
    fp_zero(&a->real);
    fp_zero(&a->img);
}

void fp2_copy(fp2_t* a, fp2_t* a_copy)
{
    fp_copy(&a->real, &a_copy->real);
    fp_copy(&a->img, &a_copy->img);
}  



void fp2_copy_masked(fp2_t* a, fp2_t* b, uint64_t mask)
{
    fp_copy_masked(&a->real, &b->real, mask);
    fp_copy_masked(&a->img, &b->img, mask);
}


int fp2_greater_equ_pos(fp2_t* a, fp2_t* b)
{
    f2p_t a_real_sq;
    f2p_t a_img_sq;

    f2p_t b_real_sq;
    f2p_t b_img_sq;

    f2p_t size_a;
    f2p_t size_b;

    fp_mul(&a->real, &a->real, &a_real_sq);
    fp_mul(&a->img, &a->img, &a_img_sq);

    fp_mul(&b->real, &b->real, &a_real_sq);
    fp_mul(&b->img, &b->img, &a_img_sq);

    f2p_add(&a_real_sq, &a_img_sq, &size_a);
    f2p_add(&b_real_sq, &b_img_sq, &size_b);

    return f2p_greater_equ_pos(&size_a, &size_b);
}