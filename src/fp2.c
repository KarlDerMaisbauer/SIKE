
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