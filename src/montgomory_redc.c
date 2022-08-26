

#include"montgomory_redc.h"
#include"fp.h"



void gcd(fp_t* m, fp_t* n, fp_t* m_mul, fp_t* n_mul)
{
    fp_t old_r;
    fp_t r_;
    fp_t old_a;
    fp_t a;
    fp_t old_b;
    fp_t b;
    fp_copy(n, &old_r);
    fp_copy(m, &r_);

    fp_zero(&old_a);
    fp_zero(&a);
    a[0] = 1;


    fp_zero(&old_b);
    fp_zero(&b);
    old_b[0] = 1;

    fp_t zero;
    fp_zero(&zero);

    fp_t quotient;
    fp_t temp;

    f2p_t sub_big;
    fp_t subtractor;

    fp_zero(&subtractor);
    

    while(!fp_equ(&r_, &zero))
    {
        fp_zero(&quotient);
        fp_div(&old_r, &r_, &quotient);
        fp_copy(&old_r, &temp);
        fp_copy(&r_, &old_r);
        fp_mul(&quotient, &r_, &sub_big);
        f2p_to_fp(&sub_big, &subtractor);
        fp_sub(&temp, &subtractor, &r_);

        fp_zero(&temp);
        fp_copy(&old_a, &temp);
        fp_copy(&a, &old_a);
        fp_mul(&quotient, &a, &sub_big);
        f2p_to_fp(&sub_big, &subtractor);
        fp_sub(&temp, &subtractor, &a);


        fp_zero(&temp);
        fp_copy(&old_b, &temp);
        fp_copy(&b, &old_b);
        fp_mul(&quotient, &b, &sub_big);
        f2p_to_fp(&sub_big, &subtractor);
        fp_sub(&temp, &subtractor, &b); 
    }
    fp_copy(&old_a, m_mul);
    fp_copy(&old_b, n_mul);
}

void init(fp_t* mod)
{
    r[0] = 1;
    while(fp_greater_pos(mod, &r))
    {
        fp_lshift(&r, 1);
    }
    f2p_t r2_part;
    fp_mul(&r, &r, &r2_part);
    f2p_mod(&r2_part, mod, &r2);
    fp_t n_strich_pre;
    gcd(&r, mod, &r_minus, &n_strich_pre);
    fp_additive_inverse(&n_strich_pre, &n_strich); 
}

void REDC(fp_t* T, fp_t* mod, fp_t* res)
{
    fp_t tmr;
    f2p_t tmrmn;
    fp_t m;
    fp_mod_2k(T, &r, &tmr);
    fp_mul(&tmr, &n_strich, &tmrmn);
    //f2p_to_fp(&tmrmn, &tmr);
    f2p_mod(&tmrmn, &r, &m);


    f2p_t MN;
    fp_t mn;
    fp_t tmn;
    fp_mul(&m, mod, &MN);
    f2p_to_fp(&MN, &mn);
    fp_add(&mn, T, &tmn);
    fp_t t;


    int shift = get_shift(&r);
    fp_rshift(&t, shift);

    if( fp_greater_equ_pos(&t, mod))
        fp_copy(&t, res);
    else
        fp_sub(&t, mod, res);
}

void REDCL(f2p_t* T, fp_t* mod, fp_t* res)
{
    fp_t tmr;       // T mod R
    f2p_t tmrmn;    // (T mod R) * n_strich
    fp_t m;

    f2p_mod_2k(T, &r, &tmr);              // T mod R     
    fp_mul(&tmr, &n_strich, &tmrmn);       // (T mod R) * n_strich
    f2p_mod_2k(&tmrmn, &r, &m);              // m <- (T mod R) * n_strich mod R

    f2p_t mn;       //m*n
    f2p_t mnpt;     //m*n+T
    f2p_t t;        //m*n+T/R

    f2p_t R;
    f2p_zero(&R);
    for(int i = 0; i < WORDS; i++)
    {
        R[i] = r[i];
    }
    fp_mul(&m, mod, &mn);
    f2p_add(&mn, T, &mnpt);
    f2p_div_2k(&mnpt, &R, &t);

    fp_t t_small;
    for(int i = 0; i < WORDS; i++)
    {
        t_small[i] = t[i];
    }
    if(fp_greater_equ_pos(&t_small, mod))
    {
        fp_sub(&t_small, mod, res);
    }
    else
    {
        fp_copy(&t_small, res);
    }
}

void MODMUL(fp_t* a, fp_t* b, fp_t* mod, fp_t* res)
{
    fp_t a_mon;
    fp_t b_mon;
    f2p_t a_mult;
    f2p_t b_mult;

    fp_mul(a, &r2, &a_mult);
    fp_mul(b, &r2, &b_mult);

    REDCL(&a_mult, mod, &a_mon);
    REDCL(&b_mult, mod, &b_mon);

    f2p_t ab; // a*b in montgomory form
    fp_t ab_2;
    fp_mul(&a_mon, &b_mon, &ab);
    REDCL(&ab, mod, &ab_2);
    f2p_zero(&ab);
    for(int i = 0; i < WORDS; i++)
    {
        ab[i] = ab_2[i];
    }
    REDCL(&ab, mod, res);

}

void ModMul(fp_t* a, fp_t* b, fp_t* mod, fp_t* res)
{
    fp_t a_mon;
    fp_t b_mon;
    f2p_t a_mult;
    f2p_t b_mult;

    fp_mul(a, &r, &a_mult);
    fp_mul(b, &r, &b_mult);

    f2p_mod(&a_mult, mod, &a_mon);
    f2p_mod(&b_mult, mod, &b_mon);

    fp_t mon_pro_1;
    fp_t one;
    fp_zero(&one);
    one[0] = 1;
    MonPro(&a_mon, &b_mon, mod, &mon_pro_1);
    MonPro(&mon_pro_1, &one, mod, res);
}


//helperfunctions
void MonPro(fp_t* a, fp_t* b, fp_t* mod, fp_t* res)
{
    f2p_t ab;
    fp_t t;
    fp_mul(a, b, &ab);
    f2p_mod_2k(&ab, &r, &t);
    f2p_t tn_strich;
    fp_t m;
    fp_mul(&t, &n_strich, &tn_strich);
    f2p_mod_2k(&tn_strich, &r, &m);
    f2p_t mn;
    fp_mul(&m, mod, &mn);
    f2p_t addal;
    f2p_add(&mn, &ab, &addal);
    f2p_t u_l;
    fp_t u;
    f2p_t R;
    f2p_zero(&R);
    for(int i = 0; i < WORDS; i++)
    {
        R[i] = r[i];
    }
    f2p_div_2k(&addal, &R, &u_l);
    f2p_to_fp(&u_l, &u);
    if(fp_greater_equ_pos(&u, mod))
    {
        fp_sub(&u, mod, res);
    }
    else
    {
        fp_copy(&u, res);
    }
}
int get_shift(fp_t* a)
{
    int ctr = 0;
    fp_t a_in;
    fp_copy(a, &a_in);
    while(a_in[0] != 1)
    {
        ctr++;
        fp_rshift(&a_in, 1);
    }
    return ctr;
}
