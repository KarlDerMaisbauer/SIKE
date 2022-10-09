

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

    
    fp2_t one;
    f2p2_t one2;
    fp2_t R2;
    fp_zero(&R2.img);
    fp_copy(&r2, &R2.real);

    fp_zero(&one.real);
    fp_zero(&one.img);
    one.real[0] = 1;
    fp2_mul(&one, &R2, &one2);

    REDCL2(&one2, mod, &one_mont);
}

void REDC(fp_t* T, fp_t* mod, fp_t* res)
{
    fp_t tmr;
    f2p_t tmrmn;
    fp_t m;

    fp_mod_2k(T, &r, &tmr);             // tmr <- T mod r
    fp_mul(&tmr, &n_strich, &tmrmn);    // tmrmn <- tmr * n_strich
    f2p_mod_2k(&tmrmn, &r, &m);         // m <- tmrmn mod r



    f2p_t mN;
    f2p_t T2;
    f2p_t t;

    fp_to_f2p(T, &T2);

    fp_mul(&m, mod, &mN);              // mN <- mod * m
    f2p_add(&mN, &T2, &t);             // t <- T + mN



    

    f2p_t mod2;
    fp_to_f2p(mod, &mod2);
    int shift = get_shift(&r);
    f2p_rshift(&t, shift);            // t <- t / R (= (T + mN) / R)

    if( f2p_greater_equ_pos(&t, &mod2))
        f2p_to_fp(&t, res);
    else
    {
        fp_t t2;
        f2p_to_fp(&t, &t2);
        fp_sub(&t2, mod, res);
    }
        
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

void REDCL2(f2p2_t* T, fp_t* mod, fp2_t* res)
{
    fp2_t tmr;       // T mod R
    f2p2_t tmrmn;   // (T mod R) * n_strich
    fp2_t m;

    f2p_mod_2k(&T->real, &r, &tmr.real);              // T mod R 
    f2p_mod_2k(&T->img, &r, &tmr.img);                // T mod R  


    fp_mul(&tmr.real, &n_strich, &tmrmn.real);       // (T mod R) * n_strich
    fp_mul(&tmr.img, &n_strich, &tmrmn.img);       // (T mod R) * n_strich

    f2p_mod_2k(&tmrmn.real, &r, &m.real);              // m <- (T mod R) * n_strich mod R
    f2p_mod_2k(&tmrmn.img, &r, &m.img);              // m <- (T mod R) * n_strich mod R

    f2p2_t mn;       //m*n
    f2p2_t mnpt;     //m*n+T
    f2p2_t t;        //m*n+T/R

    f2p_t R;
    f2p_zero(&R);
    for(int i = 0; i < WORDS; i++)
    {
        R[i] = r[i];
    }


    fp_mul(&m.real, mod, &mn.real);         // mn <- m * N
    fp_mul(&m.img, mod, &mn.img);           // mn <- m * N


    f2p_add(&mn.real, &T->real, &mnpt.real);    // mnpt <- m * N + T
    f2p_add(&mn.img, &T->img, &mnpt.img);       // mnpt <- m * N + T

    f2p_div_2k(&mnpt.real, &R, &t.real);
    f2p_div_2k(&mnpt.img, &R, &t.img);


    fp2_t t_small;
    f2p_to_fp(&t.real, &t_small.real);
    f2p_to_fp(&t.img, &t_small.img);

    fp2_t mod2;

    fp_copy(mod, &mod2.real);
    fp_zero(&mod2.img);

    if(fp2_greater_equ_pos(&t_small, &mod2))
    {
        fp_sub(&t_small.real, mod, &res->real);
        fp_copy(&t_small.img, &res->img);
    }
    else
    {
        fp2_copy(&t_small, res);
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
