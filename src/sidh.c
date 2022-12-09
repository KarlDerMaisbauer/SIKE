

#include "sidh.h"
#include "montgomery_curve.h"
#include "montgomory_redc.h"
#include "point_arithmetic.h"
#include "isogeny_strat.h"

//public params  points need to be in montgomory form !!!!!!!!!!!!!!!

void isogen_2(public_params_t params, fp_t* secret_key, fp_t* mod, public_key_t* p2, isogeny_strat_t* strategy_e2)
{
    mont_curve_t A;
    mont_curve_t A_plus
    // 1. -----------------------------
    fp2_copy(&six_mont, &A.A);
    fp2_copy(&one_mont, &A.C);

    fp2_copy(&eight_mont, &A_plus.A);
    fp2_copy(&four_mont, &A_plus.C);
    // 1. -----------------------------

    // 2. -----------------------------
    proj_point_t X1;
    proj_point_t X2;
    proj_point_t X3;

    X1.X = params.X_P3;
    X1.Z = one_mont;

    X2.X = params.X_Q3;
    X2.Z = one_mont;

    X3.X = params.X_R3;
    X3.Z = one_mont;
    // 2. -----------------------------
    proj_point_t Xs;
    
    Ladder3p(&params.xP2, &params.xQ2, &params.xR2, secret_key, &Xs, &A, mod); // 3.

    proj_point_t opt_input[3];
    opt_input[0] = X1;
    opt_input[1] = X2;
    opt_input[2] = X3;
    mont_curve_t A_plus_new;

    // need variable for strategy
    e_2_iso(&A_24, &Xs, opt_input, 3, mod, strategy_e2, &A24_new);          // 4.
    fp2_copy(&A24_new, &A24);


    // public key 
    fp2_t X1_temp;
    fp2_t X2_temp;
    fp2_t X3_temp;
    fp2_div_mont(&opt_input[0].X , &opt_input[0].Z, mod, &X1_temp);
    fp2_div_mont(&opt_input[1].X , &opt_input[1].Z, mod, &X2_temp);
    fp2_div_mont(&opt_input[2].X , &opt_input[2].Z, mod, &X3_temp);
    f2p2_t X1_L_temp;
    f2p2_t X2_L_temp;
    f2p2_t X3_L_temp;
    fp2_to_f2p2(&x1_temp, &X1_L_temp);
    fp2_to_f2p2(&x2_temp, &X2_L_temp);
    fp2_to_f2p2(&x3_temp, &X3_L_temp);
    REDCL2(&X1_L_temp,  mod, &p2->X1);
    REDCL2(&X2_L_temp,  mod, &p2->X2);
    REDCL2(&X3_L_temp,  mod, &p2->X3);
}

void isogen_3(public_params_t params, fp_t* secret_key, fp_t* mod, public_key_t* p3, isogeny_strat_t* strategy_e3)
{
    mont_curve_t A;
    mont_curve_t A_plus

    // 1. -----------------------------
    fp2_copy(&six_mont, &A.A);
    fp2_copy(&one_mont, &A.C);

    fp2_copy(&eight_mont, &A_plus.A);
    fp2_copy(&four_mont, &A_plus.C);
    // 1. -----------------------------

    // 2. -----------------------------
    proj_point_t X1;
    proj_point_t X2;
    proj_point_t X3;

    X1.X = params.X_P2;
    X1.Z = one_mont;

    X2.X = params.X_Q2;
    X2.Z = one_mont;

    X3.X = params.X_R2;
    X3.Z = one_mont;
    // 2. -----------------------------
    proj_point_t Xs;
    
    Ladder3p(&params.xP3, &params.xQ3, &params.xR3, secret_key, &Xs, &A, mod);  // 3.

    proj_point_t opt_input[3];
    opt_input[0] = X1;
    opt_input[1] = X2;
    opt_input[2] = X3;
    mont_curve_t A_plus_new;

    // need variable for strategy
    e_3_iso(&A_24, &Xs, opt_input, 3, mod, strategy_e3, &A24_new);
    fp2_copy(&A24_new, &A24);


    // public key 
    fp2_t X1_temp;
    fp2_t X2_temp;
    fp2_t X3_temp;
    fp2_div_mont(&opt_input[0].X , &opt_input[0].Z, mod, &X1_temp);
    fp2_div_mont(&opt_input[1].X , &opt_input[1].Z, mod, &X2_temp);
    fp2_div_mont(&opt_input[2].X , &opt_input[2].Z, mod, &X3_temp);
    f2p2_t X1_L_temp;
    f2p2_t X2_L_temp;
    f2p2_t X3_L_temp;
    fp2_to_f2p2(&x1_temp, &X1_L_temp);
    fp2_to_f2p2(&x2_temp, &X2_L_temp);
    fp2_to_f2p2(&x3_temp, &X3_L_temp);
    REDCL2(&X1_L_temp,  mod, &p2->X1);
    REDCL2(&X2_L_temp,  mod, &p2->X2);
    REDCL2(&X3_L_temp,  mod, &p2->X3);
}










void isoex_2(fp_t* secret_key, public_key_t* public_key, fp2_t* j_inv, fp_t* mod, isogeny_strat_t* strategy_e2)
{
    mont_curve_t A;
    get_A(&public_key->X1, &public_key->X2, &public_key->X3, mod, &A.A);
    fp2_copy(&one_mont, &A.C);

    proj_point_t Xs;
    Ladder3p(&public_key->X1, &public_key->X2, &public_key->X3, secret_key, &Xs, &A, mod);

    mont_curve_t A24;
    fp2_add(&A.A, &two_mont, &A24.A);
    A24.C = four_mont;

    mont_curve_t A24_New;
    e_2_iso(&A24, &Xs, 0, 0, mod, strategy_e2, &A24_new);

    fp2_t A24_4;
    fp2_t C24_2;
    curve_copy(&A24_new, &A24);


    fp2_mul_mont(&A24.A, &four_mont, mod, &A24_4);
    fp2_mul_mont(&A24.C, &two_mont,  mod, &C24_2);

    fp2_sub(&A24_4, C24_2, &A.A);
    fp2_copy(&A.C, &A24.C);

    fp2_t j_inv_mont;
    f2p2_t j_inv_mont_L;

    jinvariant(&A.A, &A.C, mod, &j_inv_mont);
    fp2_to_f2p2(&j_inv_mont, &j_inv_mont_L);

    REDCL2(&j_inv_mont_L, mod, j_inv);
}

void isoex_3(fp_t* secret_key, public_key_t* public_key, fp2_t* j_inv, fp_t* mod, isogeny_strat_t* strategy_e3)
{
    mont_curve_t A;
    get_A(&public_key->X1, &public_key->X2, &public_key->X3, mod, &A.A);
    fp2_copy(&one_mont, &A.C);

    proj_point_t Xs;
    Ladder3p(&public_key->X1, &public_key->X2, &public_key->X3, secret_key, &Xs, &A, mod);



    mont_curve_t A24;
    fp2_add(&A.A, &two_mont, &A24.A);
    fp2_sub(&A.A, &two_mont, &A24.C);

    mont_curve_t A24_New;
    e_3_iso(&A24, &Xs, 0, 0, mod, strategy_e3, &A24_new);
    curve_copy(&A24_new, &A24);

    fp2_t A24_p;        // A24^- + A24^+


    fp2_add(&A24.C, &A24.A, &A24_p);
    fp2_sub(&A24.A, &A24.C, &A.C);        // A24^+ - A24^-

    fp2_mul_mont(&A24_p, two_mont, mod, &A.A)                // 2 * (A24^- + A24^+)


    fp2_t j_inv_mont;
    f2p2_t j_inv_mont_L;

    jinvariant(&A.A, &A.C, mod, &j_inv_mont);
    fp2_to_f2p2(&j_inv_mont, &j_inv_mont_L);

    REDCL2(&j_inv_mont_L, mod, j_inv);
}