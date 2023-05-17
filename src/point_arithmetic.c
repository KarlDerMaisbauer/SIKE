


#include"point_arithmetic.h"
#include"point_arithmetic_helper.h"
#include"fp2.h"
#include"fp_helper.h"
#include"montgomory_redc.h"



void xDBL(proj_point_t* P, mont_curve_t* A24, fp_t* mod, proj_point_t* P2)
{

    fp2_t t0temp;
    fp2_t t1temp;
    fp2_t t0;
    fp2_t t1;


    fp2_sub(&P->X, &P->Z, &t0temp);                     // 1.
    fp2_add(&P->X, &P->Z, &t1temp);                     // 2.

    fp2_mul_mont(&t0temp, &t0temp, mod, &t0);           // 3.
    fp2_mul_mont(&t1temp, &t1temp, mod, &t1);           // 4.

    fp2_t Z2p_temp;

    fp2_mul_mont(&A24->C, &t0, mod, &Z2p_temp);         // 5.
    fp2_mul_mont(&Z2p_temp, &t1, mod, &P2->X);          // 6. X_[2]p set

    fp2_t t1minus;

    fp2_sub(&t1, &t0, &t1minus);                        // 7.
    fp2_mul_mont(&A24->A, &t1minus, mod, &t0);          // 8.

    fp2_t Z2p_temp_2;

    fp2_add(&Z2p_temp, &t0, &Z2p_temp_2);          // 9.
    fp2_mul_mont(&Z2p_temp_2, &t1minus, mod, &P2->Z);   // 10.
}

void xDBL_no_redc(proj_point_t* P, mont_curve_t* A24, fp_t* mod, proj_point_t* P2)
{

    fp2_t t0temp;
    fp2_t t1temp;
    fp2_t t0;
    fp2_t t1;
    f2p2_t t0_l;
    f2p2_t t1_l;


    fp2_sub(&P->X, &P->Z, &t0temp);                     // 1.
    fp2_add(&P->X, &P->Z, &t1temp);                     // 2.

    fp2_mul(&t0temp, &t0temp, &t0_l);           // 3.
    fp2_mul(&t1temp, &t1temp, &t1_l);           // 4.
    f2p_mod(&t0_l.real, mod, &t0.real);
    f2p_mod(&t0_l.img, mod, &t0.img);
    f2p_mod(&t1_l.real, mod, &t1.real);
    f2p_mod(&t1_l.img, mod, &t1.img);

    fp2_t Z2p_temp;
    f2p2_t Z2p_temp_l;

    fp2_mul(&A24->C, &t0, &Z2p_temp_l);         // 5.
    f2p_mod(&Z2p_temp_l.real, mod, &Z2p_temp.real);
    f2p_mod(&Z2p_temp_l.img, mod, &Z2p_temp.img);
    fp2_mul(&Z2p_temp, &t1, &Z2p_temp_l);
    //f2p2_mod(&Z2p_temp_l, mod, &P2->X);          // 6. X_[2]p set
    f2p_mod(&Z2p_temp_l.real, mod, &P2->X.real);
    f2p_mod(&Z2p_temp_l.img, mod, &P2->X.img);

    fp2_t t1minus;

    fp2_sub(&t1, &t0, &t1minus);                        // 7.
    fp2_mul(&A24->A, &t1minus, &t0_l);          // 8.
    f2p_mod(&t0_l.real, mod, &t0.real);
    f2p_mod(&t0_l.img, mod, &t0.img);

    fp2_t Z2p_temp_2;

    fp2_add(&Z2p_temp, &t0, &Z2p_temp_2);          // 9.
    fp2_mul(&Z2p_temp_2, &t1minus, &Z2p_temp_l); //&P2->Z);   // 10.
    f2p_mod(&Z2p_temp_l.real, mod, &P2->Z.real);
    f2p_mod(&Z2p_temp_l.img, mod, &P2->Z.img);
}

void xDBLe(proj_point_t* P, mont_curve_t* A24, int64_t e, fp_t* mod, proj_point_t* P2e)
{
    proj_point_t x_inner;
    proj_pt_copy(P, &x_inner);

    for(int64_t i = 1; i <= e; i++)
    {
        proj_point_t x_part;
        proj_pt_zero(&x_part);
        xDBL(&x_inner, A24, mod, &x_part);
        proj_pt_copy(&x_part, &x_inner);
        
    }
    proj_pt_copy(&x_inner, P2e);
}

void xDBLADD(proj_point_t* P, proj_point_t* Q, proj_point_t* PmQ, fp_t* mod, mont_curve_t* A24, proj_point_t* P2, proj_point_t* PpQ)
{
    // copy variables
    fp2_t c1;



    fp2_t t0;
    fp2_t t1;
    fp2_add(&P->X, &P->Z, &t0);                     // 1. t0 <- Xp + Zp
    fp2_sub(&P->X, &P->Z, &t1);                     // 2. t1 <- Xp - Zp

    fp2_t t2;
    fp2_mul_mont(&t0, &t0, mod, &P2->X);            // 3.  X_[2]P <-t0^2   
    fp2_sub(&Q->X, &Q->Z, &t2);                     // 4.  t2 <- Xq - Zq

    fp2_add(&Q->X, &Q->Z, &PpQ->X);                 // 5.  X_P+Q <- Xq + Zq
    fp2_copy(&t0, &c1);                             // c1 = t0 for 6.
    fp2_mul_mont(&c1, &t2, mod, &t0);               // 6.  t0 <- t0 * t2


    fp2_mul_mont(&t1, &t1, mod, &P2->Z);            // 7.  Z_[2]p <- t2^2
    fp2_copy(&t1, &c1);                             // c1 = t1 for 8
    fp2_mul_mont(&c1, &PpQ->X, mod, &t1);           // 8.  t1 <- t1 * X_R+Q


    fp2_sub(&P2->X, &P2->Z, &t2);                   // 9.  t2 <- X_[2]p - Z_[2]p
    fp2_copy(&P2->X, &c1);                          // c1 = X_[2]p for 10.
    fp2_mul_mont(&c1, &P2->Z, mod, &P2->X);         // 10. X_[2]p <- X_[2]p * Z_[2]p ret val for X_[2]p set

    fp2_mul_mont(&A24->A, &t2, mod, &PpQ->X);       // 11. X_P+Q <- a^+_24 * t2 
    fp2_sub(&t0, &t1, &PpQ->Z);                     // 12. Z_P+Q <- t0 - t1 

    fp2_copy(&P2->Z, &c1);                           // c1 = Z_[2]p for 13.
    fp2_add(&PpQ->X, &c1, &P2->Z);                  // 13. Z_[2]p <- X_P+Q + Z_[2]p 
    fp2_add(&t0, &t1, &PpQ->X);                     // 14. X_P+Q <- t0 + t1

    fp2_copy(&P2->Z, &c1);                          // c1 = Z_[2]P for 15.
    fp2_mul_mont(&c1, &t2, mod, &P2->Z);            // 15. Z_[2]p <- Z[2]p * t2  ret val for Z_[2]p set
    fp2_copy(&PpQ->Z, &c1);                          // c1 = Z_P+Q for 16.
    fp2_mul_mont(&c1, &c1, mod, &PpQ->Z);           // 16. Z_P+Q <- Z_P+Q^2 

    fp2_copy(&PpQ->X, &c1);                         // c1 = X_P+Q for 17
    fp2_mul_mont(&c1, &c1, mod, &PpQ->X);           // 17. X_P+Q <- X_P+Q^2 (XPpQ_sq = X_P+Q)
    fp2_copy(&PpQ->Z, &c1);                         // c1 = Z_P+Q for 18
    fp2_mul_mont(&PmQ->X, &c1, mod, &PpQ->Z);       // 18. Z_P+Q <- X_Q-P * Z_P+Q ret val for PpQ set

    fp2_copy(&PpQ->X, &c1);                         // c1 = X_P+Q for 19
    fp2_mul_mont(&PmQ->Z, &c1, mod, &PpQ->X);  // 19. X_P+Q <- Z_Q-P * X_P+Q ret val for PpQ set
}

void xTPL(proj_point_t* P, mont_curve_t* A24, fp_t* mod, proj_point_t* P3)
{
    // copy variables
    fp2_t c1;



    fp2_t t0;
    fp2_t t2;
    fp2_sub(&P->X, &P->Z, &t0);                     // 1.  t0 <- X_P - Z_P
    fp2_mul_mont(&t0, &t0, mod, &t2);               // 2.  t2 <- t_0^2

    fp2_t t1;
    fp2_t t3;
    fp2_add(&P->X, &P->Z, &t1);                     // 3.  t1 <- X_P + Z_P
    fp2_mul_mont(&t1, &t1, mod, &t3);               // 4.  t3 <- t_1^2

    fp2_t t4;
    fp2_copy(&t0, &c1);                             // c1 = t0 for 6.
    fp2_add(&t0, &t1, &t4);                         // 5.  t4 <- t1 + t0
    fp2_sub(&t1, &c1, &t0);                         // 6.  t0 <- t1 - t0

    
    fp2_mul_mont(&t4, &t4, mod, &t1);               // 7.  t1 <- t4^2
    fp2_copy(&t1, &c1);                             // c1 = t1 for 8.
    fp2_sub(&c1, &t3, &t1);                         // 8.  t1 <- t1 - t3

    fp2_t t5;
    fp2_copy(&t1, &c1);                             // c1 = t1 for 9.
    fp2_sub(&c1, &t2, &t1);                         // 9.  t1 <- t1 -t1
    fp2_mul_mont(&t3, &A24->A, mod, &t5);           // 10. t5 <- t3 * A24^+

    fp2_t t6;
    fp2_copy(&t3, &c1);                             // c1 = t3 for 11.
    fp2_mul_mont(&c1, &t5, mod, &t3);               // 11. t3 <- t5 * t3
    fp2_mul_mont(&t2, &A24->C, mod, &t6);           // 12. t6 <- t2 * A24^-

    fp2_copy(&t2, &c1);                             // c1 = t2 for 13.
    fp2_mul_mont(&c1, &t6, mod, &t2);               // 13. t2 <- t2 * t6
    fp2_copy(&t3, &c1);                             // c1 = t3 for 14.
    fp2_sub(&t2, &c1, &t3);                         // 14. t3 <- t2 - t3

    fp2_copy(&t1, &c1);                             // c1 = t1 for 16
    fp2_sub(&t5, &t6, &t2);                         // 15. t2 <- t5 - t6
    fp2_mul_mont(&t2, &c1, mod, &t1);               // 16. t1 <- t2 * t1

    fp2_add(&t3, &t1, &t2);                         // 17. t2 <- t3 + t1
    fp2_copy(&t2, &c1);                             // c1 = t2 for 18.
    fp2_mul_mont(&c1, &c1, mod, &t2);               // 18. t2 <- t2^2

    fp2_copy(&t1, &c1);                             // c1 = t1 for 20.
    fp2_mul_mont(&t2, &t4, mod, &P3->X);            // 19. X_[3]P <- t2 * t4
                                                    // retval for X_[3]P set
    fp2_sub(&t3, &c1, &t1);                         // 20. t1 <- t3 - t1

    fp2_copy(&t1, &c1);                             // c1 = t1 for 21.
    fp2_mul_mont(&c1, &c1, mod, & t1);              // 21. t1 <- t1^2
    fp2_mul_mont(&t1, &t0, mod, &P3->Z);            // Z_[3]P <- t1 * t0
}   

void xTPL_no_redc(proj_point_t* P, mont_curve_t* A24, fp_t* mod, proj_point_t* P3)
{
    // copy variables
    fp2_t c1;
    f2p2_t mul_tmp;
    f2p2_zero(&mul_tmp);



    fp2_t t0;
    fp2_t t2;
    fp2_sub(&P->X, &P->Z, &t0);                     // 1.  t0 <- X_P - Z_P
    fp2_mul(&t0, &t0, &mul_tmp);               // 2.  t2 <- t_0^2
    f2p2_mod(&mul_tmp, mod, &t2);

    fp2_t t1;
    fp2_t t3;
    fp2_add(&P->X, &P->Z, &t1);                     // 3.  t1 <- X_P + Z_P
    fp2_mul(&t1, &t1, &mul_tmp);               // 4.  t3 <- t_1^2
    f2p2_mod(&mul_tmp, mod, &t3);

    fp2_t t4;
    fp2_copy(&t0, &c1);                             // c1 = t0 for 6.
    fp2_add(&t0, &t1, &t4);                         // 5.  t4 <- t1 + t0
    fp2_sub(&t1, &c1, &t0);                         // 6.  t0 <- t1 - t0

                                // c1 = t1 for 8.
    fp2_mul(&t4, &t4, &mul_tmp);               // 7.  t1 <- t4^2
    f2p2_mod(&mul_tmp, mod, &t1);
    fp2_copy(&t1, &c1); 
    fp2_sub(&c1, &t3, &t1);                         // 8.  t1 <- t1 - t3

    fp2_t t5;
    fp2_copy(&t1, &c1);                             // c1 = t1 for 9.
    fp2_sub(&c1, &t2, &t1);                         // 9.  t1 <- t1 -t1
    fp2_mul(&t3, &A24->A, &mul_tmp);           // 10. t5 <- t3 * A24^+
    f2p2_mod(&mul_tmp, mod, &t5);

    fp2_t t6;
    fp2_copy(&t3, &c1);                             // c1 = t3 for 11.
    fp2_mul(&c1, &t5, &mul_tmp);               // 11. t3 <- t5 * t3
    f2p2_mod(&mul_tmp, mod, &t3);
    fp2_mul(&t2, &A24->C, &mul_tmp);           // 12. t6 <- t2 * A24^-
    f2p2_mod(&mul_tmp, mod, &t6);

    fp2_copy(&t2, &c1);                             // c1 = t2 for 13.
    fp2_mul(&c1, &t6, &mul_tmp);               // 13. t2 <- t2 * t6
    f2p2_mod(&mul_tmp, mod, &t2);
    fp2_copy(&t3, &c1);                             // c1 = t3 for 14.
    fp2_sub(&t2, &c1, &t3);                         // 14. t3 <- t2 - t3

    fp2_copy(&t1, &c1);                             // c1 = t1 for 16
    fp2_sub(&t5, &t6, &t2);                         // 15. t2 <- t5 - t6
    fp2_mul(&t2, &c1, &mul_tmp);               // 16. t1 <- t2 * t1
    f2p2_mod(&mul_tmp, mod, &t1);

    fp2_add(&t3, &t1, &t2);                         // 17. t2 <- t3 + t1
    fp2_copy(&t2, &c1);                             // c1 = t2 for 18.
    fp2_mul(&c1, &c1, &mul_tmp);               // 18. t2 <- t2^2
    f2p2_mod(&mul_tmp, mod, &t2);

    fp2_copy(&t1, &c1);                             // c1 = t1 for 20.
    fp2_mul(&t2, &t4, &mul_tmp);            // 19. X_[3]P <- t2 * t4
                                                    // retval for X_[3]P set
    f2p2_mod(&mul_tmp, mod, &P3->X)   ;                                             
    fp2_sub(&t3, &c1, &t1);                         // 20. t1 <- t3 - t1

    fp2_copy(&t1, &c1);                             // c1 = t1 for 21.
    fp2_mul(&c1, &c1, &mul_tmp);              // 21. t1 <- t1^2
    f2p2_mod(&mul_tmp, mod, &t1);
    fp2_mul(&t1, &t0, &mul_tmp);            // Z_[3]P <- t1 * t0
    f2p2_mod(&mul_tmp, mod, &P3->Z);
}   

void xTPLe(proj_point_t* P, mont_curve_t* A24, int64_t e, fp_t* mod, proj_point_t* P3e)
{
    proj_point_t x_inner;
    proj_pt_copy(P, &x_inner);
    for(int64_t i = 1; i <= e; i++)
    {
        proj_point_t x_part;
        proj_pt_zero(&x_part);
        xTPL(&x_inner, A24, mod, &x_part);
        proj_pt_copy(&x_part, &x_inner);
        
    }
    proj_pt_copy(&x_inner, P3e);
}

void Ladder3p(fp2_t* P, fp2_t* Q, fp2_t* QmP, fp_t* m, proj_point_t* PpmQ, mont_curve_t* A, fp_t* mod)
{
    int len = fp_get_len(m);

    proj_point_t P0;
    proj_point_t P1;
    proj_point_t P2;

    mont_curve_t a24;
    fp2_copy(&one_mont, &a24.C);
    fp2_t temp;
    fp2_zero(&temp);
    fp2_add(&A->A, &two_mont, &temp);
    fp2_div_mont(&temp, &four_mont, mod, &a24.A);
    fp_t m_temp;
    fp_copy(m, &m_temp);

    //P0.X = *P;
    fp2_copy(Q, &P0.X);
    //P0.Z = one_mont;
    fp2_copy(&one_mont, &P0.Z);

    //P1.X = *Q;
    //P1.Z = one_mont;
    fp2_copy(P, &P1.X);
    fp2_copy(&one_mont, &P1.Z);

    //P2.X = *QmP;
    fp2_copy(QmP, &P2.X);
    //P2.Z = one_mont;
    fp2_copy(&one_mont, &P2.Z);

    for(int i =  0; i < len; i++)
    {
        //uint64_t condition = 1 & (((*m)[(i/64)]) >> (i % 64));
        uint64_t condition = 1 & m_temp[0];
        fp_rshift(&m_temp, 1); 

        //LADDER_inner_part(&P0, &P1, &P2, &a24, condition, mod);
        proj_point_t P0_temp;
        proj_point_t b_out;
        if(condition)
        {
            xDBLADD(&P0, &P1, &P2, mod, &a24, &P0_temp, &b_out);
            proj_pt_copy(&P0_temp, &P0);
            proj_pt_copy(&b_out, &P1);
        }
        else
        {
            xDBLADD(&P0, &P2, &P1, mod, &a24, &P0_temp, &b_out);
            proj_pt_copy(&P0_temp, &P0);
            proj_pt_copy(&b_out, &P2);
        }

    }
    //fp2_copy(&P1.X, &PpmQ->X);
    //fp2_copy(&P1.Z, &PpmQ->Z);
    proj_pt_copy(&P1, PpmQ);
}