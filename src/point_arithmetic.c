


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

    fp2_add(&Z2p_temp, &t1minus, &Z2p_temp_2);          // 9.
    fp2_mul_mont(&Z2p_temp_2, &t1minus, mod, &P2->Z);   // 10.
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
    fp2_t t0;
    fp2_t t1;
    fp2_add(&P->X, &P->Z, &t0);                    // 1. t0 <- Xp + Zp
    fp2_sub(&P->X, &P->Z, &t1);                    // 2. t1 <- Xp - Zp

    fp2_t t0sq;
    fp2_t t2;
    fp2_mul_mont(&t0, &t0, mod, &t0sq);             // 3.  t0sq <-t0^2   (t0sq = X2p)
    fp2_sub(&Q->X, &Q->Z, &t2);                     // 4.  t2 <- Xq - Zq

    fp2_t X_QpP;
    fp2_t t02;
    fp2_add(&Q->X, &Q->Z, &X_QpP);                  // 5.  X_P+Q <- Xq + Zq (X_QpP = X_P+Q)
    fp2_mul_mont(&t0, &t2, mod, &t02);              // 6.  t0 <- t0 * t2 (t02 = t0)

    fp2_t Z2p;
    fp2_t t12;
    fp2_mul_mont(&t2, &t2, mod, &Z2p);              // 7.  Z_[2]p <- t2^2 (Z2p = Z_[2]p)
    fp2_mul_mont(&t1, &X_QpP, mod, &t12);           // 8.  t1 <- t1 * X_R+Q (t12 = t1)


    fp2_sub(&t0sq, &Z2p, &t2);                      // 9.  t2 <- X_[2]p - Z_[2]p
    fp2_mul_mont(&t0sq, &Z2p, mod, &P2->X);         // 10. X_[2]p <- X_[2]p * Z_[2]p ret val for X_[2]p set

    fp2_t Z_PpQ;
    fp2_mul_mont(&A24->A, &t2, mod, &X_QpP);        // 11. X_P+Q <- a^+_24 * t2 (X_QpP = X_P+Q)
    fp2_sub(&t02, &t12, &Z_PpQ);                    // 12. Z_P+Q <- t0 - t1 (Z_PpQ = Z_P+Q)

    fp2_t Z2p_2;
    fp2_add(&X_QpP, &Z2p, &Z2p_2);                  // 13. Z_[2]p <- X_P+Q + Z_[2]p (Z2p_2 = Z_[2]p)
    fp2_add(&t02, &t12, &X_QpP);                    // 14. X_P+Q <- t0 + t1 (X_QpP = X_P+Q )

    fp2_t ZPpQ_sq;
    fp2_mul_mont(&Z2p_2, &t2, mod, &P2->Z);         // 15. Z_[2]p <- Z[2]p * t2  ret val for Z_[2]p set
    fp2_mul_mont(&Z_PpQ, &Z_PpQ, mod, &ZPpQ_sq);    // 16. Z_P+Q <- Z_P+Q^2 (ZPpQ_sq = Z_P+Q)

    fp2_t XPpQ_sq;
    fp2_mul_mont(&X_QpP, &X_QpP, mod, &XPpQ_sq);    // 17. X_P+Q <- X_P+Q^2 (XPpQ_sq = X_P+Q)
    fp2_mul_mont(&PmQ->X, &ZPpQ_sq, mod, &PpQ->Z);  // 18. Z_P+Q <- X_Q-P * Z_P+Q ret val for PpQ set

    fp2_mul_mont(&PmQ->Z, &XPpQ_sq, mod, &PpQ->X);  // 19. X_P+Q <- Z_Q-P * X_P+Q ret val for PpQ set
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

    fp2_copy(&t1, &c1);                             // c1 = t1 for 8.
    fp2_mul_mont(&t4, &t4, mod, &t1);               // 7.  t1 <- t4^2
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

void Ladder3p(fp2_t* P, fp2_t* Q, fp2_t* QmP, fp_t* m, fp2_t* PpmQ, mont_curve_t* A, fp_t* mod)
{
   int len = fp_get_len(m);

    proj_point_t P0;
    proj_point_t P1;
    proj_point_t P2;

    mont_curve_t a24;
    fp2_copy(&one_mont, &a24.C);

    fp2_t a24_val;
    fp2_zero(&a24_val);
    fp2_add(&A->A, &four_mont, &a24_val);
    fp2_copy(&a24_val, &a24.A);


    P0.X = *P;
    P0.Z = one_mont;

    P1.X = *Q;
    P1.Z = one_mont;

    P2.X = *QmP;
    P2.Z = one_mont;

        

    for(int i =  0; i < len; i++)
    {
        uint64_t condition = 1 & ((*m)[(i/64)]) << (i % 64);
        LADDER_inner_part(&P0, &P1, &P2, &a24, condition, mod);
    }
    fp2_copy(&P1.X, PpmQ);

}