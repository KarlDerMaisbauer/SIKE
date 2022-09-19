

#include"montgomery_curve.h"

#include"fp2.h"
#include"fp.h"
#include"fp_helper.h"

void jinvariant(fp2_t* A, fp2_t* C, fp_t* mod, fp2_t* j)
{
    // copy vars
    fp2_t c1;

    fp2_t t1;
    fp2_mulm(A, A, mod, j);                     // 1.  j <- A^2
    fp2_mulm(C, C, mod, &t1);                   // 2.  t1 <- C^2

    fp2_t t0;
    fp2_addm(&t1, &t1, mod, &t0);               // 3.  t0 <- t1 + t1
    fp2_copy(&t0, &c1);                         // c1 = t0 for 4.
    fp2_subm(j, &c1, mod, &t0);                 // 4.  t0 <- j - t0

    fp2_copy(&t0, &c1);                         // c1 = t0 for 5.
    fp2_subm(&c1, &t1, mod, &t0);               // 5.  t0 <- t0 - t1
    fp2_subm(&t0, &t1, mod, j);                 // 6.  j <- t0 - t1

    fp2_copy(&t1, &c1);                         // c1 = t1 for 7.
    fp2_mulm(&c1, &c1, mod, &t1);               // 7.  t1 <- t1^2
    fp2_copy(j, &c1);                           // c1 = j for 8.
    fp2_mulm(&c1, &t1, mod, j);                 // 8.  j <- j * t1

    fp2_copy(&t0, &c1);                         // c1 = t0 for 9.
    fp2_addm(&c1, &c1, mod, &t0);               // 9.  t0 <- t0 + t0
    fp2_copy(&t0, &c1);                         //c1 = t0 for 10.
    fp2_addm(&c1, &c1, mod, &t0);               // 10. t0 <- t0 + t0

    fp2_mulm(&t0, &t0, mod, &t1);               // 11. t1 <- t0^2
    fp2_copy(&t0, &c1);                         // c1 = t0 for 12.
    fp2_mulm(&c1, &t1, mod, &t0);               // 12. t0 <- t0 * t1

    fp2_copy(&t0, &c1);                         // c1 = t0 for 13.
    fp2_addm(&c1, &c1, mod, &t0);               // 13. t0 <- t0 + t0
    fp2_copy(&t0, &c1);                         // c1 = t0 for 14.
    fp2_addm(&c1, &c1, mod, &t0);               // 14. t0 <- t0 + t0

    fp2_copy(j, &c1);                           // c1 = j for 15.
    fp2_mult_inv(&c1, mod, j);                  // 15. j <- 1/j
    fp2_copy(j, &c1);                           // c1 = j for 16.
    fp2_mulm(&c1, &t0, mod, j);                 // 16. j <- j * t0
}

void get_A(fp2_t* xp, fp2_t* xq, fp2_t* xqp, fp_t* mod, fp2_t* A)
{
    // copy vars
    fp2_t c1;

    fp2_t one;
    fp_zero(&one.real);
    fp_zero(&one.img);
    one.real[0] = 1;

    fp2_t t0;
    fp2_t t1;
    fp2_addm(xp, xq, mod, &t1);                 // 1.  t1 <- xp + xq
    fp2_mulm(xp, xq, mod, &t0);                 // 2.  t0 <- xp * xq

    fp2_mulm(xqp, &t1, mod, A);                 // 3.  A <- xqp * t1
    fp2_copy(A, &c1);                           // c1 = A for 4.
    fp2_addm(&c1, &t0, mod, A);                 // 4.  A <- A + t0

    fp2_copy(&t0, &c1);                         // c1 = t0 for 5.
    fp2_mulm(&c1, xqp, mod, &t0);               // 5.  t0 <- t0 * xqp
    fp2_copy(A, &c1);                           // c1 = A for 6.
    fp2_subm(&c1, &one, mod, A);                // 6. A <- A - 1

    fp2_copy(&t0, &c1);                         // c1 = t0 for 7.
    fp2_addm(&c1, &c1, mod, &t0);               // 7.  t0 <- t0 + t0
    fp2_copy(&t1, &c1);                         // c1 = t1 for 8.
    fp2_addm(&c1, xqp, mod, &t1);               // 8.  t1 <- t1 + xqp

    fp2_copy(&t0 , &c1);                        // c1 = t0 for 9.
    fp2_addm(&c1, &c1, mod, &t0);               // 9.  t0 <- t0 + t0
    fp2_copy(A, &c1);                           // c1 = A for 10.
    fp2_mulm(&c1, &c1, mod, A);                 // 10. A <- A^2

    fp2_copy(&t0, &c1);                         // c1 = t0 for 11.
    fp2_mult_inv(&c1, mod, &t0);                // 11. t0 <- 1/t0
    fp2_copy(A, &c1);                           // c1 = A for 12.
    fp2_mulm(&c1, &t0, mod, A);                 // A <- A * t0

    fp2_copy(A, &c1);                           // c1 = A for 13.
    fp2_subm(&c1, &t1, mod, A);                 // A <- A - t1
}

void iso_2_curve(proj_point_t* P2, fp_t* mod, mont_curve_t* A24)
{
    // copy vars
    fp2_t c1;

    fp2_mulm(&P2->X, &P2->X, mod, &A24->A);     // 1.  A_24^+ <- X_P2^2
    fp2_mulm(&P2->Z, &P2->Z, mod, &A24->C);     // 2.  C_24^+ <- Z_P2^2

    fp2_copy(&A24->A, &c1);                     // c1 = A_24^+ for 3.
    fp2_subm(&A24->C, &c1, mod, &A24->A);       // 3. A_24^+ <- C_24 - A_24^+
}

void iso_2_eval(proj_point_t* P2,proj_point_t* Q, fp_t* mod, proj_point_t* Q_strich)
{
    // copy vars
    fp2_t c1;

    fp2_t t0;
    fp2_t t1;
    fp2_addm(&P2->X, &P2->Z, mod, &t0);         // 1.  t0 <- X_p2 + Z_P2
    fp2_subm(&P2->X, &P2->Z, mod, &t1);         // 2.  t1 <- X_p2 - Z_P2


    fp2_t t2;
    fp2_t t3;
    fp2_addm(&Q->X, &Q->Z, mod, &t2);           // 3.  t2 <- X_Q + Z_Q
    fp2_subm(&Q->X, &Q->Z, mod, &t3);           // 4.  t3 <- X_Q - Z_Q

    fp2_copy(&t0, &c1);                         // c1 = t0 for 5.
    fp2_mulm(&c1, &t3, mod, &t0);               // 5.  t0 <- t0 * t3
    fp2_copy(&t1, &c1);                         // c1 = t1 for 6.
    fp2_mulm(&c1, &t2, mod, &t1);               // 6.  t1  <- t1 * t2

    fp2_addm(&t0, &t1, mod, &t2);               // 7.  t2 <- t0 + t1   
    fp2_subm(&t0, &t1, mod, &t3);               // 8.  t3 <- t0 - t1   

    fp2_mulm(&Q->X, &t2, mod, &Q_strich->X);    // 9.  X_Q' <- X_Q * t2
    fp2_mulm(&Q->Z, &t3, mod, &Q_strich->Z);    // 10. Z_Q' <- Z_Q * t3
}

void iso_4_curve(proj_point_t* P4, fp_t* mod, mont_curve_t* A24, fp2_t* K1, fp2_t* K2, fp2_t* K3)
{
    // copy vars
    fp2_t c1;

    fp2_subm(&P4->X, &P4->Z, mod, K2);          // 1.  K_2 <- X_P4 - Z_P4   K_2 set
    fp2_addm(&P4->X, &P4->Z, mod, K3);          // 2.  K_3 <- X_P4 + Z_P4   K_3 set

    fp2_mulm(&P4->Z, &P4->Z, mod, K1);          // 3.  K_1 <- Z_p4^2
    fp2_copy(K1, &c1);                          // c1 = K1 for 4.
    fp2_addm(&c1, &c1, mod, K1);                // 4.  K1  <- K1 + K1

    fp2_mulm(K1, K1, mod, &A24->C);             // 5.  C24 <- K1^2          C24 set
    fp2_copy(K1, &c1);                          // c1 = K1 for 6.
    fp2_addm(&c1, &c1, mod, K1);                // 6.  K1 <- K1 + K1        K1 set

    fp2_mulm(&P4->X, &P4->X, mod, &A24->A);     // 7.  A_24 <- X_P4^2
    fp2_copy(&A24->A, &c1);                     // c1 = A24 for 8.
    fp2_addm(&c1, &c1, mod, &A24->A);           // 8.  A_24 <- A_24 + A_24

    fp2_copy(&A24->A, &c1);                     // c1 = A24 for 9.
    fp2_mulm(&c1, &c1, mod, &A24->A);           // 9.  A24 <- A24^2
}

void iso_4_eval(fp2_t* K1, fp2_t* K2, fp2_t* K3, proj_point_t* Q, fp_t* mod, proj_point_t* Q_strich)
{
    // copy vars
    fp2_t c1;

    // additional var
    fp2_t X_Q;                                  // for X_Q
    fp2_t Z_Q;                                  // for Z_Q


    fp2_t t0;
    fp2_t t1;
    fp2_addm(&Q->X, &Q->Z, mod, &t0);           // 1.  t0 <- X_Q + Z_Q
    fp2_subm(&Q->X, &Q->Z, mod, &t1);           // 2.  t0 <- X_Q - Z_Q

    fp2_mulm(&t0, K2, mod, &X_Q);               // 3.  X_Q <- t0 * K_2
    fp2_mulm(&t1, K3, mod, &Z_Q);               // 4.  Z_Q <- t1 * K_3

    fp2_copy(&t0, &c1);                         // c1 = t0 for 5.
    fp2_mulm(&c1, &t1, mod, &t0);               // 5.  t0 <- t0 * t1
    fp2_copy(&t0, &c1);                         // c1 = t0 for 6.
    fp2_mulm(&c1, K1, mod, &t0);                // 6.  t0 <- t0 * K1

    fp2_addm(&X_Q, &Z_Q, mod, &t1);             // 7.  t1 <- X_Q + Z_Q
    fp2_copy(&Z_Q, &c1);                        // c1 = Z_Q for 8.
    fp2_subm(&X_Q, &c1, mod, &Z_Q);             // 8.  Z_Q <- X_Q - Z_Q

    fp2_copy(&t1, &c1);                         // c1 = t1 for 9.
    fp2_mulm(&c1, &c1, mod, &t1);               // 9.  t1 <- t1^2
    fp2_copy(&Z_Q, &c1);                        // c1 = Z_Q for 10.
    fp2_mulm(&c1, &c1, mod, &Z_Q);              // 10. Z_Q <- Z_Q^2

    fp2_addm(&t0, &t1, mod, &X_Q);              // 11. X_Q <- t0 + t1
    fp2_copy(&t0, &c1);                         // c1 = t0 for 12.
    fp2_subm(&Z_Q, &c1, mod, &t0);              // 12. t0 <- Z_Q - t0

    fp2_mulm(&X_Q, &t1, mod, &Q_strich->X);     // 13. X_Q' <- X_Q * t1
    fp2_mulm(&Z_Q, &t0, mod, &Q_strich->Z);     // 14. Z_Q' <- Z_Q * t0
}

void iso_3_curve(proj_point_t* P3, fp_t* mod, mont_curve_t* A24, fp2_t* K1, fp2_t* K2)
{
    // copy vars
    fp2_t c1;

    fp2_t t0;
    fp2_subm(&P3->X, &P3->Z, mod, K1);          // 1.  K1 <- X_P3 - Z_P3        K1 set
    fp2_mulm(K1, K1, mod, &t0);                 // 2.  t0 <- K1^2

    fp2_t t1;
    fp2_addm(&P3->X, &P3->Z, mod, K2);          // 3.  K2 <- X_P3 + Z_P3        K2 set
    fp2_mulm(K2, K2, mod, &t1);                 // 4.  t1 <- K2^2

    fp2_t t2;
    fp2_t t3;
    fp2_addm(&t0, &t1, mod, &t2);               // 5.  t2 <- t0 + t1
    fp2_addm(K1, K2, mod, &t3);                 // 6.  t3 <- K1 + K2

    fp2_copy(&t3, &c1);                         // c1 = t3 for 7.
    fp2_mulm(&c1, &c1, mod, &t3);               // 7.  t3 <- t3^2
    fp2_copy(&t3, &c1);                         // c1 = t3 for 8.
    fp2_subm(&c1, &t2, mod, &t3);               // 8.  t3 <- t3 - t2

    fp2_addm(&t1, &t3, mod, &t2);               // 9. t2 <- t1 + t3
    fp2_copy(&t3, &c1);                         // c1 = t3 for 10.
    fp2_addm(&c1, &t0, mod, &t3);               // 10. t3 <- t3 + t0

    fp2_t t4;
    fp2_addm(&t3, &t0, mod, &t4);               // 11. t4 <- t3 + t0
    fp2_copy(&t4, &c1);                         // c1 = t4 for 12.
    fp2_addm(&c1, &c1, mod, &t4);               // 12. t4 <- t4 + t4

    fp2_copy(&t4, &c1);                         // c1 = t4 for 13.
    fp2_addm(&t1, &c1, mod, &t4);               // 13. t4 <- t1 + t4
    fp2_mulm(&t2, &t4, mod, &A24->C);           // 14. A_24^- <- t2 * t4        A_24^- set

    fp2_addm(&t1, &t2, mod, &t4);               // 15. t4 <- t1 + t2
    fp2_copy(&t4, &c1);                         // c1 = t4 for 16.
    fp2_addm(&c1, &c1, mod, &t4);               // 16. t4 <- t4 + t4

    fp2_copy(&t4, &c1);                         // c1 = t4 for 17.
    fp2_addm(&t0, &c1, mod, &t4);               // 17. t4 <- t0 + t4
    fp2_mulm(&t3, &t4, mod, &A24->A);           // 18. A_24^+ <- t3 * t4

}

void iso_3_eval(fp2_t* K1, fp2_t* K2, proj_point_t* Q, fp_t* mod, proj_point_t* Q_strich)
{
    // copy vars
    fp2_t c1;

    fp2_t t0;
    fp2_t t1;
    fp2_addm(&Q->X, &Q->Z, mod, &t0);           // 1.  t0 <- X_Q + Z_Q
    fp2_subm(&Q->X, &Q->Z, mod, &t1);           // 2.  t1 <- X_Q - Z_Q

    fp2_copy(&t0, &c1);                         // c1 = t0 for 3.
    fp2_mulm(K1, &c1, mod, &t0);                // 3.  t0 <- K1 * t0
    fp2_copy(&t1, &c1);                         // c1 = t1 for 4.
    fp2_mulm(K2, &c1, mod, &t1);                // 4. t1 <- K2 * t1

    fp2_t t2;
    fp2_addm(&t0, &t1, mod, &t2);               // 5.  t2 <- t0 + t1
    fp2_copy(&t0, &c1);                         // c1 = t0 for 6.
    fp2_subm(&t1, &c1, mod, &t0);               // 6.  t0 <- t1 - t0

    fp2_copy(&t2, &c1);                         // c1 = t2 for 7.
    fp2_mulm(&c1, &c1, mod, &t2);               // 7.  t2 <- t2^2
    fp2_copy(&t0, &c1);                         // c1 = t0 for 8.
    fp2_mulm(&c1, &c1, mod, &t0);               // 8.  t0 <- t0^2
    fp2_mulm(&Q->X, &t2, mod, &Q_strich->X);    // 9.  X_Q' <- X_Q * t2     X_Q' set
    fp2_mulm(&Q->Z, &t0, mod, &Q_strich->Z);    // 10  Z_Q' <- Z_Q * t0     Z_Q' set
}



void proj_pt_copy_masked(proj_point_t* a, proj_point_t* b, uint64_t mask)
{
    fp2_copy_masked(&a->X, &b->X, mask);
    fp2_copy_masked(&a->Z, &b->Z, mask);
}