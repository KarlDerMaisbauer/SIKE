

#include"montgomery_curve.h"
#include"montgomory_redc.h"

#include"fp2.h"
#include"fp.h"
#include"fp_helper.h"
#include"isogeny_strat.h"
#include"point_arithmetic.h"

void jinvariant(fp2_t* A, fp2_t* C, fp_t* mod, fp2_t* j)
{
    // copy vars
    fp2_t c1;

    fp2_t t1;
    fp2_mul_mont(A, A, mod, j);                 // 1.  j <- A^2
    fp2_mul_mont(C, C, mod, &t1);               // 2.  t1 <- C^2

    fp2_t t0;
    fp2_add(&t1, &t1, &t0);                     // 3.  t0 <- t1 + t1
    fp2_copy(&t0, &c1);                         // c1 = t0 for 4.
    fp2_sub(j, &c1, &t0);                       // 4.  t0 <- j - t0

    fp2_copy(&t0, &c1);                         // c1 = t0 for 5.
    fp2_sub(&c1, &t1, &t0);                     // 5.  t0 <- t0 - t1
    fp2_sub(&t0, &t1, j);                       // 6.  j <- t0 - t1

    fp2_copy(&t1, &c1);                         // c1 = t1 for 7.
    fp2_mul_mont(&c1, &c1, mod, &t1);           // 7.  t1 <- t1^2
    fp2_copy(j, &c1);                           // c1 = j for 8.
    fp2_mul_mont(&c1, &t1, mod, j);             // 8.  j <- j * t1

    fp2_copy(&t0, &c1);                         // c1 = t0 for 9.
    fp2_add(&c1, &c1, &t0);                     // 9.  t0 <- t0 + t0
    fp2_copy(&t0, &c1);                         //c1 = t0 for 10.
    fp2_add(&c1, &c1, &t0);                     // 10. t0 <- t0 + t0

    fp2_mul_mont(&t0, &t0, mod, &t1);           // 11. t1 <- t0^2
    fp2_copy(&t0, &c1);                         // c1 = t0 for 12.
    fp2_mul_mont(&c1, &t1, mod, &t0);           // 12. t0 <- t0 * t1

    fp2_copy(&t0, &c1);                         // c1 = t0 for 13.
    fp2_add(&c1, &c1, &t0);                     // 13. t0 <- t0 + t0
    fp2_copy(&t0, &c1);                         // c1 = t0 for 14.
    fp2_add(&c1, &c1, &t0);                     // 14. t0 <- t0 + t0

    fp2_copy(j, &c1);                           // c1 = j for 15.
    fp2_mult_inv(&c1, mod, j);                  // 15. j <- 1/j
    fp2_copy(j, &c1);                           // c1 = j for 16.
    fp2_mul_mont(&c1, &t0, mod, j);             // 16. j <- j * t0
}

void get_A(fp2_t* xp, fp2_t* xq, fp2_t* xqp, fp_t* mod, fp2_t* A)
{
    // copy vars
    fp2_t c1;


    fp2_t t0;
    fp2_t t1;
    fp2_add(xp, xq, &t1);                       // 1.  t1 <- xp + xq
    fp2_mul_mont(xp, xq, mod, &t0);             // 2.  t0 <- xp * xq

    fp2_mul_mont(xqp, &t1, mod, A);             // 3.  A <- xqp * t1
    fp2_copy(A, &c1);                           // c1 = A for 4.
    fp2_add(&c1, &t0, A);                       // 4.  A <- A + t0

    fp2_copy(&t0, &c1);                         // c1 = t0 for 5.
    fp2_mul_mont(&c1, xqp, mod, &t0);           // 5.  t0 <- t0 * xqp
    fp2_copy(A, &c1);                           // c1 = A for 6.
    fp2_sub(&c1, &one_mont, A);                      // 6. A <- A - 1

    fp2_copy(&t0, &c1);                         // c1 = t0 for 7.
    fp2_add(&c1, &c1, &t0);                     // 7.  t0 <- t0 + t0
    fp2_copy(&t1, &c1);                         // c1 = t1 for 8.
    fp2_add(&c1, xqp, &t1);                     // 8.  t1 <- t1 + xqp

    fp2_copy(&t0 , &c1);                        // c1 = t0 for 9.
    fp2_add(&c1, &c1, &t0);                     // 9.  t0 <- t0 + t0
    fp2_copy(A, &c1);                           // c1 = A for 10.
    fp2_mul_mont(&c1, &c1, mod, A);             // 10. A <- A^2

    fp2_copy(&t0, &c1);                         // c1 = t0 for 11.
    fp2_mult_inv(&c1, mod, &t0);                // 11. t0 <- 1/t0
    fp2_copy(A, &c1);                           // c1 = A for 12.
    fp2_mul_mont(&c1, &t0, mod, A);             // A <- A * t0

    fp2_copy(A, &c1);                           // c1 = A for 13.
    fp2_sub(&c1, &t1, A);                       // A <- A - t1
}

void iso_2_curve(proj_point_t* P2, fp_t* mod, mont_curve_t* A24)
{
    // copy vars
    fp2_t c1;

    fp2_mul_mont(&P2->X, &P2->X, mod, &A24->A); // 1.  A_24^+ <- X_P2^2
    fp2_mul_mont(&P2->Z, &P2->Z, mod, &A24->C); // 2.  C_24^+ <- Z_P2^2

    fp2_copy(&A24->A, &c1);                     // c1 = A_24^+ for 3.
    fp2_sub(&A24->C, &c1, &A24->A);             // 3. A_24^+ <- C_24 - A_24^+
}

void iso_2_eval(proj_point_t* P2,proj_point_t* Q, fp_t* mod, proj_point_t* Q_strich)
{
    // copy vars
    fp2_t c1;

    fp2_t t0;
    fp2_t t1;
    fp2_add(&P2->X, &P2->Z, &t0);               // 1.  t0 <- X_p2 + Z_P2
    fp2_sub(&P2->X, &P2->Z, &t1);               // 2.  t1 <- X_p2 - Z_P2


    fp2_t t2;
    fp2_t t3;
    fp2_add(&Q->X, &Q->Z, &t2);                 // 3.  t2 <- X_Q + Z_Q
    fp2_sub(&Q->X, &Q->Z, &t3);                 // 4.  t3 <- X_Q - Z_Q

    fp2_copy(&t0, &c1);                         // c1 = t0 for 5.
    fp2_mul_mont(&c1, &t3, mod, &t0);           // 5.  t0 <- t0 * t3
    fp2_copy(&t1, &c1);                         // c1 = t1 for 6.
    fp2_mul_mont(&c1, &t2, mod, &t1);           // 6.  t1  <- t1 * t2

    fp2_add(&t0, &t1, &t2);                     // 7.  t2 <- t0 + t1   
    fp2_sub(&t0, &t1, &t3);                     // 8.  t3 <- t0 - t1   

    fp2_mul_mont(&Q->X, &t2, mod, &Q_strich->X);// 9.  X_Q' <- X_Q * t2
    fp2_mul_mont(&Q->Z, &t3, mod, &Q_strich->Z);// 10. Z_Q' <- Z_Q * t3
}

void iso_4_curve(proj_point_t* P4, fp_t* mod, mont_curve_t* A24, fp2_t* K1, fp2_t* K2, fp2_t* K3)
{
    // copy vars
    fp2_t c1;

    fp2_sub(&P4->X, &P4->Z, K2);                // 1.  K_2 <- X_P4 - Z_P4   K_2 set
    fp2_add(&P4->X, &P4->Z, K3);                // 2.  K_3 <- X_P4 + Z_P4   K_3 set

    fp2_mul_mont(&P4->Z, &P4->Z, mod, K1);      // 3.  K_1 <- Z_p4^2
    fp2_copy(K1, &c1);                          // c1 = K1 for 4.
    fp2_add(&c1, &c1, K1);                      // 4.  K1  <- K1 + K1

    fp2_mul_mont(K1, K1, mod, &A24->C);         // 5.  C24 <- K1^2          C24 set
    fp2_copy(K1, &c1);                          // c1 = K1 for 6.
    fp2_add(&c1, &c1, K1);                      // 6.  K1 <- K1 + K1        K1 set

    fp2_mul_mont(&P4->X, &P4->X, mod, &A24->A); // 7.  A_24 <- X_P4^2
    fp2_copy(&A24->A, &c1);                     // c1 = A24 for 8.
    fp2_add(&c1, &c1, &A24->A);                 // 8.  A_24 <- A_24 + A_24

    fp2_copy(&A24->A, &c1);                     // c1 = A24 for 9.
    fp2_mul_mont(&c1, &c1, mod, &A24->A);       // 9.  A24 <- A24^2
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
    fp2_add(&Q->X, &Q->Z, &t0);                 // 1.  t0 <- X_Q + Z_Q
    fp2_sub(&Q->X, &Q->Z, &t1);                 // 2.  t0 <- X_Q - Z_Q

    fp2_mul_mont(&t0, K2, mod, &X_Q);           // 3.  X_Q <- t0 * K_2
    fp2_mul_mont(&t1, K3, mod, &Z_Q);           // 4.  Z_Q <- t1 * K_3

    fp2_copy(&t0, &c1);                         // c1 = t0 for 5.
    fp2_mul_mont(&c1, &t1, mod, &t0);           // 5.  t0 <- t0 * t1
    fp2_copy(&t0, &c1);                         // c1 = t0 for 6.
    fp2_mul_mont(&c1, K1, mod, &t0);            // 6.  t0 <- t0 * K1

    fp2_add(&X_Q, &Z_Q, &t1);                   // 7.  t1 <- X_Q + Z_Q
    fp2_copy(&Z_Q, &c1);                        // c1 = Z_Q for 8.
    fp2_sub(&X_Q, &c1, &Z_Q);                   // 8.  Z_Q <- X_Q - Z_Q

    fp2_copy(&t1, &c1);                         // c1 = t1 for 9.
    fp2_mul_mont(&c1, &c1, mod, &t1);           // 9.  t1 <- t1^2
    fp2_copy(&Z_Q, &c1);                        // c1 = Z_Q for 10.
    fp2_mul_mont(&c1, &c1, mod, &Z_Q);          // 10. Z_Q <- Z_Q^2

    fp2_add(&t0, &t1, &X_Q);                    // 11. X_Q <- t0 + t1
    fp2_copy(&t0, &c1);                         // c1 = t0 for 12.
    fp2_sub(&Z_Q, &c1, &t0);                    // 12. t0 <- Z_Q - t0

    fp2_mul_mont(&X_Q, &t1, mod, &Q_strich->X); // 13. X_Q' <- X_Q * t1
    fp2_mul_mont(&Z_Q, &t0, mod, &Q_strich->Z); // 14. Z_Q' <- Z_Q * t0
}

void iso_3_curve(proj_point_t* P3, fp_t* mod, mont_curve_t* A24, fp2_t* K1, fp2_t* K2)
{
    // copy vars
    fp2_t c1;

    fp2_t t0;
    fp2_sub(&P3->X, &P3->Z, K1);                // 1.  K1 <- X_P3 - Z_P3        K1 set
    fp2_mul_mont(K1, K1, mod, &t0);             // 2.  t0 <- K1^2

    fp2_t t1;
    fp2_add(&P3->X, &P3->Z, K2);                // 3.  K2 <- X_P3 + Z_P3        K2 set
    fp2_mul_mont(K2, K2, mod, &t1);             // 4.  t1 <- K2^2

    fp2_t t2;
    fp2_t t3;
    fp2_add(&t0, &t1, &t2);                     // 5.  t2 <- t0 + t1
    fp2_add(K1, K2, &t3);                       // 6.  t3 <- K1 + K2

    fp2_copy(&t3, &c1);                         // c1 = t3 for 7.
    fp2_mul_mont(&c1, &c1, mod, &t3);           // 7.  t3 <- t3^2
    fp2_copy(&t3, &c1);                         // c1 = t3 for 8.
    fp2_sub(&c1, &t2, &t3);                     // 8.  t3 <- t3 - t2

    fp2_add(&t1, &t3, &t2);                     // 9. t2 <- t1 + t3
    fp2_copy(&t3, &c1);                         // c1 = t3 for 10.
    fp2_add(&c1, &t0, &t3);                     // 10. t3 <- t3 + t0

    fp2_t t4;
    fp2_add(&t3, &t0, &t4);                     // 11. t4 <- t3 + t0
    fp2_copy(&t4, &c1);                         // c1 = t4 for 12.
    fp2_add(&c1, &c1, &t4);                     // 12. t4 <- t4 + t4

    fp2_copy(&t4, &c1);                         // c1 = t4 for 13.
    fp2_add(&t1, &c1, &t4);                     // 13. t4 <- t1 + t4
    fp2_mul_mont(&t2, &t4, mod, &A24->C);       // 14. A_24^- <- t2 * t4        A_24^- set

    fp2_add(&t1, &t2, &t4);                     // 15. t4 <- t1 + t2
    fp2_copy(&t4, &c1);                         // c1 = t4 for 16.
    fp2_add(&c1, &c1, &t4);                     // 16. t4 <- t4 + t4

    fp2_copy(&t4, &c1);                         // c1 = t4 for 17.
    fp2_add(&t0, &c1, &t4);                     // 17. t4 <- t0 + t4
    fp2_mul_mont(&t3, &t4, mod, &A24->A);       // 18. A_24^+ <- t3 * t4

}

void iso_3_eval(fp2_t* K1, fp2_t* K2, proj_point_t* Q, fp_t* mod, proj_point_t* Q_strich)
{
    // copy vars
    fp2_t c1;

    fp2_t t0;
    fp2_t t1;
    fp2_add(&Q->X, &Q->Z, &t0);                 // 1.  t0 <- X_Q + Z_Q
    fp2_sub(&Q->X, &Q->Z, &t1);                 // 2.  t1 <- X_Q - Z_Q

    fp2_copy(&t0, &c1);                         // c1 = t0 for 3.
    fp2_mul_mont(K1, &c1, mod, &t0);            // 3.  t0 <- K1 * t0
    fp2_copy(&t1, &c1);                         // c1 = t1 for 4.
    fp2_mul_mont(K2, &c1, mod, &t1);            // 4. t1 <- K2 * t1

    fp2_t t2;
    fp2_add(&t0, &t1, &t2);                     // 5.  t2 <- t0 + t1
    fp2_copy(&t0, &c1);                         // c1 = t0 for 6.
    fp2_sub(&t1, &c1, &t0);                     // 6.  t0 <- t1 - t0

    fp2_copy(&t2, &c1);                         // c1 = t2 for 7.
    fp2_mul_mont(&c1, &c1, mod, &t2);           // 7.  t2 <- t2^2
    fp2_copy(&t0, &c1);                         // c1 = t0 for 8.
    fp2_mul_mont(&c1, &c1, mod, &t0);           // 8.  t0 <- t0^2
    fp2_mul_mont(&Q->X, &t2, mod, &Q_strich->X);// 9.  X_Q' <- X_Q * t2     X_Q' set
    fp2_mul_mont(&Q->Z, &t0, mod, &Q_strich->Z);// 10  Z_Q' <- Z_Q * t0     Z_Q' set
}







void e_2_iso(mont_curve_t* A24p, proj_point_t* S, proj_point_t* opt_input, int64_t opt_input_size, fp_t* mod, isogeny_strat_t strategy ,mont_curve_t* A24p_new)
{
    proj_point_t queue[strategy.queue_len];
    int8_t next_free_pos = 0;
    int64_t strat_size = strategy.strat_len;

    curve_copy(A24p, A24p_new);
    // Push first Point S
    queue[next_free_pos] = *S;
    next_free_pos++;


    for(int64_t i = 0; i < strat_size; i++)
    {
        int64_t curr_strat_val = strategy.strategy[i];
        if(curr_strat_val == 0)
        {
            fp2_t K1;
            fp2_t K2;
            fp2_t K3;

            // simulate take point out of queue !!!
            next_free_pos--;
            iso_4_curve(&(queue[next_free_pos-1]), mod, A24p_new, &K1, &K2, &K3);

            //update queue
            for(int64_t pos = 0; pos < next_free_pos; pos++)
            {
                proj_point_t update_pt = queue[pos];
                iso_4_eval(&K1, &K2, &K3, &update_pt, mod, &(queue[pos]));
            }
            
            // update optional input
            for(int64_t pos = 0; pos < opt_input_size; pos++)
            {
               proj_point_t update_pt = opt_input[pos]; 
               iso_4_eval(&K1, &K2, &K3, &update_pt, mod, &(opt_input[pos]));
            }
        }
        else
        {
            proj_point_t new_val;
            xDBLe(&(queue[next_free_pos-1]), A24p, 2*curr_strat_val, mod, &new_val);

            // add new point
            queue[next_free_pos] = new_val;
            next_free_pos++;
        }
    }
}

void e_3_iso(mont_curve_t* A24, proj_point_t* S, proj_point_t* opt_input, int64_t opt_input_size, fp_t* mod, isogeny_strat_t strategy ,mont_curve_t* A24_new)
{
    proj_point_t queue[strategy.queue_len];
    int8_t next_free_pos = 0;
    int64_t strat_size = strategy.strat_len;

    curve_copy(A24, A24_new);
    // Push first Point S
    queue[next_free_pos] = *S;
    next_free_pos++;


    for(int64_t i = 0; i < strat_size; i++)
    {
        int64_t curr_strat_val = strategy.strategy[i];
        if(curr_strat_val == 0)
        {
            fp2_t K1;
            fp2_t K2;

            // simulate take point out of queue !!!
            next_free_pos--;
            iso_3_curve(&(queue[next_free_pos-1]), mod, A24_new, &K1, &K2);

            //update queue
            for(int64_t pos = 0; pos < next_free_pos; pos++)
            {
                proj_point_t update_pt = queue[pos];
                iso_3_eval(&K1, &K2, &update_pt, mod, &(queue[pos]));
            }
            
            // update optional input
            for(int64_t pos = 0; pos < opt_input_size; pos++)
            {
               proj_point_t update_pt = opt_input[pos]; 
               iso_3_eval(&K1, &K2, &update_pt, mod, &(opt_input[pos]));
            }
        }
        else
        {
            proj_point_t new_val;
            xTPLe(&(queue[next_free_pos-1]), A24, curr_strat_val, mod, &new_val);

            // add new point
            queue[next_free_pos] = new_val;
            next_free_pos++;
        }
    }
}







void proj_pt_copy_masked(proj_point_t* a, proj_point_t* b, uint64_t mask)
{
    fp2_copy_masked(&a->X, &b->X, mask);
    fp2_copy_masked(&a->Z, &b->Z, mask);
}


void curve_copy(mont_curve_t* a, mont_curve_t* b)
{
    fp2_copy(&a->A, &b->A);
    fp2_copy(&a->C, &b->C);
}

