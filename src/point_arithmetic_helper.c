

#include"point_arithmetic_helper.h"
#include"point_arithmetic.h"
#include"fp_helper.h"


void LADDER_inner_part(proj_point_t* P0, proj_point_t* P1, proj_point_t* P2, mont_curve_t* a24, uint64_t condition, fp_t* mod)
{
    // a serves a substitue for the secont input of xDBLADD
    proj_point_t a_in;

    // a serves a substitue for the  third input of xDBLADD
    proj_point_t b_in;
    proj_point_t b_out;

    uint64_t pos_mask = 0;
    
    for(int i = 0; i < 64; i++)
    {
        pos_mask = pos_mask | (condition << i);
    }
    uint64_t neg_mask = ~pos_mask;

    // second input definition
    proj_pt_copy_masked(P1, &a_in, pos_mask);
    proj_pt_copy_masked(P2, &a_in, neg_mask);

    // third input definition
    proj_pt_copy_masked(P2, &b_in, pos_mask);
    proj_pt_copy_masked(P1, &b_in, neg_mask);


    xDBLADD(P0, &a_in, &b_in, mod, a24, P0, &b_out);

    // define output
    proj_pt_copy_masked(&b_out, P1, pos_mask);
    proj_pt_copy_masked(&b_out, P2, neg_mask);
}