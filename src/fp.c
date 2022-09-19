

#include"fp.h"
#include"fp_helper.h"
#include"montgomory_redc.h"
#include"sike_encodings.h"


void fp_init(const char* val, fp_t* mod, fp_t* a)
{
    ostofp(val, mod, a);
}

void fp_addm(fp_t* a, fp_t* b, fp_t* mod, fp_t* res)
{
    fp_t res_part;
    fp_add(a, b, &res_part);
    fp_mod(&res_part, mod, res);
}


void fp_subm(fp_t* a, fp_t* b, fp_t*mod, fp_t *res)
{
    fp_t res_part;
    fp_sub(a, b, &res_part);
    fp_mod(&res_part, mod, res);
}


void fp_mulm(fp_t* a, fp_t* b, fp_t*mod, fp_t* res)
{
    MODMUL(a, b, mod, res);
}

void fp_divm(fp_t* a, fp_t* b, fp_t*mod, fp_t* res)
{
    fp_t mod_inv;
    fp_t b_inv;
    gcd(b, mod, &b_inv, &mod_inv);
    f2p_t ab;
    fp_mul(a, &b_inv, &ab);
    f2p_mod(&ab, mod, res);
}