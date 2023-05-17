


#ifndef SIDH_H_
#define SIDH_H_

#include "public_params.h"
#include "public_key.h"
#include "montgomory_redc.h"
#include "isogeny_strat.h"


void isogen_2(public_params_t* params, fp_t* secret_key, fp_t* mod, public_key_t* p2, const isogeny_strat_t* strategy_e2);

void isogen_3(public_params_t* params, fp_t* secret_key, fp_t* mod, public_key_t* p3, const isogeny_strat_t* strategy_e3);

void isoex_2(fp_t* secret_key, public_key_t* public_key, fp2_t* j_inv, fp_t* mod, const isogeny_strat_t* strategy_e2);

void isoex_3(fp_t* secret_key, public_key_t* public_key, fp2_t* j_inv, fp_t* mod, const isogeny_strat_t* strategy_e3);


#endif