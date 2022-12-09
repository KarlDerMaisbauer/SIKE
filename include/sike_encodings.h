

#ifndef SIKE_ENCODINGS_H_
#define SIKE_ENCODINGS_H_




#include<stdint.h>

#include"fp.h"
#include"fp2.h"
#include"params.h"
#include"public_params.h"


void ostoi(const char*  value, fp_t* res);

int ostofp(const char* value, fp_t* mod, fp_t* res);

int ostofp2(const char* valuea0, const char* valuea1, fp_t* mod, fp2_t* res);

char* itoos(fp_t* value);

int64_t get_e(const char* value);

void params_translate(const params_t* params_encoded, public_params_t* params_translated);

void params_translate_redec(const params_t* params_encoded, public_params_t* params_translated);
#endif