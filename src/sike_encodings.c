

#include<string.h>

#include"sike_encodings.h"
#include"sike_encodings_helper.h"
#include"fp_helper.h"
#include"params.h"
#include"public_params.h"
#include"montgomory_redc.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


void ostoi(const char*  value, fp_t* res)
{
    int length = strlen(value) - 16;
    int ctr = 0;
    while(length > 0)
    {
        (*res)[ctr] = c16_to_hex(value, length);
        length -= 16;
        ctr++;
    }
    (*res)[ctr] = c_end_to_hex(value, length);
}

int ostofp(const char* value, fp_t* mod, fp_t* res)
{
    ostoi(value, res);
    return fp_smaller_zero(res) || fp_greater_equ_pos(res, mod);
}

int ostofp2(const char* valuea0, const char* valuea1, fp_t* mod, fp2_t* res)
{
    if(ostofp(valuea0, mod, &(res->real)))
        return 1;
    if(ostofp(valuea1, mod, &(res->img)))
        return 1;
    return 0;
}

char* itoos(fp_t* value)
{
    uint8_t mask = 0xF;
    int startpoint = 0;
    int start_pos = 0;
    char* number;
    for(int i = WORDS - 1; i >= 0; i--)
    {
        if((*value)[i] == 0)
            continue;
        for(int pos = 15; pos >= 0; pos--)
        {
            uint8_t hexval = ((*value)[i] >> (4*pos)) & mask;
            if(hexval != 0)
            {
                int size = (i*16 + pos + 3)*sizeof(char);
                number = (char*)malloc(size);
                startpoint = i;
                start_pos = pos;
                number[0] = '0';
                number[1] = 'x';
                number[size-1] = '\0';
                break;

            }
        }
        break;
    }
    
    int ctr = 2;
    for(int i = start_pos; i >= 0; i--)
    {
        uint8_t hexval = ((*value)[startpoint] >> (4*i)) & mask;
        number[ctr] = itohex(hexval);
        ctr++;

    }
    for(int i = startpoint-1; i >= 0; i--)
    {
        for(int hexpos = 0; hexpos < 16; hexpos++)
        {
            uint8_t hexval = ((*value)[i] >> (4*(15-hexpos))) & mask;
            int position = 16*((startpoint-1)-i)+hexpos+ctr;
            number[position] = itohex(hexval);
        }
    }
    return number;
}

int64_t get_e(const char* value)
{
    return strtoll(value, 0, 16);
}

void params_translate(const params_t* params_encoded, public_params_t* params_translated)
{
    //zero the translated params
    fp_zero(&params_translated->e2);
    fp_zero(&params_translated->e3);
    fp_zero(&params_translated->p);
    fp2_zero(&params_translated->xP2);
    fp2_zero(&params_translated->xQ2);
    fp2_zero(&params_translated->xR2);
    fp2_zero(&params_translated->xP3);
    fp2_zero(&params_translated->xQ3);
    fp2_zero(&params_translated->xR3);

    fp_t modtemp;
    fp_zero(&modtemp);
    ostoi(params_encoded->p, &modtemp);
    assert(!fp_smaller_zero(&modtemp));
    fp_copy(&modtemp, &params_translated->p);


    // e2
    ostofp(params_encoded->e2, &modtemp, &params_translated->e2);

    //e3
    ostofp(params_encoded->e3, &modtemp, &params_translated->e3);

    // xQ2
    ostofp(params_encoded->xQ20, &modtemp, &params_translated->xQ2.real);
    ostofp(params_encoded->xQ21, &modtemp, &params_translated->xQ2.img);

    // xP2
    ostofp(params_encoded->xP20, &modtemp, &params_translated->xP2.real);
    ostofp(params_encoded->xP21, &modtemp, &params_translated->xP2.img);

    // xR2
    ostofp(params_encoded->xR20, &modtemp, &params_translated->xR2.real);
    ostofp(params_encoded->xR21, &modtemp, &params_translated->xR2.img);


    // xQ3
    ostofp(params_encoded->xQ30, &modtemp, &params_translated->xQ3.real);
    ostofp(params_encoded->xQ31, &modtemp, &params_translated->xQ3.img);

    // xP3
    ostofp(params_encoded->xP30, &modtemp, &params_translated->xP3.real);
    ostofp(params_encoded->xP31, &modtemp, &params_translated->xP3.img);

    // xR3
    ostofp(params_encoded->xR30, &modtemp, &params_translated->xR3.real);
    ostofp(params_encoded->xR31, &modtemp, &params_translated->xR3.img);


}

void params_translate_redec(const params_t* params_encoded, public_params_t* params_translated)
{
    fp_zero(&params_translated->e2);
    fp_zero(&params_translated->e3);
    fp_zero(&params_translated->p);
    fp2_zero(&params_translated->xP2);
    fp2_zero(&params_translated->xQ2);
    fp2_zero(&params_translated->xR2);
    fp2_zero(&params_translated->xP3);
    fp2_zero(&params_translated->xQ3);
    fp2_zero(&params_translated->xR3);

    fp_t xrtemp;
    fp_t xitemp;
    fp_t modtemp;
    fp_zero(&xrtemp);
    fp_zero(&xitemp);
    fp_zero(&modtemp);
    ostoi(params_encoded->p, &modtemp);
    assert(!fp_smaller_zero(&modtemp));
    fp_copy(&modtemp, &params_translated->p);
    init(&params_translated->p);

    fp2_t temp;
    fp2_zero(&temp);
    // do they need ranslation?
    // e2
    ostofp(params_encoded->e2, &modtemp, &params_translated->e2);

    //e3
    ostofp(params_encoded->e3, &modtemp, &params_translated->e3);

    // xQ2
    ostofp(params_encoded->xQ20, &modtemp, &temp.real);
    ostofp(params_encoded->xQ21, &modtemp, &temp.img);

    fp2_to_mont(&temp, &modtemp, &params_translated->xQ2);
    fp2_zero(&temp);

    // xP2
    ostofp(params_encoded->xP20, &modtemp, &temp.real);
    ostofp(params_encoded->xP21, &modtemp, &temp.img);

    fp2_to_mont(&temp, &modtemp, &params_translated->xP2);
    fp2_zero(&temp);

    // xR2
    ostofp(params_encoded->xR20, &modtemp, &temp.real);
    ostofp(params_encoded->xR21, &modtemp, &temp.img);

    fp2_to_mont(&temp, &modtemp, &params_translated->xR2);
    fp2_zero(&temp);


    // xQ3
    ostofp(params_encoded->xQ30, &modtemp, &temp.real);
    ostofp(params_encoded->xQ31, &modtemp, &temp.img);

    fp2_to_mont(&temp, &modtemp, &params_translated->xQ3);
    fp2_zero(&temp);

    // xP3
    ostofp(params_encoded->xP30, &modtemp, &temp.real);
    ostofp(params_encoded->xP31, &modtemp, &temp.img);

    fp2_to_mont(&temp, &modtemp, &params_translated->xP3);
    fp2_zero(&temp);

    // xR3
    ostofp(params_encoded->xR30, &modtemp, &temp.real);
    ostofp(params_encoded->xR31, &modtemp, &temp.img);

    fp2_to_mont(&temp, &modtemp, &params_translated->xR3);
    fp2_zero(&temp);
}
