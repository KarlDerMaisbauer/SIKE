

#include<string.h>

#include"sike_encodings.h"
#include"sike_encodings_helper.h"
#include"fp_helper.h"

#include <stdio.h>
#include <stdlib.h>


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
