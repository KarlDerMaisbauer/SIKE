




#include"sike_encodings_helper.h"




uint64_t hex16_to_int(char* str, int16_t pos)
{
    int16_t stop = pos + 16;
    int8_t ctr = 0;
    char* val = "0000000000000000";
    for(; pos > pos_new; pos++)
    {
        val[ctr] = str[pos];
        ctr++
    }
    return strtol(val, 0, 16);
}