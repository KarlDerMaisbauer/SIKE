




#include<string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include"sike_encodings_helper.h"


uint64_t c16_to_hex(const char* val, int pos)
{
    char value[] = "0000000000000000";
    for(int i = 15; i >= 0; i--)
    {
        value[i] = val[pos+i];
    }
    return (uint64_t)strtoull(value, 0, 16);
}

uint64_t c_end_to_hex(const char* val, int pos)
{
    char value[] = "0000000000000000";
    for(int i = 15; i >= 0 && pos + i >= 0; i--)
    {
        if(val[pos + i] == 'x')
            break;
        value[i] = val[pos+i];
    }
    return strtoull(value, 0, 16);
}


char itohex(uint8_t value)
{
    switch(value) {
        case 0: return '0';
        case 1: return '1';
        case 2: return '2';
        case 3: return '3';
        case 4: return '4';
        case 5: return '5';
        case 6: return '6';
        case 7: return '7';
        case 8: return '8';
        case 9: return '9';
        case 10: return 'A';
        case 11: return 'B';
        case 12: return 'C';
        case 13: return 'D';
        case 14: return 'E';
        case 15: return 'F';
    }
    return 'R';
}