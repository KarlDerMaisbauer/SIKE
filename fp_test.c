

#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include "fp.h"
#include "fp_helper.h"
#include "fp2.h"
#include "montgomory_redc.h"
#include "montgomery_curve.h"

#include "sike_encodings_helper.h"
#include "sike_encodings.h"
#include "params.h"
#include "point_arithmetic.h"
#include "isogeny_strat.h"
#include "public_params.h"
#include "public_key.h"
#include "sidh.h"


#include <pthread.h>


#define NUM_TESTS 40
#define NUM_TESTS_COSTLY 10
#define NUM_LADDER_TESTS 5
#define NUM_E_ISO_TESTS 1

#define S434 434
#define S503 503
#define S610 610
#define S751 751

#define MODE S610

#define BITS MODE


#define TRUE 1
#define FALSE 0


#define FPTESTS FALSE
#define FP2TESTS FALSE
#define REDCTESTS FALSE
#define POINTARTESTS TRUE
#define XDBLE FALSE
#define XTPLE FALSE
#define LADDER3P TRUE
#define MONTCURVETESTS FALSE
#define MONTCURVETESTSREST FALSE
#define E_2_ISO FALSE
#define E_3_ISO FALSE
#define KEYEXCHANGE FALSE 

#if MODE == S434
    #define STRAT_E2 p434_e2
    #define STRAT_E3 p434_e3
    #define PARAMS SIKEp434
#elif MODE == S503
    #define STRAT_E2 p503_e2
    #define STRAT_E3 p503_e3
    #define PARAMS SIKEp503
#elif MODE == S610
    #define STRAT_E2 p610_e2
    #define STRAT_E3 p610_e3
    #define PARAMS SIKEp610
#elif MODE == S751
    #define STRAT_E2 p751_e2
    #define STRAT_E3 p751_e3
    #define PARAMS SIKEp751
#endif

void rand_fill_fp_pos(fp_t* a);
void rand_fill_fp_neg(fp_t* a);
void rand_fill_fp_2_input(fp_t* a, fp_t* b, int i);
void rand_fill_fp_divisor(fp_t* a);
void rand_fill_fp_divisor_pos_neg(fp_t* a, int i);
void rand_2k_num(fp_t* a);
void rand_fp_smaller_mod(fp_t* a, fp_t* mod);



void rand_fill_f2p_pos(f2p_t* a);
void rand_fill_f2p_neg(f2p_t* a);
void rand_fill_f2p_2_input(f2p_t* a, f2p_t* b, int i);
void rand_fill_f2p_divisor(f2p_t* a);
void rand_2k_num_f2p(f2p_t* a);
void rand_fill_f2p_divisor_pos_neg(f2p_t* a, int i);
void rand_fp2_smaller_mod(fp2_t* a, fp_t* mod);


void rand_fill_fp2(fp2_t* a);
void rand_fill_fp2_2_input(fp2_t* a, fp2_t* b, int i);


void rand_fill_fp_for_f2p_pos(f2p_t* a);
void rand_fill_fp_for_f2p_neg(f2p_t* a);
void rand_fill_fp_for_f2p_2_input(f2p_t* a, f2p_t* b, int i);

// projpoint fkts
void rand_fill_proj_pt(proj_point_t* p, fp_t* mod);


void* alice_pub_calc(void* params);

void* bob_pub_calc(void* params);

void* alice_sec_calc(void* params);

void* bob_sec_calc(void* params);


typedef struct{
    fp_t* secret;
    public_params_t* parameters;
    public_key_t* calced_key;
}pub_calc_t;

typedef struct{
    fp_t* secret;
    public_params_t* parameters;
    public_key_t* pub_key;
    fp2_t* jinv;
}sec_calc_t;


int main(void){
    printf("Mode %d\n", MODE);
    //srand(time(NULL));
    if(!KEYEXCHANGE)
    {
        //srand(time(NULL));
        fp_t mont_mod;
        fp_zero(&mont_mod);
        if(MODE == S434)
        {

            mont_mod[0] = 0xFFFFFFFFFFFFFFFF;
            mont_mod[1] = 0xFFFFFFFFFFFFFFFF;
            mont_mod[2] = 0xFFFFFFFFFFFFFFFF;
            mont_mod[3] = 0xFDC1767AE2FFFFFF;
            mont_mod[4] = 0x7BC65C783158AEA3;
            mont_mod[5] = 0x6CFC5FD681C52056;
            mont_mod[6] = 0x0002341F27177344;
        }
        else if(MODE == S503)
        {
            
            mont_mod[0] = 0xFFFFFFFFFFFFFFFF;
            mont_mod[1] = 0xFFFFFFFFFFFFFFFF;
            mont_mod[2] = 0xFFFFFFFFFFFFFFFF;
            mont_mod[3] = 0xABFFFFFFFFFFFFFF;
            mont_mod[4] = 0x13085BDA2211E7A0;
            mont_mod[5] = 0x1B9BF6C87B7E7DAF;
            mont_mod[6] = 0x6045C6BDDA77A4D0;
            mont_mod[7] = 0x004066F541811E1E;
        }
        else if(MODE == S610)
        {
                
            /*1e5cdaf9ce34d6910574427c2119c2d27005a09e1531d1b4e7690233f1c6d5c59a48d2cf7b6df187c561c02b6c0a137ee291942c0caa527720cf83aa0e7ae
            mont_mod[0] = 0xFFFFFFFFFFFFFFFF;
            mont_mod[1] = 0xFFFFFFFFFFFFFFFF;
            mont_mod[2] = 0xFFFFFFFFFFFFFFFF;
            mont_mod[3] = 0xFFFFFFFFFFFFFFFF;
            mont_mod[4] = 0x6E01FFFFFFFFFFFF;
            mont_mod[5] = 0xB1784DE8AA5AB02E;
            mont_mod[6] = 0x9AE7BF45048FF9AB;
            mont_mod[7] = 0xB255B2FA10C4252A;
            mont_mod[8] = 0x819010C251E7D88C;
            mont_mod[9] = 0x000000027BF6A768;*/
            
            public_params_t params_simple;
            params_translate(&PARAMS, &params_simple);
            fp_copy(&params_simple.p, &mont_mod);
            fp_print(&mont_mod);
            LINE

        
        }
        else if(MODE == S751)
        {
            mont_mod[0]  = 0xFFFFFFFFFFFFFFFF;
            mont_mod[1]  = 0xFFFFFFFFFFFFFFFF;
            mont_mod[2]  = 0xFFFFFFFFFFFFFFFF;
            mont_mod[3]  = 0xFFFFFFFFFFFFFFFF;
            mont_mod[4]  = 0xFFFFFFFFFFFFFFFF;
            mont_mod[5]  = 0xEEAFFFFFFFFFFFFF;
            mont_mod[6]  = 0xE3EC968549F878A8;
            mont_mod[7]  = 0xDA959B1A13F7CC76;
            mont_mod[8]  = 0x084E9867D6EBE876;
            mont_mod[9]  = 0x8562B5045CB25748;
            mont_mod[10] = 0x0E12909F97BADC66;
            mont_mod[11] = 0x00006FE5D541F71C;
        }
        else
        {
            printf("No valid mode\n");
            return -1;
        }
        init(&mont_mod);
        printf("Mode %d\n", MODE);
        //init(&mont_mod);

        FILE* file = fopen("erg.txt", "w+");
        fprintf(file, "%d\n", MODE);

        if(FPTESTS)
        {
            // test fp_Add
            printf("Generating Testdata for fp_add\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp_t add1;
                fp_t add2;
                fp_t res;
                fp_zero(&add1);
                fp_zero(&add2);
                fp_zero(&res);
                rand_fill_fp_2_input(&add1, &add2, i);
                fp_add(&add1, &add2, &res);
                fprintf(file, "+\n");
                fp_fprintf(file, &add1);
                fp_fprintf(file, &add2);
                fp_fprintf(file, &res);
            }


            // test fp_sub
            printf("Generating Testdata for fp_sub\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp_t add1;
                fp_t add2;
                fp_t res;

                fp_zero(&add1);
                fp_zero(&add2);
                fp_zero(&res);
                rand_fill_fp_2_input(&add1, &add2, i);
                fp_sub(&add1, &add2, &res);
                fprintf(file, "-\n");
                fp_fprintf(file, &add1);
                fp_fprintf(file, &add2);
                fp_fprintf(file, &res);
            }


            // test f2p_Add
            printf("Generating Testdata for f2p_add\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                f2p_t add1;
                f2p_t add2;
                f2p_t res;

                f2p_zero(&add1);
                f2p_zero(&add2);
                f2p_zero(&res);
                rand_fill_f2p_2_input(&add1, &add2, i);
                f2p_add(&add1, &add2, &res);
                fprintf(file, "++\n");
                f2p_fprintf(file, &add1);
                f2p_fprintf(file, &add2);
                f2p_fprintf(file, &res);
            }

            // test f2p_sub
            printf("Generating Testdata for f2p_sub\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                f2p_t add1;
                f2p_t add2;
                f2p_t res;

                f2p_zero(&add1);
                f2p_zero(&add2);
                f2p_zero(&res);
                rand_fill_f2p_2_input(&add1, &add2, i);
                f2p_sub(&add1, &add2, &res);
                fprintf(file, "--\n");
                f2p_fprintf(file, &add1);
                f2p_fprintf(file, &add2);
                f2p_fprintf(file, &res);
            }

            //test for fp_mul
            printf("Generating Testdata for fp_mul\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp_t add1;
                fp_t add2;
                f2p_t res;

                fp_zero(&add1);
                fp_zero(&add2);
                f2p_zero(&res);
                //rand_fill_fp_pos(&add1);
                rand_fill_fp_2_input(&add1, &add2, i);
                //add1[1] = 0x345;
                //add1[0] = 0x8765434;
                //add2[2] = 87654;
                //add2[0] = 876543;
                
                fp_mul(&add1, &add2, &res);
                //uint64_t add3 = (uint64_t)((uint32_t)rand());
                //add3 = add3;
                
                //fp_mul_small(&add1, add3, &res);
                fprintf(file, "*\n");
                fp_fprintf(file, &add1);
                fp_fprintf(file, &add2);
                //fprintf(file, "0x0000000000%lx\n", add3);
            
                f2p_fprintf(file, &res);
                //f2p_print(&res);
                //LINE
            }
            //test for fp_div
            printf("Generating Testdata for fp_div\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp_t add1;
                fp_t add2;
                fp_t res;

                fp_zero(&add1);
                fp_zero(&add2);
                fp_zero(&res);
                if((i % 4)>= 2)
                {
                    rand_fill_fp_pos(&add1);
                }
                else
                {
                    rand_fill_fp_neg(&add1);
                }
                rand_fill_fp_divisor_pos_neg(&add2, i);
                fp_div(&add1, &add2, &res);
                fprintf(file, "/\n");
                fp_fprintf(file, &add1);
                fp_fprintf(file, &add2);
                fp_fprintf(file, &res);
            }
            /*
            Former Not working division example
            fp_t add1;
            fp_t add2;
            fp_t res;

            fp_zero(&add1);
            fp_zero(&add2);
            fp_zero(&res);
            
            add1[0] = 0x7891315c7ce72dad;
            add1[1] = 0x1ec19fd501ccecab;
            add1[2] = 0x70aaa34664304d85;
            add1[3] = 0x5ebb87a36ddb716d;
            add1[4] = 0x689c387722f72bac;
            add1[5] = 0x2a20de9952c2d153;
            add1[6] = 0x6ef489235c3d6841;
            add1[7] = 0x1ff353bd963bc;

            
            add2[0] = 0x1395f1f561a6c9f6;
            add2[1] = 0x2da872d643e2480d;
            add2[2] = 0x124f05a12281efc9;
            add2[3] = 0x727051661f275425;
            add2[4] = 0x470d908225186830;
            add2[5] = 0x342f9bee5dd7e18f;
            add2[6] = 0x148209ba5b5ecbf3;

            

            fprintf(file, "/\n");
            fp_div(&add1, &add2, &res);
            fp_fprintf(file, &add1);
            fp_fprintf(file, &add2);
            fp_fprintf(file, &res);*/
            /*
            //test for fp_mod
            printf("Generating Testdata for fp_mod\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp_t add1;
                fp_t mod;
                fp_t res;

                fp_zero(&add1);
                fp_zero(&mod);
                fp_zero(&res);
                switch(i % 4)
                {
                    case 0:
                        rand_fill_fp_pos(&add1);
                        rand_fill_fp_pos(&mod);
                        break;
                    case 1:
                        rand_fill_fp_neg(&add1);
                        rand_fill_fp_pos(&mod);
                        break;
                    case 2:
                        rand_fill_fp_pos(&add1);
                        rand_fill_fp_divisor(&mod);
                        break;
                    case 3:
                        rand_fill_fp_neg(&add1);
                        rand_fill_fp_divisor(&mod);
                        break;
                    default:
                        assert(0 && "fp_mod invalid configuration");
                }
                fp_mod(&add1, &mod, &res);
                fprintf(file, "%%\n");
                fp_fprintf(file, &add1);
                fp_fprintf(file, &mod);
                fp_fprintf(file, &res);
            }*/

            /*
            //test for fp_mod_2k
            printf("Generating Testdata for fp_mod_2k\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp_t add1;
                fp_t mod;
                fp_t res;
                fp_zero(&add1);
                fp_zero(&res);
                fp_zero(&mod);
                if(i % 2)
                {
                    rand_fill_fp_pos(&add1); 
                }
                else
                {
                    rand_fill_fp_neg(&add1);
                }
                rand_2k_num(&mod);
                fp_mod_2k(&add1, &mod, &res);
                fprintf(file, "%%2k\n");
                fp_fprintf(file, &add1);
                fp_fprintf(file, &mod);
                fp_fprintf(file, &res);
            }

            // test for fp_addm
            printf("Generating Testdata for fp_addm\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp_t add1;
                fp_t add2;
                fp_t mod;
                fp_t res;
                fp_zero(&add1);
                fp_zero(&add2);
                fp_zero(&mod);
                fp_zero(&res);
                rand_fill_fp_2_input(&add1, &add2, i);
                rand_fill_fp_divisor(&mod);
                fp_addm(&add1, &add2, &mod, &res);
                fprintf(file, "+%%\n");
                fp_fprintf(file, &add1);
                fp_fprintf(file, &add2);
                fp_fprintf(file, &mod);
                fp_fprintf(file, &res);
            }

            // test for fp_subm
            printf("Generating Testdata for fp_subm\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp_t add1;
                fp_t add2;
                fp_t mod;
                fp_t res;
                fp_zero(&add1);
                fp_zero(&add2);
                fp_zero(&mod);
                fp_zero(&res);
                rand_fill_fp_2_input(&add1, &add2, i);
                rand_fill_fp_divisor(&mod);
                fp_subm(&add1, &add2, &mod, &res);
                fprintf(file, "-%%\n");
                fp_fprintf(file, &add1);
                fp_fprintf(file, &add2);
                fp_fprintf(file, &mod);
                fp_fprintf(file, &res);
            }

            //test for f2p_mul
            printf("Generating Testdata for f2p_mul\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                f2p_t add1;
                f2p_t add2;
                f2p_t res;
                f2p_zero(&add1);
                f2p_zero(&add2);
                f2p_zero(&res);
                rand_fill_fp_for_f2p_2_input(&add1, &add2, i);
                f2p_mul(&add1, &add2, &res);
                fprintf(file, "**\n");
                f2p_fprintf(file, &add1);
                f2p_fprintf(file, &add2);
                f2p_fprintf(file, &res);
            }


            //test for f2p_div
            printf("Generating Testdata for f2p_div\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                f2p_t add1;
                f2p_t add2;
                f2p_t res;

                f2p_zero(&add1);
                f2p_zero(&add2);
                f2p_zero(&res);
                if((i % 4)>= 2)
                {
                    rand_fill_f2p_pos(&add1);
                }
                else
                {
                    rand_fill_f2p_neg(&add1);
                }
                rand_fill_f2p_divisor_pos_neg(&add2, i);
                f2p_div(&add1, &add2, &res);
                fprintf(file, "//\n");
                f2p_fprintf(file, &add1);
                f2p_fprintf(file, &add2);
                f2p_fprintf(file, &res);
            }




            //test for f2p_mod
            printf("Generating Testdata for f2p_mod\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                f2p_t add1;
                fp_t mod;
                fp_t res;

                f2p_zero(&add1);
                fp_zero(&mod);
                fp_zero(&res);
                switch(i % 4)
                {
                    case 0:
                        rand_fill_f2p_pos(&add1);
                        rand_fill_fp_pos(&mod);
                        break;
                    case 1:
                        rand_fill_f2p_neg(&add1);
                        rand_fill_fp_pos(&mod);
                        break;
                    case 2:
                        rand_fill_f2p_pos(&add1);
                        rand_fill_fp_divisor(&mod);
                        break;
                    case 3:
                        rand_fill_f2p_neg(&add1);
                        rand_fill_fp_divisor(&mod);
                        break;
                    default:
                        assert(0 && "f2p_mod invalid configuration");
                }
                f2p_mod(&add1, &mod, &res);
                fprintf(file, "%%%%\n");
                f2p_fprintf(file, &add1);
                fp_fprintf(file, &mod);
                fp_fprintf(file, &res);
            }




            //test for f2p_div_2k
            printf("Generating Testdata for f2p_div_2k\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                f2p_t add1;
                f2p_t add2;
                f2p_t res;

                f2p_zero(&add1);
                f2p_zero(&add2);
                f2p_zero(&res);
                if(i % 2)
                {
                    rand_fill_f2p_pos(&add1); 
                }
                else
                {
                    rand_fill_f2p_neg(&add1);
                }
                rand_2k_num_f2p(&add2);
                f2p_div_2k(&add1, &add2, &res);
                fprintf(file, "//2k\n");
                f2p_fprintf(file, &add1);
                f2p_fprintf(file, &add2);
                f2p_fprintf(file, &res);
            }


            //test for f2p_mod_2k
            printf("Generating Testdata for f2p_mod_2k\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                f2p_t add1;
                fp_t mod;
                fp_t res;

                f2p_zero(&add1);
                fp_zero(&mod);
                fp_zero(&res);
                if(i % 2)
                    rand_fill_f2p_pos(&add1);
                else
                    rand_fill_f2p_neg(&add1);
                rand_2k_num(&mod);
                f2p_mod_2k(&add1, &mod, &res);
                fprintf(file, "%%%%2k\n");
                f2p_fprintf(file, &add1);
                fp_fprintf(file, &mod);
                fp_fprintf(file, &res);
            }

            // test fp_get_len
            printf("Generating Testdata for fp_get_len\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp_t add1;
                int64_t res;
                fp_zero(&add1);
                rand_fill_fp_pos(&add1);
                res = fp_get_len(&add1);
                fprintf(file, "fp_get_len\n");
                fp_fprintf(file, &add1);
                fprintf(file, "%ld\n", res);
            }*/
        }

        if(FP2TESTS)
        {
            // test fp2_add
            printf("Generating Testdata for fp2_add\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp2_t add1;
                fp2_t add2;
                fp2_t res;

                fp2_zero(&add1);
                fp2_zero(&add2);
                fp2_zero(&res);
                rand_fill_fp2_2_input(&add1, &add2, i);
                fp2_add(&add1, &add2, &res);
                fprintf(file, "+fp2\n");
                fp2_fprintf(file, &add1);
                fp2_fprintf(file, &add2);
                fp2_fprintf(file, &res);
            }


            // test fp2_sub
            printf("Generating Testdata for fp2_sub\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp2_t add1;
                fp2_t add2;
                fp2_t res;

                fp2_zero(&add1);
                fp2_zero(&add2);
                fp2_zero(&res);
                rand_fill_fp2_2_input(&add1, &add2, i);
                fp2_sub(&add1, &add2, &res);
                fprintf(file, "-fp2\n");
                fp2_fprintf(file, &add1);
                fp2_fprintf(file, &add2);
                fp2_fprintf(file, &res);
            }


            // test fp2_mul
            printf("Generating Testdata for fp2_mul\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp2_t add1;
                fp2_t add2;
                f2p2_t res;
                fp2_zero(&add1);
                fp2_zero(&add2);
                f2p2_zero(&res);
                rand_fill_fp2_2_input(&add1, &add2, i);
                fp2_mul(&add1, &add2, &res);
                fprintf(file, "*fp2\n");
                fp2_fprintf(file, &add1);
                fp2_fprintf(file, &add2);
                f2p2_fprintf(file, &res);
            }
        }

        if(REDCTESTS)
        {
            printf("Generating Testdata for montgomory reduction\n");
            
            
            printf("with fp\n");

            


            printf("Generating Testdata for REDCL\n");


            printf("REDCL sanitiy\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp_t add1;
                fp_t add2;
                fp_t res;
                f2p_t mul1;
                fp_zero(&add1);
                fp_zero(&res);
                //rand_fill_fp_2_input(&add1, &add2, i);
                rand_fp_smaller_mod(&add1, &mont_mod);
                //add1[0] = 34265234;
                fp_zero(&add2);
                f2p_zero(&mul1);
                fp_mul(&add1, &r, &mul1);
                REDCL(&mul1, &mont_mod, &add2);
                //fp_to_f2p(&add2, &mul1);
                //REDCL(&mul1, &mont_mod, &res);

                fprintf(file, "REDECL-Sanity\n");
                fp_fprintf(file, &add1);
                fp_fprintf(file, &mont_mod);
                fp_fprintf(file, &add2);
            }




            printf("addition\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp_t add1;
                fp_t add2;
                fp_t res;
                f2p_t mul1;
                f2p_t mul2;
                fp_zero(&add1);
                fp_zero(&add2);
                fp_zero(&res);
                //rand_fill_fp_2_input(&add1, &add2, i);
                rand_fp_smaller_mod(&add1, &mont_mod);
                rand_fp_smaller_mod(&add2, &mont_mod);
                fp_t add1_redc;
                fp_t add2_redc;
                fp_t res_redc;
                f2p_t res_redc_l;
                fp_zero(&add1_redc);
                fp_zero(&add2_redc);
                fp_zero(&res_redc);
                f2p_zero(&res_redc_l);

                fp_mul(&add1, &r2, &mul1);
                fp_mul(&add2, &r2, &mul2);
                REDCL(&mul1, &mont_mod, &add1_redc);
                REDCL(&mul2, &mont_mod, &add2_redc);

                fp_add(&add1_redc, &add2_redc, &res_redc);
                fp_to_f2p(&res_redc, &res_redc_l);

                REDCL(&res_redc_l, &mont_mod, &res);


                fprintf(file, "REDECL+\n");
                fp_fprintf(file, &add1);
                fp_fprintf(file, &add2);
                fp_fprintf(file, &mont_mod);
                fp_fprintf(file, &res);
            }

            printf("subtraction\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp_t add1;
                fp_t add2;
                fp_t res;
                f2p_t mul1;
                f2p_t mul2;
                fp_zero(&add1);
                fp_zero(&add2);
                fp_zero(&res);
                //rand_fill_fp_2_input(&add1, &add2, i);
                rand_fp_smaller_mod(&add1, &mont_mod);
                rand_fp_smaller_mod(&add2, &mont_mod);
                fp_t add1_redc;
                fp_t add2_redc;
                fp_t res_redc;
                f2p_t res_redc_l;
                fp_zero(&add1_redc);
                fp_zero(&add2_redc);
                fp_zero(&res_redc);
                f2p_zero(&res_redc_l);

                fp_mul(&add1, &r2, &mul1);
                fp_mul(&add2, &r2, &mul2);
                REDCL(&mul1, &mont_mod, &add1_redc);
                REDCL(&mul2, &mont_mod, &add2_redc);

                fp_sub(&add1_redc, &add2_redc, &res_redc);
                fp_to_f2p(&res_redc, &res_redc_l);
                //f2p_t c1;
                //f2p_copy(&res_redc_l, &c1);
                //f2p_mod(&c1, &mont_mod, &res_redc);
                //fp_to_f2p(&res_redc, &res_redc_l);

                REDCL(&res_redc_l, &mont_mod, &res);


                fprintf(file, "REDECL-\n");
                fp_fprintf(file, &add1);
                fp_fprintf(file, &add2);
                fp_fprintf(file, &mont_mod);
                fp_fprintf(file, &res);
            }


            printf("Generating Testdata for MODMUL\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp_t add1;
                fp_t add2;
                fp_t res;

                fp_zero(&add1);
                fp_zero(&add2);
                fp_zero(&res);
                rand_fp_smaller_mod(&add1, &mont_mod);
                rand_fp_smaller_mod(&add2, &mont_mod);
                MODMUL(&add1, &add2, &mont_mod, &res);
                fprintf(file, "MODMUL\n");
                fp_fprintf(file, &add1);
                fp_fprintf(file, &add2);
                fp_fprintf(file, &mont_mod);
                fp_fprintf(file, &res);
            }

            printf("Generating Testdata for REDCL2\n");
            printf("REDCL2 sanitiy\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp2_t add1;
                fp2_t add2;
                fp2_t rrr;
                fp2_t res;
                f2p2_t mul1;
                fp2_zero(&rrr);
                fp2_zero(&add1);
                fp2_zero(&res);
                //rand_fill_fp_2_input(&add1, &add2, i);
                fp_copy(&r, &(rrr.real));
                //fp_copy(&r, &(rrr.img));
                rand_fp2_smaller_mod(&add1, &mont_mod);
                //add1[0] = 34265234;
                fp2_zero(&add2);
                f2p2_zero(&mul1);
                fp2_mul(&add1, &rrr, &mul1);
                REDCL2(&mul1, &mont_mod, &add2);
                //fp_to_f2p(&add2, &mul1);
                //REDCL(&mul1, &mont_mod, &res);

                fprintf(file, "REDECL2-Sanity\n");
                fp2_fprintf(file, &add1);
                fp_fprintf(file, &mont_mod);
                fp2_fprintf(file, &add2);
            }



            printf("addition\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp2_t add1;
                fp2_t add2;
                fp2_t res;
                f2p2_t mul1;
                f2p2_t mul2;
                fp2_zero(&add1);
                fp2_zero(&add2);
                fp2_zero(&res);
                //rand_fill_fp_2_input(&add1, &add2, i);
                rand_fp2_smaller_mod(&add1, &mont_mod);
                rand_fp2_smaller_mod(&add2, &mont_mod);
                fp2_t add1_redc;
                fp2_t add2_redc;
                fp2_t res_redc;
                f2p2_t res_redc_l;
                fp2_zero(&add1_redc);
                fp2_zero(&add2_redc);
                fp2_zero(&res_redc);
                f2p2_zero(&res_redc_l);

            
                fp2_mul(&add1, &r22, &mul1);
                fp2_mul(&add2, &r22, &mul2);
                REDCL2(&mul1, &mont_mod, &add1_redc);
                REDCL2(&mul2, &mont_mod, &add2_redc);

                fp2_add(&add1_redc, &add2_redc, &res_redc);
                fp2_to_f2p2(&res_redc, &res_redc_l);

                REDCL2(&res_redc_l, &mont_mod, &res);


                fprintf(file, "REDECL2+\n");
                fp2_fprintf(file, &add1);
                fp2_fprintf(file, &add2);
                fp_fprintf(file, &mont_mod);
                fp2_fprintf(file, &res);
            }


            printf("mult_inv\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp2_t add1;
                fp2_t res;
                f2p2_t mul1;
                fp2_zero(&add1);
                fp2_zero(&res);
                rand_fp2_smaller_mod(&add1, &mont_mod);
                fp2_t add1_redc;
                fp2_t res_redc;
                fp2_zero(&add1_redc);
                fp2_zero(&res_redc);

            
                fp2_mul(&add1, &r22, &mul1);;
                REDCL2(&mul1, &mont_mod, &add1_redc);

                fp2_mult_inv_mont(&add1_redc, &mont_mod, &res_redc);

                REDC2(&res_redc, &mont_mod, &res);


                fprintf(file, "Mont_mult_inv\n");
                fp2_fprintf(file, &add1);
                fp_fprintf(file, &mont_mod);
                fp2_fprintf(file, &res);
            }


            printf("multiplication\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp2_t add1;
                fp2_t add2;
                fp2_t res;
                f2p2_t mul1;
                f2p2_t mul2;
                fp2_zero(&add1);
                fp2_zero(&add2);
                fp2_zero(&res);
                //rand_fill_fp_2_input(&add1, &add2, i);
                rand_fp2_smaller_mod(&add1, &mont_mod);
                rand_fp2_smaller_mod(&add2, &mont_mod);
                fp2_t add1_redc;
                fp2_t add2_redc;
                fp2_t res_redc;
                f2p2_t res_redc_l;
                fp2_zero(&add1_redc);
                fp2_zero(&add2_redc);
                fp2_zero(&res_redc);
                f2p2_zero(&res_redc_l);

            
                fp2_mul(&add1, &r22, &mul1);
                fp2_mul(&add2, &r22, &mul2);
                REDCL2(&mul1, &mont_mod, &add1_redc);
                REDCL2(&mul2, &mont_mod, &add2_redc);

                fp2_mul(&add1_redc, &add2_redc, &res_redc_l);
                //fp2_to_f2p2(&res_redc, &res_redc_l);

                REDCL2(&res_redc_l, &mont_mod, &res_redc);
                f2p2_zero(&res_redc_l);
                fp2_to_f2p2(&res_redc, &res_redc_l);
                REDCL2(&res_redc_l, &mont_mod, &res);


                fprintf(file, "REDECL2*\n");
                fp2_fprintf(file, &add1);
                fp2_fprintf(file, &add2);
                fp_fprintf(file, &mont_mod);
                fp2_fprintf(file, &res);
            }


            printf("fp2_mul_mont\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                fp2_t add1;
                fp2_t add2;
                fp2_t res;
                f2p2_t mul1;
                f2p2_t mul2;
                fp2_zero(&add1);
                fp2_zero(&add2);
                fp2_zero(&res);
                //rand_fill_fp_2_input(&add1, &add2, i);
                rand_fp2_smaller_mod(&add1, &mont_mod);
                rand_fp2_smaller_mod(&add2, &mont_mod);
                fp2_t add1_redc;
                fp2_t add2_redc;
                fp2_t res_redc;
                f2p2_t res_redc_l;
                fp2_zero(&add1_redc);
                fp2_zero(&add2_redc);
                fp2_zero(&res_redc);
                f2p2_zero(&res_redc_l);

            
                fp2_mul(&add1, &r22, &mul1);
                fp2_mul(&add2, &r22, &mul2);
                REDCL2(&mul1, &mont_mod, &add1_redc);
                REDCL2(&mul2, &mont_mod, &add2_redc);

                fp2_mul_mont(&add1_redc, &add2_redc, &mont_mod, &res_redc);
                //fp2_to_f2p2(&res_redc, &res_redc_l);
                fp2_to_f2p2(&res_redc, &res_redc_l);
                REDCL2(&res_redc_l, &mont_mod, &res);



                fprintf(file, "fp2_mul_mont\n");
                fp2_fprintf(file, &add1);
                fp2_fprintf(file, &add2);
                fp_fprintf(file, &mont_mod);
                fp2_fprintf(file, &res);
            }
        }

        mont_curve_t A;
        mont_curve_t A_redc;
        mont_curve_t A24p;
        mont_curve_t A24p_redc;
        mont_curve_t A24m;
        mont_curve_t A24m_redc;
        
        fp2_t C2;
        fp2_zero(&C2);
        proj_pt_zero((proj_point_t*)&A_redc);
        proj_pt_zero((proj_point_t*)&A);
        proj_pt_zero((proj_point_t*)&A24p);
        proj_pt_zero((proj_point_t*)&A24p_redc);
        proj_pt_zero((proj_point_t*)&A24m);
        proj_pt_zero((proj_point_t*)&A24m_redc);

        A.A.real[0] = 6;
        A.C.real[0] = 1;
        fp2_copy(&six_mont, &A_redc.A);
        fp2_copy(&one_mont, &A_redc.C);
        fp2_add(&A.C, &A.C, &C2);
        fp2_add(&A.A, &C2, &A24p.A);
        fp2_add(&C2, &C2, &A24p.C);

        fp2_copy(&A24p.A, &A24m.A);
        fp2_sub(&A.A, &C2, &A24m.C);

        f2p2_t temporary;
        f2p2_zero(&temporary);

        fp2_mul(&(A24p.A), &r22, &temporary);
        REDCL2(&temporary, &mont_mod, &(A24p_redc.A));
        fp2_mul(&(A24p.C), &r22, &temporary);
        REDCL2(&temporary, &mont_mod, &(A24p_redc.C));

        fp2_copy(&A24p_redc.A, &A24m_redc.A);
        fp2_mul(&(A24m.C), &r22, &temporary);
        REDCL2(&temporary, &mont_mod, &(A24m_redc.C));
        
        if(POINTARTESTS)
        {
            printf("Generating Testdata for Point Arithmetic\n");
            printf("Generating Testdata for xDBL\n");
            


            for(int i = 0; i < NUM_TESTS; i ++)
            {
                proj_point_t add1;
                f2p2_t temp;
                proj_point_t add1_redc;
                f2p2_zero(&temp);
                proj_pt_zero(&add1);
                proj_pt_zero(&add1_redc);
                rand_fill_proj_pt(&add1, &mont_mod);
                fp2_mul(&(add1.X), &r22, &temp);
                REDCL2(&temp, &mont_mod, &(add1_redc.X));
                fp2_mul(&(add1.Z), &r22, &temp);
                REDCL2(&temp, &mont_mod, &(add1_redc.Z));
                proj_point_t res;
                proj_pt_zero(&res);
                xDBL(&add1_redc, &A24p_redc, &mont_mod, &res);
                fprintf(file, "xDBL\n");
                proj_pt_fprintf(file, &add1);
                proj_pt_fprintf(file, (proj_point_t*)&A24p);
                fp_fprintf(file, &mont_mod);
                fp2_t X_red;
                fp2_t Y_red;
                fp2_zero(&X_red);
                fp2_zero(&Y_red);
                f2p2_t t;
                fp2_to_f2p2(&(res.X), &t);
                REDCL2(&t, &mont_mod, &X_red);
                fp2_to_f2p2(&(res.Z), &t);
                REDCL2(&t, &mont_mod, &Y_red);
                fp2_fprintf(file, &X_red);
                fp2_fprintf(file, &Y_red);
            }

            if(XDBLE)
            {
                printf("Generating Testdata for XDBLe\n");
                for(int i = 0; i < NUM_TESTS_COSTLY; i ++)
                {
                    printf("xDBLe %2dth test of %d\n", i + 1, NUM_TESTS_COSTLY);
                    proj_point_t add1;
                    f2p2_t temp;
                    proj_point_t add1_redc;
                    int e = rand() % 217;
                    f2p2_zero(&temp);
                    proj_pt_zero(&add1);
                    proj_pt_zero(&add1_redc);
                    rand_fill_proj_pt(&add1, &mont_mod);
                    fp2_mul(&(add1.X), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(add1_redc.X));
                    fp2_mul(&(add1.Z), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(add1_redc.Z));
                    proj_point_t res;
                    proj_pt_zero(&res);
                    xDBLe(&add1_redc, &A24p_redc, e, &mont_mod, &res);
                    fprintf(file, "xDBLe\n");
                    proj_pt_fprintf(file, &add1);
                    proj_pt_fprintf(file, (proj_point_t*)&A24p);
                    fprintf(file, "%d\n", e);
                    fp_fprintf(file, &mont_mod);
                    fp2_t X_red;
                    fp2_t Y_red;
                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    f2p2_t t;
                    fp2_to_f2p2(&(res.X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    fp2_to_f2p2(&(res.Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);
                }
            }
            
            printf("Generating Testdata for xDBLADD\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                proj_point_t add1;
                proj_point_t add2;
                proj_point_t add3;
                f2p2_t temp;
                proj_point_t add1_redc;
                proj_point_t add2_redc;
                proj_point_t add3_redc;
                f2p2_zero(&temp);
                proj_pt_zero(&add1);
                proj_pt_zero(&add1_redc);
                proj_pt_zero(&add2);
                proj_pt_zero(&add2_redc);
                proj_pt_zero(&add3);
                proj_pt_zero(&add3_redc);
                rand_fill_proj_pt(&add1, &mont_mod);
                rand_fill_proj_pt(&add2, &mont_mod);
                rand_fill_proj_pt(&add3, &mont_mod);
                fp2_mul(&(add1.X), &r22, &temp);
                REDCL2(&temp, &mont_mod, &(add1_redc.X));
                fp2_mul(&(add1.Z), &r22, &temp);
                REDCL2(&temp, &mont_mod, &(add1_redc.Z));

                fp2_mul(&(add2.X), &r22, &temp);
                REDCL2(&temp, &mont_mod, &(add2_redc.X));
                fp2_mul(&(add2.Z), &r22, &temp);
                REDCL2(&temp, &mont_mod, &(add2_redc.Z));

                fp2_mul(&(add3.X), &r22, &temp);
                REDCL2(&temp, &mont_mod, &(add3_redc.X));
                fp2_mul(&(add3.Z), &r22, &temp);
                REDCL2(&temp, &mont_mod, &(add3_redc.Z));


                proj_point_t res1;
                proj_point_t res2;
                proj_pt_zero(&res1);
                proj_pt_zero(&res2);
                xDBLADD(&add1_redc, &add2_redc, &add3_redc, &mont_mod, &A24p_redc, &res1, &res2);

                fprintf(file, "xDBLADD\n");
                proj_pt_fprintf(file, &add1);
                proj_pt_fprintf(file, &add2);
                proj_pt_fprintf(file, &add3);
                proj_pt_fprintf(file, (proj_point_t*)&A24p);
                fp_fprintf(file, &mont_mod);


                fp2_t X_red;
                fp2_t Y_red;
                fp2_zero(&X_red);
                fp2_zero(&Y_red);
                f2p2_t t;
                fp2_to_f2p2(&(res1.X), &t);
                REDCL2(&t, &mont_mod, &X_red);
                fp2_to_f2p2(&(res1.Z), &t);
                REDCL2(&t, &mont_mod, &Y_red);
                fp2_fprintf(file, &X_red);
                fp2_fprintf(file, &Y_red);

                f2p2_zero(&t);
                fp2_to_f2p2(&(res2.X), &t);
                REDCL2(&t, &mont_mod, &X_red);
                fp2_to_f2p2(&(res2.Z), &t);
                REDCL2(&t, &mont_mod, &Y_red);
                fp2_fprintf(file, &X_red);
                fp2_fprintf(file, &Y_red);
            }

            printf("Generating Testdata for xTPL\n");
            for(int i = 0; i < NUM_TESTS; i ++)
            {
                proj_point_t add1;
                f2p2_t temp;
                proj_point_t add1_redc;
                f2p2_zero(&temp);
                proj_pt_zero(&add1);
                proj_pt_zero(&add1_redc);
                rand_fill_proj_pt(&add1, &mont_mod);
                fp2_mul(&(add1.X), &r22, &temp);
                REDCL2(&temp, &mont_mod, &(add1_redc.X));
                fp2_mul(&(add1.Z), &r22, &temp);
                REDCL2(&temp, &mont_mod, &(add1_redc.Z));
                proj_point_t res;
                proj_pt_zero(&res);
                xTPL(&add1_redc, &A24m_redc, &mont_mod, &res);
                fprintf(file, "xTPL\n");
                proj_pt_fprintf(file, &add1);
                proj_pt_fprintf(file, (proj_point_t*)&A24m);
                fp_fprintf(file, &mont_mod);
                fp2_t X_red;
                fp2_t Y_red;
                fp2_zero(&X_red);
                fp2_zero(&Y_red);
                f2p2_t t;
                fp2_to_f2p2(&(res.X), &t);
                REDCL2(&t, &mont_mod, &X_red);
                fp2_to_f2p2(&(res.Z), &t);
                REDCL2(&t, &mont_mod, &Y_red);
                fp2_fprintf(file, &X_red);
                fp2_fprintf(file, &Y_red);
            }

            if(XTPLE)
            {
                printf("Generating Testdata for xTPLe\n");
                for(int i = 0; i < NUM_TESTS_COSTLY; i ++)
                {
                    printf("xTPLe %2dth test of %d\n", i + 1, NUM_TESTS_COSTLY);
                    proj_point_t add1;
                    f2p2_t temp;
                    proj_point_t add1_redc;
                    int e = rand() % 217;
                    f2p2_zero(&temp);
                    proj_pt_zero(&add1);
                    proj_pt_zero(&add1_redc);
                    rand_fill_proj_pt(&add1, &mont_mod);
                    fp2_mul(&(add1.X), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(add1_redc.X));
                    fp2_mul(&(add1.Z), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(add1_redc.Z));
                    proj_point_t res;
                    proj_pt_zero(&res);
                    xTPLe(&add1_redc, &A24p_redc, e, &mont_mod, &res);
                    fprintf(file, "xTPLe\n");
                    proj_pt_fprintf(file, &add1);
                    proj_pt_fprintf(file, (proj_point_t*)&A24p);
                    fprintf(file, "%d\n", e);
                    fp_fprintf(file, &mont_mod);
                    fp2_t X_red;
                    fp2_t Y_red;
                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    f2p2_t t;
                    fp2_to_f2p2(&(res.X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    fp2_to_f2p2(&(res.Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);
                }
            }

            if(LADDER3P)
            {
                printf("Generating Testdata for ladder3p\n");
                for(int i = 0; i < NUM_LADDER_TESTS; i ++)
                {
                    char* post = i == 0 ? "st" : i == 1 ? "nd" : "th";
                    printf("ladder3p %2d%s test of %d\n", i + 1, post, NUM_LADDER_TESTS);
                    proj_point_t add1;
                    proj_point_t add2;
                    proj_point_t add3;
                    f2p2_t temp;
                    proj_point_t add1_redc;
                    proj_point_t add2_redc;
                    proj_point_t add3_redc;
                    f2p2_zero(&temp);
                    proj_pt_zero(&add1);
                    proj_pt_zero(&add1_redc);
                    proj_pt_zero(&add2);
                    proj_pt_zero(&add2_redc);
                    proj_pt_zero(&add3);
                    proj_pt_zero(&add3_redc);
                    rand_fill_proj_pt(&add1, &mont_mod);
                    rand_fill_proj_pt(&add2, &mont_mod);
                    rand_fill_proj_pt(&add3, &mont_mod);
                    fp2_mul(&(add1.X), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(add1_redc.X));
                    //fp2_mul(&(add1.Z), &r22, &temp);
                    //REDCL2(&temp, &mont_mod, &(add1_redc.Z));

                    fp2_mul(&(add2.X), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(add2_redc.X));
                    //fp2_mul(&(add2.Z), &r22, &temp);
                    //REDCL2(&temp, &mont_mod, &(add2_redc.Z));

                    fp2_mul(&(add3.X), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(add3_redc.X));
                    //fp2_mul(&(add3.Z), &r22, &temp);
                    //REDCL2(&temp, &mont_mod, &(add3_redc.Z));

                    fp_t secret_key;
                    rand_fill_fp_pos(&secret_key);
                    //fp_zero(&secret_key);
                    //secret_key[0] = i+1;
                    proj_point_t res1;
                    proj_pt_zero(&res1);
                    

                    //xDBLADD(&add1_redc, &add2_redc, &add3_redc, &mont_mod, &A24p_redc, &res1, &res2);

                    Ladder3p(&add1_redc.X, &add2_redc.X, &add3_redc.X, &secret_key, &res1, &A24p_redc, &mont_mod);

                    fprintf(file, "ladder3p\n");
                    fp2_fprintf(file, &add1.X);
                    fp2_fprintf(file, &add2.X);
                    fp2_fprintf(file, &add3.X);
                    fp_fprintf(file, &secret_key);
                    proj_pt_fprintf(file, (proj_point_t*)&A24p);
                    fp_fprintf(file, &mont_mod);


                    fp2_t X_red;
                    fp2_t Y_red;
                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    f2p2_t t;
                    fp2_to_f2p2(&(res1.X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    fp2_to_f2p2(&(res1.Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);
                }
            }
        }

        if(MONTCURVETESTS)
        {
            printf("Generating Testdata for for Montgomery curve functions\n");

            if(MONTCURVETESTSREST)
            {
                printf("Generating Testdata for for jinvariant\n");
                for(int i = 0; i < NUM_TESTS; i ++)
                {
                    fp2_t AA;
                    fp2_t CC;
                    fp2_t j;
                    f2p2_t mul1;
                    f2p2_t mul2;
                    fp2_zero(&AA);
                    fp2_zero(&CC);
                    fp2_zero(&j);

                    rand_fp2_smaller_mod(&AA, &mont_mod);
                    rand_fp2_smaller_mod(&CC, &mont_mod);
                    fp2_t A_redc;
                    fp2_t C_redc;
                    fp2_t j_redc;
                    f2p2_t j_redc_l;
                    fp2_zero(&A_redc);
                    fp2_zero(&C_redc);
                    fp2_zero(&j_redc);
                    f2p2_zero(&j_redc_l);

                
                    fp2_mul(&AA, &r22, &mul1);
                    fp2_mul(&CC, &r22, &mul2);
                    REDCL2(&mul1, &mont_mod, &A_redc);
                    REDCL2(&mul2, &mont_mod, &C_redc);

                    jinvariant(&A_redc, &C_redc, &mont_mod, &j_redc);

                    fp2_to_f2p2(&j_redc, &j_redc_l);
                    fp2_zero(&j);
                    REDC2(&j_redc, &mont_mod, &j);



                    fprintf(file, "jinvariant\n");
                    fp2_fprintf(file, &AA);
                    fp2_fprintf(file, &CC);
                    fp_fprintf(file, &mont_mod);
                    fp2_fprintf(file, &j);
                }

                printf("Generating Testdata for 2_iso_curve\n");
                for(int i = 0; i < NUM_TESTS; i ++)
                {
                    proj_point_t add1;
                    f2p2_t temp;
                    proj_point_t add1_redc;
                    f2p2_zero(&temp);
                    proj_pt_zero(&add1);
                    proj_pt_zero(&add1_redc);
                    rand_fill_proj_pt(&add1, &mont_mod);
                    fp2_mul(&(add1.X), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(add1_redc.X));
                    fp2_mul(&(add1.Z), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(add1_redc.Z));
                    proj_point_t res;
                    proj_pt_zero(&res);
                    iso_2_curve(&add1_redc, &mont_mod, (mont_curve_t*)&res);
                    fprintf(file, "iso_2_curve\n");
                    proj_pt_fprintf(file, &add1);
                    fp_fprintf(file, &mont_mod);
                    fp2_t X_red;
                    fp2_t Y_red;
                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    f2p2_t t;
                    fp2_to_f2p2(&(res.X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    fp2_to_f2p2(&(res.Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);
                    //fp2_fprintf(file, &res.X);
                    //fp2_fprintf(file, &res.Z);
                }

                printf("Generating Testdata for 2_iso_eval\n");
                for(int i = 0; i < NUM_TESTS; i ++)
                {
                    proj_point_t add1;
                    f2p2_t temp;
                    proj_point_t add1_redc;
                    proj_point_t add2;
                    proj_point_t add2_redc;
                    f2p2_zero(&temp);
                    proj_pt_zero(&add1);
                    proj_pt_zero(&add1_redc);
                    proj_pt_zero(&add2);
                    proj_pt_zero(&add2_redc);
                    rand_fill_proj_pt(&add1, &mont_mod);
                    rand_fill_proj_pt(&add2, &mont_mod);
                    fp2_mul(&(add1.X), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(add1_redc.X));
                    fp2_mul(&(add1.Z), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(add1_redc.Z));

                    fp2_mul(&(add2.X), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(add2_redc.X));
                    fp2_mul(&(add2.Z), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(add2_redc.Z));
                    proj_point_t res;
                    proj_pt_zero(&res);
                    iso_2_eval(&add1_redc, &add2_redc, &mont_mod, &res);
                    fprintf(file, "2_iso_eval\n");
                    proj_pt_fprintf(file, &add1);
                    proj_pt_fprintf(file, &add2);
                    fp_fprintf(file, &mont_mod);
                    fp2_t X_red;
                    fp2_t Y_red;
                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    f2p2_t t;
                    f2p2_zero(&t);
                    fp2_to_f2p2(&(res.X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    f2p2_zero(&t);
                    fp2_to_f2p2(&(res.Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);
                }

                printf("Generating Testdata for 4_iso_curve\n");
                for(int i = 0; i < NUM_TESTS; i ++)
                {
                    f2p2_t temp;
                    f2p2_zero(&temp);

                    proj_point_t P4;
                    proj_point_t P4_redc;
                    proj_pt_zero(&P4);
                    proj_pt_zero(&P4_redc);
                    rand_fill_proj_pt(&P4, &mont_mod);
                    fp2_mul(&(P4.X), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(P4_redc.X));
                    fp2_mul(&(P4.Z), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(P4_redc.Z));

                    fp2_t K1;
                    fp2_t K1_redc;

                    fp2_t K2;
                    fp2_t K2_redc;

                    fp2_t K3;
                    fp2_t K3_redc;

                    fp2_zero(&K1);
                    fp2_zero(&K1_redc);
                    fp2_zero(&K2);
                    fp2_zero(&K2_redc);
                    fp2_zero(&K3);
                    fp2_zero(&K3_redc);


                    proj_point_t res;
                    proj_pt_zero(&res);
                    iso_4_curve(&P4_redc, &mont_mod, (mont_curve_t*)&res, &K1_redc, &K2_redc, &K3_redc);
                    REDC2(&K1_redc, &mont_mod, &K1);
                    REDC2(&K2_redc, &mont_mod, &K2);
                    REDC2(&K3_redc, &mont_mod, &K3);

                    fprintf(file, "iso_4_curve\n");
                    proj_pt_fprintf(file, &P4);
                    fp2_fprintf(file, &K1);
                    fp2_fprintf(file, &K2);
                    fp2_fprintf(file, &K3);
                    fp_fprintf(file, &mont_mod);
                    fp2_t X_red;
                    fp2_t Y_red;
                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    f2p2_t t;
                    fp2_to_f2p2(&(res.X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    fp2_to_f2p2(&(res.Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);
                }


                printf("Generating Testdata for 4_iso_eval\n");
                for(int i = 0; i < NUM_TESTS; i ++)
                {
                    f2p2_t temp;
                    f2p2_zero(&temp);

                    proj_point_t Q;
                    proj_point_t Q_redc;
                    proj_pt_zero(&Q);
                    proj_pt_zero(&Q_redc);
                    rand_fill_proj_pt(&Q, &mont_mod);
                    fp2_mul(&(Q.X), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(Q_redc.X));
                    fp2_mul(&(Q.Z), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(Q_redc.Z));

                    fp2_t K1;
                    fp2_t K1_redc;
                    fp2_zero(&K1);
                    fp2_zero(&K1_redc);
                    rand_fp2_smaller_mod(&K1, &mont_mod);
                    fp2_mul(&(K1), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(K1_redc));


                    fp2_t K2;
                    fp2_t K2_redc;
                    fp2_zero(&K2);
                    fp2_zero(&K2_redc);
                    rand_fp2_smaller_mod(&K2, &mont_mod);
                    fp2_mul(&(K2), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(K2_redc));


                    fp2_t K3;
                    fp2_t K3_redc;
                    fp2_zero(&K3);
                    fp2_zero(&K3_redc);
                    rand_fp2_smaller_mod(&K3, &mont_mod);
                    fp2_mul(&(K3), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(K3_redc));

                    proj_point_t res;
                    proj_pt_zero(&res);
                    iso_4_eval(&K1_redc, &K2_redc, &K3_redc, &Q_redc, &mont_mod, &res);

                    fprintf(file, "iso_4_eval\n");
                    proj_pt_fprintf(file, &Q);
                    fp2_fprintf(file, &K1);
                    fp2_fprintf(file, &K2);
                    fp2_fprintf(file, &K3);
                    fp_fprintf(file, &mont_mod);
                    fp2_t X_red;
                    fp2_t Y_red;
                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    f2p2_t t;
                    fp2_to_f2p2(&(res.X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    fp2_to_f2p2(&(res.Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);
                }

                printf("Generating Testdata for 3_iso_curve\n");
                for(int i = 0; i < NUM_TESTS; i ++)
                {
                    f2p2_t temp;
                    f2p2_zero(&temp);

                    proj_point_t P3;
                    proj_point_t P3_redc;
                    proj_pt_zero(&P3);
                    proj_pt_zero(&P3_redc);
                    rand_fill_proj_pt(&P3, &mont_mod);
                    fp2_mul(&(P3.X), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(P3_redc.X));
                    fp2_mul(&(P3.Z), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(P3_redc.Z));

                    fp2_t K1;
                    fp2_t K1_redc;

                    fp2_t K2;
                    fp2_t K2_redc;

                    fp2_zero(&K1);
                    fp2_zero(&K1_redc);
                    fp2_zero(&K2);
                    fp2_zero(&K2_redc);


                    proj_point_t res;
                    proj_pt_zero(&res);
                    iso_3_curve(&P3_redc, &mont_mod, (mont_curve_t*)&res, &K1_redc, &K2_redc);
                    REDC2(&K1_redc, &mont_mod, &K1);
                    REDC2(&K2_redc, &mont_mod, &K2);

                    fprintf(file, "iso_3_curve\n");
                    proj_pt_fprintf(file, &P3);
                    fp2_fprintf(file, &K1);
                    fp2_fprintf(file, &K2);
                    fp_fprintf(file, &mont_mod);
                    fp2_t X_red;
                    fp2_t Y_red;
                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    f2p2_t t;
                    fp2_to_f2p2(&(res.X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    fp2_to_f2p2(&(res.Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);
                }

                printf("Generating Testdata for 3_iso_eval\n");
                for(int i = 0; i < NUM_TESTS; i ++)
                {
                    f2p2_t temp;
                    f2p2_zero(&temp);

                    proj_point_t Q;
                    proj_point_t Q_redc;
                    proj_pt_zero(&Q);
                    proj_pt_zero(&Q_redc);
                    rand_fill_proj_pt(&Q, &mont_mod);
                    fp2_mul(&(Q.X), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(Q_redc.X));
                    fp2_mul(&(Q.Z), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(Q_redc.Z));

                    fp2_t K1;
                    fp2_t K1_redc;
                    fp2_zero(&K1);
                    fp2_zero(&K1_redc);
                    rand_fp2_smaller_mod(&K1, &mont_mod);
                    fp2_mul(&(K1), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(K1_redc));


                    fp2_t K2;
                    fp2_t K2_redc;
                    fp2_zero(&K2);
                    fp2_zero(&K2_redc);
                    rand_fp2_smaller_mod(&K2, &mont_mod);
                    fp2_mul(&(K2), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(K2_redc));


                    proj_point_t res;
                    proj_pt_zero(&res);
                    iso_3_eval(&K1_redc, &K2_redc, &Q_redc, &mont_mod, &res);

                    fprintf(file, "iso_3_eval\n");
                    proj_pt_fprintf(file, &Q);
                    fp2_fprintf(file, &K1);
                    fp2_fprintf(file, &K2);
                    fp_fprintf(file, &mont_mod);
                    fp2_t X_red;
                    fp2_t Y_red;
                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    f2p2_t t;
                    fp2_to_f2p2(&(res.X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    fp2_to_f2p2(&(res.Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);
                }

                printf("Generating Testdata for get_A\n");
                for(int i = 0; i < NUM_TESTS; i ++)
                {
                    fp2_t xp;
                    fp2_t xq;
                    fp2_t xpq;
                    fp2_t xp_redc;
                    fp2_t xq_redc;
                    fp2_t xpq_redc;

                    fp2_t A;
                    fp2_t A_non_redc;

                    rand_fill_fp2(&xp);
                    rand_fill_fp2(&xq);
                    rand_fill_fp2(&xpq);

                    fp2_to_mont(&xp, &mont_mod, &xp_redc);
                    fp2_to_mont(&xq, &mont_mod, &xq_redc);
                    fp2_to_mont(&xpq, &mont_mod, &xpq_redc);

                    get_A(&xp_redc, &xq_redc, &xpq_redc, &mont_mod, &A);

                    fp2_from_mont(&A, &mont_mod, &A_non_redc);

                    fprintf(file, "get_A\n");
                    fp2_fprintf(file, &xp);
                    fp2_fprintf(file, &xq);
                    fp2_fprintf(file, &xpq);

                    fp2_fprintf(file, &A_non_redc);
                }

            }

            if(E_2_ISO)
            {
                printf("Generating Testdata for e_2_iso\n");
                for(int i = 0; i < NUM_E_ISO_TESTS; i ++)
                {
                    char* post = i == 0 ? "st" : i == 1 ? "nd" : "th";
                    printf("e_2_iso %2d%s test of %d\n", i + 1, post, NUM_E_ISO_TESTS);
                    proj_point_t opt_input[3];

                    proj_pt_zero(&(opt_input[0]));
                    proj_pt_zero(&(opt_input[1]));
                    proj_pt_zero(&(opt_input[2]));



                    f2p2_t temp;
                    f2p2_zero(&temp);

                    proj_point_t Q;
                    proj_point_t Q_redc;
                    proj_pt_zero(&Q);
                    proj_pt_zero(&Q_redc);
                    rand_fill_proj_pt(&Q, &mont_mod);
                    fp2_mul(&(Q.X), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(Q_redc.X));
                    fp2_mul(&(Q.Z), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(Q_redc.Z));


                    proj_point_t res;
                    proj_pt_zero(&res);
                    e_2_iso(&A24p_redc, &Q_redc, opt_input, 3, &mont_mod, STRAT_E2, (mont_curve_t*)&res);

                    fprintf(file, "e_2_iso\n");
                    proj_pt_fprintf(file, &Q);
                    proj_pt_fprintf(file, (proj_point_t*)&A24p);
                    fp_fprintf(file, &mont_mod);




                    fp2_t X_red;
                    fp2_t Y_red;
                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    f2p2_t t;
                    fp2_to_f2p2(&(res.X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    fp2_to_f2p2(&(res.Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);

                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    fp2_to_f2p2(&(opt_input[0].X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    fp2_to_f2p2(&(opt_input[0].Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);

                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    fp2_to_f2p2(&(opt_input[1].X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    fp2_to_f2p2(&(opt_input[1].Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);

                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    fp2_to_f2p2(&(opt_input[2].X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    fp2_to_f2p2(&(opt_input[2].Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);

                }
            }

            if(E_3_ISO)
            {
                printf("Generating Testdata for e_3_iso\n");
                for(int i = 0; i < NUM_E_ISO_TESTS; i ++)
                {
                    char* post = i == 0 ? "st" : i == 1 ? "nd" : "th";
                    printf("e_3_iso %2d%s test of %d\n", i + 1, post, NUM_E_ISO_TESTS);
                    proj_point_t opt_input[3];

                    proj_pt_zero(&(opt_input[0]));
                    proj_pt_zero(&(opt_input[1]));
                    proj_pt_zero(&(opt_input[2]));



                    f2p2_t temp;
                    f2p2_zero(&temp);

                    proj_point_t Q;
                    proj_point_t Q_redc;
                    proj_pt_zero(&Q);
                    proj_pt_zero(&Q_redc);
                    rand_fill_proj_pt(&Q, &mont_mod);
                    fp2_mul(&(Q.X), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(Q_redc.X));
                    fp2_mul(&(Q.Z), &r22, &temp);
                    REDCL2(&temp, &mont_mod, &(Q_redc.Z));


                    proj_point_t res;
                    proj_pt_zero(&res);
                    e_3_iso(&A24p_redc, &Q_redc, opt_input, 3, &mont_mod, STRAT_E3, (mont_curve_t*)&res);

                    fprintf(file, "e_3_iso\n");
                    proj_pt_fprintf(file, &Q);
                    proj_pt_fprintf(file, (proj_point_t*)&A24p);
                    fp_fprintf(file, &mont_mod);




                    fp2_t X_red;
                    fp2_t Y_red;
                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    f2p2_t t;
                    fp2_to_f2p2(&(res.X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    fp2_to_f2p2(&(res.Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);
                    /*
                    fp2_print(&(res.X));
                    LINE
                    fp2_print(&(res.Z));
                    LINE
                    fp2_print(&X_red);
                    LINE
                    fp2_print(&Y_red);
                    LINE
                    LINE
                    LINE
                    */
                    

                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    fp2_to_f2p2(&(opt_input[0].X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    fp2_to_f2p2(&(opt_input[0].Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);

                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    fp2_to_f2p2(&(opt_input[1].X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    fp2_to_f2p2(&(opt_input[1].Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);

                    fp2_zero(&X_red);
                    fp2_zero(&Y_red);
                    fp2_to_f2p2(&(opt_input[2].X), &t);
                    REDCL2(&t, &mont_mod, &X_red);
                    fp2_to_f2p2(&(opt_input[2].Z), &t);
                    REDCL2(&t, &mont_mod, &Y_red);
                    fp2_fprintf(file, &X_red);
                    fp2_fprintf(file, &Y_red);

                }
            }


        }

        fclose(file);
        return 0;
    }
    else
    {
        printf("Simple key exchange start\n");
        public_params_t params;
        
        params_translate_redec(&PARAMS, &params);
        printf("test\n\n\n\n\n");
        printf("size of p = %ld\n",fp_get_len(&(params.p)));
        printf("p is:\n");
        fp_print(&params.p);
        LINE
        
        fp_t secret_Alice;
        fp_t secret_Bob;

        public_key_t pup_Alice;
        public_key_t pup_Bob;

        fp_zero(&secret_Alice);
        fp_zero(&secret_Bob);
        fp_t temp_Sec;
        fp_zero(&temp_Sec);
        rand_fill_fp_pos(&temp_Sec);//, &params.p);
        //rand_fp_smaller_mod(&temp_Sec, &params.p);
        printf("Alice pre mod:\n");
        fp_print(&temp_Sec);
        LINE LINE
        fp_copy(&temp_Sec, &secret_Alice);
        fp_zero(&temp_Sec);
        rand_fp_smaller_mod(&temp_Sec, &params.p);
        //rand_fill_fp_pos(&temp_Sec);//, &params.p);
        //temp_Sec[WORDS-1] = 0;//temp_Sec[WORDS-1] >> 1;
        //temp_Sec[WORDS-2] = 0;
        //temp_Sec[WORDS-3] = 0;
        printf("Bob pre mod:\n");
        fp_print(&temp_Sec);
        LINE LINE
        fp_copy(&temp_Sec, &secret_Bob);

        printf("Alices secret is:\n");
        fp_print(&secret_Alice);
        LINE
        printf("Bobs secret is:\n");
        fp_print(&secret_Bob);
        LINE
        pthread_t alice;
        pthread_t bob;

        pub_calc_t pub_calc_alice = {.secret = &secret_Alice, 
                                    .parameters = &params,
                                    .calced_key = &pup_Alice};
        
        pub_calc_t pub_calc_bob =   {.secret = &secret_Bob, 
                                    .parameters = &params,
                                    .calced_key = &pup_Bob};


        int fail = pthread_create(&alice, NULL, &alice_pub_calc, (void*)&pub_calc_alice);
        if(fail)
        {
            printf("alice can not calculate\n");
            return -1;
        }
        fail = pthread_create(&bob, NULL, &bob_pub_calc, (void*)&pub_calc_bob);
        if(fail)
        {
            printf("bob can not calculate\n");
            return -1;
        }
        printf("Starting calculating public keys\n");
        pthread_join(alice, NULL);
        pthread_join(bob, NULL);

        // print public
        {
            printf("Alices pub key is:\n");
            fp2_t temp;
            fp2_zero(&temp);
            fp2_from_mont(&pup_Alice.x1, &params.p, &temp);
            fp2_print(&temp);
            LINE
            fp2_zero(&temp);
            fp2_from_mont(&pup_Alice.x2, &params.p, &temp);
            fp2_print(&temp);
            LINE
            fp2_zero(&temp);
            fp2_from_mont(&pup_Alice.x3, &params.p, &temp);
            fp2_print(&temp);
            LINE LINE LINE
            printf("Bobs pub key is:\n");
            fp2_zero(&temp);
            fp2_from_mont(&pup_Bob.x1, &params.p, &temp);
            fp2_print(&temp);
            LINE
            fp2_zero(&temp);
            fp2_from_mont(&pup_Bob.x2, &params.p, &temp);
            fp2_print(&temp);
            LINE
            fp2_zero(&temp);
            fp2_from_mont(&pup_Bob.x3, &params.p, &temp);
            fp2_print(&temp);
            LINE LINE LINE
        }

        fp2_t shared_Alice;
        fp2_t shared_Bob;

        sec_calc_t sec_calc_alice = {.secret = &secret_Alice, 
                                    .parameters = &params,
                                    .pub_key = &pup_Bob,
                                    .jinv = &shared_Alice};

        sec_calc_t sec_calc_bob = {.secret = &secret_Bob, 
                                    .parameters = &params,
                                    .pub_key = &pup_Alice,
                                    .jinv = &shared_Bob};

        fail = pthread_create(&alice, NULL, &alice_sec_calc, (void*)&sec_calc_alice);
        if(fail)
        {
            printf("alice can not calculate\n");
            return -1;
        }
        fail = pthread_create(&bob, NULL, &bob_sec_calc, (void*)&sec_calc_bob);
        if(fail)
        {
            printf("bob can not calculate\n");
            return -1;
        }                        
        pthread_join(alice, NULL);
        pthread_join(bob, NULL);
        /*
        printf("Alices secret is:\n");
        fp_print(&secret_Alice);
        LINE
        


        //Alice public key
        printf("Calculating Shared key of Alice\n");
        isogen_2(&params, &secret_Alice, &(params.p), &pup_Alice, &STRAT_E2);

        printf("Alices pub key is:\n");
        fp2_t temp;
        fp2_zero(&temp);
        fp2_from_mont(&pup_Alice.x1, &params.p, &temp);
        fp2_print(&temp);
        LINE
        fp2_zero(&temp);
        fp2_from_mont(&pup_Alice.x2, &params.p, &temp);
        fp2_print(&temp);
        LINE
        fp2_zero(&temp);
        fp2_from_mont(&pup_Alice.x3, &params.p, &temp);
        fp2_print(&temp);
        LINE LINE LINE


        //Bob public key
        printf("Bobs secret is:\n");
        fp_print(&secret_Bob);
        LINE
        printf("Calculating Shared key of Bob\n");
        isogen_3(&params, &secret_Bob, &(params.p), &pup_Bob, &STRAT_E3);
        fp2_zero(&temp);
        fp2_from_mont(&pup_Bob.x1, &params.p, &temp);
        fp2_print(&temp);
        LINE
        fp2_zero(&temp);
        fp2_from_mont(&pup_Bob.x2, &params.p, &temp);
        fp2_print(&temp);
        LINE
        fp2_zero(&temp);
        fp2_from_mont(&pup_Bob.x3, &params.p, &temp);
        fp2_print(&temp);
        LINE LINE LINE*/

        

        //shared secred Alice
        //printf("Calculating Shared secret of Alice\n");
        //isoex_2(&secret_Alice, &pup_Bob, &shared_Alice, &(params.p), &STRAT_E2);

        //shared secred Bob
        //printf("Calculating Shared secret of Bob\n");
        //isoex_3(&secret_Bob, &pup_Alice, &shared_Bob, &(params.p), &STRAT_E3);


        //fp2_t Alice_jinv;
        //fp2_t Bob_jinv;

        //fp2_from_mont(&shared_Alice, &(params.p), &Alice_jinv);
        //fp2_from_mont(&shared_Bob,   &(params.p), &Bob_jinv);

        printf("Shared secred Alice got:\n");
        fp2_print(&shared_Alice);
        LINE
        LINE
        printf("Shared secred Bob got:\n");
        fp2_print(&shared_Bob);
        LINE

        if(fp2_equ(&shared_Alice, &shared_Bob))
        {
            printf("#####################\n");
            printf("it's the same picture\n");
            printf("#####################\n");
        }
        else
        {
            printf("#####################\n");
            printf("not the same\n");
            printf("#####################\n");
        }

        return 0;
    }
}

void rand_fill_fp_pos(fp_t* a)
{
    int bits = 0;
    int i = 0;
    while(bits < BITS)
    {
        uint64_t lower  = (uint64_t)rand();
        uint64_t higher = (uint64_t)rand();
        (*a)[i] = ((higher << 32) | lower);
        bits += 64;
        i++;
    }
    //bits -= 64;
    int remaining_bits = BITS - bits + 64;
    uint64_t lower  = (uint64_t)rand();
    uint64_t higher = (uint64_t)rand();
    (*a)[i] = ((higher << 32) | lower) >> (64 - remaining_bits);
}

void rand_fill_fp_neg(fp_t* a)
{
    rand_fill_fp_pos(a);
    
    for(int i = WORDS - 1; i >= 0; i--)
    {
        uint64_t mask = (((uint64_t)1) << 63);   
        while(!((*a)[i] & mask) && (*a)[i] != 0xffffffffffffffffULL)
        {
            (*a)[i] =  (*a)[i] | mask;
            mask >>= 1;
        }
    }
}

void rand_fill_fp_2_input(fp_t* a, fp_t* b, int i)
{
    switch(i % 4){
            case 0:
                rand_fill_fp_pos(a);
                rand_fill_fp_pos(b);  
                break;
            case 1:
                rand_fill_fp_neg(a);
                rand_fill_fp_pos(b);  
                break;
            case 2:
                rand_fill_fp_pos(a);  
                rand_fill_fp_neg(b);
                break;
            case 3:
                rand_fill_fp_neg(a);  
                rand_fill_fp_neg(b);
                break;
            default:
                assert(0 && "no valid number fill\n");
                break;
        }
}

void rand_fill_fp_divisor(fp_t* a)
{
    int bits = 0;
    int i = 0;
    while(bits < BITS)
    {

        uint64_t lower  = (uint64_t)rand();
        uint64_t higher = (uint64_t)rand();
        (*a)[i] = ((higher << 32) | lower);
        bits += 64;
        i++;
    }
    //bits -= 64;
    if(!(rand() % 2))
    {
        int remaining_bits = BITS - bits + 64;
        uint64_t lower  = (uint64_t)rand();
        uint64_t higher = (uint64_t)rand();
        (*a)[i] = ((higher << 32) | lower) >> (64 - remaining_bits);
    }
}

void rand_fill_fp_divisor_pos_neg(fp_t* a, int i)
{
    fp_t a_pre;
    fp_zero(&a_pre);
    rand_fill_fp_divisor(&a_pre);
    if(i % 2)
    {
        fp_additive_inverse(&a_pre, a);
    }
    else
    {
        fp_copy(&a_pre, a);
    }
}

void rand_fill_f2p_pos(f2p_t* a)
{
    int bits = 0;
    int i = 0;
    while(bits < 2*BITS)
    {
        uint64_t lower  = (uint64_t)rand();
        uint64_t higher = (uint64_t)rand();
        (*a)[i] = ((higher << 32) | lower);
        bits += 64;
        i++;
    }
    //bits -= 64;
    int remaining_bits = BITS - bits + 64;
    uint64_t lower  = (uint64_t)rand();
    uint64_t higher = (uint64_t)rand();
    (*a)[i] = ((higher << 32) | lower) >> (64 - remaining_bits);
}

void rand_fill_f2p_neg(f2p_t* a)
{
    rand_fill_f2p_pos(a);
    
    for(int i = DWORDS - 1; i >= 0; i--)
    {
        uint64_t mask = (((uint64_t)1) << 63);   
        while(!((*a)[i] & mask) && (*a)[i] != 0xffffffffffffffffULL)
        {
            (*a)[i] =  (*a)[i] | mask;
            mask >>= 1;
        }
    }
}

void rand_fill_f2p_2_input(f2p_t* a, f2p_t* b, int i)
{
    switch(i % 4){
            case 0:
                rand_fill_f2p_pos(a);
                rand_fill_f2p_pos(b);  
                break;
            case 1:
                rand_fill_f2p_neg(a);
                rand_fill_f2p_pos(b);  
                break;
            case 2:
                rand_fill_f2p_pos(a);  
                rand_fill_f2p_neg(b);
                break;
            case 3:
                rand_fill_f2p_neg(a);  
                rand_fill_f2p_neg(b);
                break;
            default:
                assert(0 && "no valid number fill\n");
                break;
        }
}

void rand_2k_num(fp_t* a)
{
    int shift = rand();
    shift = shift % (WORDS * 64);
    fp_zero(a);
    (*a)[0] = 1;
    fp_lshift(a, shift);
}

void rand_2k_num_f2p(f2p_t* a)
{
    int shift = rand();
    shift = shift % (2 * WORDS * 64);
    f2p_zero(a);
    (*a)[0] = 1;
    f2p_lshift(a, shift);
}

void rand_fill_fp_for_f2p_pos(f2p_t* a)
{
    int bits = 0;
    int i = 0;
    while(bits < BITS)
    {
        uint64_t lower  = (uint64_t)rand();
        uint64_t higher = (uint64_t)rand();
        (*a)[i] = ((higher << 32) | lower);
        bits += 64;
        i++;
    }
    //bits -= 64;
    int remaining_bits = BITS - bits + 64;
    uint64_t lower  = (uint64_t)rand();
    uint64_t higher = (uint64_t)rand();
    (*a)[i] = ((higher << 32) | lower) >> (64 - remaining_bits);
}

void rand_fill_fp_for_f2p_neg(f2p_t* a)
{
    rand_fill_fp_for_f2p_pos(a);
    
    for(int i = DWORDS - 1; i >= 0; i--)
    {
        uint64_t mask = (((uint64_t)1) << 63);   
        while(!((*a)[i] & mask) && (*a)[i] != 0xffffffffffffffffULL)
        {
            (*a)[i] =  (*a)[i] | mask;
            mask >>= 1;
        }
    }
}

void rand_fill_fp_for_f2p_2_input(f2p_t* a, f2p_t* b, int i)
{
    switch(i % 4)
    {
        case 0:
            rand_fill_fp_for_f2p_pos(a);
            rand_fill_fp_for_f2p_pos(b);
            break;
        case 1:
            rand_fill_fp_for_f2p_neg(a);
            rand_fill_fp_for_f2p_pos(b);
            break;
        case 2:
            rand_fill_fp_for_f2p_pos(a);
            rand_fill_fp_for_f2p_neg(b);
            break;
        case 3:
            rand_fill_fp_for_f2p_neg(a);
            rand_fill_fp_for_f2p_neg(b);
            break;
        default:
            assert(0 && "f2p_rand fill invalid configuration");
    }
}

void rand_fill_f2p_divisor(f2p_t* a)
{
   int bits = 0;
    int i = 0;
    while(bits < 2*BITS)
    {

        uint64_t lower  = (uint64_t)rand();
        uint64_t higher = (uint64_t)rand();
        (*a)[i] = ((higher << 32) | lower);
        bits += 64;
        i++;
    }
    //bits -= 64;
    if(!(rand() % 2))
    {
        int remaining_bits = BITS - bits + 64;
        uint64_t lower  = (uint64_t)rand();
        uint64_t higher = (uint64_t)rand();
        (*a)[i] = ((higher << 32) | lower) >> (64 - remaining_bits);
    } 
}

void rand_fill_f2p_divisor_pos_neg(f2p_t* a, int i)
{
    f2p_t a_pre;
    f2p_zero(&a_pre);
    rand_fill_f2p_divisor(&a_pre);
    if(i % 2)
    {
        f2p_additive_inverse(&a_pre, a);
    }
    else
    {
        f2p_copy(&a_pre, a);
    }
}

void rand_fill_fp2(fp2_t* a)
{
    rand_fill_fp_pos(&(a->real));
    rand_fill_fp_pos(&(a->img));


}

void rand_fill_fp2_2_input(fp2_t* a, fp2_t* b, int i)
{
    switch(i % 16)
    {
        case 0:
            rand_fill_fp_pos(&(a->real));
            rand_fill_fp_pos(&(a->img));  
            rand_fill_fp_pos(&(b->real));
            rand_fill_fp_pos(&(b->img));
            break;
        case 1:
            rand_fill_fp_neg(&(a->real));
            rand_fill_fp_pos(&(a->img));  
            rand_fill_fp_pos(&(b->real));
            rand_fill_fp_pos(&(b->img));
            break;
        case 2:
            rand_fill_fp_pos(&(a->real));
            rand_fill_fp_neg(&(a->img));  
            rand_fill_fp_pos(&(b->real));
            rand_fill_fp_pos(&(b->img));
            break;
        case 3:
            rand_fill_fp_neg(&(a->real));
            rand_fill_fp_neg(&(a->img));  
            rand_fill_fp_pos(&(b->real));
            rand_fill_fp_pos(&(b->img));
            break;
        case 4:
            rand_fill_fp_pos(&(a->real));
            rand_fill_fp_pos(&(a->img));  
            rand_fill_fp_neg(&(b->real));
            rand_fill_fp_pos(&(b->img));
            break;
        case 5:
            rand_fill_fp_neg(&(a->real));
            rand_fill_fp_pos(&(a->img));  
            rand_fill_fp_neg(&(b->real));
            rand_fill_fp_pos(&(b->img));
            break;
        case 6:
            rand_fill_fp_pos(&(a->real));
            rand_fill_fp_neg(&(a->img));  
            rand_fill_fp_neg(&(b->real));
            rand_fill_fp_pos(&(b->img));
            break;
        case 7:
            rand_fill_fp_neg(&(a->real));
            rand_fill_fp_neg(&(a->img));  
            rand_fill_fp_neg(&(b->real));
            rand_fill_fp_pos(&(b->img));
            break;
        case 8:
            rand_fill_fp_pos(&(a->real));
            rand_fill_fp_pos(&(a->img));  
            rand_fill_fp_pos(&(b->real));
            rand_fill_fp_neg(&(b->img));
            break;
        case 9:
            rand_fill_fp_neg(&(a->real));
            rand_fill_fp_pos(&(a->img));  
            rand_fill_fp_pos(&(b->real));
            rand_fill_fp_neg(&(b->img));
            break;
        case 10:
            rand_fill_fp_pos(&(a->real));
            rand_fill_fp_neg(&(a->img));  
            rand_fill_fp_pos(&(b->real));
            rand_fill_fp_neg(&(b->img));
            break;
        case 11:
            rand_fill_fp_neg(&(a->real));
            rand_fill_fp_neg(&(a->img));  
            rand_fill_fp_pos(&(b->real));
            rand_fill_fp_neg(&(b->img));
            break;
        case 12:
            rand_fill_fp_pos(&(a->real));
            rand_fill_fp_pos(&(a->img));  
            rand_fill_fp_neg(&(b->real));
            rand_fill_fp_neg(&(b->img));
            break;
        case 13:
            rand_fill_fp_neg(&(a->real));
            rand_fill_fp_pos(&(a->img));  
            rand_fill_fp_neg(&(b->real));
            rand_fill_fp_neg(&(b->img));
            break;
        case 14:
            rand_fill_fp_pos(&(a->real));
            rand_fill_fp_neg(&(a->img));  
            rand_fill_fp_neg(&(b->real));
            rand_fill_fp_neg(&(b->img));
            break;
        case 15:
            rand_fill_fp_neg(&(a->real));
            rand_fill_fp_neg(&(a->img));  
            rand_fill_fp_neg(&(b->real));
            rand_fill_fp_neg(&(b->img));
            break;
        default:
            assert(0 && "invalid input configuration\n");
    }
}

void rand_fp_smaller_mod(fp_t* a, fp_t* mod)
{
    rand_fill_fp_pos(a);
    fp_t mask;
    fp_t one;
    fp_zero(&one);
    fp_zero(&mask);
    one[0] = 1;
    int i = WORDS - 1;
    for(; i >= 0;)
    {
        if((*mod)[i] == 0)
        {
            mask[i] = 0xffffffffffffffffULL;
            i--;
            continue;
        }
        break;
    }
    uint64_t m = (uint64_t)1 << 63;
    //printf("mask:     %lx\n", m);
    while(m > 0)
    {
        //printf("1\n");
        if(!(m & (*mod)[i]))
        {
            //printf("mod:      %lx\n",(*mod)[i]);
            mask[i] |= m;
            m = m >> 1;
            //printf("new_mask: %lx\n", m);
            continue;
        }
        mask[i] |= m;
        break;
    }

    fp_t a_pre;
    for(int i = 0; i < WORDS; i++)
    {
        a_pre[i] = (*a)[i] & (~mask[i]);
    }
    fp_sub(&a_pre, &one, a);
}

void rand_fp2_smaller_mod(fp2_t* a, fp_t* mod)
{
    rand_fp_smaller_mod(&(a->real), mod);
    rand_fp_smaller_mod(&(a->img),  mod);
}

void rand_fill_proj_pt(proj_point_t* p, fp_t* mod)
{
    rand_fp2_smaller_mod(&(p->X), mod);
    rand_fp2_smaller_mod(&(p->Z), mod);
}

void* alice_pub_calc(void* params)
{
    pub_calc_t* par = (pub_calc_t*)params;

    isogen_2((par->parameters), (par->secret), &(par->parameters->p), (par->calced_key), &STRAT_E2);

    return NULL;
}


void* bob_pub_calc(void* params)
{
    pub_calc_t* par = (pub_calc_t*)params;

    isogen_3((par->parameters), (par->secret), &(par->parameters->p), (par->calced_key), &STRAT_E3);

    return NULL;
}


void* alice_sec_calc(void* params)
{
    sec_calc_t* par = (sec_calc_t*)params;
    isoex_2(par->secret, par->pub_key, par->jinv, &(par->parameters->p), &STRAT_E2);
    return NULL;

}

void* bob_sec_calc(void* params)
{
    sec_calc_t* par = (sec_calc_t*)params;
    isoex_3(par->secret, par->pub_key, par->jinv, &(par->parameters->p), &STRAT_E3);
    return NULL;

}