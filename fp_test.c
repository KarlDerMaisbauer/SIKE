

#include <stdio.h>
#include <stdint.h>

#include "fp.h"
#include "fp_helper.h"
#include "montgomory_redc.h"
#include "montgomery_curve.h"

#include "sike_encodings_helper.h"
#include "sike_encodings.h"
#include "params.h"
#include "point_arithmetic.h"


#include <limits.h>
int main(void){
    //b = 122;
    //printf("%d\n",b);
    /*
    fp_t e       ;//= fp_init();
    fp_t p       ;//= fp_init();
    fp_t er      ;//= fp_init();
    fp_t er_inv  ;//= fp_init();
    fp_t er_sub  ;//= fp_init();
    fp_t three   ;//= fp_init();
    fp_t one     ;//= fp_init();
    fp_t two     ;
    fp_t div_erg ;//= fp_init();
    f2p_t res    ;
    //fp_t* a       = fp_init();
    //fp_t* mod     = fp_init();
    //fp_t* mod_erg = fp_init();
    for(int i = 0; i < WORDS; i++)
    {
        e[i] = 0xFFFFFFFFFFFFFFFF;
        p[i] = 0x111111111111111F;
    }
    fp_zero(&one);
    fp_zero(&two);
    fp_zero(&three);
    one[2]   = 1;
    two[2]   = 0x0000000000000003;
    three[2] = 3;

    //add
    printf("Addition:\n");
    fp_add(&e, &p, &er);
    fp_print(&e);
    printf(" +\n");
    fp_print(&p);
    printf(" = \n");
    fp_print(&er);
    printf(" \n\n\n");

    printf("additive invers:\n");
    fp_print(&er);
    printf("\n");
    fp_additive_inverse(&er, &er_inv);
    fp_print(&er_inv);
    printf(" \n\n\n");
*/
    // sub
/*
    fp_t mone;
    fp_t six;
    fp_t sub_res;
    fp_zero(&six);
    fp_zero(&sub_res);
    six[0] = 6;
    for(int i = 0; i < WORDS; i++)
    {
        mone[i] = 0xFFFFFFFFFFFFFFFF;
    }
    printf("Substraction:\n");
    fp_sub(&mone, &six, &sub_res);
    fp_print(&mone);
    printf(" -\n");
    fp_print(&six);
    printf(" =\n");
    fp_print(&sub_res);
    printf(" \n\n\n");
*/
    
/*
    // mul
    printf("Multiplication:\n");
    fp_mul(&two, &e, &res);
    fp_print(&two);
    printf(" *\n");
    fp_print(&e);
    printf(" =\n");
    f2p_print(&res);
    printf(" \n\n\n");

    //left shift
    printf("left shift large:\n");
    f2p_print(&res);
    printf("\n");
    f2p_lshift(&res, 10);
    f2p_print(&res);
    printf(" \n\n");

    printf("left shift small:\n");
    fp_print(&er);
    printf("\n");
    fp_lshift(&er, 4);
    fp_print(&er);
    printf(" \n\n\n");
*/
    //div
  /*  
    printf("Division:\n");
    fp_t three   ;//= fp_init();
    fp_t one     ;//= fp_init();
    fp_t div_erg ;//= fp_init();

    fp_zero(&one);
    fp_zero(&three);
    one[0] = 17;
    three[0]=3663;
    fp_div(&three, &one, &div_erg);
    fp_print(&three);
    printf(" /\n");
    fp_print(&one);
    printf(" =\n");
    fp_print(&div_erg);
    printf(" \n\n");


    //mod
    printf("Modulo:\n");
    fp_t a;
    fp_t mod;
    fp_t mod_erg;
    fp_zero(&a);
    fp_zero(&mod);
    fp_zero(&mod_erg);
    a[8] = 2;
    a[0] = 1;
    mod[5] = 2;

    fp_mod(&a, &mod, &mod_erg);
    fp_print(&a);
    printf(" %%\n");
    fp_print(&mod);
    printf(" =\n");
    fp_print(&mod_erg);
    printf("\n");*/
    //int64_t dssd = small_div(-433, 11);
    //printf("%d\n",dssd);

        

    //gcd
/*
    printf("GCD:\n");
    fp_t nine;
    fp_t seven;

    fp_zero(&nine);
    fp_zero(&seven);

    nine[0] = 12;
    seven[0] = 17;

    fp_t ninemul;
    fp_t sevenmul;


    fp_zero(&ninemul);
    fp_zero(&sevenmul);
    printf("\n");
    gcd(&seven, &nine, &sevenmul, &ninemul);

    f2p_t sevmul;
    f2p_t ninmul;
    f2p_t erg;
    fp_mul(&seven, &sevenmul, &sevmul); 
    fp_mul(&nine, &ninemul, &ninmul); 
    f2p_add(&sevmul, &ninmul, &erg);
    fp_print(&ninemul);
    printf("\n");
    fp_t mod_ver;
    fp_mod(&ninemul, &seven, &mod_ver);
    fp_print(&ninemul);
    printf("%%\n");
    fp_print(&sevenmul);
    printf("=\n");
    fp_print(&mod_ver);
    printf("\n");
    printf("erg\n");
    f2p_print(&erg);
    printf("\n");
    //init(&seven);
    printf("\n\n\n");*/
/*

    // k2 div
    printf("2k div:\n");
    f2p_t A;
    f2p_t div;
    f2p_t divkres;
    f2p_zero(&A);
    f2p_zero(&div);
    f2p_zero(&divkres);
    div[5] = 1;
    A[2] = 2345;
    A[6] = 324;
    f2p_div_2k(&A, &div, &divkres);
    f2p_print(&A);
    printf("/\n");
    f2p_print(&div);
    printf("=\n");
    f2p_print(&divkres);
    printf("\n\n\n");



    // k2 mod
    printf("2k mod:\n");
    fp_t modularednum;
    fp_t modulo;
    fp_t modulo_res;

    fp_zero(&modulo);
    modulo[0] = 1;
    fp_lshift(&modulo, 97);
    modularednum[0] = 352456464475364234;

    fp_mod_2k(&modularednum, &modulo, &modulo_res);
    fp_print(&modularednum);
    printf("%%\n");
    fp_print(&modulo);
    printf("=\n");
    fp_print(&modulo_res);
    printf("\n\n\n");
    fp_t modnum_inv;
    fp_additive_inverse(&modularednum, &modnum_inv);
    fp_mod_2k(&modnum_inv, &modulo, &modulo_res);
    fp_print(&modnum_inv);
    printf("%%\n");
    fp_print(&modulo);
    printf("=\n");
    fp_print(&modulo_res);
    printf("\n\n\n");
    */

/*
    // montomory
    printf("montgomory mod:\n");

    fp_t threee;
    fp_t five;
    fp_t mont_mod;
    fp_t mont_res;
    fp_zero(&threee);
    fp_zero(&five);
    fp_zero(&mont_mod);

    mont_mod[0] = 19531;
    threee[0] = 13451;
    five[0] = 95784;

    init(&mont_mod);


    MODMUL(&threee, &five, &mont_mod, &mont_res);

    fp_print(&threee);
    printf("*\n");
    fp_print(&five);
    printf("%%\n");
    fp_print(&mont_mod);
    printf("=\n");
    fp_print(&mont_res);
    printf("\n\n\n");
    */

/*
    printf("divm:\n");

    fp_t ad;
    fp_t a_inv;
    fp_t mod_div;
    fp_t a_inv_pre;
    fp_t mod_inv;
    fp_zero(&ad);
    fp_zero(&a_inv);
    fp_zero(&mod_div);

    ad[0] = 12;
    mod_div[0] = 17;
    gcd(&a, &mod_div, &a_inv_pre, &mod_inv);
    fp_mod(&a_inv_pre, &mod_div, &a_inv);
    fp_t res_divm;
    fp_divm(&ad, &a_inv, &mod_div, &res_divm);
    fp_print(&ad);
    printf("/\n");
    fp_print(&a_inv);
    printf("%%\n");
    fp_print(&mod_div);
    printf("=\n");
    fp_print(&res_divm);
    printf("\n\n\n");*/

/*
    printf("16char to hex:\n");

    const char* strtest = "3563463634663456";
    uint64_t strerg = c16_to_hex(strtest, 0);
    printf("should : %16s\n", strtest);
    printf("is     : %16lx\n",strerg);

    printf("\n\n\n");
    printf("end char to hex:\n");

    const char* strtestend = "0x35634663456";
    uint64_t strergend = c_end_to_hex(strtestend, 0);
    printf("should : %16s\n", strtestend);
    printf("is     : %16lx\n",strergend);
    printf("\n\n\n");

    printf("ostoi:\n");
    const char* readnumhex = "0x2341F271773446CFC5FD681C520567BC65C783158AEA3FDC1767AE2FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
    fp_t readnum;
    fp_zero(&readnum);
    ostoi(readnumhex, &readnum);

    printf("sould:\n");
    printf("%192s\n", readnumhex);

    for(int i = WORDS-1; i >= 0; i--)
    {
        printf("%16lx", readnum[i]);
    }

    printf("\n\n\n");


    printf("ostofp:\n");

    fp_t faulty;
    fp_t faulty_tol;
    fp_zero(&faulty);
    fp_zero(&faulty_tol);
    fp_t ostomod;
    fp_zero(&ostomod);
    ostomod[5] = 234;
    const char* faultystr = "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
    const char* faultystr_tol = "0x0000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";

    if(ostofp(faultystr, &ostomod, &faulty))
        printf("should be and is faulty\n");
    else
        printf("should be but is not faulty\n");
    fp_print(&faulty);
    printf("\n\n\n");
    if(ostofp(faultystr_tol, &ostomod, &faulty_tol))
        printf("should be and is faulty\n");
    else
        printf("should be but is not faulty\n");
    fp_print(&faulty_tol);
    printf("\n\n\n");
    


    printf("itoos\n");
    fp_t hexnum;
    fp_zero(&hexnum);
    ostofp(readnumhex, &faulty, &hexnum);
    char* test = itoos(&hexnum);
    for(int i = WORDS-1; i >= 0; i--)
    {
        printf("%16lx", hexnum[i]);
    }
    printf("\n%192s\n",test);*/

    fp_t mod;
    fp_zero(&mod);
    ostoi(SIKEp434.p, &mod);
    //fp_print(&mod);
    //printf("bef init\n");
    init(&mod);
    //printf("aft init\n");
    mont_curve_t A;
    curve_zero(&A);
    fp2_copy(&six_mont, &A.A);
    fp2_copy(&one_mont, &A.C);

    public_params_t test;
    //printf("\n");
    params_translate(&SIKEp434, &test);

    public_params_t translated;

    params_translate_redec(&SIKEp434, &translated);

    mont_curve_t A_new;
    fp2_t erg;

    fp_t one;
    fp_zero(&one);
    one[0] = 1;

    Ladder3p(&translated.xP2, &translated.xQ2, &translated.xR2, &one, &erg, &A_new, &mod);
    fp2_t erg_red;

    REDC2(&erg, &mod, &erg_red);
    
    fp_print(&erg_red.real);
    printf("\n");
    fp_print(&erg_red.img);
    /*
    printf("\np:\n");
    fp_print(&test.p);
    printf("\n");

    printf("\ne2:\n");
    fp_print(&test.e2);
    printf("\n");

    printf("\ne3:\n");
    fp_print(&test.e3);
    printf("\n");

    //xQ2
    printf("\nQ2:\n");
    printf("\n");
    fp_print(&test.xQ2.real);
    printf("\n");

    fp_print(&test.xQ2.img);
    printf("\n");


    //xP2
    printf("\nP2:\n");
    printf("\n");
    fp_print(&test.xP2.real);
    printf("\n");

    fp_print(&test.xP2.img);
    printf("\n");

    //xR2
    printf("\nR2:\n");
    printf("\n");
    fp_print(&test.xR2.real);
    printf("\n");

    fp_print(&test.xR2.img);
    printf("\n");

    printf("montgomory reduced params\n");*/
    
    return 0;
}


