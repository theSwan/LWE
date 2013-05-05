#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gmp.h>
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "libbgv.h"

double bgv_get_dvn()
{
	return dvn;   
}

void gen_q(fmpz_t q, long len)
{
        mpz_t tmp, hold;
        mpz_init(tmp);
        mpz_init(hold);
        mpz_set_ui(hold, 1);
        mpz_mul_2exp(tmp, hold, len);
        mpz_nextprime(hold, tmp);
        fmpz_set_mpz(q, hold);
        mpz_clear(tmp);
        mpz_clear(hold);
}

void hcrypt_random(mpz_t tmp)
{
	FILE *fp;
	fp = fopen("/dev/urandom", "rb");
        int len = 9;
	if (fp) {
		unsigned char *bytes;
		bytes = (unsigned char *) malloc (len * sizeof(unsigned char));
              
                if (fread(bytes, sizeof(unsigned char), len, fp)) {
                        mpz_import(tmp, len, 1, sizeof(unsigned char), 0, 0, bytes);
                }
                
                else {
                        printf("file read error\n");
                }
                
		fclose(fp);
		free(bytes);
	}
        
        else {
                printf("random number generation error\n");
        }
}

void guassian_vec(fmpz_mat_t mat, long len)
{
        double tdvn = bgv_get_dvn();
	long a = (long)ceil(-10*tdvn);
	long b = (long)floor(+10*tdvn);
	long x, i;
	double p;
	mpz_t randseed;
	mpz_init(randseed);
	hcrypt_random(randseed);
	unsigned long int useed = mpz_get_ui(randseed);
	srand(useed);
	for( i = 0 ; i < len ; i++) {
                x = rand()%(b - a) + a;
                fmpz_set_si(fmpz_mat_entry(mat, i, 0), x);
	}
	mpz_clear(randseed);
}


void unif_mat(fmpz_mat_t mat, fmpz_t q, long row, long col)
{
        long i, x, j;
	flint_rand_t state;
        flint_randinit(state);
        
	for( i = 0 ; i < row; i++) {
                for( j = 0 ; j < col; j++ ){
                        fmpz_randtest_mod_signed(fmpz_mat_entry(mat, i, j), state, q);
                }
	}
	flint_randclear(state);
}

void fmpz_smod(fmpz_t num, fmpz_t q)
{
        fmpz_mod(num, num, q);
        fmpz_t tmp;
        fmpz_init(tmp);
        fmpz_cdiv_q_si(tmp, q, 2);
        if(fmpz_cmp(num, tmp) > 0) {
                fmpz_sub(num, num, q);
        }
        fmpz_clear(tmp);
}

void fmpz_smod_ui(fmpz_t num, long mod)
{
        fmpz_mod_ui(num, num, mod);
        long tmp = mod/2 + 1;
        if(fmpz_cmp_si(num, tmp) > 0){
                fmpz_sub_ui(num, num, mod);
        }
}
