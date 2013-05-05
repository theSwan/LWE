#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gmp.h>
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "libbgv.h"

param_node_t *e_setup(long miu, long lamda, long b, param_node_t *param)
{
        param = param_node_init(param);
	gen_q(param->q, miu);
	fmpz_t tmp;
	fmpz_init(tmp);
	fmpz_fdiv_q_si(tmp, param->q, bgv_get_bigb());
	long prod, mi;
	prod = lamda * fmpz_flog_ui(tmp, 2);
        fmpz_set_si(tmp, prod);
        mi = fmpz_flog_ui(tmp, 2);
        
        prod = pow(2, mi);
	if(b == 1) {
		param->n = prod;
	}  /* LWE */
	else {
		param->n = 1;
	} /* RLWE */
        
	param->bign = ceil((2 * param->n + 1) * fmpz_flog_ui(param->q, 2));
        
	fmpz_clear(tmp);
	return param;
}


void e_skeygen(fmpz_mat_t sk, param_node_t *param)
{
        guassian_vec(sk, (param->n + 1));
        fmpz_set_si(fmpz_mat_entry(sk, 0, 0), 1);
        long i;
        for( i = 1 ; i <= param->n ; i++ ) {
                fmpz_smod(fmpz_mat_entry(sk, i, 0), param->q);
        }
}

void e_pkeygen(fmpz_mat_t pk, param_node_t *param, fmpz_mat_t sk)
{
        fmpz_mat_t ppk, ee, bb, ss, tmp, tmp1;
        fmpz_mat_init(ppk, param->bign, param->n);
        fmpz_mat_init(ee, param->bign, 1);
        fmpz_mat_init(bb, param->bign, 1);
        fmpz_mat_init(ss, param->n, 1);
        
        long i, j;
        for( i = 0 ; i < param->n ; i++ ) {
                fmpz_set(fmpz_mat_entry(ss, i, 0), fmpz_mat_entry(sk, i+1, 0 ));
        }
        guassian_vec(ee, param->bign);
        unif_mat(ppk, param->q, param->bign, param->n);
        
        fmpz_mat_init(tmp, param->bign, 1);
        fmpz_mat_init(tmp1, param->bign, 1);
        fmpz_mat_mul(tmp, ppk, ss);
        fmpz_mat_scalar_mul_si(tmp1, ee, 2);
        fmpz_mat_add(bb, tmp, tmp1);
        for( i = 0 ; i < param->bign ; i++ ) {
                fmpz_set(fmpz_mat_entry(pk, i, 0), fmpz_mat_entry(bb, i, 0));
                fmpz_smod(fmpz_mat_entry(pk, i, 0), param->q);
        }
        for( i = 0 ; i < param->bign ; i++ ) {
                for( j = 1; j <= param->n; j++ ){
                        fmpz_neg(fmpz_mat_entry(pk, i, j), fmpz_mat_entry(ppk, i, j-1));
                        fmpz_smod(fmpz_mat_entry(pk, i, j), param->q);
                }
        }
        
        fmpz_mat_clear(tmp);
        fmpz_mat_clear(tmp1);
        fmpz_mat_clear(ee);
        fmpz_mat_clear(ss);
        fmpz_mat_clear(bb);
        fmpz_mat_clear(ppk);
}

void e_encrypt(fmpz_mat_t ct, param_node_t *param, fmpz_mat_t pk, fmpz_t ms)
{
        long i, j;
        
        fmpz_mat_t mm, rr, tmp, tmp1;
        fmpz_mat_init(mm, 1 + param->n, 1);
        fmpz_mat_init(rr, param->bign, 1);
        fmpz_mat_init(tmp, 1 + param->n, 1);
        fmpz_mat_init(tmp1, 1 + (param->n), param->bign);
        for( i = 0 ; i < 1 + param->n ; i++ ) {
                for( j = 0; j < param->bign; j++ ){
                        fmpz_set(fmpz_mat_entry(tmp1, i, j), fmpz_mat_entry(pk, j, i));
                }
        }
        fmpz_mat_zero(mm);
        fmpz_set(fmpz_mat_entry(mm, 0, 0), ms);
        
        fmpz_t t;
        fmpz_init_set_ui(t,2);
        unif_mat(rr, t, param->bign, 1);
        fmpz_mat_mul(tmp, tmp1, rr);
        fmpz_mat_add(ct, mm, tmp);
        
        for( i = 0; i < param->n + 1 ; i++) {
                fmpz_smod(fmpz_mat_entry(ct, i, 0), param->q);
        }
        
        fmpz_clear(t);
        fmpz_mat_clear(tmp);
        fmpz_mat_clear(tmp1);
        fmpz_mat_clear(mm);
        fmpz_mat_clear(rr);
}

void e_decrypt(fmpz_t ms, param_node_t *param, fmpz_mat_t sk, fmpz_mat_t ct)
{
        fmpz_t tmp;
        fmpz_init(tmp);
        fmpz_zero(ms);
        
        long i;
        
        for( i = 0 ; i < param->n + 1 ; i++) {
                fmpz_mul(tmp, fmpz_mat_entry(ct, i, 0), fmpz_mat_entry(sk, i, 0));
                fmpz_add(ms, ms, tmp);
        }
      //  fmpz_print(ms);
      //  printf("\n");
        fmpz_smod(ms, param->q);
       // fmpz_print(ms);
       // printf("\n");
        fmpz_smod_ui(ms, 2);
        fmpz_clear(tmp);
}
