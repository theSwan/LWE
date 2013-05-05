#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gmp.h>
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "libbgv.h"


long bgv_get_bigb()
{
        return bigb;
}
param_node_t *param_node_init(param_node_t *pnt)
{
	pnt = (param_node_t *)malloc(sizeof(param_node_t));
	pnt->next = NULL;
	pnt->n = pnt->bign = 0;
	fmpz_init(pnt->q);
	return pnt;
}

void powers(fmpz_mat_t po, fmpz_mat_t x, fmpz_t qq, long xrow)
{
	long len = fmpz_clog_ui(qq, 2);
	long qrow = xrow * len;
	long i, j;
	for( i = 0 ; i < xrow ; i++) {
		fmpz_set(fmpz_mat_entry(po, i, 0), fmpz_mat_entry(x, i, 0));
	}
	for( i = 1 ; i < len ; i++) {
		for( j = 0 ; j < xrow ; j++) {
			fmpz_mul_si(fmpz_mat_entry(po, j + i * xrow, 0), fmpz_mat_entry(po, j + (i-1)*xrow, 0), 2);
			fmpz_smod(fmpz_mat_entry(po, j + i * xrow, 0), qq);
		}
	}
}

void bitdecomp(fmpz_mat_t bits, fmpz_mat_t x, fmpz_t qq, long xrow)
{
        fmpz_t t;
        fmpz_init_set_ui(t,2);
	long len = fmpz_clog(qq, t);
	long i, j, k;
        fmpz_t hold, negval;
        fmpz_init(hold);
        fmpz_init(negval);
        for(i = 0; i < xrow; i++) {
                fmpz_set(hold, fmpz_mat_entry(x, i, 0));
                k = 0;
                if(fmpz_cmp_si(hold, 0)<0) {
                        fmpz_neg(hold, hold);
                        while ( !fmpz_is_zero(hold) ) {
                                fmpz_mod(negval, hold, t);
                                j = i + k * xrow;
                                fmpz_neg(fmpz_mat_entry(bits, j, 0), negval);
                                fmpz_tdiv_q(hold, hold, t);
                                k++;
                        }
                }
                else {
                        while ( !fmpz_is_zero(hold) ) {
                                j = i + k * xrow;
                                fmpz_mod(fmpz_mat_entry(bits, j, 0), hold, t);
                                fmpz_tdiv_q(hold, hold, t);
                                k++;
                        }
                }
        }
        
        fmpz_clear(t);
        fmpz_clear(negval);
}

void vec_tensor(fmpz_mat_t tensor, fmpz_mat_t x, fmpz_t qq, long xrow)
{
        long i, j;
        for( i = 0 ; i < xrow ; i++ ) {
                for( j = 0 ; j < xrow ; j++ ){
                        fmpz_mul(fmpz_mat_entry(tensor,j+i*xrow,0), fmpz_mat_entry(x,i,0), fmpz_mat_entry(x,j,0));
        		fmpz_smod(fmpz_mat_entry(tensor,j+i*xrow,0), qq);
                }
        }
}


void switchkeygen(fmpz_mat_t mapb, fmpz_mat_t s1, fmpz_mat_t s2, fmpz_t qq, long n1, long n2)
{
	fmpz_mat_t sp1;
	param_node_t *param;
	param = (param_node_t *)malloc(sizeof(param_node_t));
	long i;
	param->n = n2 - 1;
	param->bign = n1 * fmpz_clog_ui(qq, 2);
	fmpz_init_set(param->q, qq);
	param->next = NULL;
	e_pkeygen(mapb, param, s2);
	long len = fmpz_clog_ui(qq, 2);
	long qrow = n1 * len;
	powers(sp1, s1, qq, n1);
	for( i = 0 ; i < param->bign ; i++) {
		fmpz_add(fmpz_mat_entry(mapb, i, 0), fmpz_mat_entry(mapb, i, 0), fmpz_mat_entry(sp1, i, 0));
                fmpz_smod(fmpz_mat_entry(mapb, i, 0), qq);
        }
	fmpz_mat_clear(sp1);
	free(param);
}

void scale(fmpz_mat_t c2, fmpz_mat_t c1, fmpz_t qq, fmpz_t pp, long row)
{
        long i;
        
        fmpz_t coeff, tmp, tmp1, tmp2, tmp3;
        fmpz_init(coeff);
        fmpz_init(tmp);
        fmpz_init(tmp1);
        fmpz_init(tmp2);
        fmpz_init(tmp3);
        int iseven, flag;
        for( i = 0 ; i < row ; i++ ) {
                fmpz_set(tmp, fmpz_mat_entry(c1, i, 0));
                fmpz_mod_ui(tmp1, tmp, 2);   /* tmp1 = base = tmp % r */
                flag = fmpz_get_si(tmp1);
                fmpz_mul(coeff, tmp, pp);
                fmpz_fdiv_q(tmp2, coeff, qq);
                //fmpz_print(tmp2);
                fmpz_mod_ui(tmp3, tmp2, 2);
                iseven = fmpz_get_si(tmp3);
              //  printf(" %d\n", iseven);
                if(flag != iseven) {
                        if(fmpz_cmp_si(tmp2, 0) > 0) {
                                fmpz_sub_ui(tmp2, tmp2, 1);
                        }
                        else if(fmpz_cmp_si(tmp2, 0) < 0) {
                                fmpz_add_ui(tmp2, tmp2, 1);
                        }
                }
                
                fmpz_set(fmpz_mat_entry(c2, i, 0), tmp2);
        }
        fmpz_clear(coeff);
        fmpz_clear(tmp);
        fmpz_clear(tmp1);
        fmpz_clear(tmp2);
        fmpz_clear(tmp3);
}

void switchkey(fmpz_mat_t c3, fmpz_mat_t mapb, fmpz_mat_t c1, fmpz_t qq, long c1row, long bcol)
{
	fmpz_mat_t bd, bdt;
	long len = fmpz_clog_ui(qq, 2);
	long qrow = c1row * len;
	fmpz_mat_init(bd, qrow, 1);
	bitdecomp(bd, c1, qq, c1row);
	long i, j;
	
	fmpz_poly_mat_init(bdt, 1, qrow);
        for( i = 0 ; i < qrow ; i++ ) {
                fmpz_set(fmpz_mat_entry(bdt, 0, i), fmpz_mat_entry(bd, i, 0));
	}
	fmpz_mat_mul(c3, bdt, mapb);
	for( i = 0 ; i < bcol ; i++ ) {
                fmpz_smod(fmpz_mat_entry(c3, 0, i), qq);
	}
	fmpz_mat_clear(bd);
	fmpz_mat_clear(bdt);
}

void hcrypt_bgv_refresh(fmpz_mat_t c3, fmpz_mat_t c, fmpz_mat_t map, fmpz_t qq, fmpz_t pp, long crow, long mapcol)
{        
        fmpz_mat_t c1;
        long len = fmpz_clog_ui(qq, 2);
	long row = crow * len;
	fmpz_mat_init(c1, row, 1);
        powers(c1, c, qq, crow);
        fmpz_mat_t c2;
        fmpz_mat_init(c2, row, 1);
        scale(c2, c1, qq, pp, row);
        switchkey(c3, map, c2, pp, row, mapcol);
        
	fmpz_mat_clear(c1);
	fmpz_mat_clear(c2);
}
