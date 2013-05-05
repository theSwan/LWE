#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gmp.h>
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "libbgv.h"

char str[100];

int main(int argc, char *args[])
{
	FILE *fp;
        
        if((fp = fopen(args[1], "r")) == NULL)
        {
                printf("file read error\n");
                exit(0);
        }
        
        fmpz_t tmp;
        fmpz_init(tmp);
        fgets(str, 100, fp);
        fmpz_set_str(tmp, str, 10);
        
        long lev, d, i, j;
        lev = fmpz_get_si(tmp);
        fgets(str, 100, fp);
        fmpz_set_str(tmp, str, 10);
        d = fmpz_get_si(tmp);
        
        param_node_t *hp, *pr, *ps;
        hp = param_node_init(hp);
        
        ps = hp;
        for( i = 0 ; i <= lev; i++ ) {
                pr = param_node_init(pr);
                fgets(str, 100, fp);
                fmpz_set_str(pr->q, str, 10);
                fgets(str, 100, fp);
                fmpz_set_str(tmp, str, 10);
                pr->n = fmpz_get_si(tmp);
                fgets(str, 100, fp);
                fmpz_set_str(tmp, str, 10);
                pr->bign = fmpz_get_si(tmp);
                ps->next = pr;
                ps = pr;
        }
        ps->next = NULL;
        
        fclose(fp);
        
        hp = hp->next;
        
        long high, low, ctlev, k, row1, row2, row;
        
        fmpz_set_str(tmp, args[4] ,10);
        high = fmpz_get_si(tmp);
        
        fmpz_set_str(tmp, args[5],10);
        low = fmpz_get_si(tmp);
        ctlev = low - 1;
        
        long l = lev;
        
        while( l > high ) {
                hp = hp->next;
                l--;
        }
        fmpz_mat_t levhigh, levlow;
        row1 = hp->n + 1;
        fmpz_mat_init(levhigh, row1, 1);
        int flag;
        if((fp = fopen(args[2], "r")) != NULL) {
                flag = fmpz_mat_fread(fp, levhigh);
        }
        fclose(fp);
        while( l > low ) {
                hp = hp->next;
                l--;
        }
        row2 = 1 + hp->n;
        fmpz_mat_init(levlow, row2, 1);
        if((fp = fopen(args[3], "r")) != NULL) {
                flag = fmpz_mat_fread(fp, levlow);
        }
        fclose(fp);
        l = lev;
        while( l > high ) {
                hp = hp->next;
                l--;
        }
        l--;
        fmpz_mat_t ctmp, tmpp;
        char s[]="pkb0.txt";
        while( l >= low ) {
                row = (1+hp->n)*(1+hp->n);
                fmpz_mat_init(ctmp, 1 + hp->next->n, 1);
                fmpz_mat_init(tmpp, row, 1);
                fmpz_mat_zero(tmpp);
                for( k = 0 ; k < 1+hp->n; k++) {
                        fmpz_set(fmpz_mat_entry(tmpp, k ,0), fmpz_mat_entry(levhigh, k, 0));
                }
                fmpz_mat_t pkb;
                s[3] = '0'+l-1;
                if((fp = fopen(s, "r")) != NULL) {
                        flag = fmpz_mat_fread(fp, pkb);
                }
                fclose(fp);
                hcrypt_bgv_refresh(ctmp, tmpp, pkb, hp->q, hp->next->q, row, 1 + hp->next->n);
                fmpz_mat_swap(ctmp, levhigh);
                fmpz_mat_clear(ctmp);
                fmpz_mat_clear(tmpp);
                fmpz_mat_clear(pkb);
                hp = hp->next;
                l--;
        }
        fmpz_mat_t c3, ct;
        row = row1 * row2;
        fmpz_mat_init(c3, row, 1);
        
        for( i = 0 ; i < row1 ; i++ ) {
                for( j = 0 ; j < row2 ; j++ ){
                        fmpz_mul(fmpz_mat_entry(c3, j + i * row1, 0), fmpz_mat_entry(levhigh, i, 0), fmpz_mat_entry(levlow, j, 0));
                        fmpz_smod(fmpz_mat_entry(c3, j + i * row1, 0), hp->q);
                }
        }
        
        row1 = 1 + hp->next->n;
        fmpz_mat_t pkb;
        s[3] = '0'+l-1;
        if((fp = fopen(s, "r")) != NULL) {
                flag = fmpz_mat_fread(fp, pkb);
        }
        fclose(fp);
        
        fmpz_mat_init(ct, row1, 1);
        hcrypt_bgv_refresh(ct, c3, pkb, hp->q, hp->next->q, row, row1);
        fmpz_mat_print(ct);
        
        fmpz_mat_clear(levhigh);
        fmpz_mat_clear(levlow);
        fmpz_mat_clear(c3);
        fmpz_mat_clear(ct);
	return 0;
}
