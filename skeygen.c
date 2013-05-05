#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gmp.h>
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "libbgv.h"

int main(int argc, char *args[]) /*setup.txt:lev,d, {q, n, bign} */
{
        FILE *fp;
        
        if((fp = fopen(args[1], "r")) == NULL)
        {
                printf("file read error\n");
                exit(0);
        }
        
        char str[100];
        
        fmpz_t tmp;
        fmpz_init(tmp);
        
        fgets(str, 100, fp);
        fmpz_set_str(tmp, str, 10);
        
        long lev, d, i, j, row;
        lev = fmpz_get_si(tmp);
        fgets(str, 100, fp);
        fmpz_set_str(tmp, str, 10);
        d = fmpz_get_si(tmp);
        FILE *skf;
        
        char name[]="sk1.txt";
        for( i = lev ; i >= 0; i-- ) {
                param_node_t *h;
                fmpz_mat_t sk;
                h = param_node_init(h);
                fgets(str, 100, fp);
                fmpz_set_str(h->q, str, 10);
                fgets(str, 100, fp);
                fmpz_set_str(tmp, str, 10);
                h->n = fmpz_get_si(tmp);
                fgets(str, 100, fp);
                fmpz_set_str(tmp, str, 10);
                h->bign = fmpz_get_si(tmp);
                row = 1 + h->n;
                fmpz_mat_init(sk, row, 1);
                e_skeygen(sk, h);
                name[2] = '0' + i;
                if((skf = fopen(name, "w")) == NULL)
                {
                        printf("file open error\n");
                        exit(0);
                }
                int flag = fmpz_mat_fprint(skf, sk);
                if(flag < 0) {
                        printf("file write error\n");
                        exit(0);
                }
                fclose(skf);
                fmpz_mat_clear(sk);
        }
        fclose(fp);
        fmpz_clear(tmp);
        return 0;
}