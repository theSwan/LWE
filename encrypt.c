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

int main(int argc, char *args[])/* setup.txt, ms, pk.txt{row, col, pka, row, col, pkb} */
{
        FILE *fp;
        
        if((fp = fopen(args[1], "r")) == NULL)
        {
                printf("setup file read error\n");
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
        
        param_node_t *param;
        param = param_node_init(param);
        fgets(str, 100, fp);
        fmpz_set_str(param->q, str, 10);
        fgets(str, 100, fp);
        fmpz_set_str(tmp, str, 10);
        param->n = fmpz_get_si(tmp);
        fgets(str, 100, fp);
        fmpz_set_str(tmp, str, 10);
        param->bign = fmpz_get_si(tmp);
        param->next = NULL;

        fclose(fp);
        
        fmpz_t ms;
        fmpz_init(ms);
        fmpz_set_str(ms, args[2], 10);
        
        char name[]="pka0.txt";
        name[3]='0'+lev;
        
        int flag;
        
        fmpz_mat_t pka;
        fmpz_mat_init(pka, param->bign, 1 + (param->n) );
        if((fp = fopen(name, "r")) != NULL)
        {
                flag = fmpz_mat_fread(fp, pka);
        }
        fclose(fp);
        
        printf("%ld\n", lev);
        
        fmpz_mat_t ct;
        
        fmpz_mat_init(ct, 1 + param->n, 1);
        
        e_encrypt(ct, param, pka, ms);
        fmpz_mat_print(ct);
        
        fmpz_mat_clear(pka);
        fmpz_mat_clear(ct);
        fmpz_clear(ms);
        free(param);
        return 0;
}
/*output: ct->lev, row, col, ct->text */