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

int main(int argc, char *args[])/* setup.txt, ct.txt{ct->lev, ct->row, ct->col, ct}, lev*/
{
        long lev;
        fmpz_t tmp;
        fmpz_init(tmp);
        fmpz_set_str(tmp, args[3], 10);
        lev = fmpz_get_si(tmp);
        
        FILE *fp;
        
        if((fp = fopen(args[1], "r")) == NULL)
        {
                printf("file read error\n");
                exit(0);
        }
        
        fgets(str, 100, fp);
        fmpz_set_str(tmp, str, 10);
        
        long level, i, j;
        level = fmpz_get_si(tmp);
        fgets(str, 100, fp);
        
        param_node_t *pr;
        pr = param_node_init(pr);
        
        for( i = level ; i >= lev; i-- ) {
                fgets(str, 100, fp);
                fmpz_set_str(pr->q, str, 10);
                fgets(str, 100, fp);
                fmpz_set_str(tmp, str, 10);
                pr->n = fmpz_get_si(tmp);
                fgets(str, 100, fp);
                fmpz_set_str(tmp, str, 10);
                pr->bign = fmpz_get_si(tmp);
        }
        pr->next = NULL;
        fclose(fp);
        
        fmpz_mat_t ct;
        fmpz_mat_init(ct, 1 + pr->n, 1);
        int flag;
        if((fp = fopen(args[2], "r")) != NULL)
        {
                flag = fmpz_mat_fread(fp, ct);
        }
        fclose(fp);
        
        char name[]="sk0.txt";
        name[2]='0'+lev;
        
        fmpz_mat_t sk;
        fmpz_mat_init(sk, 1 + pr->n, 1);
        if((fp = fopen(name, "r")) != NULL)
        {
                flag = fmpz_mat_fread(fp, sk);
        }
        fclose(fp);
        
        fmpz_t ms;
        fmpz_init(ms);
        e_decrypt(ms, pr, sk, ct);
        fmpz_print(ms);
        printf("\n");
        
        fmpz_mat_clear(sk);
        fmpz_mat_clear(ct);
        fmpz_clear(ms);
        free(pr);
        return 0;
        
}