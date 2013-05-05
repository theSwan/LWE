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
int main(int argc, char *args[])/*setup.txt*/
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
        
        param_node_t *ph, *pr, *ps, *param;
        ph = param_node_init(ph);

        ps = ph;
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
        param = ph->next;
        pr = param->next;
        
        fclose(fp);

        fmpz_mat_t ss, sr;
        fmpz_mat_init(ss, 1 + param->n, 1);
        fmpz_mat_init(sr, 1 + pr->n, 1);
        char name[]="sk0.txt";
        char aname[]="pka0.txt";
        char bname[]="pkb0.txt";
        name[2]='0'+lev;
        aname[3]=name[2];
        bname[3]=aname[3];
        int flag;
        if((fp = fopen(name, "r")) != NULL)
        {
                flag = fmpz_mat_fread(fp, ss);
        }
        fclose(fp);
        name[2]='0'+lev-1;
        if((fp = fopen(name, "r")) != NULL)
        {
                flag = fmpz_mat_fread(fp, sr);
        }
        fclose(fp);
        fmpz_mat_t pka,pkb;
        fmpz_mat_init(pkb, 1, 1);
        fmpz_mat_init(pka, param->bign, 1 + (param->n));
        e_pkeygen(pka, param, ss);
        if((fp = fopen(aname, "w")) != NULL) {
                flag = fmpz_mat_fprint(fp, pka);
        }
        fclose(fp);
        if((fp = fopen(bname, "w")) != NULL) {
                flag = fmpz_mat_fprint(fp, pkb);
        }
        fclose(fp);
        fmpz_mat_clear(pka);
        fmpz_mat_clear(pkb);
                
        fmpz_mat_t s1, tensor;
	long row1, row2, row3, row4, row5, len, llog;
        
        for( i = lev-1 ; i >= 0 ; i-- ) {
                aname[3]='0'+i;
                bname[3]=aname[3];
                llog = fmpz_clog_ui(pr->q, 2);
		fmpz_mat_init(pka, pr->bign, 1 + (pr->n));

                e_pkeygen(pka, pr, sr);
                
                if((fp = fopen(aname, "w")) != NULL) {
                        flag = fmpz_mat_fprint(fp, pka);
                }
                fclose(fp);
                
                row1 = param->n + 1;
        	row2 = row1 * row1;
        	fmpz_mat_init(tensor, row2, 1);
                vec_tensor(tensor, ss, param->q, row1);
                
                len = fmpz_clog_ui(param->q, 2);
		row3 = row2 * len;
                fmpz_mat_init(s1, row3, 1);
                bitdecomp(s1, tensor, param->q, row2);
                
                row4 = row3 * llog;
		row5 = 1 + pr->n;
		fmpz_mat_init(pkb, row4, row5);
                
                switchkeygen(pkb, s1, sr, pr->q, row3, row5);
                
                fmpz_mat_clear(s1);
                fmpz_mat_clear(tensor);
                fmpz_mat_clear(ss);
         
                if((fp = fopen(bname, "w")) != NULL) {
                        flag = fmpz_mat_fprint(fp, pkb);
                }
                fclose(fp);
                
                if(i != 0) {
                        pr = pr->next;
                        param = param->next;
                        name[2]='0'+i-1;
                        fmpz_mat_init(ss, 1 + pr->n, 1);
                        if((fp = fopen(name, "r")) != NULL) {
                                flag = fmpz_mat_fread(fp, ss);
                        }
                        fclose(fp);
                        fmpz_mat_swap(ss, sr);
                }
                fmpz_mat_clear(pka);
                fmpz_mat_clear(pkb);
        }
        
        fmpz_mat_clear(sr);
        return 0;
}