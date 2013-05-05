#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gmp.h>
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"

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

int main()
{
        fmpz_t t, q;
        fmpz_init(t);
        fmpz_init(q);
        fmpz_set_si(q, 3);
        fmpz_set_si(t, 1000);
        fmpz_smod_ui(t, 7);
        fmpz_print(t);
        printf("\n");
        return 0;
}