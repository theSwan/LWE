#ifndef LIBBGV_H
#define LIBBGV_H
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gmp.h>
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"

#ifdef __cplusplus
extern "C" {
#endif

        typedef struct param_node_t {
                fmpz_t q;
                long n;
                long bign;
                struct param_node_t *next;
        } param_node_t;

        static const double pi = 3.1415926;
        static const double dvn = 1.0;
        static const long bigb = 20;

        double bgv_get_dvn();
        long bgv_get_bigb();
        void hcrypt_random(mpz_t r);
        void gen_q(fmpz_t q, long len);
        void guassian_vec(fmpz_mat_t mat, long len);
        void unif_mat(fmpz_mat_t mat, fmpz_t q, long row, long col);
        void fmpz_smod(fmpz_t num, fmpz_t q);
        void fmpz_smod_ui(fmpz_t num, long mod);
                
        param_node_t *param_node_init(param_node_t *pnt);
        void powers(fmpz_mat_t po, fmpz_mat_t x, fmpz_t qq, long xrow);
        void bitdecomp(fmpz_mat_t bits, fmpz_mat_t x, fmpz_t qq, long xrow);
        void vec_tensor(fmpz_mat_t tensor, fmpz_mat_t x, fmpz_t qq, long xrow);
        void switchkeygen(fmpz_mat_t mapb, fmpz_mat_t s1, fmpz_mat_t s2, fmpz_t qq, long n1, long n2);
        void scale(fmpz_mat_t c2, fmpz_mat_t c1, fmpz_t qq, fmpz_t pp, long row);
        void switchkey(fmpz_mat_t c3, fmpz_mat_t mapb, fmpz_mat_t c1, fmpz_t qq, long c1row, long bcol);
        void hcrypt_bgv_refresh(fmpz_mat_t c3, fmpz_mat_t c, fmpz_mat_t map, fmpz_t qq, fmpz_t pp, long crow, long mapcol);
        param_node_t *e_setup(long miu, long lamda, long b, param_node_t *param);
        void e_skeygen(fmpz_mat_t sk, param_node_t *param);
        void e_pkeygen(fmpz_mat_t pk, param_node_t *param, fmpz_mat_t sk);
        void e_encrypt(fmpz_mat_t ct, param_node_t *param, fmpz_mat_t pk, fmpz_t ms);
        void e_decrypt(fmpz_t ms, param_node_t *param, fmpz_mat_t sk, fmpz_mat_t ct);
        
#ifdef __cplusplus
}
#endif
#endif
