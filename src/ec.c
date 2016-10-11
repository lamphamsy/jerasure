/* *
 * Copyright (c) 2016, Scality
 * All rights reserved.
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <gf_complete.h>
#include "galois.h"
#include "jerasure.h"
#include "cauchy.h"
#include "ec.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

int *gfDensity;     // binary density of GF elements
int *gfList;        // list of GF elements with an order of increasing binary
                    //  density
int *gfListInv;     // inversed elements of gfList
int *pOnes;
int gfCard;
int alpha;

void show(int *mat, int m, int k) {
    int i, j;
    // show pcm
    for (i = 0; i < m; i++) {
        for (j = 0; j < k; j++) {
            printf("%2d ", mat[i*k + j]);
        }
        printf("\n");
    }
}

int compDensity(const void * elem1, const void * elem2) {
    return ( gfDensity[*(int*)elem1] - gfDensity[*(int*)elem2] );
}

int compOnes(const void * elem1, const void * elem2) {
    return ( pOnes[*(int*)elem1] - pOnes[*(int*)elem2] );
}

void init(int w) {
    int i;
    int n_ones;

    gfCard = 1 << w;

    // alloc mem
    gfDensity = talloc(int, gfCard);
    gfList = talloc(int, gfCard - 1);
    gfListInv = talloc(int, gfCard - 1);

    // gen gfDensity
    for (i = 0; i < gfCard - 1; i++) {
        n_ones = cauchy_n_ones(i + 1, w);
        gfDensity[i + 1] = n_ones;
        gfList[i] = i + 1;
    }
    // sort gfList
    qsort(gfList, gfCard - 1, sizeof(* gfList), compDensity);

    // gen gfListInv
    for (i = 0; i < gfCard - 1; i++) {
        gfListInv[i] = galois_single_divide(1, gfList[i], w);
    }
}

int *ec_pcm(int k, int m, int w) {
    int *pcm = NULL;
    int *res = NULL;
    int *x = NULL;
    int *y = NULL;
    int *bList = NULL;
    int *rList = NULL;
    int i, j, u, v;
    int a, b;
    int idx;

    if (w < 30 && (1 << w) < m) return NULL;
    if (w < 30 && (1 << w) < k) return NULL;

    init(w);
    // get alpha
    alpha = gfListInv[k];

    // printf("SHOW gfList\n");
    // show(gfList, 1, gfCard - 1);
    // printf("SHOW gfListInv\n");
    // show(gfListInv, 1, gfCard - 1);
    // printf("ALPHA %d\n", alpha);

    u = gfCard - k;
    pcm = talloc(int, u*k);
    if (pcm == NULL) { goto end; }
    res = talloc(int, m*k);
    if (res == NULL) { goto end; }
    x = talloc(int, u);
    if (x == NULL) { goto end; }
    y = talloc(int, k);
    if (y == NULL) { goto end; }
    bList = talloc(int, gfCard);
    if (bList == NULL) { goto end; }
    pOnes = talloc(int, u);
    if (pOnes == NULL) { goto end; }
    rList = talloc(int, u);
    if (rList == NULL) { goto end; }

    // init bList as all elements;
    for (i = 0; i < gfCard; i++) {
        bList[i] = 1;
    }

    // init pOnes as zeros
    for (i = 0; i < u; i++) {
        pOnes[i] = 0;
    }
    // init rList
    for (i = 0; i < u; i++) {
        rList[i] = i;
    }

    // init x[0] & x[1]
    x[0] = 2;
    x[1] = 1;
    bList[x[0]] = 0;
    bList[x[1]] = 0;
    // gen y
    for (i = 0; i < k; i++) {
        a = x[0] ^ galois_single_multiply(x[1],
            galois_single_multiply(alpha, gfList[i], w), w);
        b = 1 ^ galois_single_multiply(alpha, gfList[i], w);
        y[i] = galois_single_divide(a, b, w);
        bList[y[i]] = 0;
    }
    // last elements of x are other elements
    idx = 2;
    for (i = 0; i < gfCard; i++) {
        if (bList[i] == 1) {
            x[idx] = i;
            idx++;
        }
    }
    // printf("SHOW X\n");
    // show(x, 1, u);
    // printf("SHOW Y\n");
    // show(y, 1, k);

    // gen pcm
    for (i = 0; i < u; i++) {
        for (j = 0; j < k; j++) {
            pcm[i * k + j] = galois_single_divide(1, x[i] ^ y[j], w);
        }
    }
    // printf("SHOW PCM\n");
    // show(pcm, u, k);

    // normalize 1st row by multiplying elements in jth column by (x[0] + y[j])
    for (j = 0; j < k; j++) {
        a = x[0] ^ y[j];
        for (i = 0; i < u; i++) {
            pcm[i * k + j] = galois_single_multiply(a, pcm[i * k + j], w);
        }
        pOnes[0] += cauchy_n_ones(pcm[j], w);
    }

    // normalize 2nd row by dividing its elements by alpha
    for (j = k; j < 2 * k; j++) {
        pcm[j] = galois_single_divide(pcm[j], alpha, w);
        pOnes[1] += cauchy_n_ones(pcm[j], w);
    }

    // printf("OPT 2 ROWs\n");
    // show(pcm, u, k);

    // optimize next rows
    for (i = 2; i < u; i++) {
        // calculate n_ones of its elements
        a = 0;
        for (j = 0; j < k; j++) {
            a += cauchy_n_ones(pcm[i*k + j], w);
        }

        // optimize
        idx = 0;
        for (v = 2; v < gfCard; v++) {
            b = 0;
            for (j = 0; j < k; j++) {
                b += cauchy_n_ones(
                        galois_single_multiply(v, pcm[i*k + j], w), w);
            }
            if (b < a) {
                a = b;
                idx = v;
            }
        }
        // multiply idx to ith row
        if (idx > 0) {
            for (j = 0; j < k; j++) {
                pcm[i*k + j] = galois_single_multiply(pcm[i*k + j], idx, w);
            }
        }
        for (j = 0; j < k; j++) {
            pOnes[i] += cauchy_n_ones(pcm[i*k + j], w);
        }
    }

    // printf("OPT PCM\n");
    // show(pcm, u, k);
    //
    // printf("pOnes\n");
    // show(pOnes, 1, u);
    //
    // printf("rList\n");
    // show(rList, 1, u);

    // sort rList according to pOnes
    qsort(rList, u, sizeof(*rList), compOnes);

    // printf("rList sorted\n");
    // show(rList, 1, u);

    // result matrix is formed from m rows whose indices are first m elements
    //  from pOnes
    for (i = 0; i < m; i++) {
        for (j = 0; j < k; j++) {
            res[i*k + j] = pcm[rList[i]*k + j];
        }
    }

    // printf("RESULT\n");
    // show(res, m, k);

end:
    if (pcm) {
        free(pcm);
    }
    if (x) {
        free(x);
    }
    if (y) {
        free(y);
    }
    if (bList) {
        free(bList);
    }
    if (pOnes) {
        free(pOnes);
    }
    if (rList) {
        free(rList);
    }
    if (gfDensity) {
        free(gfDensity);
    }
    if (gfList) {
        free(gfList);
    }
    if (gfListInv) {
        free(gfListInv);
    }

    return res;
}
