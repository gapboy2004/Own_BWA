#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include "bwt.h"

#define MATCH     2
#define MISMATCH -3
#define GAP_OPEN -4
#define GAP_EXT  -1

#define MAX(a,b) ((a)>(b)?(a):(b))

#define STOP 0
#define FROM_DIAG 1
#define FROM_LEFT 2
#define FROM_UP   3

SWResult sw_extend(const char *read, int qlen,
                   const char *ref,  int rlen)
{
    /* --- allocate --- */
    int **H  = malloc((qlen+1) * sizeof(int*));
    int **E  = malloc((qlen+1) * sizeof(int*));
    int **F  = malloc((qlen+1) * sizeof(int*));
    int **tb = malloc((qlen+1) * sizeof(int*));
    for(int i = 0; i <= qlen; i++){
        H[i]  = calloc(rlen+1, sizeof(int));
        E[i]  = calloc(rlen+1, sizeof(int));
        F[i]  = calloc(rlen+1, sizeof(int));
        tb[i] = calloc(rlen+1, sizeof(int));
    }

    /* --- init row 0 และ col 0 (semi-global) --- */
    // for(int j = 1; j <= rlen; j++){
    //     H[0][j]  = GAP_OPEN + (j-1) * GAP_EXT;
    //     E[0][j]  = H[0][j];
    //     tb[0][j] = FROM_LEFT;
    // }
    // for(int i = 1; i <= qlen; i++){
    //     H[i][0]  = GAP_OPEN + (i-1) * GAP_EXT;
    //     F[i][0]  = H[i][0];
    //     tb[i][0] = FROM_UP;
    // }

    for(int j = 1; j <= rlen; j++){
        H[0][j]  = 0;
        F[0][j] = -99;
        E[0][j]  = -99;
        tb[0][j] = -99;
    }
    for(int i = 1; i <= qlen; i++){
        H[i][0]  = 0;
        F[i][0]  = -99;
        E[i][0]  = -99;
        tb[i][0] = -99;
    }

    int best_score = INT_MIN, best_i = qlen, best_j = rlen;

    /* --- fill --- */
    for(int i = 1; i <= qlen; i++){
        for(int j = 1; j <= rlen; j++){

            E[i][j] = MAX(H[i][j-1] + GAP_OPEN,
                          E[i][j-1] + GAP_EXT);

            F[i][j] = MAX(H[i-1][j] + GAP_OPEN,
                          F[i-1][j] + GAP_EXT);

            int sub  = (read[i-1] == ref[j-1]) ? MATCH : MISMATCH;
            int diag = H[i-1][j-1] + sub;

            H[i][j] = MAX(MAX(MAX(diag, E[i][j]), F[i][j]),0);  /* ไม่มี 0 */

            if (H[i][j] == 0)
                tb[i][j] = STOP;
            else if (H[i][j] == diag)
                tb[i][j] = FROM_DIAG;
            else if (H[i][j] == E[i][j])
                tb[i][j] = FROM_LEFT;
            else
                tb[i][j] = FROM_UP;

            if(H[i][j] > best_score){
                best_score = H[i][j];
                best_i = i; best_j = j;
            }
        }
    }
    printf("H matrix:\n");
    for(int i = 0; i <= qlen; i++){
        for(int j = 0; j <= rlen; j++){

            printf("%d \t", H[i][j]);
        }
        printf("\n");
    }

    printf("E matrix:\n");
    for(int i = 0; i <= qlen; i++){
        for(int j = 0; j <= rlen; j++){

            printf("%d \t", E[i][j]);
        }
        printf("\n");
    }
    printf("F matrix:\n");
    for(int i = 0; i <= qlen; i++){
        for(int j = 0; j <= rlen; j++){

            printf("%d \t", F[i][j]);
        }
        printf("\n");
    }
    printf("tb matrix:\n");
    for(int i = 0; i <= qlen; i++){
        for(int j = 0; j <= rlen; j++){

            printf("%d \t", tb[i][j]);
        }
        printf("\n");
    }

    /* --- traceback --- */
    char raw_cigar[512];
    int  clen = 0;
    int  i = best_i, j = best_j;

    while(i > 0 && j > 0){
        if(tb[i][j] == FROM_DIAG){
            raw_cigar[clen++] = (read[i-1]==ref[j-1]) ? 'M' : 'X';
            i--; j--;
        } else if(tb[i][j] == FROM_LEFT){
            raw_cigar[clen++] = 'D';
            j--;
        } else {
            raw_cigar[clen++] = 'I';
            i--;
        }
    }

    /* --- build result --- */
    SWResult res;
    res.score = best_score;
    res.qbeg  = i;   res.qend = best_i;
    res.rbeg  = j;   res.rend = best_j;

    /* compress cigar */
    int ci = 0;
    for(int k = clen-1; k >= 0;){
        char op  = raw_cigar[k];
        int  cnt = 0;
        while(k >= 0 && raw_cigar[k] == op){ cnt++; k--; }
        ci += sprintf(res.cigar + ci, "%d%c", cnt, op);
    }
    res.cigar[ci] = '\0';

    /* --- free --- */
    for(int k = 0; k <= qlen; k++){
        free(H[k]); free(E[k]); free(F[k]); free(tb[k]);
    }
    free(H); free(E); free(F); free(tb);

    printf("Best score = %d\n", best_score);
printf("Best cell  = (%d,%d)\n", best_i, best_j);
    return res;
}