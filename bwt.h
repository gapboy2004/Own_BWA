#ifndef BWA_BWT_H
#define BWA_BWT_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ── alphabet ──────────────────────────────────────────────
   encode: A=0, C=1, G=2, T=3, $=4  (sentinel always smallest)
   ─────────────────────────────────────────────────────────── */
#define ALPHA_SIZE  5
#define SENTINEL    '$'

/* TODO: implement encode_base()
   input : char c  ('A','C','G','T','$')
   output: uint8_t  0-4                  */
uint8_t encode_base(char c);

/* TODO: implement decode_base()
   input : uint8_t b  0-4
   output: char                          */
char decode_base(uint8_t b);


/* ── suffix array ──────────────────────────────────────────
   sa[i] = starting position of i-th smallest suffix
   ─────────────────────────────────────────────────────────── */
typedef struct {
    uint32_t *sa;   /* suffix array, length n+1          */
    uint64_t  n;    /* length of original string (no $)  */
} SuffixArray;

/* TODO: build_suffix_array()
   naive O(n^2 log n) is fine for Phase 1
   ─────────────────────────────────────────────────────────── */
SuffixArray *build_suffix_array(const char *text, uint64_t n);
void         sa_free(SuffixArray *sa);


/* ── BWT ───────────────────────────────────────────────────
   bwt[i] = text[ sa[i] - 1 ]  (wrap: if sa[i]==0 → '$')
   ─────────────────────────────────────────────────────────── */
typedef struct {
    uint8_t  *bwt;      /* BWT string, length n+1, encoded   */
    uint64_t  n;        /* length (including sentinel)        */
    uint64_t  eof_pos;  /* position of '$' in BWT            */
} BWT;

/* TODO: build_bwt_from_sa()
   given text + suffix array → BWT string
   ─────────────────────────────────────────────────────────── */
BWT *build_bwt_from_sa(const char *text, const SuffixArray *sa);
void bwt_free(BWT *bwt);


/* ── FM-index ──────────────────────────────────────────────
   C[c]    = number of bases in BWT that are < c
   Occ[c][i] = occurrences of c in bwt[0..i-1]
   ─────────────────────────────────────────────────────────── */
typedef struct {
    BWT      *bwt;
    uint64_t  C[ALPHA_SIZE];        /* cumulative counts      */
    uint64_t **Occ;                 /* Occ[alpha][pos]        */
} FMIndex;

/* TODO: build_fm_index()
   compute C[] and Occ[][] from BWT
   ─────────────────────────────────────────────────────────── */
FMIndex *build_fm_index(BWT *bwt);
void     fm_free(FMIndex *fm);

/* TODO: lf_map()
   LF(i) = C[bwt[i]] + Occ[bwt[i]][i]
   ─────────────────────────────────────────────────────────── */
uint64_t lf_map(const FMIndex *fm, uint64_t i);

#endif /* BWA_BWT_H */

typedef struct {
    uint64_t qbeg;   /* ตำแหน่งเริ่มใน read  */
    uint64_t qend;   /* ตำแหน่งสิ้นสุดใน read */
    uint64_t lo;     /* SA range เริ่ม        */
    uint64_t hi;     /* SA range สิ้นสุด      */
} SMEM;

typedef struct {
    uint64_t rpos;   /* ตำแหน่งใน genome (จาก SA) */
    uint64_t qbeg;   /* ตำแหน่งเริ่มใน read       */
    uint64_t qend;   /* ตำแหน่งสิ้นสุดใน read     */
    int      score;  /* chain score                */
    int      prev;   /* index ของ SMEM ก่อนหน้า   */
} ChainSeed;

typedef struct {
    int  score;
    int  qbeg, qend;
    int  rbeg, rend;
    char cigar[256];
} SWResult;

int find_smems(const FMIndex *fm, 
               const char *read, 
               uint64_t m, 
               SMEM *smems);

int chain_seeds(const uint32_t *SA, 
                SMEM *smems, 
                int n_smems, 
                ChainSeed *seeds);

/* ── SW Extension ────────────────────────────────────────── */
// (คงเดิมไว้ที่นี่ได้ แต่แนะนำให้แยกไฟล์ในอนาคต)
SWResult sw_extend(const char *read, int qlen, 
                   const char *ref, int rlen);

/* ── Helper functions ที่ขาดใน Header ────────────────────── */
void lf_roundtrip(const FMIndex *fm);
void backward_search(const FMIndex *fm, const char *pattern, uint64_t m, uint64_t *lo, uint64_t *hi);