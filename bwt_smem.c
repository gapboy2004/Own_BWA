#include "bwt.h"
#include <stdio.h>

/* ค้นหา pattern ใน genome โดยไม่ต้อง scan ทั้ง genome */
void backward_search(const FMIndex *fm,
                     const char *pattern,
                     uint64_t m,   /* ความยาว pattern */
                     uint64_t *lo, /* output: range เริ่ม */
                     uint64_t *hi)
{ /* output: range สิ้นสุด */
    *lo = 0;
    *hi = fm->bwt->n - 1;

    for (int64_t i = m - 1; i >= 0; i--)
    {
        uint8_t c = encode_base(pattern[i]);
        *lo = fm->C[c] + fm->Occ[c][*lo];
        *hi = fm->C[c] + fm->Occ[c][*hi + 1] - 1;
        if (*lo > *hi)
        {
            /* pattern ไม่พบใน genome */
            *lo = *hi = UINT64_MAX;
            return;
        }
    }
}

int find_smems(const FMIndex *fm,
               const char *read,
               uint64_t m,
               SMEM *smems)
{

    int n_smems = 0;
    uint64_t i = 0; /* position ปัจจุบันใน read */

    while (i < m)
    {
        uint64_t lo = 0;
        uint64_t hi = fm->bwt->n - 1;
        uint64_t j = i; /* ขยายไปขวาจาก i */
        uint64_t last_lo = lo, last_hi = hi;

        /* ขยาย pattern ไปขวาจนกว่า range จะ empty */
        while (j < m)
        {
            uint8_t c = encode_base(read[j]);

            if (c >= ALPHA_SIZE) // testcase: ถ้าเจอ base ที่ไม่รู้จัก (เช่น 'N') ให้หยุดขยาย
            {          
                j++;   
                break; 
            }
            
            uint64_t nlo = fm->C[c] + fm->Occ[c][lo];
            uint64_t nhi = fm->C[c] + fm->Occ[c][hi + 1] - 1;

            if (nlo > nhi)
                break; /* ขยายต่อไม่ได้ → หยุด */

            last_lo = nlo;
            last_hi = nhi;
            j++;
        }

        /* บันทึก SMEM ถ้า match อย่างน้อย 1 ตัว */
        if (j > i)
        {
            smems[n_smems].qbeg = i;
            smems[n_smems].qend = j;
            smems[n_smems].lo = last_lo;
            smems[n_smems].hi = last_hi;
            n_smems++;

            i = j; /* กระโดดไป position ถัดจาก SMEM */
        }
        else
        {
            i++; /* ตัวนี้ match ไม่ได้เลย ข้ามไป */
        }
    }

    return n_smems;
}

#define MAX_GAP 100
#define GAP_PENALTY 1

int chain_seeds(const uint32_t *SA,
                SMEM *smems,
                int n_smems,
                ChainSeed *seeds)
{

    int n_seeds = 0;

    /* Step 1: expand SMEM → seed แต่ละ hit */
    for (int i = 0; i < n_smems; i++)
    {
        for (uint64_t r = smems[i].lo; r <= smems[i].hi; r++)
        {
            seeds[n_seeds].rpos = SA[r];
            seeds[n_seeds].qbeg = smems[i].qbeg;
            seeds[n_seeds].qend = smems[i].qend;
            seeds[n_seeds].score = (int)(smems[i].qend - smems[i].qbeg);
            seeds[n_seeds].prev = -1;
            n_seeds++;
        }
    }

    /* Step 2: เรียง seed ตาม rpos */
    /* naive bubble sort — Phase 1 */
    for (int i = 0; i < n_seeds - 1; i++)
    {
        for (int j = i + 1; j < n_seeds; j++)
        {
            if (seeds[j].rpos < seeds[i].rpos)
            {
                ChainSeed tmp = seeds[i];
                seeds[i] = seeds[j];
                seeds[j] = tmp;
            }
        }
    }

    /* Step 3: DP หา chain score */
    for (int i = 1; i < n_seeds; i++)
    {
        for (int j = 0; j < i; j++)
        {

            /* เงื่อนไข chain ได้ */
            if (seeds[j].qbeg >= seeds[i].qbeg)
                continue; /* read ต้องไม่สลับ  */
            if (seeds[j].rpos >= seeds[i].rpos)
                continue; /* genome ต้องไม่สลับ */

            uint64_t gap = seeds[i].rpos - seeds[j].rpos;
            if (gap > MAX_GAP)
                continue; /* ไม่ไกลเกินไป      */

            int len_i = (int)(seeds[i].qend - seeds[i].qbeg);
            int penalty = (int)gap * GAP_PENALTY;
            int new_score = seeds[j].score + len_i - penalty;

            if (new_score > seeds[i].score)
            {
                seeds[i].score = new_score;
                seeds[i].prev = j;
            }
        }
    }

    return n_seeds;
}