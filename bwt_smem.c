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
               SMEM *smems) {

    int n_smems = 0;
    uint64_t i  = 0;   /* position ปัจจุบันใน read */

    while (i < m) {
        uint64_t lo = 0;
        uint64_t hi = fm->bwt->n - 1;
        uint64_t j  = i;   /* ขยายไปขวาจาก i */
        uint64_t last_lo = lo, last_hi = hi;

        /* ขยาย pattern ไปขวาจนกว่า range จะ empty */
        while (j < m) {
            uint8_t c    = encode_base(read[j]);
            uint64_t nlo = fm->C[c] + fm->Occ[c][lo];
            uint64_t nhi = fm->C[c] + fm->Occ[c][hi + 1] - 1;

            if (nlo > nhi) break;   /* ขยายต่อไม่ได้ → หยุด */

            last_lo = nlo;
            last_hi = nhi;
            j++;
        }

        /* บันทึก SMEM ถ้า match อย่างน้อย 1 ตัว */
        if (j > i) {
            smems[n_smems].qbeg = i;
            smems[n_smems].qend = j;
            smems[n_smems].lo   = last_lo;
            smems[n_smems].hi   = last_hi;
            n_smems++;

            i = j;   /* กระโดดไป position ถัดจาก SMEM */
        } else {
            i++;     /* ตัวนี้ match ไม่ได้เลย ข้ามไป */
        }
    }

    return n_smems;
}