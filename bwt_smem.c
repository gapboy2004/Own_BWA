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