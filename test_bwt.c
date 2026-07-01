#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h> // ใช้สำหรับ PRIu64 เพื่อความปลอดภัยของ Format
#include "bwt.h"

int main(void)
{
    const char *text = "GATTC$";
    uint64_t n = strlen(text) - 1;

    /* 1. Build Index ตามลำดับ */
    SuffixArray *sa = build_suffix_array(text, n);
    BWT *bwt = build_bwt_from_sa(text, sa);
    FMIndex *fm = build_fm_index(bwt);

    /* 2. พิมพ์ Suffix Array */
    printf("i\tSA[i]\tSuffix\n");
    printf("─────────────────────────\n");
    for (uint64_t i = 0; i <= n; i++)
        printf("%" PRIu64 "\t%u\t%s\n", i, sa->sa[i], text + sa->sa[i]);

    /* 3. ทดสอบการค้นหา (SMEMs) */
    const char *read = "GATC";
    uint64_t m = strlen(read);
    
    // ใช้ขนาดที่เหมาะสม หรือจัดการแบบ Dynamic ถ้า Read ยาวขึ้น
    SMEM smems[64];
    ChainSeed seeds[256];

    int n_smems = find_smems(fm, read, m, smems);
    int n_seeds = chain_seeds(sa->sa, smems, n_smems, seeds);

    if (n_seeds > 0) {
        /* หา seed ที่มี score สูงสุด */
        int best = 0;
        for (int i = 1; i < n_seeds; i++)
            if (seeds[i].score > seeds[best].score)
                best = i;

        /* Traceback หาจุดเริ่ม Alignment ที่แท้จริง */
        int idx = best;
        int first = best;
        while (seeds[idx].prev != -1) {
            idx = seeds[idx].prev;
            first = idx;
        }

        uint64_t alignment_start = (uint64_t)seeds[first].rpos - seeds[first].qbeg;

        printf("\nBest chain (score=%d):\n", seeds[best].score);
        printf("Alignment starts at genome pos: %" PRIu64 "\n", alignment_start);
        
        /* แสดงรายการ Seed ใน Chain */
        idx = best;
        while (idx != -1) {
            printf("  read[%" PRIu64 ",%" PRIu64 ") → ref pos %" PRIu64 "\n",
                   seeds[idx].qbeg, seeds[idx].qend, seeds[idx].rpos);
            idx = seeds[idx].prev;
        }
    } else {
        printf("\nNo seeds found.\n");
    }

    /* 4. Smith-Waterman Extension */
    printf("\n--- Smith-Waterman ---\n");
    const char *ref = "GATTC";
    SWResult r = sw_extend(read, (int)m, ref, (int)strlen(ref));

    printf("score : %d\n", r.score);
    printf("query : [%d,%d)\n", r.qbeg, r.qend);
    printf("ref   : [%d,%d)\n", r.rbeg, r.rend);
    printf("cigar : %s\n", r.cigar);

    /* 5. Cleanup ตามลำดับการสร้าง (ย้อนหลัง) */
    fm_free(fm);
    bwt_free(bwt);
    sa_free(sa);

    return 0;
}