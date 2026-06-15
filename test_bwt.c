/* test_bwt.c */
#include <stdio.h>
#include <string.h>
#include "bwt.h"

int main(void)
{
    const char *text = "ATCGATCG$";
    uint64_t n = strlen(text) - 1; /* 6 ไม่นับ $ */

    SuffixArray *sa = build_suffix_array(text, n);

    /* พิมพ์ SA แต่ละแถว */
    printf("i\tsa[i]\tsuffix\n");
    printf("─────────────────────────\n");
    for (uint64_t i = 0; i <= n; i++)
        printf("%llu\t%u\t%s\n", (unsigned long long)i, sa->sa[i], text + sa->sa[i]);

    BWT *bwt = build_bwt_from_sa(text, sa);
    /* พิมพ์ BWT */
    printf("\nbwt: ");
    for (uint64_t i = 0; i < bwt->n; i++)
        printf("%c", decode_base(bwt->bwt[i]));
    printf("\n");

    /* พิมพ์ FM Index */
    FMIndex *fm = build_fm_index(bwt);

    // printf("\nFM Index:\n");
    // printf("C:\n");
    // for (int c = 0; c < ALPHA_SIZE; c++)
    //      printf("C[%c] = %llu\n", decode_base(c), (unsigned long long)fm->C[c]);
    // printf("\nOcc:\n");
    // for (int c = 0; c < ALPHA_SIZE; c++) {
    //     printf("Occ[%c]: ", decode_base(c));
    //     for (uint64_t i = 0; i <= bwt->n; i++)
    //         printf("%llu ", (unsigned long long)fm->Occ[c][i]);
    //     printf("\n");
    // }

    /* ทดสอบ LF-mapping */
    printf("\nLF-mapping:\n");
    for (uint64_t i = 0; i < bwt->n; i++)
    {
        uint64_t lf_i = lf_map(fm, i);
        printf("LF(%llu) = %llu\n", (unsigned long long)i, (unsigned long long)lf_i);
    }
    /* ทดสอบ LF-mapping แบบกลับด้าน */
    printf("\nLF-roundtrip:\n");
    lf_roundtrip(fm);

    printf("i\tbwt[i]\tLF(i)\n");
    for (uint64_t i = 0; i < fm->bwt->n; i++)
        printf("%llu\t%c\t%llu\n", i, decode_base(fm->bwt->bwt[i]), lf_map(fm, i));

    uint64_t lo, hi;
    backward_search(fm, "ATCG", 4, &lo, &hi);
    printf("ATCG found at SA rows [%llu, %llu]\n", lo, hi);
    sa_free(sa);
    bwt_free(bwt);
    fm_free(fm);
    return 0;
}
