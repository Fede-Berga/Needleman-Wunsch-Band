#include "sw.h"
#define BUZZ_SIZE 1024
#define DEBUG

void read_text(char* s1, char* s2) {
    printf("Inserisci la stringa s1 : ");
    scanf("%s", s1);

    printf("Inserisci la stringa s2 : ");
    scanf("%s", s2);

    #ifdef DEBUG
        printf("\ns1 : %s\n", s1);
        printf("\ns2 : %s\n", s2);
    #endif
}

int score(char c1, char c2){
    if (c1 == c2 == '*')
    {
        #ifdef DEBUG
            printf("c1 == c2 == '*'\n");
        #endif
        return 1;
    }

    if (c1 == '*' || c2 == '*')
    {   
        #ifdef DEBUG
            printf("c1 == '*' || c2 == '*'\n");
        #endif
        return -4;
    }
    
    if (c1 == c2)
    {
        #ifdef DEBUG
            printf("c1 == c2\n");
        #endif
        for (uint64_t i = 0; i < SCORING_LEN; i++)
        {
            if (scoring_lookup[i] == c1)
            {
                return scoring_matrix[i][i];
            }
        }  
    }

    #ifdef DEBUG
            printf("c1 != c2\n");
    #endif
    
    unsigned int i_c1, i_c2;

    for (uint64_t i = 0; i < SCORING_LEN; i++)
    {
        if (scoring_lookup[i] == c1)
        {
            i_c1 = i;
        }
    } 

    for (uint64_t i = 0; i < SCORING_LEN; i++)
    {
        if (scoring_lookup[i] == c2)
        {
            i_c2 = i;
        }
    }

    printf("%d, %d\n", i_c1, i_c2);

    return scoring_matrix[i_c1][i_c2];
}

range_s range(band_s band, uint64_t y) {
    range_s r = {0, 0};

    r.sx = y - band.extra;

    if (y < band.extra) 
        r.sx = 0;

    r.dx = y + (band.base - 1) + band.extra;

    if (r.dx > band.l2) 
        r.dx = band.l2;
    
    return r;
}

uint64_t conv(band_s band, uint64_t x, uint64_t y) {
    assert(y <= band.l1);

    assert(x >= y || y - x <= band.extra);
    assert(x <= y || x - y <= band.base - 1 + band.extra);

    uint64_t width = band.base + 2 * band.extra;
    
    return (width * y) + (x - y + band.extra);
}

static void band_align(uint64_t base, uint64_t extra, char* s1, char* s2) {
    uint64_t l1 = strlen(s1);
    uint64_t l2 = strlen(s2);
    uint64_t band = base + 2 * extra;

    cell_s* m = (cell_s*)malloc(band * (l1 + 1) * sizeof(cell_s));
    assert(m != NULL && "Cannot allocate matrix\n");

    /* boundary conditions */
    m[0].cell = 0;
    for (uint64_t x = 1; x < band; x++)
        m[x].cell = 0;
    band_s b = {
        .l1 = l1,
        .l2 = l2,
        .base = base,
        .extra = extra
    };

    for (uint64_t y = 1; y <= l1; y++) {
        range_s r = range(b, y);
        {
            cell_s* tp = m + conv(b, r.sx, y);
            tp->cell = 0;
            tp->prev_x = 0;
            tp->prev_y = 0;
            tp = m + conv(b, r.dx, y);
            tp->cell = 0;
            tp->prev_x = 0;
            tp->prev_y = 0;
        }
        for (uint64_t x = r.sx + 1; x <= r.dx; x++) {
            #ifdef DEBUG
                printf("StepB: %ld, %ld\n", x, y);
            #endif
            cell_s* tp = m + conv(b, x, y);
            tp->cell = m[conv(b, x - 1, y - 1)].cell;
            tp->prev_x = x -1;
            tp->prev_y = y -1;
            if (s1[y - 1] == s2[x - 1]) ++(tp->cell);
            if (m[conv(b, x - 1, y)].cell > tp->cell) {
                tp->cell = m[conv(b, x - 1, y)].cell;
                tp->prev_x = x - 1;
                tp->prev_y = y;
            }
            if (x < r.dx && m[conv(b, x, y - 1)].cell > tp->cell) {
                tp->cell = m[conv(b, x, y - 1)].cell;
                tp->prev_x = x;
                tp->prev_y = y - 1;
            }
            #ifdef DEBUG
                printf("StepB result: %ld, %ld, %ld\n", m[conv(b, x, y)].cell,
                m[conv(b, x, y)].prev_x, m[conv(b, x, y)].prev_y);
            #endif
        }
    }
}



int main(int argc, char const *argv[])
{
    char s1[BUZZ_SIZE];
    char s2[BUZZ_SIZE];

    read_text(s1, s2);

    to_upper_case(s1);
    to_upper_case(s2);
    
    for (size_t i = 0; i < strlen(s1); i++)
    {
        printf("score : %d\n", score(s1[i], s2[i]));
    }
    

    return 0;
}
