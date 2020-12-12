#include "sw.h"
#define BUZZ_SIZE 1024
#define DEBUG

void read_text(char* s1, char* s2) {
    printf("Inserisci la stringa s1 : ");
    scanf("%s", s1);

    printf("Inserisci la stringa s2 : ");
    scanf("%s", s2);

    #ifdef DEBUG
        printf("s1 : %s\n", s1);
        printf("s2 : %s\n", s2);
    #endif
}

int score(char c1, char c2){
    #ifdef DEBUG
        printf("c1 : %c, c2: %c\n", c1, c2);
    #endif
    if (c1 == c2 == '*')
    {
        return 1;
    }

    if (c1 == '*' || c2 == '*')
    {   
        return -LGP;
    }
    
    if (c1 == c2)
    {
        for (uint64_t i = 0; i < SCORING_LEN; i++)
        {
            if (scoring_lookup[i] == c1)
            {
                return scoring_matrix[i][i];
            }
        }  
    }
    
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

    return scoring_matrix[i_c1][i_c2];
}

bool string_analizer(char * string){ //Controllo che i carateri prsenti in string siano solo quelli ammessi
    bool correct;
    size_t j;
    for (size_t i = 0; i < strlen(string); i++)
    {
        correct = false;
        j = 0;
        while (i < SCORING_LEN && !correct)
        {
            if (scoring_lookup[j] == string[i])
            {
                correct = true;
            }
            ++j;
        }
        if (!correct)
        {
            return false;
        } 
    }
    return true;
}

int max(int match, int indel_s1, int indel_s2){
    int max;
    if (match > indel_s1)
    {
        if (match > indel_s2)
        {
            max = match;
        }
        else //indel_s2 >= match
        {
            max = indel_s2;
        }  
    }
    else //indel_s1 >= match
    {
        if (indel_s1 > indel_s2)
        {
            max = indel_s1;
        }
        else //indel_s2 >= indel_s1
        {
            max = indel_s2;
        } 
    }

    #ifdef DEBUG
        printf("max: %d\n", max);
    #endif
    
    return max;   
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
    assert(x <= y || x - y <= band.base + band.extra);

    uint64_t width = band.base + 2 * band.extra;
    
    return (width * y) + (x - y + band.extra);
}

static int band_align(uint64_t base, uint64_t extra, char* s1, char* s2) {
    uint64_t l1 = strlen(s1);
    uint64_t l2 = strlen(s2);
    uint64_t band = base + 2 * extra;
    int match, indel_s1, indel_s2;

    cell_s* m = (cell_s*)malloc(band * (l1 + 1) * sizeof(cell_s));
    assert(m != NULL && "Cannot allocate matrix\n");

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
                printf("\nStepB: %ld, %ld\n", x, y);
            #endif
            cell_s* tp = m + conv(b, x, y);
            
            if ((x - 1) == 0 && (y - 1) == 0)
            {
                #ifdef DEBUG
                    printf("(x - 1) == (y - 1) == 0\n");
                #endif
                match = 0;
                indel_s1 = -LGP;
                indel_s2 = -LGP;
            }
            else if ((x - 1) == 0)
            {
                #ifdef DEBUG
                    printf("(x - 1) == 0\n");
                #endif
                match = -LGP * (y - 1);
                indel_s1 = m[conv(b, x, y - 1)].cell;
                indel_s2 = -LGP * y;
            }
            else if ((y - 1) == 0)
            {   
                #ifdef DEBUG
                    printf("(y - 1) == 0\n");
                #endif
                match = -LGP * (x - 1);
                indel_s1 = -LGP * x;
                indel_s2 = m[conv(b, x - 1, y)].cell;
            }
            else
            {
                #ifdef DEBUG
                    printf("else\n");
                #endif
                match = m[conv(b, x - 1, y - 1)].cell;
                indel_s1 = m[conv(b, x, y - 1)].cell ;
                indel_s2 = m[conv(b, x - 1, y)].cell; 
            }

            match += score(s1[y - 1], s2[x - 1]);
            indel_s1 += score('*', s2[x - 1]);
            indel_s2 += score(s1[y - 1], '*'); 

            tp->cell = max(match, indel_s1, indel_s2);

            #ifdef DEBUG
                printf("result: match : %d, indel_s1 : %d, indel_s2 : %d\n", match, indel_s1, indel_s2);    
            #endif

            if (match == tp->cell)
            {
                tp->prev_x = x - 1;
                tp->prev_y = y - 1;
            }
            else if (indel_s1 == tp->cell)
            {
                tp->prev_x = x;
                tp->prev_y = y - 1;
            }
            else //if (indel_s2 == tp->cell)
            {
                tp->prev_x = x - 1;
                tp->prev_y = y;
            }
            
            #ifdef DEBUG
                printf("result: %d, %ld, %ld\n", m[conv(b, x, y)].cell, m[conv(b, x, y)].prev_x, m[conv(b, x, y)].prev_y);
            #endif
        }
    }

    uint64_t x = l2;
    uint64_t y = l1;

    #ifdef DEBUG
        printf("StepC: Reconstruction\n");
    #endif

    for (cell_s c = m[conv(b, x, y)]; x != 0 && y != 0; x = c.prev_x, y = c.prev_y, c = m[conv(b, x, y)])
    {
        #ifdef DEBUG
            printf("StepC: %ld,%ld,%d,%ld,%ld=%c%c\n", x, y, c.cell, c.prev_x, c.prev_y, s1[y-1], s2[x-1]);
        #endif
        range_s r = range(b, y);
        if (x == r.sx && r.sx > 0 || x == r.dx && r.dx < l2)
        {
            printf("Hai raggiunto il confine della banda\n");
            return -1;
        }
        
    }

    x = l2;
    y = l1;

    uint64_t i = 0;

    char allignment_s1[BUZZ_SIZE] = "\0";
    char allignment_s2[BUZZ_SIZE] = "\0";

    for (cell_s c = m[conv(b, x, y)]; x != 0 && y != 0; x = c.prev_x, y = c.prev_y, c = m[conv(b, x, y)])
    {
        #ifdef DEBUG
            printf("%ld, %ld\n", x, y);
        #endif
        if (x != c.prev_x && y != c.prev_y)
        {
            allignment_s1[i] = s1[y - 1];
            allignment_s2[i] = s2[x - 1];
        }
        else if(y == c.prev_y)
        {
            allignment_s1[i] = '-';
            allignment_s2[i] = s2[x - 1];
        }
        else
        {
            allignment_s2[i] = '-';
            allignment_s1[i] = s1[y - 1];
        }
        i++;
    }

    if (x == 0)
    {
        while (y != 0)
        {
            #ifdef DEBUG
                printf("%ld, %ld\n", x, y);
            #endif
            allignment_s2[i] = '-';
            allignment_s1[i] = s1[y - 1];
            i++;
            y--;
        }
    }

    if (y == 0)
    {
        while (x != 0)
        {
            #ifdef DEBUG
                printf("%ld, %ld\n", x, y);
            #endif
            allignment_s1[i] = '-';
            allignment_s2[i] = s2[x - 1];
            i++;
            x--;
        }
    }  
    
    printf("\nRicostruzione avvenuta\n\n");

    allignment_s1[i] = '\0';
    allignment_s2[i] = '\0';

    reverse(allignment_s1, 0, strlen(allignment_s1) - 1);
    reverse(allignment_s2, 0, strlen(allignment_s2) - 1);

    printf("%s\n", allignment_s1);
    printf("%s\n", allignment_s2);
    printf("score : %d\n", m[conv(b, l2, l1)].cell);

    return 0;
}

void swap(char *str1, char *str2) 
{ 
    char temp[BUZZ_SIZE]; 
    strcpy(temp, str1);
    strcpy(str1, str2);
    strcpy(str2, temp);
}

int main(int argc, char const *argv[])
{
    char s1[BUZZ_SIZE];
    char s2[BUZZ_SIZE];

    read_text(s1, s2);

    to_upper_case(s1);
    to_upper_case(s2);
    
    if (strlen(s2) < strlen(s1)) 
    {
        swap(s1, s2); 
        #ifdef DEBUG
            printf("len1, len2 : %ld %ld\n", strlen(s1), strlen(s2));
            printf("s1, s2: %s,%s\n", s1, s2);
        #endif
    }

    uint64_t l1 = strlen(s1);
    uint64_t l2 = strlen(s2);
    uint64_t base = l2 - l1 + 1;
    uint64_t extra = 1;
    band_align(base, 0, s1, s2);
    return 0;

    #ifdef DEBUG
        printf("Instance lengths: %ld,%ld\n", l1, l2);
    #endif

    int ret_val = -1;

    for(;base + 2 * extra <= l2 + 1 && ret_val == -1; extra *= 2) 
    {
        printf("\n\nStep:\nbase : %ld, extra : %ld\n", base, extra);

        ret_val = band_align(base, extra, s1, s2);
    }
    if (ret_val == -1)
    {
        band_align(base, l1, s1, s2);
    }
    

    return 0;
}
