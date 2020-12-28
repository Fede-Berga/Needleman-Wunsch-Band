#ifndef SW_H
#define SW_H
#define SCORING_LEN 24
#define LGP 4
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>

/**
 * Definisce le informazioni della matrice di programmazione dinamica.
 */
typedef struct cell_s {
        int cell; ///< Soluzione al sottoproblema corrisondente.
        uint64_t prev_x; ///< Link per la ricostruzione in x.
        uint64_t prev_y; ///< Link per la ricostruzione in y.
} cell_s; 

/**
 * Definisce i limiti destro e sinistro della banda.
 */
typedef struct range_s {
        uint64_t sx; ///< Limite sinistro della banda.
        uint64_t dx; ///< Limite destro della banda.
} range_s;

/**
 * Definisce informazioni sulla banda.
 */
typedef struct band_s {
        uint64_t l1; ///< Lunghezza della prima stringa.
        uint64_t l2; ///< Lunghezza della seconda stringa.
        uint64_t base; ///< Ampiezza minima della banda.
        uint64_t extra; ///< Ampliamento della banda.
} band_s;

/**
 * Matrice di score.
 */
const int scoring_matrix[SCORING_LEN][SCORING_LEN] = {
    { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -LGP},
    {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -LGP},
    {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -LGP},
    {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -LGP},
    { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -LGP},
    {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -LGP},
    {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -LGP},
    { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -LGP},
    {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -LGP},
    {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -LGP},
    {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -LGP},
    {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -LGP},
    {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -LGP},
    {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -LGP},
    {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -LGP},
    { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -LGP},
    { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -LGP},
    {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -LGP},
    {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -LGP},
    { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -LGP},
    {-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -LGP},
    {-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -LGP},
    { 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -LGP},
    {-LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP, -LGP,  LGP}
};

/**
 * Chiavi della matrice di score.
 */
char scoring_lookup[SCORING_LEN] = {'A',  'R',  'N',  'D',  'C',  'Q',  'E',  'G',  'H',  'I',  'L',  'K',  'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V',  'B',  'Z',  'X',  '*'};

/**
 * Effetua l'analisi della stringa in imput.
 * In particolare controlla che i caratteri nella string siano
 * inclusi in scoring_lookup.
 * 
 * @param string Stringa da analizzare
 * 
 * @return true se la stringa rispetta i parametri stabiliti.
 */
bool string_analizer(char * string);
bool do_string_analizer(char * string, uint64_t i);

/**
 * Converte una stringa lower case in una stringa upper case.
 * 
 * @param string Stringa da convertire.
 * 
 * @post La stringa è convertita.
 */
void to_upper_case(char * string);
void do_to_upper_case(char * string, uint64_t i);

/**
 * Effettua il reverse di una stringa.
 * 
 * @param string La stringa da invertire.
 * 
 * @post La stringa è invertita.
 */
void reverse(char * string);
void do_reverse(char * string, int begin, int end);

/**
 * Legge due stringhe.
 * 
 * @param s1 Prima stringa da leggere.
 * 
 * @param s2 Seconda stringa da leggere.
 * 
 * @post Le stringhe sono lette.
 */
void read_text(char * s1, char * s2);

/**
 * Restituisce il puteggio relativo ai due caratteri, contenuto nella matrice di score.
 * 
 * @param c1 Primo carattere
 * @param c2 Secondo carattere.
 */
int score(char c1, char c2);

/**
 * Restituisce il massimo tra tre interi.
 * 
 * @param match Primo numero.
 * 
 * @param indel_s1 Secondo numero.
 * 
 * @param indel_s2 Terzo numero.
 * 
 * @return Massimo tra i tre numeri.
 */
int max(int match, int indel_s1, int indel_s2);

/**
 * Per una quota y restutuisce i limiti destro e sinistro della banda
 * 
 * @param band Informazioni sulla banda.
 * 
 * @param y Quota.
 * 
 * @return limiti destro e sinistro della banda.
 */
range_s range(band_s band, uint64_t y);

/**
 * Linearizza le cordinate della banda.
 * 
 * @param band Informazioni sulla banda.
 * 
 * @param x Ascissa non linearizzata.
 * 
 * @param y Ordinata non linearizzata.
 * 
 * @return Coordinata linearizzata.
 */
uint64_t conv(band_s band, uint64_t x, uint64_t y);

/**
 * Effettua l'allineamento di due stringhe.
 * 
 * @param base Numero minimo di operazioni per allineare le due stringhe.
 * 
 * @param extra Ampliamento della banda.
 * 
 * @param s1 Prima stringa.
 * 
 * @param s2 Seconda stringa.
 * 
 * @return -1 se la Banda non è adeguatamnte ampia, 0 altrimenti.
 */
static int band_align(uint64_t base, uint64_t extra, char* s1, char* s2, char * allignment_s1, char * allignment_s2);

/**
 * Inverte due stringhe.
 * 
 * @param s1 Prima stringa.
 * 
 * @param s2 Seconda stringa.
 */
void swap(char *str1, char *str2);

#endif