#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define intsize 64
typedef long long unsigned int typ;

typedef struct cell{
  typ * genotype;
  int * mut_pos;
  int num_mut;
}cell;

typedef struct generation{
  cell * cells; /* all the cells in this generation */
  int cell_num;
}generation;

char * bitfunc(typ n);

int frequency(int k);
