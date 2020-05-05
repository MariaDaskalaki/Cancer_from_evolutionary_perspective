#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>

#define intsize 64
typedef long long unsigned int typ;

typ global_fitfactor;
long double global_growth_rate;


char * bitfunc(typ n);

void sorting(typ* array, int size);

void freq_func(typ* array, int size);

unsigned int countSetBits(typ n);

void setbits_positions(int N, typ* pop,int total_num);

int* individual_setbits_positions(int n, typ pop);

double* sum_array(double* array, int size);

double find_min(double* array, int size);

int distinct_elements(typ* array, int size);


