#include "pop_sim_libraries.h"
#include "pop_sim_pop_models.h"

int mutate(typ input, typ *pop, int generation);

int binomial(gsl_rng * r, double mu, int N);

typ fitfactors_func(unsigned k);

float weight_selection_position(unsigned k);

float weight_selection_position1(float* array,gsl_rng* r,int size,int* array1);

void com_fitness_func2(float* array,gsl_rng* r,int N, double* fitness, double* com_fitness, typ *pop, int selected_pos);

void com_fitness_func1(int N, double *fitness, double *com_fitness,int k, typ *pop,int fitness_weight);

void selection_func(float* array,gsl_rng* r,int N1, int N2,typ* pop1, typ* pop2, double *fitness, double* com_fitness,int selected_pos);

float* weight_position_array_func(gsl_rng*r,int size);

int * population_change(int initialN, int generations,long double growth_rate,gsl_rng* r);

void print_genome_without(int N, typ* pop, int generation);

void mutation_process(int number_of_mutations, long long int positions, typ* pop, int generation);

int print_genome_with(int N, typ* pop);

void population_realloc(int N, typ* pop1, typ* pop2);

void fixed_positions(int N, typ* pop);


