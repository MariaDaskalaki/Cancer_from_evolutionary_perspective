# Cancer_from_evolutionary_perspective

Forward time simulations of cancer cells in discrete time of generations, based on two populations(previous cell population and current cell population).

The initial cell population consists of healthy cells.

The code is written in C and takes two command line parameters: initial number of cells and number of generations.

A "Frequency.txt" , a "Frequency_mutations.txt" and a "Population_size.txt" file are derived from the C-code. 

Plots of txt files are derived from the corresponding python codes.

The number of mutations for each generation is a random number from binomial distribution, inside the "pool" of total number of cells for the generation * genome_positions(intsize=64) with mutation rate mu=10^-4. The mutation rate of normal cells is between 10^-8 to 10^-6, while the mutation rate of cancer cells is obviously higher. (int binomial)

Mutations happen randomly in any random individual(cell) of the population in any random position. Each position of each individual can be mutated more than once. In code the mutation is represented by a "magic" number mask as the right-shift operator of 1 in the position of mutation, procedure based on bitwise operations. (int mutate)

With the functions unsigned int countSetBits, void setbits_positions and int* individual_setbits_positions we get respectively the number of total setbits(1s) for each individual(cell) and comulatively for each generation, so both inherited and new mutations are found, the positions that have been mutated for each generation and the mutated positions for each individual of each generation. 

The final attempt for finding a ftiness function was to give gravity and "mutation significancy" in certain positions of the genome. The number of positions that matter is random. Each genome position is controlled by a "fit_factor" based again on bitwise operations. The scores are given in the function typ fitfactors_func which returns the total "fit_factor_score", a bitwise number which represents the positions that matter. This number will serve as control in the function void com_fitness_func2, in order to check whether we have mutation in position that matters or not, consistently whether we have SMALL fitness (1.0) or LARGER fitness (described above).

As far as the interactions between genome positions are concerned:

1) Each genome position out of 64 is evalueated with either positive or negative number "weight" taken from gaussian distribution according to the law of great numbers, to which biological systems obey. (float* weight_position_array_func). 
2) As a result we have weight_position_array. (weight_position_array_func)
3) If we do not proceed futher each position is considered independet, which means we have no interactions between them (simple linear model).
4) With negative weights we take into account the negative interactions between positions. 
5) Each element of  new_weight_position_array is derived from the formula/equation 2*weight[i]+sum of the weights of all previous and next elements.
6) As a result, each position is affected mostly by itself (thats why coefficient 2) and from all previous and next positions. 
7) So now, we have multivariate linear model. 
8) These new scores will be the initial fitness scores.(new_weight_position_array)
9) However, due to the choice of gaussian distribution negative ftiness scores will arise. As fitness we define the probability of giving offsprings, so negative fitness scores are not possible. We accept values only >=0.
10) Thats why we had segemetation fault in our first attempt.
11) The first thought was to remove the negative "weights" but it didnt seem right and the next was to replae fitness values <0 with equal to 0 using if. If I would do this my LARGER fitness would always be less than my SMALL fitness=1.0 for positions that do not matter. So, in the end there would be fixation of positions that do not matter.
12) After a lot of search I found that roulette wheel selection algorithm keeps working for negative fitnesses after normalization of all fitness scores.
13) In this way, we could keep the negative weights and the final fitness scores would have positive values without changing the difference between them (which is important for selection).
14) So the norm_new_weight_position_array was created, each element of whichis derives from new_weight_position_array-minimum_element_of_new_weight_position_array. ( weight_selection_position1,find_min)

In function void com_fitness_func2 each individual of each generation is evaluated with a certain fitness score given as the sum of the weights of its positions that have been mutated.

For each generation the kids are identical to their parents become parents of the next generation, based on the idea of the identical division of the cell cycle. The probability for an individual to become parent is proportional to its fitness.
In function void selection_func we have implemeted the Rhoulette wheel selection algorithm: 
Calculate S = the sum of a finesses (comulative_fitness_array)
For each individual of the second population,generate a random number between 0 and S.
Starting from the top of the first population, keep adding the finesses to the partial sum P (com_fitness[k]), till P<S.
The individual for which P exceeds S is the chosen individual to become parent.
The selection_func also returns the growth_rate of each generation as the number of proliferative cells/ total number of cells(prolferative and non proliferative), with the help of the function int distinct_elements.

void print_genome_without , int print_genome_with are for printing each genotype as a binary string, with the help of char * bitfunc.

There have been two implementations based on two different cancer evolution approaches:
1) The Wright-Fisher model for non constant population size with 4 different cancer growth population models:
  population_array[i]=exponential_linear(initialN,i,growth_rate,2); // for 10000 exp until 3 generation.
  population_array[i]=generalized_logistic(initialN,i,(double)pow((double)10,(double)4),growth_rate,0.25); //capacity~10^12 check paper Estimating tumor growth rates in vivo, for the parameters.
  population_array[i]=Gompertz_model(initialN,i,growth_rate,gsl_ran_lognormal(r,-2.9,0.71)); // check paper Estimating tumor growth rates in vivo, for parameters.
  population_array[i]=von_bertalanffy(initialN,i,growth_rate,gsl_ran_lognormal(r,-2.9,0.71));
 For each population model there is a threshold of 100000 cells for each generation.
 2) The branching process for stohastic non constant population size based on Poisson distribution:
  int * population_change_via_branching(gsl_rng* r, int initialN, int generations, double lambda, double ChangingFactor)

Main function consists of the following parts:
1) Construction of population array based on certain population model
2) Loop of ngenerations with biological procecess:
  int numberOfMutations=binomial(r,mu,size)
  print_genome_without(size,pop_prev,i)
  mutation_process(numberOfMutations,total_pop_pos,pop_prev,i)
  int total_setbits_num=print_genome_with(size,pop_prev)
  freq_func(pop_prev,size)
  selection_func(weight_position_array,r,size,next_size, pop_prev, pop_cur, fitness, com_fitness,selected_pos)
  setbits_positions(size,pop_prev,total_setbits_num)
  population_realloc(next_size, pop_prev, pop_cur) ( the kids become parents)
 3) At the end of generations loop the function void fixed_positions prints the number of positions==mutations that have been fixed and the number of position==mutations that fave been fixed from those that matter.


The two main codes are pop_sim_interactions_with_gauss_weights.c AND pop_sim_branching.c
The runs are for 100 intilial population and 1500 generations.
The population reaches up to 100000 (for example in Gompertz_model).
As the populations size increases the fixation becomes more difficult even if the mutation rate is interrupted after 300 generations.



