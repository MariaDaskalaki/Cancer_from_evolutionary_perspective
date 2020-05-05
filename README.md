# Cancer_from_evolutionary_perspective

Forward time simulations of cancer cells in discrete time of generations, based on two populations(previous cell population and current cell population).

The initial cell population consists of healthy cells.

The code is written in C and takes two command line parameters: initial number of cells and number of generations.

For the implementation pop_sim_interactions_with_gauss_weights.c: 
  1) pop_sim_libraries.h must be at the current directory
  2) gcc -o pop_sim_interactions_with_gauss_weights pop_sim_interactions_with_gauss_weughts.c -g -lgsl -lgslcblas -lm
  3) ./pop_sim_interactions_with_gauss_weights <number of cells> <numberof generations>
  
For the implementation pop_sim_CA.c:
  1) pop_sim_libraries_CA.h must be at the current directory
  2) gcc -o pop_sim_CA pop_sim_CA.c -g -lgsl -lgslcblas -lm
  3) ./pop_sim_CA <number of cells> <number of generations>

 From pop_sim_interactions_with_gauss_weights.c a "Frequency.txt" , a "Frequency_mutations.txt" and a "Population_size.txt" file are derived from the C-code. 

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

POP_SIM_CA: Implementing Gompertz model with one-dimensional cellular automata

1) A one dimensional cellular automaton is proposed for implementing the Gompertz model to formulate tumor growth dynamics, which is called Gompertz CA.

2) First a discrete model is derived from standard Gompertz model. Then, stohastically evolutionary rules for the Gompertz CA are developped based on fitness (genome-position dependent) and the discrete Gompertz model. The Gopertz CA is therefore randomized by the stohastically evolutionary rules, so that it can be represented with a stohastic process or a stohastic differential equation.

3) A growing tumor is regarded as a self-organized complex evolutionary system. Tumors, in terms of complex systems, can be analogous to cellular automata and tumor cells to automaton cells.

4) The standard conitnious Gompertz model:
The Gompertz model is the best-known mathematical equation for modeling tumor growth:

V(t)=Voexp(A/B(1-exp(-Bt), V(t=0)=Vo

5) The discrete Gompertz model:

dV(t)/dt=AVoexp(-Bt)exp(A/B(1-exp(-Bt))

Let t--->0 (small enough) then 

dV(t)/dt=(V(T+to)-V(t))/To, where t now is a discrete time, t=kTo (k=1,2,...) and To is the discrete time interval.

So, V((k+1)To)=V(kTo)+AVoToexp(A/B(1-exp(-BkTo))), discrete Gompertz model (undelies the Gompertz CA as the one of the two evolutionary stohastic rules).

6) The Compertz CA is constructed with a one dimensional cellular automaton:
                  TGM=<t,Cells,CellSpace,Neighbors,Rules>
   The Gompertz CA defines its cell neighborhood according to the Moore neighborhood: Ngh(i)={Cell(i-1), Cell(i+1)}.
   The cell neighborhood is a field of evolutionary rules that regulates the states of Cell(i)
   
7) Gompertz CA cells are divided into normal cells and abnormal(tumorous cells), which are represented with different states of Gompertz CA cells: Si(t) in {0,1}, where 0 indicates that cell(i) is normal at discrete time t and 1 indicates an abnormal Gompertz CA cell (tumorous).

8) The evolutionary rules of Gompertz CA:
The evolution of a given cell depends on the state of the two neighbors and the fitness of cell(i). A normal CA cell probably at discrete time t probably evolves into a tumorous CA cell at the next time t+1.

The state of any Gompertz CA cell is not reversible: 0 can become 1 but 1 can not bacome 0. Gompertz CA cells take ecolutionary computing simoultaneously. At time t=0, the state if any Gompertz CA cell expept those that become parents ( so increased fitness) is zero.

9) Rules: 

Si(t+To)=si(t), if si(t)=1 or pi(t)<=P and si(t)+1, otherwise

pi(t)=po, if si-1(t)=0 and si+1(t)=0 and individuals_condition=1(see 13) else pi(t)=0, AND 2po if si-1(t)=1 and si+1(t)=1, AND po(1-ΔV(t)/ΔVmax), otherwise

where pi(t) is called the evolutionary probability of Cell(i) at discrete time t
      P is the probable threshold belonging to [0,1] and generated with a uniform distribution (probability function)
      po is the fitness of Cell(i) at discrete time t
      ΔV(t) and ΔVmax are the ideal increment and the probable maximum increment in the volume of the tumor simulated by the Gompertz model CA at discrete time t.
      
      ΔV(t)=V(t+To)-V(t)=AVoToexp(-Bt)exp(A/B(1-exp(-Bt)))
      
      ΔVmax=constant= Gompertz capacity
      
10) Rules Explanation:
The abnormal cell (Si(t)=1) will give abnormal ( 0--->1 but not 1---->0)
The healthy cell can probably give abnormal
If we have a healthy cell and a neighborhood of healthy cells then it can become abnormal with probability proportional to its fitness.
If we have a healthy cell and a neighborhhod of abnormal cells then it will become abnormal with probability propotrional to two times of its fitness.
If we have a healthy cell and a neighborhood of one healthy and one abnormal cell then it will become abnormal with probability proporional to its fitness and its available space according to Gompertz capacity.

11. In order to take into account the spatial heterogeneity due to cell-environment interactions we have placed the cells in a 1D array of size N (cur_index_array), so now each offspring has a certain position close to its parent position. In order to find the position of each offspring we have sorted the offsprings array according to the parent index array. (typ* sorting_1).

12. The spatial heterogeneity is implied by different fitness values according to site. As we get closer to the centre of the array fitness values decrease in order to show the environmental stress due to low oxygen ( hypoxia) or overpopulation ( necrosis), the opposite happens as we get closer to the edges of the array. We used  values from bimodal distribution from a linear function decreasing for the first half of the array and increasing for the other half. The maximun of the function is 10 and the minimum is 1 (double* fitness_site_func). As a result each individual has a fitness according to its position(fitness_individuals_site_array).

13. If a cell is neutral (undergone 0 mutations) with individuals_condition=0 ( neutrality) it means that its fitness value  depends only on its site. If a cell has undeergone mutations in beneficial position ( non neutral/mutant) with individuals_condition=1 then its fitness value is the sum of fitnesses of the positions of mutations and the fitness of its site. 

14. In conclusion, in this implementation we have taken into account the genome position interactions,   the cell- to cell interactions, resulting in modelling the invasion of tumorous cells into a healthy tissue and the environmnet -to- cell interactions resulting in delaying of fixation time (as we have a low connectivity graph implementation)  in comparison with the space free Gompertz.

15. Derived .txt files: Frequencies_CA.txt ( genotypes:frequencies), Frequencies_CA_mutant.txt(mutant genotypes:frequencies),
Frequencies_mutations_CA.txt(all genotypes mutations:frequencies), Frequencies_mutations_CA_mutant.txt(mutant genotypes mutations:frequencies), Number of mutants.txt (Generation '\t' Number of mutants).

Both of the implementations are in the directory Comparison.



