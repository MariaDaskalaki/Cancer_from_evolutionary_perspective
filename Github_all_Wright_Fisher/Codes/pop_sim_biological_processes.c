#include "pop_sim_biological_processes.h"

int mutate(typ input, typ *pop, int generation)
{
  
  typ k;
  //  if (generation == 10) // για να συμβει σίγουρα η μετάλλαξη στη θέση 0 στη γενιά 10.                                                                                                                   
  // k = 0;                                                                                                                                                                                                 
  // else                                                                                                                                                                                                   
  k= rand() % (input);

  int individual = k/intsize;
  int position = k % intsize;

  typ mask = (typ)1 << position;
  pop[individual] = pop[individual] ^ mask;
  fprintf(stderr, "k:%llu, %llu\n", k, mask);
  fprintf(stderr, "The mutation has happened in the %d individual of the current population and in POSITION %d\n",individual,position);


  return position;
}

int binomial(gsl_rng * r, double mu, int N){
long long int total_pop_pos=N*intsize;
int total_num_mut = gsl_ran_binomial (r, mu, total_pop_pos);
  // printf("Total Number of Mutations is %d\n", total_num_mut);                                                                                                                                            

// printf ("\n");                                                                                                                                                                                           

  return total_num_mut;
}


typ fitfactors_func(unsigned k){
  typ* fitfactor_array=(typ*)calloc(intsize,sizeof(typ));
  typ* fitfactor_array_final=(typ*)calloc(intsize,sizeof(typ));
  typ fitfactor=0;
  //printf("The fitfactors masks are\n");                                                                                                                                                                  \
                                                                                                                                                                                                            
  for (unsigned i=0;i<intsize;i++){
    fitfactor_array[i]=(typ)1 << i;
    fitfactor=fitfactor | (fitfactor_array[i]);
    fitfactor_array_final[i]=fitfactor;
    //printf("The fitfactor is %llu\n",fitfactor_array_final[i]);	\
                                                                                                                                                                                                            
    //printf("The %d fitfactor is %llu\n",k,fitfactor_array_final[k]);                                                                                                                                     \
                                                                                                                                                                                                            

  }
  //printf("The %d fitfactor is %llu\n",k,fitfactor_array_final[k-1]);                                                                                                                                      

  typ ret = fitfactor_array_final[k-1];
  free(fitfactor_array);
  free(fitfactor_array_final);
  return ret;
}

float weight_selection_position(unsigned k){
  float* weight_position_array=(float*)calloc(intsize,sizeof(float));
  float* com_weight_position_array=(float*)calloc(intsize,sizeof(float));
  for (int i=0;i<intsize;i++){
    weight_position_array[i]=(float)rand()/(float)RAND_MAX * 10;
    com_weight_position_array[0]=weight_position_array[0];
    if (i>0){
      com_weight_position_array[i]=weight_position_array[i]+com_weight_position_array[i-1];
    }
    //printf("The weights for each position are\n");                                                                                                                                                        
    //printf("%f\n",weight_position_array[i]);                                                                                                                                                              
    //printf("The com_weights are\n");                                                                                                                                                                      
    //printf("%f\n",com_weight_position_array[i]);                                                                                                                                                          
  }
  //float max=weight_position_array[0];                                                                                                                                                                     
  //int selected_pos=0;                                                                                                                                                                                     
  //for (int i=0;i<intsize;i++){                                                                                                                                                                            
  //if (weight_position_array[i]>max){                                                                                                                                                                      
  //max=weight_position_array[i];                                                                                                                                                                           
  //selected_pos=i;                                                                                                                                                                                         
  //}                                                                                                                                                                                                       
  //}                                                                                                                                                                                                       
  //printf("The max value is %f in position %d\n",max,selected_pos);                                                                                                                                        
  free(weight_position_array);
  float x = com_weight_position_array[k];
  free(com_weight_position_array);
  return x;
}

float weight_selection_position1(float* array,gsl_rng* r,int size,int* array1){
  
  float* new_weight_position_array=(float*)calloc(intsize,sizeof(float));
  float* norm_new_weight_position_array=(float*)calloc(intsize,sizeof(float));
  float com_weight_position;

  new_weight_position_array=sum_array(array,intsize);


  float min_fitness=find_min(new_weight_position_array,intsize);

  for (int i=0; i<intsize; i++){
    norm_new_weight_position_array[i]=(float)new_weight_position_array[i]-(float)min_fitness;
  }
  
  for (int i=0; i<size; i++){
    com_weight_position+=norm_new_weight_position_array[array1[i]];
  }


  float x=com_weight_position;

  free(new_weight_position_array);
  free(norm_new_weight_position_array);
 
  return x;
}


void com_fitness_func2(float* array,gsl_rng* r,int N, double* fitness, double* com_fitness, typ *pop, int selected_pos){
  // int k=rand()%64 +1;                                                                                                                                                                                    
  //printf("The positions that matter are %d\n",k);                                                                                                                                                          
  typ fitfactor=fitfactors_func(selected_pos);
  //int* each_individual_mutations;
  for (int i=0;i<N;i++){
    fprintf(stderr,"%llu %llu %llu\n",pop[i],fitfactor,pop[i] &fitfactor);
    if ((typ)(pop[i] & fitfactor) == (typ)0){
      fprintf(stderr,"SMALL FIT\n");
      fitness[i]=1.;
    }
    else{
      int num_of_individuals_muts=countSetBits(pop[i]);
      //fprintf(stderr,"the total number of mutations for individual %d is %d\n",i,num_of_individuals_muts);
      int* individual_array= individual_setbits_positions(num_of_individuals_muts,pop[i]);
      fprintf(stderr,"LARGER FIT\n");
      fitness[i]= weight_selection_position1(array,r,num_of_individuals_muts,individual_array); // υποτίθεται είναι με τις αλληλεπιδράσεις (πολυμεταβλητό γραμμικό μοντέλο), ενώ το weight_selection_position είναι με το απλό γραμμικό μοντέλο όπου κάθε θέση είναι ανεξάρτητη.
    
    }
    fprintf(stderr, "Fitness of %d is  %f\n", i, fitness[i]);
    com_fitness[0]=fitness[0];
    if (i>0){
      com_fitness[i]=fitness[i]+com_fitness[i-1];
      fprintf(stderr, "Com_fitness of %d is %f\n",i,com_fitness[i]);
    }
  }
  global_fitfactor = fitfactor;
}


void com_fitness_func1(int N, double *fitness, double *com_fitness,int k, typ *pop,int fitness_weight){

  typ fitFactor = (typ)1 << k;

  for (int i=0;i<N;i++){
    fprintf(stderr,"%llu %llu %llu\n", pop[i], fitFactor, pop[i] & fitFactor);
    if( (typ)(pop[i] & fitFactor) == (typ)0 )
      {
        fprintf(stderr,"SMALL FIT\n");
        fitness[i] = 1.;
      }
    else
      {

        fprintf(stderr,"LARGER FIT\n");
        fitness[i] = fitness_weight;

      }

    fprintf(stderr, "Fitness of %d is %f\n", i, fitness[i]);

    com_fitness[0]=fitness[0];
    if (i>0){
      com_fitness[i]=fitness[i]+com_fitness[i-1];
    }
  }
}

void selection_func(float* array,gsl_rng* r,int N1, int N2,typ* pop1, typ* pop2, double *fitness, double* com_fitness,int selected_pos){

  typ k;
  typ* parent_array=(typ*)calloc(N2,sizeof(typ));                                                                                                                                                         
  //com_fitness_func(N, fitness, com_fitness);                                                                                                                                                              
  //com_fitness_func1(N, fitness, com_fitness, weight_selection_position(), pop1,2000);                                                                                                                     
  com_fitness_func2(array,r,N1,fitness,com_fitness,pop1,selected_pos);

  fprintf(stderr,"Total fitness is %f\n", com_fitness[N1-1]);
  //double rand_num=(((double)rand() / (double)RAND_MAX) * com_fitness[N1-1]);
  //int sum=0;
  for (int i=0;i<N2;i++){
    double r=(double)rand()/RAND_MAX * com_fitness[N1-1];
    //double r=(double)rand() / ((double)RAND_MAX + 1) * com_fitness[N1-1];
      for (k=0;k<N1;k++){
        
	if (com_fitness[k]>r) break;
      }

    typ parent = k;
    assert(parent<N1);
    pop2[i] = pop1[parent];
    fprintf(stderr,"The individual %d of next population has as parent the individual %llu of the current population\n",i, parent);
    parent_array[i]=parent;                                                                                                                                                                               
    printf("%llu\n",parent_array[i]);                                                                                                                                                                    
  }

  unsigned int distinct=distinct_elements(parent_array,N2);                                                                                                                                               
  long double growth_rate=(long double)distinct/(long double)N1;                                                                                                                                          
  fprintf(stderr,"The number of proliferative individuals(parents) are %d\n",distinct);                                                                                                                   
  fprintf(stderr,"The growth rate id %Lf\n",growth_rate);                                                                                                                                                 
  //global_growth_rate=growth_rate;


  free(parent_array);                                                                                                                                                                                     
}


float* weight_position_array_func(gsl_rng*r,int size){

  float* weight_position_array1=(float*)calloc(size,sizeof(float));

  for (int i=0;i<size;i++){
    weight_position_array1[i]=gsl_ran_gaussian(r,1);
    // weight_position_array[i]=(float)rand()/(float)RAND_MAX*10-11;                                                                                                                                        
  }
  return weight_position_array1;
  free(weight_position_array1);
}


int * population_change(int initialN, int generations,long double growth_rate,gsl_rng* r) // growth_rate=ln(2)/Tc, where Tc=doubling time of cells and here it is considered to be 1 generation.
{

    int* population_array=(int*)calloc(generations,sizeof(int));
    population_array[0]=initialN;
    for (int i = 1; i < generations; ++i)
    {
      //population_array[i]=rand()%1000 + initialN; //rand() % N + A will give you a value in the range [A, A + N], not [A, N]
      //population_array[i]=exponential_linear(initialN,i,growth_rate,2); // for 10000 exp until 3 generation.
      population_array[i]=generalized_logistic(initialN,i,(double)pow((double)10,(double)4),growth_rate,0.25); //capacity~10^12 check paper Estimating tumor growth rates in vivo, for the parameters.
      //population_array[i]=generalized_logistic(initialN,i,10000,growth_rate,0.3); // check paper Cellular interactions constrain tumor growth, for parameters.
      //population_array[i]=Gompertz_model(initialN,i,growth_rate,gsl_ran_lognormal(r,-2.9,0.71)); // check paper Estimating tumor growth rates in vivo, for parameters.
      //population_array[i]=von_bertalanffy(initialN,i,4,2); // check paper Cellular interactions constrain tumor growth, for parameters.
      //population_array[i]=von_bertalanffy(initialN,i,growth_rate,gsl_ran_lognormal(r,-2.9,0.71));
      

    }

    return population_array;
}


void print_genome_without(int N, typ* pop, int generation){

  for(int k = 0; k < N; ++k){
      char* s=bitfunc(pop[k]);
      fprintf(stderr,"%d  -- The genome (WITHOUT) of individual %d of the current population is %s\n", generation, k,s);
      free(s);
    }
}


void mutation_process(int number_of_mutations, long long int positions, typ* pop, int generation){
  
  typ* position_of_mutations_array=(typ*)calloc(number_of_mutations,sizeof(typ));
    for (int j=0;j<number_of_mutations;j++){
      position_of_mutations_array[j]=mutate(positions, pop, generation); // οι θέσεις που έγιναν μεταλλάξεις μπαινουν σε ενα array
      //printf("Mutation_array is %d\n",position_of_mutations_array[j]);                                                                                                                                    
      //mutate(total_pop_pos, pop_prev, i);                                                                                                                                                                 

    }

    free(position_of_mutations_array);
}


int print_genome_with(int N, typ* pop){
  int x=0;
  for(int k = 0; k < N; ++k){
      char* s=bitfunc(pop[k]);
      int num_of_setbits=countSetBits(pop[k]);
      x+=num_of_setbits;
      fprintf(stderr,"The genome (WITH) of individual %d of the current population is %s\n",k,s);
      fprintf(stderr,"The individual %d has %d number of setbits/mutations\n",k,num_of_setbits);
      individual_setbits_positions(num_of_setbits,pop[k]);
      free(s);
  }

  if (x==0){
    x=1;
  }
 
 
  fprintf(stderr,"The total number of setbits for this generation are %d\n",x);

  return x;
}


void population_realloc(int N, typ* pop1, typ* pop2){

  for (int k=0; k<N; k++){
      pop1[k]=pop2[k];
      // char* s=bitfunc(pop_cur[k]); /                                                                                                                                                                     
      // printf("The genome (without the possible upcoming mutation) of individual %d of the current population is %s\n",k,s); /                                                                            
    }
}

  

void fixed_positions(int N, typ* pop){

  for (int i=0;i<N;i++){
    int setBits = countSetBits(pop[i]);
    int commonBits = countSetBits(pop[i] & global_fitfactor); // οι κοινες θεσεις ( άσσοι) μεταξύ των άσσεων του fitfactor και του μεταλλαγμένου γονιδιώματος).                                         
    fprintf(stderr,"Individual:%d\t Number of positions that have been fixed (frome those that matter):%d\t Number of positions that have contibuted (from those that matter):%d\n", i, setBits, commonBits\
);

  }
}
