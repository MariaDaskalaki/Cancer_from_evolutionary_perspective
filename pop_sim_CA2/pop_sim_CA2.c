#include "pop_sim_libraries.h"


typ global_fitfactor;
long double global_growth_rate;


char * bitfunc(typ n){
  int nbits = intsize;
  //char *s=malloc(nbits+1);                                                                                                                                                                                
  char * s=(char*)calloc(nbits + 1,sizeof(char));
  s[nbits]='\0'; // The \0 is treated as NULL character. The applications of '\0' lies in determining the end of string. Here it says that the end of the string will be the nbits=32th digit.              
  //assert(n > 0);                                                                                                                                                                                          
  //unsigned int u=*(unsigned int*)&n; // store integer n in a positive(as we have unsigned type) pointer variable of integer u. & is a misc operator in c and &a returns the actual address of the variabl\
e a.                                                                                                                                                                                                        
  int i;
  typ mask = (typ)1 << (nbits-1); // << left shift operator. The left shift and right shift operators should not be used for negative numbers. The result of is undefined behaviour if any of the operands \
is a negative number. Thats why we use unsigned type. If the number is shifted more than the size of integer the behaviour is undefined. For example 1 << 33 is undefined if integers are stored using 32 b\
its. Here we have stored integers of sixe nbits=intsize=define 32. >> Binary Left Shift Operator. The left operands value (here is 1) is moved left by the number of bits specified by the right operand ( \
here is 32). So here its like filling the 32 positions from right to left.                                                                                                                                  
  for (i=0;i<nbits;i++,mask>>=1){ // mask>>=1 is the same as mask=mask>>1 so the left from above 1 is moved right by 1 bit.                                                                                 
    s[i]=((n&mask)!=0)+'0'; // & Binary AND Operator copies a bit to the result if exists in both operands. Here we check if 1 (which is !=0) is in both u and mask we place and what is left is filled up \
with 0.                                                                                                                                                                                                     
  }
  // printf("%llu, %llu, %s\n", n, mask, s);                                                                                                                                                                

  return s;
}


void sorting(typ* array, int size){
  //fprintf(stderr,"The elements in ascending order are\n");                                                                                                                                               \
                                                                                                                                                                                                            
  for (int i=0;i<size;i++){
    for (int j=i+1;j<size;j++){
      if (array[i]>array[j]){
        typ temp=array[i];
        array[i]=array[j];
        array[j]=temp;
      }
    }
    if (array[i]!=0){
      FILE *f1=fopen("Frequencies.txt","a");
      fprintf(f1,"%Lf\t",(long double)array[i]/(long double)size);
      fclose(f1);

    }
  }
}



void freq_func(typ* array, int size){
  typ* freq=(typ*)calloc(size,sizeof(typ));
  for(int i=0; i<size; i++){
        freq[i] = -1;
    }
  for(int i=0; i<size; i++){
    int count=1;
    for(int j=i+1; j<size; j++){
      if(array[i]==array[j]){
        count++;
        freq[j] = 0;
      }
    }
    if(freq[i] != 0){
      freq[i] = count;
    }
  }

  for(int i=0; i<size; i++){

    if(freq[i] != 0){
      
	char * s = bitfunc(array[i]);
	fprintf(stderr, "%llu:The genome %s has frequency %f which means that occurs %llu times\n",array[i], s, (double)freq[i]/(double)size,freq[i]);
	FILE *f1=fopen("Frequencies_CA.txt","a");                                                                                                                                                             
	fprintf(f1, "%llu:%Lf\t",array[i],(long double)freq[i]/(long double)size);
	fclose(f1);
	free(s);
    
     
    }
  }


  //sorting(freq,size);                                                                                                                                                                                     
  free(freq);
}

void freq_func_mutant(typ* array, int size){
  typ* freq=(typ*)calloc(size,sizeof(typ));
  for(int i=0; i<size; i++){
        freq[i] = -1;
    }
  for(int i=0; i<size; i++){
    int count=1;
    for(int j=i+1; j<size; j++){
      if(array[i]==array[j]){
        count++;
        freq[j] = 0;
      }
    }
    if(freq[i] != 0){
      freq[i] = count;
    }
  }

  for(int i=0; i<size; i++){

    if(freq[i] != 0){
      
	char * s = bitfunc(array[i]);
	fprintf(stderr, "%llu:The mutant genome %s has frequency %f which means that occurs %llu times\n",array[i], s, (double)freq[i]/(double)size,freq[i]);
	FILE *f3=fopen("Frequencies_CA_mutant.txt","a");                                                                                                                                                             
	fprintf(f3, "%llu:%Lf\t",array[i],(long double)freq[i]/(long double)size);
	fclose(f3);
	free(s);
    
     
    }
  }


  //sorting(freq,size);                                                                                                                                                                                     
  free(freq);
}




typ Gompertz_discrete_model(int initial, int t, long double growth_rate, long double decay_rate){
  long double c;
  c=(long double)growth_rate/(long double)decay_rate;
  typ newSize=(growth_rate*initial)*(exp(-decay_rate*t))*(exp(c*(1-exp(-decay_rate*t))));

  if (newSize > 10000){
    newSize=10000;
  }
  return newSize;
}




int mutate(typ input, typ *pop, int generation)
{

  typ k;
  //  if (generation == 10) // για να συμβει σίγουρα η μετάλλαξη στη θέση 0 στη γενιά 10.                                                                                                                  \
                                                                                                                                                                                                            
  // k = 0;                                                                                                                                                                                                \
                                                                                                                                                                                                            
  // else                                                                                                                                                                                                  \
                                                                                                                                                                                                            
  k= rand() % (input);

  typ individual = k/intsize;
  int position = k % intsize;

  typ mask = (typ)1 << position;
  pop[individual] = pop[individual] ^ mask;
  fprintf(stderr, "k:%llu, %llu\n", k, mask);
  fprintf(stderr, "The mutation has happened in the %llu individual of the current population and in POSITION %d\n",individual,position);


  return position;
}





int binomial(gsl_rng * r, double mu, int N){
long long int total_pop_pos=N*intsize;
int total_num_mut = gsl_ran_binomial (r, mu, total_pop_pos);
  // printf("Total Number of Mutations is %d\n", total_num_mut);                                                                                                                                            

// printf ("\n");                                                                                                                                                                                           

  return total_num_mut;
}


unsigned int countSetBits(typ n){
  unsigned int count=0;
  while(n){
    count+=n & 1;
    n >>=1;
  }
  return count;
}



void setbits_positions(int N, typ* pop,int total_num){
  int* position_array=(int*)calloc(64,sizeof(int));
  for (int i=0;i<N;i++){
    unsigned int setbit_pos=0;
    while(pop[i]){
      if((pop[i]&1)!=0){
        fprintf(stderr, "The setbit position is %d\n",setbit_pos);
        position_array[setbit_pos]+=1;
      }
      setbit_pos++;
      pop[i]>>=1;
    }
  }
  for (int j=0;j<64;j++){
    FILE* f2=fopen("Frequencies_mutations_CA.txt","a");
    fprintf(stderr,"The mutation in position %d has occured %d times\n",j,position_array[j]);
    fprintf(f2,"%d:%f\t",j,(double)position_array[j]/((double)total_num));   
    fclose(f2);
  }
  free(position_array);
}

void setbits_positions_mutant(int N, typ* pop,int total_num){
  int* position_array=(int*)calloc(64,sizeof(int));
  for (int i=0;i<N;i++){
    unsigned int setbit_pos=0;
    while(pop[i]){
      if((pop[i]&1)!=0){
        fprintf(stderr, "The setbit position is %d\n",setbit_pos);
        position_array[setbit_pos]+=1;
      }
      setbit_pos++;
      pop[i]>>=1;
    }
  }
  for (int j=0;j<64;j++){
    FILE* f4=fopen("Frequencies_mutations_CA_mutant.txt","a");
    fprintf(stderr,"The mutant mutation in position %d has occured %d times\n",j,position_array[j]);
    fprintf(f4,"%d:%f\t",j,(double)position_array[j]/((double)total_num));   
    fclose(f4);
  }
  free(position_array);
}

int* individual_setbits_positions(int n, typ pop){
  int* individual_muts_array=(int*)calloc(n,sizeof(int));
  unsigned int setbit_pos=0;
  int i=0;
  while (pop){
    if (pop&1==1){
      //fprintf(stderr,"The mut is in position %d\n",setbit_pos);
      individual_muts_array[i]=setbit_pos;
      i++;
    }
    setbit_pos++;
    pop>>=1;
  }

  for (int i=0;i<n;i++){
    //fprintf(stderr,"The mutation_array is %d\n",individual_muts_array[i]);
  }
  return individual_muts_array;
  free(individual_muts_array);
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

double* sum_array(double* array, int size){
  float previous_sum=0;
  float after_sum=0;
  double* sum_array1=(double*)calloc(size,sizeof(double));

  for (int i=0; i<size; i++){
    sum_array1[i]=previous_sum;
    previous_sum+=array[i];
  }

  for (int i=size-1; i>=0; i--){
    sum_array1[i]+=after_sum;
    after_sum+=array[i];
  }

  for (int i=0; i<size; i++){
    sum_array1[i]=2*array[i]+sum_array1[i];
    //fprintf(stderr, "%f\n", sum_array[i]);
  }
  
  return sum_array1;
  free(sum_array1);
}

double find_min(double* array, int size){
  double min_elem;
  min_elem=array[0];
  for (int i=0; i<size; i++){
    if (min_elem>array[i]){
      min_elem=array[i];
    }
  }
  //fprintf(stderr,"The minimum element is %f\n",min_elem);

  return min_elem;
}

//double find_max(double* array, int size){
//double max_elem;
//max_elem=array[0];
//for (int i=0; i<size; i++){
//  if (max_elem<array[i]){
//    max_elem=array[i];
//  }
//}
//return max_elem;
//}


double weight_selection_position1(double* array,gsl_rng* r,int size,int* array1){
  
  double* new_weight_position_array=(double*)calloc(intsize,sizeof(double));
  double* norm_new_weight_position_array=(double*)calloc(intsize,sizeof(double));
  double com_weight_position;

  new_weight_position_array=sum_array(array,intsize);


  double min_fitness=find_min(new_weight_position_array,intsize);

  for (int i=0; i<intsize; i++){
    norm_new_weight_position_array[i]=(float)new_weight_position_array[i]-(float)min_fitness+3;
    //fprintf(stderr,"The norm_new_weight for position %d is %f\n",i,norm_new_weight_position_array[i]);
  }
  
  for (int i=0; i<size; i++){
    com_weight_position+=norm_new_weight_position_array[array1[i]];
  }


  double x=com_weight_position;

  free(new_weight_position_array);
  free(norm_new_weight_position_array);
 
  return com_weight_position;
}
  


void com_fitness_func2(double* array,gsl_rng* r,int N, double* fitness, double* com_fitness, typ *pop, int selected_pos){
  // int k=rand()%64 +1;                                                                                                                                                                                    
  //printf("The positions that matter are %d\n",k);                                                                                                                                                          
  typ fitfactor=fitfactors_func(selected_pos);
  //int* each_individual_mutations;
  for (int i=0;i<N;i++){
    fprintf(stderr,"%llu %llu %llu\n",pop[i],fitfactor,pop[i] &fitfactor);
    if ((typ)(pop[i] & fitfactor) == (typ)0){
      fprintf(stderr,"SMALL FIT\n");
      //fitness[i]=find_min(initial_fitness_array,intsize)-1;
      fitness[i]=2.;
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
  //for (int i=0; i<N; i++){
  //fitness[i]=fitness[i]/com_fitness[N-1];
  //fprintf(stderr,"The new fitness for %d is %f\n",i,fitness[i]);
  //}

  global_fitfactor = fitfactor;
}

double* softmax(double *input, int input_len){
  // assert (input != NULL);                                                                                                                                                                                
  // assert (input_len != 0);                                                                                                                                                                               
    int i;
    double m;
    /* Find maximum value from input array */
    m = input[0];
    for (i = 1; i < input_len; i++) {
        if (input[i] > m) {
            m = input[i];
        }
    }

    double sum = 0;
    for (i = 0; i < input_len; i++) {
        sum += exp(input[i]-m);
    }

    for (i = 0; i < input_len; i++) {
        input[i] = exp(input[i] - m - log(sum));
        //fprintf(stderr,"The new element is %f\n",input[i]);

    }
    return input;
}


  

void selection_func(double* array,gsl_rng* r,int N,typ* pop1, typ* pop2, double *fitness, double* com_fitness,int selected_pos){

  typ k;
  typ* parent_array=(typ*)calloc(N,sizeof(typ));
  //com_fitness_func(N, fitness, com_fitness);                                                                                                                                                             \
                                                                                                                                                                                                            
  //com_fitness_func1(N, fitness, com_fitness, weight_selection_position(), pop1,2000);                                                                                                                    \
                                                                                                                                                                                                            
  com_fitness_func2(array,r,N,fitness,com_fitness,pop1,selected_pos);

  fprintf(stderr,"Total fitness is %f\n", com_fitness[N-1]);

  for (int i=0;i<N;i++){
    double r=(double)rand()/RAND_MAX * com_fitness[N-1];
    //double r=(double)rand() / ((double)RAND_MAX + 1) * com_fitness[N1-1];
      for (k=0;k<N;k++){
        
	if (com_fitness[k]>r) break;
      }

    typ parent = k;
    assert(parent<N);
    pop2[i] = pop1[parent];
    fprintf(stderr,"The individual %d of next population has as parent the individual %llu of the current population\n",i, parent);
    parent_array[i]=parent;
  }

  fitness=softmax(fitness,N);

  for (int i=0; i<N; i++){
    //fitness[i]=fitness[i]/com_fitness[N-1];
    fprintf(stderr,"The new fitness for %d is %f\n",i,fitness[i]);
  }

  //double rand_num=(((double)rand() / (double)RAND_MAX) * com_fitness[N1-1]);                                                                                                                              
  //int sum=0;                                                                                                                                                                                              

  //unsigned int distinct=distinct_elements(parent_array,N2);
  //long double growth_rate=(long double)distinct/(long double)N1;
  //fprintf(stderr,"The number of proliferative individuals(parents) are %d\n",distinct);
  //fprintf(stderr,"The growth rate id %Lf\n",growth_rate);
  //global_growth_rate=growth_rate;                                                                                                                                                                         


  free(parent_array);
}

  

int update_state(int* state, int size, double evolutionary_probability, double* fitness, double* probability, int generation, long double growth_rate, long double decay_rate, typ capacity){
  int count=0;
  for (int i=0; i<size; i++){
    if (state[i-1]==0 && state[i+1]==0){
      probability[i]=fitness[i];
    }
    else if (state[i-1]==1 && state[i+1]==1){
      probability[i]=2*fitness[i];
    }
    else {
      probability[i]=fitness[i]*(1-(Gompertz_discrete_model(1,generation,growth_rate,decay_rate)/capacity));
    }
    if (state[i]==1 || probability[i]<=evolutionary_probability){
      state[i]=state[i];
    }
    else {
      state[i]=state[i]+1;
    }
    if (state[i]==1){
      count++;
    }
    fprintf(stderr,"The new state array for individual %d is %d\n",i,state[i]);
  }
  FILE* f5=fopen("Number of mutants.txt","a");
  fprintf(stderr,"The total number of 1_states is %d\n",count);
  fprintf(f5, "%d\t",count);
  fclose(f5);
  return count;
}

void find_1_states_positions(int* state,int size, int num_of_1s, typ* pop, int total_num){
  int* states_1_array=(int*)calloc(num_of_1s,sizeof(int));
  typ* mutant_genots_array=(typ*)calloc(num_of_1s,sizeof(typ));
  int index=0;
  int j=0;
  for (int i=0; i<size; i++){
    if (state[i]==1){
      states_1_array[j]=index;
      j++;
	}
    index++;
  }
  for (int i=0; i<num_of_1s; i++){
    fprintf(stderr,"The 1 is in position %d\n",states_1_array[i]);
    mutant_genots_array[i]=pop[states_1_array[i]];
    fprintf(stderr,"The mutant genotype is %llu\n",mutant_genots_array[i]);
  }

  freq_func_mutant(mutant_genots_array,num_of_1s);

  setbits_positions_mutant(num_of_1s, mutant_genots_array, total_num);

  free(mutant_genots_array);
  //return mutant_genots_array;
      
}
  



double* weight_position_array_func(gsl_rng*r,int size){

  double* weight_position_array1=(double*)calloc(size,sizeof(double));

  for (int i=0;i<size;i++){
    weight_position_array1[i]=gsl_ran_gaussian(r,1);
    // weight_position_array[i]=(float)rand()/(float)RAND_MAX*10-11;                                                                                                                                       \
                                                                                                                                                                                                            
  }
  return weight_position_array1;
  free(weight_position_array1);
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
      //printf("Mutation_array is %d\n",position_of_mutations_array[j]);                                                                                                                                   \
                                                                                                                                                                                                            
      //mutate(total_pop_pos, pop_prev, i);                                                                                                                                                                \
                                                                                                                                                                                                            

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
      // char* s=bitfunc(pop_cur[k]); /                                                                                                                                                                    \
                                                                                                                                                                                                            
      // printf("The genome (without the possible upcoming mutation) of individual %d of the current population is %s\n",k,s); /                                                                           \
                                                                                                                                                                                                            
    }
}



void fixed_positions(int N, typ* pop){

  for (int i=0;i<N;i++){
    int setBits = countSetBits(pop[i]);
    int commonBits = countSetBits(pop[i] & global_fitfactor); // οι κοινες θεσεις ( άσσοι) μεταξύ των άσσεων του fitfactor και του μεταλλαγμένου γονιδιώματος).                                             
    fprintf(stderr,"Individual:%d\t Number of positions that have been fixed (frome those that matter):%d\t Number of positions that have contibuted (from those that matter):%d\n", i, setBits, commonBits\
\
);

  }
}


int main(int argc, char* argv[]){
  global_fitfactor = 0;
  int N=atoi(argv[1]);
  int generations=atoi(argv[2]);
  int initialN=N; // για Ν σταθερό έτσι όπως έχουμε φτιάξει τον κώδικα ακολουθεί μοντέλο Wright-Fisher.
  int states_1_num;
  

  typ* pop_prev=(typ*)calloc(N,sizeof(typ));
  typ* pop_cur=(typ*)calloc(N,sizeof(typ));
  int* state_array=(int*)calloc(N,sizeof(int));
  double* fitness=(double*)calloc(N,sizeof(double));
  double* com_fitness=(double*)calloc(N,sizeof(double));
  double* probability_array=(double*)calloc(N,sizeof(double));

  double mu=(double)pow((double)10,(double)-2); // for cancer cells, while for normal cells is 10^-8 or 10^-6;                                                                                              
  fprintf(stderr,"The mutation rate is %f\n",mu);
  long double growth_rate=log(2);
  int selected_pos=rand()%64+1;
  fprintf(stderr,"The positions that matter are %d\n",selected_pos);


  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  gsl_rng_set(r, time(NULL));

  double evolutionary_probability=gsl_ran_flat(r,0,1);
  //double evolutionary_probability=0.01;
  fprintf(stderr,"The evolutionary probability is %f\n",evolutionary_probability);

  double* weight_position_array=weight_position_array_func(r,intsize);

  for (int i=0; i<generations; i++){
    if (i>300){
      mu=0;
    }

    long long int total_pop_pos=N*intsize;


    int numberOfMutations=binomial(r,mu,N);
    fprintf(stderr,"The %d generation which has undergone %d  mutations\n",i,numberOfMutations);
    //fitness_chars_func(N);                                                                                                                                                                                

    print_genome_without(N,pop_prev,i);


    mutation_process(numberOfMutations,total_pop_pos,pop_prev,i);


    int total_setbits_num=print_genome_with(N,pop_prev);


    fprintf(stderr,"The genome frequencies are:\n");


    FILE* f1=fopen("Frequencies_CA.txt","a");
    fprintf(f1,"Generation %d\t",i);
    fclose(f1);

    freq_func(pop_prev,N);// τρέχει σωστά και βγάζει συχνότητα γοιδιωμάτων μετα την μετάλλαξη.                                                                                                          \
                                                                                                                                                                                                            

    f1=fopen("Frequencies_CA.txt","a");
    fprintf(f1,"\n");
    fclose(f1);

    if (i==0){
      com_fitness_func2(weight_position_array,r,N,fitness,com_fitness,pop_prev,selected_pos);
      int count=0;
      for (int j=0; j<N; j++){
	if (fitness[j]!=2){
	  state_array[j]=1;
	  count++;
	}
	fprintf(stderr,"The initial state array in generation %d for individual %d is %d\n",i,j,state_array[j]);
      }

      FILE* f5=fopen("Number of mutants.txt","a");
      fprintf(f5,"Generation %d\t",i);
      fclose(f5);

      f5=fopen("Number of mutants.txt","a");
      fprintf(stderr,"The initial number of 1_states are %d\n",count);
      fprintf(f5, "%d\t",count);
      fclose(f5);

      f5=fopen("Number of mutants.txt","a");
      fprintf(f5,"\n");
      fclose(f5);
    
      
      FILE* f3=fopen("Frequencies_CA_mutant.txt","a");
      fprintf(f3,"Generation %d\t",i);
      fclose(f3);

      FILE* f4=fopen("Frequencies_mutations_CA_mutant.txt","a");
      fprintf(f4,"Generation %d\t",i);
      fclose(f4);

      
      find_1_states_positions(state_array, N, count,pop_prev, total_setbits_num);

      f3=fopen("Frequencies_CA_mutant.txt","a");
      fprintf(f3,"\n");
      fclose(f3);

      f4=fopen("Frequencies_mutations_CA_mutant.txt","a");
      fprintf(f4,"\n");
      fclose(f4);
    }

    selection_func(weight_position_array,r,N, pop_prev, pop_cur, fitness, com_fitness,selected_pos);

    if (i>0){

      FILE* f5=fopen("Number of mutants.txt","a");
      fprintf(f5,"Generation %d\t",i);
      fclose(f5);


      states_1_num=update_state(state_array, N, evolutionary_probability, fitness, probability_array, i, growth_rate,  0.2, N);

      f5=fopen("Number of mutants.txt","a");
      fprintf(f5,"\n");
      fclose(f5);
    

    fprintf(stderr,"The mutant genome frequencies are:\n");

    FILE* f3=fopen("Frequencies_CA_mutant.txt","a");
    fprintf(f3,"Generation %d\t",i);
    fclose(f3);

    FILE* f4=fopen("Frequencies_mutations_CA_mutant.txt","a");
    fprintf(f4,"Generation %d\t",i);
    fclose(f4);

    find_1_states_positions(state_array, N, states_1_num,pop_prev, total_setbits_num);

    f3=fopen("Frequencies_CA_mutant.txt","a");
    fprintf(f3,"\n");
    fclose(f3);

    f4=fopen("Frequencies_mutations_CA_mutant.txt","a");
    fprintf(f4,"\n");
    fclose(f4);

    }


    
    FILE* f2=fopen("Frequencies_mutations_CA.txt","a");
    fprintf(f2,"Generation %d\t",i);
    fclose(f2);


    setbits_positions(N,pop_prev,total_setbits_num);

    f2=fopen("Frequencies_mutations_CA.txt","a");
    fprintf(f2,"\n");
    fclose(f2);

    population_realloc(N, pop_prev, pop_cur);

  }

  fixed_positions(N,pop_cur);

  free(pop_prev);
  free(pop_cur);
  free(fitness);
  free(com_fitness);
  gsl_rng_free(r);
  free(weight_position_array);
  free(state_array);
  free(probability_array);

  return 0;
}






    

