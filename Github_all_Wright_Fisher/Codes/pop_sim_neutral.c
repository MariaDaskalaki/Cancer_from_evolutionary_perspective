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
	FILE *f1=fopen("Frequencies.txt","a");                                                                                                                                                             
	fprintf(f1, "%llu:%Lf\t",array[i],(long double)freq[i]/(long double)size);
	fclose(f1);
	free(s);
    
     
    }
  }


  //sorting(freq,size);                                                                                                                                                                                     
  free(freq);
}



typ exponential_linear(int initial, int t, long double growth_rate, double linear_const){
  long double t_met=((double)1/growth_rate)*log10((double)linear_const/(double)(growth_rate*initial));
  fprintf(stderr,"The change time is %Lf\n",t_met);
  //printf("%Lf\n",t_met);                                                                                                                                                                                 \
                                                                                                                                                                                                            
  if (t<=t_met){
    typ newSize=initial*exp(growth_rate*t);
    //printf("%llu\n",newSize);                                                                                                                                                                             
    return newSize;
  }
  else{
    typ newSize=linear_const*t+initial;
    //printf("%llu\n",newSize);                                                                                                                                                                             
    return newSize;
  }

} // τρέχει σωστά η συνάρτηση του εκθετικού, αλλά πρεπει στον κανονικό κώδικα να βρω το growth_rate και μία ρεαλιστική γραμμική σταθερά.



typ generalized_logistic(int initial,int t, long double capacity, long double growth_rate, double n_const){
  double numerator=initial*capacity;
  //printf("%f\n",numerator);                                                                                                                                                                               
  double denominator_power=(double)1/(double)n_const;
  //printf("%f\n",denominator_power);                                                                                                                                                                       
  double denominator=(double)pow((double)pow(initial,n_const)+((double)pow(capacity,n_const)-(double)pow(initial,n_const))*exp(-growth_rate*n_const*t),denominator_power);
  //printf("%f\n",denominator);                                                                                                                                                                             
  typ newSize=(double)numerator/(double)denominator;
  //printf("%llu\n",newSize);                                                                                                                                                                               
  return newSize;
} /* τρέχει σωστά */



typ Gompertz_model(int initial,int t, long double growth_rate, long double decay_rate){
  //long double c=(long double)growth_rate/(long double)decay_rate;
  double c=27.631; //log(10^12)
  typ newSize=initial*(exp(c*(1-exp(-decay_rate*t))));
  //printf("%llu\n",newSize);
  if (newSize > 100000){
    newSize=100000;
  }
  return newSize;
} /* τρέχει σωστά.*/



typ dcc_model(int initial, int t, long double growth_rate, int b_const){
  long double initial_capacity=1;
  long double praxis=(long double)powl(t,((double)5/(double)3));
  // printf("%Lf\n",praxis);                                                                                                                                                                               \
                                                                                                                                                                                                            
  long double capacity=initial_capacity*exp((((long double)(3*b_const*praxis))/(long double)5));
  //printf("%Lf\n",capacity);                                                                                                                                                                              \
                                                                                                                                                                                                            
  long double praxis3=1-exp(growth_rate*t);
  //printf("%Lf\n",praxis3);                                                                                                                                                                               \
                                                                                                                                                                                                            
  long double praxis2=(long double)powl((long double)capacity,(long double)praxis3); // οσο μεγαλώνει η βάση με αρνητικό εκθέτη τείνει στο 0.                                                              \
                                                                                                                                                                                                            
  //printf("%Lf\n",praxis2);                                                                                                                                                                               \
                                                                                                                                                                                                            
  typ newSize=initial*exp(-growth_rate*t)*((long double)powl(capacity,(1-exp(growth_rate*t))));
  //printf("%llu\n",newSize);                                                                                                                                                                              \
                                                                                                                                                                                                            
  return newSize;
} // δεν τρέχει σωστά.



typ von_bertalanffy(int initial, int t , long double growth_rate, long double decay_rate){ // for decay_rate=0 we have power law model
  double gamma=(double)2/(double)3; // "second type growth" von_bertalanffy model                                                                     \
                                                                                                                                                                                                            
  long double exponential=exp(-decay_rate*(1-gamma)*t);
  long double power=(long double)1/(long double)(1-gamma);
  long double praxis=(long double)growth_rate/(long double)decay_rate+((long double)pow(initial,1-gamma)-((long double)growth_rate/(long double)decay_rate))*exponential;
  typ newSize=(long double)powl((long double)praxis,(long double)power);
  //printf("%llu\n",newSize);                                                                                                                                                                              \
                                                                                                                                                                                                            
  return newSize;
} // τρέχει σωστά                                                                                                                                                                                          \


int mutate(typ input, typ *pop, int generation)
{
  
  typ k;
  //  if (generation == 10) // για να συμβει σίγουρα η μετάλλαξη στη θέση 0 στη γενιά 10.                                                                                                                   
  // k = 0;                                                                                                                                                                                                 
  // else                                                                                                                                                                                                   
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
    FILE* f2=fopen("Frequencies_mutations.txt","a");
    fprintf(stderr,"The mutation in position %d has occured %d times\n",j,position_array[j]);
    fprintf(f2,"%d:%f\t",j,(double)position_array[j]/((double)total_num));   
    fclose(f2);
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

float* sum_array(float* array, int size){
  float previous_sum=0;
  float after_sum=0;
  float* sum_array1=(float*)calloc(size,sizeof(float));

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

float find_min(float* array, int size){
  float min_elem;
  min_elem=array[0];
  for (int i=0; i<size; i++){
    if (min_elem>array[i]){
      min_elem=array[i];
    }
  }
  //fprintf(stderr,"The minimum element is %f\n",min_elem);

  return min_elem;
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

void com_fitness_func_neutral(int N, double* fitness, double* com_fitness){
  for (int i=0; i<N; i++){
    fitness[i]=1.;
    com_fitness[0]=fitness[0];
    if (i>0){
      com_fitness[i]=fitness[i]+com_fitness[i-1];
    }
    fprintf(stderr, "NEUTRAL FIT of individual %d is %f\n",i,fitness[i]);
    fprintf(stderr, "Com_fitness of individual %d is %f\n",i,com_fitness[i]);
  }
}


  

  


void com_fitness_func2(float* array,gsl_rng* r,int N, double* fitness, double* com_fitness, typ *pop, int selected_pos){
  // int k=rand()%64 +1;                                                                                                                                                                                    
  //printf("The positions that matter are %d\n",k);

  int* individual_array;
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
      individual_array= individual_setbits_positions(num_of_individuals_muts,pop[i]);
      fprintf(stderr,"LARGER FIT\n");
      fitness[i]= 1.; // υποτίθεται είναι με τις αλληλεπιδράσεις (πολυμεταβλητό γραμμικό μοντέλο), ενώ το weight_selection_position είναι με το απλό γραμμικό μοντέλο όπου κάθε θέση είναι ανεξάρτητη.
    
    }
    fprintf(stderr, "Fitness of %d is  %f\n", i, fitness[i]);
    com_fitness[0]=fitness[0];
    if (i>0){
      com_fitness[i]=fitness[i]+com_fitness[i-1];
      fprintf(stderr, "Com_fitness of %d is %f\n",i,com_fitness[i]);
    }
  }
  global_fitfactor = fitfactor;
  free(individual_array);
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


int distinct_elements(typ* array, int size) {
    int count = 1;

    // Pick all elements one by one                                                                                                                                                                        \
                                                                                                                                                                                                            
    for (int i = 1; i < size; i++){
        int j = 0;
        for (j = 0; j < i; j++){
            if (array[i] == array[j])
                break;
        }

        if (i == j){
            count++;
        }
    }
    return count;
}


void selection_func(float* array,gsl_rng* r,int N1, int N2,typ* pop1, typ* pop2, double *fitness, double* com_fitness,int selected_pos){

  typ k;
  typ* parent_array=(typ*)calloc(N2,sizeof(typ));                                                                                                                                                         
  //com_fitness_func(N, fitness, com_fitness);                                                                                                                                                              
  //com_fitness_func1(N, fitness, com_fitness, weight_selection_position(), pop1,2000);                                                                                                                     
  //com_fitness_func2(array,r,N1,fitness,com_fitness,pop1,selected_pos);

  com_fitness_func_neutral(N1,fitness,com_fitness);

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
      //population_array[i]=rand()%20 + initialN; //rand() % N + A will give you a value in the range [A, A + N], not [A, N]
      population_array[i]=exponential_linear(initialN,i,growth_rate,2); // for 10000 exp until 3 generation.
      //population_array[i]=generalized_logistic(initialN,i,(double)pow((double)10,(double)4),growth_rate,0.25); //capacity~10^12 check paper Estimating tumor growth rates in vivo, for the parameters.
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
      //individual_setbits_positions(num_of_setbits,pop[k]);
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


  
  


int main(int argc, char* argv[]){
  global_fitfactor = 0;
  int N=atoi(argv[1]);
  int generations=atoi(argv[2]);
  int initialN=N; // για Ν σταθερό έτσι όπως έχουμε φτιάξει τον κώδικα ακολουθεί μοντέλο Wright-Fisher.                                                                                                     
  //int generations;                                                                                                                                                                                        
  // printf("Enter the number N of cells\n");                                                                                                                                                               
  //scanf("%d",&N);                                                                                                                                                                                         
  //printf("Enter the number of generations\n");                                                                                                                                                            
  //scanf("%d",&generations);                                                                                                                                                                               
  typ* pop_prev=(typ*)calloc(N,sizeof(typ));
  typ* pop_cur=(typ*)calloc(N,sizeof(typ));
  double* fitness=(double*)calloc(N,sizeof(double));
  double* com_fitness=(double*)calloc(N,sizeof(double));
  int* population_array;

  double mu=(double)pow((double)10,(double)-4); // for cancer cells, while for normal cells is 10^-8 or 10^-6;
  fprintf(stderr,"The mutation rate is %f\n",mu);
  long double growth_rate=log(2);
  int selected_pos=rand()%64+1;
  fprintf(stderr,"The positions that matter are %d\n",selected_pos);
  //int genChange=1;                                                                                                                                                                                        
  //int newSize = initialN = N;                                                                                                                                                                             
  int size;
  int next_size;



  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  gsl_rng_set(r, time(NULL));

  //float* weight_position_array=weight_position_array_func(r,intsize);

  population_array=population_change(initialN,generations,growth_rate,r);



  for (int i=0;i<generations;i++){
    if (i>300){ // αν διακοπεί ο μεταλλακτικός ρυθμός περιμένουμε εγκαθίδρυση                                                                                                                               
      mu=0;
    }



                                                                                                                                                                                        

      size=population_array[i];
      next_size=population_array[i+1];

      if (next_size>size){                                                                                                                                                                                   
        pop_cur=(typ*)realloc(pop_cur,sizeof(typ)*next_size);
        pop_prev=(typ*)realloc(pop_prev,sizeof(typ)*next_size);
        fitness=(double*)realloc(fitness,sizeof(typ)*next_size);
        com_fitness=(double*)realloc(com_fitness,sizeof(typ)*next_size);
        assert(pop_cur!=NULL);
        assert(pop_prev!=NULL);
        assert(fitness!=NULL);
        assert(com_fitness!=NULL);
      }

    long long int total_pop_pos=size*intsize;


    int numberOfMutations=binomial(r,mu,size);
    fprintf(stderr,"The %d generation which has undergone %d  mutations\n",i,numberOfMutations);
    //fitness_chars_func(N);

    print_genome_without(size,pop_prev,i);


    mutation_process(numberOfMutations,total_pop_pos,pop_prev,i);


    int total_setbits_num=print_genome_with(size,pop_prev);



    fprintf(stderr,"The genome frequencies are:\n");


    FILE* f1=fopen("Frequencies.txt","a");
    fprintf(f1,"Generation %d\t",i);
    fclose(f1);

    freq_func(pop_prev,size);// τρέχει σωστά και βγάζει συχνότητα γοιδιωμάτων μετα την μετάλλαξη.                                                                                                           

    f1=fopen("Frequencies.txt","a");
    fprintf(f1,"\n");
    fclose(f1);


    selection_func(weight_position_array,r,size,next_size, pop_prev, pop_cur, fitness, com_fitness,selected_pos);

    FILE* f2=fopen("Frequencies_mutations.txt","a");
    fprintf(f2,"Generation %d\t",i);
    fclose(f2);


    setbits_positions(size,pop_prev,total_setbits_num);

    f2=fopen("Frequencies_mutations.txt","a");
    fprintf(f2,"\n");
    fclose(f2);

    population_realloc(next_size, pop_prev, pop_cur);

  }

  fixed_positions(size,pop_cur);

  size=next_size;

  for (int i=0;i<generations;i++){
    //fprintf(stderr,"%d\n",population_array[i]);
    FILE* f3=fopen("Population_size.txt","a");
    fprintf(f3,"%d\n",population_array[i]);
    fclose(f3);
  }

  free(population_array);
  free(fitness);
  free(com_fitness);
  free(pop_prev);
  free(pop_cur);
  gsl_rng_free(r);
  //free(weight_position_array);

  return 0;
}
