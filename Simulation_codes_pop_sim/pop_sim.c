#include "pop_sim_libraries_CA.h"


typ global_fitfactor;
long double global_growth_rate;
double* global_fitness_site_array;

#define MAXVOL 10000

time_t t;


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

int freq_func(typ* array, int size, FILE* f1, int generation){
  typ* freq=(typ*)calloc(size,sizeof(typ));
  int fixation_param=0;
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
	//fprintf(stderr, "%llu:The genome %s has frequency %f which means that occurs %llu times and population size %d\n",array[i], s, (double)freq[i]/(double)size,freq[i],size);
	fprintf(f1, "%llu:%Lf\t",array[i],(long double)freq[i]/(long double)size);
	if (generation>2*size){
	  if ((long double)freq[i]/(long double)size==1){
	    fixation_param=1;
	  }
	  else{
	    fixation_param=0;
	  }
	}
	free(s);


    }
  }
  free(freq);
  return fixation_param;
}


void freq_func_mutant(typ* array, int size, FILE* f3){
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
	//fprintf(stderr, "%llu:The mutant genome %s has frequency %f which means that occurs %llu times and CAPopsize is %d\n",array[i], s, (double)freq[i]/(double)size,freq[i],size);
	fprintf(f3, "%llu:%Lf\t",array[i],(long double)freq[i]/(long double)size);
	free(s);


    }
  }
  free(freq);
}

long double Gompertz_discrete_model(int initial_cancer_cells, int t, long double growth_rate, long double decay_rate,double radius){ /* it is used for the calculation of the population size of mutant cells */
  long double c;
  c=(double)growth_rate/(double)decay_rate;
  double radius_power=pow(radius,(double)3);
  double volume=(double)(4*3.14*radius_power)/(double)3;
  double initial_tum_volume=initial_cancer_cells*volume;
  //double newSize2=(growth_rate*initial_tum_volume*1.0)*(exp(-decay_rate*t))*(exp(c*(1-exp(-decay_rate*t))));
  //double newSize2=growth_rate*initial_tum_volume*exp(-decay_rate*t)*exp(c*(1-exp(-decay_rate*t)))+growth_rate*initial_tum_volume*1.0*exp(c*(1-exp(-decay_rate*t)));
  double newSize2=growth_rate*initial_tum_volume*exp(-decay_rate*t)*exp(c*(1-exp(-decay_rate*t)));

  if(newSize2<0){
    newSize2=1;
  }

  //assert(newSize2 > 0);

  if (newSize2 > MAXVOL){
    newSize2=MAXVOL;
  }

  return newSize2;
}

int mutate(typ input, typ *pop, int generation){

  typ k;

  k= rand() % (input);

  typ individual = k/intsize;
  int position = k % intsize;

  typ mask = (typ)1 << position;
  pop[individual] = pop[individual] | mask; // ^ xor bitwise operator (if a mutation rehits it can convert 1 to 0 as 1 and 1 returns 0) and | or bitwise operator ( if a mutation rehits it remains 1 as 1 and 1 returns 1).
  //fprintf(stderr, "k:%llu, %llu\n", k, mask);
  //fprintf(stderr, "The mutation has happened in the %llu individual of the current population and in POSITION %d\n",individual,position);


  return position;
}

int binomial(gsl_rng * r, double mu, int N){
  long long int total_pop_pos=N*intsize;
  int total_num_mut = gsl_ran_binomial (r, mu, total_pop_pos);

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


void setbits_positions(int N, typ* pop,int total_num,FILE* f2){ //finds the 1 in the genome number ( positions that are mutated) 
  int* position_array=(int*)calloc(64,sizeof(int));
  for (int i=0;i<N;i++){
    unsigned int setbit_pos=0;
    while(pop[i]){
      if((pop[i]&1)!=0){
        //fprintf(stderr, "The setbit position is %d\n",setbit_pos);
        position_array[setbit_pos]+=1;
      }
      setbit_pos++;
      pop[i]>>=1;
    }
  }
  for (int j=0;j<64;j++){
    //fprintf(stderr,"The mutation in position %d has occured %d times\n",j,position_array[j]);
    //fprintf(f2,"%d:%f\t",j,(double)position_array[j]/((double)total_num));
    fprintf(f2,"%d:%f\t",j,(double)position_array[j]/(double)N);
  }
  free(position_array);
}


void setbits_positions_mutant(int N, typ* pop,int total_num, FILE* f4){ // finds the 1 in the mutants number ( positions in mutants that are mutated) 
  int* position_array=(int*)calloc(64,sizeof(int));
  for (int i=0;i<N;i++){
    unsigned int setbit_pos=0;
    while(pop[i]){
      if((pop[i]&1)!=0){
        //fprintf(stderr, "The setbit position is %d\n",setbit_pos);
        position_array[setbit_pos]+=1;
      }
      setbit_pos++;
      pop[i]>>=1;
    }
  }
  for (int j=0;j<64;j++){
    //fprintf(stderr,"The mutant mutation in position %d has occured %d times\n",j,position_array[j]);
    //fprintf(f4,"%d:%f\t",j,(double)position_array[j]/((double)total_num));
    fprintf(f4,"%d:%f\t",j,(double)position_array[j]/(double)N);
  }
  free(position_array);
}


int* individual_setbits_positions(int n, typ pop){
  int* individual_muts_array=(int*)calloc(n,sizeof(int));
  unsigned int setbit_pos=0;
  int i=0;
  while (pop){
    if (pop&1==1){
      individual_muts_array[i]=setbit_pos;
      i++;
    }
    setbit_pos++;
    pop>>=1;
  }
  return individual_muts_array;
}


typ fitfactors_func(unsigned k){
  typ* fitfactor_array=(typ*)calloc(intsize,sizeof(typ));
  typ* fitfactor_array_final=(typ*)calloc(intsize,sizeof(typ));
  typ fitfactor=0;

  for (unsigned i=0;i<intsize;i++){
    fitfactor_array[i]=(typ)1 << i;
    fitfactor=fitfactor | (fitfactor_array[i]);
    fitfactor_array_final[i]=fitfactor;
  }

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
  }

  return sum_array1;
}

double find_min(double* array, int size){
  double min_elem;
  min_elem=array[0];
  for (int i=0; i<size; i++){
    if (min_elem>array[i]){
      min_elem=array[i];
    }
  }

  return min_elem;
}

double weight_selection_position1(double* array,gsl_rng* r,int size,int* array1){

  //double* new_weight_position_array=(double*)calloc(intsize,sizeof(double));
  double* norm_new_weight_position_array=(double*)calloc(intsize,sizeof(double));
  double com_weight_position=0.;
  //double* com_weight_position=(double*)calloc(intsize,sizeof(double));

  double* new_weight_position_array=sum_array(array,intsize);


  double min_fitness=find_min(new_weight_position_array,intsize);

  for (int i=0; i<intsize; i++){
    //norm_new_weight_position_array[i]=(float)new_weight_position_array[i]-(float)min_fitness+3;
    norm_new_weight_position_array[i]=(float)new_weight_position_array[i]-(float)min_fitness;
    //norm_new_weight_position_array[i]=1;
    assert(norm_new_weight_position_array[i] >= 0);
  }

  for (int i=0; i<size; i++){
    com_weight_position+=norm_new_weight_position_array[array1[i]];
    //fprintf(stderr,"heyyyy---%f\n",com_weight_position);
  }


  double x=com_weight_position;

  free(new_weight_position_array);
  free(norm_new_weight_position_array);

  //return norm_new_weight_position_array;

  return com_weight_position;
}

typ* sorting_1(typ* index_array, typ* array, int N){
  int temp;
  typ* new_index=(typ*)calloc(sizeof(typ),N);

  for(int i = 0 ; i < N ; i++)
    {
        for(int j = i+1 ; j < N ; j++)
        {
            if(index_array[i] > index_array[j])
            {
                // Swap elements in first array
                temp = index_array[j] ;
                index_array[j] = index_array[i] ;
                index_array[i] = temp ;


                // Swap elements in second array here itself
                temp = array[j] ;
                array[j] = array[i] ;
                array[i] = temp ;
            }
        }
	new_index[array[i]]=i;
    }
  for (int i=0; i<N; i++){
    //fprintf(stderr,"The sorted_index_array is %llu\n",index_array[i]);
  }
  for (int i=0; i<N; i++){
    //fprintf(stderr,"The sorted_array is %llu\n",array[i]);
  }
  for (int i=0; i<N; i++){
    //fprintf(stderr,"The %d individual that is the offspring of the individual with sorted_site %llu is now at site %llu\n",i,index_array[i],new_index[i]);
  }
  return new_index;
}


double* fitness_site_func_linear(int N){
  double* fitness_site_array=(double*)calloc(N,sizeof(double));

  FILE* f8=fopen("Linear_fitness_site_distribution.txt", "w");

  for (int i=0; i<=N/2; i++){
    fitness_site_array[i]=(-((double)18/(double)N)*i)+10;
  }
  for (int i=N/2+1; i<N; i++){
    //fitness_site_array[i]= (18.0/N) * i - 8.0; //(((double)18/(double)(N-2))*i)+(1-((double)(18*N)/(double)(2*(N-2))));
    //fitness_site_array[i]=((double)18/(double)(2-N))*i+10+18*((double)(N-1)/(double)(2-N));
    fitness_site_array[i]=(((double)18/(double)(N-2))*i)+(1-((double)(18*N)/(double)(2*(N-2))));
  }
  for (int i=0; i<N; i++){
    //fprintf(stderr,"The fitness_site_array is %f\n",fitness_site_array[i]);
    fprintf(f8,"%d\t%f\n",i,fitness_site_array[i]);
  }

  fclose(f8);

  return fitness_site_array;
}

double* fitness_site_func_alter(int N){
  double* fitness_site_array=(double*)calloc(N,sizeof(double));

  FILE* f8=fopen("Alternative_fitness_site_distribution.txt","w");

  for (int i=0; i<N; i++){
    fitness_site_array[i]=sin(i)*sin(i);
    fprintf(f8,"%d\t%f\n",i,fitness_site_array[i]);
  }
  fclose(f8);

  return(fitness_site_array);
}

double* alter_fitness(int N, typ* index_array){
  double* fitness_individuals_site_array=(double*)calloc(N,sizeof(double));

  for (int i=0; i<N; i++){
    fitness_individuals_site_array[i]=global_fitness_site_array[index_array[i]];
  }
  return fitness_individuals_site_array;
}

double* linear_fitness(int N, typ* index_array){
  double* fitness_individuals_site_array=(double*)calloc(N,sizeof(double));

  for (int i=0; i<N; i++){
    fitness_individuals_site_array[i]=global_fitness_site_array[index_array[i]];
    //fprintf(stderr,"The individual %d of the current population with position %llu has fitness %f\n",i,index_array[i],fitness_individuals_site_array[i]);
  }

  return fitness_individuals_site_array;
}

double* fitness_site_func5(int N,gsl_rng* r){ /* spatial fitness values from random beta bimodal distribution */
  double* fitness_site_array=(double*)calloc(N,sizeof(double));
  int* array=(int*)calloc(N,sizeof(int));
  double* linearized_array=(double*)calloc(N,sizeof(double));
  int input_start=0;
  int input_end=N-1;
  double output_start=0.01;
  double output_end=0.99;
  double slope = 1.0 * (output_end - output_start) / (input_end - input_start);

  FILE* f8=fopen("Beta_distribution.txt","w");

  for (int i=0; i<N; i++){
    array[i]=i;
  }

  for (int i=0; i<N; i++){
    linearized_array[i]=output_start + slope * (array[i] - input_start);
    fitness_site_array[i]=gsl_ran_beta_pdf(linearized_array[i],0.5,0.4);
    //fprintf(stderr,"The fitness in position %d is %f\n",i,fitness_site_array[i]);
    fprintf(f8, "%d\t%f\n",i,fitness_site_array[i]);
  }

  fclose(f8);

  free(array);
  free(linearized_array);
  return fitness_site_array;

}

double* beta_fitness(int N, typ* index_array){
  double* fitness_individuals_site_array=(double*)calloc(N,sizeof(double));

  for (int i=0; i<N; i++){
    fitness_individuals_site_array[i]=global_fitness_site_array[index_array[i]];
    //fprintf(stderr,"The individual %d of the current population with position %llu has fitness %f\n",i,index_array[i],fitness_individuals_site_array[i]);
  }

  return fitness_individuals_site_array;
}

double* fitness_site_func4(int N, typ* index_array){ /* spatial fitness values from linear function with minimum 0 in the middle and maximum 2 at the edges */
  double* fitness_site_array=(double*)calloc(N,sizeof(double));
  double* fitness_individuals_site_array=(double*)calloc(N,sizeof(double));

  for (int i=0; i<N/2; i++){
    fitness_site_array[i]=(-((double)(4*i)/(double)N))+2;
  }
  for (int i=N/2; i=N; i++){
    fitness_site_array[i]=((double)(4*i)/(double)(N-2))-((double)(4*N)/(double)(2*(N-2)));
  }
  for (int i=0; i<N; i++){
    fitness_individuals_site_array[i]=fitness_site_array[index_array[i]];
  }

  free(fitness_site_array);
  return fitness_individuals_site_array;
}

double* fitness_site_func_constant(int N, double c){
  double* fitness_site_array=(double*)calloc(N,sizeof(double));

  FILE* f8=fopen("Constant_fitness_site_distribution.txt", "w");

  for (int i=0; i<N; i++){
    fitness_site_array[i]=c;
    fprintf(f8, "%d\t%f\n", i, fitness_site_array[i]);
  }

  fclose(f8);

  return fitness_site_array;
}


double* constant_fitness(int N, typ* index_array){
  double* fitness_individuals_site_array=(double*)calloc(N,sizeof(double));

  for (int i=0; i<N; i++){
    fitness_individuals_site_array[i]=global_fitness_site_array[index_array[i]];
    //fprintf(stderr,"The individual %d of the current population with position %llu has fitness %f\n",i,index_array[i],fitness_individuals_site_array[i]);
  }

  return fitness_individuals_site_array;
}


double* fitness_site_func2(int N, typ* index_array){ /* spatial fitness values from parabolic function with minimum 0 in the middle */
  double* fitness_site_array=(double*)calloc(N,sizeof(double));
  double* fitness_individuals_site_array=(double*)calloc(N,sizeof(double));
  int half=(double)N/(double)2;

  for (int i=0; i<N; i++){
    fitness_site_array[i]=(double)pow((double)(i-half),(double)2);
  }
  for (int i=0; i<N; i++){
    fitness_individuals_site_array[i]=fitness_site_array[index_array[i]];
  }

  free(fitness_site_array);
  return fitness_individuals_site_array;
}

double* fitness_site_func3(int N, typ* index_array, gsl_rng* r){ /* spatial fitness values from random gauss distribution */
  double* fitness_site_array=(double*)calloc(N,sizeof(double));
  double* fitness_individuals_site_array=(double*)calloc(N,sizeof(double));

  for (int i=0; i<N; i++){
    fitness_site_array[i]=gsl_ran_gaussian(r,1);
  }
  for (int i=0; i<N; i++){
    fitness_individuals_site_array[i]=fitness_site_array[index_array[i]];
  }

  free(fitness_site_array);
  return fitness_individuals_site_array;
}


int* fitness_func2(double* array,gsl_rng* r,int N, double* fitness, typ *pop, int selected_pos, typ* index_array, int* individuals_condition,char* user_fitness_site,double constant_param,FILE* f7,char* cell_cancer_model){

  double* fitness_individuals;

  if (!strcmp(user_fitness_site,"linear")){
    fitness_individuals=linear_fitness(N,index_array);
  }

  if(!strcmp(user_fitness_site,"constant")){
    fitness_individuals=constant_fitness(N,index_array);
  }

  if(!strcmp(user_fitness_site,"beta")){
    fitness_individuals=beta_fitness(N,index_array);
  }

  if (!strcmp(user_fitness_site,"alter")){
    fitness_individuals=alter_fitness(N,index_array);
  }

  typ fitfactor=fitfactors_func(selected_pos);
  double fitfactor_thres;
  double random_number;

  if (!strcmp(cell_cancer_model,"neutral")){
    for (int i=0; i<N; i++){
      fitness[i]=fitness_individuals[i];
      individuals_condition[i]=0;
    }
  }
  else{

  if (strcmp(cell_cancer_model,"deterministic")){

    fprintf(f7,"Fitfactor %llu\n",fitfactor);
  }

  if (strcmp(cell_cancer_model,"neutral")){

  for (int i=0;i<N;i++){
    //fprintf(stderr,"%llu %llu %llu\n",pop[i],fitfactor,pop[i] &fitfactor);
    if ((typ)(pop[i] & fitfactor) == (typ)0){
      //fprintf(stderr,"SMALL FIT\n");
      fitness[i]=fitness_individuals[i];
      if(fitness[i] < 0){
	//fprintf(stderr, "fitness of %d is negative: %f\n", i, fitness[i]);
	assert(fitness[i] >= 0);
      }
      //fitness[i]=2.;
      individuals_condition[i]=0;
    }
    else{
      int num_of_individuals_muts=countSetBits(pop[i]);
      int* individual_array= individual_setbits_positions(num_of_individuals_muts,pop[i]);
      fitness[i]=weight_selection_position1(array,r,num_of_individuals_muts,individual_array)+fitness_individuals[i];// υποτίθεται είναι με τις αλληλεπιδράσεις (πολυμεταβλητό γραμμικό μοντέλο), ενώ το weight_selection_position είναι με το απλό γραμμικό μοντέλο όπου κάθε θέση είναι ανεξάρτητη.
      //fprintf(stderr,"cancer----%f+%f\n",weight_selection_position1(array,r,num_of_individuals_muts,individual_array),fitness_individuals[i]);
      if (strcmp(cell_cancer_model,"deterministic")){
	//fitfactor_thres=((double)(pop[i] & fitfactor)/(double)fitfactor)*100;
	fitfactor_thres=(double)(pop[i] & fitfactor)/(double)fitfactor;
	random_number=((double)rand()/(RAND_MAX + 1.0));
	fprintf(f7,"Fitfactor_thres %f and random_number %f\n", fitfactor_thres, random_number);
	if (random_number<fitfactor_thres){
	  //fprintf(stderr,"LARGER FIT and cancerous\n");
	  individuals_condition[i]=1;
	}
	else{
	  //fprintf(stderr,"LARGER FIT but no cancerous\n");
	  individuals_condition[i]=0;
	}
      }
      else{
	//fprintf(stderr,"LARGER FIT and definately cancerous\n");
	individuals_condition[i]=1;
      }
      free(individual_array);

    }

    //fprintf(stderr, "Fitness by site and genotype of %d is  %f and is %d\n", i, fitness[i],individuals_condition[i]);
  }

  global_fitfactor = fitfactor;
  }
  }

  free(fitness_individuals);
  return individuals_condition;
}

double* softmax(double *input, int input_len){
  // assert (input != NULL);
  // assert (input_len != 0);

    double* softmax_array=(double*)calloc(input_len,sizeof(double));
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
        softmax_array[i] = exp(input[i] - m - log(sum));
        //fprintf(stderr,"The new fitness_probs is %f\n",softmax_array[i]);

    }
    return softmax_array;
}


double find_max(double* array, int size){
  double max_elem;
  max_elem=array[0];
  for (int i=0; i<size; i++){
    if (max_elem<array[i]){
      max_elem=array[i];
    }
  }
  return max_elem;
}

double* linear_conversion(double* array, int size){
  int new_min=0;
  int new_max=1;
  double old_min=find_min(array,size);
  double old_max=find_max(array,size);
  double sum=0;

  for (int i=0; i<size; i++){
    array[i]=(( (double)(array[i] - old_min) / (double)(old_max - old_min) ) * (new_max - new_min)) + new_min;
  }
  return array;
}

double* logit(double* array, int size){
  for (int i=0; i<size; i++){
    array[i]=((double)exp(array[i]))/((double)(1 + exp(array[i])));
    //fprintf(stderr,"The logit is %f\n",array[i]);
  }
  return array;
}

typ* selection_func(int N1, int N2 ,typ* pop1, typ* pop2, double *fitness, double* com_fitness, typ* index_array,int* state_due_to_parents_array,int* state_array,char* inherited_state){

  int k;
  typ* parent_array=(typ*)calloc(N2,sizeof(typ));
  typ* index_parent_array=(typ*)calloc(N2,sizeof(typ));
  int children_from_cancerous=0;
  int children_from_healthy=0;
  //double odds_ratio=0;

  for (int i=0;i<N2;i++){
    double r= ((double)rand()/(RAND_MAX + 1.0)) * com_fitness[N1-1];
    //fprintf(stderr, "r: %f, rand: %f, com_fitness: %f\n", r, r, com_fitness[N1-1]);
    assert(r < com_fitness[N1-1]);
    //double r=(double)rand() / ((double)RAND_MAX + 1) * com_fitness[N1-1];
      for (k=0;k<N1;k++){
	if (com_fitness[k]>r) break;
      }

    int  parent = k;
    assert(parent<N1);
    pop2[i] = pop1[parent];
    parent_array[i]=parent;
    index_parent_array[i]=index_array[parent_array[i]];
    if (!strcmp(inherited_state,"inherited")){
      state_due_to_parents_array[i]=state_array[parent_array[i]];
    }
    else{
      state_due_to_parents_array[i]=0;
    }
    //fprintf(stderr,"The individual %d of the next population has as parent the individual %llu of the current population with position %llu\n",i,parent_array[i],index_parent_array[i]);
    //fprintf(stderr,"The parent_array is %llu\n",parent_array[i]);
    //fprintf(stderr,"The cell %d of the next population has state due to its parents %d\n",i,state_due_to_parents_array[i]);
  }

  //fprintf(stderr,"%d children are offsprings of cancerous cells and %d children are offsprings of healthy cells\n",children_from_cancerous,children_from_healthy);

  free(parent_array);
  return index_parent_array;
}


void update_state_neutral(double* array,gsl_rng* r,typ* pop1,int selected_pos,typ* index_array,int* individuals_condition,char* user_fitness_site,double constant_param,char* cell_cancer_model,FILE* f7,double* fitness, double* com_fitness, int size){
  
  fitness_func2(array,r,size,fitness,pop1,selected_pos,index_array,individuals_condition,user_fitness_site,constant_param,f7,cell_cancer_model);
  com_fitness[0]=fitness[0];
  for (int i=1; i<size; i++){
    com_fitness[i]=fitness[i]+com_fitness[i-1];
  }
}
  
  


int update_state_linear_time_dep(int* state1, int* state2, int size, double evolutionary_probability, double* fitness, double* com_fitness, double* probability, double* array, int* individuals_condition,gsl_rng* r, typ* pop1, int selected_pos, typ* index_array, char* user_fitness_site, double constant_param,int generation, int generations, FILE* f6, FILE* f5,FILE* f7,char* cell_cancer_model){

  fitness_func2(array,r,size,fitness,pop1,selected_pos, index_array, individuals_condition,user_fitness_site,constant_param,f7,cell_cancer_model);
  int count_due_to_genotype=0;

  for (int i=0; i<size; i++){
    if (individuals_condition[i]!=0 || state1[i]==1){
      state2[i]=1;
      fitness[i]=fitness[i]+10;
      count_due_to_genotype++;
    }
    else{
      state2[i]=0;
    }
    //fprintf(stderr,"The state due to genotype or due to parents is %d\n",state2[i]);
  }

 fprintf(f6, "%d\t",count_due_to_genotype);

 int count=0;
 int prev;
 int next;
 typ* temp=(typ*)calloc(size,sizeof(typ));

 for (int i=0; i<size; i++){
    temp[index_array[i]]=i;
    //fprintf(stderr,"The position of cell %d is %llu\n",i,index_array[i]);
  }

 for (int i=0; i<size; i++){
    prev=index_array[i]-1;
    if (prev<0){
      prev=size-1;
    }
    next=index_array[i]+1;
    if (next==size){
      next=0;
    }
    //fprintf(stderr,"Cell %d with site %llu has neighbor from left the site %d so cell %llu and neighbor from right the site %d so cell %llu\n",i,index_array[i],prev,temp[prev],next,temp[next]);

    if (state2[i]==0){
      if (state2[temp[prev]]==0 && state2[temp[next]]==0){
	  probability[i]=0;
	  //fprintf(stderr,"The cell %d has healthy neighbors the cell %llu and the cell %llu and no mutations with probability %f\n",i,temp[prev],temp[next],probability[i]);
      }
      else if (state2[temp[prev]]==1 && state2[temp[next]]==1){
        probability[i]=1;
        //fprintf(stderr,"The cell %d has mutant neighbors the cell %llu and the cell %llu and probability %f\n",i,temp[prev],temp[next],probability[i]);
      }
      else {
        probability[i]=generation/generations;
        //fprintf(stderr,"The cell %d has one healthy the cell %llu and one mutant the cell %llu neighbor and probability %f\n",i,temp[prev],temp[next],probability[i]);
      }
      if (probability[i]<=evolutionary_probability){
	state2[i]=0;
      }
      else {
      //state2[i]=state1[i]+1;
	state2[i]=1;
	fitness[i]=fitness[i]+10;
      //fitness[i]=((double)(fitness[prev]+fitness[next])/(double)2);
      }
    }



    if (state2[i]==1){
      count++;
      }

    //fprintf(stderr,"The total fitenss for selection for individual %d is %f\n",i,fitness[i]);
    //fprintf(stderr,"The new state array for individual %d is %d\n",i,state2[i]);

    com_fitness[0]=fitness[0];
    if (i>0){
      com_fitness[i]=fitness[i]+com_fitness[i-1];
      //fprintf(stderr, "Com_fitness of %d is %f\n",i,com_fitness[i]);
    }
 }

 //fprintf(stderr,"The total number of 1 states is %d\n",count);

 fprintf(f5, "%d\t",count);

 free(temp);
 return count;
}


int update_state_no_neighbors(int* state1, int* state2, int size, double* fitness, double* com_fitness, double* array, int* individuals_condition, gsl_rng* r, typ* pop1, int selected_pos, typ* index_array, char* user_fitness_site, double constant_param, FILE* f6, FILE* f5, FILE* f7,char* cell_cancer_model){

  fitness_func2(array, r, size, fitness, pop1, selected_pos, index_array, individuals_condition, user_fitness_site, constant_param,f7,cell_cancer_model);
  int count=0;

  for (int i=0; i<size; i++){
    if (individuals_condition[i]!=0 || state1[i]==1){
      state2[i]=1;
      fitness[i]=fitness[i]+10;
      count++;
    }
    else{
      state2[i]=0;
    }
    //fprintf(stderr, "The cell %d has state %d due to genotype or due to its parents\n", i, state2[i]);
  }

  //fprintf(stderr,"The total number of 1 states is %d\n",count);

  com_fitness[0]=fitness[0];
  for (int k=1; k<size; k++){
    com_fitness[k]=fitness[k]+com_fitness[k-1];
  }
  fprintf(f5, "%d\t",count);
  fprintf(f6, "%d\t",count);

  return count;
}


int update_state_linear_time_ind(int* state1, int*state2, int size, double evolutionary_probability, double* fitness, double* com_fitness, double* probability,double* array,int* individuals_condition, gsl_rng* r, typ* pop1, int selected_pos, typ* index_array, char* user_fitness_site, double constant_param, FILE* f6, FILE* f5, FILE* f7,char* cell_cancer_model){

  fitness_func2(array,r,size,fitness,pop1,selected_pos, index_array, individuals_condition,user_fitness_site,constant_param,f7,cell_cancer_model);
  int count_due_to_genotype=0;
  int count=0;
  int prev;
  int next;
  typ* temp=(typ*)calloc(size,sizeof(typ));

  for (int i=0; i<size; i++){
    temp[index_array[i]]=i;
    //fprintf(stderr,"The position of cell %d is %llu\n",i,index_array[i]);                                                                                                                                 
  }

  for (int i=0; i<size; i++){
    prev=index_array[i]-1;
    if (prev<0){
      prev=size-1;
    }
    next=index_array[i]+1;
    if (next==size){
      next=0;
    }
  }




  for (int i=0; i<size; i++){
      if (individuals_condition[i]!=0 || state1[i]==1){
      state2[i]=1;
      fitness[i]=fitness[i]+10;
      count_due_to_genotype++;
      count++;
      }
      else{
	 if (state2[temp[prev]]==0 && state2[temp[next]]==0){
         probability[i]=0;
         //fprintf(stderr,"The cell %d has healthy neighbors the cell %llu and the cell %llu and no mutations with probability %f\n",i,temp[prev],temp[next],probability[i]);                               
       }
       else if (state2[temp[prev]]==1 && state2[temp[next]]==1){
         probability[i]=1;
         //fprintf(stderr,"The cell %d has mutant neighbors the cell %llu and the cell %llu and probability %f\n",i,temp[prev],temp[next],probability[i]);                                                  
       }
       else {
         probability[i]=0.5;
         //fprintf(stderr,"The cell %d has one healthy the cell %llu and one mutant the cell %llu neighbor and probability %f\n",i,temp[prev],temp[next],probability[i]);                                   
       }

       if (probability[i]<=evolutionary_probability){
         state2[i]=0;
       }
       else {
      //state2[i]=state1[i]+1;                                                                                                                                                                              
         state2[i]=1;
	 count++;
         fitness[i]=fitness[i]+10;

       }
      //fprintf(stderr,"The state due to genotype or due to parents is %d\n", state2[i]);
      }
  }

  fprintf(f6, "%d\t",count_due_to_genotype);


    //fprintf(stderr,"The total fitenss for selection for individual %d is %f\n",i,fitness[i]);
    //fprintf(stderr,"The new state array for individual %d is %d\n",i,state2[i]);
  for (int i=0; i<size; i++){
    com_fitness[0]=fitness[0];
    if (i>0){
      com_fitness[i]=fitness[i]+com_fitness[i-1];
      //fprintf(stderr, "Com_fitness of %d is %f\n",i,com_fitness[i]);
    }
  }

  //fprintf(stderr,"The total number of 1 states is %d\n", count);

  fprintf(f5, "%d\t",count);

  free(temp);
  return count;
}


int update_state_no_genotype_effect(int* state1,int* state2, int size, double evolutionary_probability, double* fitness,double* com_fitness, double* probability, int generation, long double growth_rate, long double decay_rate, double radius,int* individuals_condition, double* array, gsl_rng* r, typ* pop1, int selected_pos, typ* index_array,char* user_fitness_site,double constant_param,int initial_cancer_cells, FILE* f5,FILE* f7,char* cell_cancer_model){

  fitness_func2(array,r,size,fitness,pop1,selected_pos, index_array, individuals_condition,user_fitness_site,constant_param,f7,cell_cancer_model);
  int count=0;

  if (generation==1){
      state1[size/2]=1;
      //state1[0]=1;
      //state1[size-1]=0;
    }

  for (int i=0; i<size; i++){
    if (state1[i]==1){
      state2[i]=1;
      fitness[i]=fitness[i]+10;
    }
    else{
      state2[i]=0;
    }
    //fprintf(stderr,"The cell %d has state due to its parents %d\n",i,state2[i]);
  }

  int prev;
  int next;
  typ* temp=(typ*)calloc(size,sizeof(typ));
  double radius_pow=pow(radius,(double)3);
  double prob_max_tum_increment=(double)(2*4*3.14*radius_pow)/(double)3;
  double Gompertz_param=(double)(Gompertz_discrete_model(initial_cancer_cells,generation,growth_rate,decay_rate,radius))/prob_max_tum_increment;
  //fprintf(stderr,"The Compertz_param for generations %d is %f\n",generation,Gompertz_param);
  assert(Gompertz_param<1);


  for (int i=0; i<size; i++){
    temp[index_array[i]]=i;
    //fprintf(stderr,"The position of cell %d is %llu\n",i,index_array[i]);
  }
  for (int i=0; i<size; i++){
    prev=index_array[i]-1;
    if (prev<0){
      prev=size-1;
    }
    next=index_array[i]+1;
    if (next==size){
      next=0;
    }
    //fprintf(stderr,"Cell %d with site %llu has neighbor from left the site %d so cell %llu and neighbor from right the site %d so cell %llu\n",i,index_array[i],prev,temp[prev],next,temp[next]);

    if (state2[i]==0){
      if (state2[temp[prev]]==0 && state2[temp[next]]==0){
	probability[i]=0;
	//fprintf(stderr,"The cell %d has healthy neighbors the cell %llu and the cell %llu and no mutations with probability %f\n",i,temp[prev],temp[next],probability[i]);
      }

      else if (state2[temp[prev]]==1 && state2[temp[next]]==1){
	probability[i]=1;
	//fprintf(stderr,"The cell %d has mutant neighbors the cell %llu and the cell %llu and probability %f\n",i,temp[prev],temp[next],probability[i]);
      }
      else {
	probability[i]=Gompertz_param;
	//fprintf(stderr,"The cell %d has one healthy the cell %llu and one mutant the cell %llu neighbor and probability %f\n",i,temp[prev],temp[next],probability[i]);
      }

      if (probability[i]<=evolutionary_probability){
	state2[i]=0;
      }
      else {
      //state2[i]=state1[i]+1;
	state2[i]=1;
      //fitness[i]=((double)(fitness[prev]+fitness[next])/(double)2);
	fitness[i]=fitness[i]+10;
      }
    }

    if (state2[i]==1){
      count++;
    }
    //fprintf(stderr,"The total fitenss for selection for individual %d is %f\n",i,fitness[i]);
    //fprintf(stderr,"The new state array for individual %d is %d\n",i,state2[i]);
    com_fitness[0]=fitness[0];
    if (i>0){
      com_fitness[i]=fitness[i]+com_fitness[i-1];
      //fprintf(stderr, "Com_fitness of %d is %f\n",i,com_fitness[i]);
    }
  }


  //fprintf(stderr,"The total number of 1_states is %d\n",count);
  fprintf(f5, "%d\t",count);

  free(temp);
  return count;
}


int update_state_non_linear_gom_time_dep(int* state1,int* state2, int size, double evolutionary_probability, double* fitness,double* com_fitness, double* probability, int generation, long double growth_rate, long double decay_rate, double radius,int* individuals_condition, double* array, gsl_rng* r, typ* pop1, int selected_pos, typ* index_array,char* user_fitness_site,double constant_param,int initial_cancer_cells,int* state_array,FILE* f6, FILE* f5, FILE* f7,char* cell_cancer_model){

  //double* fitness_probs=(double*)calloc(size,sizeof(double));

  fitness_func2(array,r,size,fitness,pop1,selected_pos, index_array, individuals_condition,user_fitness_site,constant_param,f7,cell_cancer_model);
  int count_due_to_genotype=0;
  int count=0;
  int prev;
  int next;
  typ* temp=(typ*)calloc(size,sizeof(typ));
  double radius_pow=pow(radius,(double)3);
  double prob_max_tum_increment=(double)(2*4*3.14*radius_pow)/(double)3;
  double Gompertz_param=(double)(Gompertz_discrete_model(initial_cancer_cells,generation,growth_rate,decay_rate,radius))/prob_max_tum_increment;
  //fprintf(stderr,"The Compertz_param for generations %d is %f\n",generation,Gompertz_param);                                                                                                              
  assert(Gompertz_param<1);

  for (int i=0; i<size; i++){
    temp[index_array[i]]=i;
    //fprintf(stderr,"The position of cell %d is %llu\n",i,index_array[i]);                                                                                                                                 
  }
  for (int i=0; i<size; i++){
    prev=index_array[i]-1;
    if (prev<0){
      prev=size-1;
    }
    next=index_array[i]+1;
    if (next==size){
      next=0;
    }
    //fprintf(stderr,"Cell %d with site %llu has neighbor from left the site %d so cell %llu and neighbor from right the site %d so cell %llu\n",i,index_array[i],prev,temp[prev],next,temp[next]);
  }




  for (int i=0; i<size; i++){
    if (individuals_condition[i]!=0 || state1[i]==1){
      state2[i]=1;
      fitness[i]=fitness[i]+10;
      count_due_to_genotype++;
      count++;
    }
    else{
      if (state2[temp[prev]]==0 && state2[temp[next]]==0){
        probability[i]=0;
        //fprintf(stderr,"The cell %d has healthy neighbors the cell %llu and the cell %llu and no mutations with probability %f\n",i,temp[prev],temp[next],probability[i]);                                
      }

      else if (state2[temp[prev]]==1 && state2[temp[next]]==1){
        probability[i]=1;
        //fprintf(stderr,"The cell %d has mutant neighbors the cell %llu and the cell %llu and probability %f\n",i,temp[prev],temp[next],probability[i]);                                                   
      }
      else {
        probability[i]=Gompertz_param;
        //fprintf(stderr,"The cell %d has one healthy the cell %llu and one mutant the cell %llu neighbor and probability %f\n",i,temp[prev],temp[next],probability[i]);                                    
      }

      if (probability[i]<=evolutionary_probability){
        state2[i]=state1[i];
      }
      else {
      //state2[i]=state1[i]+1;                                                                                                                                                                              
        state2[i]=1;
      //fitness[i]=((double)(fitness[prev]+fitness[next])/(double)2);                                                                                                                                       
        fitness[i]=fitness[i]+10;
	count++;
      }

    }
    //fprintf(stderr,"The state of cell %d due to its genotype or due to its parents is %d\n",i,state2[i]);
  }


  fprintf(f6, "%d\t",count_due_to_genotype);

  for (int i=0; i<size; i++){

    com_fitness[0]=fitness[0];
    if (i>0){
      com_fitness[i]=fitness[i]+com_fitness[i-1];
      if(com_fitness[i] < 0 || fitness[i] < 0){
	//fprintf(stderr, "com: %f, fit: %f\n", com_fitness[i], fitness[i]);
	assert(com_fitness[i] > 0);
	assert(fitness[i] > 0);
      }
    }//fprintf(stderr, "Com_fitness of %d is %f\n",i,com_fitness[i]);
  }

  //fitness_probs=softmax(fitness,size);


  //fprintf(stderr,"The total number of 1_states is %d\n",count);
  fprintf(f5, "%d\t",count);

  free(temp);
  return count;
}

void find_1_states_positions(int* state,int size, int num_of_1s, typ* pop, int total_num,FILE* f3, FILE* f4){
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

  //  assert(j == num_of_1s);

  for (int i=0; i<num_of_1s; i++){
    //fprintf(stderr,"The 1 is in position %d\n",states_1_array[i]);
    mutant_genots_array[i]=pop[states_1_array[i]];
    //fprintf(stderr,"The mutant genotype is %llu\n",mutant_genots_array[i]);
  }

  freq_func_mutant(mutant_genots_array,num_of_1s,f3);
  setbits_positions_mutant(num_of_1s, mutant_genots_array, total_num,f4);

  free(states_1_array);
  free(mutant_genots_array);

}

double* weight_position_array_func(gsl_rng*p_random,int size, FILE* f10){

  double* weight_position_array1=(double*)calloc(size,sizeof(double));

  for (int i=0;i<size;i++){
    weight_position_array1[i]=gsl_ran_gaussian(p_random,1);
    //weight_position_array1[i]=1;
    fprintf(f10,"%d\t%f\n",i,weight_position_array1[i]);
    // weight_position_array[i]=(float)rand()/(float)RAND_MAX*10-11;                                                                                                                                       \

  }
  return weight_position_array1;
}


void population_change_via_branching(gsl_rng* r, int initialN, int generations, double lambda, double ChangingFactor, int* population_array)
{

  //int* population_array=(int*)calloc(generations,sizeof(int));
    population_array[0]=initialN;
    for (int i = 1; i < generations; ++i)
    {
      lambda=lambda*ChangingFactor;
      if (lambda<0){
	lambda=0;
      }
      int x=0;
      if (population_array[i-1]>0){
	for (int j=0; j<population_array[i-1]; j++){
	  x+=gsl_ran_poisson(r, 1+lambda);
	}
	population_array[i]=x;
      }


    }

    //return population_array;
}

void print_genome_without(int N, typ* pop, int generation){

  for(int k = 0; k < N; ++k){
      char* s=bitfunc(pop[k]);
      //fprintf(stderr,"%d  -- The genome (WITHOUT) of individual %d of the current population is %s\n", generation, k,s);
      free(s);
    }
}


void mutation_process(int number_of_mutations, long long int positions, typ* pop, int generation){

  typ* position_of_mutations_array=(typ*)calloc(number_of_mutations,sizeof(typ));
    for (int j=0;j<number_of_mutations;j++){
      position_of_mutations_array[j]=mutate(positions, pop, generation); // οι θέσεις που έγιναν μεταλλάξεις μπαινουν σε ενα array
    }

    free(position_of_mutations_array);
}


int print_genome_with(int N, typ* pop){
  int x=0;
  for(int k = 0; k < N; ++k){
      char* s=bitfunc(pop[k]);
      int num_of_setbits=countSetBits(pop[k]);
      x+=num_of_setbits;
      //fprintf(stderr,"The genome (WITH) of individual %d of the current population is %s\n",k,s);
      //fprintf(stderr,"The individual %d has %d number of setbits/mutations\n",k,num_of_setbits);
      //individual_setbits_positions(num_of_setbits,pop[k]);
      free(s);
  }

  if (x==0){
    x=1;
  }


  //fprintf(stderr,"The total number of setbits for this generation are %d\n",x);

  return x;
}

void population_realloc(int N, typ* pop1, typ* pop2){

  for (int k=0; k<N; k++){
      pop1[k]=pop2[k];

    }
}

void common_positions(typ n, FILE* f7){
  fprintf(f7, "The common positions are:");
  unsigned int common_pos=0;
  while(n){
    if((n&1)!=0){
      fprintf(f7, "%d\t",common_pos);
    }
    common_pos++;
    n>>=1;
  }
  fprintf(f7,"\n");
}

void common_fixed_positions(typ n, FILE* f7){
  fprintf(f7, "The common fixed positions are:");
  unsigned int common_pos=0;
  while(n){
    if((n&1)!=0){
      fprintf(f7, "%d\t",common_pos);
    }
    common_pos++;
    n>>=1;
  }
  fprintf(f7,"\n");
}


void fixed_positions(int N, typ* pop,FILE* f7){
  int setBits=0;
  int commonBits=0;
  typ final_pop=0;

  for (int i=0;i<N;i++){
    final_pop=pop[i];
    setBits = countSetBits(pop[i]);
    commonBits = countSetBits(pop[i] & global_fitfactor); // οι κοινες θεσεις ( άσσοι) μεταξύ των άσσεων του fitfactor και του μεταλλαγμένου γονιδιώματος).
    //fprintf(stderr,"Individual:%d\t Number of positions that have been fixed (from all the positions):%d\t Number of positions that have contibuted (from those that matter):%d\n", i, setBits, commonBits\
\
);

  }
  fprintf(f7,"Number of positions that have been fixed(from all the positions):%d\t Number of positions that have contirbuted (from those that matter):%d\n",setBits,commonBits);
  common_positions(final_pop & global_fitfactor, f7);
  common_fixed_positions(final_pop,f7);
}


void average_fitness(int N, double* fitness, int generation, FILE* f8){
  double* sum_fitness=(double*)calloc(N,sizeof(double));
  double* dif_avg=(double*)calloc(N,sizeof(double));
  double* sum_dif_avg=(double*)calloc(N,sizeof(double));
  double average_fitness=0;
  double variance=0;
  sum_fitness[0]=fitness[0];
  int i;
  for (i=1; i<N; i++)
    sum_fitness[i]=fitness[i]+sum_fitness[i-1];

  average_fitness=(double)sum_fitness[N-1]/(double)N;

  for (int i=0; i<N; i++){
    dif_avg[i]=pow((fitness[i]-average_fitness), (double)2);
  }
  sum_dif_avg[0]=dif_avg[0];

  for (int i=1; i<N; i++){
    sum_dif_avg[i]=dif_avg[i]+sum_dif_avg[i-1];
  }

  variance=(double)sum_dif_avg[N-1]/(double)(N-1);
  fprintf(f8,"Generation %d\t%f\t%f\n",generation,average_fitness,variance);

  free(sum_fitness);
  free(dif_avg);
  free(sum_dif_avg);
}


int main(int argc, char* argv[]){
  time_t start, end;
  time(&start);
  global_fitfactor = 0;
  int N;
  int generations;
  double constant_param;
  int fixation=0;
  char user_fitness_site[100];
  char user_neighbor_effect[100];
  char cell_cancer_model[100];
  char inherited_state[100];
  char population_model[100];
  user_fitness_site[0] = '\0';
  user_neighbor_effect[0] = '\0';
  cell_cancer_model[0]  = '\0';
  inherited_state[0] = '\0';
  population_model[0] = '\0';

  int cmd_i;
  int pop_size;
  int pop_next_size;
  unsigned seed;
  //int* population_array;

  for (cmd_i=1; cmd_i<argc; cmd_i++){
    if( (!strcmp(argv[cmd_i], "-N" ) ) ){
      N=atoi(argv[++cmd_i]);
      continue;                                                                                                                                                                                          \

    }
    if ((!strcmp(argv[cmd_i], "-gens"))){
      generations=atoi(argv[++cmd_i]);
      continue;
      //exit(0);                                                                                                                                                                                           \

    }
    if ((!strcmp(argv[cmd_i], "-beta"))){
      memcpy(user_fitness_site, "beta", strlen("beta")+1);
      continue;
      //global_fitness_site_array=fitness_site_func5(N,r);
    }
    if ((!strcmp(argv[cmd_i], "-linear"))){
      memcpy(user_fitness_site, "linear", strlen("linear")+1);
      continue;
    }
    if ((!strcmp(argv[cmd_i], "-constant"))){
      memcpy(user_fitness_site, "constant", strlen("constant")+1);
      constant_param=atof(argv[++cmd_i]);
      continue;
    }
    if ((!strcmp(argv[cmd_i], "-alter"))){
      memcpy(user_fitness_site, "alter", strlen("alter")+1);
      continue;
    }

    if ((!strcmp(argv[cmd_i], "-linear_time_ind"))){
      memcpy(user_neighbor_effect, "linear_time_ind", strlen("linear_time_ind")+1);
      continue;
    }
    if ((!strcmp(argv[cmd_i], "-linear_time_dep"))){
      memcpy(user_neighbor_effect, "linear_time_dep", strlen("linear_time_dep")+1);
      continue;
    }
    if ((!strcmp(argv[cmd_i], "-gom_time_dep"))){
      memcpy(user_neighbor_effect, "gom_time_dep", strlen("gom_time_dep")+1);
      continue;
    }
    if ((!strcmp(argv[cmd_i], "-no_genotype_effect"))){
      memcpy(user_neighbor_effect, "no_genotype_effect", strlen("no_genotype_effect")+1);
      continue;
    }
    if ((!strcmp(argv[cmd_i],"-no_neighbors_effect"))){
      memcpy(user_neighbor_effect, "no_neighbors_effect", strlen("no_neighbors_effect")+1);
      continue;
    }

    if ((!strcmp(argv[cmd_i],"-neutral_nei"))){
      memcpy(user_neighbor_effect, "neutral", strlen("neutral")+1);
      continue;
    }

    if ((!strcmp(argv[cmd_i],"-deterministic"))){
      memcpy(cell_cancer_model, "deterministic", strlen("deterministic")+1);
      continue;
    }

    if ((!strcmp(argv[cmd_i],"-neutral_ca"))){
      memcpy(cell_cancer_model, "neutral", strlen("neutral")+1);
      continue;
    }

    if ((!strcmp(argv[cmd_i],"-stochastic"))){
      memcpy(cell_cancer_model, "stochastic", strlen("stochastic")+1);
      continue;
    }


    if ((!strcmp(argv[cmd_i],"-inherited"))){
      memcpy(inherited_state, "inherited", strlen("inherited")+1);
      continue;
    }

    if ((!strcmp(argv[cmd_i], "-neutral_st"))){
      memcpy(inherited_state, "neutral", strlen("neutral")+1);
      continue;
    }

    if ((!strcmp(argv[cmd_i], "-non_herited"))){
      memcpy(inherited_state, "non_herited", strlen("non_herited")+1);
      continue;
    }


    if ((!strcmp(argv[cmd_i], "-branching"))){
      memcpy(population_model, "branching", strlen("branching")+1);
      continue;
    }

    if ((!strcmp(argv[cmd_i], "-non_branching"))){
      memcpy(population_model, "non_branching", strlen("non_branching")+1);
      continue;
    }

    if ((!strcmp(argv[cmd_i], "-seed"))){
      seed=atoi(argv[++cmd_i]);
      continue;
    }

    fprintf(stderr,"Argument %s does not exist\n", argv[cmd_i]);
    exit(-1);
  }

  FILE* f1=fopen("Frequencies_CA.txt","w");
  FILE* f2=fopen("Frequencies_mutations_CA.txt","w");
  FILE* f3=fopen("Frequencies_CA_mutant.txt","w");
  FILE* f4=fopen ("Frequencies_mutations_CA_mutant.txt","w");
  FILE* f5=fopen("Number of mutants.txt","w");
  FILE* f6=fopen("Number_of_mutants_genotype.txt","w");
  FILE* f7=fopen("Positions_that_matter.txt","w");
  FILE* f8=fopen("Average_fitness.txt","w");
  FILE* f9=fopen("Population_size.txt","w");
  FILE* f10=fopen("Position_weights.txt", "w");
  
  int initialN=N; // για Ν σταθερό έτσι όπως έχουμε φτιάξει τον κώδικα ακολουθεί μοντέλο Wright-Fisher.
  int states_1_num;
  double radius=37.5;
  int initial_cancer_cells;
  //int total_children_from_cancerous;
  //int total_children_from_healthy;
  //unsigned seed=3;
  //unsigned seed = 25; //(unsigned) time(&t);
  //unsigned seed=time(&t);
  unsigned seed_p_random=3;

  srand(seed);


  typ* pop_prev=(typ*)calloc(N,sizeof(typ));
  typ* pop_cur=(typ*)calloc(N,sizeof(typ));
  int* state_array=(int*)calloc(N,sizeof(int));
  int* state_due_to_parents_array=(int*)calloc(N,sizeof(int));
  double* fitness=(double*)calloc(N,sizeof(double));
  double* com_fitness=(double*)calloc(N,sizeof(double));
  double* probability_array=(double*)calloc(N,sizeof(double));
  typ* offspring_array=(typ*)calloc(N,sizeof(typ));
  typ* initial_index_array=(typ*)calloc(N,sizeof(typ));
  int* individuals_condition=(int*)calloc(N,sizeof(int));
  int* population_array=(int*)calloc(generations+1,sizeof(int));

  double mu=((double)pow((double)10,(double)-4)); // for cancer cells, while for normal cells is 10^-8 or 10^-6;
  //fprintf(stderr,"The mutation rate is %f\n",mu);
  long double growth_rate=log(2);
  //int selected_pos=rand()%64+1;
  int selected_pos=15; // I keep it constant in order to find if there is a certain motif of finding specific mutations.
  //fprintf(stderr,"The positions that matter are %d\n",selected_pos);


  fprintf(f7,"The positions that matter are %d\n", selected_pos);

  fprintf(f5,"The total number of initial healthy cells are %d\n",initialN);


  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng * p_random;
  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  p_random=gsl_rng_alloc(T);

  gsl_rng_set(r, seed);
  gsl_rng_set(p_random,seed_p_random);

  double evolutionary_probability=gsl_ran_flat(r,0,1);
  //double evolutionary_probability=0.89;
  //double evolutionary_probability=0.01;
  //fprintf(stderr,"The evolutionary probability is %f\n",evolutionary_probability);

  fprintf(f7,"The evolutionary probability is %f\n", evolutionary_probability);

  double* weight_position_array=weight_position_array_func(p_random,intsize, f10);

  if (!strcmp(population_model,"branching")){
    //population_array=population_change_via_branching(r,initialN,generations,2,0.8);
    population_change_via_branching(r,initialN,generations,2,0.8, population_array);
    for (int i=0; i<generations; i++){
      fprintf (stderr, "%d\n", population_array[i]);
    }
  }

  else {
    if (!strcmp(user_fitness_site,"beta")){
      global_fitness_site_array=fitness_site_func5(N,r);
    }

    if (!strcmp(user_fitness_site,"linear")){
      global_fitness_site_array=fitness_site_func_linear(N);
    }

    if (!strcmp(user_fitness_site,"constant")){
      global_fitness_site_array=fitness_site_func_constant(N,constant_param);
    }
    if(!strcmp(user_fitness_site,"alter")){
      global_fitness_site_array=fitness_site_func_alter(N);
    }
  }

  for (int i=0; i<initialN; i++){
    //offspring_array[i]=i;
    initial_index_array[i]=i;
  }

  for (int i=0; i<generations; i++){
    //if (i>300){
    //mu=0;
    //}

    if (!strcmp(population_model,"branching")){
      pop_size=population_array[i];
      pop_next_size=population_array[i+1];

      if (pop_next_size>pop_size){
	pop_cur=(typ*)realloc(pop_cur,sizeof(typ)*pop_next_size);
	pop_prev=(typ*)realloc(pop_prev,sizeof(typ)*pop_next_size);
	fitness=(double*)realloc(fitness,sizeof(typ)*pop_next_size);
	com_fitness=(double*)realloc(com_fitness,sizeof(typ)*pop_next_size);
	state_array=(int*)realloc(state_array,sizeof(int)*pop_next_size);
	state_due_to_parents_array=(int*)realloc(state_due_to_parents_array,sizeof(int)*pop_next_size);
	probability_array=(double*)realloc(probability_array,sizeof(double)*pop_next_size);
	offspring_array=(typ*)realloc(offspring_array,sizeof(typ)*pop_next_size);
	initial_index_array=(typ*)realloc(initial_index_array,sizeof(typ)*pop_next_size);
	individuals_condition=(int*)realloc(individuals_condition,sizeof(int)*pop_next_size);
	assert(pop_cur!=NULL);
	assert(pop_prev!=NULL);
	assert(fitness!=NULL);
	assert(com_fitness!=NULL);
      }

      if (!strcmp(user_fitness_site, "linear")){
	global_fitness_site_array=fitness_site_func_linear(pop_size);
      }

      if (!strcmp(user_fitness_site,"beta")){
	global_fitness_site_array=fitness_site_func5(pop_size,r);
      }

      if (!strcmp(user_fitness_site,"constant")){
	global_fitness_site_array=fitness_site_func_constant(pop_size,constant_param);
      }

      if(!strcmp(user_fitness_site,"alter")){
	global_fitness_site_array=fitness_site_func_alter(pop_size);
      }

    }

    else {
      pop_size=pop_next_size=initialN=N;
    }

    fprintf(f9,"Generation %d\t\%d\n",i,pop_size);

    for (int i=0; i<pop_next_size; i++){
      offspring_array[i]=i;
    }

    long long int total_pop_pos=pop_size*intsize;


    int numberOfMutations=binomial(r,mu,pop_size);
    //fprintf(stderr,"The %d generation which has undergone %d  mutations\n",i,numberOfMutations);

    print_genome_without(pop_size,pop_prev,i);


    mutation_process(numberOfMutations,total_pop_pos,pop_prev,i);


    int total_setbits_num=print_genome_with(pop_size,pop_prev);


    //fprintf(stderr,"The genome frequencies are:\n");

    fprintf(f1,"Generation %d\t",i);

    fixation=freq_func(pop_prev,pop_size,f1,i);// τρέχει σωστά και βγάζει συχνότητα γοιδιωμάτων μετα την μετάλλαξη.


    fprintf(f1,"\n");

    if (i==0){
      fitness_func2(weight_position_array,r,pop_size,fitness,pop_prev,selected_pos, initial_index_array, individuals_condition,user_fitness_site,constant_param,f7,cell_cancer_model);
      com_fitness[0]=fitness[0];

      for (int k=1; k<pop_size; k++){
	//fprintf(stderr, "k: %d, popsize: %d, nextpop: %d\n", k, pop_size, pop_next_size);
	com_fitness[k]=fitness[k]+com_fitness[k-1];
	//fprintf(stderr,"The com fitness of individual %d is %f\n",k,com_fitness[k]);
      }
      int count=0;
      //state_array[pop_size/2]=1;

      for (int j=0; j<pop_size; j++){
	if (individuals_condition[j]!=0){
	  state_array[j]=1;
	  count++;
	}
	//fprintf(stderr,"The initial state array in generation %d for individual %d is %d\n",i,j,state_array[j]);
      }
      initial_cancer_cells=count;

      fprintf(f5,"Generation %d\t",i);

      //fprintf(stderr,"The initial number of 1_states are %d\n",count);
      fprintf(f5, "%d\t",count);
      fprintf(f5,"\n");

      if (strcmp(user_neighbor_effect,"no_genotype_effect") || strcmp(user_neighbor_effect,"neutral_nei")){
	fprintf(f6, "Generation %d\t",i);
	fprintf(f6, "%d\t",count);
	fprintf(f6,"\n");
      }

      fprintf(f3,"Generation %d\t",i);
      fprintf(f4,"Generation %d\t",i);
      find_1_states_positions(state_array, pop_size, count,pop_prev, total_setbits_num,f3,f4);
      fprintf(f3,"\n");
      fprintf(f4,"\n");
    }

    if (i>0){

    fprintf(f5,"Generation %d\t",i);

    if (!strcmp(cell_cancer_model,"neutral")){
      update_state_neutral(weight_position_array,r,pop_prev,selected_pos,initial_index_array,individuals_condition,user_fitness_site,constant_param,cell_cancer_model,f7,fitness,com_fitness,pop_size);
    }

    if (strcmp(user_neighbor_effect,"no_genotype_effect") || strcmp(user_neighbor_effect,"neutral_nei")){

      fprintf(f6,"Generation %d\t",i);
    }

    if (!strcmp(user_neighbor_effect,"gom_time_dep")){


	states_1_num=update_state_non_linear_gom_time_dep(state_due_to_parents_array,state_array, pop_size, evolutionary_probability, fitness,com_fitness, probability_array, i, growth_rate,  0.2, radius,individuals_condition,weight_position_array,r,pop_prev,selected_pos,initial_index_array,user_fitness_site,constant_param,initial_cancer_cells,state_array,f6,f5,f7,cell_cancer_model);
      }

    if (!strcmp(user_neighbor_effect,"linear_time_ind")){

	states_1_num=update_state_linear_time_ind(state_due_to_parents_array, state_array, pop_size, evolutionary_probability, fitness, com_fitness, probability_array,weight_position_array,individuals_condition,r, pop_prev, selected_pos, initial_index_array, user_fitness_site, constant_param,f6,f5,f7,cell_cancer_model);
      }

    if (!strcmp(user_neighbor_effect,"linear_time_dep")){

	states_1_num=update_state_linear_time_dep(state_due_to_parents_array,state_array, pop_size, evolutionary_probability, fitness, com_fitness, probability_array, weight_position_array, individuals_condition,r,pop_prev, selected_pos, initial_index_array, user_fitness_site, constant_param, i,generations,f6,f5,f7,cell_cancer_model);
      }

    if (!strcmp(user_neighbor_effect,"no_genotype_effect")){
	states_1_num= update_state_no_genotype_effect(state_due_to_parents_array,state_array,pop_size, evolutionary_probability, fitness, com_fitness, probability_array, i,growth_rate,0.2,radius, individuals_condition, weight_position_array, r, pop_prev, selected_pos, initial_index_array, user_fitness_site, constant_param, initial_cancer_cells,f5,f7,cell_cancer_model);
      }

    if (!strcmp(user_neighbor_effect,"no_neighbors_effect")){
        states_1_num=update_state_no_neighbors(state_due_to_parents_array, state_array, pop_size, fitness, com_fitness, weight_position_array, individuals_condition, r, pop_prev, selected_pos, initial_index_array, user_fitness_site, constant_param, f6, f5,f7,cell_cancer_model);

      }

      fprintf(f5,"\n");

      if (strcmp(user_neighbor_effect,"no_genotype_effect") || strcmp(user_neighbor_effect,"neutral")){
	fprintf(f6,"\n");
      }

      if (strcmp(user_neighbor_effect,"neutral") || strcmp(cell_cancer_model,"neutral") || strcmp(inherited_state,"neutral")){

      //fprintf(stderr,"The mutant genome frequencies are:\n");

      fprintf(f3,"Generation %d\t",i);

      fprintf(f4,"Generation %d\t",i);

      find_1_states_positions(state_array, pop_size, states_1_num,pop_prev, total_setbits_num,f3,f4);

      fprintf(f3,"\n");
      fprintf(f4,"\n");
      }

    }

    typ* index_parents_array= selection_func(pop_size, pop_next_size, pop_prev, pop_cur, fitness, com_fitness,initial_index_array,state_due_to_parents_array,state_array,inherited_state);

    typ* cur_index_array=sorting_1(index_parents_array,offspring_array,pop_next_size);

      for (int j=0; j<pop_next_size; j++){
      initial_index_array[j]=cur_index_array[j];
      offspring_array[j]=j;
      //fprintf(stderr,"The new initial index is %llu\n",initial_index_array[j]);
      }

      fprintf(f2,"Generation %d\t",i);

      setbits_positions(pop_size,pop_prev,total_setbits_num,f2);

      fprintf(f2,"\n");

      average_fitness(pop_size,fitness,i,f8);

      population_realloc(pop_next_size, pop_prev, pop_cur);

      free(cur_index_array);
      free(index_parents_array);
      if (!strcmp(population_model,"branching")){
	free(global_fitness_site_array);
      }


      if (fixation==1){
	break;
      }

  }

  fixed_positions(pop_size,pop_cur,f7);

  pop_size=pop_next_size;

  //for (int i=0;i<generations;i++){
  //fprintf(f9,"%d\n",population_array[i]);
  //}

  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(f4);
  fclose(f5);
  fclose(f6);
  fclose(f7);
  fclose(f8);
  fclose(f9);
  fclose(f10);

  time(&end);


  double time_taken=end-start;
  fprintf(stderr,"The time taken for the programm to run is %f sec\n",time_taken);

  free(pop_prev);
  free(pop_cur);
  free(fitness);
  free(com_fitness);
  gsl_rng_free(r);
  gsl_rng_free(p_random);
  free(weight_position_array);
  free(state_array);
  free(state_due_to_parents_array);
  free(probability_array);
  free(individuals_condition);
  free(offspring_array);
  free(initial_index_array);
  //free(global_fitness_site_array);
  free(population_array);
  if (strcmp(population_model,"branching")){
    free (global_fitness_site_array);
  }

  return 0;
}
