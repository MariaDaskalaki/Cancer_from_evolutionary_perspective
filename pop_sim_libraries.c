#include "pop_sim_libraries.h"

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




