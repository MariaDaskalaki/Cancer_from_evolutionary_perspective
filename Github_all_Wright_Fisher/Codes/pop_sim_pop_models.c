#include "pop_sim_pop_models"

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
  //printf("%llu\n",newSize);                                                                                                                                                                              \
                                                                                                                                                                                                            
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
