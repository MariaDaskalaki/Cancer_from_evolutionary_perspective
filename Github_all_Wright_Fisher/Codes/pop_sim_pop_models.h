#include "pop_sim_libraries.h"


typ exponential_linear(int initial, int t, long double growth_rate, double linear_const);

typ generalized_logistic(int initial,int t, long double capacity, long double growth_rate, double n_const);

typ Gompertz_model(int initial,int t, long double growth_rate, long double decay_rate);

typ dcc_model(int initial, int t, long double growth_rate, int b_const);

typ von_bertalanffy(int initial, int t , long double growth_rate, long double decay_rate);



