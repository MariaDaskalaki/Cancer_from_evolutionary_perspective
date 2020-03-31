# Cancer_from_evolutionary_perspective

Forward time simulations of cancer cells in discrete time of generations, based on two populations(previous cell population and current cell population).

The initial cell population consists of healthy cells.

The code is written in C and takes two command line parameters: initial number of cells and number of generations.

A "Frequency.txt" and a "Frequency_mutations.txt" file are derived from the C-code. 

Plots of txt files are derived from the corresponding python codes.

The number of mutations for each generation is a random number from binomial distribution, inside the "pool" of total number of cells for the generation * genome_positions(intsize=64) with mutation rate mu=10^-4. The mutation rate of normal cells is between 10^-8 to 10^-6, while the mutation rate of cancer cells is obviously higher. (int binomial)

Mutations happen randomly in any random individual(cell) of the population in any random position. Each position of each individual can be mutated more than once. In code the mutation is represented by a "magic" number mask as the right-shift operator of 1 in the position of mutation, procedure based on bitwise operations. (int mutate)

With the functions unsigned int countSetBits, void setbits_positions and int* individual_setbits_positions we get respectively the number of total setbits(1s) for each individual(cell) and comulatively for each generation, so both inherited and new mutations are found, the positions that have been mutated for each generation and the mutated positions for each individual of each generation. 

The final attempt for finding a ftiness function was to give gravity and "mutation significancy" in certain positions of the genome. The number of positions that matter is random. Each genome position is controlled by a "fit_factor" based again on bitwise operations. The scores are given in the function typ fitfactors_func which returns the total "fit_factor_score", a bitwise number which represents the positions that matter. This number will serve as control in the function void com_fitness_func2, in order to check whether we have mutation in position that matters or not, consistently whether we have SMALL fitness (1.0) or LARGER fitness (described above).

Όσον αφορά το αρχείο pop_sim_interactions_with_gauss_weights.c

1) Eβαλα βαρη σε καθε θεση απο τις 64 απο την γκαουσιανη κατανομη ( κανονικη) γτ συμφωνα με τον νομο των μεγαλων αριθμων ολα τεινουν ν ακολουθουν κανονικη κατανομη και αυτο ισχυει και για βιολογικα συστηματα.

2) Eτσι εφτιαξα το weight_position_array. (weight_position_array_func)

3) Αλλά με αυτο θεωρειτει η καθε θεση ανεξαρτητη ( δλδ δεν επηρεαζεται η μια απο την αλλη), απλό γραμμικό μοντέλο.

4) Εφοσον ειναι απο κανονικη κατανομη αυτα τα βαρη εχουν και θετικες και αρνητικες τιμες γτ ηθελα να πιασω και το ενδεχομενο μια θεση να επηρεαζει αρνητικα την αλλη.

5) Μετα εφτιαξα το new_weight_position_array το οποιο καθε στοιχειο του βγαινει σαν συνδυασμος 2* το αντιστοιχο βαρος + το αθροισμα ολων των προηγουμενων βαρων και ολων των επομενων. (weight_selection_position1, sum_array))

6) Άρα καθε θεση πια επηρεαζεται περισσοτερο απο τον εαυτο της ( λογικο) γι αυτο και ο συντελεστης 2* αλλα και απ ολες τις προηγουμενες και απ ολες τις επομενες.

7) Άρα εχουμε παει σε πολυμεταβλητο γραμμικο μοντελο.

8) Αυτα τωρα τα νεα σκορ θα ειναι τα fitness scores.

9) Ωστοσο λογω της κανονικης κατανομης θα προκυψουν και αρνητικα βαρη αρα και αρνητικα scores....και αυτο ειναι θεμα γτ το selection δεν μπορει να λειτουργησει με αρνητικα fitnesses, μόνο με θετικα γτ ουσιαστικα fitness ειναι η ικανοτητα να δωσεις παιδια....βασικα ειναι μια πιθανοτητα να δωσεις παιδια....αρα δεν μπορει να ειναι αρνητικη μονο μεγαλυτερη ή ιση του 0 και επισης η δινεις παιδια μεγαλυτερο απο 0 η δε δινεις 0.

10) Άρα γι αυτο μου εκανε abort. Άρα το selection λειτουργουσε σωστα.

11) Η μια εκδοχη ηταν ν αφαιρεσω τ αρνητικα βαρη αλλά δεν ήθελα γτ ήθελα και η δεύτερη εκδοχή ήταν  οπου βρισκω fitness <0  με μια if να το θετω 0, όμως σε αυτό υπήρχε το θέμα ότι οταν το small fitness το βαζω 1 για τις περιοχες που δεν μετρανε τοσο πολυ θα βγει μικροτερο απο το 1 οποτε εντελει θα εκγαθιδτυθουν περιοχες που δεν μετράνε.

12) Μετά από πολύ ψάξιμο βρήκα ότι το roulette wheel selection algorithm λειτουργεί και για αρνητικά fitness scores εφόσον πρώτα κανονικοποιηθούν όλα τα fitnesses.

13) Kαι τοτε ναι θα διατηρουσα και τ αρνητικα βαρη και εντελει ολα τα fitness scores θα επαιρναν τιμες απο το 0 και πανω αλλα εντελει η διαφορα ( που αυτη ουσιαστικα μας ενδιαφερει για το selection) θα παραμεινει ιδια.

14) και ετσι δημιουργησα το norm_new_weight_position_array καθε στοιχειο του οποιου ειναι το new_weight_position_array-minimum_element_of_new_weight_position_array. ( weight_selection_position1,find_min)

In function void com_fitness_func2 each individual of each generation is evaluated with a certain fitness score given as the sum of the weights of its positions that have been mutated.

For each generation the kids are identical to their parents become parents of the next generation, based on the idea of the identical division of the cell cycle. The probability for an individual to become parent is proportional to its fitness.
In function void selection_func we have implemeted the Rhoulette wheel selection algorithm: 
Calculate S = the sum of a finesses (comulative_fitness_array)
For each individual of the second population,generate a random number between 0 and S.
Starting from the top of the first population, keep adding the finesses to the partial sum P (com_fitness[k]), till P<S.
The individual for which P exceeds S is the chosen individual to become parent.
The selection_func also returns the growth_rate of each generation as the number of proliferative cells/ total number of cells(prolferative and non proliferative), with the help of the function int distinct_elements.

void print_genome_without , int print_genome_with are for printing each genotype as a binary string, with the help of char * bitfunc.

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




Να σημειωθεί ότι όλα τα txt και png αρχεία έχουν παραχθεί από τον κώδικα pop_sim_interactions_with_gauss_weights.c καθώς θεωρώ ότι μέχρι στιγμής αυτός ειναι ο σωστότερος κώδικας.
