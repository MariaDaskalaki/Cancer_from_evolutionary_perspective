# Cancer_from_evolutionary_perspective

Forward time simulations of cancer cells in discrete time of generations, based on two populations(previous cell population and current cell population).

The initial cell population consists of healthy cells.

The code is written in C and takes two command line parameters: initial number of cells and number of generations.

A "Frequency.txt" and a "Frequency_mutations.txt" file are derived from the C-code. 

Plots of txt files are derived from the corresponding python codes.


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
