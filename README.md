To run the program

**1. fill-up parameter values in the file *parameters* as follows**

    *b1:            ceRNA_1 synthesis rate 
    *b2:            ceRNA_2 synthesis rate
    *beta:          miRNA synthesis rate
    *d1:            ceRNA_1 degradation rate
    *d2:            ceRNA_2 degradation rate
    *delta:         miRNA degradation rate
    *sigma1:		complex_1 stoichiometric decay rate
    *sigma2:		complex_2 stoichiometric  decay rate
    *k1pl:          complex_1 association rate
    *k2pl:          complex_2 association rate
    *k1min:         complex_1 dissociation rate 
    *k2min:         complex_2 dissociation rate 
    *kappa1:		complex_1 catalytic decay rate
    *kappa2:		complex_2 catalytic decay rate
    *n1:            TF_1 occupancy
    *n2:            TF_2 occupancy
    *nmu:           TF_mu occupancy 
    *h:             number of cooperatively binding TF molecules
    *kin:           TF binding rate	
    *kout:          TF unbinding rate
    *fmax:          maximum number of input TF molecules
    *fmin:          minimum number of input TF molecules
    *m2min:         minimum size of the target population
    *m2max:         maximum size of the target population
    *channel:       assign 1 for the miRNA-channel and 2 for the TF-channel
    *step:          advancement of TF concentrations 
    *MC_steps:      number of Monte-Carlo steps
    *MC_eq:         number of Monte-Carlo steps to equilibrate the system



**2. make sure that you have g++ compiler**

**3. run the bash script *bash.sh***

