#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "gillespie.h"
#include "infcalc.h"


const float Nmax = 200000;

///////// initialisation

float step = 1;

long MCsteps = 20000;
long MCeq = 1000;


// calculating m2 average and variance
void generateDataFromUD( float * params, int write, char * filename){
    
    
    FILE * pFile;
    FILE * pFile2;
    
    
    pFile = fopen (filename, "w");
    FILE * myAvFile = fopen ("steady_state", "w");
    
    params[28] = params[21];
    
    
    while(params[28] <= params[20] + 0.1){
        
        if( Nmax > m1 + m2 + c1 + c2 + mu){
            
            Tstop = MCeq;
            gillespie( params, 0, pFile, NULL); // equilibrating the system
            
            
            Tstop = MCsteps; // number of monte carlo steps
            gillespie( params, 1, pFile, myAvFile); // recording the steady state
            
        }
        params[28]  += step;
        
    }
    fclose(myAvFile);
    
    fclose(pFile);
    
    
}


// processing the optimal input
int generateDataFromGivenDist(float * params, int write, char * filename, int stopTime){
    
    
    char s [50];
    float fm1q = 0;
    float probabiliti = 0;
    
    FILE * pFile1 = fopen ("prob_notNormalised", "r");
    FILE * pFile = fopen ("alll", "w");
    
    
    int i = 0, u = 0;
    int m2maxx = 0;
    float normalizationFactor = 0;
    
    if (pFile1 != NULL)
        
        while (fgetc(pFile1) != EOF)
            
        {
            fscanf(pFile1, "%s\t%f\t%f\n", s, &params[28], &probabiliti );
            i++;
            if( Nmax >   m1 + m2 + c1 + c2 + mu && i >= step * u && i>1){
                //   printf("t%f\t%f\n", params[26], probabiliti);
                
                Tstop = MCeq;
                
                gillespie(params, 0, pFile, NULL); // equilibrating the system
                
                normalizationFactor += probabiliti;
                
                if(m2 > m2maxx) m2maxx = m2;
                
                Tstop = probabiliti * stopTime;
                gillespie(params, 1, pFile, NULL);  // recording the steady state
                u++;
            }
        }
    
    
    fclose(pFile);
    
    fclose(pFile1);
    
    return m2maxx;
    
    
    
}

/////////////////////////-------------- main


int main(){
    
    float fm10 = 1;
    srand (time(NULL));
    
    
    char filename [10] = "alll1";
    char meanfilename [15] = "steady_state";
    
    
    float b1=0,  b2=0,  bmu=0,  d1=0,  d2=0,  dmu=0,  sigma1=0,  sigma2=0,  k1pl=0,  k2pl=0,  k1min=0, k2min=0, kappa1=0, kappa2=0,  n1var=0,  n2var=0,  nmuvar=0,  h=0,  kin=0, kout=0, tfmax=0, tfmin=0, m2min=0, m2max=0, channel=0, fvar=0;
    
    
    float params [29] = { b1,  b2,  bmu,  d1,  d2,
        dmu,  sigma1,  sigma2,  k1pl,  k2pl,
        k1min, k2min, kappa1, kappa2,  n1var,
        n2var,  nmuvar,  h,  kin, kout,
        tfmax, tfmin, m2min, m2max, channel,
        step, MCsteps, MCeq, fvar};
    
    
    float maxprob = 0;
    
    char filename1 [10] = "alll";
    char s[20];
    
    FILE * paramsFile;
    
    paramsFile = fopen("parameters","r"); // read mode
    
    int i= 0;
    
    if (paramsFile != NULL){
        
        while (fgetc(paramsFile) != EOF)
        {
            
            fscanf(paramsFile, "%s\t%f\n", s, &params[i]);
            i ++;
            
        }
    }
    
    
    
    step = params[25];
    tfmax = params[20];
    MCsteps = (long) params[26];
    MCeq = (long) params[27];
    
    maxbound = 0;
    minbound = 10000;
    minVariance2 = 10000;
    
    
    generateDataFromUD(params, 0, filename); // generating the output
    maxprob = calculateOptInf(params[20], 2, meanfilename); // calculating optimal input
    
    
    m2max = generateDataFromGivenDist(params, 1, filename1, MCsteps/maxprob); // processing optimal input distribution
    
    printf("MI\t\t Delta m2 \t m2min \t\t variance(m2min)\n"); // displaying results
    
    calculateInf2( 0, b1/d1, 0, m2max,  1,  tfmax, filename1, step); // calculating MI
    
    
    printf("\t%f\t%f\t%f\n", maxbound - minbound, minbound, minVariance2 ); // displaying results
    
    
    
    return 0;
    
}
