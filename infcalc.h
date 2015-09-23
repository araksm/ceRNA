#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <float.h>
#include <vector>
#include <iostream>     // std::cout


using std::vector;

#define _USE_MATH_DEFINES
float cycles = 500;




inline float log2mine(const float x){
    if(x>0)
        return  log(x) / log(2);
    else
        return 0;
}



// calculating optimal input distribution
float calculateOptInf(float lim, int oneortwo, char * meanfilename){
    
    FILE * pFile = fopen (meanfilename, "r");
    
    float fm1 = 0;
    
    float maxprob = 0;
    
    char s [50];
    
    float meanM2 = 0;
    float varianceM2 = 0;
    
    float meanM1 = 0;
    float varianceM1 = 0;
    
    float meanMu = 0;
    float varianceMu = 0;
    
    long meanN = 0;
    
    float meanC1 = 0;
    float meanC2 = 0;
    
    float inf1 = 0;
    float inf2 = 0;
    
    float mean1current = 0;
    float mean2current = 0;
    float f1mcurrent = 0;
    float var1current = 0;
    float var2current = 0;
    float var12current = 0;
    float var12 = 0;
    float temp1 = 0;
    float temp2 = 0;
    int i = 0;
    FILE * pFile1;
    
    pFile1 = fopen("prob_notNormalised", "w");
    
    
    if (pFile != NULL)
        while (fgetc(pFile) != EOF && fm1 < lim)
        {
            i++;
            
            fscanf(pFile, "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", s, &meanM1, &varianceM1, &meanM2, &varianceM2, &meanMu, &varianceMu, &meanC1, &meanC2, &fm1);
            
            if(i>0 && fm1 - f1mcurrent > 0 && varianceM1 >= 0 && varianceM2 > 0 )	{
                
                temp1 = (meanM1 - mean1current)/ (fm1 - f1mcurrent);
                temp1 =  temp1 * temp1;
                
                inf1+= sqrt(2* temp1/(varianceM1 + var1current));
                
                
                if(oneortwo == 1) {
                    fprintf(pFile1, "optProb:\t%f\t%f\n", fm1, sqrt(2* temp1/(varianceM1 + var1current)));
                    if(sqrt(2* temp1/(varianceM1 + var1current)) > maxprob && i > 1)
                    {
                        maxprob = sqrt(2* temp1/(varianceM1 + var1current));
                    }
                }
                temp2 = (meanM2 - mean2current)/ (fm1 - f1mcurrent);
                temp2 =  temp2 * temp2;
                
                inf2 += sqrt(2* temp2/(varianceM2 + var2current));
                
                if(oneortwo == 2) {
                    fprintf(pFile1, "optProb:\t%f\t%f\n", fm1, sqrt(2* temp2/(varianceM2 + var2current)));
                    
                    if(sqrt(2* temp2/(varianceM2 + var2current)) > maxprob && i>1)
                    {
                        maxprob = sqrt(2* temp2/(varianceM2 + var2current));
                    }
                }
                
                
            }
            
            mean1current = meanM1;
            mean2current = meanM2;
            var1current = varianceM1;
            var2current = varianceM2;
            var12current = var12;
            f1mcurrent = fm1;
            
        }
    inf1 /= sqrt(2* M_PI*exp(1));
    inf1 = log2mine(inf1);
    
    
    inf2 /= sqrt(2* M_PI*exp(1));
    inf2 = log2mine(inf2);
    
    
    // theoretical prediction
    //  if(oneortwo == 2){
    //     if(inf2> 0){
    //   printf("optinf2 : %f \t", inf2);
    ////     }
    //    else{
    //   printf("optinf2 : %f \t", 0.0);
    //    }
    // }
    
    fclose(pFile);
    fclose(pFile1);
    
    return  maxprob;
    
}


//calculating mutual information
float calculateInf2(long m1min, long m1max, long m2min, long m2max,   float varmin, float varmax, char * filename, float step){
    
    unsigned long mu = 0, m1 = 0, m2 = 0, mut = 0,  count = 0, length = 0, b1norm = 0, jointnorm1 = 0, jointnorm2 = 0, jointnorm12 = 0, numM1 = 0, numM2 = 0, numM12 = 0, c1 = 0, c2 = 0;
    float var = 0;
    int t = 0;
    
    float timeMu, time;
    float sumMu = 0;
    float inf1 = 0;
    float inf2 = 0;
    float inf3 = 0;
    
    sumMu = 0;
    
    unsigned long breaksB = (unsigned long) (varmax - varmin)/step;
    unsigned long breaksM2 = (unsigned long) (m2max - m2min)/step;
    
    
    vector<float> countm2 (breaksM2+1);
    
    vector<vector<float> > m2count;
    
    vector<float> b1count(breaksB+1);
    
    // Set up sizes. (HEIGHT x WIDTH)
    m2count.resize(breaksB+1);
    for (int i = 0; i <= breaksB; ++i)
        m2count[i].resize(breaksM2+1);
    
    
    float db1m = step;//(float) (b1max - b1min)/breaksB;
    float dm2 = step;//(float) (m2max - m2min)/breaksM;
    
    float dm2f = step; // (float) (m2max - m2min)/breaks;
    
    FILE * pFile;
    
    char s [50];
    
    long p = 0, q = 0;
    
    
    
    for(p = 0; p <= breaksM2; p++){
        countm2[p] = 0;
        
    }
    
    for(p = 0; p <= breaksB; p++){
        b1count[p] = 0;
        
        for(q = 0; q <= breaksM2; q++){
            m2count[p][q] = 0;
            
        }
    }
    /////////////////////////// inf by mu
    
    pFile = fopen ("alll", "r");
    
    if (pFile != NULL){
        
        //   printf("kuku%i\t", var<varmax );
        while (fgetc(pFile) != EOF && var<varmax)
        {
            
            fscanf(pFile, "%s\t%f\t%li\t%li\t%li\t%li\t%li\t%f\n", s, &timeMu, &m1, &m2, &mu, &c1, &c2, &var);
            
            //    printf("\t%li\n", m2 );
            
            
            for(p = 0; p <= breaksM2; p++ ) {
                
                if(m2 >= m2min + p*dm2 && m2 < m2min + (p+1)*dm2)
                {
                    countm2[p] ++;
                }
                
            }
            
            for(p = 0; p <= breaksB; p++ ) {
                
                
                if(var >= varmin + p*db1m && var < varmin + (p+1)*db1m ) {
                    
                    b1count[p] ++;
                    
                    for(q = 0; q <= breaksM2; q++){
                        if(m2 >= m2min + q*dm2 && m2 < m2min + (q+1)*dm2)
                        {
                            m2count[p][q] ++;
                        }
                    }
                    
                    
                }
                
            }
        }
    }
    
    
    fclose (pFile);
    
    
    for(p = 0; p <= breaksB; p++){
        b1norm += b1count[p];
    }
    
    
    
    for(p = 0; p <= breaksM2; p++){
        numM2 += countm2[p];
    }
    
    
    
    for(p = 0; p <= breaksB; p++){
        for(q = 0; q <= breaksM2; q++){
            jointnorm2 += m2count[p][q];
        }
    }
    
    inf2 = 0;
    for(p = 0; p <= breaksB; p++){
        for(q = 0; q <= breaksM2; q++){
            inf2 += (float) m2count[p][q] * (log2mine(m2count[p][q]) - log2mine(jointnorm2) - log2mine ( b1count[p]) + log2mine(b1norm) -  log2mine (countm2[q]) + log2mine(numM2))/jointnorm2;
            // printf("inf2: %f\t", inf2);
            
        }
    }
    
    
    
    
    if(inf2> 0)
    {
        printf("%f ", inf2);
    }
    else
    {
        printf("%f \t", 0.0);
    }
    
    return inf2;
    
}






