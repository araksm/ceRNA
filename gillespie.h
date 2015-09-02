//
//  gillespie.h
//
//
//  Created by Araks  Martirosyan on 12/12/13.
//
//

#ifndef ____gillespie__
#define ____gillespie__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "mt64.h"
#include "rando3.h"

#endif /* defined(____gillespie__) */

float Tstop = 30000;
float dt = 1;
long m1 = 0, m2 = 0, mu = 0, c1 = 0, c2 = 0, p2 = 0;


float meanMu = 0;
float meanM1 = 0;
float meanM2 = 0;
float maxbound = 0;
float minbound = 10000;
float minVariance2 = 10000;


//Gillespie Algorithm

int gillespie ( float * params, int write, FILE * pFile, FILE* myAvfile){
    
    float b1,  b2,  bmu,  d1,  d2,  dmu,  sigma1,  sigma2,  k1pl,  k2pl,  k1min, k2min, kappa1, kappa2,  n1var,  n2var,  nmuvar,  h,  kin, kout, tfmax, tfmin, m2min, m2max, channel, step, fvar;
    
    // reading parameters
    b1 = params[0];
    b2 = params[1];
    bmu = params[2];
    d1 = params[3];
    d2 = params[4];
    dmu = params[5] ;
    sigma1 = params[6] ;
    sigma2 = params[7] ;
    k1pl = params[8] ;
    k2pl = params[9] ;
    k1min = params[10] ;
    k2min = params[11] ;
    kappa1 = params[12] ;
    kappa2 = params[13] ;
    n1var = params[14] ;
    n2var = params[15] ;
    nmuvar = params[16] ;
    h = params[17] ;
    kin = params[18] ;
    kout = params[19] ;
    tfmax = params[20];
    tfmin = params[21];
    m2min = params[22];
    m2max = params[23];
    channel = params[24];

    fvar = params[26];
    
    init_genrand64(7266447313870364031); // initializing random number generator
    
    
    
    float varianceM2 = 0;
    float varianceM1 = 0;
    float varianceM12 = 0;
    float varianceMu = 0;
    float meanC1 = 0;
    float meanC2 = 0;
    float fmwrite = 0;
    unsigned long meanN = 0;
    
    
    
    
    float epsilon = DBL_EPSILON + 0.000000001;
    int stop = 0, i = 0;
    float  a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a19, a20, omega;

    float nmu0 = 1;
    float nm10 = 1;
    float nm20 = 1;
    
    // choosing channel: 1 - miRNA channel, 2 - TF channel
    
    if(channel < 2){
        nm10 = (kin * pow(fvar, h) ) / (kin *pow(fvar, h) + kout);
        nm20 = n2var;
    }
    
    else
        if (channel > 1)  {
            nm10 = n1var;
            nm20 = (kin *pow(fvar, h) ) / (kin  * pow(fvar, h)  + kout);
        }
    
    
    nmu0 = nmuvar;
    
    
    float u = 0.0;
    float timeT = 0.01;
    float timeTflow = 0.0;
    float r1 = 0.0;
    
    // calculating the probabilities of events
    
    a1 = b1 * nm10 ;
    
    a2 = d1 * m1 ;
    a3 = k1pl * m1 * mu ;
    a4 = k1min * c1 ;
    
    a5 = b2 * nm20 ;
    a6 = d2 * m2 ;
    a7 = k2pl * m2 * mu ;
    a8 = k2min * c2 ;
    
    a9 = bmu * nmu0 ;
    a10 = dmu * mu ;
    a11 = kappa1 * c1 ;
    a12 = kappa2 * c2 ;
    
    a19 = sigma1 * c1 ;
    a20 = sigma2 * c2 ;
    
    omega = a1 + a2 + a3 + a4 + a5 + a6 + a7+ a8 + a9 + a10 + a11 + a12 + a19 + a20;
    
    
    while ( timeTflow < Tstop && !( fabs(omega) < epsilon )) { // start while loop
        
        //random number in [1, 0]
        r1 = genrand64_real1();
        
        //random number in [1, 0]
        u = genrand64_real1();
        
        // next reaction's time period
        timeT = (float) -log(u)/ omega;
        
        
        // choosing next reaction
        r1 = r1*omega;
        
        
        // ceRNA1
        if(0 <= r1  && r1  <= a1) {
            m1 ++;
            
        }
        
        else if ( r1   <= a1 + a2){
            
            m1 -- ;
        }
        //unbinding
        else if ( r1   <=  a1 + a2 + a3 ){
            
            m1 -- ;
            mu -- ;
            c1 ++ ;
            
        }
        else if (  r1   <= a1 + a2 + a3 + a4 ){
            
            m1 ++ ;
            mu ++ ;
            c1 -- ;
            
        }
        // ceRNA2
        else if (  r1   <= a1 + a2 + a3 + a4 + a5 ) {
          
            m2 ++ ;
            
        }
        
        else if (  r1   <= a1 + a2 + a3 + a4 + a5 + a6){
            
            m2 -- ;
            
        }
        // binding
        else if (   r1   <= a1 + a2 + a3 + a4 + a5 + a6 + a7){
           
            m2 -- ;
            mu -- ;
            c2 ++ ;
            
            
        }
        
        else if ( r1   <= a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 ){
            
            m2 ++;
            mu ++;
            c2 --;
            
        }
        // miRNA
        else if( r1   <= a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 ) {
            
            mu ++;
            
        }
        else if ( r1   <= a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 ){
            
            mu -- ;
            
        }
        // recycling
        else if(  r1   <= a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 ) {
            
            mu ++ ;
            c1 --;
            
        }
        else if ( r1  <= a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 ){
            
            mu ++;
            c2 --;
        }
        
        //c degradation
        else if( r1   <= a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a19 ) {
            
            c1 -- ;
        }
        
        else if ( r1   <= omega ){
            
            c2 -- ;
            
        }
        
        timeTflow += timeT;
        
        // recording
        if(timeTflow > dt * i && m1>=0 && m2>=0 && mu>=0 && c1>=0 && c2>=0 ) {
            if(write > 0) fprintf( pFile, "tau:\t%d\t%li\t%li\t%li\t%li\t%li\t%f\n", i, m1, m2, mu, c1, c2, fvar);
            
            meanM2 +=  m2;
            varianceM2 += m2*m2;
            fmwrite = fvar;
            meanM1 +=  m1;
            varianceM1 += m1*m1;
            meanMu +=  mu;
            varianceMu += mu *mu;
            varianceM12 += m1*m2;
            meanC1 += c1;
            meanC2 += c2;
            meanN ++;
            
            i++;
        }
        
        
        a2 = d1 * m1 ;
        a3 = k1pl * m1 * mu ;
        a4 = k1min * c1 ;
        
        a6 = d2 * m2 ;
        a7 = k2pl * m2 * mu ;
        a8 = k2min * c2 ;
        
        a10 = dmu * mu ;
        a11 = kappa1 * c1 ;
        a12 = kappa2 * c2 ;
        
        a19 = sigma1 * c1 ;
        a20 = sigma2 * c2 ;
        
        omega = a1 + a2 + a3 + a4 + a5 + a6 + a7+ a8 + a9 + a10 + a11 + a12 + a19 + a20;
    } //end of  while loop
    
    
    // calculating average
    meanM1 = meanM1/meanN;
    meanM2 = meanM2/meanN;
    meanMu = meanMu/meanN;
    meanC1 = meanC1/meanN;
    meanC2 = meanC2/meanN;
    
    //calculating variance
    varianceM2 = varianceM2/meanN - meanM2 * meanM2;
    varianceM1 = varianceM1/meanN - meanM1 * meanM1;
    varianceMu = varianceMu/meanN - meanMu * meanMu;
    varianceM12 = varianceM12/meanN - meanM1*meanM2;
    
    // writing to the file steady_state
    
    if(
       meanM1 >= 0 && varianceM1 >= 0 && meanM2 >= m2min  && meanM2 <= m2max  && varianceM2 >= 0 &&
       myAvfile != NULL ){
        
        
        fprintf( myAvfile, "state:\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", meanM1, varianceM1, meanM2, varianceM2, meanMu, varianceMu, meanC1, meanC2, fvar);
    
        // calculating m2min, variancem2 (m2min) and m2max
        if(meanM2 < minbound) {
            minVariance2 = varianceM2;
            minbound = meanM2;
        };
            if(meanM2 > maxbound) maxbound = meanM2;
        
        
    }
    
    

    return 0;
    
    
}
