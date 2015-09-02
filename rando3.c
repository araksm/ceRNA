#include "mt64.h"
#include <math.h>

int rand_pm1()
{
  double r;
  r=genrand64_real1();
  if(r>0.5)
        return 1;
  else
    return (-1);
}

double drand_exp(double lam){
  double r;
  
  r = genrand64_real3();
  return -log(r)/lam;
}

/* --- Normal extraction --- */

double randG(double Mean,double Sigma)
{
  double v1,v2,gauss,r;
  do {
    v1=2.*genrand64_real2()-1.;
    v2=2.*genrand64_real2()-1.;
    r=(double)(v1*v1+v2*v2);
  } while (r>1);
  gauss=v1*sqrt(-2*log(r)/r);
  return (gauss*Sigma+Mean);
}







