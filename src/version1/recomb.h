#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#include "global.h"
#include "individual.h"

/* Routine for real polynomial mutation of an T */
void realmutation(vector<double> &ind, double rate)
{
    long double rnd, delta1, delta2, mut_pow, deltaq;
    long double y, yl, yu, val, xy;
    long double eta_m = etam;

	int id_rnd = rand()%nvar;

    for (int j=0; j<nvar; j++)
    {
        if (rnd_uni(&rnd_uni_init)<=rate)
        {
            y  = ind[j];
            yl = vlowBound[j];
            yu = vuppBound[j];
            delta1 = (y-yl)/(yu-yl);
            delta2 = (yu-y)/(yu-yl);
            rnd = rnd_uni(&rnd_uni_init);
            mut_pow = 1.0/(eta_m+1.0);
            if (rnd <= 0.5)
            {
                xy = 1.0-delta1;
                val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta_m+1.0)));
                deltaq =  pow(val,mut_pow) - 1.0;
            }
            else
            {
                xy = 1.0-delta2;
                val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta_m+1.0)));
                deltaq = 1.0 - (pow(val,mut_pow));
            }
            y = y + deltaq*(yu-yl);
            if (y<yl)
                y = yl;
            if (y>yu)
                y = yu;
            ind[j] = y;
        }
    }
    return;
}

void diff_evo_xoverA(vector<double> &ind0, vector<double> &ind1, vector<double> &ind2, vector<double> &ind3, vector<double> &child, double CR, double F)
{
	// Check Whether the cross-over to be performed
	/*Loop over no of variables*/
	int idx_rnd = rand()%(nvar);//int(rnd_uni(&rnd_uni_init)*nvar*nInd);
	for(int n=0;n<nvar;n++)
	{
	  double rnd = rnd_uni(&rnd_uni_init);
	  if(rnd<CR||n==idx_rnd)
		  child[n] = ind1[n] + F*(ind2[n] - ind3[n]);
	  else
		  child[n] = ind0[n];

	  if(child[n]<vlowBound[n]){
 	       child[n] = ind0[n];//vlowBound[n] + rnd*(ind0.x_var[n] - vlowBound[n]);
	  }
	  if(child[n]>vuppBound[n]){ 
	        child[n] = ind0[n];//vuppBound[n] - rnd*(vuppBound[n] - ind0.x_var[n]);
	  }
	  if(child[n]<vlowBound[n]) child[n] = vlowBound[n];
	  if(child[n]>vuppBound[n]) child[n] = vuppBound[n];
	}
}
#endif
