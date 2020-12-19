#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#include "global.h"
#include "individual.h"

/* Routine for real polynomial mutation of an T */
void realmutation(CIndividual &child, double rate)
{
    long double rnd, delta1, delta2, mut_pow, deltaq;
    long double y, yl, yu, val, xy;
    long double eta_m = etam;

    for(int k =0; k < nInd;k++)
    {
	if(!child.changed[k])continue;
	int id_rnd = rand()%nvar;
        for (int j=0; j<nvar; j++)
        {
            if (rnd_uni<=rate)
            {
                y  = child.x_var[k][j];
                yl = vlowBound[j];
                yu = vuppBound[j];
                delta1 = (y-yl)/(yu-yl);
                delta2 = (yu-y)/(yu-yl);
                rnd = rnd_uni;
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
                child.x_var[k][j] = y;
            }
        }
    }
}
void diff_evo_xoverA(CIndividual &ind0, CIndividual &ind1, CIndividual &ind2, CIndividual &ind3, CIndividual &child, double CR, double F, vector<bool> &changed)
{
	child.x_var= ind0.x_var;
	child.y_obj= ind0.y_obj;
	// Check Whether the cross-over to be performed
	/*Loop over no of variables*/
	int idx_rnd = rand()%(nvar*nInd);//int(rnd_uni(&rnd_uni_init)*nvar*nInd);
        for(int i = 0; i < nInd; i++)
        {
          for(int j = 0; j < nInd; j++)
          {
	         cost_1[i*nInd+j] = -distance_obj(ind0.y_obj[i], ind1.y_obj[j]);
	         cost_2[i*nInd+j] = -distance_obj(ind0.y_obj[i], ind2.y_obj[j]);
	         cost_3[i*nInd+j] = -distance_obj(ind0.y_obj[i], ind3.y_obj[j]);
          }
        } 
       KM.hungarian(cost_1, asg_1);
       KM.hungarian(cost_2, asg_2);
       KM.hungarian(cost_3, asg_3);
	for(int n=0;n<nvar*nInd; n++)
	{
	  double rnd = rnd_uni;
	  if(rnd<CR||n==idx_rnd)
	  {
              child.x_var[n/nvar][n%nvar] = ind1.x_var[asg_1[n/nvar]][n%nvar] + F*(ind2.x_var[asg_2[n/nvar]][n%nvar] - ind3.x_var[asg_3[n/nvar]][n%nvar]);
	      if(child.x_var[n/nvar][n%nvar]<vlowBound[n%nvar])
 	           child.x_var[n/nvar][n%nvar] = ind0.x_var[n/nvar][n%nvar];//vlowBound[n] + rnd*(ind0.x_var[n] - vlowBound[n]);
	      if(child.x_var[n/nvar][n%nvar]>vuppBound[n%nvar])
	            child.x_var[n/nvar][n%nvar] = ind0.x_var[n/nvar][n%nvar];//vuppBound[n] - rnd*(vuppBound[n] - ind0.x_var[n]);
	      if(child.x_var[n/nvar][n%nvar]<vlowBound[n%nvar]) child.x_var[n/nvar][n%nvar] = vlowBound[n%nvar];
	      if(child.x_var[n/nvar][n%nvar]>vuppBound[n%nvar]) child.x_var[n/nvar][n%nvar] = vuppBound[n%nvar];
	       changed[n/nvar]=true;
	  }
	}
}
void diff_evo_xoverA_exp(CIndividual &ind0, CIndividual &ind1, CIndividual &ind2, CIndividual &ind3, CIndividual &child, double CR, double F)
{
	child.x_var= ind0.x_var;
	child.y_obj= ind0.y_obj;

	// Check Whether the cross-over to be performed
	/*Loop over no of variables*/
	int n = rand()%(nvar*nInd);//int(rnd_uni(&rnd_uni_init)*nvar*nInd);
        for(int i = 0; i < nInd; i++)
        {
          for(int j = 0; j < nInd; j++)
          {
	         cost_1[i*nInd+j] = -distance_obj(ind0.y_obj[i], ind1.y_obj[j]);
	         cost_2[i*nInd+j] = -distance_obj(ind0.y_obj[i], ind2.y_obj[j]);
	         cost_3[i*nInd+j] = -distance_obj(ind0.y_obj[i], ind3.y_obj[j]);
          }
        } 
       KM.hungarian(cost_1, asg_1);
       KM.hungarian(cost_2, asg_2);
       KM.hungarian(cost_3, asg_3);

	int cont =0;
	do{
         child.x_var[n/nvar][n%nvar] = ind1.x_var[asg_1[n/nvar]][n%nvar] + F*(ind2.x_var[asg_2[n/nvar]][n%nvar] - ind3.x_var[asg_3[n/nvar]][n%nvar]);
	  if(child.x_var[n/nvar][n%nvar]<vlowBound[n%nvar])
 	       child.x_var[n/nvar][n%nvar] = ind0.x_var[n/nvar][n%nvar];
	  if(child.x_var[n/nvar][n%nvar]>vuppBound[n%nvar])
	        child.x_var[n/nvar][n%nvar] = ind0.x_var[n/nvar][n%nvar];
	  if(child.x_var[n/nvar][n%nvar]<vlowBound[n%nvar]) child.x_var[n/nvar][n%nvar] = vlowBound[n%nvar];
	  if(child.x_var[n/nvar][n%nvar]>vuppBound[n%nvar]) child.x_var[n/nvar][n%nvar] = vuppBound[n%nvar];
	   child.changed[n/nvar]=true;
	   n++;
	   n %= (nvar*nInd);
	   cont++;
	}
        while(rnd_uni < CR && cont < (nvar*nInd) );
}
#endif
