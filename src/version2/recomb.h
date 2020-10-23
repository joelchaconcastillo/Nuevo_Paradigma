#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#include "global.h"
#include "individual.h"

/* Routine for real polynomial mutation of an T */
void realmutation(CIndividual &ind, double rate)
{
    long double rnd, delta1, delta2, mut_pow, deltaq;
    long double y, yl, yu, val, xy;
    long double eta_m = etam;

	int id_rnd = int(rnd_uni(&rnd_uni_init)*nvar*nInd);

    for (int j=0; j<nvar*nInd; j++)
    {
        if (rnd_uni(&rnd_uni_init)<=rate)
        {
            y  = ind.x_var[j];
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
            ind.x_var[j] = y;
        }
    }
    return;
}


void diff_evo_xoverA(CIndividual &ind0, CIndividual &ind1, CIndividual &ind2, CIndividual &ind3, CIndividual &child, double CR, double F)
{
	// Check Whether the cross-over to be performed
	/*Loop over no of variables*/
	//int idx_rnd = rand()%(nvar*nInd);//int(rnd_uni(&rnd_uni_init)*nvar*nInd);
   for(int ids=0; ids<nInd;ids++)
   {
	int idx_rnd = rand()%(nvar);//int(rnd_uni(&rnd_uni_init)*nvar*nInd);
   //	int ids=rand()%nInd;
	for(int n=0;n<nvar;n++)
	{
	  double rnd = rnd_uni(&rnd_uni_init);
	  if(rnd<CR||n==idx_rnd)
		  child.x_var[ids*nvar + n] = ind1.x_var[ids*nvar + n] + F*(ind2.x_var[ids*nvar + n] - ind3.x_var[ids*nvar + n]);
	  else
		  child.x_var[ids*nvar + n] = ind0.x_var[ids*nvar + n];

	if(child.x_var[ids*nvar + n]<vlowBound[ids*nvar + n]){
 	       child.x_var[ids*nvar + n] = ind0.x_var[ids*nvar + n];//vlowBound[n] + rnd*(ind0.x_var[n] - vlowBound[n]);
	  }
	  if(child.x_var[n]>vuppBound[n]){ 
	        child.x_var[ids*nvar + n] = ind0.x_var[ids*nvar + n];//vuppBound[n] - rnd*(vuppBound[n] - ind0.x_var[n]);
	  }
	  if(child.x_var[ids*nvar + n]<vlowBound[ids*nvar + n]) child.x_var[ids*nvar + n] = vlowBound[ids*nvar + n];
	  if(child.x_var[ids*nvar + n]>vuppBound[ids*nvar + n]) child.x_var[ids*nvar + n] = vuppBound[ids*nvar + n];
	}
   }
}

//diff_evo_xoverB
void diff_evo_xoverB(CIndividual &ind0, CIndividual &ind1, CIndividual &ind2, CIndividual &child, double rate)
{
	int idx_rnd = int(rnd_uni(&rnd_uni_init)*nvar);


    double CR   =  (rnd_uni(&rnd_uni_init)<0.5)?0.2:1.0;
	for(int n=0;n<nvar;n++)
	{
	  /*Selected Two Parents*/

	  // strategy one 
	  // child.x_var[n] = ind0.x_var[n] + rate*(ind2.x_var[n] - ind1.x_var[n]);
	  
	  //*
	  // strategy two

	  double rnd1 = rnd_uni(&rnd_uni_init);
	  //double CR   = 1.0;
	  if(rnd1<CR||n==idx_rnd)
		  child.x_var[n] = ind0.x_var[n] + (rate)*(ind2.x_var[n] - ind1.x_var[n]);
	  else
		  child.x_var[n] = ind0.x_var[n];
	  //*/

	  // handle the boundary voilation
	  if(child.x_var[n]<vlowBound[n]){
	          double rnd = rnd_uni(&rnd_uni_init);
//	          double rnd =-0.1+1.2*rnd_uni(&rnd_uni_init);
 	       //child.x_var[n] = vlowBound[n] + rnd*(vuppBound[n] - vlowBound[n]);
 	        child.x_var[n] = ind0.x_var[n];// vlowBound[n] + rnd*(ind0.x_var[n] - vlowBound[n]);
	  }
	  if(child.x_var[n]>vuppBound[n]){ 
	          double rnd = rnd_uni(&rnd_uni_init);
	          //double rnd =-0.1+1.2*rnd_uni(&rnd_uni_init);
 	       // child.x_var[n] = vlowBound[n] + rnd*(vuppBound[n] - vlowBound[n]);
	        child.x_var[n] = ind0.x_var[n];//vuppBound[n] - rnd*(vuppBound[n] - ind0.x_var[n]);
	  }
	  if(child.x_var[n]<vlowBound[n]) child.x_var[n] = vlowBound[n];
	  if(child.x_var[n]>vuppBound[n]) child.x_var[n] = vuppBound[n];
	}
}
void diff_evo_xover2B(CIndividual &ind0, CIndividual &ind1, CIndividual &ind2, CIndividual &child, double rate, CIndividual &best)
{
	int idx_rnd = int(rnd_uni(&rnd_uni_init)*nvar);


    double CR   =  (rnd_uni(&rnd_uni_init)<0.5)?0.2:1.0;
	for(int n=0;n<nvar;n++)
	{
	  /*Selected Two Parents*/

	  // strategy one 
	  // child.x_var[n] = ind0.x_var[n] + rate*(ind2.x_var[n] - ind1.x_var[n]);
	  
	  //*
	  // strategy two

	  double rnd1 = rnd_uni(&rnd_uni_init);
	  //double CR   = 1.0;
	  if(rnd1<CR||n==idx_rnd)
		  child.x_var[n] = ind0.x_var[n] + (rate)*0.5*(ind2.x_var[n] - ind1.x_var[n]) + rate*0.5*(ind0.x_var[n] - best.x_var[n]);
	  else
		  child.x_var[n] = ind0.x_var[n];
	  //*/

	  // handle the boundary voilation
	  if(child.x_var[n]<vlowBound[n]){
	          double rnd = rnd_uni(&rnd_uni_init);
//	          double rnd =-0.1+1.2*rnd_uni(&rnd_uni_init);
 	       //child.x_var[n] = vlowBound[n] + rnd*(vuppBound[n] - vlowBound[n]);
 	        child.x_var[n] = ind0.x_var[n];// vlowBound[n] + rnd*(ind0.x_var[n] - vlowBound[n]);
	  }
	  if(child.x_var[n]>vuppBound[n]){ 
	          double rnd = rnd_uni(&rnd_uni_init);
	          //double rnd =-0.1+1.2*rnd_uni(&rnd_uni_init);
 	       // child.x_var[n] = vlowBound[n] + rnd*(vuppBound[n] - vlowBound[n]);
	        child.x_var[n] = ind0.x_var[n];//vuppBound[n] - rnd*(vuppBound[n] - ind0.x_var[n]);
	  }
	  if(child.x_var[n]<vlowBound[n]) child.x_var[n] = vlowBound[n];
	  if(child.x_var[n]>vuppBound[n]) child.x_var[n] = vuppBound[n];
	}
}
void diff_evo_xoverC(CIndividual &ind0, CIndividual &ind1, CIndividual &ind2, vector<double> &xdiff,  CIndividual &child,  double rate)
{
      double rnd = rnd_uni(&rnd_uni_init), rnd2 = rnd_uni(&rnd_uni_init);
	  for(int n=0;n<nvar;n++)
	  {
		  /*Selected Two Parents*/
		  
		  if(rnd<1)
		      child.x_var[n] = ind0.x_var[n] + rate*(ind2.x_var[n] - ind1.x_var[n]);
		  else
			  child.x_var[n] = ind0.x_var[n] + rnd2*xdiff[n];
	
		  if(child.x_var[n]<vlowBound[n]) child.x_var[n] = vlowBound[n];
		  if(child.x_var[n]>vuppBound[n]) child.x_var[n] = vuppBound[n];
	  }
}



#endif
