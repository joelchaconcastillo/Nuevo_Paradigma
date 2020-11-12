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
void diff_evo_xoverA(CIndividual &ind0, CIndividual &ind1, CIndividual &ind2, CIndividual &ind3, CIndividual &child, double CR, double F)
{
	// Check Whether the cross-over to be performed
	/*Loop over no of variables*/
	int idx_rnd = rand()%(nvar*nInd);//int(rnd_uni(&rnd_uni_init)*nvar*nInd);

	for(int n=0;n<nvar*nInd;n++)
	{
	  double rnd = rnd_uni(&rnd_uni_init);
	  if(rnd<CR||n==idx_rnd)
		  child.x_var[n/nvar][n%nvar] = ind1.x_var[n/nvar][n%nvar] + F*(ind2.x_var[n/nvar][n%nvar] - ind3.x_var[n/nvar][n%nvar]);
	  else
		  child.x_var[n/nvar][n%nvar] = ind0.x_var[n/nvar][n%nvar];
	  if(child.x_var[n/nvar][n%nvar]<vlowBound[n%nvar]){
 	       child.x_var[n/nvar][n%nvar] = ind0.x_var[n/nvar][n%nvar];//vlowBound[n] + rnd*(ind0.x_var[n] - vlowBound[n]);
	  }
	  if(child.x_var[n/nvar][n%nvar]>vuppBound[n%nvar]){ 
	        child.x_var[n/nvar][n%nvar] = ind0.x_var[n/nvar][n%nvar];//vuppBound[n] - rnd*(vuppBound[n] - ind0.x_var[n]);
	  }
	  if(child.x_var[n/nvar][n%nvar]<vlowBound[n%nvar]) child.x_var[n/nvar][n%nvar] = vlowBound[n%nvar];
	  if(child.x_var[n/nvar][n%nvar]>vuppBound[n%nvar]) child.x_var[n/nvar][n%nvar] = vuppBound[n%nvar];
	}
}
void diff_evo_xoverA_exp(CIndividual &ind0, CIndividual &ind1, CIndividual &ind2, CIndividual &ind3, CIndividual &child, double CR, double F, vector<bool> &changed)
{
	// Check Whether the cross-over to be performed
	/*Loop over no of variables*/
	int n = rand()%(nvar*nInd);//int(rnd_uni(&rnd_uni_init)*nvar*nInd);
	child.x_var= ind0.x_var;
	child.y_obj= ind0.y_obj;
        vector<int> asg_1(nInd), asg_2(nInd), asg_3(nInd);
	vector<vector<double> > dist_matrix_1(nInd, vector<double> (nInd)), dist_matrix_2(nInd, vector<double> (nInd)), dist_matrix_3(nInd, vector<double> (nInd));


        for(int i = 0; i < nInd; i++)
        {
          for(int j = i; j < nInd; j++)
          {
    	    if(i==j) dist_matrix_1[i][j]=dist_matrix_2[i][j]= dist_matrix_3[i][j]=DBL_MAX;
    	     else
	     {
	         dist_matrix_1[j][i] = dist_matrix_1[i][j] = distance_obj(ind0.y_obj[i], ind1.y_obj[j]);
	         dist_matrix_2[j][i] = dist_matrix_2[i][j] = distance_obj(ind0.y_obj[i], ind2.y_obj[j]);
	         dist_matrix_3[j][i] = dist_matrix_3[i][j] = distance_obj(ind0.y_obj[i], ind3.y_obj[j]);
	     }
          }
        } 
        KuhnMunkres(asg_1, dist_matrix_1);
        KuhnMunkres(asg_2, dist_matrix_2);
        KuhnMunkres(asg_3, dist_matrix_3);

	int cont =0;
	do{
         child.x_var[n/nvar][n%nvar] = ind1.x_var[asg_1[n/nvar]][n%nvar] + F*(ind2.x_var[asg_2[n/nvar]][n%nvar] - ind3.x_var[asg_3[n/nvar]][n%nvar]);
	  if(child.x_var[n/nvar][n%nvar]<vlowBound[n%nvar])
 	       child.x_var[n/nvar][n%nvar] = ind0.x_var[n/nvar][n%nvar];//vlowBound[n] + rnd*(ind0.x_var[n] - vlowBound[n]);
	  if(child.x_var[n/nvar][n%nvar]>vuppBound[n%nvar])
	        child.x_var[n/nvar][n%nvar] = ind0.x_var[n/nvar][n%nvar];//vuppBound[n] - rnd*(vuppBound[n] - ind0.x_var[n]);
	  if(child.x_var[n/nvar][n%nvar]<vlowBound[n%nvar]) child.x_var[n/nvar][n%nvar] = vlowBound[n%nvar];
	  if(child.x_var[n/nvar][n%nvar]>vuppBound[n%nvar]) child.x_var[n/nvar][n%nvar] = vuppBound[n%nvar];
	   changed[n/nvar]=true;
	   n++;
	   n %= (nvar*nInd);
	   cont++;
	}
        while(rnd_uni(&rnd_uni_init) < CR && cont < (nvar*nInd) );
}
void diff_evo_xoverA_exp(CIndividual &ind0, CIndividual &ind1, CIndividual &ind2, CIndividual &ind3, CIndividual &child, double CR, double F)
{
	// Check Whether the cross-over to be performed
	/*Loop over no of variables*/
	int n = rand()%(nvar*nInd);//int(rnd_uni(&rnd_uni_init)*nvar*nInd);
	child.x_var = ind0.x_var;
	int cont =0;
	do{
         child.x_var[n/nvar][n%nvar] = ind1.x_var[n/nvar][n%nvar] + F*(ind2.x_var[n/nvar][n%nvar] - ind3.x_var[n/nvar][n%nvar]);
	  if(child.x_var[n/nvar][n%nvar]<vlowBound[n%nvar])
 	       child.x_var[n/nvar][n%nvar] = ind0.x_var[n/nvar][n%nvar];//vlowBound[n] + rnd*(ind0.x_var[n] - vlowBound[n]);
	  if(child.x_var[n/nvar][n%nvar]>vuppBound[n%nvar])
	        child.x_var[n/nvar][n%nvar] = ind0.x_var[n/nvar][n%nvar];//vuppBound[n] - rnd*(vuppBound[n] - ind0.x_var[n]);
	  if(child.x_var[n/nvar][n%nvar]<vlowBound[n%nvar]) child.x_var[n/nvar][n%nvar] = vlowBound[n%nvar];
	  if(child.x_var[n/nvar][n%nvar]>vuppBound[n%nvar]) child.x_var[n/nvar][n%nvar] = vuppBound[n%nvar];
	   n++;
	   n %= (nvar*nInd);
	   cont++;
	}
        while(rnd_uni(&rnd_uni_init) < CR && cont < (nvar*nInd) );
}
double dist_obj(vector<double> &a, vector<double> &b)
{
	double d = 0.0;
	for(int i = 0; i < a.size(); i++) d += (a[i]-b[i])*(a[i]-b[i]);
	return d;
}
void diff_evo_xoverA_knn(CIndividual &ind0, CIndividual &ind1, CIndividual &ind2, CIndividual &ind3, CIndividual &child, double CR, double F)
{
       vector<int> s1(nInd),s2(nInd),s3(nInd);
       for(int k = 0; k < nInd; k++)
       {
	       vector<pair<double, int> > mind(3, make_pair(DBL_MAX, -1));
	   for(int i = 0; i < nInd; i++)
	   {
		   if(mind[0].first > dist_obj(ind0.y_obj[k], ind1.y_obj[i])) mind[0]=make_pair( dist_obj(ind0.y_obj[k], ind1.y_obj[i]), i);
		   if(mind[1].first > dist_obj(ind0.y_obj[k], ind2.y_obj[i])) mind[1]=make_pair( dist_obj(ind0.y_obj[k], ind2.y_obj[i]), i);
		   if(mind[2].first > dist_obj(ind0.y_obj[k], ind3.y_obj[i])) mind[2]=make_pair( dist_obj(ind0.y_obj[k], ind3.y_obj[i]), i);
	   }
	  // random_shuffle(mind.begin(), mind.end());
       s1[k]=mind[0].second;
       s2[k]=mind[1].second;
       s3[k]=mind[2].second;

       } 

	// Check Whether the cross-over to be performed
	/*Loop over no of variables*/
	int idx_rnd = rand()%(nvar*nInd);//int(rnd_uni(&rnd_uni_init)*nvar*nInd);

	for(int n=0;n<nvar*nInd;n++)
	{
	  double rnd = rnd_uni(&rnd_uni_init);
	  if(rnd<CR||n==idx_rnd)
		  child.x_var[n/nvar][n%nvar] = ind1.x_var[s1[n/nvar]][n%nvar] + F*(ind2.x_var[s2[n/nvar]][n%nvar] - ind3.x_var[s3[n/nvar]][n%nvar]);
	  else
		  child.x_var[n/nvar][n%nvar] = ind0.x_var[n/nvar][n%nvar];
	  if(child.x_var[n/nvar][n%nvar]<vlowBound[n%nvar]){
 	       child.x_var[n/nvar][n%nvar] = ind0.x_var[n/nvar][n%nvar];//vlowBound[n] + rnd*(ind0.x_var[n] - vlowBound[n]);
	  }
	  if(child.x_var[n/nvar][n%nvar]>vuppBound[n%nvar]){ 
	        child.x_var[n/nvar][n%nvar] = ind0.x_var[n/nvar][n%nvar];//vuppBound[n] - rnd*(vuppBound[n] - ind0.x_var[n]);
	  }
	  if(child.x_var[n/nvar][n%nvar]<vlowBound[n%nvar]) child.x_var[n/nvar][n%nvar] = vlowBound[n%nvar];
	  if(child.x_var[n/nvar][n%nvar]>vuppBound[n%nvar]) child.x_var[n/nvar][n%nvar] = vuppBound[n%nvar];
	}
}
#endif
