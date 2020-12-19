#ifndef __GLOBAL_H_
#define __GLOBAL_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
#include <cassert>
#include <algorithm>
#define MAX_VAR 2000
#include "Kuhn_Munkres.h"
#define rnd_uni (rand()/double(RAND_MAX))
using namespace std;

int     nvar,  nobj;                    //  the number of variables and objectives
int param_k, param_l;

double  vlowBound[MAX_VAR] ,   vuppBound[MAX_VAR];   //  lower and upper bounds of variables

char    strTestInstance[256];
char    strpath[800];
vector <double> idealpoint;
int		etax    = 20, 	etam    = 50;   // distribution indexes of crossover and mutation
double  realx,  realm,  realb = 0.9;    // crossover, mutation, selection probabilities
double Di, Df, CR, F; // distance available of the hypersphere...
int nInd, nPop, nOffspring, nWeight;
long long max_nfes;
double *namda;
////Hungarian information
double *cost_1, *cost_2, *cost_3;
int *asg_1, *asg_2, *asg_3;
Hungarian KM;
struct strIndividual
{
    vector<vector<double> > x_var, y_obj;
    vector<double>  fitness, changed;
};
bool operator<(const vector<double> &y_obj1, const vector<double> &y_obj2)
{
    for(int n=0; n<nobj; n++)
    {
        if(y_obj2[n]<y_obj1[n]) return false;
    }
    if(y_obj1==y_obj2) return false;
    return true;
}
bool operator<<(const vector<double> &y_obj1, const vector<double> &y_obj2)
{
    for(int n=0; n< y_obj1.size(); n++)
    {
	if(y_obj1[n]<y_obj2[n])
	 return true; 
        else if(y_obj1[n]>y_obj2[n])
	 return false;
    }
    return true;
}
double distance_obj(vector<double> &a, vector<double> &b)
{
   double dist =0.0;
   for(int i = 0; i < a.size(); i++)
     dist += ((a[i]-b[i])*(a[i]-b[i]));
   return dist;
}
double fitnessfunction(vector <double> &y_obj, double *namda)
{
	double max_fun = -1.0e+30;
	for(int n=0; n<nobj; n++)
	{
		double feval, diff = fabs(y_obj[n] - idealpoint[n]);
		if(namda[n]==0) feval = diff/0.0001;
		else feval = diff/namda[n];
		max_fun = max(max_fun, feval);
	}
	return max_fun;;
}
vector<set<int> > non_dominated_sorting(vector<vector<double> > &y_obj)
{
  vector<vector<int> > domin_to(y_obj.size());
  vector<int> times_dominated(y_obj.size(), 0);
  vector<set<int> > fronts(1);
   int current_rank = 0;
   for(int pidx1=0; pidx1 < y_obj.size(); pidx1++)
   {
      for(int pidx2=0; pidx2 < y_obj.size(); pidx2++)
      {
	if(pidx1 == pidx2) continue;
        if( y_obj[pidx1] < y_obj[pidx2]) domin_to[pidx1].push_back(pidx2);
 	else if( y_obj[pidx2] < y_obj[pidx1]) times_dominated[pidx1]++;
      }
      if( times_dominated[pidx1] == 0)
      {
         fronts[current_rank].insert(pidx1);
      }
  }
  //ranking.... 
  set<int> next_front;
  while(true)
  {
     for(auto i:fronts[current_rank])
     {
       for(auto j : domin_to[i])
	{
	   times_dominated[j]--;
	   if(times_dominated[j] == 0)
	   {
	      next_front.insert(j);
	   }
	}
     }
    if(next_front.empty())break;
    fronts.push_back(next_front);
    current_rank++;
    next_front.clear();
  }
  return fronts;
}
void eval_R2(vector<vector<double> > &y_obj, vector<double> fitness)
{
  vector<set<int> > fronts = non_dominated_sorting(y_obj);
  fitness.assign(nInd, 0);
  for(int r = 0; r < fronts.size(); r++)
  {
     for(int w = 0; w < nWeight; w++)
     {
       double minv = DBL_MAX;
       for(auto k:fronts[r])
       {
           minv = min(minv, fitnessfunction(y_obj[k], &namda[w*nobj]));
       } 
       fitness[r] += minv;
     }
  }
}
#endif
