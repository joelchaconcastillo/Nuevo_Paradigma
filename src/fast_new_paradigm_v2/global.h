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
double cost_1[MAX_VAR][MAX_VAR], cost_2[MAX_VAR][MAX_VAR], cost_3[MAX_VAR][MAX_VAR];
int asg_1[MAX_VAR], asg_2[MAX_VAR], asg_3[MAX_VAR];
////
double distance_obj(vector<double> &a, vector<double> &b)
{
   double dist =0.0;
   for(int i = 0; i < a.size(); i++)
     dist += ((a[i]-b[i])*(a[i]-b[i]));
   return dist;
}
double fitnessfunction(vector <double> &y_obj, double *namda)
{
    // ASF
	double max_fun = -1.0e+30;
        double sum = 0.0;
	for(int n=0; n<nobj; n++)
	{
		double diff = fabs(y_obj[n] - idealpoint[n]);
		double feval;
		if(namda[n]==0) 
			//feval = 0.0001*diff;
//			feval = 0.0001*diff;
			feval = diff/0.0001;
		else
			//feval = diff*namda[n];
			feval = diff/namda[n];
	        sum +=feval;
		if(feval>max_fun) max_fun = feval;
	}
	return max_fun;// + 0.1*sum;
}
#endif
