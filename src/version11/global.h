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
#include "Kuhn_Munkres.h"

using namespace std;

#include "random.h"


//------------- Parameters in test instance ------------------

int     nvar,  nobj;                    //  the number of variables and objectives
int param_k, param_l;

double  vlowBound[3000] ,   vuppBound[3000];   //  lower and upper bounds of variables

char    strTestInstance[256];
char    strpath[800];


//------------- Parameters in random number ------------------
int     seed    = 177;
long    rnd_uni_init;        


//------------- Parameters in MOEA/D -------------------------

vector <double> idealpoint;
double          scale[100];  


int		etax    = 20, 	etam    = 50;   // distribution indexes of crossover and mutation

double  realx,  realm,  realb = 0.9;    // crossover, mutation, selection probabilities

double Di, Df, CR, F; // distance available of the hypersphere...
int nInd, nPop, nOffspring, nWeight;
long long max_nfes;

double distance_obj(vector<double> &a, vector<double> &b)
{
   double dist =0.0;
   for(int i = 0; i < a.size(); i++)
     dist += (a[i]-b[i])*(a[i]-b[i]);
   return dist;
}
#endif
