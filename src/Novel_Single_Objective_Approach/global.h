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

int nvar, nobj, param_k, param_l, etax = 20, etam = 50, nInd, nPop, nOffspring, nWeight, NP;
struct strindividual
{
   double *x_var, *y_obj, *ff;//fitness of each front
   bool *changed;
   int *f, *sf;//fronts, size-fronts
};
double *x_var, *y_obj, *idealpoint, *fitness;
bool *changed;
int *fronts, *size_fronts, *domin_to, *times_dominated, *size_domin_to, *hypermat_assig;
double *namda, *costs, *memo_dist;

set<pair<double, int> > *w_set;
double * contribution_R2;
bool *rejected;

long long max_nfes;
char    strTestInstance[256], strpath[800];
double  vlowBound[MAX_VAR] ,   vuppBound[MAX_VAR], realx,  realm,  realb = 0.9, Di, Df, CR, F;   //  lower and upper bounds of variables
Hungarian KM;

#endif
