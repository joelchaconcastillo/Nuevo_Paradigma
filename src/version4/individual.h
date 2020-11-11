#ifndef __INDIVIDUAL_H_
#define __INDIVIDUAL_H_

#include "global.h"
#include "problem.h"
class CIndividual{
    public:
	CIndividual();
	virtual ~CIndividual();

	vector<vector<double> > x_var;
	vector<vector<double> > y_obj;
        vector<bool> modified;
	int    count;

	void   rnd_init();
	void   obj_eval();
	void   show_objective();
	void   show_variable();
	vector<double>  fitness;
    bool   operator<(const CIndividual &ind2);
//    bool   operator<(const vector<double> &y_obj1, const vector<double> &y_obj2);
    bool   operator<<(const CIndividual &ind2);
    bool   operator==(const CIndividual &ind2);
    void   operator=(const CIndividual &ind2);
};
CIndividual::CIndividual()
{
	x_var.assign(nInd, vector<double>(nvar, 0));
        y_obj.assign(nInd, vector<double>(nobj, 0));
        modified.assign(nInd, false);
	fitness.assign(nInd, 0);
}
CIndividual::~CIndividual()
{
}
void CIndividual::rnd_init()
{
    for(int k= 0; k < nInd; k++)
    {
       for(int n = 0; n<nvar; n++)
       {
           x_var[k][n] = vlowBound[n] + rnd_uni(&rnd_uni_init)*(vuppBound[n] - vlowBound[n]);    
       }
    }
}
void CIndividual::obj_eval()
{
	vector<int> v;
	priority_queue<pair<double, int> > pq;
        for(int k = 0;k <nInd; k++)
	{
	   v.push_back(k);
	   if(!modified[k])continue;
	   if(!strcmp("UF1", strTestInstance))  CEC09_F1(y_obj[k], x_var[k]);
	   if(!strcmp("UF2", strTestInstance))  CEC09_F2(y_obj[k], x_var[k]);
	   if(!strcmp("UF3", strTestInstance))  CEC09_F3(y_obj[k], x_var[k]);
	   if(!strcmp("UF4", strTestInstance))  CEC09_F4(y_obj[k], x_var[k]);
	   if(!strcmp("UF5", strTestInstance))  CEC09_F5(y_obj[k], x_var[k]);
	   if(!strcmp("UF6", strTestInstance))  CEC09_F6(y_obj[k], x_var[k]);
	   if(!strcmp("UF7", strTestInstance))  CEC09_F7(y_obj[k], x_var[k]);
	   if(!strcmp("UF8", strTestInstance))  CEC09_F8(y_obj[k], x_var[k]);
	   if(!strcmp("UF9", strTestInstance))  CEC09_F9(y_obj[k], x_var[k]);
	   if(!strcmp("UF10", strTestInstance)) CEC09_F10(y_obj[k], x_var[k]);

	   //WFG test instances....
	   if(!strcmp("WFG1", strTestInstance))  wfg1(y_obj[k], x_var[k]);
	   if(!strcmp("WFG2", strTestInstance))  wfg2(y_obj[k], x_var[k]);
	   if(!strcmp("WFG3", strTestInstance))  wfg3(y_obj[k], x_var[k]);
	   if(!strcmp("WFG4", strTestInstance))  wfg4(y_obj[k], x_var[k]);
	   if(!strcmp("WFG5", strTestInstance))  wfg5(y_obj[k], x_var[k]);
	   if(!strcmp("WFG6", strTestInstance))  wfg6(y_obj[k], x_var[k]);
	   if(!strcmp("WFG7", strTestInstance))  wfg7(y_obj[k], x_var[k]);
	   if(!strcmp("WFG8", strTestInstance))  wfg8(y_obj[k], x_var[k]);
	   if(!strcmp("WFG9", strTestInstance))  wfg9(y_obj[k], x_var[k]);

           //DTLZ test instances....
	   if(!strcmp("DTLZ1", strTestInstance))  dtlz1(y_obj[k], x_var[k]);
	   if(!strcmp("DTLZ2", strTestInstance))  dtlz2(y_obj[k], x_var[k]);
	   if(!strcmp("DTLZ3", strTestInstance))  dtlz3(y_obj[k], x_var[k]);
	   if(!strcmp("DTLZ4", strTestInstance))  dtlz4(y_obj[k], x_var[k]);
	   if(!strcmp("DTLZ5", strTestInstance))  dtlz5(y_obj[k], x_var[k]);
	   if(!strcmp("DTLZ6", strTestInstance))  dtlz6(y_obj[k], x_var[k]);
	   if(!strcmp("DTLZ7", strTestInstance))  dtlz7(y_obj[k], x_var[k]);
//	   pq.push(make_pair(-y_obj[k][0], k));
        }
   vector<vector<double> > x=x_var, y=y_obj;
	random_shuffle(v.begin(), v.end());
        for(int k = 0;k <nInd; k++)
	{
///		int v = pq.top().second;
///		pq.pop();
///		x_var[k] = x[v];
///		y_obj[k] = y[v];
	   x_var[k]=x[v[k]];
	   y_obj[k]=y[v[k]];
	}
}
void CIndividual::show_objective()
{
    for(int n=0; n<nobj*nInd; n++)
		printf("%f ",y_obj[n/nInd][n%nobj]);
	printf("\n");
}
void CIndividual::show_variable()
{
    for(int k = 0; k < nInd; k++)
    for(int n=0; n<nvar; n++)
		printf("%f ",x_var[k][n]);
	printf("\n");
}
void CIndividual::operator=(const CIndividual &ind2)
{
	x_var = ind2.x_var;
	y_obj = ind2.y_obj;
	fitness = ind2.fitness;
	modified = ind2.modified;
}
bool CIndividual::operator<(const CIndividual &ind2)
{
    bool dominated = true;
    for(int n=0; n<nobj; n++)
    {
        if(ind2.y_obj[n]<y_obj[n]) return false;
    }
    if(ind2.y_obj==y_obj) return false;
    return dominated;
}
bool CIndividual::operator<<(const CIndividual &ind2)
{
	bool dominated = true;
//    for(int n=0; n<nobj; n++)
//	{
//		if(ind2.y_obj[n]<y_obj[n]  - 0.0001) return false;
//	}
//	if(ind2.y_obj==y_obj) return false;
	return dominated;
}
bool CIndividual::operator==(const CIndividual &ind2)
{
	if(ind2.y_obj==y_obj) return true;
	else return false;
}
bool operator<(const vector<double> &y_obj1, const vector<double> &y_obj2)
{
    bool dominated = true;
    for(int n=0; n<nobj; n++)
    {
        if(y_obj2[n]<y_obj1[n]) return false;
    }
    if(y_obj1==y_obj2) return false;
    return dominated;
}
bool operator<<(const vector<double> &y_obj1, const vector<double> &y_obj2)
{
    for(int n=0; n< y_obj1.size(); n++)
    {
//        if( fabs(y_obj1[n]-y_obj2[n]) < 10)
//	 return true;
	if(y_obj1[n]<y_obj2[n])
	 return true; 
        else if(y_obj1[n]>y_obj2[n])
	 return false;
    }
    return true;
}
#endif

