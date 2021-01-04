#ifndef __EVOLUTION_H_
#define __EVOLUTION_H_

#include <cfloat>
#include <set>
#include <queue>
#include <map>
#include <unordered_set>
#include <iomanip>
#include <iterator>
#include "global.h"
#include "problem.h"

class CMOEAD
{

   public:
	CMOEAD();
	~CMOEAD();

	void init_population();                
	void obj_eval(double *x_var, double *y_obj);
	inline void update_reference(double *point){for(int n=0; n<nobj; n++) if(point[n]<idealpoint[n]) idealpoint[n] = point[n];}
	void replacement_phase();
	void evol_population();                                    
	void exec_emo(int run);
	void save_front(char savefilename[4024]); 
	void save_pos(char savefilename[4024]);
	void update_parameterD();
        void update_archive();
	double distance_var( int a, int b);
	inline int* pointer_hyp(int a, int b){  return hypermat_assig +a*(nPop+nOffspring)*nInd + b*nInd; }
	inline double* pointer_dist(int a, int b){ if(a > b) swap(a, b);  return memo_dist + a*(nPop+nOffspring) + b; }

	void realmutation(double *x_var, double rate);
        void diff_evo_xoverA_exp(double *ind0, double *ind1, double *ind2, double *ind3, double *child, double CR, double F, int *asg_1, int *asg_2, int *asg_3, bool *childchanged);
	double distance_obj(double *a, double *b);
	double fitnessfunction(double *y_obj, double *namda);
	void non_dominated_sorting(double *y_obj, int *f, int *sf, int npoints);
        bool compareR2(int a, int b);
        void eval_R2(double *y_obj, int *front, int size_front, double &fitness);
	bool dominate(double *y_obj1, double *y_obj2);
	void add_archive(double *y_obj, double *x_var);
   private:
        
	int* child_idx, *parent_idx, *f_archive, *sf_archive; 
        strindividual *pool;
	double *archive_y, *archive_x, *cpy_archive_y, *cpy_archive_x; //archive...
	long long nfes;          //  the number of function evluations
	double	D;	//Current minimum distance
	int n_archive, max_archive, size_archive;
};
CMOEAD::CMOEAD()
{

}

void CMOEAD::update_parameterD()
{
      double TElapsed = nfes, TEnd = max_nfes;
      D = Di - Di * (TElapsed / (TEnd*Df));
      D = max(D, 0.0);
}
double CMOEAD::distance_var(int a, int b)
{
   double *distab=pointer_dist(a,b);
   if(*distab > 0.0) return *distab;
   int *asg_ab = pointer_hyp(a,b), *asg_ba = pointer_hyp(b,a);
   if(*asg_ab==-1)
   {
     for(int i = 0; i < nInd; i++)
      for(int j = 0; j < nInd; j++)
       costs[i*nInd+j] = -distance_obj(pool[a].y_obj+i*nobj, pool[b].y_obj+j*nobj);
     KM.hungarian(costs, asg_ab);
     for(int i = 0; i < nInd; i++) asg_ba[asg_ab[i]] = i;
   }
   double dist = 0.0;
   for(int i = 0; i < nInd; i++)
   {
      for(int j = 0; j < nvar; j++)
      {
         double factor = (pool[a].x_var[nvar*i+j]-pool[b].x_var[asg_ab[i]*nvar+j])/(vuppBound[j]-vlowBound[j]);
         dist += factor*factor;
      }
   }
   (*distab) = sqrt(dist);
   return sqrt(dist);
}
CMOEAD::~CMOEAD()
{
   delete [] pool;
   delete [] x_var;
   delete [] y_obj;
   delete [] idealpoint;
   delete [] fitness;
   delete [] changed;

   delete [] fronts; 
   delete [] size_fronts;
   delete [] domin_to;
   delete [] times_dominated;
   delete [] size_domin_to;

   delete [] namda;
   delete [] costs;
   delete [] hypermat_assig;
   delete [] memo_dist;
   
   delete [] parent_idx;
   delete [] child_idx;
   
   delete [] archive_y;
   delete [] cpy_archive_y;
   delete [] archive_x;
   delete [] cpy_archive_x;
   delete [] f_archive;
   delete [] sf_archive;
   delete [] contribution_R2;
   delete [] w_set;
   delete [] rejected;
}
void CMOEAD::init_population()
{
    NP = nOffspring + nPop;
    char filename[1024];
    size_archive = 100;
    max_archive = 2*size_archive;
    n_archive =  0;
    // Read weight vectors from a data file
    sprintf(filename,"%s/ParameterSetting/Weight/W%dD_%d.dat", strpath, nobj, nWeight);
    std::ifstream readf(filename);
    KM.init(nInd);
    pool = new strindividual[NP];
    parent_idx = new int[nPop];
    child_idx = new int[nOffspring];
    x_var = new double[nvar*nInd*NP];
    y_obj = new double[nobj*nInd*NP];
    idealpoint = new double[nobj]; 
    fitness = new double[NP*nInd];
    changed = new bool[NP*nInd];

    fronts = new int[NP*nInd*nInd];
    size_fronts = new int[NP*nInd];

    domin_to = new int[max(nInd*nInd, max_archive*max_archive)];
    times_dominated = new int[max(nInd, max_archive)];
    size_domin_to = new int[max(nInd, max_archive)];

    hypermat_assig = new int[NP*NP*nInd];

    namda = new double[nobj*nWeight];
    costs = new double[nInd*nInd*3];
    memo_dist = new double[NP*NP];
   
    archive_x = new double[max_archive*nvar];
    cpy_archive_x = new double[max_archive*nvar];
    archive_y = new double[max_archive*nobj];
    cpy_archive_y = new double[max_archive*nobj];
    f_archive =  new int[max_archive*max_archive];
    sf_archive = new int[max_archive];
    contribution_R2 = new double[max_archive];
    w_set = new set<pair<double, int> > [nWeight];
    rejected = new bool[max_archive];
    

    memset(hypermat_assig, -1, sizeof(int)*NP*NP*nInd);
    fill_n(fitness, NP*nInd, -1.0);
    fill_n(memo_dist, NP*NP, -1);
    fill_n(idealpoint, nobj, 1000000000);
    fill_n(archive_y, max_archive*nobj, 1000);
    fill_n(archive_x, max_archive*nvar, 0);
    memset(changed, false, sizeof(bool)*nInd*NP);
    // Load weight vectors
    for(int i=0; i< nWeight; i++)
	for(int j=0; j<nobj; j++)
	 readf>>namda[i*nobj + j];
    
    for(int i=0; i < NP; i++)
    {
	pool[i].x_var = x_var + i*nvar*nInd;
	pool[i].y_obj = y_obj + i*nobj*nInd;
        pool[i].changed = changed + i*nInd;
        pool[i].f = fronts + i*nInd*nInd;
        pool[i].ff = fitness + i*nInd;
        pool[i].sf = size_fronts + i*nInd;

        for(int k= 0; k < nInd; k++)
	{  
	   for(int n = 0; n<nvar; n++) pool[i].x_var[k*nvar + n] = vlowBound[n] + rnd_uni*(vuppBound[n] - vlowBound[n]);     
	   obj_eval(pool[i].x_var + k*nvar, pool[i].y_obj+ k*nobj);
	   update_reference(pool[i].y_obj+ k*nobj);
	   add_archive(pool[i].y_obj+ k*nobj, pool[i].x_var + k*nvar);
	}
        fill_n(pool[i].ff, nInd, -1.0);
        non_dominated_sorting(pool[i].y_obj, pool[i].f, pool[i].sf, nInd);
 	
	if( i < nPop)
           parent_idx[i] = i;
	else
	   child_idx[i-nPop] = i;
	nfes +=nInd;
     } 
     readf.close();
}
void CMOEAD::evol_population()
{
   for(int i = 0; i < nOffspring; i++)
   {
      int idx1=parent_idx[rand()% nPop], idx2=parent_idx[rand()%nPop], idx3=parent_idx[rand()%nPop], idx_target = parent_idx[i], c_idx = child_idx[i];
      while(idx1 == idx_target) idx1=parent_idx[rand()%nPop];
      while(idx2 == idx1 || idx2 == idx_target) idx2=parent_idx[rand()%nPop];
      while(idx3 == idx2 || idx3 == idx1 || idx3 == idx_target) idx3=parent_idx[rand()%nPop];
      for(int m = 0; m < nobj*nInd; m++)   pool[c_idx].y_obj[m] = pool[idx_target].y_obj[m];
      for(int n = 0; n < nvar*nInd; n++)   pool[c_idx].x_var[n] = pool[idx_target].x_var[n];

      int *asg_1  = pointer_hyp(idx_target, idx1) , *asg_2 = pointer_hyp(idx_target, idx2), *asg_3 = pointer_hyp(idx_target, idx3);
//      if( *asg_1 == -1 || *asg_2 == -1 || *asg_3 == -1)
//      {
//          for(int ii = 0; ii < nInd; ii++)
//          {
//         	for(int jj = 0; jj < nInd; jj++)
//         	{
//        	   if(*asg_1 == -1) costs[ii*nInd+jj] =-distance_obj(pool[idx_target].y_obj+ii*nobj, pool[idx1].y_obj+nobj*jj);
//        	   if(*asg_2 == -1) costs[nInd*nInd+ii*nInd+jj] =-distance_obj(pool[idx_target].y_obj+ii*nobj, pool[idx2].y_obj+nobj*jj);
//        	   if(*asg_3 == -1) costs[2*nInd*nInd + ii*nInd+jj] =-distance_obj(pool[idx_target].y_obj + ii*nobj, pool[idx3].y_obj+nobj*jj);
//         	}
//          }
//      }
//      if(*asg_1 == -1) KM.hungarian(costs, asg_1);
//      if(*asg_2 == -1) KM.hungarian(costs + nInd*nInd, asg_2);
//      if(*asg_3 == -1) KM.hungarian(costs + 2*nInd*nInd, asg_3);
      diff_evo_xoverA_exp(pool[idx_target].x_var, pool[idx1].x_var, pool[idx2].x_var, pool[idx3].x_var, pool[c_idx].x_var, CR, F, asg_1, asg_2, asg_3, pool[c_idx].changed);

      for(int k = 0; k < nInd; k++)
      {	
	   if(!pool[c_idx].changed[k]) continue;
           realmutation(pool[c_idx].x_var+k*nvar, 1.0/nvar);
           obj_eval(pool[c_idx].x_var+k*nvar, pool[c_idx].y_obj+k*nobj);
           update_reference(pool[c_idx].y_obj+k*nobj); 
	   pool[c_idx].changed[k]=false;
	   add_archive(pool[c_idx].y_obj+ k*nobj, pool[c_idx].x_var + k*nvar);
     	   nfes++;
      }
      non_dominated_sorting(pool[c_idx].y_obj, pool[c_idx].f, pool[c_idx].sf, nInd);
      
     for(int s = 0; s< nPop; s++)
     {
	int p_idx = parent_idx[s];
      *(pointer_hyp(p_idx, c_idx))= *(pointer_hyp(c_idx, p_idx)) =-1, *(pointer_dist(c_idx, p_idx))=-1;
     }
   }
   replacement_phase();
}
void CMOEAD::exec_emo(int run)
{
        char filename1[5024];
        char filename2[5024];
	srand(run);
	//initialization
	nfes      = 0;
	init_population();
        replacement_phase();

	sprintf(filename1,"%s/POS/POS_R2_EMOA_%s_RUN%d_seed_%d_nobj_%d_nvar_%d_DI_%lf_DF_%lf_CR_%lf_F_%lf",strpath, strTestInstance,run, run, nobj, nvar, Di/sqrt(nvar*nInd), Df, CR, F);
	sprintf(filename2,"%s/POF/POF_R2_EMOA_%s_RUN%d_seed_%d_nobj_%d_nvar_%d_DI_%lf_DF_%lf_CR_%lf_F_%lf",strpath, strTestInstance,run, run, nobj, nvar, Di/sqrt(nvar*nInd), Df, CR, F);
        long long current = nfes;
	long long accumulator = 0, bef = nfes;
	//save_pos(filename1);
        save_front(filename2);
	while(nfes<max_nfes)
	{
		update_parameterD();
		accumulator += nfes - bef ;
                if(accumulator > 0.1*(max_nfes)  )
		{
	           accumulator -= 0.1*(max_nfes);
		   save_pos(filename1);
		   save_front(filename2);
		}
		bef=nfes;
		evol_population();
	}

        update_archive();
	save_pos(filename1);
	save_front(filename2);
}
void CMOEAD::save_front(char saveFilename[4024])
{

    std::fstream fout;
    fout.open(saveFilename,fstream::app|fstream::out );
    for(int n=0; n < nPop; n++)
    {
       for(int i = 0; i < nInd; i++)
       {
          for(int k=0;k<nobj;k++)
             fout<<pool[parent_idx[n]].y_obj[i*nobj + k]<<"  ";
          for(int k=0;k<nobj;k++)
             fout<<pool[child_idx[n]].y_obj[i*nobj +k]<<"  ";
          fout<<"\n";
      }
    }
    for(int n=0; n < size_archive; n++)
    {
          for(int k=0;k<nobj;k++)
             fout<<archive_y[n*nobj + k]<<"  ";
          fout<<"\n";
    }
    fout.close();
}
void CMOEAD::save_pos(char saveFilename[4024])
{
   std::fstream fout; //fout.open(saveFilename,std::ios::out);
   fout.open(saveFilename, fstream::app|fstream::out);
   fout.close();
}
void CMOEAD::replacement_phase()
{
  fill_n(fitness, NP*nInd, -1.0);
  auto compare_l = [&](const int &a, const int &b)->bool{return compareR2(a, b);};
  unordered_set<int> penalized, survivors;
  priority_queue<int, vector<int>, decltype(compare_l)> candidates(compare_l);
  for(int i = 0; i < NP; i++) candidates.push(i);
  while(!candidates.empty() && survivors.size() < nPop)
  {
     int idx = candidates.top(); candidates.pop();
     bool flagIsSurvivor=true;
     for(auto s:survivors) 
	if(distance_var(s, idx) < D){flagIsSurvivor = false; break;}
     if(flagIsSurvivor)
	survivors.insert(idx);
     else penalized.insert(idx);
  }
  vector<double> mindist(NP, DBL_MAX);
  for(auto p:penalized)
    for(auto s:survivors)
       mindist[p] = min(mindist[p], distance_var(s, p));
  while(survivors.size() < nPop)
  {
     pair<double, int> max_dcn(-1, -1);
     for(auto p:penalized) 
	if(max_dcn.first < mindist[p]) max_dcn = make_pair(mindist[p], p);
     //update max_dcn..
     survivors.insert(max_dcn.second);
     penalized.erase(max_dcn.second);
     for(auto p:penalized) mindist[p] = min(mindist[p], distance_var(max_dcn.second, p));
  } 
  vector<bool> checked(NP, false);
  int i = 0;
  for(auto s:survivors)
  {
     parent_idx[i++] = s;
     checked[s] = true;
  }
  for(int i = 0, j=0; i < NP; i++) if(!checked[i]) child_idx[j++] = i;
}
void CMOEAD::obj_eval(double *x_var, double *y_obj)
{
   if(!strcmp("UF1", strTestInstance))  CEC09_F1(y_obj, x_var);
  else if(!strcmp("UF2", strTestInstance))  CEC09_F2(y_obj, x_var);
  else if(!strcmp("UF3", strTestInstance))  CEC09_F3(y_obj, x_var);
  else if(!strcmp("UF4", strTestInstance))  CEC09_F4(y_obj, x_var);
  else if(!strcmp("UF5", strTestInstance))  CEC09_F5(y_obj, x_var);
  else if(!strcmp("UF6", strTestInstance))  CEC09_F6(y_obj, x_var);
  else if(!strcmp("UF7", strTestInstance))  CEC09_F7(y_obj, x_var);
  else if(!strcmp("UF8", strTestInstance))  CEC09_F8(y_obj, x_var);
  else if(!strcmp("UF9", strTestInstance))  CEC09_F9(y_obj, x_var);
  else if(!strcmp("UF10", strTestInstance)) CEC09_F10(y_obj, x_var);
  else if(!strcmp("WFG1", strTestInstance))  wfg1(y_obj, x_var);
  else if(!strcmp("WFG2", strTestInstance))  wfg2(y_obj, x_var);
  else if(!strcmp("WFG3", strTestInstance))  wfg3(y_obj, x_var);
  else if(!strcmp("WFG4", strTestInstance))  wfg4(y_obj, x_var);
  else if(!strcmp("WFG5", strTestInstance))  wfg5(y_obj, x_var);
  else if(!strcmp("WFG6", strTestInstance))  wfg6(y_obj, x_var);
  else if(!strcmp("WFG7", strTestInstance))  wfg7(y_obj, x_var);
  else if(!strcmp("WFG8", strTestInstance))  wfg8(y_obj, x_var);
  else if(!strcmp("WFG9", strTestInstance))  wfg9(y_obj, x_var);
  else if(!strcmp("DTLZ1", strTestInstance))  dtlz1(y_obj, x_var);
  else if(!strcmp("DTLZ2", strTestInstance))  dtlz2(y_obj, x_var);
  else if(!strcmp("DTLZ3", strTestInstance))  dtlz3(y_obj, x_var);
  else if(!strcmp("DTLZ4", strTestInstance))  dtlz4(y_obj, x_var);
  else if(!strcmp("DTLZ5", strTestInstance))  dtlz5(y_obj, x_var);
  else if(!strcmp("DTLZ6", strTestInstance))  dtlz6(y_obj, x_var);
  else if(!strcmp("DTLZ7", strTestInstance))  dtlz7(y_obj, x_var);
}


void CMOEAD::realmutation(double *x_var, double rate)
{
    long double rnd, delta1, delta2, mut_pow, deltaq;
    long double y, yl, yu, val, xy;
    long double eta_m = etam;
    int id_rnd = rand()%nvar;
    for (int j=0; j<nvar; j++)
    {
        //if (rnd_uni<=rate || id_rnd==j)
        if (rnd_uni<=rate)
        {
            y  = x_var[j];
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
            x_var[j] = y;
        }
    }
}
void CMOEAD::diff_evo_xoverA_exp(double *ind0, double *ind1, double *ind2, double *ind3, double *child, double CR, double F, int *asg_1, int *asg_2, int *asg_3, bool *childchanged)
{
	int n = rand()%(nvar*nInd);
	int cont =0;
	do{
         int i = n/nvar, j = n%nvar;
         child[n] = ind1[asg_1[i]*nvar + j] + F*(ind2[asg_2[i]*nvar + j] - ind3[asg_3[i]*nvar+j]);
	  if(child[n] < vlowBound[j]) child[n] = ind0[n];
	  if(child[n]>vuppBound[j]) child[n] = ind0[n];
	  if(child[n] < vlowBound[j]) child[n] = vlowBound[n];
	  if(child[n]>vuppBound[j]) child[n] = vuppBound[n];

	   childchanged[i]=true;
	   n++;
	   n %= (nvar*nInd);
	   cont++;
	}
        while(rnd_uni < CR && cont < (nvar*nInd) );
}
double CMOEAD::distance_obj(double *a, double *b)
{
   double dist =0.0;
   for(int i = 0; i < nobj; i++)
     dist += ((a[i]-b[i])*(a[i]-b[i]));
   return dist;
}
double CMOEAD::fitnessfunction(double *y_obj, double *namda)
{
	double max_fun = -1.0e+30;
	for(int n=0; n<nobj; n++)
	{
		double feval, diff = fabs(y_obj[n] - idealpoint[n]);
		if(namda[n]==0) feval = diff/0.0001;
		else feval = diff/namda[n];
		max_fun = max(max_fun, feval);
	}
	return max_fun;
}
void CMOEAD::non_dominated_sorting(double *y_obj, int *f, int *sf, int npoints) 
{
  memset(times_dominated, 0, sizeof(int)*npoints);
  memset(size_domin_to, 0, sizeof(int)*npoints);
  memset(sf,0, sizeof(int)*npoints);
  for(int pidx1=0; pidx1 < npoints; pidx1++)
  {
    for(int pidx2=0; pidx2 < npoints; pidx2++)
    {
      if(pidx1 == pidx2) continue;
      if( dominate(y_obj + pidx1*nobj, y_obj + pidx2*nobj)) domin_to[pidx1*npoints + size_domin_to[pidx1]] = pidx2, size_domin_to[pidx1]++;
      else if( dominate(y_obj + pidx2*nobj, y_obj + pidx1*nobj) ) times_dominated[pidx1]++;
    }
    if(times_dominated[pidx1] == 0) f[(*sf)++] = pidx1;
 }
  while(*sf !=0)
  {
     for(int i = 0; i < *sf; i++)
     {
        for(int j = 0, idxp=*(f+i); j < size_domin_to[idxp]; j++)
	{
           int idxc = domin_to[idxp*npoints+j];
	   times_dominated[idxc]--;
	   if(times_dominated[idxc] == 0){ f[npoints+*(sf+1)] = idxc; (*(sf+1))++;}
	}
     }
     f+=npoints;
     sf++;
  }
}
bool CMOEAD::dominate(double *y_obj1, double *y_obj2)
{
   int equal = 0;
   for(int i = 0; i < nobj; i++)
	if(y_obj2[i]<y_obj1[i]) return false;
	else if(y_obj2[i] == y_obj1[i]) equal++;
   if(equal == nobj) return false;
   return true;
}
bool CMOEAD::compareR2(int a, int b)
{
   for(int i = 0; i < nInd; i++)
   {
        eval_R2(pool[a].y_obj, pool[a].f+i*nInd, pool[a].sf[i], pool[a].ff[i]);
        eval_R2(pool[b].y_obj, pool[b].f+i*nInd, pool[b].sf[i], pool[b].ff[i]);
      if(pool[a].ff[i] > pool[b].ff[i]) return true;
      else if(pool[b].ff[i] > pool[a].ff[i]) return false;
   }
   return false;
}
void CMOEAD::eval_R2(double *y_obj, int *front, int size_front, double &fitness)
{
     if(size_front == 0) {fitness = 0.0; return;} 
     if( fitness > 0 ) return;
     fitness = 0.0;
     for(int w = 0; w < nWeight; w++)
     {
       double minv = DBL_MAX;
       for(int i = 0; i < size_front; i++)
           minv = min(minv, fitnessfunction(&y_obj[front[i]*nobj], &namda[w*nobj]));
       fitness += minv;
     }
}
void CMOEAD::add_archive(double *y_obj, double *x_var)
{
   for(int m = 0; m < nobj; m++) archive_y[n_archive*nobj + m ] = y_obj[m];
   for(int n = 0; n < nvar; n++) archive_x[n_archive*nvar + n ] = x_var[n];
   n_archive++;
   if( n_archive == max_archive) update_archive();
}
void CMOEAD::update_archive()
{
  memset(contribution_R2, 0, sizeof(double)*max_archive);
  memset(rejected, true, sizeof(bool)*max_archive);
  non_dominated_sorting(archive_y, f_archive, sf_archive, n_archive);  
  if( *sf_archive < size_archive)
  {
     n_archive = 0;
     for(int i = 0; i < *sf_archive; i++)
     {
       int idx = f_archive[i];
       for(int m = 0; m < nobj; m++) cpy_archive_y[n_archive*nobj + m] = archive_y[idx*nobj+m];
       for(int n = 0; n < nvar; n++) cpy_archive_x[n_archive*nvar + n] = archive_x[idx*nvar+n];
       n_archive++;
     }
     for(int i = 0; i < n_archive; i++)
     {
       for(int m = 0; m < nobj; m++) archive_y[i*nobj + m] = cpy_archive_y[i*nobj+m];
       for(int n = 0; n < nvar; n++) archive_x[i*nvar + n] = cpy_archive_x[i*nvar+n];
     }
    return;
  }
  for(int  i = 0; i < nWeight; i++) w_set[i].clear();

  for(int w_idx = 0; w_idx < nWeight; w_idx++)
  {
      for(int i = 0; i < *sf_archive; i++)
      {
	 int idx = f_archive[i];
         w_set[w_idx].insert(make_pair(fitnessfunction(archive_y + idx*nobj, &namda[w_idx*nobj]), idx));
	 rejected[idx]=false;
      }
      contribution_R2[w_set[w_idx].begin()->second] += (next(w_set[w_idx].begin(), 1)->first - next(w_set[w_idx].begin(), 0)->first);
  }
  n_archive = *sf_archive;
  while(n_archive > size_archive)
  {
      pair<double, int> min_info(10000000, -1);
      //take the worst contribution-individual..                   
      for(int i = 0; i < *sf_archive; i++)
      {
	 int idx = f_archive[i];
         if(rejected[idx])continue;
         if(min_info.first > contribution_R2[idx])
           min_info = make_pair(contribution_R2[idx], idx);
      }
     //update contributions... 
     memset(contribution_R2, 0, sizeof(double)*max_archive);
     for(int w_idx = 0; w_idx < nWeight; w_idx++)
     {
        w_set[w_idx].erase(make_pair(fitnessfunction(archive_y + min_info.second*nobj, namda + w_idx*nobj), min_info.second));
        contribution_R2[w_set[w_idx].begin()->second] += (next(w_set[w_idx].begin(), 1)->first - next(w_set[w_idx].begin(), 0)->first);
     }
     rejected[min_info.second]=true;
     n_archive--;
  }
  n_archive = 0;
  for(int i = 0; i < *sf_archive; i++)
  {
     int idx = f_archive[i];
     if(rejected[idx]) continue;
    for(int m = 0; m < nobj; m++) cpy_archive_y[n_archive*nobj + m] = archive_y[idx*nobj+m];
    for(int n = 0; n < nvar; n++) cpy_archive_x[n_archive*nvar + n] = archive_x[idx*nvar+n];
    n_archive++;
  }
  for(int i = 0; i < n_archive; i++)
  {
    for(int m = 0; m < nobj; m++) archive_y[i*nobj + m] = cpy_archive_y[i*nobj+m];
    for(int n = 0; n < nvar; n++) archive_x[i*nvar + n] = cpy_archive_x[i*nvar+n];
  }
}
#endif
