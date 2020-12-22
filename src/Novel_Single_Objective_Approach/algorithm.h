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
#include "recomb.h"
#include "problem.h"

class CMOEAD
{

   public:
	CMOEAD();
	virtual ~CMOEAD();

	void init_population();                
	void obj_eval(vector<double> &x_var, vector<double> &y_obj);
	bool update_reference(vector<double> &point); 
	void replacement_phase();
	void evol_population();                                    
	void exec_emo(int run);
	void save_front(char savefilename[4024]); 
	void save_pos(char savefilename[4024]);
	void update_parameterD();
	void update_external_file(vector<vector<double> > &archive);
	double distance_var( int a, int b);
	//inline int* pointer_hyp(int a, int b){ if(a > b) swap(a, b);  return hypermat_assig +a*(nPop+nOffspring)*nInd + b*nInd; }
	inline int* pointer_hyp(int a, int b){  return hypermat_assig +a*(nPop+nOffspring)*nInd + b*nInd; }
	inline double* pointer_dist(int a, int b){ if(a > b) swap(a, b);  return memo_dist + a*(nPop+nOffspring) + b; }
  	inline void get_cost(double *cost, vector<vector<double> > &set1, vector<vector<double> > &set2)
        {
	   for(int i = 0;i < nInd*nInd;i++) cost[i] = -distance_obj(set1[i/nInd], set2[i%nInd]);
        }
	vector <strIndividual> pool;
   private:
        
	vector<int> child_idx, parent_idx, inv_parent_idx; 
	vector<vector<double> > R2_pop;     // weight vector
	// algorithm parameters
	long long nfes;          //  the number of function evluations
	double	D;	//Current minimum distance
	int n_archive;
};
CMOEAD::CMOEAD()
{

}
CMOEAD::~CMOEAD()
{
   delete[] namda;
   delete[] cost_1;
   delete[] cost_2;
   delete[] cost_3;
//   delete[] asg_1; 
//   delete[] asg_2;
//   delete[] asg_3;
   delete[] hypermat_assig;
   delete[] memo_dist;
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
   int *asg_1 = pointer_hyp(a,b);
   if(*asg_1==-1)
   {
     for(int i = 0; i < nInd; i++)
      for(int j = 0; j < nInd; j++)
       cost_1[i*nInd+j] = -distance_obj(pool[a].y_obj[i], pool[b].y_obj[j]);
     KM.hungarian(cost_1, asg_1);
   }
   double dist = 0.0;
   for(int i = 0; i < nInd; i++)
   {
      for(int j = 0; j < nvar; j++)
      {
         double factor = (pool[a].x_var[i][j]-pool[b].x_var[asg_1[i]][j])/(vuppBound[j]-vlowBound[j]);
         dist += factor*factor;
      }
   }
   (*distab) = sqrt(dist);
   return sqrt(dist);
}
void CMOEAD::init_population()
{
    idealpoint = vector<double>(nobj, 1.0e+30);
    char filename[1024];
    // Read weight vectors from a data file
    sprintf(filename,"%s/ParameterSetting/Weight/W%dD_%d.dat", strpath, nobj, nWeight);
    std::ifstream readf(filename);
    KM.init(nInd);
    namda = new double[nobj*nWeight];
    cost_1 = new double[nInd*nInd];
    cost_2 = new double[nInd*nInd];
    cost_3 = new double[nInd*nInd];
    hypermat_assig = new int[(nPop+nOffspring)*(nPop+nOffspring)*nInd];
    memset(hypermat_assig, -1, sizeof(int)*(nPop+nOffspring)*(nPop+nOffspring)*nInd);
    memo_dist = new double[(nPop+nOffspring)*(nPop+nOffspring)];
    for(int i = 0; i < (nPop+nOffspring)*(nPop+nOffspring); i++) memo_dist[i]=-1;

    n_archive=100;
    // Load weight vectors
    for(int i=0; i< nWeight; i++)
	for(int j=0; j<nobj; j++)
	 readf>>namda[i*nobj + j];
    for(int i=0; i< nPop+nOffspring; i++)
    {
        strIndividual ind;
        ind.x_var.assign(nInd, vector<double>(nvar, 0));
        ind.y_obj.assign(nInd, vector<double>(nobj, 0));
	ind.changed.assign(nInd, false);
	// random initialization..
        for(int k= 0; k < nInd; k++)
	{  
	   for(int n = 0; n<nvar; n++)
             ind.x_var[k][n] = vlowBound[n] + rnd_uni*(vuppBound[n] - vlowBound[n]);     
	    obj_eval(ind.x_var[k], ind.y_obj[k]), update_reference(ind.y_obj[k]), R2_pop.push_back(ind.y_obj[k]);
	}

        ind.fronts = non_dominated_sorting(ind.y_obj);
        for(int k = 0; k < ind.fronts.size(); k++) eval_R2(ind, k);
	pool.push_back(ind); 
	if( i < nPop)
	   parent_idx.push_back(i);
	else
	   child_idx.push_back(i);
	nfes +=nInd;
     }
     update_external_file(R2_pop);
     readf.close();
}
bool CMOEAD::update_reference(vector<double> &point)
{
  bool updated = false;
  for(int n=0; n<nobj; n++)
     if(point[n]<idealpoint[n])
        idealpoint[n] = point[n], updated=true;
  return updated;
}
void CMOEAD::evol_population()
{
   for(int i = 0; i < nOffspring; i++)
   {
      int idx1=parent_idx[rand()% nPop], idx2=parent_idx[rand()%nPop], idx3=parent_idx[rand()%nPop], idx_target = parent_idx[i];
      while(idx1 == idx_target) idx1=parent_idx[rand()%nPop];
      while(idx2 == idx1 || idx2 == idx_target) idx2=parent_idx[rand()%nPop];
      while(idx3 == idx2 || idx3 == idx1 || idx3 == idx_target) idx3=parent_idx[rand()%nPop];

      strIndividual &child = pool[child_idx[i]], &ind0 = pool[idx_target], &ind1 = pool[idx1], &ind2 = pool[idx2], &ind3 = pool[idx3];
      child = ind0;
      int *asg_1  = pointer_hyp(idx_target, idx1) , *asg_2 =pointer_hyp(idx_target, idx2), *asg_3 = pointer_hyp(idx_target, idx3);
      if(*asg_1 == -1) get_cost(cost_1, ind0.y_obj, ind1.y_obj), KM.hungarian(cost_1, asg_1);
      if(*asg_2 == -1) get_cost(cost_2, ind0.y_obj, ind1.y_obj), KM.hungarian(cost_2, asg_2);
      if(*asg_3 == -1) get_cost(cost_3, ind0.y_obj, ind1.y_obj), KM.hungarian(cost_3, asg_3);

      diff_evo_xoverA_exp(ind0, ind1, ind2, ind3, child, CR, F, asg_1, asg_2, asg_3);

      for(int k = 0; k < nInd; k++)
      {	
	   if(!child.changed[k])continue;
           realmutation(child.x_var[k], 1.0/nvar);
           obj_eval(child.x_var[k], child.y_obj[k]);
           update_reference(child.y_obj[k]); 
     	   R2_pop.push_back(child.y_obj[k]);
     	   nfes++;
      }
      child.changed.assign(nInd, false);
      child.fronts = non_dominated_sorting(child.y_obj);

      if(R2_pop.size() >= 2*n_archive) 
      update_external_file(R2_pop);
   }
   //updating all R2 contributions cuz reference vector has changed..
   for(int idx1 = 0; idx1 < pool.size(); idx1++)
   {
      strIndividual &ind = pool[idx1];
      ind.fitness.clear();
//      for(int k = 0; k < ind.fronts.size(); k++) eval_R2(ind, k);
      for(auto idx2:child_idx) *(pointer_hyp(idx1, idx2))=-1, *(pointer_dist(idx1, idx2))=-1;
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
	//	   save_pos(filename1);
		   save_front(filename2);
	//	cout << nfes/(double)max_nfes << " ... percent.."<<endl; 
		}
		bef=nfes;
		evol_population();
	}

        update_external_file(R2_pop);
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
             fout<<pool[parent_idx[n]].y_obj[i][k]<<"  ";
          for(int k=0;k<nobj;k++)
             fout<<pool[child_idx[n]].y_obj[i][k]<<"  ";

          fout<<"\n";
      }
    }
    for(int n=0; n < n_archive; n++)
    {
          for(int k=0;k<nobj;k++)
             fout<<R2_pop[n][k]<<"  ";
          fout<<"\n";
    }
    fout.close();
}
void CMOEAD::save_pos(char saveFilename[4024])
{
//   std::fstream fout; //fout.open(saveFilename,std::ios::out);
//   fout.open(saveFilename, fstream::app|fstream::out);
//   for(int n=0; n<nPop; n++)
//   {
//      for(int k=0;k<nvar;k++)
//         fout<<pool[parent_idx[n]].x_var[k] << "  ";
//      for(int k=0;k<nvar;k++)
//   	 fout<<pool[parent_idx[n]].x_var[k]/(vuppBound[k]-vlowBound[k]) << "  "; //fout<<population[n].indiv.x_var[k]<< fixed << setprecision(30) << "  ";
//   	fout<<"\n";
//   }
//   fout.close();
}
void CMOEAD::replacement_phase()
{
  auto compare_l = [&](const int &a, const int &b)->bool
  {
     strIndividual &ind_a = pool[a], &ind_b = pool[b];
//	return (ind_a.fitness>ind_b.fitness);
     int rank = 0, sf1 = ind_a.fronts.size(), sf2=ind_b.fronts.size();
     do{
	if(ind_a.fitness.size() <= rank) {eval_R2(ind_a, rank); continue;}
	if(ind_b.fitness.size() <= rank) {eval_R2(ind_b, rank); continue;}	
        if(ind_a.fitness[rank] > ind_b.fitness[rank]) return true;
        else if(ind_b.fitness[rank] > ind_a.fitness[rank]) return false;
        rank++;	
     }while(rank<nInd && rank < sf1 && rank <sf2);
     return true;
  };
  unordered_set<int> penalized, survivors;
  priority_queue<int, vector<int>, decltype(compare_l)> candidates(compare_l);
  for(int i = 0; i < pool.size(); i++) candidates.push(i);
  while(!candidates.empty() && survivors.size() < nPop)
  {
     int idx = candidates.top(); candidates.pop();
     bool flagIsSurvivor=true;
     for(auto s:survivors) 
	if(distance_var(s, idx) <= D){flagIsSurvivor = false; break;}
     if(flagIsSurvivor)
	survivors.insert(idx);
     else penalized.insert(idx);
  }
  vector<double> mindist(pool.size(), DBL_MAX);
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
  vector<bool> checked(pool.size(), false);
  int i = 0;
  for(auto s:survivors)
  {
     parent_idx[i++] = s;
     checked[s] = true;
  }
  for(int i = 0, j=0; i < pool.size(); i++) if(!checked[i]) child_idx[j++] = i;
}
void CMOEAD::obj_eval(vector<double> &x_var, vector<double> &y_obj)
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
//void CMOEAD::update_external_file(vector<vector<double> > &archive)
//{
//  vector<int> multiset_R2((int)archive.size());
//  for(int i = 0 ; i < multiset_R2.size(); i++) multiset_R2[i]=i;
//
//  vector<double> contribution_R2(multiset_R2.size(), 0);
//  vector< vector<double> > fitness_table(nWeight, vector<double>(multiset_R2.size()));
//  vector< set<pair<double, int> > > w_set(nWeight);
//  for(int w_idx = 0; w_idx < nWeight; w_idx++)
//  {
//      for(auto idx:multiset_R2)
//      {
//          double gx = fitnessfunction(archive[idx], &namda[w_idx*nobj]);//atof(sz);
//         fitness_table[w_idx][idx] = gx;
//         w_set[w_idx].insert(make_pair(gx, idx));
//      }
//      contribution_R2[w_set[w_idx].begin()->second] += (next(w_set[w_idx].begin(), 1)->first - next(w_set[w_idx].begin(), 0)->first);
//  }
//  while(multiset_R2.size() > n_archive)
//  {
//      pair<double, int> min_info(10000000, -1);
//      //take the worst contribution-individual..                   
//      for(int idx = 0; idx < multiset_R2.size(); idx++)
//      {
//         if(min_info.first > contribution_R2[multiset_R2[idx]])
//           min_info = make_pair(contribution_R2[multiset_R2[idx]], idx);
//      }
//     //update contributions... 
//     contribution_R2.assign(archive.size(), 0.0);
//     for(int w_idx = 0; w_idx < nWeight; w_idx++)
//     {
//        w_set[w_idx].erase(make_pair(fitness_table[w_idx][multiset_R2[min_info.second]], multiset_R2[min_info.second]));
//
//        contribution_R2[w_set[w_idx].begin()->second] += (next(w_set[w_idx].begin(), 1)->first - next(w_set[w_idx].begin(), 0)->first);
//     }
//     iter_swap(multiset_R2.begin()+min_info.second, multiset_R2.end()-1);
//     multiset_R2.pop_back();
//  }
//  vector<vector<double > > tmp = archive;
//  archive.resize(n_archive);
//  for(int i = 0; i < n_archive; i++) archive[i]=tmp[multiset_R2[i]];
//}
void CMOEAD::update_external_file(vector<vector<double> > &archive)
{
  unordered_set<int> selected = non_dominated_sorting(archive)[0];

  vector<double> contribution_R2(archive.size(), 0.0);
  vector< set<pair<double, int> > > w_set(nWeight);
  for(int w_idx = 0; w_idx < nWeight; w_idx++)
  {
      for(auto idx:selected)
         w_set[w_idx].insert(make_pair(fitnessfunction(archive[idx], &namda[w_idx*nobj]), idx));
      contribution_R2[w_set[w_idx].begin()->second] += (next(w_set[w_idx].begin(), 1)->first - next(w_set[w_idx].begin(), 0)->first);
  }
  while(selected.size() > n_archive)
  {
      pair<double, int> min_info(10000000, -1);
      //take the worst contribution-individual..                   
      for(auto idx:selected)
      {
         if(min_info.first > contribution_R2[idx])
           min_info = make_pair(contribution_R2[idx], idx);
      }
     //update contributions... 
     contribution_R2.assign(archive.size(), 0.0);
     for(int w_idx = 0; w_idx < nWeight; w_idx++)
     {
        w_set[w_idx].erase(make_pair(fitnessfunction(archive[min_info.second], &namda[w_idx*nobj]), min_info.second));
        contribution_R2[w_set[w_idx].begin()->second] += (next(w_set[w_idx].begin(), 1)->first - next(w_set[w_idx].begin(), 0)->first);
     }
     selected.erase(min_info.second);
  }
  vector<vector<double> > tmp = archive;
  archive.clear();
  for(auto idx:selected) archive.push_back(tmp[idx]);
}
#endif
