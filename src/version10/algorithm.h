#ifndef __EVOLUTION_H_
#define __EVOLUTION_H_

#include <cfloat>
#include <set>
#include <queue>
#include <map>
#include <unordered_set>
#include <iomanip>
#include "global.h"
#include "recomb.h"
#include "common.h"
#include "individual.h"

class CMOEAD
{

   public:
	CMOEAD();
	virtual ~CMOEAD();

	void init_population();                
	void update_reference(CIndividual &ind); 
	void replacement_phase();
	void evol_population();                                    
	void exec_emo(int run);
	void save_front(char savefilename[4024]); 
	void save_pos(char savefilename[4024]);
	void update_parameterD();
	void update_external_file(vector<vector<double> > &archive);
	double distance_var( int a, int b);
	void dominance_information(vector<unordered_set<int> > &Sp, vector<int> &Np, unordered_set<int> &candidates_front, vector<vector<double> >  &y_obj);
	void eval_R2(CIndividual &indiv);
	vector<set<int> > non_dominated_sorting(vector<vector<double> > &y_obj);

   private:
        struct compare
        {
	   bool operator()(const pair<vector<double>, int> &a, const pair<vector<double>, int> &b){return b.first<<a.first; }
        };
	vector <CIndividual> pool;
	vector<int> child_idx, parent_idx, inv_parent_idx;
	vector<vector<double> > namda, R2_pop;     // weight vector


        vector<vector<double> > dist_matrix;
        vector<int> assignaments;
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

}
void CMOEAD::update_parameterD()
{
      double TElapsed = nfes, TEnd = max_nfes;
      D = Di - Di * (TElapsed / (TEnd*Df));
      D = max(D, 0.0);
}
double CMOEAD::distance_var(int a, int b)
{
///distance matrix..
    for(int i = 0; i < nInd; i++)
    {
      for(int j = 0; j < nInd; j++)
      {
	if(i==j) dist_matrix[i][j]=DBL_MAX;
	else dist_matrix[i][j] = distance_obj(pool[a].y_obj[i], pool[b].y_obj[j]);
      }
    }
    KuhnMunkres(assignaments, dist_matrix);
   double dist = 0 ;
   for(int i = 0; i < nInd; i++)
   {
      for(int j = 0; j < nvar; j++)
      {
         double factor = (pool[a].x_var[i][j]-pool[b].x_var[assignaments[i]][j])/(vuppBound[j]-vlowBound[j]);
         dist += factor*factor;
      }
   }
   return sqrt(dist);
}
void CMOEAD::init_population()
{
    idealpoint = vector<double>(nobj, 1.0e+30);
    char filename[1024];
    // Read weight vectors from a data file
    sprintf(filename,"%s/ParameterSetting/Weight/W%dD_%d.dat", strpath, nobj, nWeight);
    std::ifstream readf(filename);
    namda.resize(nWeight, vector<double> (nobj, 0.0));
    dist_matrix.assign(nInd, vector<double> (nInd, 0.0));
    assignaments.assign(nInd, 0);
    n_archive=100;
//    R2_pop.assign(n_archive, vector<double> (nobj, 1000000000));
    // Load weight vectors
    for(int i=0; i< nWeight; i++)
	for(int j=0; j<nobj; j++)
	 readf>>namda[i][j];

    for(int i=0; i< nPop+nOffspring; i++)
    {
	       CIndividual ind;
		// Randomize and evaluate solution
		ind.rnd_init();
		ind.obj_eval();
		// Initialize the reference point
		update_reference(ind);
		eval_R2(ind);
		// Save in the population
		pool.push_back(ind); 
		if( i < nPop)
		{
		   parent_idx.push_back(i);
		}
		else
		   child_idx.push_back(i);
		nfes +=nInd;
     }

     for(int i = 0; i < parent_idx.size(); i++)
     for(int j = 0; j < pool[parent_idx[i]].y_obj.size(); j++) R2_pop.push_back(pool[parent_idx[i]].y_obj[j]);

     for(int i = 0; i < child_idx.size(); i++)
     for(int j = 0; j < pool[child_idx[i]].y_obj.size(); j++) R2_pop.push_back(pool[child_idx[i]].y_obj[j]);

     update_external_file(R2_pop);
     readf.close();
}
void CMOEAD::update_reference(CIndividual &ind)
{
   for(int k = 0; k <nInd; k++)
   {
      for(int n=0; n<nobj; n++)
      {
         if(ind.y_obj[k][n]<idealpoint[n])
         {
            idealpoint[n] = ind.y_obj[k][n];
         }
      }
  }
}
void CMOEAD::evol_population()
{
   for(int i = 0; i < nOffspring; i++)
   {
      int idx_target = i;
      int idx1=rand()% nPop, idx2=rand()%nPop, idx3=rand()%nPop;
      while(idx1 == idx_target) idx1=rand()%nPop;
      while(idx2 == idx1 || idx2 == idx_target) idx2=rand()%nPop;
      while(idx3 == idx2 || idx3 == idx1 || idx3 == idx_target) idx3=rand()%nPop;
	vector<bool> changed(nInd, false);
      // produce a child solution
      CIndividual &child = pool[child_idx[i]];
   //   diff_evo_xoverA_knn(pool[parent_idx[idx_target]], pool[parent_idx[idx1]], pool[parent_idx[idx2]], pool[parent_idx[idx3]], child, CR, F);
      diff_evo_xoverA_exp(pool[parent_idx[idx_target]], pool[parent_idx[idx1]], pool[parent_idx[idx2]], pool[parent_idx[idx3]], child, CR, F, changed);
     // diff_evo_xoverA(pool[parent_idx[idx_target]], pool[parent_idx[idx1]], pool[parent_idx[idx2]], pool[parent_idx[idx3]], child, CR, F);
  //        realmutation(child.x_var[rand()%nInd], 1.0/(nvar));

      for(int k = 0; k < nInd; k++)
      {
         if(changed[k])
         {
	   nfes++;
           realmutation(child.x_var[k], 1.0/(nvar));
	 }
      }
      child.obj_eval(changed);
      for(int k = 0; k < nInd; k++)
      {
         if(changed[k])
	 R2_pop.push_back(child.y_obj[k]);
      }
      update_reference(child); //O(M)
   }
 
  if(R2_pop.size() >= 200) 
   update_external_file(R2_pop);
   replacement_phase();
}
void CMOEAD::exec_emo(int run)
{
        char filename1[5024];
        char filename2[5024];
	seed = run;
	seed = (seed + 23)%1377;
	rnd_uni_init = -(long)seed;

	//initialization
	nfes      = 0;
	init_population();

	sprintf(filename1,"%s/POS/POS_R2_EMOA_%s_RUN%d_seed_%d_nobj_%d_nvar_%d_DI_%lf_DF_%lf_CR_%lf_F_%lf",strpath, strTestInstance,run, seed, nobj, nvar, Di/sqrt(nvar*nInd), Df, CR, F);
	sprintf(filename2,"%s/POF/POF_R2_EMOA_%s_RUN%d_seed_%d_nobj_%d_nvar_%d_DI_%lf_DF_%lf_CR_%lf_F_%lf",strpath, strTestInstance,run, seed, nobj, nvar, Di/sqrt(nvar*nInd), Df, CR, F);
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
		}
		bef=nfes;
		evol_population();
//	        nfes += nOffspring*nInd;
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
//    for(int n=0; n < nPop; n++)
//    {
//       for(int i = 0; i < nInd; i++)
//       {
//          for(int k=0;k<nobj;k++)
//             fout<<pool[child_idx[n]].y_obj[i][k]<<"  ";
//          fout<<"\n";
//      }
//    }

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
  unordered_set<int> penalized, survivors;
  priority_queue<pair<vector<double>, int>, vector<pair<vector<double>, int>>, compare> candidates;
  for(int i = 0; i < pool.size(); i++)
  {
    eval_R2(pool[i]);
    candidates.push(make_pair(pool[i].fitness, i));
  }
  while(!candidates.empty() && survivors.size() < nPop)
  {
     int idx = candidates.top().second; candidates.pop();
     bool flagIsSurvivor=true;
     for(auto s:survivors) 
	if(distance_var(s, idx) < D){flagIsSurvivor = false; break;}
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
void CMOEAD::eval_R2(CIndividual &indiv)
{
  vector<set<int> > fronts = non_dominated_sorting(indiv.y_obj);
  indiv.fitness.assign(nInd, 0);
  //adding earch rank...
  for(int r = 0; r < fronts.size(); r++)
  {
     for(int w = 0; w < nWeight; w++)
     {
       double minv = DBL_MAX;
       for(auto k:fronts[r])
       {
           minv = min(minv, fitnessfunction(indiv.y_obj[k], namda[w]));
       } 
       indiv.fitness[r] += minv;
     }
  }
}
vector<set<int> > CMOEAD::non_dominated_sorting(vector<vector<double> > &y_obj)
{
  //dominance ranking..
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
void CMOEAD::update_external_file(vector<vector<double> > &archive)
{
  vector<int> multiset_R2((int)archive.size());
  for(int i = 0 ; i < multiset_R2.size(); i++) multiset_R2[i]=i;

  vector<double> contribution_R2(multiset_R2.size(), 0);
  vector< vector<double> > fitness_table(nWeight, vector<double>(multiset_R2.size()));
  vector< set<pair<double, int> > > w_set(nWeight);
  for(int w_idx = 0; w_idx < nWeight; w_idx++)
  {
      for(auto idx:multiset_R2)
      {
             double gx = fitnessfunction(archive[idx], namda[w_idx]);//atof(sz);
         fitness_table[w_idx][idx] = gx;
         w_set[w_idx].insert(make_pair(gx, idx));
      }
      contribution_R2[w_set[w_idx].begin()->second] += (next(w_set[w_idx].begin(), 1)->first - next(w_set[w_idx].begin(), 0)->first);
  }
  while(multiset_R2.size() > n_archive)
  {
      pair<double, int> min_info(10000000, -1);
      //take the worst contribution-individual..                   
      for(int idx = 0; idx < multiset_R2.size(); idx++)
      {
         if(min_info.first > contribution_R2[multiset_R2[idx]])
           min_info = make_pair(contribution_R2[multiset_R2[idx]], idx);
      }
     //update contributions... 
     contribution_R2.assign(archive.size(), 0.0);
     for(int w_idx = 0; w_idx < nWeight; w_idx++)
     {
        w_set[w_idx].erase(make_pair(fitness_table[w_idx][multiset_R2[min_info.second]], multiset_R2[min_info.second]));

        contribution_R2[w_set[w_idx].begin()->second] += (next(w_set[w_idx].begin(), 1)->first - next(w_set[w_idx].begin(), 0)->first);
     }
     iter_swap(multiset_R2.begin()+min_info.second, multiset_R2.end()-1);
     multiset_R2.pop_back();
  }
  vector<vector<double > > tmp = archive;
  archive.resize(n_archive);
  for(int i = 0; i < n_archive; i++) archive[i]=tmp[multiset_R2[i]];
}
#endif
