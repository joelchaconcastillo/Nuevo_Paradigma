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
	void replacement_phase(CIndividual &parent,CIndividual &child);
	void evol_population();                                    
	void exec_emo(int run);
	void save_front(char savefilename[4024]); 
	void save_pos(char savefilename[4024]);
	void update_parameterD();
	double distance_var( vector<double> &a, vector<double> &b);
	void penalization(unordered_set<int> &candidates, unordered_set<int> &penalized, vector<double> &distances, unordered_set<int> &candidates_front, vector<unordered_set<int> > &Sp, vector<int> &Np);
	void dominance_information(vector<unordered_set<int> > &Sp, vector<int> &Np, unordered_set<int> &candidates_front, vector<vector<double> >  &y_obj);
	void R2_contribution_subset(unordered_set<int> &candidates, unordered_set<int> &candidates_front, vector<int> &survivors, unordered_set<int> &survivors_front, unordered_set<int> &penalized, vector<set<pair<double, int> > > &survivors_weight, vector<double> &distances, vector<vector<double> > &table_fitness, vector<vector<double> > &x_var);
	void update_lowest_front(unordered_set<int> &candidates_front, unordered_set<int> &survivors_front, vector<unordered_set<int> > &Sp, vector<int> &Np);

	vector<set<int> > non_dominated_sorting(vector<vector<double> > &y_obj);

   private:
	vector <CIndividual> pool;
	vector<int> child_idx, parent_idx, inv_parent_idx;
	vector<vector<double> > namda;     // weight vector
	// algorithm parameters
	long long nfes;          //  the number of function evluations
	double	D;	//Current minimum distance

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
double CMOEAD::distance_var( vector<double> &a, vector<double> &b)
{
   double dist = 0 ;
   for(int i = 0; i < a.size(); i++)
   {
      double factor = (a[i]-b[i])/(vuppBound[i]-vlowBound[i]);
      dist += factor*factor;
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
		// Save in the population
		pool.push_back(ind); 
		if( i < nPop)
		{
		   parent_idx.push_back(i);
		}
		else
		   child_idx.push_back(i);
		nfes++;
     }
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
      // produce a child solution
      CIndividual &child = pool[child_idx[i]];
      for(int k = 0; k < nInd; k++)
      {
         diff_evo_xoverA(pool[parent_idx[idx_target]].x_var[k], pool[parent_idx[idx1]].x_var[k], pool[parent_idx[idx2]].x_var[k], pool[parent_idx[idx3]].x_var[k], child.x_var[k], CR, F);
      // apply polynomial mutation
          realmutation(child.x_var[k], 1.0/(nvar));
      }
      child.obj_eval();
      update_reference(child); //O(M)
   }
   for(int i = 0; i < nOffspring; i++)
   replacement_phase(pool[parent_idx[i]], pool[child_idx[i]]);
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

	sprintf(filename1,"%s/POS/POS_R2_EMOA_%s_RUN%d_seed_%d_nobj_%d_nvar_%d_DI_%lf_DF_%lf_CR_%lf_F_%lf",strpath, strTestInstance,run, seed, nobj, nvar, Di/sqrt(nvar), Df, CR, F);
	sprintf(filename2,"%s/POF/POF_R2_EMOA_%s_RUN%d_seed_%d_nobj_%d_nvar_%d_DI_%lf_DF_%lf_CR_%lf_F_%lf",strpath, strTestInstance,run, seed, nobj, nvar, Di/sqrt(nvar), Df, CR, F);
        long long current = nfes;
	long long accumulator = 0, bef = nfes;
	//save_pos(filename1);
        save_front(filename2);
	while(nfes<max_nfes)
	{
		update_parameterD();
		evol_population();
		accumulator += nfes - bef ;
                if(accumulator > 0.01*(max_nfes)  )
		{
	           accumulator -= 0.01*(max_nfes);
	//	   save_pos(filename1);
		   save_front(filename2);
		}
		bef=nfes;
	        nfes += nOffspring*nInd;
	}
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
          fout<<"\n";
      }
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
void CMOEAD::dominance_information(vector<unordered_set<int> > &Sp, vector<int> &Np, unordered_set<int> &candidates_front, vector<vector<double> >  &y_obj)
{
   for(int pidx1=0; pidx1< y_obj.size(); pidx1++)
   {
      for(int pidx2=0; pidx2< y_obj.size(); pidx2++)
      {
	if(pidx1 == pidx2) continue;
        if( y_obj[pidx1] < y_obj[pidx2]) Sp[pidx1].insert(pidx2);
 	else if( y_obj[pidx2] < y_obj[pidx1]) Np[pidx1]++;
      }
      if( Np[pidx1] == 0) candidates_front.insert(pidx1);
   }
}
void CMOEAD::penalization(unordered_set<int> &candidates, unordered_set<int> &penalized, vector<double> &distances, unordered_set<int> &candidates_front, vector<unordered_set<int> > &Sp, vector<int> &Np)
{
  unordered_set<int> tmp = candidates;
  for(auto c_idx:tmp)
  {
     if( distances[c_idx] <= D ) 
     {
	penalized.insert(c_idx);
        for(auto idx:Sp[c_idx])
	{
	  Np[idx]--; 
	  if(Np[idx]==0 && distances[idx] > D) candidates_front.insert(idx);
	}
        for(int idx=0; idx<Sp.size(); idx++) Sp[idx].erase(c_idx);
        candidates_front.erase(c_idx);
	candidates.erase(c_idx);
     }
  } 
}
void CMOEAD::update_lowest_front(unordered_set<int> &candidates_front, unordered_set<int> &survivors_front, vector<unordered_set<int> > &Sp, vector<int> &Np)
{
  for(auto s_idx:survivors_front)
  {
	for(auto idx:Sp[s_idx])
        {
           Np[idx]--;
           if(Np[idx] == 0) candidates_front.insert(idx);
        }
        Np[s_idx]--;
  }
  survivors_front.clear();
}
void CMOEAD::R2_contribution_subset(unordered_set<int> &candidates, unordered_set<int> &candidates_front, vector<int> &survivors, unordered_set<int> &survivors_front, unordered_set<int> &penalized, vector<set<pair<double, int> > > &survivors_weight, vector<double> &distances, vector<vector<double> > &table_fitness, vector<vector<double> > &x_var)
{
  pair<double, int> max_contribution(-DBL_MAX, -1);
  if( survivors_front.empty())
  {
     survivors_weight.assign(nWeight, set<pair<double, int>>());
     for(auto c_idx:candidates_front)
     {
	double sum = 0;
        for(int w_id=0; w_id < nWeight; w_id++)
        {
           sum -= table_fitness[w_id][c_idx];
        }
	if(max_contribution.first < sum) max_contribution = make_pair(sum, c_idx);
     }
    max_contribution.first *=-1; //this is not necessary..
  }
  else
  {
     for(auto c_f:candidates_front)
     {
       double Total_contribution = 0.0;
       for(int w_id = 0; w_id < nWeight; w_id++)
          Total_contribution += max(0.0, survivors_weight[w_id].begin()->first - table_fitness[w_id][c_f]);
       if(Total_contribution > max_contribution.first) max_contribution = make_pair(Total_contribution, c_f);
     }
  }
  survivors_front.insert(max_contribution.second);
  survivors.push_back(max_contribution.second);
  candidates.erase(max_contribution.second);
  candidates_front.erase(max_contribution.second);
  for(int w_id=0; w_id < nWeight; w_id++)  survivors_weight[w_id].insert(make_pair(table_fitness[w_id][max_contribution.second], max_contribution.second));
  for(auto idx_c:candidates) distances[idx_c] = min(distances[idx_c], distance_var(x_var[idx_c], x_var[max_contribution.second]));
  for(auto idx_p:penalized) distances[idx_p] = min(distances[idx_p], distance_var(x_var[idx_p], x_var[max_contribution.second]));
}
void CMOEAD::replacement_phase(CIndividual &parent, CIndividual &child)
{
   vector<vector<double> > x_var, y_obj, table_fitness;
   vector<double> distances;
   vector<unordered_set<int> > Sp;
    vector<int> Np;
   unordered_set<int> candidates, penalized, survivors_front, candidates_front;
   vector<int> survivors;
   vector<set<pair<double, int> > > survivors_weight;

 
     for(int k = 0,idx=0; k < parent.y_obj.size(); k++)
     {
        y_obj.push_back(parent.y_obj[k]);
	x_var.push_back(parent.x_var[k]);
	candidates.insert(idx++);
        y_obj.push_back(child.y_obj[k]);
	x_var.push_back(child.x_var[k]);
	candidates.insert(idx++);
     }
   distances.assign((int)y_obj.size(), DBL_MAX);
   table_fitness.assign(nWeight, vector<double> ((int)y_obj.size())); 
   survivors_weight.assign(nWeight, set<pair<double, int> > ());
   Np.assign((int)y_obj.size(), 0);
   Sp.assign((int)y_obj.size(), unordered_set<int>());
   for(int w_id = 0; w_id < nWeight; w_id++)
     for(auto c_idx:candidates)
        table_fitness[w_id][c_idx] = fitnessfunction(y_obj[c_idx], namda[w_id]);

   dominance_information(Sp, Np, candidates_front, y_obj);

  R2_contribution_subset(candidates, candidates_front, survivors, survivors_front, penalized, survivors_weight, distances, table_fitness, x_var);
  for(auto s_idx:survivors) for(auto c_idx:candidates) distances[c_idx] = min(distances[c_idx], distance_var(x_var[s_idx], x_var[c_idx]));

      while(survivors.size() < nInd)
      {
         //penalize nearest individuals
           if(!candidates.empty()) 
    	penalization(candidates, penalized, distances, candidates_front, Sp, Np);
         if( candidates.empty()) break;
         if(candidates_front.empty())
           update_lowest_front(candidates_front, survivors_front, Sp, Np);
         else 
            R2_contribution_subset(candidates, candidates_front, survivors, survivors_front, penalized, survivors_weight, distances, table_fitness, x_var);
    
      }
      while( survivors.size() < nInd)
      {
          pair<double, int> max_spread = make_pair(-DBL_MAX, -1);
          for( auto p_idx : penalized)
          {
             if( distances[p_idx] > max_spread.first) max_spread = make_pair(distances[p_idx], p_idx);
          }
          survivors.push_back(max_spread.second);
          penalized.erase(max_spread.second);
          for( auto p_idx : penalized)
          {
             distances[p_idx] = min(distances[p_idx], distance_var(x_var[p_idx], x_var[max_spread.second]));
          }
      }
  int idx = 0;
  for(auto s_idx:survivors)
  {
    parent.y_obj[idx] = y_obj[s_idx];
    parent.x_var[idx] = x_var[s_idx];
    idx++; 
  } 
}
#endif
