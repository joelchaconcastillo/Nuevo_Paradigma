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
	double distance_var( vector<double> &a, vector<double> &b);
	void penalization(unordered_set<int> &candidates, unordered_set<int> &penalized, vector<double> &distances, unordered_set<int> &candidates_front);
        void dominance_information(); 
        void update_fronts(); 
	void pick_penalized(unordered_set<int> &penalized, unordered_set<int> &survivors, vector<double> &distances);
	void R2_contribution_subset(unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &survivors, unordered_set<int> &survivors_front, unordered_set<int> &penalized, vector<set<pair<double, int> > > &survivors_weight, vector<double> &distances);
	void update_lowest_front(unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &survivors_front);
	void diversity_information(unordered_set<int> &survivors, unordered_set<int> &penalized, vector<double> &distances);

	vector<set<int> > non_dominated_sorting(vector<vector<double> > &y_obj);
        void table_fitness_information();
	void eval_R2(CIndividual &indiv);
	void updateTarget(CIndividual &target, CIndividual &trial);

   private:
	vector <CIndividual> pool;
	vector<int> child_idx, parent_idx, inv_parent_idx;
	vector<vector<double> > namda;     // weight vector
        vector<int> Np,Rp;//rank
	vector<unordered_set<int> > Sp;//dominated indexes and inverse
        vector<vector<int> > fronts ;
        vector<vector<double > > table_fitness;
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

//    table_fitness.assign(nWeight, vector<double> (nPop+nOffspring));

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
    for(int i=0; i< nPop+nOffspring; i++) eval_R2(pool[i]);
     readf.close( );
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
      diff_evo_xoverA(pool[parent_idx[idx_target]], pool[parent_idx[idx1]], pool[parent_idx[idx2]], pool[parent_idx[idx3]], child, CR, F);
      // apply polynomial mutation
   //   realmutation(child, 1.0/(nvar*nInd));
      child.obj_eval();
      update_reference(child); //O(M)

//      updateTarget(pool[parent_idx[idx_target]], child);
      eval_R2(child);
      if( child.fitness < pool[parent_idx[idx_target]].fitness) pool[parent_idx[idx_target]] =  child;
	//cout << child.fitness<<" __ " << pool[parent_idx[idx_target]].fitness <<endl;
   }
//   replacement_phase();
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
	//	update_parameterD();
		evol_population();
		accumulator += nfes - bef ;
      //          if(accumulator > 0.01*(max_nfes)  )
		{
	           accumulator -= 0.01*(max_nfes);
	//	   save_pos(filename1);
		   save_front(filename2);
		}
		bef=nfes;
	        nfes += nOffspring;
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
   std::fstream fout; //fout.open(saveFilename,std::ios::out);
   fout.open(saveFilename, fstream::app|fstream::out);
   for(int n=0; n<nPop; n++)
   {
      for(int k=0;k<nvar;k++)
         fout<<pool[parent_idx[n]].x_var[k] << "  ";
      for(int k=0;k<nvar;k++)
   	 fout<<pool[parent_idx[n]].x_var[k]/(vuppBound[k]-vlowBound[k]) << "  "; //fout<<population[n].indiv.x_var[k]<< fixed << setprecision(30) << "  ";
   	fout<<"\n";
   }
   fout.close();
}
void CMOEAD::update_fronts()
{
   vector<int> current_Np = Np;
   fronts.assign(1, vector<int>());
   int rank = 0;
   for(auto idx:parent_idx)
   {
      if(current_Np[idx]==0)
      {
        fronts[rank].push_back(idx);
        Rp[idx] = rank;
      }
   }
   while(true)
   {
      vector<int> next_front;
      for(auto idx:fronts[rank])
      {
	for(auto idx_dominated:Sp[idx])
        {
	  current_Np[idx_dominated]--;
          if(current_Np[idx_dominated]  == 0) 
          {
	     next_front.push_back(idx_dominated);
	     Rp[idx_dominated] = rank+1;
          }
        }
      }
      if(next_front.empty()) break;
      fronts.push_back(next_front);
      rank++;
   }
}
void CMOEAD::dominance_information()
{
   Sp.assign(nPop+nOffspring, unordered_set<int>());
   Np.assign(nPop+nOffspring, 0);
   Rp.assign(nPop+nOffspring, 0);
   fronts.assign(1, vector<int>());
   int rank = 0;
   for(int pidx1=0; pidx1< nPop+nOffspring; pidx1++)
   {
      for(int pidx2=0; pidx2< nPop+nOffspring; pidx2++)
      {
	if(pidx1 == pidx2) continue;
        if( pool[pidx1] < pool[pidx2]) Sp[pidx1].insert(pidx2);
 	else if( pool[pidx2] < pool[pidx1]) Np[pidx1]++;
      }
      if( Np[pidx1] == 0)
      {
         fronts[rank].push_back(pidx1);
	 Rp[pidx1]=rank;
      }
   }
}
void CMOEAD::table_fitness_information()
{
  //getting R2-contribution of each individual...
//  for(int w_id = 0; w_id < nWeight; w_id++)
//  {
//     for(int idx = 0; idx < pool.size(); idx++)
//     {
//        table_fitness[w_id][idx] = fitnessfunction(pool[idx].y_obj, namda[w_id]);
//     }
//  }
}
void CMOEAD::diversity_information(unordered_set<int> &survivors, unordered_set<int> &candidates, vector<double> &distances)
{
   for(auto s_idx:survivors)
   {
      for(auto c_idx:candidates)
      {
	distances[c_idx] = min(distances[c_idx], distance_var(pool[s_idx].x_var, pool[c_idx].x_var));
      }
   } 
}
void CMOEAD::penalization(unordered_set<int> &candidates, unordered_set<int> &penalized, vector<double> &distances, unordered_set<int> &candidates_front)
{
  if(candidates.empty()) return;
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
void CMOEAD::pick_penalized(unordered_set<int> &penalized, unordered_set<int> &survivors, vector<double> &distances)
{
   pair<double, int> max_spread = make_pair(-DBL_MAX, -1);
   for( auto p_idx : penalized)
   {
      if( distances[p_idx] > max_spread.first) max_spread = make_pair(distances[p_idx], p_idx);
   }
   survivors.insert(max_spread.second);
   penalized.erase(max_spread.second);
   for( auto p_idx : penalized)
   {
      distances[p_idx] = min(distances[p_idx], distance_var( pool[p_idx].x_var, pool[max_spread.second].x_var));
   }
}
void CMOEAD::update_lowest_front(unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &survivors_front)
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
void CMOEAD::R2_contribution_subset(unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &survivors, unordered_set<int> &survivors_front, unordered_set<int> &penalized, vector<set<pair<double, int> > > &survivors_weight, vector<double> &distances)
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
  survivors.insert(max_contribution.second);
  candidates.erase(max_contribution.second);
  candidates_front.erase(max_contribution.second);
  for(int w_id=0; w_id < nWeight; w_id++)  survivors_weight[w_id].insert(make_pair(table_fitness[w_id][max_contribution.second], max_contribution.second));
  for(auto idx_c:candidates) distances[idx_c] = min(distances[idx_c], distance_var(pool[idx_c].x_var, pool[max_contribution.second].x_var));
  for(auto idx_p:penalized) distances[idx_p] = min(distances[idx_p], distance_var(pool[idx_p].x_var, pool[max_contribution.second].x_var));
}
void CMOEAD::replacement_phase()
{

  unordered_set<int> survivors, candidates, penalized, survivors_front, candidates_front;
  vector<double> distances((int)pool.size(), DBL_MAX);
  vector<set<pair<double, int> > > survivors_weight;
  for(int idx = 0; idx < pool.size(); idx++) candidates.insert(idx); 
  table_fitness_information();
  dominance_information(); 

   for(auto c_idx:candidates) if( Np[c_idx] == 0) candidates_front.insert(c_idx);

  R2_contribution_subset(candidates, candidates_front, survivors, survivors_front, penalized, survivors_weight, distances);
  diversity_information(survivors, candidates, distances);

  while( survivors.size() < nPop)
  {  
     penalization(candidates, penalized, distances, candidates_front);

     if( candidates.empty()) break;

     if(candidates_front.empty())
       update_lowest_front(candidates, candidates_front, survivors_front);
     else 
        R2_contribution_subset(candidates, candidates_front, survivors, survivors_front, penalized, survivors_weight, distances);
  }

  while( survivors.size() < nPop )
  {
       pick_penalized(penalized, survivors, distances);
  }

  int idx = 0; 
  vector<bool> setted( pool.size(), false);
  for(auto ite_s = survivors.begin() ; ite_s != survivors.end(); ite_s++, idx++)
  {
    parent_idx[idx] = *ite_s;
    setted[*ite_s] = true;
  }
  idx = 0;
  for(int i = 0; i < pool.size();  i++) if(!setted[i]) child_idx[idx++]=i;


}
void CMOEAD::eval_R2(CIndividual &indiv)
{
  vector<set<int> > fronts = non_dominated_sorting(indiv.y_obj);
  double totalsum = 0.0; 
  //adding earch rank...
//  for(int r = 0; r < fronts.size(); r++)
  {
	int r = 0;// fronts.size()-1;
     for(int w = 0; w < nWeight; w++)
     {
       double minv = DBL_MAX;
       for(auto k:fronts[r])
       {
           minv = min(minv, fitnessfunction(indiv.y_obj[k], namda[w]));
       } 
       totalsum += minv;
     }
  }
  indiv.fitness = totalsum;
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
void CMOEAD::updateTarget(CIndividual &target, CIndividual &trial)
{
   //merging...
   vector<double> x_var;
   vector<vector<double> > y_obj;
   for(auto x_value:target.x_var) x_var.push_back(x_value);
   for(auto x_value:trial.x_var) x_var.push_back(x_value);
   for(auto y_vec:target.y_obj) y_obj.push_back(y_vec);
   for(auto y_vec:trial.y_obj) y_obj.push_back(y_vec);

   vector<set<int> > fronts = non_dominated_sorting(y_obj); 
   set<int> selected_points;
   int rank = 0, total_size=0;
   //selecting taking into account the rank..
   while(selected_points.size() < nInd)
   {
     for(auto idx:fronts[rank]) selected_points.insert(idx);
     rank++;
   } 
   rank--; 

   vector< set<pair<double, int> > > sorted_fitness(nWeight);
   vector<double> contribution(y_obj.size(), 0.0);
   for(int w = 0; w < nWeight; w++)
   { 
      for(auto idx:fronts[rank])
      {	
        sorted_fitness[w].insert(make_pair(fitnessfunction(y_obj[idx], namda[w]), idx));
      }
      contribution[sorted_fitness[w].begin()->second] += next(sorted_fitness[w].begin(),1)->first - sorted_fitness[w].begin()->first;
   }
   //removing the worst individuals based in the R2-indicator..
   while(selected_points.size() > nInd)
   {
     pair<double, int> min_c(DBL_MAX, -1);
      for(auto idx:fronts[rank]) if(min_c.first > contribution[idx] ) min_c = make_pair(contribution[idx], idx);
     //remove it..
     selected_points.erase(min_c.second);
     fronts[rank].erase(min_c.second);
     //update it..
     contribution.assign(y_obj.size(), 0.0);
     for(int w = 0; w < nWeight; w++)
     { 
      sorted_fitness[w].erase( make_pair(fitnessfunction(y_obj[min_c.second], namda[w]), min_c.second) );
      contribution[sorted_fitness[w].begin()->second] += next(sorted_fitness[w].begin(),1)->first - sorted_fitness[w].begin()->first;
     } 
   }    
  //update target..
  int ite = 0;
  for(auto idx:selected_points)
  {
    target.y_obj[ite] = y_obj[idx];
    for(int i = 0; i < nvar; i++) target.x_var[ite*nvar + i] = x_var[idx*nvar + i];
    ite++;
  }
}
#endif
