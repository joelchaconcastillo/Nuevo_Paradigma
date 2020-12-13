/*
  Author: Joel Chacón Castillo
  Description: 
    This version is taken from https://brc2.com/the-algorithm-workshop/
*/
#include <bits/stdc++.h>
#define EPSILON 1e-4
#define ASTRID 0
#define PRIME 1
#define NONE 2
#define ROW_COVERED 0
#define COL_COVERED 1
using namespace std;
void print(vector< vector<double> > &costs, vector< vector< bool> > &covered_astrid, vector<vector<int> > &state)
{
  int  N = costs.size();
  cout << " ";
  for(int j = 0; j < N; j++) if(covered_astrid[COL_COVERED][j]) cout << "#";else cout <<" ";
  cout <<endl;
  for(int i =0; i < N ; i++)
  {
    if(covered_astrid[ROW_COVERED][i]) cout << "#"; else cout <<" ";
    for(int j = 0; j < N; j++)
    {
	cout << costs[i][j];
	if(state[i][j]==ASTRID) cout << "*";
	else if(state[i][j]==PRIME) cout << "'";
	else cout << " ";
    } 
    cout <<endl;
  }
}
void step_1_2(vector<vector<double> > &costs, vector<vector<int> > &state)
{
  int N = costs.size();
  vector< vector< bool> > covered_astrid(2, vector<bool> (N, false));
  //steps one and two..
  for(int i = 0; i < N ; i++)
  {
    double minv =DBL_MAX ;
    for(int j = 0; j < N; j++) minv = min(minv, costs[i][j]);
    for(int j = 0; j < N; j++)
    {
       costs[i][j] -=minv;
      if(fabs(costs[i][j]) < EPSILON && !covered_astrid[ROW_COVERED][i] && !covered_astrid[COL_COVERED][j])
      state[i][j]=ASTRID, covered_astrid[ROW_COVERED][i]=true, covered_astrid[COL_COVERED][j]=true;
    }
  }
}
int step_3(vector< vector< bool> > &covered_astrid, vector< vector<int> > &state)
{
  int N = state.size();
  int count_cols = 0;
  for(int i = 0; i < N ; i++)//cols..
  {
    if(covered_astrid[COL_COVERED][i]){count_cols++; continue;}
    for(int j = 0; j < N; j++)//rows..
    {
       if(state[j][i] == ASTRID)
       {
        covered_astrid[COL_COVERED][i] = true;
	count_cols++;
	break; //go to the next
       }
    }
  }
  return count_cols;
}
int step_4(pair<int, int> &rc, vector<vector<double> > &costs, vector<vector<int> > &state, vector< vector< bool> > &covered_astrid)
{

  int N = costs.size();
  for(int i = 0; i < N && rc.first == -1; i++) //find an uncovered zero entry
  {
      if(!covered_astrid[ROW_COVERED][i])
      for(int j = 0; j < N && rc.second == -1; j++)
        if( !covered_astrid[COL_COVERED][j] &&fabs(costs[i][j]) < EPSILON) rc = make_pair(i, j);
  }
  if( rc.first != -1)
  {
   state[rc.first][rc.second] = PRIME;
    for(int k=0; k < N; k++)
      if(state[rc.first][k] == ASTRID)
      {
	covered_astrid[ROW_COVERED][rc.first]=true, covered_astrid[COL_COVERED][k]=false; 
	return 4;
      }
    return 5;
   }
  return 6;
}
void step_5(pair<int, int> &rc, vector<vector<int> > &state, vector< vector< bool> > &covered_astrid)
{
  int N = state.size();
  queue<pair<int, int> > path;
  path.push(rc);
  while(true)
  {
    rc.first = -1;
    for(int i=0; i < N; i++) if(state[i][rc.second]==ASTRID){ rc.first=i; break;}
    if(rc.first==-1) break;
    path.push(rc);
    for(int i=0; i < N; i++) if(state[rc.first][i]==PRIME){ rc.second=i; break;}
    path.push(rc);
  }
  while(!path.empty())
  {
     state[path.front().first][path.front().second] = 1-state[path.front().first][path.front().second];
     path.pop();
  }
  covered_astrid[ROW_COVERED].assign(N, false);
  covered_astrid[COL_COVERED].assign(N, false);
  for(int i = 0; i < N*N ; i++) if(state[i/N][i%N]==PRIME) state[i/N][i%N]=NONE;
}
void step_6(vector< vector<double> > &costs, vector< vector< bool> > &covered_astrid)
{

   int N = costs.size();
   double minv = DBL_MAX;
   for(int i = 0; i < N; i++)
   {
      if(covered_astrid[ROW_COVERED][i]) continue;
      for(int j = 0; j < N; j++)
      {
         if(!covered_astrid[COL_COVERED][j])
          minv = min(costs[i][j], minv);
      }
   }

   for(int i = 0; i < N; i++)
   {
      for(int j = 0; j < N; j++)
      {
	if(covered_astrid[ROW_COVERED][i]) costs[i][j] +=minv;
	if(!covered_astrid[COL_COVERED][j]) costs[i][j] -=minv;
      }
   }
}
void KuhnMunkres(vector<int> &assignament, vector<vector<double> > &costs)
{
  int N = (int)costs.size(), cols_covered_count;
  vector<vector<int> > state(N, vector<int>(N, NONE));// 0: astrid,1: prime, 2: none
  vector<vector<bool> > covered_astrid(4, vector<bool>(N, false)); //0: row_astrid, 1:col_astrid, 2:row_covered, 3:col_covered

  step_1_2(costs, state);
  cols_covered_count = step_3(covered_astrid, state);
  while(cols_covered_count < N)
  {
     cols_covered_count = 0; 
     pair<int, int> rc(-1, -1); //row, col
     int nxt_step = step_4(rc, costs, state, covered_astrid);

     if(nxt_step == 5)
     {
      step_5(rc, state, covered_astrid);
      cols_covered_count = step_3(covered_astrid, state);
	cout << cols_covered_count <<endl;
     }
     else if(nxt_step == 6)
     {
        step_6(costs, covered_astrid);
     }
  }
  print(costs, covered_astrid, state);
/////////////////////////////////////////////////////////////////////////
 for(int i = 0; i < N; i++)
   for(int j = 0; j < N; j++)
      if(state[i][j] == ASTRID) assignament[i] = j;
}
int main()
{
  int N = 3;
  vector<int> assignament(N, -1);
  vector<vector<double> > costs(N, vector<double> (N, 0));
  costs[0] = {1, 2, 3};
  costs[1] = {2, 4, 6};
  costs[2] = {3, 6, 9};

  KuhnMunkres(assignament, costs);

  for(int i = 0; i < N ; i++)
    cout << i << " "<< assignament[i]<<endl;
  return 0;
}
