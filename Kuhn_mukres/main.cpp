/*
  Author: Joel Chacón Castillo
  Description: 
    This version is taken from https://brc2.com/the-algorithm-workshop/
*/
#include <bits/stdc++.h>
#include "Kuhn_Munkres.h"
int main()
{
  int N = 3;
  vector<int> assignament(N, -1);
  vector<vector<double> > costs(N, vector<double> (N, 0));
  costs[0] = {1, 2, 3};
  costs[1] = {2, 4, 6};
  costs[2] = {3, 6, 9};

  KuhnMunkres(assignament, costs);
  cout << "Assignament.."<<endl;
  for(int i = 0; i < N ; i++)
    cout << i << " "<< assignament[i]<<endl;
  return 0;
}
