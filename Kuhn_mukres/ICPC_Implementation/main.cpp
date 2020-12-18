#include <bits/stdc++.h>
#include "km.hpp"
using namespace std;
int main()
{
 n = 3; 
 cost[0][0] =-3;
 cost[0][1] =-9;
 cost[0][2] =-500;

 cost[1][0] =-4;
 cost[1][1] =-3;
 cost[1][2] =-9;

 cost[2][0] =-48;
 cost[2][1] =-22;
 cost[2][2] =-91;
 cout << hungarian()<<endl;
 for(int i = 0; i < n; i++) cout << xy[i]<<endl;

 return 0;
}
