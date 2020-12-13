#include <bits/stdc++.h>
#include "km.hpp"
using namespace std;
int main()
{
 n = 3; 
 cost[0][0] =-1;
 cost[0][1] =-2;
 cost[0][2] =-3;

 cost[1][0] =-2;
 cost[1][1] =-4;
 cost[1][2] =-6;

 cost[2][0] =-3;
 cost[2][1] =-6;
 cost[2][2] =-9;
 cout << hungarian()<<endl;
 for(int i = 0; i < n; i++) cout << xy[i]<<endl;

 return 0;
}
