#ifndef KUHN_MUNKRES
#define KUHN_MUNKRES
using namespace std;
class Hungarian
{
  private:
   const double EPS = 1e-9;
   const double INF = 1e14;
   double *lx, *ly, *slack, *cost;
   int n, max_match, *xy, *yx, *slackx, *prev2;
   bool *S, *T; 
   int *q;


  public:
   Hungarian();
   ~Hungarian();
   void hungarian(double *cost, int *assig);

   void init(int n_nodes);
   void add_to_tree(int x, int prevx);
   void augment();
   void init_labels();
   void update_labels();
};
Hungarian::Hungarian()
{
}
void Hungarian::init(int n_nodes)
{
   n = n_nodes;
   lx =  new double[n_nodes];
   ly =  new double[n_nodes];
   slack =  new double[n_nodes];
//   xy = new int[n_nodes];
   yx = new int[n_nodes];
   slackx = new int[n_nodes];
   prev2 = new int[n_nodes];
   S = new bool[n_nodes];
   T = new bool[n_nodes];
   q = new int[n];
}
Hungarian::~Hungarian()
{
  delete []lx;
  delete []ly;
  delete []slack;
 // delete []xy;
  delete []yx;
  delete []slackx;
  delete []prev2;
  delete []S;
  delete []T;
  delete []q;
}
void Hungarian::hungarian(double *cost_p, int *assig)
{
        xy = assig;
	cost = cost_p;
   	max_match = 0, memset(xy, -1, sizeof(int)*n); 
   	memset(yx, -1, sizeof(int)*n), init_labels(), augment(); //steps 1-3
}
void Hungarian::add_to_tree(int x, int prevx) 
{
   	S[x] = true, prev2[x] = prevx; 
   	for(int y=0; y < n; y++) if (lx[x] + ly[y] - cost[x*n+y] < slack[y] - EPS)
   		slack[y] = lx[x] + ly[y] - cost[x*n+y], slackx[y] = x;
}
void Hungarian::update_labels(){
   	double delta = INF; 
   	for(int y = 0; y < n; y++) if (!T[y]) delta = std::min(delta, slack[y]);
   	for(int x = 0 ; x < n; x++) if (S[x]) lx[x] -= delta;
   	for(int y = 0; y < n; y++) if (T[y]) ly[y] += delta; else slack[y] -= delta;
}
void Hungarian::init_labels()
{
           memset(lx, 0, sizeof(double)*n);
           memset(ly, 0, sizeof(double)*n);
   	for(int x = 0; x < n; x++) for(int y = 0; y < n; y++) lx[x] = std::max(lx[x], cost[x*n+y]);
}
void Hungarian::augment() 
{
   	if (max_match == n) return; 
   	int x, y, root, wr = 0, rd = 0; 
   	memset(S, false, sizeof(bool)*n), memset(T, false, sizeof(bool)*n); 
   	memset(prev2, -1, sizeof(int)*n); 
   	for(int x = 0; x < n; x++) if (xy[x] == -1){
   		q[wr++] = root = x, prev2[x] = -2;
   		S[x] = true; break; }
   	for(int y = 0; y < n; y++) slack[y] = lx[root] + ly[y] - cost[root*n+y], slackx[y] = root;
   	while (true){
   		while (rd < wr){
   			x = q[rd++];
   			for (y = 0; y < n; y++) if (cost[x*n+y] == lx[x] + ly[y] && !T[y]){
   				if (yx[y] == -1) break; T[y] = true; 
   				q[wr++] = yx[y], add_to_tree(yx[y], x); }
   			if (y < n) break; }
   		if (y < n) break; 
   		update_labels(), wr = rd = 0; 
   		for (y = 0; y < n; y++) if (!T[y] && slack[y] == 0){
   			if (yx[y] == -1){x = slackx[y]; break;}
   			else{
   				T[y] = true; 
   				if (!S[yx[y]]) q[wr++] = yx[y], add_to_tree(yx[y], slackx[y]);
   			}}
   		if (y < n) break; }
   	if (y < n){
   		max_match++; 
   		for (int cx = x, cy = y, ty; cx != -2; cx = prev2[cx], cy = ty)
   			ty = xy[cx], yx[cy] = cx, xy[cx] = cy;
   		augment(); }
}
#endif
