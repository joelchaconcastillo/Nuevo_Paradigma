const double EPS = 1e-9;
const double INF = 1e14;
#define N 502 
//Dado un grafo bipartito completo con costos no negativos, encuentra el matching perfecto de minimo costo.
double cost[N][N], lx[N], ly[N], slack[N]; //llenar: cost=matriz de adyacencia
int n, max_match, xy[N], yx[N], slackx[N],prev2[N];//n=cantidad de nodos
bool S[N], T[N]; //sets S and T in algorithm
void add_to_tree(int x, int prevx) {
	S[x] = true, prev2[x] = prevx; 
	for(int y=0; y < n; y++) if (lx[x] + ly[y] - cost[x][y] < slack[y] - EPS)
		slack[y] = lx[x] + ly[y] - cost[x][y], slackx[y] = x;
}
void update_labels(){
	double delta = INF; 
	for(int y = 0; y < n; y++) if (!T[y]) delta = std::min(delta, slack[y]);
	for(int x = 0 ; x < n; x++) if (S[x]) lx[x] -= delta;
	for(int y = 0; y < n; y++) if (T[y]) ly[y] += delta; else slack[y] -= delta;
}
void init_labels(){
        memset(lx, 0, sizeof lx);
        memset(ly, 0, sizeof ly);
	for(int x = 0; x < n; x++) for(int y = 0; y < n; y++) lx[x] = std::max(lx[x], cost[x][y]);
}
void augment() {
	if (max_match == n) return; 
	int x, y, root, q[N], wr = 0, rd = 0; 
	memset(S, false, sizeof(S)), memset(T, false, sizeof(T)); 
	memset(prev2, -1, sizeof(prev2)); 
	for(int x = 0; x < n; x++) if (xy[x] == -1){
		q[wr++] = root = x, prev2[x] = -2;
		S[x] = true; break; }
	for(int y = 0; y < n; y++) slack[y] = lx[root] + ly[y] - cost[root][y], slackx[y] = root;
	while (true){
		while (rd < wr){
			x = q[rd++];
			for (y = 0; y < n; y++) if (cost[x][y] == lx[x] + ly[y] && !T[y]){
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
double hungarian(){
	double ret = 0; max_match = 0, memset(xy, -1, sizeof(xy)); 
	memset(yx, -1, sizeof(yx)), init_labels(), augment(); //steps 1-3
	for(int x = 0; x < n; x++)
  	{
	 ret += cost[x][xy[x]]; 
	}
	return ret;
}
