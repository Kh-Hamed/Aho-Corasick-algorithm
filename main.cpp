// C++ program for implementation of Aho Corasick algorithm
// for string matching
using namespace std;
#include <bits/stdc++.h>

// Max number of states in the matching machine.
// Should be equal to the sum of the length of all keywords.
const long long int MAXS = 550;

// Maximum number of characters in input alphabet
const long long int MAXC = 26;

// OUTPUT FUNCTION IS IMPLEMENTED USING out[]
// Bit i in this mask is one if the word with index i
// appears when the machine enters this state.
double out[MAXS];
long long int V=1;
// FAILURE FUNCTION IS IMPLEMENTED USING f[]
long long int f[MAXS];

// GOTO FUNCTION (OR TRIE) IS IMPLEMENTED USING g[][]
long long int g[MAXS][MAXC];
long long int adjacent[MAXS][MAXS];
// Builds the string matching machine.
// arr - array of words. The index of each keyword is important:
//		 "out[state] & (1 << i)" is > 0 if we just found word[i]
//		 in the text.
// Returns the number of states that the built machine has.
// States are numbered 0 up to the return value - 1, inclusive.
long long int buildMatchingMachine(char arr[][MAXS],long long int k,double pt[])
{
	// Initialize all values in output function as 0.
	memset(out, 0, sizeof out);

	// Initialize all values in goto function as -1.
	memset(g, -1, sizeof g);

	// Initially, we just have the 0 state
	long long int states = 1;
	long long int i=0,j=0;
	long long int hold=0;
	long long int fl=0;
	long long int ch=0;
	long long int state;
	long long int failure=0;
	long long int t=0;

	// Construct values for goto function, i.e., fill g[][]
	// This is same as building a Trie for arr[]
	for (i = 0; i < k; ++i)
	{
		//const string &word = arr[i];
		long long int currentState = 0;
		char word[MAXS];
		t=0;
		while(arr[i][t]!=0){
			word[t]=arr[i][t];
			t=t+1;
		}
		

		// Insert all characters of current word in arr[]
		for ( j = 0; j < t; ++j)
		{
			ch = word[j] - 'a';

			// Allocate a new node (create a new state) if a
			// node for ch doesn't exist.
			if (g[currentState][ch] == -1){
				g[currentState][ch] = states++;
				adjacent[currentState][g[currentState][ch]]=1;
				V+=1;
				
			}
			currentState = g[currentState][ch];
		}

		// Add current word in output function
		if(out[currentState] !=0)
		out[currentState] += pt[i];
		else
		out[currentState] = pt[i];
	}
	out[0]=0;

	// For all characters which don't have an edge from
	// root (or state 0) in Trie, add a goto edge to state
	// 0 itself
	for ( ch = 0; ch < MAXC; ++ch){
		if (g[0][ch] == -1){
			g[0][ch] = 0;
		}
	}
	adjacent[0][0]=1;
			

	// Now, let's build the failure function

	// Initialize values in fail function
	memset(f, -1, sizeof f);

	// Failure function is computed in breadth first order
	// using a queue
	queue<long long int> q;

	// Iterate over every possible input
	for ( ch = 0; ch < MAXC; ++ch)
	{
		// All nodes of depth 1 have failure function value
		// as 0. For example, in above diagram we move to 0
		// from states 1 and 3.
		if (g[0][ch] != 0)
		{	
			//tp=g[0][ch];
			f[g[0][ch]] = 0;
			q.push(g[0][ch]);
		}
	}

	// Now queue has states 1 and 3
	while (q.size())
	{
		// Remove the front state from queue
		state = q.front();
		q.pop();

		// For the removed state, find failure function for
		// all those characters for which goto function is
		// not defined.
		for (ch = 0; ch < MAXC; ++ch)
		{
			// If goto function is defined for character 'ch'
			// and 'state'
			if (g[state][ch] != -1)
			{
				// Find failure state of removed state
				failure = f[state];

				// Find the deepest node labeled by proper
				// suffix of string from root to current
				// state.
				while (g[failure][ch] == -1)
					failure = f[failure];

				failure = g[failure][ch];
				f[g[state][ch]] = failure;

				// Merge output values
				out[g[state][ch]] += out[failure];

				// Insert the next level node (of Trie) in Queue
				q.push(g[state][ch]);
			}
		}
	}

for ( ch = 0; ch < MAXC; ++ch)
{
	if (g[0][ch] != 0)
	{
		for ( i = 0; i < MAXC; ++i)
		{
			hold=g[0][ch];
			if(g[hold][i]==-1){
				fl = f[hold];
				while (g[fl][i] == -1){
						fl = f[fl];
				}
				adjacent[hold][g[fl][i]]=1;
			}
				// Insert the next level node (of Trie) in Queue
		}
			q.push(g[0][ch]);
	}
}

while (q.size())
	{
		// Remove the front state from queue
		state = q.front();
		q.pop();
		// For the removed state, find failure function for
		// all those characters for which goto function is
		// not defined.
		
		for ( ch = 0; ch < MAXC; ++ch)
		{
			if (g[state][ch] != -1){	
				for ( i = 0; i < MAXC; ++i)
				{
					hold=g[state][ch];
				if(g[hold][i]==-1){
					fl = f[hold];
					while (g[fl][i] == -1){
						fl = f[fl];
					}
					adjacent[hold][g[fl][i]]=1;
				}
				// Insert the next level node (of Trie) in Queue
				}
			q.push(g[state][ch]);
			
			}	
		}
	}
	return states;
}

struct Edge {
    long long int src, dest;
	double weight;
};
 
// a structure to represent a connected, directed and
// weighted graph
struct Graph {
    // V-> Number of vertices, E-> Number of edges
    long long int V, E;
 
    // graph is represented as an array of edges.
    struct Edge* edge;
};
 
// Creates a graph with V vertices and E edges
struct Graph* createGraph(long long int V,long long int E)
{
    struct Graph* graph = new Graph;
    graph->V = V;
    graph->E = E;
    graph->edge = new Edge[graph->E];
    return graph;
}
 
// The main function that finds shortest distances
// from src to all other vertices using Bellman-
// Ford algorithm.  The function also detects
// negative weight cycle
bool isNegCycleBellmanFord(struct Graph* graph,
                           long long int src)
{
    long long int V = graph->V;
    long long int E = graph->E;
    double dist[V];
    long long int i=0;
    long long int j =0;
 
    // Step 1: Initialize distances from src
    // to all other vertices as INFINITE
    for (long long int i = 0; i < V; i++)
        dist[i] = DBL_MAX;
    dist[src] = 0;
 
    // Step 2: Relax all edges |V| - 1 times.
    // A simple shortest path from src to any
    // other vertex can have at-most |V| - 1
    // edges
    for ( i = 0; i <= V - 1; i++) {
        for (j = 0; j < E; j++) {
            long long  int u = graph->edge[j].src;
            long long int v = graph->edge[j].dest;
            double weight = graph->edge[j].weight;
            if (dist[u] + weight < dist[v])
                dist[v] = dist[u] + weight;
        }
    }
 
    // Step 3: check for negative-weight cycles.
    // The above step guarantees shortest distances
    // if graph doesn't contain negative weight cycle.
    // If we get a shorter path, then there
    // is a cycle.
    for ( i = 0; i < E; i++) {
        long long int u = graph->edge[i].src;
        long long int v = graph->edge[i].dest;
        double weight = graph->edge[i].weight;
        if (dist[u] + weight < dist[v])
            return true;
    }
 
    return false;
}

// Driver program to test above
int main()
{
	memset(adjacent, 0, sizeof adjacent);
	double point[MAXS];
	long long int n=0,i=0,k=0;
	long long int E=0;
	long long int u=0,cnt=1;
	double min=DBL_MAX;
	double max=-DBL_MAX;
	double initial=0;
	scanf("%lld",&n); 
	char arr[n][MAXS];
	memset(arr, 0, sizeof arr);
	for( i=0;i<n;++i){
		scanf("%lf",&point[i]);
	}
	for( i=0;i<n;++i){
		cin >> arr[i];
	}
	buildMatchingMachine(arr, n,point);
	for(u=0;u<V;++u){
		if(min>=out[u])
			min=out[u];
		if(max<=out[u])
			max=out[u];
	}
	for(u=0;u<V;++u){
		for(k=0;k<V;++k){
			if(adjacent[u][k]!=0){
				E+=1;
			}
		}
	}
	double up_ban=max;
	double low_ban=min;
	double hold=0;
	struct Graph* graph = createGraph(V, E);
	i=0;
	initial=(up_ban+low_ban)/2;
	hold=initial;
	for(u=0;u<V;++u){
		for(k=0;k<V;++k){
			if(adjacent[u][k]!=0){
				graph->edge[i].src = u;
	    		graph->edge[i].dest = k;
	    		graph->edge[i].weight = initial -out[k];
	    		i+=1;
    		}
		}
	}
	while(cnt<=50)
	{
	i=0;
	initial=(up_ban+low_ban)/2;
	for(u=0;u<V;++u){
		for(k=0;k<V;++k){
			if(adjacent[u][k]!=0){
	    		graph->edge[i].weight = initial -out[k];
	    		i+=1;
    		}
		}
	}
	if (isNegCycleBellmanFord(graph,0)){
		low_ban=(up_ban+low_ban)/2;
	}
	else{
		up_ban=(up_ban+low_ban)/2;
	}
	cnt+=1;
}
printf("%0.6f",(up_ban+low_ban)/2);		
	return 0;
}
