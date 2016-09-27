/*******************************************************************************
 * MC658 - Projeto e Análise de Algoritmos III - 2s2016
 * Prof: Flavio Keidi Miyazawa
 * PED: Mauro Henrique Mulati
 * Usa ideias e código de Rafael Arakaki e Flávio Keidi Miyazawa 
 ******************************************************************************/

/*******************************************************************************
 * EDITE ESTE ARQUIVO APENAS ONDE INDICADO
 * DIGITE SEU RA: 123153
 * SUBMETA SOMENTE ESTE ARQUIVO
 ******************************************************************************/

#include <iostream>
#include <float.h>
#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include "tsp_bt_bnb.h"

bool is_circuit(TSP_Data &tsp, vector<Node> r, double &r_weight) {
	// Number of nodes
	if(tsp.NNodes != (int)r.size()){
		return false;
	}

	// Maps a boolean y to each vertex of the graph g. Initialized as false.
	NodeBoolMap y(tsp.g);
	for(NodeIt o(tsp.g); o!=INVALID; ++o){
		y[o] = false;
	}
	
	// Maps a boolean x to each edge of the graph g. Initialized as false.
	EdgeBoolMap x(tsp.g);
	for(EdgeIt e(tsp.g); e!=INVALID; ++e){  // We set every edge out of the the solution
		x[e] = false;
	}
	
	int nedges = 0;
	
	Node ni = INVALID;
	Node nj = INVALID;
	Node ad = INVALID;
	
	// Lets check each edge inferred by the sequence of nodes.
	for(auto i = r.begin(); i != r.end(); ++i){
		ni = *i;
		// cerr << endl;
		// cerr << "ni: " << tsp.vname[ni] << endl;
		// cerr << "  y[ni]: " << y[ni] << endl;
		
		// If ni is the last element, lets examine the edge (last,first).
		if( next(i) == r.end()){
			nj = *(r.begin());
		}
		else{
			nj = *next(i);
		}
		
		// cerr << "nj: " << tsp.vname[nj] << endl;
		// cerr << "  y[nj]: " << y[nj] << endl;
		
		// Check if nj is in the adjacents of ni, i.e., if exists edge (ni,nj) in the graph. We have to verify all the adjacent edges, since we are in a kind of lists of adjacency structure.
		for(IncEdgeIt e(tsp.g, ni); e != INVALID; ++e){
			ad = tsp.g.v(e);
			if( ad == ni ){
				ad = tsp.g.u(e);
			}
			// cerr << "  ad: " << tsp.vname[ad] << endl;
			
			// cerr << "  x[e]: " << x[e] << "  y[ad]: " << y[ad] << endl;
			
			// If there is an adjacent of ni, which equals nj.
			// Also, it verifies subcycles
			if( ad == nj ){
				if( !x[e] ){  // It is not yet in the solution
					x[e] = true;   // The edge is now verified by the solution
					r_weight += tsp.weight[e];
					nedges++;
					// cerr << "  OK  x[e]: " << x[e];
				}
				else{  // It is already in the solution
					return false;
				}
				
				if( !y[ad] ){  // It is not yet in the solution
					y[ad] = true;  // The adjacent node is now verified by the solution
					// cerr << "  y[ad]: " << y[ad] << endl;
					
				}
				else{  // It is already in the solution
					return false;
				}
			}
		}
	}
	
	// cerr << "nedges: " << nedges << endl;
	// cerr << "tsp.NNodes: " << tsp.NNodes << endl;	
	
	// Number of edges
	if(nedges != tsp.NNodes){
		return false;
	}

	return true;
}

bool vector_contains(vector<Node> r, Node n) {
	Node ni = INVALID;
	for(auto i = r.begin(); i != r.end(); ++i){
		ni = *i;
		if (ni == n) {
			return true;
		}
	}
	return false;
}

void tsp_bt(TSP_Data &tsp, vector<Node> r, int maxTime, clock_t beginExec, bool &timedOut, bool &optimalSolution) {
	clock_t now = clock();
	double elapsed_time = (double) (now-beginExec) / CLOCKS_PER_SEC;
	if (elapsed_time > maxTime) {
		timedOut = true;
		return;
	}

	double new_r_weight = 0.0;

	bool r_is_circuit = is_circuit(tsp, r, new_r_weight);

	if(r_is_circuit){
		tsp.BestCircuit = r;
		tsp.BestCircuitValue = new_r_weight;
		return;
	}

	if(new_r_weight >= tsp.BestCircuitValue){
		r.pop_back();
		return;
	}

	for(NodeIt v(tsp.g); v != INVALID; ++v){
		if(!vector_contains(r, v) && !timedOut) {
			r.push_back(v);
			bool x = false;
			tsp_bt(tsp, r, maxTime, beginExec, timedOut, x);
		}
	}
	if (!timedOut) {
		optimalSolution = true;	
	}
}

bool bt(TSP_Data &tsp, int maxTime)
/*******************************************************************************
 * SUBSTITUIA O CONTEÚDO DESTE MÉTODO POR SUA IMPLEMENTAÇÃO DE BACKTRACKING.
 * ENTRETANTO, NÃO ALTERE A ASSINATURA DO MÉTODO.
 ******************************************************************************/
{
	clock_t beginExec = clock();

	greedy(tsp, maxTime);

	clock_t now = clock();
	double elapsed_time = (double) (now-beginExec) / CLOCKS_PER_SEC;
	cerr << endl;
	cerr << "greedyTime : " << elapsed_time;
	cerr << endl;

	bool optimalSolution = false;
	bool timedOut = false;

	tsp_bt(tsp, vector<Node>(), maxTime, beginExec, timedOut, optimalSolution);

	cerr << endl;
	cerr << "timedOut : " << timedOut << " | optimalSolution : " << optimalSolution << "  ";
	cerr << endl;

	return optimalSolution;

}

bool bnb(TSP_Data &tsp,  int maxTime)
/*******************************************************************************
 * SUBSTITUIA O CONTEÚDO DESTE MÉTODO POR SUA IMPLEMENTAÇÃO DE BRANCH AND BOUND.
 * ENTRETANTO, NÃO ALTERE A ASSINATURA DO MÉTODO.
 ******************************************************************************/
{
	return greedy(tsp, maxTime);
}

bool greedy(TSP_Data &tsp,  int maxTime)
/*******************************************************************************
 * Algoritmo guloso para o TSP.
 ******************************************************************************/
{
	// Maps a boolean x to each edge of the graph g. Initialized as false.
	EdgeBoolMap x(tsp.g);
	for(EdgeIt e(tsp.g); e!=INVALID; ++e){  // We set every edge out of the the solution
		x[e] = false;
	}
	
	// NodeBoolMap y(tsp.g, false);  // Maps a boolean y to each vertex of the graph g. Initialized as false.
	NodeBoolMap y(tsp.g);
	for(NodeIt o(tsp.g); o!=INVALID; ++o){
		y[o] = false;
	}
	
	double cost = 0.0;
	int nedges = 0;
	
	// cerr << endl;
	// for(NodeIt o(tsp.g); o != INVALID; ++o){
	// 	cerr << tsp.g.id(o) << ":" << y[o] << "  ";
	// }
	// cerr << endl;

	NodeIt nit(tsp.g);
	Node n = nit;
	Node f = nit;
	
	while(nedges != tsp.NNodes){
		// Put the vertex in the solution
		y[n] = true;
		tsp.BestCircuit.push_back(n);
		
		// cerr << "n: " << tsp.g.id(n) << "  y[n]: " << y[n] << endl;
		
		double wmin = DBL_MAX;  // min weight
		IncEdgeIt emin = INVALID;  // min inc edge of n
		Node nmin = INVALID;
		
		IncEdgeIt e(tsp.g, n);
		Node op = INVALID;
		
		// cerr << "wmin: " << wmin << endl;
		
		for(; e != INVALID; ++e){

			
			op = tsp.g.v(e);
			if( op == n ){
				op = tsp.g.u(e);
			}
			
			// cerr << "   (" << tsp.g.id(tsp.g.u(e)) << ", " << tsp.g.id(tsp.g.v(e)) << ")  x: " << x[e] << "  c: " << tsp.weight[e] << " op: " << tsp.g.id(op) << " y[op] " << y[op] << endl;
			
			if( ! y[ op ] ){        // The expression in [] returns the "destin" vertex of edge e
				if( tsp.weight[e] < wmin ){
					wmin = tsp.weight[e];
					emin = e;
					nmin = op;
				}
			}
		}
		
		if( wmin < DBL_MAX ){  // If got some edge
			// cerr << "wmin: " << wmin << endl;
			x[emin] = true;  // Puts the edge e in the solution: this data will be visible outside this function
			nedges++;
			cost += wmin;
			n = nmin;
			// cerr << "new n: " << tsp.g.id(n) << endl;
		}
		else{
			cout << "Error: could not found a minimum weight value." << endl;
			exit(1);
		}
		
		if( nedges == tsp.NNodes - 1 ){
			y[f] = false;
		}
		
		// cerr << "nedges: " << nedges << endl;
		// cerr << endl;
	}
	
	if( nedges > 0 ){
		tsp.BestCircuitValue = cost;
	}

	return false;
}
