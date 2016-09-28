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

bool vector_contains(vector<Node> r, Node n) {
	for(int i = 0; i < r.size(); i++){
		if (r[i] == n) {
			return true;
		}
	}
	return false;
}


bool is_EdgeBoolMap_true(TSP_Data &tsp, EdgeBoolMap &x) {
	for(EdgeIt o(tsp.g); o!=INVALID; ++o){
		if (!x[o]) {
			return false;
		}
	}
	return true;
}

bool is_NodeBoolMap_true(TSP_Data &tsp, NodeBoolMap &y) {
	for(NodeIt o(tsp.g); o!=INVALID; ++o){
		if (!y[o]) {
			return false;
		}
	}
	return true;
}

void tsp_bt(TSP_Data &tsp, vector<Node> circuit, double &circuitWeightSoFar, double lastWeightAdded, int maxTime, clock_t beginExec, bool &timedOut) {
	clock_t now = clock();
	double elapsed_time = (double) (now-beginExec) / CLOCKS_PER_SEC;
	if (elapsed_time > maxTime) {
		timedOut = true;
		return;
	}

	if((int)circuit.size() == tsp.NNodes){
		double weightFromLastToFirst = 666666.6;

		Node lastNode = circuit.back();
		Node firstNode = circuit.front();
		for (IncEdgeIt e(tsp.g, lastNode); e!=INVALID; ++e) {
			Node ad = tsp.g.v(e);
			if( ad == lastNode ){
				ad = tsp.g.u(e);
			}
			if(ad == firstNode ){
				weightFromLastToFirst = tsp.weight[e];
			}
		}

		if ((circuitWeightSoFar + weightFromLastToFirst) < tsp.BestCircuitValue) {
			tsp.BestCircuit = circuit;
			tsp.BestCircuitValue = circuitWeightSoFar + weightFromLastToFirst;
		}
		return;
	}

	if(circuitWeightSoFar >= tsp.BestCircuitValue){
		circuit.pop_back();
		circuitWeightSoFar -= lastWeightAdded;
		return;
	}

	Node lastNode = circuit.back();
	for (IncEdgeIt e(tsp.g, lastNode); e!=INVALID; ++e) {
		Node ad = tsp.g.v(e);
		if( ad == lastNode ){
			ad = tsp.g.u(e);
		}
		if( !vector_contains(circuit, ad) ){
			circuit.push_back(ad);
			circuitWeightSoFar += tsp.weight[e];
			
			tsp_bt(tsp, circuit, circuitWeightSoFar, tsp.weight[e], maxTime, beginExec, timedOut);
		}
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

	bool timedOut = false;

	vector<Node> circuit = vector<Node>();
	Node firstNode = INVALID;
	string firstNodeStr = ""; 
	for (NodeIt v(tsp.g); v!=INVALID; ++v) {
		if (firstNodeStr.empty() || firstNodeStr > tsp.vname[v]) {
			firstNodeStr = tsp.vname[v];
			firstNode = v;
		}
	}
	circuit.push_back(firstNode);

	double x = 0.0;
	tsp_bt(tsp, circuit, x, 0.0, maxTime, beginExec, timedOut);

	return !timedOut;
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
