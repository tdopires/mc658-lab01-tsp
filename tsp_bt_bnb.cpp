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


/************************************************************************************************************************************************************/

string indent(int depth) {
	stringstream res;
	for (int j = 0; j < depth; j++) {
		res << "    ";
	}
	return res.str();
}

int find_set_index(vector<vector<Node>> sets, Node node) {
	for (int i = 0; i < sets.size(); i++) {
		for (int j = 0; j < sets[i].size(); j++) {
			if (sets[i][j] == node) {
				return i;
			}
		}
	}
	return -1;
}

bool isEdgeBoolMap_AllTrue(TSP_Data &tsp, EdgeBoolMap &x) {
	for(EdgeIt e(tsp.g); e!=INVALID; ++e){
		if (!x[e]) {
			return false;
		}
	}
	return true;
}

vector<Edge> orderEdgesByWeight(TSP_Data &tsp) {
	vector<Edge> orderedEdges = vector<Edge>();

	EdgeBoolMap x(tsp.g);
	for(EdgeIt e(tsp.g); e!=INVALID; ++e){
		x[e] = false;
	}

	while (!isEdgeBoolMap_AllTrue(tsp, x)) {
		Edge minEdge;
		double minWeight = DBL_MAX;
		for(EdgeIt e(tsp.g); e!=INVALID; ++e){
			if (tsp.weight[e] < minWeight && !x[e]) {
				minWeight = tsp.weight[e];
				minEdge = e;
			}
		}
		x[minEdge] = true;
		orderedEdges.push_back(minEdge);
	}
	return orderedEdges;
}

// Algoritmo de Kruskal
double calculate_mst_weight(TSP_Data &tsp, vector<Node> nodesToIgnore, vector<Edge> orderedEdges, int maxTime, clock_t beginExec, bool &timedOut) {
	double weight = 0.0;
	vector<Edge> result = vector<Edge>();
	vector<vector<Node>> sets = vector<vector<Node>>();
	
	for (NodeIt v(tsp.g); v!=INVALID; ++v) {
		if ( vector_contains(nodesToIgnore, v) )
			continue;

		vector<Node> newset = vector<Node>();
		newset.push_back(v);
		sets.push_back(newset);

		if (((double)((clock() - beginExec) / CLOCKS_PER_SEC)) > maxTime) {
			timedOut = true;
			return weight;
		}
	}

	for (int i = 0; i < orderedEdges.size(); i++) {
		Edge e = orderedEdges[i];
		Node u = tsp.g.u(e);
		Node v = tsp.g.v(e);
		if ( vector_contains(nodesToIgnore, u) || vector_contains(nodesToIgnore, v) )
			continue;

		int u_index = find_set_index(sets, u);
		int v_index = find_set_index(sets, v);
		if ( u_index != v_index ) {
			result.push_back(e);
			weight += tsp.weight[e];

			sets[u_index].insert(sets[u_index].end(), sets[v_index].begin(), sets[v_index].end());
			sets[v_index] = vector<Node>();
		}
		if (((double)((clock() - beginExec) / CLOCKS_PER_SEC)) > maxTime) {
			timedOut = true;
			return weight;
		}
	}
	
	return weight;
}

double smallest_weight_on_edges(TSP_Data &tsp, Node node) {
	double smallest_weight = DBL_MAX;
	for (IncEdgeIt e(tsp.g, node); e!=INVALID; ++e) {
		Node ad = tsp.g.v(e);
		if( ad == node ) {
			ad = tsp.g.u(e);
		}

		if (tsp.weight[e] < smallest_weight) {
			smallest_weight = tsp.weight[e];
		}
	}
	return smallest_weight;
}
double smallest_weight_on_edges_excluding(TSP_Data &tsp, Node node, Node toExclude) {
	double smallest_weight = DBL_MAX;
	for (IncEdgeIt e(tsp.g, node); e!=INVALID; ++e) {
		Node ad = tsp.g.v(e);
		if( ad == node ) {
			ad = tsp.g.u(e);
		}
		if (ad == toExclude) {
			continue;
		}

		if (tsp.weight[e] < smallest_weight) {
			smallest_weight = tsp.weight[e];
		}
	}
	return smallest_weight;
}

double calculate_weight_between_circuit_and_rest(TSP_Data &tsp, vector<Node> circuit) {
	if ((int)circuit.size() == 1) {
		Node firstNode = circuit[0];
		return smallest_weight_on_edges(tsp, firstNode);
	}

	Node firstNode = circuit[0];
	Node secondNode = circuit[1];
	double w1 = smallest_weight_on_edges_excluding(tsp, firstNode, secondNode);

	Node beforeLastNode = circuit[((int)circuit.size())-2];
	Node lastNode = circuit[((int)circuit.size())-1];
	double w2 = smallest_weight_on_edges_excluding(tsp, beforeLastNode, lastNode);
	return w1 + w2;

}


void tsp_bnb(TSP_Data &tsp, vector<Edge> orderedEdges, double weightSoFar, Node node, vector<Node> circuit, int depth, int maxTime, clock_t beginExec, bool &timedOut) {
	clock_t now = clock();
	double elapsed_time = (double) (now-beginExec) / CLOCKS_PER_SEC;
	if (elapsed_time > maxTime) {
		timedOut = true;
		return;
	}
	
	if ((int)circuit.size() == tsp.NNodes ) { // is a leaf and we have a circuit
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

		if ((weightSoFar + weightFromLastToFirst) <= tsp.BestCircuitValue) {
			// cerr << endl;
			// cerr << indent(depth) << "( graphNode: " << tsp.vname[node]; 
			// cerr << " || weightSoFar: " << weightSoFar;
			// cerr << " || weightFromLastToFirst: " << weightFromLastToFirst;
			// cerr << " =====> NEW BEST CIRCUIT VALUE: " << (weightSoFar + weightFromLastToFirst) << " ) ";

			tsp.BestCircuit = circuit;
			tsp.BestCircuitValue = weightSoFar + weightFromLastToFirst;
		}

	} else {
		// cerr << endl;
		// cerr << indent(depth) << "( node: " << tsp.vname[node]; 
		// cerr << " || weightSoFar: " << weightSoFar;
		// cerr << endl;


		for (IncEdgeIt e(tsp.g, node); e!=INVALID; ++e) {
			Node ad = tsp.g.v(e);
			if( ad == node ){
				ad = tsp.g.u(e);
			}

			if ( !vector_contains(circuit, ad) ) {

				vector<Node> oldCircuit = circuit;
				circuit.push_back(ad);

				double weightChildMst = calculate_mst_weight(tsp, circuit, orderedEdges, maxTime, beginExec, timedOut);

				double weightcircuitBetweenRest = calculate_weight_between_circuit_and_rest(tsp, circuit);

				// cerr << endl;
				// cerr << indent(depth) << " => going to node: " << tsp.vname[ad] << " mst do resto " << weightChildMst; 
				
				double newWeightSoFar = weightSoFar + tsp.weight[e];

				if ( (newWeightSoFar + weightChildMst + weightcircuitBetweenRest) < tsp.BestCircuitValue ) {

					tsp_bnb(tsp, orderedEdges, newWeightSoFar, ad, circuit, depth + 1, maxTime, beginExec, timedOut);

				}

				circuit = oldCircuit;
			}

			if (((double)((clock() - beginExec) / CLOCKS_PER_SEC)) > maxTime) {
				timedOut = true;
				return;
			}
		}

		// cerr << " ) ";

	}
}

bool bnb(TSP_Data &tsp,  int maxTime)
/*******************************************************************************
 * SUBSTITUIA O CONTEÚDO DESTE MÉTODO POR SUA IMPLEMENTAÇÃO DE BRANCH AND BOUND.
 * ENTRETANTO, NÃO ALTERE A ASSINATURA DO MÉTODO.
 ******************************************************************************/
{
	clock_t beginExec = clock();

	greedy(tsp, maxTime);

	bool timedOut = false;

	Node firstNode = INVALID;
	string firstNodeStr = ""; 
	for (NodeIt v(tsp.g); v!=INVALID; ++v) {
		if (firstNodeStr.empty() || firstNodeStr > tsp.vname[v]) {
			firstNodeStr = tsp.vname[v];
			firstNode = v;
		}
	}

	vector<Node> circuit = vector<Node>();
	circuit.push_back(firstNode);

	vector<Edge> orderedEdges = orderEdgesByWeight(tsp);
	tsp_bnb(tsp, orderedEdges, 0.0, firstNode, circuit, 0, maxTime, beginExec, timedOut);

	return !timedOut;

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
