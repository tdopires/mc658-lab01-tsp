/*******************************************************************************
 * MC658 - Projeto e Análise de Algoritmos III - 2s2016
 * Prof: Flavio Keidi Miyazawa
 * PED: Mauro Henrique Mulati
 * Usa ideias e código de Rafael Arakaki e Flávio Keidi Miyazawa 
 ******************************************************************************/

/*******************************************************************************
 * ATENÇÃO: NÃO ALTERE ESTE ARQUIVO
 ******************************************************************************/

#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include "myutils.h"
#include "tsp.h"
#include "tsp_bt_bnb.h"

using namespace lemon;

//------------------------------------------------------------------------------
TSP_Data::TSP_Data(ListGraph &graph,
                   NodeStringMap &nodename,
                   NodePosMap &posicaox,
                   NodePosMap &posicaoy,
                   EdgeValueMap &eweight):
                   g(graph),
                   vname(nodename),
                   ename(graph),
                   vcolor(graph),
                   ecolor(graph),
                   weight(eweight),
                   posx(posicaox),
                   posy(posicaoy),
                   AdjMat(graph,eweight,MY_INF),
                   BestCircuit()  // BestCircuit(countEdges(graph))
{
	NNodes=countNodes(this->g);
	NEdges=countEdges(this->g);
	BestCircuitValue = DBL_MAX;
	// max_perturb2opt_it = 3000;  // default value
}
//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	// Variables to be obtained from parameters
	int    exec    = 0;      // 0: not set, 1: BT, 2: BNB
	int    maxTime = 0;      // 0: not set
   bool verbose   = false;
   string inputFile_name;   // Input graph file
   string outputFile_name;  // Output sol file
	
	// Reading program arguments
   for(int i = 1; i < argc; i++){
      const string arg(argv[i]);
      string next;
      if((i+1) < argc){
         next = string(argv[i+1]);
		}
      else{
         next = string("");
		}
		
      if( exec != 0 && (arg.find("-k") == 0 || arg.find("-a") == 0 ) ){
         cout << "Erro ao ler parametro \"" << arg << "\", somente pode haver um parametro de modo de execucao." << endl;
         showUsage();
         exit(1);
      }
      else if( arg.find("-k") == 0 ){
         exec = 1;
      }
      else if( arg.find("-a") == 0 ){
         exec = 2;
      }
      else if( arg.find("-v") == 0 ){
        verbose = true;
      }
      else if( arg.find("-t") == 0 && next.size() > 0){
         maxTime = atoi(next.c_str()); i++; continue;
      }
      else if( arg.find("-i") == 0 && next.size() > 0){
         inputFile_name = next; i++; continue;
      }
      else if( arg.find("-o") == 0 && next.size() > 0){
         outputFile_name = next; i++; continue;
      }
      else{
         cout << "Parametro invalido: \"" << arg << "\"" << " (ou faltando argumento)" << endl;
         showUsage();
         exit(1);
      }
   }

   // Required parameters
   if( exec == 0 ){
      cout << "Nenhum modo de execucao selecionado dentre: -k ou -a" << endl;
      showUsage(); 
		exit(1);
   }
   
   if( inputFile_name.size() < 1 ){
      cout << ((inputFile_name.size() < 1)? "nome do arq de grafo, ":"") 
			  << endl;
      showUsage(); 
		exit(1);
   }

   if( outputFile_name.size() < 1 ){
      cout << ((outputFile_name.size() < 1)? "nome do arq de saida.":"") 
			  << endl;
      showUsage(); 
		exit(1);
   }
   
   if( maxTime == 0 ){
        maxTime = 30;   // Default of 30s
   }

	// int seed=1;     // mhmulati
	// srand48(seed);  // mhmulati
	
	// Variables that represent the weighted undirected graph of the asymmetric tsp
	ListGraph     g;
	EdgeValueMap  weight(g);
	NodeStringMap vname(g);
	NodePosMap    posx(g), 
	              posy(g);

	// Read the graph from the imput file
	if ( ! ReadListGraph(inputFile_name, g, vname, weight, posx, posy) ) {
	  cout << "Erro na leitura do arquivo de entrada " << inputFile_name << endl;
	  exit(1);
	}
	
	// Init the graph data structure
	TSP_Data tsp(g, vname, posx, posy, weight);
	
	// cerr << tspInstanceAsString(tsp) << endl;
	
   double elapsed_time = DBL_MAX;
   clock_t before = clock();
   bool foundOptimalSolution = false;
	
	switch(exec){
		case 1:{
			foundOptimalSolution = bt(tsp, maxTime);
			break;
		}
		case 2:{
			foundOptimalSolution = bnb(tsp, maxTime);
			break;
		}
	}
	
   clock_t after = clock();
   elapsed_time = (double) (after-before) / CLOCKS_PER_SEC;

	// Verificar se é de fato uma solução para a instância
	bool isSol = checkSol(tsp);
	
	// Imprimir a solucao em arquivo
   // writeOutputFile(outputFile_name, inputFile_name, tsp, elapsed_time, maxTime, exec, foundOptimalSolution);
	ofstream myfile;
   myfile.open(outputFile_name);
	myfile << resultAsString(tsp, exec, inputFile_name, elapsed_time, maxTime, foundOptimalSolution, isSol) << flush;
	myfile.close();
	
  // Imprime a solução na tela
	cout << endl;
   if(verbose){
		cout << resultAsString(tsp, exec, inputFile_name, elapsed_time, maxTime, foundOptimalSolution, isSol) << flush;
		if( tsp.BestCircuitValue < DBL_MAX && isSol){
			ViewTspCircuit(tsp);
		}
		cout << endl;
	}
	
	return !foundOptimalSolution;
}
//------------------------------------------------------------------------------
void ViewTspCircuit(TSP_Data &tsp)
{
	ListGraph h;
	NodeStringMap h_vname(h);  // node names
	NodeNodeMap g2h(tsp.g);  // maps a node of g to a node of h
	NodePosMap h_posx(h);
	NodePosMap h_posy(h);
	NodeColorMap vcolor(h);   // color of the vertices
	EdgeColorMap acolor(h);  // color of edges
	EdgeStringMap aname(h);  // name of edges
	
	for(NodeIt v(tsp.g); v != INVALID; ++v){
		Node hv;
		hv = h.addNode();
		g2h[v] = hv;
		h_posx[hv] = tsp.posx[v];
		h_posy[hv] = tsp.posy[v];
		h_vname[hv] = tsp.vname[v];
		vcolor[hv] = BLUE;
	}
	
	for(int i = 0; i < tsp.NNodes; i++){
		Node u, v;
		Edge a;
		u = tsp.BestCircuit[i]; 
		v = tsp.BestCircuit[(i+1) % tsp.NNodes]; 
		a = h.addEdge(g2h[u] , g2h[v]);
		aname[a] = "";
		acolor[a] = BLUE;
	}
	
	ViewListGraph(h,h_vname,aname,h_posx,h_posy,vcolor,acolor,"TSP Circuit with cost "+DoubleToString(tsp.BestCircuitValue));
}
//------------------------------------------------------------------------------
void showUsage()
// Usage information
{
	cout << "Usage:"
	     << "./tsp <modo_operacao> (um dentre: -k backtracking, -a branch_and_bound) -t <tempo_max_em_segundos> {-v: mostra solução na tela}"
	     << "-i <nome_arquivo_entrada> -o <nome_arquivo_saida>"
	     << endl;
}
//------------------------------------------------------------------------------
/* void writeOutputFile(string outputfile, 
                     string graphname, 
                     TSP_Data &tsp, 
                     double elapsed_time, 
                     int max_time, 
                     int exec, 
                     bool opt)
{
   ofstream myfile;
   myfile.open(outputfile);
	if(exec == 1){
		myfile << "BACKTRACKING" << '\t'<< opt << '\t' << graphname << '\t';
		if( tsp.BestCircuitValue < DBL_MAX ){
			myfile << tsp.BestCircuitValue << '\t' << elapsed_time << '\t' << max_time << endl;
		}
		else{
			myfile << "NAO ENCONTROU SOLUCAO VIAVEL" << endl;
		}
	}
	else{
		myfile << "BRANCH AND BOUND" << '\t'<< opt  << '\t' << graphname << '\t';
		if( tsp.BestCircuitValue < DBL_MAX ){
			myfile << tsp.BestCircuitValue << '\t' << elapsed_time << '\t' << max_time << endl;
		}
		else{
			myfile << "NAO ENCONTROU SOLUCAO VIAVEL" << endl;
		}
	}
   
   myfile.close();
} */
//------------------------------------------------------------------------------
bool checkSol(TSP_Data &tsp)
{
	// cerr << tspInstanceAsString(tsp) << endl;
	
	// cerr << tspSolutionAsString(tsp) << endl;
	
	// Number of nodes
	if(tsp.NNodes != (int)tsp.BestCircuit.size()){
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
	
	double cost = 0.0;
	int nedges = 0;
	
	Node ni = INVALID;
	Node nj = INVALID;
	Node ad = INVALID;
	
	// Lets check each edge inferred by the sequence of nodes.
	for(auto i = tsp.BestCircuit.begin(); i != tsp.BestCircuit.end(); ++i){
		ni = *i;
		// cerr << endl;
		// cerr << "ni: " << tsp.vname[ni] << endl;
		// cerr << "  y[ni]: " << y[ni] << endl;
		
		// If ni is the last element, lets examine the edge (last,first).
		if( next(i) == tsp.BestCircuit.end()){
			nj = *(tsp.BestCircuit.begin());
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
					cost += tsp.weight[e];
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
	
	// cerr << "cost: " << cost << endl;
	// cerr << "tsp.BestCircuitValue: " << tsp.BestCircuitValue << endl;	
	
	// Cost
	if(cost < tsp.BestCircuitValue - MY_EPS || cost > tsp.BestCircuitValue + MY_EPS){
		return false;
	}
	
	return true;
}
//------------------------------------------------------------------------------
string tspInstanceAsString(TSP_Data &tsp)
{
	stringstream tspInstanceAsString;
	tspInstanceAsString << "TSP Instance" << endl;
	tspInstanceAsString << "n: " << tsp.NNodes << "  m: " << tsp.NEdges << endl;
	tspInstanceAsString << "Nodes:";
	for(NodeIt v(tsp.g); v != INVALID; ++v){
		tspInstanceAsString << " " << tsp.vname[v];
	}
	tspInstanceAsString << endl;
	 
	tspInstanceAsString << "Edges:";
	for(EdgeIt e(tsp.g); e != INVALID; ++e){
		tspInstanceAsString << "(" << tsp.vname[tsp.g.u(e)] << "," << tsp.vname[tsp.g.v(e)] << "):" << tsp.weight[e] << " ";
	}
	tspInstanceAsString << endl;
	
	return tspInstanceAsString.str();
}
//------------------------------------------------------------------------------
string tspSolutionAsString(TSP_Data &tsp)
{
	stringstream tspSolutionAsString;
	tspSolutionAsString << "BestCircuitValue  : " << tsp.BestCircuitValue << endl;
	tspSolutionAsString << "BestCircuit       :";
	for(auto i = tsp.BestCircuit.begin(); i != tsp.BestCircuit.end(); ++i){
		tspSolutionAsString << " " << tsp.vname[*i];
	}
	tspSolutionAsString << endl;
	return tspSolutionAsString.str();
}
//------------------------------------------------------------------------------
string resultAsString(TSP_Data &tsp, int exec, string inputFile_name, double elapsed_time, int maxTime, bool foundOptimalSolution, bool isSol)
{
	stringstream resultAsString;

	resultAsString << "Algorithm         : " << (exec == 1?"BACKTRACKING":"BRANCH AND BOUND") << endl;
	resultAsString << "Instance          : " << inputFile_name << endl;
	resultAsString << "ElapsedTime       : " << elapsed_time << "s" << endl;
	resultAsString << "TimeLimit         : " << maxTime << "s" <<  endl;
	if( tsp.BestCircuitValue < DBL_MAX && !isSol){
		resultAsString << "Solution returned as feasible by the student's program, but it has not passed the checking of feasibility." << endl;
	}
	else if( tsp.BestCircuitValue < DBL_MAX && isSol){
		resultAsString << "Guaranteed Optimal: " <<  (foundOptimalSolution?"Yes":"No") << endl;
		resultAsString << tspSolutionAsString(tsp);
	}
	else{
		resultAsString << "Have not found feasible solution." << endl;
	}

	return resultAsString.str();
}
//------------------------------------------------------------------------------
