/*******************************************************************************
 * MC658 - Projeto e Análise de Algoritmos III - 2s2016
 * Prof: Flavio Keidi Miyazawa
 * PED: Mauro Henrique Mulati
 * Usa ideias e código de Rafael Arakaki e Flávio Keidi Miyazawa 
 ******************************************************************************/

#ifndef TSP_H
#define TSP_H

/* ATENÇÃO: NÃO ALTERE ESTE ARQUIVO */

/* 
 * This is the type used to obtain the pointer to the problem data. This pointer
 * is stored in the branch and cut tree. And when we define separation routines,
 * we can recover the pointer and access the problem data again.
 */
class TSP_Data {
public:
	TSP_Data(ListGraph &graph,
	         NodeStringMap &nodename,
	         NodePosMap &posicaox,
	         NodePosMap &posy,
	         EdgeValueMap &eweight);
	ListGraph &g;
	int NNodes,NEdges;
	// int max_perturb2opt_it; // maximum number of iterations for heuristic Perturb2OPT
	NodeStringMap &vname;
	EdgeStringMap ename;
	NodeColorMap vcolor;
	EdgeColorMap ecolor;
	EdgeValueMap &weight;
	NodePosMap &posx;
	NodePosMap &posy;
	AdjacencyMatrix AdjMat; // adjacency matrix
	vector<Node> BestCircuit; // vector containing the best circuit found
	double BestCircuitValue;
};

// Usage information
void showUsage();
void ViewTspCircuit(TSP_Data &tsp);

bool checkSol(TSP_Data &tsp);
string tspInstanceAsString(TSP_Data &tsp);
string tspSolutionAsString(TSP_Data &tsp);
string resultAsString(TSP_Data &tsp, int exec, string inputFile_name, double elapsed_time, int maxTime, bool foundOptimalSolution, bool isSol);

// void writeOutputFile(string outputFile_name, string inputFile_name, TSP_Data &tsp, double elapsed_time, int maxTime, int exec, bool opt);

// Global variables reference
// extern bool verbose; //0=not set, 1=print solution after optimization

#endif
