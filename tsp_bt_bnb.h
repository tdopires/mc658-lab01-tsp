/*******************************************************************************
 * MC658 - Projeto e Análise de Algoritmos III - 2s2016
 * Prof: Flavio Keidi Miyazawa
 * PED: Mauro Henrique Mulati
 * Usa ideias e código de Rafael Arakaki e Flávio Keidi Miyazawa 
 ******************************************************************************/

/*******************************************************************************
 * ATENÇÃO: NÃO ALTERE ESTE ARQUIVO
 ******************************************************************************/

#ifndef TSP_BT_BNB_H
#define TSP_BT_BNB_H

#include "tsp.h"

bool bt(TSP_Data &tsp, int maxTime);
bool bnb(TSP_Data &tsp, int maxTime);
bool greedy(TSP_Data &tsp, int maxTime);

#endif
