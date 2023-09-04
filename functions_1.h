/*
 * functions_1.h
 *
 *  Created on: Oct 26, 2012
 *      Author: ketandat
 */

#ifndef FUNCTIONS_1_H_
#define FUNCTIONS_1_H_

#include "functions_u.h"
#include "structures_1.h"
#include "structures_2.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

Layout* readLayoutFromFile(std::string, int& );
void writeLayoutToFile(const string _layoutfile, int efcount, int nfcount, Barrier* blist, Barrier* nflist, IOPoint* iolist);

void getGridlinecoord( std::vector<Barrier> , std::vector<IOPoint> , std::vector<float> &, std::vector<float> &, float , float );
void getVGridlinecoord(std::vector<Barrier>, std::vector<IOPoint>, std::vector<float> &, float, float);
void getHGridlinecoord(std::vector<Barrier>, std::vector<IOPoint>, std::vector<float> &, float, float);
void createFlowMatrix(float*, Node* , Layout*, float*, float*, int , int );
void swap(int *, int *);
void permute(vector<int*>&, int *, int, int);
void removeDuplicates(vector<float>&);

#endif /* FUNCTIONS_1_H_ */


