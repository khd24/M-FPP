/*
 * functions_2.h
 *
 *  Created on: Oct 26, 2012
 *      Author: ketandat
 */

#ifndef FUNCTIONS_2_H_
#define FUNCTIONS_2_H_

#include <float.h>
#include <vector>
#include "functions_u.h"
#include "structures_1.h"
#include "structures_2.h"

using namespace std;

void createNetwork(Node*, vector<float>, vector<float>);
void createCells(Cell* , Node* , vector<float> , vector<float>);
void getBarrierCellOverlap(Node*, Cell*, vector<Barrier>, vector<float>, vector<float>);
void getIONodeOverlap(Node*, vector<IOPoint> &, vector<float> , vector<float>);

#endif /* FUNCTIONS_2_H_ */
