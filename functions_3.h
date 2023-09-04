/*
 * functions_3.h
 *
 *  Created on: Oct 26, 2012
 *      Author: ketandat
 */

#ifndef FUNCTIONS_3_H_
#define FUNCTIONS_3_H_

#include "functions_u.h"
#include "functions_1.h"
#include "functions_2.h"
#include "functions_4.h"
#include <set>
#include <queue>

using namespace std;

void bottomRightFill_dom(Node*, Cell* , vector<float>& , float , float);
void topRightFill_dom(Node*, Cell* , vector<float>& , float , float);
void bottomLeftFill_dom(Node*, Cell* , vector<float>& , float , float);
void topLeftFill_dom(Node*, Cell*, vector<float>& , float , float);
void generateCandidates(ObjectiveFunction *, int* , int , int , bool &, float *, vector<Barrier> &, vector<Barrier> &, vector<IOPoint> &, float* , float* , float *, float , float );
void getSingleNFCoordinates(vector<float> &, vector<Barrier> &, vector<IOPoint> &, float , float , float , float );

void generateCandidates3(ObjectiveFunction *, int* , int , int , bool &, float *, vector<Barrier> &, vector<Barrier> &, vector<IOPoint> &, float* , float* , float *, float , float , vector<set<NFCandidate, CompareByPosition>> &, vector<set<float, FloatCompare>> &, vector<set<float, FloatCompare>> &);
void addNewCandidates(int , vector<NFCandidate> &, vector<set<NFCandidate, CompareByPosition>> &, vector<Barrier> , vector<Barrier> , vector<IOPoint> , float , float ); 
bool isBoundaryPoint(int, float *, int , vector<set<float, FloatCompare>> &, vector<set<float, FloatCompare>> &);

#endif /* FUNCTIONS_3_H_ */
