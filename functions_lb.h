/*
 * functions_lb.h
 *
 *  Created on: Jun 13, 2020
 *      Author: ketandat, aanandan
 */

/*
deal with the data from problemsset
for # of NF= N

creating:
1. cost function of placing NF at any given node in the layout with regards to EF (N by N^2 matrix)
2. distance function of one node in the layout to all other node in the layout, used to extract cost function of placing NF at any given node in the layout with regards to NF (N^2 by N^2 matrix)
3. partitioning the facility with the layout location, return the index of the location (1-D array with size N)


*/

#ifndef FUNCTIONS_LB_H_
#define FUNCTIONS_LB_H_

#include <vector>
#include <omp.h>
#include "functions_4.h"
#include "structures_1.h"
#include "structures_2.h"
#include "functions_u.h"
#include "functions_3.h"

using namespace std;

void updateFlowMatrices_bnb(float* efef2, float* efnf2, float* nfnf2, float* efef, float* efnf, float* nfnf, int efcount, int nfcount, vector<int> placed, vector<int> to_be_placed);
float computeFiniteLB(vector<Barrier> blist, vector<IOPoint> ilist, vector<Barrier> nlist, float* efefflow, float* efnfflow, float* nfnfflow, float lwidth, float lheight, float* times, int lbclass);
void allPairSP(float *, Node *, int , int );
void constructXCosts(float * , float *, float* , vector<IOPoint> , int , int , int , int );
void getEFInteraction(float &, vector<IOPoint> , float *, float* , int );
//void solveQSAP(float &, float *, float *, float *, int , int );

#endif /* FUNCITONS_LB_H_ */