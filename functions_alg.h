/*
 * functions_opt.h
 *
 *  Created on: Oct 18, 2012
 *      Author: ketandat
 */

#ifndef FUNCTIONS_ALG_H_
#define FUNCTIONS_ALG_H_

#include <omp.h>
#include <vector>
#include <set>
 //#include <gurobi_c++.h>
#include "functions_lb.h"
#include "functions_3.h"
#include "structures_1.h"
#include "structures_2.h"
#include "CandidateRepository.h"
#include "BBController.h"

using namespace std;

void executeCornerOptimalHeuristic_allPerm(Layout *, vector<int*> , float* , ObjectiveFunction *);
void executeCornerOptimalHeuristic_givenPerm(Layout *, int* , float* times, ObjectiveFunction *);

void executeComprehensiveOptimal(Layout *, vector<int*> , float* , ObjectiveFunction *);

void executeTrueComprehensiveOptimal(Layout* layout, vector<int*> permutations, float* times, ObjectiveFunction* optimal, float maxGap, int maxTime);

void generatePermutation_ALDEP(Layout *, int* );

void lb_QAP(void);
void lb_finite (Layout *, vector<int*> , float* , ObjectiveFunction *, bool, bool);  // Implement with Rectilinear approximation and/or Dijkstra 

void updateFlowMatrices(float *, float *, float *, float *, float*, float *, int, int, vector<int>, vector<int>);
void getUniquePermutations(vector<int*> &, int , int );
long long int integer_keygen(set<int> , int );
bool vector_isequal(vector<long long int> , vector<long long int> );

void executeClass1(Layout *, vector<int*> &, float*, ObjectiveFunction *, int);
void executeClass2(Layout *, vector<int*> &, float* , ObjectiveFunction *, int  );

void executeRepair(Layout *, float* , ObjectiveFunction *, int , int);

void getUniqueCombinations(vector<vector<int>> &, vector<int> &, int , int , int , int, int & );
void verifyObjectiveValue(Layout *, ObjectiveFunction *);

#endif /* FUNCTIONS_ALG_H_ */
