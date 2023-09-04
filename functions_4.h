/*
 * functions_4.h
 *
 *  Created on: Oct 26, 2012
 *      Author: ketandat
 */

#ifndef FUNCTIONS_5_H_
#define FUNCTIONS_5_H_

#include <float.h>
#include <string.h>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <queue>
#include "structures_1.h"
#include "structures_2.h"
#include "functions_1.h"
#include "functions_2.h"
#include "functions_u.h"

using namespace std;

void evaluateCandidates(ObjectiveFunction *optimum, int *cur_perm, float *efefflow, float *efnfflow, float *nfnfflow, float* coord, vector<Barrier> blist, vector<Barrier> nflist, vector<IOPoint> iolist, float lwidth, float lheight, bool rectilinear_approx);
void evaluateCandidates2(ObjectiveFunction *optimum, int *cur_perm, float *efefflow, float *efnfflow, float *nfnfflow, float* coord, vector<Barrier> blist, vector<Barrier> nflist, vector<IOPoint> iolist, float lwidth, float lheight, bool rectilinear_approx);
void evaluateCandidates3(LBObjectiveFunction* optimum, float* efefflow, float* efnfflow, float* nfnfflow, vector<Barrier> blist, vector<Barrier> nflist, vector<IOPoint> iolist, float lwidth, float lheight, bool rectilinear_approx, int totalNfCount);
void executeDijkstra(Node* ,Node*, float* , int , int );
void copyNode(Node*, Node*);

void dijkstraDistances(float *, vector<IOPoint> io, Node *network,  int _hgridline_count, int _vgridline_count);
void rectilinearDistances(float *, vector<IOPoint> );

void calculateEFEFObj (float &efefobj, float *efefflow, vector<Barrier> blist,vector<IOPoint> iolist, float lwidth, float lheight, bool rectilinear_approx);

struct myclass;

#endif /* FUNCTIONS_5_H_ */
