/*
 * functions_u.h
 *
 *  Created on: Oct 26, 2012
 *      Author: ketandat
 */

#ifndef FUNCTIONS_U_H_
#define FUNCTIONS_U_H_

#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "structures_1.h"

using namespace std;

void mergesort(float* , float*, int, int );
void merge(float*,float*, int ,int, int);
int quicksearch(float* , float , int , int );
int* quickBounds(float*, float, int , int );
bool float_lt(float , float );
bool float_gt(float , float );
bool float_eq(float , float );
float getMaxValue(float*, int , float );
float getMinValue(float* , int , float);
float float_min(float, float, float, float);
string getPath(char*, int);
void printMatrix(float *, int , int , const char* );
vector<float> vectorUnion(vector<float> vec1, vector<float> vec2);
bool barrierOverlap(vector<Barrier> barriers, Barrier b);
bool barrierOverlap(Barrier b1, Barrier b2);

#endif /* FUNCTIONS_U_H_ */