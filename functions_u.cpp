/*
 * functions_u.cpp
 *
 *  Created on: Jun 13, 2020
 *      Author: ketandat
 */

#include "functions_u.h"

float getMaxValue(float* _array, int _size, float _lim) {
	float max = _lim;
	for (int i = 0; i < _size; i++) {
		if (_array[i] > max)
			max = _array[i];
	}
	return max;
}

float getMinValue(float* _array, int _size, float _lim) {
	float min = _lim;
	for (int i = 0; i < _size; i++) {
		if (_array[i] < min)
			min = _array[i];
	}
	return min;
}

float float_min(float f1, float f2, float f3, float f4)
{
	float val;

	if (!float_gt(f1, f2))
		val = f1;
	else
		val = f2;

	if (float_gt(val, f3))
		val = f3;

	if (float_gt(val, f4))
		val = f4;

	return val;
}

void mergesort(float* _array, float* _temp, int min, int max) {
	int mid;
	if (min < max) {
		mid = (min + max) / 2;
		mergesort(_array, _temp, min, mid);
		mergesort(_array, _temp, mid + 1, max);
		merge(_array, _temp, min, mid, max);
	}
}

void merge(float* _array, float* _temp, int _min, int _mid, int _max) {
	int i, j, k;

	i = _min;
	j = _mid + 1;
	k = _min;
	while (i <= _mid && j <= _max) {
		if (float_lt(_array[i], _array[j]) || float_eq(_array[i], _array[j])) {
			_temp[k] = _array[i];
			i++;
		}
		else {
			_temp[k] = _array[j];
			j++;
		}
		k++;
	}
	if (j > _max) {
		for (int z = i; z <= _mid; z++) {
			_temp[k] = _array[z];
			k++;

		}
	}
	else {
		for (int z = j; z <= _max; z++) {
			_temp[k] = _array[z];
			k++;
		}

	}
	for (int z = _min; z <= _max; z++) _array[z] = _temp[z];
}

int quicksearch(float* inarray, float value, int _min, int _max) {

	int min, mid, max;
	min = _min;
	max = _max;

	while (max - min > 1)
	{
		mid = min + ((max - min) / 2);
		if (!float_lt(value, inarray[mid]))
			min = mid;
		else if (!float_gt(value, inarray[mid]))
			max = mid;
	}

	if (float_eq(value, inarray[min]))
		return min;
	else if (float_eq(value, inarray[max]))
		return max;

	return -1;
}

int* quickBounds(float* inarray, float value, int _min, int _max)
{
	int min, mid, max;
	int* bounds = 0;
	min = _min;
	max = _max;

	while (max - min > 1)
	{
		mid = min + ((max - min) / 2);
		if (!float_lt(value, inarray[mid]))
			min = mid;
		else if (!float_gt(value, inarray[mid]))
			max = mid;
	}

	if (float_eq(value, inarray[min]))
	{
		bounds = new int[2];
		bounds[0] = min;
		bounds[1] = min;
	}

	else if (float_eq(value, inarray[max]))
	{
		bounds = new int[2];
		bounds[0] = max;
		bounds[1] = max;
	}
	else if (float_gt(value, inarray[min]) && float_lt(value, inarray[max]))
	{
		bounds = new int[2];
		bounds[0] = min;
		bounds[1] = max;
	}

	return bounds;
}

bool float_lt(float a, float b) {

	return (a < (b - 0.0001));
}

bool float_gt(float a, float b) 
{
	return (a > (b + 0.0001));
}

bool float_eq(float a, float b) 
{
	return (abs(a - b) < 0.0001);
}

string getPath(char* _name, int _rank)
{
	stringstream s;
	s << _name << _rank << ".txt";
	return s.str();
}

void printMatrix(float *matrix, int rowsize, int colsize, const char* name) {

	cout << name << endl;
	for (int i = 0; i < rowsize; i++) {
		for (int j = 0; j < colsize; j++) {
			cout << matrix[i * colsize + j] << "\t";
		}
		cout << endl;
	}
	cout << endl;
}

vector<float> vectorUnion(vector<float> vec1, vector<float> vec2) {

	vector<float> retVector;

	retVector.insert(retVector.end(), vec1.begin(), vec1.end());

	retVector.insert(retVector.end(), vec2.begin(), vec2.end());

	sort(retVector.begin(), retVector.end());

	vector<float>::iterator it = unique(retVector.begin(), retVector.end(), float_eq);
	retVector.resize(distance(retVector.begin(), it));

	return retVector;
}

bool barrierOverlap(vector<Barrier> barriers, Barrier b) {	
	
	for (auto it = barriers.begin(); it != barriers.end(); ++it) {
		
		if (barrierOverlap(*it, b)) {

			return true;
		}
	}

	return false;
}

bool barrierOverlap(Barrier b1, Barrier b2) {

	if (!float_lt(b1.getMinX(), b2.getMaxX()) || !float_lt(b2.getMinX(), b1.getMaxX()))
		return false;

	if (!float_lt(b1.getMinY(), b2.getMaxY()) || !float_lt(b2.getMinY(), b1.getMaxY()))
		return false;

	return true;
}