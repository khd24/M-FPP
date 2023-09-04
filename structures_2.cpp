/*
 * structures_2.cpp
 *
 *  Created on: Jun 13, 2020
 *      Author: ketandat
 */

#include "structures_2.h"

bool CompareByPosition::operator()(const NFCandidate &lhs, const NFCandidate &rhs) const {
	if (!float_eq(lhs.x, rhs.x))
		return float_lt(lhs.x, rhs.x);
	return float_lt(lhs.y, rhs.y);
}

bool CompareByXPosition::operator()(const NFCandidate &lhs, const NFCandidate &rhs) const {
	if (!float_eq(lhs.x, rhs.x))
		return float_gt(lhs.x, rhs.x);
	return float_gt(lhs.y, rhs.y);
}

bool CompareByYPosition::operator()(const NFCandidate &lhs, const NFCandidate &rhs) const {
	if (!float_eq(lhs.y, rhs.y))
		return float_gt(lhs.y, rhs.y);
	return float_gt(lhs.x, rhs.x);
}

bool FloatCompare::operator()(const float &lhs, const float &rhs) const {return float_lt(lhs, rhs);}