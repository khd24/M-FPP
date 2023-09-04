/*
 * structures_2.h
 *
 *  Created on: Oct 14, 2012
 *      Author: ketandat
 */

#ifndef STRUCTURES_2_H_
#define STRUCTURES_2_H_

#include "functions_u.h"

struct Node;
struct Cell;

struct Cell
{
	int key_X;
	int key_Y;
	int index;
	float left_gridline;
	float bottom_gridline;
	float right_gridline;
	float top_gridline;
	Node* top_left_node;
	Node* bottom_left_node;
	Node* bottom_right_node;
	Node* top_right_node;
	Cell* top_cell;
	Cell* bottom_cell;
	Cell* left_cell;
	Cell* right_cell;
	int barrier_Index;
	float lb;
};

struct Node
{
	int key_X;
	int key_Y;
	int index;
	float v_gridline;
	float h_gridline;
	int top_node;
	int bottom_node;
	int left_node;
	int right_node;
	Cell* bottom_right_cell;
	Cell* top_right_cell;
	Cell* top_left_cell;
	Cell* bottom_left_cell;
	float top_edge;
	float bottom_edge;
	float left_edge;
	float right_edge;
	int io_Index;
	float label;
	float obj_value;
};

struct Candidate
{
	float x;
	float y;
	float value;
	float cell_lb;
	int cell_Index;
};

struct NFCandidate
{
	float x;
	float y;
};

struct ObjectiveFunction
{
	int *perm;
	float *ptr;
	float value;
	float efefobj;
	float efnfobj;
	float nfnfobj;
	int candidatecount;
	int feasibleCandidates;
	float bestbound;
	int rank;
};

struct LBObjectiveFunction {
	float efefobj = numeric_limits<float>::max();
	float efnfobjPlaced = numeric_limits<float>::max();
	float efnfobjToBePlaced = numeric_limits<float>::max();
	float nfnfobjPlaced = numeric_limits<float>::max();
	float nfnfobjToBePlaced = numeric_limits<float>::max();
};

struct NFPriority
{
	int id;
	float interaction;
};

class NodeCompare{
public:
	inline bool operator() (Node* n1, Node* n2) const { return float_gt(n1->label, n2->label); }
};

class CandidateCompare{
public:
	bool operator() (Candidate& c1, Candidate& c2) const;
};

class CompareByPosition {
public:
	bool operator()(const NFCandidate &lhs, const NFCandidate &rhs) const;
};

class CompareByXPosition {
public:
	bool operator()(const NFCandidate &lhs, const NFCandidate &rhs) const;
};

class CompareByYPosition {
public:
	bool operator()(const NFCandidate &lhs, const NFCandidate &rhs) const;
};

class FloatCompare {
public:
	bool operator()(const float &lhs, const float &rhs) const;
};

struct NFPriorityCompare {
	inline bool operator() (const NFPriority n1, const NFPriority n2) const { return float_lt(n1.interaction, n2.interaction); }
};

enum class CoordinateType { X, Y };

enum class NodeStatus { Active, Fathomed, Infeasible, Feasible, Incumbent };

enum class BranchingStrategy { BFS, DFS, BstFS };


#endif /* STRUCTURES_2_H_ */