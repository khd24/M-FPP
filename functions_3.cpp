/*
 * functions_3.cpp
 *
 *  Created on: Jun 13, 2020
 *      Author: ketandat
 */

#include "functions_3.h"

void bottomRightFill_dom(Node* _node, Cell* _celllist, vector<float>& _nfcoord, float _nfwidth, float _nfheight)
{
	float width = 0;
	float height = 0;
	queue<Cell*> q;

	Cell* cell = _node->bottom_right_cell;

	while (cell && float_lt(width, _nfwidth))
	{
		if ((cell->barrier_Index) < 0)
		{
			width += (cell->right_gridline) - (cell->left_gridline);
			q.push(cell);
			cell = cell->right_cell;
		}
		else
			return;
	}

	while (!q.empty())
	{
		height = 0;
		cell = q.front();
		q.pop();
		while (cell && float_lt(height, _nfheight))
		{
			if ((cell->barrier_Index) < 0)
			{
				height += (cell->top_gridline) - (cell->bottom_gridline);
				cell = cell->bottom_cell;
			}
			else
				return;
		}
	}

	if (!float_lt(width, _nfwidth) && !float_lt(height, _nfheight))
	{
		_nfcoord.push_back((_node->v_gridline));
		_nfcoord.push_back((_node->h_gridline));
	}
}

void bottomLeftFill_dom(Node* _node, Cell* _celllist, vector<float>& _nfcoord, float _nfwidth, float _nfheight)
{
	float width = 0;
	float height = 0;
	queue<Cell*> q;

	Cell* cell = _node->bottom_left_cell;

	while (cell && float_lt(width, _nfwidth))
	{
		if ((cell->barrier_Index) < 0)
		{
			width += (cell->right_gridline) - (cell->left_gridline);
			q.push(cell);
			cell = cell->left_cell;
		}
		else
			return;
	}

	while (!q.empty())
	{
		height = 0;
		cell = q.front();
		q.pop();
		while (cell && float_lt(height, _nfheight))
		{
			if ((cell->barrier_Index) < 0)
			{
				height += (cell->top_gridline) - (cell->bottom_gridline);
				cell = cell->bottom_cell;
			}
			else
				return;
		}
	}

	if (!float_lt(width, _nfwidth) && !float_lt(height, _nfheight))
	{
		_nfcoord.push_back((_node->v_gridline) - _nfwidth);
		_nfcoord.push_back((_node->h_gridline));
	}
}

void topLeftFill_dom(Node* _node, Cell* _celllist, vector<float>& _nfcoord, float _nfwidth, float _nfheight)
{
	float width = 0;
	float height = 0;
	queue<Cell*> q;

	Cell* cell = _node->top_left_cell;

	while (cell && float_lt(width, _nfwidth))
	{
		if ((cell->barrier_Index) < 0)
		{
			width += (cell->right_gridline) - (cell->left_gridline);
			q.push(cell);
			cell = cell->left_cell;
		}
		else
			return;
	}

	while (!q.empty())
	{
		height = 0;
		cell = q.front();
		q.pop();
		while (cell && float_lt(height, _nfheight))
		{
			if ((cell->barrier_Index) < 0)
			{
				height += (cell->top_gridline) - (cell->bottom_gridline);
				cell = cell->top_cell;
			}
			else
				return;
		}
	}

	if (!float_lt(width, _nfwidth) && !float_lt(height, _nfheight))
	{
		_nfcoord.push_back((_node->v_gridline) - _nfwidth);
		_nfcoord.push_back((_node->h_gridline) + _nfheight);
	}
}

void topRightFill_dom(Node* _node, Cell* _celllist, vector<float>& _nfcoord, float _nfwidth, float _nfheight)
{
	float width = 0;
	float height = 0;
	queue<Cell*> q;

	Cell* cell = _node->top_right_cell;

	while (cell && float_lt(width, _nfwidth))
	{
		if ((cell->barrier_Index) < 0)
		{
			width += (cell->right_gridline) - (cell->left_gridline);
			q.push(cell);
			cell = cell->right_cell;
		}
		else
			return;
	}

	while (!q.empty())
	{
		height = 0;
		cell = q.front();
		q.pop();
		while (cell && float_lt(height, _nfheight))
		{
			if ((cell->barrier_Index) < 0)
			{
				height += (cell->top_gridline) - (cell->bottom_gridline);
				cell = cell->top_cell;
			}
			else
				return;
		}
	}

	if (!float_lt(width, _nfwidth) && !float_lt(height, _nfheight))
	{
		_nfcoord.push_back((_node->v_gridline));
		_nfcoord.push_back((_node->h_gridline) + _nfheight);
	}
}

// This procedure finds a heuristic placement which assumes corner optimality of the head NF.
void generateCandidates(ObjectiveFunction *optimal, int* permutation, int depth, int nfcount, bool &inf_flag, float *candidate, vector<Barrier> &blist, vector<Barrier> &nflist, vector<IOPoint> &iolist, float* efef, float* efnf, float *nfnf, float lwidth, float lheight) {

	if (depth == nfcount) {

		(*optimal).candidatecount++;
		evaluateCandidates(optimal, permutation, efef, efnf, nfnf, candidate, blist, nflist, iolist, lwidth, lheight, false);

	}


	else {

		int nfid = permutation[depth];
		float nf_width = nflist.at(nfid).getWidth();
		float nf_height = nflist.at(nfid).getHeight();

		vector<float> nfcoord;

		getSingleNFCoordinates(nfcoord, blist, iolist, nf_width, nf_height, lwidth, lheight);

		int candidatecount = nfcoord.size() / 2;

		for (int i = 0; i < candidatecount; i++) {

			*(candidate + 2 * nfid) = nfcoord.at(2 * i);
			*(candidate + 2 * nfid + 1) = nfcoord.at(2 * i + 1);

			Barrier newb;
			newb.setBarrier(*(candidate + 2 * nfid), *(candidate + 2 * nfid + 1), nf_width, nf_height, 1);
			blist.push_back(newb);

			generateCandidates(optimal, permutation, depth + 1, nfcount, inf_flag, candidate, blist, nflist, iolist, efef, efnf, nfnf, lwidth, lheight);

			blist.pop_back();

		}
	}
}

void getSingleNFCoordinates(vector<float> &nfcoord, vector<Barrier> &blist, vector<IOPoint> &iolist, float nf_width, float nf_height, float lwidth, float lheight) {

	vector<float> vgridlines, hgridlines;
	getGridlinecoord(blist, iolist, vgridlines, hgridlines, lwidth, lheight);

	int hgridline_count = hgridlines.size();
	int vgridline_count = vgridlines.size();

	int nodecount = (hgridline_count)*(vgridline_count);

	int ycell_count = hgridline_count - 1;
	int xcell_count = vgridline_count - 1;

	Node *network = new Node[(hgridline_count)*(vgridline_count)]; //additional space for NF nodes
	Cell *cells = new Cell[(xcell_count)*(ycell_count)]; //additional space for NF cells

	createNetwork(network, hgridlines, vgridlines);
	createCells(cells, network, hgridlines, vgridlines);
	getBarrierCellOverlap(network, cells, blist, hgridlines, vgridlines);

	for (int i = 0; i < nodecount; i++)
	{
		bottomRightFill_dom(network + i, cells, nfcoord, nf_width, nf_height);
		topRightFill_dom(network + i, cells, nfcoord, nf_width, nf_height);
		bottomLeftFill_dom(network + i, cells, nfcoord, nf_width, nf_height);
		topLeftFill_dom(network + i, cells, nfcoord, nf_width, nf_height);
	}

	delete[] network;
	delete[] cells;
}


// This procedure finds optimal placement of the NFs by finding additional candidates for head NF due to cyclic interference.
void generateCandidates3(ObjectiveFunction *optimal, int* permutation, int depth, int nfcount, bool &inf_flag, float *candidate, vector<Barrier> &blist, vector<Barrier> &nflist, vector<IOPoint> &iolist, float* efef, float* efnf, float *nfnf, float lwidth, float lheight, vector<set<NFCandidate, CompareByPosition>> &masterlist, vector<set<float, FloatCompare>> &original_qx, vector<set<float, FloatCompare>> &original_qy) {

	if (depth == nfcount) {

		int efcount = blist.size() - nfcount;

		int root_nfid = permutation[0];
		float root_nfwidth = nflist.at(root_nfid).getWidth();
		float root_nfheight = nflist.at(root_nfid).getHeight();

		vector<NFCandidate> root_cndlist;

		NFCandidate cnd;
		cnd.x = *(candidate + 2 * root_nfid);
		cnd.y = *(candidate + 2 * root_nfid + 1);

		//		if(inf_flag)
		root_cndlist.push_back(cnd);

		//		if(isBoundaryPoint(root_nfid, candidate,  nfcount, original_qx, original_qy))
		addNewCandidates(root_nfid, root_cndlist, masterlist, blist, nflist, iolist, lwidth, lheight);

		for (int i = 0; i < root_cndlist.size(); i++) {

			NFCandidate nf1cnd = root_cndlist.at(i);
			*(candidate + 2 * root_nfid) = nf1cnd.x;
			*(candidate + 2 * root_nfid + 1) = nf1cnd.y;

			blist.at(efcount).setX(nf1cnd.x);
			blist.at(efcount).setY(nf1cnd.y);

			evaluateCandidates(optimal, permutation, efef, efnf, nfnf, candidate, blist, nflist, iolist, lwidth, lheight, false);

			(*optimal).candidatecount++;

		}

		*(candidate + 2 * root_nfid) = cnd.x;
		*(candidate + 2 * root_nfid + 1) = cnd.y;
		blist.at(efcount).setX(cnd.x);
		blist.at(efcount).setY(cnd.y);
	}

	else {

		int nfid = permutation[depth];
		float nf_width = nflist.at(nfid).getWidth();
		float nf_height = nflist.at(nfid).getHeight();

		vector<float> nfcoord;

		getSingleNFCoordinates(nfcoord, blist, iolist, nf_width, nf_height, lwidth, lheight);

		int candidatecount = nfcoord.size() / 2;

		for (int i = 0; i < candidatecount; i++) {

			*(candidate + 2 * nfid) = nfcoord.at(2 * i);
			*(candidate + 2 * nfid + 1) = nfcoord.at(2 * i + 1);

			NFCandidate cnd;
			cnd.x = nfcoord.at(2 * i);
			cnd.y = nfcoord.at(2 * i + 1);
			bool flag = inf_flag || (masterlist.at(nfid).find(cnd) == masterlist.at(nfid).end());
			//			masterlist.at(nfid).insert(cnd);

			Barrier newb;
			newb.setBarrier(*(candidate + 2 * nfid), *(candidate + 2 * nfid + 1), nf_width, nf_height, 1);
			blist.push_back(newb);

			generateCandidates3(optimal, permutation, depth + 1, nfcount, flag, candidate, blist, nflist, iolist, efef, efnf, nfnf, lwidth, lheight, masterlist, original_qx, original_qy);

			blist.pop_back();

		}
	}
}

void addNewCandidates(int root_nfid, vector<NFCandidate> &root_cndlist, vector<set<NFCandidate, CompareByPosition>> &masterlist, vector<Barrier> barrierlist, vector<Barrier> nflist, vector<IOPoint> iolist, float lwidth, float lheight) {

	int nfcount = nflist.size();
	float nf_width = nflist.at(root_nfid).getWidth();
	float nf_height = nflist.at(root_nfid).getHeight();

	int efcount = barrierlist.size() - nfcount;

	vector<Barrier> blist(barrierlist);

	blist.erase(blist.begin() + efcount);

	vector<float> nfcoord;

	getSingleNFCoordinates(nfcoord, blist, iolist, nf_width, nf_height, lwidth, lheight);

	int candidatecount = nfcoord.size() / 2;

	for (int i = 0; i < candidatecount; i++) {
		NFCandidate cnd1;
		cnd1.x = nfcoord.at(2 * i);
		cnd1.y = nfcoord.at(2 * i + 1);

		//		if(masterlist.at(root_nfid).find(cnd1) == masterlist.at(root_nfid).end())
		root_cndlist.push_back(cnd1);
	}
}

bool isBoundaryPoint(int root_nfid, float *candidate, int nfcount, vector<set<float, FloatCompare>> &original_qx, vector<set<float, FloatCompare>> &original_qy) {

	for (int m = 0; m < nfcount; m++) {
		if (m != root_nfid) {

			float x = candidate[2 * m];
			float y = candidate[2 * m + 1];

			set<float, FloatCompare>::iterator itx = original_qx.at(m).find(x);
			set<float, FloatCompare>::iterator ity = original_qy.at(m).find(y);

			set<float, FloatCompare>::iterator itx_end = original_qx.at(m).end();
			set<float, FloatCompare>::iterator ity_end = original_qy.at(m).end();

			if ((itx != itx_end) || (ity != ity_end))
				return true;

		}
	}

	return false;
}