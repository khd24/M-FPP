/*
 * functions_4.cpp
 *
 *  Created on: Jun 13, 2020
 *      Author: ketandat
 */

#include "functions_4.h"

struct myclass {
	bool operator() (Node *a, Node *b) { return float_gt(a->label, b->label); }
} myobject;

void evaluateCandidates2(ObjectiveFunction *optimum, int *cur_perm, float *efefflow, float *efnfflow, float *nfnfflow, float* coord, vector<Barrier> blist, vector<Barrier> nflist, vector<IOPoint> iolist, float lwidth, float lheight, bool rectilinear_approx) {

	int efiocount = iolist.size();
	int nfiocount = nflist.size();

	vector<Barrier> barrier(blist);
	vector<IOPoint> io(iolist);

	for (int j = 0; j < nflist.size(); j++) {
		Barrier nf;
		nf.setBarrier(coord[2 * j], coord[2 * j + 1], nflist.at(j).getWidth(), nflist.at(j).getHeight(), 1.0);
		barrier.push_back(nf);

		IOPoint nfio;
		nfio.x = coord[2 * j];
		nfio.y = coord[2 * j + 1];
		nfio.isNF = true;
		io.push_back(nfio);
	}

	vector<float> vgridlines, hgridlines;

	getGridlinecoord(barrier, io, vgridlines, hgridlines, lwidth, lheight);

	int _hgridline_count = hgridlines.size();
	int _vgridline_count = vgridlines.size();
	int ycell_count = _hgridline_count - 1;
	int xcell_count = _vgridline_count - 1;

	Node *network = new Node[_hgridline_count * _vgridline_count];
	Cell *cells = new Cell[(xcell_count)*(ycell_count)]; //additional space for NF cells

	createNetwork(network, hgridlines, vgridlines);
	createCells(cells, network, hgridlines, vgridlines);
	getBarrierCellOverlap(network, cells, barrier, hgridlines, vgridlines);
	getIONodeOverlap(network, io, hgridlines, vgridlines);

	float *distances = new float[io.size() * io.size()];

	if (rectilinear_approx)
		rectilinearDistances(distances, io);
	else
		dijkstraDistances(distances, io, network, _hgridline_count, _vgridline_count);

	float obj = 0;
	float efnf = 0;
	float efef = 0;

	for (int ioid1 = 0; ioid1 < efiocount; ioid1++) {
		for (int ioid2 = ioid1 + 1; ioid2 < efiocount; ioid2++) {
			efef += distances[ioid1 * (efiocount + nfiocount) + ioid2] * efefflow[ioid1 * efiocount + ioid2];
		}

		for (int ioid2 = 0; ioid2 < nfiocount; ioid2++) {
			int tioid2 = ioid2 + efiocount;
			efnf += distances[ioid1 * (efiocount + nfiocount) + tioid2] * efnfflow[ioid1 * nfiocount + ioid2];
		}
	}

	obj += (efef + efnf);

	for (int ioid1 = 0; ioid1 < nfiocount; ioid1++) {
		for (int ioid2 = ioid1 + 1; ioid2 < nfiocount; ioid2++) {
			int tioid1 = ioid1 + efiocount;
			int tioid2 = ioid2 + efiocount;

			obj += distances[tioid1 * (efiocount + nfiocount) + tioid2] * nfnfflow[ioid1 * nfiocount + ioid2];
		}
	}

	if (obj < (*optimum).value) {
		copy(coord, coord + 2 * nflist.size(), (*optimum).ptr);
		copy(cur_perm, cur_perm + nflist.size(), (*optimum).perm);
		(*optimum).value = obj;
		(*optimum).efnfobj = efnf;
		(*optimum).efefobj = efef;
	}


	delete[] distances;
	delete[] network;
	delete[] cells;

}

void evaluateCandidates3(LBObjectiveFunction* optimum, float* efefflow, float* efnfflow, float* nfnfflow, vector<Barrier> blist, vector<Barrier> nflist, vector<IOPoint> iolist, float lwidth, float lheight, bool rectilinear_approx, int totalNfCount) {

	int efiocount = iolist.size();
	int nfiocount = nflist.size();

	vector<Barrier> barrier(blist);
	vector<IOPoint> io(iolist);

	for (auto it = nflist.begin(); it != nflist.end(); ++it) {
		barrier.push_back(*it);
		
		IOPoint nfio;
		nfio.x = it->getMinX();
		nfio.y = it->getMaxY();
		nfio.isNF = true;
		io.push_back(nfio);
	}

	float* distances = new float[io.size() * io.size()];

	if (rectilinear_approx) {
		rectilinearDistances(distances, io);
	}
	
	else {

		vector<float> vgridlines, hgridlines;

		getGridlinecoord(barrier, io, vgridlines, hgridlines, lwidth, lheight);

		int _hgridline_count = hgridlines.size();
		int _vgridline_count = vgridlines.size();
		int ycell_count = _hgridline_count - 1;
		int xcell_count = _vgridline_count - 1;

		Node* network = new Node[_hgridline_count * _vgridline_count];
		Cell* cells = new Cell[(xcell_count) * (ycell_count)]; //additional space for NF cells

		createNetwork(network, hgridlines, vgridlines);
		createCells(cells, network, hgridlines, vgridlines);
		getBarrierCellOverlap(network, cells, barrier, hgridlines, vgridlines);
		getIONodeOverlap(network, io, hgridlines, vgridlines);

		dijkstraDistances(distances, io, network, _hgridline_count, _vgridline_count);

		delete[] network;
		delete[] cells;
	}

	float efef = 0;
	float efnfPlaced = 0;
	float nfnfPlaced = 0;
	float efnfToBePlaced = 0;
	float nfnfToBePlaced = 0;

	for (int ioid1 = 0; ioid1 < efiocount; ioid1++) {
		for (int ioid2 = ioid1 + 1; ioid2 < efiocount; ioid2++) {
			efef += distances[ioid1 * (efiocount + nfiocount) + ioid2] * efefflow[ioid1 * efiocount + ioid2];
		}

		for (int ioid2 = 0; ioid2 < nfiocount - 1; ioid2++) {
			int tioid2 = ioid2 + efiocount;
			int nfid2 = nflist[ioid2].id();
			efnfPlaced += distances[ioid1 * (efiocount + nfiocount) + tioid2] * efnfflow[ioid1 * totalNfCount + nfid2];
		}

		int ioid3 = nfiocount - 1;
		int tioid3 = ioid3 + efiocount;
		int nfid3 = nflist[ioid3].id();
		efnfToBePlaced += distances[ioid1 * (efiocount + nfiocount) + tioid3] * efnfflow[ioid1 * totalNfCount + nfid3];
	}

	for (int ioid1 = 0; ioid1 < nfiocount; ioid1++) {
		int tioid1 = ioid1 + efiocount;
		int nfid1 = nflist[ioid1].id();
		
		for (int ioid2 = ioid1 + 1; ioid2 < nfiocount - 1; ioid2++) {
			int tioid2 = ioid2 + efiocount;
			int nfid2 = nflist[ioid2].id();
			nfnfPlaced += distances[tioid1 * (efiocount + nfiocount) + tioid2] * nfnfflow[nfid1 * totalNfCount + nfid2];
		}

		int ioid3 = nfiocount - 1;
		int tioid3 = ioid3 + efiocount;
		int nfid3 = nflist[ioid3].id();
		nfnfToBePlaced += distances[tioid1 * (efiocount + nfiocount) + tioid3] * nfnfflow[nfid1 * totalNfCount + nfid3];
	}

	if (efef < (*optimum).efefobj) {
		(*optimum).efefobj = efef;
	}

	if (efnfPlaced < (*optimum).efnfobjPlaced) {
		(*optimum).efnfobjPlaced = efnfPlaced;
	}

	if (nfnfPlaced < (*optimum).nfnfobjPlaced) {
		(*optimum).nfnfobjPlaced = nfnfPlaced;
	}

	if (efnfToBePlaced < (*optimum).efnfobjToBePlaced) {
		(*optimum).efnfobjToBePlaced = efnfToBePlaced;
	}

	if (nfnfToBePlaced < (*optimum).nfnfobjToBePlaced) {
		(*optimum).nfnfobjToBePlaced = nfnfToBePlaced;
	}

	delete[] distances;
}

void evaluateCandidates(ObjectiveFunction *optimum, int *cur_perm, float *efefflow, float *efnfflow, float *nfnfflow, float* coord, vector<Barrier> blist, vector<Barrier> nflist, vector<IOPoint> iolist, float lwidth, float lheight, bool rectilinear_approx) {

	int efiocount = iolist.size();
	int nfiocount = nflist.size();

	vector<Barrier> barrier(blist);
	vector<IOPoint> io(iolist);

	for (int j = 0; j < nflist.size(); j++) {
		//		Barrier nf;
		//		nf.setBarrier(coord[2 * j], coord[2 * j + 1], nflist.at(j).getWidth(), nflist.at(j).getHeight(), 1.0);
		//		barrier.push_back(nf);

		IOPoint nfio;
		nfio.x = coord[2 * j];
		nfio.y = coord[2 * j + 1];
		nfio.isNF = true;
		io.push_back(nfio);
	}

	vector<float> vgridlines, hgridlines;

	getGridlinecoord(barrier, io, vgridlines, hgridlines, lwidth, lheight);

	int _hgridline_count = hgridlines.size();
	int _vgridline_count = vgridlines.size();
	int ycell_count = _hgridline_count - 1;
	int xcell_count = _vgridline_count - 1;

	Node *network = new Node[_hgridline_count * _vgridline_count];
	Cell *cells = new Cell[(xcell_count)*(ycell_count)]; //additional space for NF cells

	createNetwork(network, hgridlines, vgridlines);
	createCells(cells, network, hgridlines, vgridlines);
	getBarrierCellOverlap(network, cells, barrier, hgridlines, vgridlines);
	getIONodeOverlap(network, io, hgridlines, vgridlines);

	float *distances = new float[io.size() * io.size()];

	if (rectilinear_approx)
		rectilinearDistances(distances, io);
	else
		dijkstraDistances(distances, io, network, _hgridline_count, _vgridline_count);

	float obj = 0;
	float efnf = 0;
	float efef = 0;

	for (int ioid1 = 0; ioid1 < efiocount; ioid1++) {
		for (int ioid2 = ioid1 + 1; ioid2 < efiocount; ioid2++) {
			efef += distances[ioid1 * (efiocount + nfiocount) + ioid2] * efefflow[ioid1 * efiocount + ioid2];
		}

		for (int ioid2 = 0; ioid2 < nfiocount; ioid2++) {
			int tioid2 = ioid2 + efiocount;
			efnf += distances[ioid1 * (efiocount + nfiocount) + tioid2] * efnfflow[ioid1 * nfiocount + ioid2];
		}
	}

	obj += (efef + efnf);

	for (int ioid1 = 0; ioid1 < nfiocount; ioid1++) {
		for (int ioid2 = ioid1 + 1; ioid2 < nfiocount; ioid2++) {
			int tioid1 = ioid1 + efiocount;
			int tioid2 = ioid2 + efiocount;

			obj += distances[tioid1 * (efiocount + nfiocount) + tioid2] * nfnfflow[ioid1 * nfiocount + ioid2];
		}
	}

	if (obj < (*optimum).value) {
		copy(coord, coord + 2 * nflist.size(), (*optimum).ptr);
		copy(cur_perm, cur_perm + nflist.size(), (*optimum).perm);
		(*optimum).value = obj;
		(*optimum).efnfobj = efnf;
		(*optimum).efefobj = efef;
	}


	delete[] distances;
	delete[] network;
	delete[] cells;

}

void dijkstraDistances(float *distances, vector<IOPoint> io, Node *network, int _hgridline_count, int _vgridline_count) {

	float *shortestpaths = new float[(_hgridline_count)*(_vgridline_count)];

	for (int id = 0; id < io.size(); id++)
	{
		IOPoint ik = io.at(id);

		fill_n(shortestpaths, (_hgridline_count)*(_vgridline_count), FLT_MAX);
		executeDijkstra(network, network + ik.nodeid, shortestpaths, (_hgridline_count), (_vgridline_count));

		for (int j = 0; j < io.size(); j++)
			distances[id * io.size() + j] = shortestpaths[io.at(j).nodeid];

	}

	delete[] shortestpaths;
}

void rectilinearDistances(float *distances, vector<IOPoint> io) {

	for (int id1 = 0; id1 < io.size(); id1++)
		for (int id2 = 0; id2 < io.size(); id2++)
			distances[id1 * io.size() + id2] = fabs(io.at(id1).x - io.at(id2).x) + fabs(io.at(id1).y - io.at(id2).y);

}

void executeDijkstra(Node* _network, Node* _source, float* _shortestpaths, int hsize, int vsize)
{

	int index, tempindex;
	bool* added_to_queue;
	bool* visited;
	float* labels;
	Node *curnode, *rchild, *tchild, *lchild, *bchild;

	added_to_queue = new bool[vsize*hsize];
	visited = new bool[vsize*hsize];
	labels = new float[vsize*hsize];

	fill_n(added_to_queue, vsize*hsize, false);
	fill_n(visited, vsize*hsize, false);
	fill_n(labels, vsize*hsize, FLT_MAX);

	priority_queue<Node*, vector<Node*>, NodeCompare> pq;

	curnode = _source;
	tempindex = curnode->index;
	added_to_queue[tempindex] = true;
	labels[tempindex] = 0;
	curnode->label = labels[tempindex];
	pq.push(curnode);

	while (!pq.empty())
	{
		curnode = pq.top();
		pq.pop();
		index = curnode->index;

		visited[index] = true;
		added_to_queue[index] = false;
		_shortestpaths[index] = curnode->label;

		if ((curnode->left_node > -1) && (curnode->left_edge < FLT_MAX - 0.0001))
		{

			lchild = _network + curnode->left_node;
			tempindex = lchild->index;

			if (!visited[tempindex] && labels[tempindex] > curnode->label + curnode->left_edge)
			{
				labels[tempindex] = curnode->label + curnode->left_edge;
				lchild->label = labels[tempindex];
				if (!added_to_queue[tempindex])
				{
					pq.push(lchild);
					added_to_queue[tempindex] = true;
				}
			}
		}

		if ((curnode->right_node > -1) && (curnode->right_edge < FLT_MAX - 0.0001))
		{
			rchild = _network + curnode->right_node;
			tempindex = rchild->index;
			if (!visited[tempindex] && labels[tempindex] > curnode->label + curnode->right_edge)
			{
				labels[tempindex] = curnode->label + curnode->right_edge;
				rchild->label = labels[tempindex];
				if (!added_to_queue[tempindex])
				{
					pq.push(rchild);
					added_to_queue[tempindex] = true;
				}
			}
		}

		if ((curnode->top_node > -1) && (curnode->top_edge < FLT_MAX - 0.0001))
		{
			tchild = _network + curnode->top_node;
			tempindex = tchild->index;

			if (!visited[tempindex] && labels[tempindex] > curnode->label + curnode->top_edge) {
				labels[tempindex] = curnode->label + curnode->top_edge;
				tchild->label = labels[tempindex];
				if (!added_to_queue[tempindex]) {
					pq.push(tchild);
					added_to_queue[tempindex] = true;
				}
			}
		}

		if ((curnode->bottom_node > -1) && (curnode->bottom_edge < FLT_MAX - 0.0001))
		{
			bchild = _network + curnode->bottom_node;
			tempindex = bchild->index;

			if (!visited[tempindex] && labels[tempindex] > curnode->label + curnode->bottom_edge) {
				labels[tempindex] = curnode->label + curnode->bottom_edge;
				bchild->label = labels[tempindex];
				if (!added_to_queue[tempindex]) {
					pq.push(bchild);
					added_to_queue[tempindex] = true;
				}
			}
		}
	}
	delete[] visited;
	delete[] added_to_queue;
	delete[] labels;
}

void copyNode(Node *n1, Node *n2) {

	n1->index = n2->index;
	n1->key_X = n2->key_X;
	n1->key_Y = n2->key_Y;

	n1->top_node = n2->top_node;
	n1->left_node = n2->left_node;
	n1->bottom_node = n2->bottom_node;
	n1->right_node = n2->right_node;

	n1->top_edge = n2->top_edge;
	n1->left_edge = n2->left_edge;
	n1->bottom_edge = n2->bottom_edge;
	n1->right_edge = n2->right_edge;
}

void calculateEFEFObj(float &efefobj, float *efefflow, vector<Barrier> blist, vector<IOPoint> iolist, float lwidth, float lheight, bool rectilinear_approx) {

	int efiocount = iolist.size();

	vector<float> vgridlines, hgridlines;

	getGridlinecoord(blist, iolist, vgridlines, hgridlines, lwidth, lheight);

	int _hgridline_count = hgridlines.size();
	int _vgridline_count = vgridlines.size();
	int ycell_count = _hgridline_count - 1;
	int xcell_count = _vgridline_count - 1;

	Node *network = new Node[_hgridline_count * _vgridline_count];
	Cell *cells = new Cell[(xcell_count)*(ycell_count)]; //additional space for NF cells

	createNetwork(network, hgridlines, vgridlines);
	createCells(cells, network, hgridlines, vgridlines);
	getBarrierCellOverlap(network, cells, blist, hgridlines, vgridlines);
	getIONodeOverlap(network, iolist, hgridlines, vgridlines);

	float *distances = new float[iolist.size() * iolist.size()];

	if (rectilinear_approx)
		rectilinearDistances(distances, iolist);
	else
		dijkstraDistances(distances, iolist, network, _hgridline_count, _vgridline_count);

	efefobj = 0;

	for (int ioid1 = 0; ioid1 < efiocount - 1; ioid1++)
		for (int ioid2 = ioid1; ioid2 < efiocount; ioid2++)
			efefobj += distances[ioid1 * (efiocount)+ioid2] * efefflow[ioid1 * efiocount + ioid2];

}