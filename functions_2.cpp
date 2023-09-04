/*
 * functions_2.cpp
 *
 *  Created on: Jun 13, 2020
 *      Author: ketandat
 */

#include "functions_2.h"

void createNetwork(Node* _out_nodelist, vector<float> _in_hgridline_coord_list, vector<float> _in_vgridline_coord_list)
{
	int _hgridline_count = _in_hgridline_coord_list.size();
	int _vgridline_count = _in_vgridline_coord_list.size();

	for (int i = 0; i < _hgridline_count; i++) {
		for (int j = 0; j < _vgridline_count; j++) {
			int k = i * _vgridline_count + j;
			_out_nodelist[k].key_X = j;
			_out_nodelist[k].key_Y = i;
			_out_nodelist[k].index = k;
			_out_nodelist[k].v_gridline = _in_vgridline_coord_list.at(j);
			_out_nodelist[k].h_gridline = _in_hgridline_coord_list.at(i);
			_out_nodelist[k].io_Index = -1;

			if (i == 0) {
				_out_nodelist[k].bottom_edge = FLT_MAX;
				_out_nodelist[k].bottom_node = -1;
			}
			else {
				_out_nodelist[k].bottom_edge = _in_hgridline_coord_list.at(i) - _in_hgridline_coord_list.at(i - 1);
				_out_nodelist[k].bottom_node = (i - 1)*_vgridline_count + j;
			}

			if (j == 0) {
				_out_nodelist[k].left_edge = FLT_MAX;
				_out_nodelist[k].left_node = -1;
			}
			else {
				_out_nodelist[k].left_edge = _in_vgridline_coord_list.at(j) - _in_vgridline_coord_list.at(j - 1);
				_out_nodelist[k].left_node = i * _vgridline_count + (j - 1);
			}

			if (i == _hgridline_count - 1) {
				_out_nodelist[k].top_edge = FLT_MAX;
				_out_nodelist[k].top_node = -1;
			}
			else {
				_out_nodelist[k].top_edge = _in_hgridline_coord_list.at(i + 1) - _in_hgridline_coord_list.at(i);
				_out_nodelist[k].top_node = (i + 1)*_vgridline_count + j;
			}

			if (j == _vgridline_count - 1) {
				_out_nodelist[k].right_edge = FLT_MAX;
				_out_nodelist[k].right_node = -1;
			}
			else {
				_out_nodelist[k].right_edge = _in_vgridline_coord_list.at(j + 1) - _in_vgridline_coord_list.at(j);
				_out_nodelist[k].right_node = i * _vgridline_count + (j + 1);
			}

			_out_nodelist[k].bottom_left_cell = 0;
			_out_nodelist[k].bottom_right_cell = 0;
			_out_nodelist[k].top_left_cell = 0;
			_out_nodelist[k].top_right_cell = 0;
		}
	}
}

void createCells(Cell* _out_celllist, Node* _in_nodelist, vector<float> _in_hgridline_coord_list, vector<float> _in_vgridline_coord_list)
{
	int _hgridline_count = _in_hgridline_coord_list.size();
	int _vgridline_count = _in_vgridline_coord_list.size();

	int xcellcount = _vgridline_count - 1;
	for (int i = 0; i < _hgridline_count - 1; i++) {
		for (int j = 0; j < _vgridline_count - 1; j++) {

			int k = i * xcellcount + j;
			_out_celllist[k].barrier_Index = -1;

			_out_celllist[k].key_X = j;
			_out_celllist[k].key_Y = i;
			_out_celllist[k].index = k;
			_out_celllist[k].left_gridline = _in_vgridline_coord_list.at(j);
			_out_celllist[k].right_gridline = _in_vgridline_coord_list.at(j + 1);
			_out_celllist[k].bottom_gridline = _in_hgridline_coord_list.at(i);
			_out_celllist[k].top_gridline = _in_hgridline_coord_list.at(i + 1);

			_in_nodelist[i*_vgridline_count + j].top_right_cell = _out_celllist + k;
			_out_celllist[k].bottom_left_node = _in_nodelist + i * _vgridline_count + j;

			_in_nodelist[i*_vgridline_count + (j + 1)].top_left_cell = _out_celllist + k;
			_out_celllist[k].bottom_right_node = _in_nodelist + i * _vgridline_count + (j + 1);

			_in_nodelist[(i + 1)*_vgridline_count + (j + 1)].bottom_left_cell = _out_celllist + k;
			_out_celllist[k].top_right_node = _in_nodelist + (i + 1)*_vgridline_count + (j + 1);

			_in_nodelist[(i + 1)*_vgridline_count + j].bottom_right_cell = _out_celllist + k;
			_out_celllist[k].top_left_node = _in_nodelist + (i + 1)*_vgridline_count + j;


			if (i == 0)
				_out_celllist[k].bottom_cell = 0;
			else
				_out_celllist[k].bottom_cell = _out_celllist + (i - 1)*xcellcount + j;

			if (i == _hgridline_count - 2)
				_out_celllist[k].top_cell = 0;
			else
				_out_celllist[k].top_cell = _out_celllist + (i + 1)*xcellcount + j;

			if (j == 0)
				_out_celllist[k].left_cell = 0;
			else
				_out_celllist[k].left_cell = _out_celllist + i * xcellcount + (j - 1);

			if (j == _vgridline_count - 2)
				_out_celllist[k].right_cell = 0;
			else
				_out_celllist[k].right_cell = _out_celllist + i * xcellcount + (j + 1);
		}
	}
}

void getBarrierCellOverlap(Node* _network, Cell* _celllist, vector<Barrier> barrierlist, vector<float> _hgridline_coord_array, vector<float> _vgridline_coord_array)
{

	int _hgridline_count = _hgridline_coord_array.size();
	int _vgridline_count = _vgridline_coord_array.size();
	int barriercount = barrierlist.size();

	float * h_gridline_array = new float[_hgridline_count];
	float * v_gridline_array = new float[_vgridline_count];

	copy(_hgridline_coord_array.begin(), _hgridline_coord_array.end(), h_gridline_array);
	copy(_vgridline_coord_array.begin(), _vgridline_coord_array.end(), v_gridline_array);

	for (int i = 0; i < barriercount; i++) {
		Barrier blist = barrierlist.at(i);
		float lvcoord = blist.getMinX();
		float rvcoord = blist.getMaxX();

		float bhcoord = blist.getMinY();
		float thcoord = blist.getMaxY();

		int xmin = quicksearch(v_gridline_array, lvcoord, 0, _vgridline_count - 1);
		int xmax = quicksearch(v_gridline_array, rvcoord, xmin, _vgridline_count - 1);
		int ymin = quicksearch(h_gridline_array, bhcoord, 0, _hgridline_count - 1);
		int ymax = quicksearch(h_gridline_array, thcoord, ymin, _hgridline_count - 1);

		for (int j = xmin; j < xmax; j++)
		{
			for (int k = ymin; k < ymax; k++)
			{
				Cell* c = _celllist + k * (_vgridline_count - 1) + j;
				c->barrier_Index = i;

				int nodeindex = k * _vgridline_count + j;
				Node* n = _network + nodeindex;
				Node* n1;

				if (k > ymin)
				{
					n->right_edge = FLT_MAX;
					if (j > xmin)
						n->left_edge = FLT_MAX;

					n1 = _network + n->right_node;
					n1->left_edge = FLT_MAX;
				}

				if (j > xmin)
				{
					n->top_edge = FLT_MAX;
					if (k > ymin)
						n->bottom_edge = FLT_MAX;

					n1 = _network + n->top_node;
					n1->bottom_edge = FLT_MAX;
				}
			}
		}
	}

	delete[] h_gridline_array;
	delete[] v_gridline_array;
}

void getIONodeOverlap(Node* _nodelist, vector<IOPoint> &ios, vector<float> _hgridline_coord_array, vector<float> _vgridline_coord_array)
{
	int iocount = ios.size();
	int _vgridline_count = _vgridline_coord_array.size();
	int _hgridline_count = _hgridline_coord_array.size();

	float * h_gridline_array = new float[_hgridline_count];
	float * v_gridline_array = new float[_vgridline_count];

	copy(_hgridline_coord_array.begin(), _hgridline_coord_array.end(), h_gridline_array);
	copy(_vgridline_coord_array.begin(), _vgridline_coord_array.end(), v_gridline_array);

	for (int i = 0; i < iocount; i++) {
		IOPoint io = ios.at(i);
		float x = io.x;
		float y = io.y;
		int xindex = quicksearch(v_gridline_array, x, 0, _vgridline_count - 1);
		int yindex = quicksearch(h_gridline_array, y, 0, _hgridline_count - 1);

		int index = (yindex >= 0 && xindex >= 0) ? yindex * _vgridline_count + xindex : -1;
		ios.at(i).nodeid = index;
	}

	delete[] h_gridline_array;
	delete[] v_gridline_array;
}