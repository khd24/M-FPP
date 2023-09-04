/*
 * functions_1.cpp
 *
 *  Created on: Jun 13, 2020
 *      Author: ketandat
 */

#include "functions_1.h"

Layout* readLayoutFromFile(std::string _layoutfile, int& _problemsize) {

	int probsize, efsize, nfsize, iosize;
	float lwidth, lheight;
	Layout* layout = 0;

	ifstream myfile(_layoutfile.c_str());

	if (!myfile) {
		std::cerr << "Error: input file not found: " << _layoutfile.c_str() << std::endl;
		exit(-1);
	}

	while (myfile.is_open() && myfile.good() && !myfile.eof()) {

		myfile >> probsize;
		myfile >> lwidth;
		myfile >> lheight;
		myfile >> efsize;
		myfile >> nfsize;

		_problemsize = probsize;
		iosize = efsize;
		layout = new Layout[probsize];
		for (int i = 0; i < probsize; i++) {

			Barrier * barriers = new Barrier[efsize];
			Barrier * nfs = new Barrier[nfsize];
			IOPoint * iops = new IOPoint[iosize];

			float *efefi = new float[iosize * iosize];
			float *efnfi = new float[iosize * nfsize];
			float *nfnfi = new float[nfsize * nfsize];

			for (int j = 0; j < efsize; j++) {
				float* a = new float[5];

				for (int k = 0; k < 4; k++)
					myfile >> a[k];

				barriers[j].setBarrier(a[0], a[1], a[2], a[3], 1.0);
				delete[] a;
			}

			for (int j = 0; j < nfsize; j++) {
				float* a = new float[5];

				for (int k = 0; k < 4; k++)
					myfile >> a[k];

				nfs[j].setBarrier(a[0], a[1], a[2], a[3], 1.0);
				//nfs[j].setBarrier(-1,-1,100,100,1.0);
				delete[] a;
			}

			for (int j = 0; j < iosize; j++) {
				float* b = new float[2];
				for (int k = 0; k < 2; k++)
					myfile >> b[k];

				iops[j].x = b[0];
				iops[j].y = b[1];
				iops[j].nodeid = -1;
				iops[j].isNF = false;

				delete[] b;
			}

			for (int j = 0; j < iosize * iosize; j++)
				myfile >> efefi[j];

			for (int j = 0; j < iosize * nfsize; j++)
				myfile >> efnfi[j];

			for (int j = 0; j < nfsize * nfsize; j++)
				myfile >> nfnfi[j];

			layout[i].width = lwidth;
			layout[i].height = lheight;
			layout[i].barriercount = efsize;
			layout[i].nfcount = nfsize;
			layout[i].iocount = iosize;
			layout[i].barrierlist = barriers;
			layout[i].nflist = nfs;
			layout[i].iolist = iops;
			layout[i].efef = efefi;
			layout[i].efnf = efnfi;
			layout[i].nfnf = nfnfi;

		}
		myfile.close();
	}
	return layout;
}

void getGridlinecoord(std::vector<Barrier> barrierlist, std::vector<IOPoint> iolist, std::vector<float> &_vgridline_coord, std::vector<float> &_hgridline_coord, float _width, float _height) {
	int efcount = barrierlist.size();
	int iocount = iolist.size();

	getHGridlinecoord(barrierlist, iolist, _hgridline_coord, _width, _height);

	getVGridlinecoord(barrierlist, iolist, _vgridline_coord, _width, _height);
}

void getHGridlinecoord(std::vector<Barrier> barrierlist, std::vector<IOPoint> iolist, std::vector<float> &_hgridline_coord, float _width, float _height) {
	int efcount = barrierlist.size();
	int iocount = iolist.size();

	for (int i = 0; i < efcount; i++) {
		_hgridline_coord.push_back(barrierlist[i].getMinY());
		_hgridline_coord.push_back(barrierlist[i].getMaxY());
	}

	for (int i = 0; i < iocount; i++) {
		_hgridline_coord.push_back(iolist[i].y);
	}

	_hgridline_coord.push_back(0);
	_hgridline_coord.push_back(_height);

	sort(_hgridline_coord.begin(), _hgridline_coord.end());

	std::vector<float>::iterator it1 = unique(_hgridline_coord.begin(), _hgridline_coord.end(), float_eq);
	_hgridline_coord.resize(distance(_hgridline_coord.begin(), it1));
}

void getVGridlinecoord(std::vector<Barrier> barrierlist, std::vector<IOPoint> iolist, std::vector<float> &_vgridline_coord, float _width, float _height) {
	int efcount = barrierlist.size();
	int iocount = iolist.size();

	for (int i = 0; i < efcount; i++) {
		_vgridline_coord.push_back(barrierlist[i].getMinX());
		_vgridline_coord.push_back(barrierlist[i].getMaxX());
	}

	for (int i = 0; i < iocount; i++) {
		_vgridline_coord.push_back(iolist[i].x);
	}

	_vgridline_coord.push_back(0);
	_vgridline_coord.push_back(_width);

	sort(_vgridline_coord.begin(), _vgridline_coord.end());

	std::vector<float>::iterator it2 = std::unique(_vgridline_coord.begin(), _vgridline_coord.end(), float_eq);
	_vgridline_coord.resize(std::distance(_vgridline_coord.begin(), it2));
}

/* swap values at two pointers */
void swap(int *x, int *y)
{
	int temp;
	temp = *x;
	*x = *y;
	*y = temp;
}

/* print permutations of string */
void permute(vector<int*> &plist, int *a, int i, int n)
{
	int j;
	if (i == n) {
		int *b = new int[n];
		copy(a, a + n, b);
		plist.push_back(b);
	}
	else
	{
		for (j = i; j < n; j++)
		{
			swap((a + i), (a + j));
			permute(plist, a, i + 1, n);
			swap((a + i), (a + j));
		}
	}
}

void writeLayoutToFile(const string _layoutfile, int efcount, int nfcount, Barrier* blist, Barrier* nflist, IOPoint* iolist) {

	ofstream layoutstream;
	layoutstream.open(_layoutfile.c_str(), std::ios_base::app);

	for (int efid = 0; efid < efcount; efid++) {
		Barrier b = blist[efid];
		layoutstream << b.getMinX() << "\t" << b.getMaxY() << "\t" << b.getWidth() << "\t" << b.getHeight() << std::endl;
	}

	for (int nfid = 0; nfid < nfcount; nfid++) {
		Barrier nf = nflist[nfid];
		layoutstream << nf.getMinX() << "\t" << nf.getMaxY() << "\t" << nf.getWidth() << "\t" << nf.getHeight() << std::endl;
	}

	for (int efid = 0; efid < efcount; efid++) {
		IOPoint io = iolist[efid];
		layoutstream << io.x << "\t" << io.y << std::endl;
	}

	layoutstream.close();

}

void removeDuplicates(vector<float>& vector) {

	sort(vector.begin(), vector.end());

	std::vector<float>::iterator it2 = unique(vector.begin(), vector.end(), float_eq);
	vector.resize(distance(vector.begin(), it2));
}