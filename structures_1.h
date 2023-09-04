/*
 * structures_1.h
 *
 *  Created on: Oct 14, 2012
 *      Author: ketandat
 */

#ifndef STRUCTURES_1_H_
#define STRUCTURES_1_H_

using namespace std;

////////////////////////////////////////////////////////////////////////////

class Barrier
{
	float x, y, width, height, congestion;
	int order_, id_;

public:

/*
 * First quadrant of XY plane.
 * Coordinates of Top-left corner of the barrier are provided as _x, _y.
 */
	void setBarrier(float _x, float _y, float _width, float _height, float _congestion);

	float getMinX(void);
	float getMaxX(void);
	float getMinY(void);
	float getMaxY(void);
	float getWidth(void);
	float getHeight(void);
	float getCongestion(void);
	int order();
	int id();

	void setWidth(float _width);
	void setHeight(float _height);
	void setX(float _x);
	void setY(float _y);
	void setOrder(int order);
	void setId(int id);
};

struct IOPoint
{
	float x;
	float y;
	int nodeid;
	bool isNF;
};

struct Layout
{
	int barriercount, iocount, nfcount;
	float width, height;
	Barrier * barrierlist;
	Barrier * nflist;
	IOPoint * iolist;
	float *efef;
	float *efnf;
	float *nfnf;
};

#endif /* STRUCTURES_1_H_ */
