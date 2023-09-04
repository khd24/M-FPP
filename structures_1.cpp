/*
 * structures_1.cpp
 *
 *  Created on: Jun 13, 2020
 *      Author: ketandat
 */

#include "structures_1.h"

void Barrier::setBarrier(float _x, float _y, float _width, float _height, float _congestion)
{
	x = _x;
	y = _y;
	width = _width;
	height = _height;
	congestion = _congestion;
	order_ = -1;
	id_ = -1;
	
}

float Barrier::getMinX(void) { return x; }
float Barrier::getMaxX(void) { return (x + width); }
float Barrier::getMinY(void) { return (y - height); }
float Barrier::getMaxY(void) { return y; }
float Barrier::getWidth(void) { return width; }
float Barrier::getHeight(void) { return height; }
float Barrier::getCongestion(void) { return congestion; }
int Barrier::order() { return order_; }
int Barrier::id() { return id_; }

void Barrier::setWidth(float _width) { width = _width; }
void Barrier::setHeight(float _height) { height = _height; }
void Barrier::setX(float _x) { x = _x; }
void Barrier::setY(float _y) { y = _y; }
void Barrier::setOrder(int order) { order_ = order; }
void Barrier::setId(int id) { id_ = id; }