/*
 *  nonoShape.cpp
 *  LaserStageRaw
 *
 *  Created by say nono on 24.08.10.
 *  Copyright 2010 www.say-nono.com. All rights reserved.
 *
 */

#include "nonoShape.h"

nonoShape::nonoShape(){
}

nonoShape::nonoShape(nonoPoint *points, int amount){
}

int nonoShape::addPoint(nonoPoint *point){
	_points.push_back(*point);
	return _points.size();
}

vector<nonoPoint> nonoShape::getPoints(){
	return _points;
}