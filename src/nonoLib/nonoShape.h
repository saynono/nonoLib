/*
 *  nonoShape.h
 *  LaserStageRaw
 *
 *  Created by say nono on 24.08.10.
 *  Copyright 2010 www.say-nono.com. All rights reserved.
 *
 */

#pragma once

#include "ofMain.h"
#include "nonoPoint.h"

class nonoShape{
	
public:
	
	nonoShape();
	nonoShape(nonoPoint *points, int amount);
	
	int addPoint(nonoPoint *point);
	vector<nonoPoint> getPoints();
	
	
	
private:
	vector<nonoPoint>		_points;
	
	
};