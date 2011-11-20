/*
 *  nonoPlane.h
 *  RoomShadowTestsNew_2
 *
 *  Created by Vincent R. on 08.01.10.
 *  Copyright 2010 www.say-nono.com. All rights reserved.
 *
 */


#ifndef _NONOPLANE
#define _NONOPLANE

#include "nonoPoint.h"
#include "nonoMath.h"
#include "nonoGeom.h"


class nonoPlane{
	
public:
	
	nonoPlane();
	nonoPlane(nonoPoint *p1,nonoPoint *p2, nonoPoint *p3);
	nonoPlane(nonoPoint *p1[],int amount);
	~nonoPlane();
	void setPoints(nonoPoint *p1[],int amount);
	bool getIntersection(nonoPoint *p1,nonoPoint *p2, nonoPoint *intersection);
	void move(double mx, double my, double mz = 0);
	void calcProps();
	
	int amountPoints;
	double a,b,c,d;
	nonoPoint normal;
	nonoPoint *points[10];
	nonoPoint center;
	bool visible;
	
};

#endif

