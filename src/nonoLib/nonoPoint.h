/*
 *  nonoPoint.h
 *  LaserShapes_01
 *
 *  Created by Vincent R. on 20.12.09.
 *  Copyright 2009 www.say-nono.com. All rights reserved.
 *
 */


#ifndef _NONOPOINT
#define _NONOPOINT

class nonoPoint{
	
public:
	float x;
	float y;
	float z;	
	
	nonoPoint();
	nonoPoint(float x, float y);
	nonoPoint(float x, float y, float z);
	nonoPoint(nonoPoint *p);
	
	void setPoints(float x_, float y_);
	void setPoints(float x_, float y_, float z_);
	
	nonoPoint cross( const nonoPoint& p ) const {
		return nonoPoint( y*p.z - z*p.y,
					   z*p.x - x*p.z,
					   x*p.y - y*p.x );
	}
	
	static nonoPoint cross( const nonoPoint& p1, const nonoPoint& p2 ) {
		return nonoPoint( p1.y*p2.z - p1.z*p2.y,
						 p1.z*p2.x - p1.x*p2.z,
						 p1.x*p2.y - p1.y*p2.x );
	}
	
	void interpolate(nonoPoint *p1, nonoPoint *p2, float percent);

//	float* toFloatArray();
	
};

#endif