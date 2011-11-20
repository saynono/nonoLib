/*
 *  nonoPoint.cpp
 *  LaserShapes_01
 *
 *  Created by Vincent R. on 20.12.09.
 *  Copyright 2009 www.say-nono.com. All rights reserved.
 *
 */

#include "nonoPoint.h"

nonoPoint::nonoPoint(){
	x = 0;
	y = 0;
	z = 0;
}


nonoPoint::nonoPoint(float x_, float y_){
	x = x_;
	y = y_;
	z = 0;
}


nonoPoint::nonoPoint(float x_, float y_, float z_){
	x = x_;
	y = y_;
	z = z_;
}

nonoPoint::nonoPoint(nonoPoint *p){
	x = p->x;
	y = p->y;
	z = p->z;
}

void nonoPoint::setPoints(float x_, float y_){
	x = x_;
	y = y_;
}

void nonoPoint::setPoints(float x_, float y_, float z_){
	x = x_;
	y = y_;
	z = z_;
}

void nonoPoint::interpolate(nonoPoint *p1, nonoPoint *p2, float percent){
	//	printf("  -> %f    -> %f\n",(p2->x - p1->x)*percent+ p1->x,percent);
	x = (p2->x - p1->x)*percent + p1->x;
	y = (p2->y - p1->y)*percent + p1->y;
	z = (p2->z - p1->z)*percent + p1->z;	
}


//float* nonoPoint::toFloatArray(){
//	float arr[3] = {(float)x,(float)y,(float)z};
//	return arr;
//}
