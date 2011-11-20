/*
 *  nonoPlane.cpp
 *  RoomShadowTestsNew_2
 *
 *  Created by Vincent R. on 08.01.10.
 *  Copyright 2010 www.say-nono.com. All rights reserved.
 *
 */

#include "nonoPlane.h"


nonoPlane::nonoPlane(){
}

nonoPlane::nonoPlane(nonoPoint *p1,nonoPoint *p2, nonoPoint *p3){
//	a = p1->y*(p2->z-p3->z) + p2->y*(p3->z-p1->z) + p3->y*(p1->z-p2->z);
//	b = p1->z*(p2->x-p3->x) + p2->z*(p3->x-p1->x) + p3->z*(p1->x-p2->x);
//	c = p1->x*(p2->y-p3->y) + p2->x*(p3->y-p1->y) + p3->x*(p1->y-p2->y);
//	d = -( p1->x*(p2->y*p3->z - p3->y*p2->z) +
//						p2->x*(p3->y*p1->z - p1->y*p3->z) +
//						p3->x*(p1->y*p2->z - p2->y*p1->z) );
////	nonoGeom::normal(p1,p2,p3,&normal);	
//	normal.x = a;
//	normal.y = b;
//	normal.z = c;	
//	nonoGeom::normalise(&normal);
	
	amountPoints = 3;
	points[0] = p1;
	points[1] = p2;
	points[2] = p3;
	visible = true;
	calcProps();
}


nonoPlane::nonoPlane(nonoPoint *ps[],int amount){
	setPoints(ps,amount);
	/*
	amountPoints = amount;
	for (int i = 0; i < amount; i++) {
		points[i] = ps[i];
	}
	
	a = points[0]->y*(points[1]->z-points[2]->z) + points[1]->y*(points[2]->z-points[0]->z) + points[2]->y*(points[0]->z-points[1]->z);
	b = points[0]->z*(points[1]->x-points[2]->x) + points[1]->z*(points[2]->x-points[0]->x) + points[2]->z*(points[0]->x-points[1]->x);
	c = points[0]->x*(points[1]->y-points[2]->y) + points[1]->x*(points[2]->y-points[0]->y) + points[2]->x*(points[0]->y-points[1]->y);
	d = -( points[0]->x*(points[1]->y*points[2]->z - points[2]->y*points[1]->z) +
		  points[1]->x*(points[2]->y*points[0]->z - points[0]->y*points[2]->z) +
		  points[2]->x*(points[0]->y*points[1]->z - points[1]->y*points[0]->z) );
	nonoGeom::normal(points[1],points[1],points[2],&normal);	
	visible = true;	
	*/
}

nonoPlane::~nonoPlane(){
}

void nonoPlane::setPoints(nonoPoint *ps[],int amount){
	
	amountPoints = amount;
	for (int i = 0; i < amount; i++) {
		points[i] = ps[i];
	}
	calcProps();
}

void nonoPlane::calcProps(){
	a = points[0]->y*(points[1]->z-points[2]->z) + points[1]->y*(points[2]->z-points[0]->z) + points[2]->y*(points[0]->z-points[1]->z);
	b = points[0]->z*(points[1]->x-points[2]->x) + points[1]->z*(points[2]->x-points[0]->x) + points[2]->z*(points[0]->x-points[1]->x);
	c = points[0]->x*(points[1]->y-points[2]->y) + points[1]->x*(points[2]->y-points[0]->y) + points[2]->x*(points[0]->y-points[1]->y);
	d = -( points[0]->x*(points[1]->y*points[2]->z - points[2]->y*points[1]->z) +
		  points[1]->x*(points[2]->y*points[0]->z - points[0]->y*points[2]->z) +
		  points[2]->x*(points[0]->y*points[1]->z - points[1]->y*points[0]->z) );
	
	double l = sqrt(a*a + b*b + c*c);
	a /= l;
	b /= l;
	c /= l;
	d /= l;	
	
	normal.x = a;
	normal.y = b;
	normal.z = c;	
	
	visible = true;
	
	float minX = std::numeric_limits<double>::max();
	float minY = std::numeric_limits<double>::max();
	float minZ = std::numeric_limits<double>::max();
	float maxX = std::numeric_limits<double>::min();
	float maxY = std::numeric_limits<double>::min();
	float maxZ = std::numeric_limits<double>::min();
	
	for (int i = 0; i < amountPoints; i++) {
		minX = std::min(points[i]->x,minX);
		minY = std::min(points[i]->y,minY);
		minZ = std::min(points[i]->z,minZ);
		maxX = std::max(points[i]->x,maxX);
		maxY = std::max(points[i]->y,maxY);
		maxZ = std::max(points[i]->z,maxZ);
	}
	
	center.x = nonoMath::interpolatedouble(minX,maxX,.5);
	center.y = nonoMath::interpolatedouble(minY,maxY,.5);
	center.z = nonoMath::interpolatedouble(minZ,maxZ,.5);
	
}

bool nonoPlane::getIntersection(nonoPoint *p1,nonoPoint *p2, nonoPoint *intersection){
//	nonoGeom::normalise(intersection);
	
	double val1 = a*p1->x + b * p1->y + c * p1->z + d;
	double val2 = a * (p1->x-p2->x) + b * (p1->y-p2->y) + c * (p1->z - p2->z);
//	if(val1 == 0) printf("POINT1 == 0\n");
//	if(val2 == 0) printf("POINT2 == 0\n");
	double u = val1/val2;
	nonoGeom::reset(intersection);
//	intersection = p1;
//	printf("   U = %f\n",u);
	nonoGeom::subtract(p2,p1,intersection);
	nonoGeom::multiply(intersection,u,intersection);
	nonoGeom::add(intersection,p1);
	return true;
}

void nonoPlane::move(double mx, double my, double mz){
	for (int i = 0; i < amountPoints; i++) {
		points[i]->x += mx;
		points[i]->y += my;
		points[i]->z += mz;
	}
	calcProps();
}