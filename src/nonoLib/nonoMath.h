/*
 *  nonoMath.h
 *  LaserShapes_01
 *
 *  Created by Vincent R. on 20.12.09.
 *  Copyright 2009 www.say-nono.com. All rights reserved.
 *
 */


#ifndef _NONOMATH
#define _NONOMATH

#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif

#include <math.h>
#include <algorithm>
#include "nonoPoint.h"

#define NONO_MIN(x,y) (x < y ? x : y)
#define NONO_MAX(x,y) (x > y ? x : y)

//#define NONO_ABS(x) (x < 0 ? -x : x)

/*
template <class T> const T& max ( const T& a, const T& b ) {
	return (b<a)?a:b;     // or: return comp(b,a)?a:b; for the comp version
}

template <class T> const T& min ( const T& a, const T& b ) {
	return (b>a)?a:b;     // or: return comp(b,a)?a:b; for the comp version
}
*/

using namespace std;

class nonoMath{

public:
	
// ANGLE STUFF

static double mySin(double val){
	return sin(toRadians(val));
}
static double myCos(double val){
	return cos(toRadians(val));
}
static double myTan(double val){
	return tan(toRadians(val));
}

static double toRadians(double val){
	return val * M_PI / 180;
}

static double toDegrees(double val){
	return val / M_PI * 180.0;
}



	
	
	// GEOMETRY
	
	static double angleRadians(double _x,double _y){
		return fmod(((atan2(-_y,_x)+(M_PI / 2)) + M_PI*2 ) , (M_PI*2)) ;
	}
	static double angle(double _x,double _y){
		return toDegrees(angleRadians(_x,_y));
	}
	static double angle(nonoPoint* p1, nonoPoint* p2){
		return toDegrees(angleRadians(- (p2->x - p1->x) , p2->y - p1->y ));
	}
	static double distance(double _x,double _y){
		return sqrt(_x*_x+_y*_y);
	}
	
	static double angleBetween2Lines(nonoPoint* center, nonoPoint* p1, nonoPoint* p2){
		
		double angle1 = toDegrees(angleRadians( p1->x - center->x , p1->y - center->y ));
		double angle2 = toDegrees(angleRadians( p2->x - center->x , p2->y - center->y ));
		if(angle1<0) angle1 += 360;
		if(angle2<0) angle2 += 360;		
		double angle = max(angle1,angle2) - min(angle1,angle2);
		if(angle>180) angle = 360-angle;		
		
		return angle;
	}
	
	static void angleLengthToPoint(double angle,double length, nonoPoint *p){
		p->x = 0;
		p->y = -length;
		rotatePoint(p,angle,p);
	}
	
	static void rotatePoint(nonoPoint *point, double angle, nonoPoint *p=NULL){
		if(p==NULL) p = point;
		angle *= -1;
		double x1 = mySin(angle+90)*point->x;
		double y1 = myCos(angle+90)*point->x;
		double x2 = mySin(angle)*point->y;
		double y2 = myCos(angle)*point->y;
		p->x = x1+x2;
		p->y = y1+y2;
	}
	 
	static void rotateCoordinates(double x,double y, double angle, nonoPoint *p){
		angle *= -1;
		double x1 = mySin(angle+90)*x;
		double y1 = myCos(angle+90)*x;
		double x2 = mySin(angle)*y;
		double y2 = myCos(angle)*y;
		p->x = x1+x2;
		p->y = y1+y2;		
//		nonoPoint &p;
//		return *(n->xew nonoPoint(x1+x2,y1+y2));
//		return *(new nonoPoint(100,200));
	}
	
	static bool getLineLineIntersection(nonoPoint *p1a, nonoPoint *p1b, nonoPoint *p2a, nonoPoint *p2b, nonoPoint *pRes){
		double xD1,yD1,xD2,yD2,xD3,yD3;  
		double dot,deg,len1,len2;  
		double segmentLen1,segmentLen2;  
		double ua,ub,div;  
		
		// calculate differences  
		xD1=p1b->x-p1a->x;  
		xD2=p2b->x-p2a->x;  
		yD1=p1b->y-p1a->y;  
		yD2=p2b->y-p2a->y;  
		xD3=p1a->x-p2a->x;  
		yD3=p1a->y-p2a->y;    
		
		// calculate the lengths of the two lines  
		len1=sqrt(xD1*xD1+yD1*yD1);  
		len2=sqrt(xD2*xD2+yD2*yD2);  
		
		// calculate angle between the two lines.  
		dot=(xD1*xD2+yD1*yD2); // dot product  
		deg=dot/(len1*len2);  
		
		// if abs(angle)==1 then the lines are parallell,  
		// so no intersection is possible  
		if(fabs(deg)==1) return false;  
		
		// find intersection Pt between two lines  
		div=yD2*xD1-xD2*yD1;
		if(div == 0) return false;
		
		ua=(xD2*yD3-yD2*xD3)/div;  
		ub=(xD1*yD3-yD1*xD3)/div;  
		
		nonoPoint pt;
		pt.x=p1a->x+ua*xD1;  
		pt.y=p1a->y+ua*yD1;  
		
		// calculate the combined length of the two segments  
		// between Pt-p1 and Pt-p2  
		xD1=pt.x-p1a->x;  
		xD2=pt.x-p1b->x;  
		yD1=pt.y-p1a->y;  
		yD2=pt.y-p1b->y;  
		segmentLen1=sqrt(xD1*xD1+yD1*yD1)+sqrt(xD2*xD2+yD2*yD2);  
		
		// calculate the combined length of the two segments  
		// between Pt-p3 and Pt-p4  
		xD1=pt.x-p2a->x;  
		xD2=pt.x-p2b->x;  
		yD1=pt.y-p2a->y;  
		yD2=pt.y-p2b->y;  
		segmentLen2=sqrt(xD1*xD1+yD1*yD1)+sqrt(xD2*xD2+yD2*yD2);  
		
		// if the lengths of both sets of segments are the same as  
		// the lenghts of the two lines the point is actually  
		// on the line segment.  
		
		// if the point isn’t on the line, return null  
		double abs1 = fabs(len1-segmentLen1);
		double abs2 = fabs(len2-segmentLen2);
		if(abs1>0.01 || abs2>0.01)  
			return false;  
		
		// return the valid intersection  
		pRes->x = pt.x;
		pRes->y = pt.y;
		
		return true;  
	}  
	
	static double distancePointToLine(nonoPoint *p, nonoPoint *lp1, nonoPoint *lp2){
		double a = angleRadians(p->x-lp1->x,p->y-lp1->y) - angleRadians(lp2->x-lp1->x,lp2->y-lp1->y);
		double d1 = distance(p->x-lp1->x,p->y-lp1->y);
		double d2 = (cos(a)*d1);
		double val = sqrt(d1*d1 - d2*d2);
		return val;
	}
	
	static double distancePointToPoint(nonoPoint *p1, nonoPoint *p2){
		double dx = p1->x-p2->x;
		double dy = p1->y-p2->y;
		double dz = p1->z-p2->z;
		double val = sqrt(dx*dx+dy*dy+dz*dz);
		return val;
	}		
	
	static void nearestPointOnLine(nonoPoint *p, nonoPoint *lp1, nonoPoint *lp2, nonoPoint *pRes){
		double a = angleRadians(p->x-lp1->x,p->y-lp1->y) - angleRadians(lp2->x-lp1->x,lp2->y-lp1->y);
		double d1 = distance(p->x-lp1->x,p->y-lp1->y);
		double d2 = (cos(a)*d1);
		double dline = distance(lp2->x-lp1->x,lp2->y-lp1->y);
		double percent = d2/dline;
		pRes->x = interpolatedouble(lp1->x,lp2->x,percent);
		pRes->y = interpolatedouble(lp1->y,lp2->y,percent);
	}
	
	
	// BEZIER STUFF
	
	static double getBezierPoint3P1D  (double f, double x1, double a1, double x2) {
		double b, e, ee, ff; 
		e = 1-f; 
		ee = e*e;
		ff = f*f;
		b = 2*f*e;
		return x2*ff+a1*b+x1*ee;
	}
	static double getBezierPoint4P1D  (double f, double x1, double a1, double a2, double x2) {
		double a, b, c, d, e, ee, ff; 
		e = 1-f; 
		ee = e*e;
		ff = f*f;
		a=f*ff;
		b=3*ff*e;
		c=3*f*ee;
		d=e*ee;			
		return (x2*a+a2*b+a1*c+x1*d);
	}
	
	static void getBezierPoint3P2D  (double f, nonoPoint *p1, nonoPoint *a1, nonoPoint *p2, nonoPoint *pRes) {
		pRes->x = getBezierPoint3P1D(f,p1->x,a1->x,p2->x);
		pRes->y = getBezierPoint3P1D(f,p1->y,a1->y,p2->y);
	}
	static void getBezierPoint4P2D  (double f, nonoPoint *p1, nonoPoint *a1, nonoPoint *a2, nonoPoint *p2, nonoPoint *pRes) {
		pRes->x = getBezierPoint4P1D(f,p1->x,a1->x,a2->x,p2->x);
		pRes->y = getBezierPoint4P1D(f,p1->y,a1->y,a2->y,p2->y);
	}
	static void getBezierPoint3P3D  (double f, nonoPoint *p1, nonoPoint *a1, nonoPoint *p2, nonoPoint *pRes) {
		pRes->x = getBezierPoint3P1D(f,p1->x,a1->x,p2->x);
		pRes->y = getBezierPoint3P1D(f,p1->y,a1->y,p2->y);
		pRes->z = getBezierPoint3P1D(f,p1->z,a1->z,p2->z);
	}

	static double getBezierLength3P2D  (nonoPoint *p1, nonoPoint *a1, nonoPoint *p2, int details = 50) {
		double ebd = 0;
		double step = 1.0/(double)details;
		nonoPoint p;
		getBezierPoint3P2D(0, p1, a1, p2,&p);
		for (double i = step; i <= 1; i+=step) {
//			double percent = i/100.0;
			nonoPoint pp;
			getBezierPoint3P2D(i, p1, a1, p2, &pp);
			ebd += distance(p.x-pp.x, p.y-pp.y);
			p.x = pp.x;
			p.y = pp.y;
		}
		return ebd;
		//			return distance(getBezierLength3P1D(p1.x,a1.x,p2.x,details),getBezierLength3P1D(p1.y,a1.y,p2.y,details));
	}
	
	/*
	static getBezierPartLength3P2D  (p1:Point, a1:Point, p2:Point,details:int = 50):Array {
		var step:double = 100/details;
		var p:Point = NonoMath.getBezierPoint3P2D(0, p1, a1, p2);
		var arr:Array = [];
		for (var i : double = step; i <= 100; i+=step) {
			var percent:double = i/100;
			var pp:Point = NonoMath.getBezierPoint3P2D(percent, p1, a1, p2);
			arr.push(NonoMath.distance(p.x-pp.x, p.y-pp.y));
			p = pp;
		}
		return arr;
		//			return distance(getBezierLength3P1D(p1.x,a1.x,p2.x,details),getBezierLength3P1D(p1.y,a1.y,p2.y,details));
	}
	*/
	
	
	static double getBezierLength3P1D  (double p1, double a1, double p2, int details = 50) {
		double result = 0;
		double p1res = getBezierPoint3P1D(0,p1,a1,p2);
		double p2res;
		double step = 1.0/details;
		for (double i = step; i <= 1; i+=step) {
			p2res = getBezierPoint3P1D(i,p1,a1,p2);
			result += fabs(p2res-p1res);
			p1res = p2res;
		}
		return result;
	}
	
	static double getBezierLength4P2D  (nonoPoint *p1, nonoPoint *a1, nonoPoint *a2, nonoPoint *p2, int details = 50) {
		return distance(getBezierLength4P1D(p1->x,a1->x,a2->x,p2->x,details),getBezierLength4P1D(p1->y,a1->y,a2->y,p2->y,details));
	}
	static double getBezierLength4P1D  (double p1, double a1, double a2, double p2, int details = 50) {
		double result = 0;
		double p1res = getBezierPoint4P1D(0,p1,a1,a2,p2);
		double p2res;
		double step = 1.0/details;
		for (double i = step; i < 1; i+=step) {
			p2res = getBezierPoint4P1D(i,p1,a1,a2,p2);
			result += fabs(p2res-p1res);
			p1res = p2res;
		}
		return result;
	}
	
	/*
	static void getBezierArcLength4P2Dintern  (points:Array, error:double, length:Array) {
		
		var len:double = 0;
		var chord:double = 0;
		len += Point.distance(Point(points[0]), Point(points[1]));
		len += Point.distance(Point(points[1]), Point(points[2]));
		len += Point.distance(Point(points[2]), Point(points[3]));
		chord = Point.distance(Point(points[0]), Point(points[3]));
		trace("(len-chord) > error : " + ((len-chord) > error) + "     -> " + (len-chord) + "			(len : " + len + "    chord : " + chord);
		if((len-chord) > error)
		{
			var ret:Array = devideBezier4P2D(.5, Point(points[0]), Point(points[1]), Point(points[2]), Point(points[3]));
			getBezierArcLength4P2Dintern([ret[0],ret[1],ret[2],ret[3]], error, length);        // try left side
			getBezierArcLength4P2Dintern([ret[4],ret[5],ret[6],ret[7]], error, length);       // try right side
			return;
		}
		length[0] = length[0] + len;
	}
	
	static double getBezierArcLength4P2D  (points:Array, double error) {
		var length:Array = [0];
		getBezierArcLength4P2Dintern(points, error, length);
		return length[0];
	}
	*/
	
	 
	static void devideBezier3P1D (double percent, double x1, double a1, double x2, double *res) {
		percent = max(0.0, (double)min(1.0, percent));
		double percent_left = percent;
		double percent_right = 1-percent;
		double q1,q2;
		q2 = getBezierPoint3P1D(percent,x1,a1,x2);
		double r0,r1;
		r0 = q2;
		q1 = (x1 * percent_right) + (a1 * percent_left);
		r1 = (a1 * percent_right) + (x2 * percent_left);
		res[0] = q1;
		res[1] = q2;
		res[2] = r1;
	}
	
	static void devideBezier3P2D (double percent, nonoPoint *p1, nonoPoint *a1, nonoPoint *p2, nonoPoint *result) {
		double *resX = new double[3];
		double *resY = new double[3];
		devideBezier3P1D(percent,p1->x,a1->x,p2->x,resX);
		devideBezier3P1D(percent,p1->y,a1->y,p2->y,resY);
		for(int i=0;i<3;i++){
			result->x = resX[i];
			result->y = resY[i];
		}
	}

	static void devideBezier3P3D (double percent, nonoPoint *p1, nonoPoint *a1, nonoPoint *p2, nonoPoint *result) {
		double *resX = new double[3];
		double *resY = new double[3];
		double *resZ = new double[3];
		devideBezier3P1D(percent,p1->x,a1->x,p2->x,resX);
		devideBezier3P1D(percent,p1->y,a1->y,p2->y,resY);
		devideBezier3P1D(percent,p1->z,a1->z,p2->z,resZ);
		for(int i=0;i<3;i++){
			result->x = resX[i];
			result->y = resY[i];
			result->z = resZ[i];
		}
	}
	
	
	 static void devideBezier4P1D (double percent, double x1, double a1, double a2, double x2, double *result) {
		percent = max(0.0, min(1.0, percent));
		double percent_left = percent;
		double percent_right = 1.0-percent;
		double q0, q1, q2, q3;
		q0 = x1;
		q3 = getBezierPoint4P1D(percent,x1,a1,a2,x2);
	 
		double r0,r1,r2,r3;
		r3 = x2;
		r0 = q3;
	 
		q1 = (q0 * percent_right) + (a1 * percent_left);
		q2 = (q1 * percent_right) + (a2 * percent_left + a1 * percent_right) * percent_left;
		r2 = (a2 * percent_right) + (r3 * percent_left);
		r1 = (r2 * percent_left) + ( a1 * percent_right + a2 * percent_left) * percent_right;
		result[0] = q1;
		result[1] = q2;
		result[2] = q3;
		result[3] = r1;
		result[4] = r2;		 
		 
	 }
	
	
	static void devideBezier4P2D  (double percent, nonoPoint *p1, nonoPoint *a1, nonoPoint *a2, nonoPoint *p2, nonoPoint *result) {
		double *resX = new double[5];
		double *resY = new double[5];		
		devideBezier4P1D(percent,p1->x,a1->x,a2->x,p2->x,resX);
		devideBezier4P1D(percent,p1->y,a1->y,a2->y,p2->y,resY);
		for(int i=0;i<5;i++){
			result[i].x = resX[i];
			result[i].y = resY[i];
		}
	}
	
	
	// UTILS
	
	static double interpolatedouble(double num1,double num2,double value){
		return num1*(1.0f-value) + num2*value;
	}
	static float interpolate(float num1,float num2,float value){
		return num1*(1.0f-value) + num2*value;
	}

	static void interpolatePoints(nonoPoint* p1,nonoPoint* p2,nonoPoint* pResult, double value){
		pResult->x = interpolate(p1->x,p2->x,value);
		pResult->y = interpolate(p1->y,p2->y,value);
		pResult->z = interpolate(p1->z,p2->z,value);		
	}
	
	/*	static interpolatedoubleArray(nums1:Array,nums2:Array,value:double):Array{
		var result:Array = [nums1.length];
		for (var i : double = 0; i < nums1.length; i++) {
			var res:double = interpolatedouble(nums1[i],nums2[i],value);
			result[i] = res;
		}
		return result;
	}
*/
};

#endif