/*
 *  nonoGeom.h
 *  LaserShapes_01
 *
 *  Created by Vincent R. on 21.12.09.
 *  Copyright 2009 www.say-nono.com. All rights reserved.
 *
 */



#ifndef _NONOGEOM
#define _NONOGEOM

#include "nonoMath.h"
#include "nonoPoint.h"
#include <vector>

class nonoGeom{
	
public:

	static void reset(nonoPoint *p){
		p->x = 0;
		p->y = 0;
		p->x = 0;		
	}
	
	
	static void multiply(nonoPoint *p, double val, nonoPoint *res = NULL){
		if(res == NULL) res = p;
		res->x = p->x * val;
		res->y = p->y * val;
		res->z = p->z * val;		
	}
	
	static void divide(nonoPoint *p, double val, nonoPoint *res = NULL){
		multiply(p,1.0/val,res);
	}
	
	
	static double getAngle(nonoPoint *p){
		return nonoMath::angle(p->x, p->y);
	}
	
	static void rotate(nonoPoint *p, double angle, nonoPoint *res = NULL){
		nonoMath::rotatePoint(p, angle,res);
		/*
		 if(res == NULL) res = p;
		 var p:Point = nonoMath::rotatePoint(p, angle,res);
		 x = p.x;
		 y = p.y;
		 */
	}
	
	static void rotateZ(nonoPoint *p, double angle, nonoPoint *res = NULL){
		if(res == NULL) res = p;
		nonoMath::rotateCoordinates(p->x,p->y, angle,res);
	}
	
	static void rotateX(nonoPoint *p, double angle, nonoPoint *res = NULL){
		if(res == NULL) res = p;
		double x = res->x;
		double y = res->y;
		double z = res->z;
		nonoMath::rotateCoordinates(p->y,p->z, angle,res);
		y = res->x;
		z = res->y;
		res->x = x;
		res->y = y;
		res->z = z;
	}
	
	static void rotateY(nonoPoint *p, double angle, nonoPoint *res = NULL){
		if(res == NULL) res = p;
		double x = res->x;
		double y = res->y;
		double z = res->z;
		nonoMath::rotateCoordinates(p->x,p->z, angle,res);
		x = res->x;
		z = res->y;
		res->x = x;
		res->y = y;
		res->z = z;
	}
		
	static double getLength(nonoPoint *p, nonoPoint *p2 = NULL){
		double l;
		if(p2 == NULL) l = sqrt(p->x*p->x + p->y*p->y + p->z*p->z);
		else{
			double x = p2->x-p->x;
			double y = p2->y-p->y;
			double z = p2->z-p->z;
			l = sqrt(x*x+y*y+z*z);
		}
		return l;
	}
	
	static void addLength(nonoPoint *p, double len, nonoPoint *res = NULL){
		if(res == NULL) res = p;
		double l = getLength(p) + len;
		double a = nonoMath::angle(p->x, p->y);
		res->x = nonoMath::mySin(a)*l;
		res->y = nonoMath::myCos(a)*l;			
	}
	 
	static void addPoint(nonoPoint *p, nonoPoint *res){
		res->x += p->x;
		res->y += p->y;
		res->z += p->z;
	}
	
	 
	static void interpolate(nonoPoint *p1, nonoPoint *p2, double value, nonoPoint *res){
		res->x = nonoMath::interpolatedouble(p1->x, p2->x, value);
		res->y = nonoMath::interpolatedouble(p1->y, p2->y, value);
		res->z = nonoMath::interpolatedouble(p1->z, p2->z, value);
	}
	 
	 
	 static void normalise(nonoPoint *p) {
		 double d = getLength(p);
		 p->x /= d;
		 p->y /= d;
		 p->z /= d;
	 }
	 	 
	static void add(nonoPoint *p, double len, nonoPoint *res = NULL){
		if(res == NULL) res = p;
		res->x = p->x + len;
		res->y = p->y + len;
		res->z = p->z + len;
	 }
	 
	static void add(nonoPoint *p, nonoPoint *p2, nonoPoint *res = NULL){
		if(res == NULL) res = p;
		res->x = p->x + p2->x;
		res->y = p->y + p2->y;
		res->z = p->z + p2->z;
	}
	
	static void subtract(nonoPoint *p, double len, nonoPoint *res = NULL){
		if(res == NULL) res = p;
		res->x = p->x - len;
		res->y = p->y - len;
		res->z = p->z - len;
	}
	
	static void subtract(nonoPoint *p, nonoPoint *p2, nonoPoint *res = NULL){
		if(res == NULL) res = p;
		res->x = p->x - p2->x;
		res->y = p->y - p2->y;
		res->z = p->z - p2->z;
	}
	
	
	
	static void normal(nonoPoint *p1, nonoPoint *p2,nonoPoint *p3, nonoPoint *res){
		
/*		var v:NonoPoint = new NonoPoint();
		var v1:NonoPoint = p2.copy();
		var v2:NonoPoint = p3.copy();
		v1.subtractMe(p1);
		v2.subtractMe(p1);
		v->x = (v1.y*v2.z) - (v1.z*v2.y);
		v->y = (v1.z*v2.x) - (v1.x*v2.z);
		v->z = (v1.x*v2.y) - (v1.y*v2.x);
*/		
		double v1x = p2->x-p1->x;
		double v1y = p2->y-p1->y;
		double v1z = p2->z-p1->z;
		double v2x = p3->x-p1->x;
		double v2y = p3->y-p1->y;
		double v2z = p3->z-p1->z;
		
		res->x = (v1y*v2z) - (v1z*v2y);
		res->y = (v1z*v2x) - (v1x*v2z);
		res->z = (v1x*v2y) - (v1y*v2x);
		
		double d = nonoGeom::getLength(res);
		nonoGeom::multiply(res,1.0/d);
	}	
	
	static void center(nonoPoint **points, int amount, nonoPoint *res){
		for (int i = 0; i < amount; i++) {
			res->x += points[i]->x;
			res->y += points[i]->y;
			res->z += points[i]->z;
		}
		nonoGeom::divide(res,amount);
	}
	
	
	// http://local.wasp.uwa.edu.au/~pbourke/geometry/planeeq/
	// http://local.wasp.uwa.edu.au/~pbourke/geometry/planeline/
	
	static bool intersectionLinePlane(nonoPoint *l1,nonoPoint *l2, nonoPoint *p1, nonoPoint *p2, nonoPoint *p3, nonoPoint *intersection){
		
		double a = p1->y*(p2->z-p3->z) + p2->y*(p3->z-p1->z) + p3->y*(p1->z-p2->z);
		double b = p1->z*(p2->x-p3->x) + p2->z*(p3->x-p1->x) + p3->z*(p1->x-p2->x);
		double c = p1->x*(p2->y-p3->y) + p2->x*(p3->y-p1->y) + p3->x*(p1->y-p2->y);
		double d = -(  (p1->x*((p2->y*p3->z) - (p3->y*p2->z))) +
			  p2->x*(p3->y*p1->z - p1->y*p3->z) +
			  p3->x*(p1->y*p2->z - p2->y*p1->z) );
		
		double l = sqrt(a*a + b*b + c*c);
		a /= l;
		b /= l;
		c /= l;
		d /= l;
		
		printf(" | \n");
		printf(" |   p1 = [x=%f,y=%f,z=%f]\n",p1->x,p1->y,p1->z);
		printf(" |   p2 = [x=%f,y=%f,z=%f]\n",p2->x,p2->y,p2->z);
		printf(" |   p3 = [x=%f,y=%f,z=%f]\n",p3->x,p3->y,p3->z);
		printf(" | \n");
		printf(" |   l1 = [x=%f,y=%f,z=%f]\n",l1->x,l1->y,l1->z);
		printf(" |   l2 = [x=%f,y=%f,z=%f]\n",l2->x,l2->y,l2->z);		
		printf(" | \n");
		printf(" |   plane = [a=%f,b=%f,c=%f,d=%f]\n",a,b,c,d);
		
		double val1 = a*l1->x + b * l1->y + c * l1->z + d;
		double val2 = a * (l1->x-l2->x) + b * (l1->y-l2->y) + c * (l1->z - l2->z);
		if(val1 == 0) return false;
		if(val2 == 0) return false;
		double u = val1/val2;

		printf(" |   u = %f\n",u);
		reset(intersection);
		nonoGeom::subtract(l2,l1,intersection);
		printf(" |   subtracted = [x=%f,y=%f,z=%f]\n",intersection->x,intersection->y,intersection->z);
		nonoGeom::multiply(intersection,u,intersection);
		printf(" |   multiplied = [x=%f,y=%f,z=%f]\n",intersection->x,intersection->y,intersection->z);
		nonoGeom::add(intersection,l1);
		
		printf(" |   intersection = [x=%f,y=%f,z=%f]\n",intersection->x,intersection->y,intersection->z);
		
		return true;
	}
	
	
	static bool intersectionLineLine(nonoPoint *p1a, nonoPoint *p1b, nonoPoint *p2a, nonoPoint *p2b, nonoPoint *pRes){
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
	
	static int pointInsidePolygon(int nPoints, nonoPoint *points, nonoPoint *pointToTest)
	{
		int i, j, c = 0;
		for (i = 0, j = nPoints-1; i < nPoints; j = i++) {
			if ( ((points[i].y > pointToTest->y) != (points[i].y > pointToTest->y)) &&
				(pointToTest->x < (points[j].x-points[i].x) * (pointToTest->y - points[i].y) / (points[j].y - points[i].y) + points[i].x) )
				c = !c;
		}
		return c;
	}	
	
	
	
	static bool pointInsidePolygonNew(nonoPoint *polygon, int N, nonoPoint *p){
	
		int counter = 0;
		int i;
		double xinters;
		nonoPoint p1,p2;
		
		p1 = polygon[0];
		for (i=1;i<=N;i++) {
			p2 = polygon[i % N];
			
//			printf(" [ %f , %f ]	 [ %f , %f ]\n",p1.x,p1.y,p2.x,p2.y);
			
			if (p->y > NONO_MIN(p1.y,p2.y)) {
				if (p->y <= NONO_MAX(p1.y,p2.y)) {
					if (p->x <= NONO_MAX(p1.x,p2.x)) {
						if (p1.y != p2.y) {
							xinters = (p->y-p1.y)*(p2.x-p1.x)/(p2.y-p1.y)+p1.x;
							if (p1.x == p2.x || p->x <= xinters)
								counter++;
						}
					}
				}
			}
			p1 = p2;
		}
		
		if (counter % 2 == 0)
			return(false);
		else
			return(true);
		
	}
	
	static double polygonArea(nonoPoint *polygon,int N)
	{
		int i,j;
		double area = 0;
		
		for (i=0;i<N;i++) {
			j = (i + 1) % N;
			area += polygon[i].x * polygon[j].y;
			area -= polygon[i].y * polygon[j].x;
		}
		
		area /= 2;
		return(area < 0 ? -area : area);
	}
	
	static double polygonArea(nonoPoint **polygon,int N)
	{
		int i,j;
		double area = 0;
		
		for (i=0;i<N;i++) {
			j = (i + 1) % N;
			area += polygon[i]->x * polygon[j]->y;
			area -= polygon[i]->y * polygon[j]->x;
		}
		
		area /= 2;
		return(area < 0 ? -area : area);
	}
	
	
//	int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy)
//	{
//		int i, j, c = 0;
//		for (i = 0, j = nvert-1; i < nvert; j = i++) {
//			if ( ((verty[i]>testy) != (verty[j]>testy)) &&
//				(testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
//				c = !c;
//		}
//		return c;
//	}	
	
//function SameSide(p1,p2, a,b)
	static bool pointOnSameSide(nonoPoint* p1, nonoPoint* p2, nonoPoint* a, nonoPoint* b){
		return true;
//		nonoPoint cp1 = nonoPoint::cross();
		
//		nonoPoint cp1 = CrossProduct(b-a, p1-a);
//		nonoPoint cp2 = CrossProduct(b-a, p2-a);
//		dot()
//		if DotProduct(cp1, cp2) >= 0 then return true
//		else return false
//			ofPoint;
	}
			
	static bool pointInTriangle(nonoPoint* pCheck, nonoPoint* p1, nonoPoint* p2, nonoPoint* p3){
		if ( pointOnSameSide(pCheck,p1, p2,p3) && pointOnSameSide(pCheck,p2, p1,p3) && pointOnSameSide(pCheck,p3, p1,p2)) return true;
		return false;
	}
	
	
	/*
	 Return whether a polygon in 2D is concave or convex
	 return 0 for incomputables eg: colinear points
	 CONVEX == 1
	 CONCAVE == -1
	 It is assumed that the polygon is simple
	 (does not intersect itself or have holes)
	 
	 http://paulbourke.net/geometry/clockwise/index.html#clockwise
	 
	 */
	
	static int _isPolygonConvex(nonoPoint *p,int n)
	{
		int i,j,k;
		int flag = 0;
		double z;
		
		if (n < 3)
			return(0);
		
		for (i=0;i<n;i++) {
			j = (i + 1) % n;
			k = (i + 2) % n;
			z  = (p[j].x - p[i].x) * (p[k].y - p[j].y);
			z -= (p[j].y - p[i].y) * (p[k].x - p[j].x);
			if (z < 0)
				flag |= 1;
			else if (z > 0)
				flag |= 2;
			if (flag == 3)
				return(-1);
		}
		if (flag != 0)
			return(1);
		else
			return(0);
	}
	
	// isLeft(): tests if a point is Left|On|Right of an infinite line.
	//    Input:  three points P0, P1, and P2
	//    Return: >0 for P2 left of the line through P0 and P1
	//            =0 for P2 on the line
	//            <0 for P2 right of the line
	//    See: the January 2001 Algorithm on Area of Triangles
	static inline float
	isLeft( nonoPoint P0, nonoPoint P1, nonoPoint P2 )
	{
		return (P1.x - P0.x)*(P2.y - P0.y) - (P2.x - P0.x)*(P1.y - P0.y);
	}
	//===================================================================
	
	
	// chainHull_2D(): Andrew's monotone chain 2D convex hull algorithm
	//     Input:  P[] = an array of 2D points
	//                   presorted by increasing x- and y-coordinates
	//             n = the number of points in P[]
	//     Output: H[] = an array of the convex hull vertices (max is n)
	//     Return: the number of points in H[]
	// http://softsurfer.com/Archive/algorithm_0109/algorithm_0109.htm
	
	static int
	chainHull_2D( nonoPoint* P, int n, nonoPoint* H )
	{
		// the output array H[] will be used as the stack
		int    bot=0, top=(-1);  // indices for bottom and top of the stack
		int    i;                // array scan index
		
		// Get the indices of points with min x-coord and min|max y-coord
		int minmin = 0, minmax;
		float xmin = P[0].x;
		for (i=1; i<n; i++)
			if (P[i].x != xmin) break;
		minmax = i-1;
		if (minmax == n-1) {       // degenerate case: all x-coords == xmin
			H[++top] = P[minmin];
			if (P[minmax].y != P[minmin].y) // a nontrivial segment
				H[++top] = P[minmax];
			H[++top] = P[minmin];           // add polygon endpoint
			return top+1;
		}
		
		// Get the indices of points with max x-coord and min|max y-coord
		int maxmin, maxmax = n-1;
		float xmax = P[n-1].x;
		for (i=n-2; i>=0; i--)
			if (P[i].x != xmax) break;
		maxmin = i+1;
		
		// Compute the lower hull on the stack H
		H[++top] = P[minmin];      // push minmin point onto stack
		i = minmax;
		while (++i <= maxmin)
		{
			// the lower line joins P[minmin] with P[maxmin]
			if (isLeft( P[minmin], P[maxmin], P[i]) >= 0 && i < maxmin)
				continue;          // ignore P[i] above or on the lower line
			
			while (top > 0)        // there are at least 2 points on the stack
			{
				// test if P[i] is left of the line at the stack top
				if (isLeft( H[top-1], H[top], P[i]) > 0)
					break;         // P[i] is a new hull vertex
				else
					top--;         // pop top point off stack
			}
			H[++top] = P[i];       // push P[i] onto stack
		}
		
		// Next, compute the upper hull on the stack H above the bottom hull
		if (maxmax != maxmin)      // if distinct xmax points
			H[++top] = P[maxmax];  // push maxmax point onto stack
		bot = top;                 // the bottom point of the upper hull stack
		i = maxmin;
		while (--i >= minmax)
		{
			// the upper line joins P[maxmax] with P[minmax]
			if (isLeft( P[maxmax], P[minmax], P[i]) >= 0 && i > minmax)
				continue;          // ignore P[i] below or on the upper line
			
			while (top > bot)    // at least 2 points on the upper stack
			{
				// test if P[i] is left of the line at the stack top
				if (isLeft( H[top-1], H[top], P[i]) > 0)
					break;         // P[i] is a new hull vertex
				else
					top--;         // pop top point off stack
			}
			H[++top] = P[i];       // push P[i] onto stack
		}
		if (minmax != minmin)
			H[++top] = P[minmin];  // push joining endpoint onto stack
		
		return top+1;
	}
	
	
//	static int
//	giftWrap( nonoPoint* P, int n, nonoPoint* H )
//	{
//		p = left-most, bottom point
//		current = p
//		do {
//			// Sentinel
//			leftmost = 1
//			for(int i=1; i<n;i++){// i from 2 to N
//				// If line current->leftmost is to the left of current->point[i]
//				if ( left( current, point[i], leftmost ) ){
//					leftmost = i
//				}
//			}
//			// Output and move on.
//			current = point[leftmost]
//		} while ( current != p );
//	}
	
	static void calculateCentroid(vector<nonoPoint*> points, nonoPoint* centroid){
		int n = points.size();
		centroid->x = 0;
		centroid->y = 0;		
		centroid->z = 0;
		for(int i=0;i<n;i++){
			centroid->x += points[i]->x;
			centroid->y += points[i]->y;
			centroid->z += points[i]->z;			
		}
		centroid->x /= (double)n;
		centroid->y /= (double)n;
		centroid->z /= (double)n;
	}
	
	static void calculateCentroid(nonoPoint* points, int n, nonoPoint* centroid){
		centroid->x = 0;
		centroid->y = 0;		
		centroid->z = 0;
		for(int i=0;i<n;i++){
			centroid->x += points[i].x;
			centroid->y += points[i].y;
			centroid->z += points[i].z;			
		}
		centroid->x /= (double)n;
		centroid->y /= (double)n;
		centroid->z /= (double)n;
	}
	
	static void calculateCentroid(nonoPoint** points, int n, nonoPoint* centroid){
		centroid->x = 0;
		centroid->y = 0;		
		centroid->z = 0;
		for(int i=0;i<n;i++){
			centroid->x += points[i]->x;
			centroid->y += points[i]->y;
			centroid->z += points[i]->z;			
		}
		centroid->x /= (double)n;
		centroid->y /= (double)n;
		centroid->z /= (double)n;
	}	
	
	static int findLeftTopPoint(nonoPoint* points, int n){
		nonoPoint centroid;
		int cp = 0;
		float closest = 1000;
		calculateCentroid(points,n,&centroid);
		for(int i=0;i<n;i++){
			float d = nonoMath::angle(&points[i],&centroid);
			if(d < 45) d += 360;
			d = d - 315.0f;
			if(d<0) d *= -1;
			if(d < closest){
				closest = d;
				cp = i;
			}
		}
		return cp;
	}
//
//	static void sortPointsClockwiseFromLeftTop(vector<nonoPoint> points){
//		
//	}
	
	static void sortPointsClockwiseFromLeftTop(vector<nonoPoint*> points, vector<nonoPoint*>* pointsOut){
		nonoPoint centroid;
		int cp = 0;
		float closest = 1000;
		calculateCentroid(points,&centroid);
		vector<PointAnglePair> pointsToSort;
		for(int i=0;i<points.size();i++){
			float d = nonoMath::angle(points[i],&centroid);
			PointAnglePair pair;
			pair.p = points[i];
			pair.angle = d;
			pointsToSort.push_back(pair);
			if(d < 45) d += 360;
			d = d - 315.0f;
			if(d<0) d *= -1;
			if(d < closest){
				closest = d;
				cp = i;
			}
		}
		
		sort(pointsToSort.begin(),pointsToSort.end(),sortPointAnglePair);
		pointsOut->clear();
		for(int i=0;i<pointsToSort.size();i++){
			pointsOut->push_back(pointsToSort[i].p);
//			printf(" [ %f , %f ]	",pointsOut->at(i)->x,pointsOut->at(i)->y);
		}
//		printf("\n");
		pointsToSort.clear();
		
	}
	
	static void __sortPointsClockwiseFromLeftTop(nonoPoint** points, int n, nonoPoint** pointsOut){
		nonoPoint centroid;
		int cp = 0;
		float closest = 1000;
		calculateCentroid(points,n,&centroid);
//		float angles[n];
		vector<PointAnglePair> pointsToSort;
		printf("|| sortPointsClockwiseFromLeftTop     ");
		for(int i=0;i<n;i++){
			float d = nonoMath::angle(points[i],&centroid);
//			angles[i] = d;
			PointAnglePair pair;
			pair.p = points[i];
			pair.angle = d;
			pointsToSort.push_back(pair);
			printf(" [ %f , %f ]	",pair.p->x,pair.p->y);
			if(d < 45) d += 360;
			d = d - 315.0f;
			if(d<0) d *= -1;
			if(d < closest){
				closest = d;
				cp = i;
			}
		}
		printf("\n");
		
		sort(pointsToSort.begin(),pointsToSort.end(),sortPointAnglePair);
		printf("|| sortPointsClockwiseFromLeftTop [%i] ",pointsToSort.size());
		for(int i=0;i<pointsToSort.size();i++){
			pointsOut[i]->x = pointsToSort[i].p->x;
			pointsOut[i]->y = pointsToSort[i].p->y;			
			printf(" [ %f , %f ]	",pointsToSort[i].p->x,pointsToSort[i].p->y);
		}
		printf("\n");
		pointsToSort.clear();		
	}
	
	static void sortPointsClockwiseFromLeftTop(nonoPoint* points, int n, nonoPoint* pointsOut){
		nonoPoint centroid;
		int cp = 0;
		float closest = 1000;
		calculateCentroid(points,n,&centroid);
		vector<PointAnglePair> pointsToSort;
		for(int i=0;i<n;i++){
			float d = nonoMath::angle(&points[i],&centroid);
			PointAnglePair pair;
			pair.p = &points[i];
			pair.angle = d;
			pointsToSort.push_back(pair);
			if(d < 45) d += 360;
			d = d - 315.0f;
			if(d<0) d *= -1;
			if(d < closest){
				closest = d;
				cp = i;
			}
		}
		
		sort(pointsToSort.begin(),pointsToSort.end(),sortPointAnglePair);
		for(int i=0;i<n;i++){
			pointsOut[i].x = pointsToSort[i].p->x;
			pointsOut[i].y = pointsToSort[i].p->y;			
		}
		pointsToSort.clear();		
	}
	
	
	
private:
	
	struct PointAnglePair {
		nonoPoint* p;
		float angle;
		int id;
	};
	
	static bool sortPointAnglePair (PointAnglePair p1, PointAnglePair p2) { 
		return (p1.angle > p2.angle); 
	}
	
};

#endif _NONOGEOM