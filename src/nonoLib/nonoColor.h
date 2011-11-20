/*
 *  nonoColor.h
 *  SVSechzehn_InfoVisuals
 *
 *  Created by say nono on 31.08.10.
 *  Copyright 2010 www.say-nono.com. All rights reserved.
 *
 */


#pragma once

#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif

#ifndef M_HEX_VALUES
#define M_HEX_VALUES = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};
#endif

#include <math.h>
#include <algorithm>
#include "nonoMath.h"
#include "ofMain.h"

/*
 template <class T> const T& max ( const T& a, const T& b ) {
 return (b<a)?a:b;     // or: return comp(b,a)?a:b; for the comp version
 }
 
 template <class T> const T& min ( const T& a, const T& b ) {
 return (b>a)?a:b;     // or: return comp(b,a)?a:b; for the comp version
 }
 */

using namespace std;

class nonoColor{
	
public:
	
	
//	static double myTan(double val){
//		return tan(toRadians(val));
//	}
	
	static string toHexString(int c) {
		char buffer [50];
		sprintf(buffer, "0x%06X",c);
		string result(buffer);
		return result;
	}
	
//	static int HextoDec(char *hex)
//	{
//		if (*hex==0) return 0;
//		return  HextoDec(hex-1)*16 +  xtod(*hex) ; 
//	}
//
//	static int xstrtoi(char *hex)      // hex string to integer
//	{
//		return HextoDec(hex+strlen(hex)-1);
//	}	
	
	
	static int getAlpha(int c) {
		return (c >> 24) & 0x000000FF;
	}
	
	static int getRed(int c) {
		return (c >> 16) & 0x000000FF;
	}
	
	static int getGreen(int c) {
		return (c >> 8) & 0x000000FF;
	}
	
	static int getBlue(int c) {
		return c & 0x000000FF;
	}
	
	static int setAlpha(int c, int clr=0) {
		return ((c & 0x000000FF) << 24) + (clr & 0xffffff);
	}
	
	static int setRed(int c, int clr=0) {
		return ((c & 0x000000FF) << 16) + (clr & 0xff00ffff);
	}
	
	static int setGreen(int c,int clr=0) {
		return ((c & 0x000000FF) << 8) + (clr & 0xffff00ff);
	}
	
	static int setBlue(int c,int clr=0) {
		return (c & 0x000000FF) + (clr & 0xffffff00);
	}
	
	static int setColors(int r,int g,int b,int a = 255) {
		int c = 0;
		c += setRed(r);
		c += setGreen(g);
		c += setBlue(b);
		c += setAlpha(a);
		return c;
	} 
	
	static int setGrey(int grey) {
		int c = 0;
		c += setRed(grey);
		c += setGreen(grey);
		c += setBlue(grey);
		return c;
	}	
	
	static int mixColors(int c1,int c2,float mix) {
		float antimix = 1 - mix;
		int r = setRed((getRed(c1) * antimix) + (getRed(c2) * mix));
		int g = setGreen((getGreen(c1) * antimix) + (getGreen(c2) * mix));
		int b = setBlue((getBlue(c1) * antimix) + (getBlue(c2) * mix));
		return (r + g + b);
	}
	
	
	
	static int HSB2RGB(float h = 1,float s = 1,float b = 1) {
		h = int(h) % 360;
		int i = (int)((int)(h / 60.0) % 6);
		float f = h / 60.0 - (int)(h / 60.0);
		float p = b * (1 - s);
		float q = b * (1 - s * f);
		float t = b * (1 - (1 - f) * s);
		
		b*=255;
		t*=255;
		p*=255;
		q*=255;
		
		switch (i) {   
			case 0: 
				return setColors(b, t, p); 
			case 1: 
				return setColors(q, b, p); 
			case 2: 
				return setColors(p, b, t); 
			case 3: 
				return setColors(p, q, b); 
			case 4: 
				return setColors(t, p, b); 
			case 5:
				return setColors(b, p, q); 
		}
		return setColors(0, 0, 0);
	}
	
//	static function RGB2HSB(c:uint):Object{
//		
//		var r:int = getRed(c);
//		var g:int = getGreen(c);
//		var b:int = getBlue(c);
//		
//		r = (r < 0)? 0 : (r>255)? 255: Math.round(r);
//		g = (g < 0)? 0 : (g>255)? 255: Math.round(g);
//		b = (b < 0)? 0 : (b>255)? 255: Math.round(b);
//		var min:Number = Math.min(r, g, b);
//		var max:Number = Math.max(r, g, b);
//		if(max==0){
//			return {h:0, s: 0, b: 0};
//		}else{
//			var s:Number = (max - min)/max * 100;
//		}
//		var bri:Number = max / 255 * 100;
//		var h:Number = _getHue(r, b, g, max, min);
//		return {h:h, s:s, b:bri};
//	}	
	
private:
	

};