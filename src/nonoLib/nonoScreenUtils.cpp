/*
 *  screenUtils.cpp
 *  futbolin09cpp
 *
 *  Created by Vincent R. on 25.07.09.
 *  Copyright 2009 www.say-nono.com. All rights reserved.
 *
 */

#include "nonoScreenUtils.h"

nonoScreenUtils::nonoScreenUtils(){
	maxScreens = 10;
	currentScreenID = 0;
	screensFound = 0;
	isActive = false;
	updateList();
}

nonoScreenUtils::~nonoScreenUtils(){
	
}

void nonoScreenUtils::setActive(bool b){
	isActive = b;
}

void nonoScreenUtils::updateList(){
	
	#if defined( TARGET_OSX )
	CGDirectDisplayID displays[maxScreens];
	CGDisplayCount numDisplays;
	CGDisplayErr err;
	
	err = CGGetActiveDisplayList(
								 maxScreens,
								 displays,
								 &numDisplays);
	if ( err != 0 )
	{
		printf("Cannot get displays (%d)\n", err);
	}
	
	screensFound = (int) numDisplays;
	for(int i=0;i<screensFound;i++){
		CGRect rect = CGDisplayBounds(displays[i]);
		CGSize size = rect.size;
		CGPoint origin = rect.origin;
		screens[i].x = (int)origin.x;
		screens[i].y = (int)origin.y;
		screens[i].w = (int)size.width;
		screens[i].h = (int)size.height;
	}
	
	#else if defined( TARGET_WIN32 )
	
	#endif
	
	
	printf("Displays found (%i)\n",screensFound);
	printf("-------------------------------\n");
	for(int i=0;i<screensFound;i++){
		printf("Display #%i		w:%i	h:%i	x:%i	y:%i\n",i,screens[i].w,screens[i].h,screens[i].x,screens[i].y);
	}
	printf("\n");
}

ScreenRect* nonoScreenUtils::nextDisplay(){
	selectDisplay(currentScreenID+1);
}

ScreenRect* nonoScreenUtils::prevDisplay(){
	selectDisplay(currentScreenID-1);
}


ScreenRect* nonoScreenUtils::selectDisplay(int id){
	
	if(screensFound==0){
		printf("ERROR: No display found for swapping.");
		return NULL;
	}
	currentScreenID = id + screensFound;
	currentScreenID %= screensFound;
	if(isActive) ofSetWindowPosition(screens[currentScreenID].x, screens[currentScreenID].y);
	return &screens[currentScreenID];
}

int nonoScreenUtils::getCurrentScreenID(){
	return currentScreenID;
}
