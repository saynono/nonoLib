/*
 *  screenUtils.h
 *  futbolin09cpp
 *
 *  Created by Vincent R. on 25.07.09.
 *  Copyright 2009 www.say-nono.com. All rights reserved.
 *
 */

#include "ofMain.h"
#include "ofAppRunner.h"

struct ScreenRect{
	int w;
	int h;
	int x;
	int y;
};

class nonoScreenUtils
{

	public:

		nonoScreenUtils();
		virtual ~nonoScreenUtils();
	
	
		void setActive(bool b);
		void updateList();
		ScreenRect* nextDisplay();
		ScreenRect* prevDisplay();
		ScreenRect* selectDisplay(int id);
		int getCurrentScreenID();

	private:
		int			maxScreens;
		ScreenRect	screens[10];
		int			currentScreenID;
		int			screensFound;
		bool		isActive;
};
