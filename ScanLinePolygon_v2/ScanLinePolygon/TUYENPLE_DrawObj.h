#pragma once

#include "TUYENPLE_RectD.h"
#include <vector>

#define CPGEOM_NA				0
#define CPGEOM_TEXT				1
#define CPGEOM_LINE				2
#define CPGEOM_ARC				3
#define CPGEOM_POINT			4
#define CPGEOM_CONTAINER		5  
#define CPGEOM_POLYGON			6 

#define DIRECTION_CLOCKWISE          1
#define DIRECTION_COUNTER_CLOCKWISE -1


class DrawObj
{

private:
	RectD m_Rect;

public:

public:
	DrawObj(void);
	virtual ~DrawObj();

	virtual DrawObj* clone() const = 0 ;
	virtual void SetRect() 
	{
		m_Rect = RectD(0.0, 0.0, 0.0, 0.0);
	};

	virtual int GetType() {
		return CPGEOM_NA;
	}

	RectD GetRect() {
		return m_Rect;
	}

	void SetRect(double x1, double y1, double x2, double y2) {
		m_Rect = RectD(x1, y1, x2, y2); // left, top, right, bottom
	}

};

typedef std::vector< DrawObj* > DrawObjArray;