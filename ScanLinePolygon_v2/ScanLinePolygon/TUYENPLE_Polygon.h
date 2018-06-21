#pragma once
#include "TUYENPLE_DrawObj.h"
#include "TUYENPLE_PointD.h"
#include "gf_types.h"

class Polygon;

typedef std::vector< Polygon* > PolygonArray;

class Polygon : public DrawObj {

private:
	DrawObjArray* m_pBoundary;
	PolygonArray* m_pCuttingArea;

public:
	Polygon(void);
	~Polygon(void);

	Polygon(const Polygon& src);
	virtual DrawObj* clone() const		{	return new Polygon(*this) ;	}

	bool IsInner(PointD Pt, const bool& bPtOnBoundaryIsInner = true  );
	int IsInner(PointD Pt, const DrawObjArray* const pBoundaryObjs);
	void ConvertToPointArray(PointDArray &, const DrawObjArray *, double dDivision = 10, int UnitType = 0/*0:degree, 1:mm */);
	gf_polygon * Make_GFPolygon(double dDivision = 10, int UnitType = 0/*0:degree, 1:mm */);

	void SetRect();

	virtual int GetType() {
		return CPGEOM_POLYGON;
	}

	DrawObjArray* getBoundary() {
		return m_pBoundary;
	}

	PolygonArray* getCuttingArea() {
		return m_pCuttingArea;
	}

	void setBoundary(DrawObjArray* b) {
		m_pBoundary = b;
	}

	void setCuttingArea(PolygonArray *p) {
		m_pCuttingArea = p;
	}
};
