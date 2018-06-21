#pragma once

#include "TUYENPLE_DrawObj.h"
#include "TUYENPLE_PointD.h"

class Line;

typedef std::vector< Line* > LineArray;

class Line : public DrawObj {

private:

	PointD m_StartPoint;
	PointD m_EndPoint;
	double m_Thickness;

public:
	Line(void);
	~Line(void);

	Line(const Line& src);
	virtual DrawObj* clone() const		{	return new Line(*this) ;	}

	void SetRect();

	virtual int GetType() {
		return CPGEOM_LINE;
	}

	void SetStartPoint(PointD p) {
		m_StartPoint = p;
	}

	void SetEndPoint(PointD p) {
		m_EndPoint = p;
	}

	PointD GetStartPoint() {
		return m_StartPoint;
	}

	PointD GetEndPoint() {
		return m_EndPoint;
	}

	void SetThickness(double d) {
		m_Thickness = d;
	}

	double GetThickness() {
		return m_Thickness;
	}

	bool Is_Pt_On_LineSeg_Math(const PointD &Pt) const;
	double GetLength() const;


};