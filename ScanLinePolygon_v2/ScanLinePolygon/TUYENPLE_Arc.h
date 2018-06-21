#pragma once

#include "TUYENPLE_DrawObj.h"
#include "TUYENPLE_Line.h"
#include "TUYENPLE_PointD.h"

class Arc : public DrawObj 
{

private:

	PointD m_StartPoint;
	PointD m_EndPoint;
	PointD m_Center;

	double m_Radius;
	double m_StartAngle;
	double m_EndAngle;
	double m_Thickness;
	int m_Direction;

public:
	Arc(void);
	~Arc(void);
	Arc(const Arc& src);
	virtual DrawObj* clone() const		{	return new Arc(*this) ;	}
	LineArray* GetLineSegArray(const double& dAngleOfSlice = 10.0);

	void SetRect();

	virtual int GetType() 
	{
		return CPGEOM_ARC;
	}

	PointD GetStartPoint() {
		return m_StartPoint;
	}

	PointD GetEndPoint() {
		return m_EndPoint;
	}

	PointD GetCenter() {
		return m_Center;
	}

	void SetStartPoint(PointD p) {
		m_StartPoint = p;
	}

	void SetEndPoint(PointD p) {
		m_EndPoint = p;
	}

	void SetCenter(PointD p) {
		m_Center = p;
	}

	void SetAngle(double s, double e) {
		m_StartAngle = s;
		m_EndAngle = e;
	}

	void SetThickness(double d) {
		m_Thickness = d;
	}

	void SetRadius(double r) {
		m_Radius = r;
	}

	void SetStartAngle(double s) {
		m_StartAngle = s;
	}

	void SetEndAngle(double e) {
		m_EndAngle = e;
	}

	double GetThickness() {
		return m_Thickness;
	}

	double GetStartAngle() {
		return m_StartAngle;
	}

	double GetEndAngle() {
		return m_EndAngle;
	}

	double GetRadius() {
		return m_Radius;
	}

	void setDirection(int dir) {
		m_Direction = dir;
	}

	int GetDirection() {
		return m_Direction;
	}

	bool IsCCW() {
		return m_Direction == DIRECTION_COUNTER_CLOCKWISE;
	}

	bool Is_Pt_On_Arc(PointD ptPoint);
	LineArray* GetLineSegArrayForCircle(const double& dAngleOfSlice);
	double GetSweepAngle();
};