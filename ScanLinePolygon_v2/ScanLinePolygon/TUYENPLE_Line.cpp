#include "TUYENPLE_Line.h"
#include "TUYENPLE_Math.h"

Line::Line() {
}

Line::~Line(void) {
}

void  Line::SetRect() {
	double left = min(m_StartPoint.GetX(), m_EndPoint.GetX());
	double top = max(m_StartPoint.GetY(), m_EndPoint.GetY());
	double right = max(m_StartPoint.GetX(), m_EndPoint.GetX());
	double bottom = min(m_StartPoint.GetY(), m_EndPoint.GetY());

	DrawObj::SetRect(left, top, right, bottom);
}

bool Line::Is_Pt_On_LineSeg_Math(const PointD &Pt) const
{
	// l => StartPt to Endpt Distance
	// l1 => ptPtOnLine to StartPt  Distance
	// l2 => ptPtOnLine to Endpt Distance
	// l3 => l1 + l2	

	double l1 = Pt.GetDistanceWith( m_StartPoint ) ;
	double l2 = Pt.GetDistanceWith( m_EndPoint ) ;
	double l3 = l1 + l2 ;

	double l = GetLength() ;

	if( fabs( l3 - l ) < MATH_TOLERANCE )
		return true ;

	return false ;
}

double Line::GetLength() const
{
	return m_StartPoint.GetDistanceWith( m_EndPoint ) ;
}

Line::Line(const Line& src)
{
	*this = src ;	
}

