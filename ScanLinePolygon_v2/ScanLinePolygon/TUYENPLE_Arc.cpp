#include "TUYENPLE_Arc.h"
#include "TUYENPLE_Math.h"
#include <math.h>


Arc::Arc(void) {
}

Arc::~Arc(void) {
}

Arc::Arc(const Arc& src)
{
	*this = src ;	
}

void  Arc::SetRect() {
	PointD StartPt = GetStartPoint() ;
	PointD EndPt = GetEndPoint() ;

	PointD ptTop, ptLeft, ptBtm, ptRight ;
	ptTop = ptLeft = ptBtm = ptRight = m_Center ;
	ptTop.m_Y += m_Radius;
	ptLeft.m_X -= m_Radius;
	ptBtm.m_Y -= m_Radius;
	ptRight.m_X += m_Radius;

	if( m_StartAngle == m_EndAngle )
	{
		DrawObj::SetRect(ptLeft.GetX(), ptTop.GetY(), ptRight.GetX(), ptBtm.GetY());
		return;
	}

	double left,top,right,bottom;
	left = min(m_StartPoint.GetX(), m_EndPoint.GetX());
	top = max(m_StartPoint.GetY(), m_EndPoint.GetY());
	right = max(m_StartPoint.GetX(), m_EndPoint.GetX());
	bottom = min(m_StartPoint.GetY(), m_EndPoint.GetY());

	if( Is_Pt_On_Arc(ptLeft) )
	{
		left = min(left, ptLeft.GetX());
		top = max(top, ptLeft.GetY());
		right = max(right, ptLeft.GetX());
		bottom = min(bottom, ptLeft.GetY());
	}
	if( Is_Pt_On_Arc(ptTop) )
	{
		left = min(left, ptTop.GetX());
		top = max(top, ptTop.GetY());
		right = max(right, ptTop.GetX());
		bottom = min(bottom, ptTop.GetY());
	}
	if( Is_Pt_On_Arc(ptRight) )
	{
		left = min(left, ptRight.GetX());
		top = max(top, ptRight.GetY());
		right = max(right, ptRight.GetX());
		bottom = min(bottom, ptRight.GetY());
	}
	if( Is_Pt_On_Arc(ptBtm) )
	{
		left = min(left, ptBtm.GetX());
		top = max(top, ptBtm.GetY());
		right = max(right, ptBtm.GetX());
		bottom = min(bottom, ptBtm.GetY());
	}

	DrawObj::SetRect(left, top, right, bottom);
	return;
}

bool Arc::Is_Pt_On_Arc(PointD ptPoint)
{
	bool bOnArc = false ;

	double dDistance = distance(ptPoint.GetX(),ptPoint.GetY(),m_Center.GetX(),m_Center.GetY());
	if(fabs(dDistance - m_Radius) >= TOLERANCE_MM)
		return bOnArc;

	//According to "Geometric Tools for Computer Graphics", page 249
	PointD ptStartPt = GetStartPoint();
	PointD ptEndPt = GetEndPoint();
	double dKrossValue;

	PointD vectorOne = PointD(ptPoint.GetX()-ptStartPt.GetX(),ptPoint.GetY()-ptStartPt.GetY());
	PointD vectorTwo = PointD(ptEndPt.GetX()-ptStartPt.GetX(),ptEndPt.GetY()-ptStartPt.GetY());

	dKrossValue = vectorOne.GetX()*vectorTwo.GetY() - vectorTwo.GetX()*vectorOne.GetY();

	if( DIRECTION_COUNTER_CLOCKWISE == m_Direction)
	{
		if(dKrossValue >= 0)
			return true;
		else
			return false;
	}
	else if( DIRECTION_CLOCKWISE == m_Direction)
	{
		if(dKrossValue <= 0)
			return true;
		else
			return false;
	}

	return bOnArc ;	
}

LineArray * Arc::GetLineSegArray(const double& dAngleOfSlice)
{
	LineArray* pLineSegArray = new LineArray;		

	int	   nSegmentCount = 0;
	double dStartAngle,dStartAngleInRad;	
	double dEndAngle;
	double dAngleDistance;
	
	PointD StartPoint;
	PointD EndPoint;
	PointD TempPoint;

	StartPoint = GetStartPoint();	EndPoint = GetEndPoint();
	dStartAngle = normalize_angle(m_StartAngle);	
	dEndAngle = normalize_angle(m_EndAngle);

	Line* pLineSeg = NULL;

	bool bIsCCW = (m_Direction == DIRECTION_COUNTER_CLOCKWISE) ;
	
	if (bIsCCW)
	{
		dAngleDistance = (dEndAngle > dStartAngle) ? (dEndAngle-dStartAngle) : (dEndAngle-dStartAngle+360.0);	
	}
	else
	{
		dAngleDistance = (dStartAngle > dEndAngle) ? (dStartAngle-dEndAngle) : (dStartAngle-dEndAngle+360.0);
	}
	
	
	if(dAngleDistance < dAngleOfSlice)
	{
		Line* pEndLineSeg = new Line;
		pEndLineSeg->SetStartPoint(StartPoint);
		pEndLineSeg->SetEndPoint(EndPoint);

		pLineSegArray->push_back(pEndLineSeg);
	
		return pLineSegArray;
	}

	while(dAngleOfSlice <= dAngleDistance)
	{
		pLineSeg = new Line;						

		pLineSeg->SetStartPoint(StartPoint);

		dStartAngle = (bIsCCW) ? (dStartAngle+dAngleOfSlice) : (dStartAngle-dAngleOfSlice);		
		
		dStartAngle = normalize_angle(dStartAngle);		

		dStartAngleInRad = DEG2RAD(dStartAngle) ;
		TempPoint.SetX( m_Center.GetX() + m_Radius*cos(dStartAngleInRad) ) ;
		TempPoint.SetY( m_Center.GetY() + m_Radius*sin(dStartAngleInRad) ) ;

		pLineSeg->SetEndPoint(TempPoint);
		
		pLineSegArray->push_back(pLineSeg);

		StartPoint = TempPoint;
		dAngleDistance -= dAngleOfSlice;
	}

	if(dAngleDistance != 0.0)
	{
		Line* pEndLineSeg = new Line;
		pEndLineSeg->SetStartPoint(TempPoint);
		pEndLineSeg->SetEndPoint(EndPoint);
		pLineSegArray->push_back(pEndLineSeg);
	}
	
	return pLineSegArray;
}

LineArray* Arc::GetLineSegArrayForCircle(const double& dAngleOfSlice)
{
	LineArray* pLineSegArray = new LineArray;		

	int	   nSegmentCount = 0;
	double dStartAngle,dStartAngleInRad;	
	double dEndAngle;
	double dAngleDistance;

	PointD StartPoint;
	PointD EndPoint;
	PointD	TempPoint;

	StartPoint = GetStartPoint();	EndPoint = GetEndPoint();
	dStartAngle = normalize_angle(m_StartAngle);	
	dEndAngle = normalize_angle(m_EndAngle);

	Line* pLineSeg = NULL;

	bool bIsCCW = (m_Direction == DIRECTION_COUNTER_CLOCKWISE) ;

	dAngleDistance = 360.0;

	while(dAngleOfSlice <= dAngleDistance)
	{
		pLineSeg = new Line;						

		pLineSeg->SetStartPoint(StartPoint);

		dStartAngle = (bIsCCW) ? (dStartAngle+dAngleOfSlice) : (dStartAngle-dAngleOfSlice);		

		dStartAngle = normalize_angle(dStartAngle);		

		dStartAngleInRad = DEG2RAD(dStartAngle) ;
		TempPoint.SetX( m_Center.GetX() + m_Radius*cos(dStartAngleInRad) ) ;
		TempPoint.SetY( m_Center.GetY() + m_Radius*sin(dStartAngleInRad) ) ;

		pLineSeg->SetEndPoint(TempPoint);

		pLineSegArray->push_back(pLineSeg);

		StartPoint = TempPoint;
		dAngleDistance -= dAngleOfSlice;
	}

	if(dAngleDistance != 0.0)
	{
		Line* pEndLineSeg = new Line;
		pEndLineSeg->SetStartPoint(TempPoint);
		pEndLineSeg->SetEndPoint(EndPoint);
		pLineSegArray->push_back(pEndLineSeg);
	}

	return pLineSegArray;
}

double Arc::GetSweepAngle()
{
	double dSweepAngle = 0.0 ;

	// calculate sweep angle
	//		refer to "determine arc sweep angle.xls" in $UDE-WORKING/Doc/
	double dStartAngle, dEndAngle ;
	dStartAngle = GetStartAngle() ;
	dEndAngle	= GetEndAngle() ;

	if( dStartAngle == dEndAngle )
		return 360.0 ;

	double dDiffAngle = dStartAngle - dEndAngle ;
	bool bNegativeDiff = dDiffAngle < 0  ;
	bool bCW = ( m_Direction == DIRECTION_CLOCKWISE ) ;
	dSweepAngle = ( bNegativeDiff ^ bCW ) ? fabs( dDiffAngle ) : ( 360.0 - fabs( dDiffAngle ) ) ;

	return dSweepAngle ;
}