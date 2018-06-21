#include "TUYENPLE_PointD.h"
#include "TUYENPLE_Math.h"

PointD::PointD(double x, double y) {
	m_X = x;
	m_Y = y;
}

PointD::~PointD(void) {
}

PointD& PointD::operator=( const PointD& rhs )
{
	if( this == &rhs )
		return *this ;

	m_X = rhs.m_X ;
	m_Y = rhs.m_Y ;

	return *this ;
}

bool PointD::operator==(const PointD& rhs) const
{
	if( NearlyEquals(m_X, m_Y, rhs.m_X, rhs.m_Y, TOLERANCE_MM * 10.0 ) )
		return true ;

	return false ;
}

bool PointD::operator!=(const PointD& rhs) const
{
	if( !( *this == rhs ) )
		return true ;

	return false ;	
}

PointD PointD::operator+(const PointD& rhs) const
{
	PointD ptRes = *this ;

	ptRes.m_X += rhs.m_X ;
	ptRes.m_Y += rhs.m_Y ;

	return ptRes ;
}

PointD& PointD::operator+=(const PointD& rhs)
{
	m_X += rhs.m_X ;	
	m_Y += rhs.m_Y ;

	return *this ;
}

PointD PointD::operator-(const PointD& rhs) const
{
	PointD ptRes = *this ;

	ptRes.m_X -= rhs.m_X ;
	ptRes.m_Y -= rhs.m_Y ;

	return ptRes ;
}

PointD& PointD::operator-=(const PointD& rhs)
{
	m_X -= rhs.m_X ;	
	m_Y -= rhs.m_Y ;

	return *this ;
}

PointD PointD::operator/( const double& dDenominator ) 
{
	PointD Result = *this ;

	Result.m_X /= dDenominator ;
	Result.m_Y /= dDenominator ;

	return Result ;
}

double PointD::GetDistanceWith(PointD Pt)  const
{
	return distance( m_X, m_Y, Pt.GetX(), Pt.GetY() ) ;
}