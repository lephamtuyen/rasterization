#include "TUYENPLE_RectD.h"

RectD::RectD( double x1 /*= 0.*/, double y1 /*= 0.*/, double x2 /*= 0.*/, double y2 /*= 0.*/ )
:	m_left(x1),
m_top(y1),
m_right(x2),
m_bottom(y2) {

}

RectD::~RectD(void) {
}

RectD& RectD::operator=( const RectD& rhs ) {
	if( this == &rhs )
		return *this ;

	m_left = rhs.m_left;
	m_top = rhs.m_top;
	m_right = rhs.m_right;
	m_bottom = rhs.m_bottom;

	return *this ;
}