#if !defined( _MathMacroFn_h_ )
#define _MathMacroFn_h_
/********************************************************************
	created:	2004/06/07
	created:	7:6:2004   14:29
	filename: 	D:\PCBX-WORKING\COMMON\MATH\MathMacroFn.h
	file path:	D:\PCBX-WORKING\COMMON\MATH
	file base:	MathMacroFn
	file ext:	h

	purpose:	simple math related function
							(Macro Functions)
*********************************************************************/

#include <math.h>

#include "MathConstant.h"

// unit convert
#define DEG2RAD(x)	( (x) * MC_PI / 180. )
#define RAD2DEG(x)	( (x) / MC_PI * 180. )

#define SMALL_VAL	0.0001//0.0001	// use for PDB classes( operator == overloading )

// calc area
#define CIRCLE_AREA( Radius )				( MC_PI * Radius * Radius )
#define PIE_AREA( Radius, AngleInRad )		( ( 0.5 * Radius * Radius * AngleInRad ) )
#define SQUARE_AREA( Width )				(  Width * Width )

// open range check 1
//		left_val < val < right_val, ( left_val, right_val )
#ifndef IS_VALUE_IN_RANGE
	#define IS_VALUE_IN_RANGE( left_val, val, right_val )			( val > left_val && val < right_val )
#endif

#define IS_VALUE_OUT_OF_RANGE( left_val, val, right_val )		( val < left_val || val > right_val )



// two point distance
//#define DISTANCE( x1, y1, x2, y2 )			( sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) ) ; )



// determine equality
inline bool NearlyEquals(double x, double y, double r)
{
	return ( (x - y) * (x - y) ) < r*r ;
}

inline bool NearlyEquals(double x1, double y1, double x2, double y2, double r)
{
	return ( ( (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) ) <  r*r ) ;
}

inline bool NearlyEquals(double x1, double y1, double z1, double x2, double y2, double z2, double r)
{
	return ( ( (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2) ) <  r*r ) ;
}



// coordinate transformation : float coord => int coord
//		ex> logical coordinate => device coordinate
//typedef void (*CoordTransformF2I)( float x1, float y1, int& x2, int& y2 ) ;


#endif // !defined( _MathMacroFn_h_ )