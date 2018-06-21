#include <math.h>

#define PI	3.1415926535897932384626433832795f

#define DEG2RAD(x)	( (x) * PI / 180. )
#define RAD2DEG(x)	( (x) / PI * 180. )

#define TOLERANCE_DEFAULT			(0.0001	    	)
#define TOLERANCE_CM				(0.0001			)
#define TOLERANCE_MM				(0.001			)
#define TOLERANCE_MICRO				(0.1				)
#define TOLERANCE_INCH				(TOLERANCE_MM / 25.4		)
#define TOLERANCE_MIL				(TOLERANCE_INCH * 1000		)
#define MATH_TOLERANCE				(0.0001		)

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) > (b) ? (b) : (a))

inline double distance( const double& x1, const double& y1,
					   const double& x2, const double& y2 )
{
	return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) ) ;
}

inline double normalize_angle( const double& dAngleInDegree )
{
	double dNormalized_Angle = dAngleInDegree ;

	if( dNormalized_Angle <= 0.0 )
		dNormalized_Angle += 360.0f ;

	dNormalized_Angle = fmod( dNormalized_Angle, 360.0 ) ;
	return dNormalized_Angle ;
}

// determine equality
inline bool NearlyEquals(double x, double y, double r)
{
	return ( (x - y) * (x - y) ) < r*r ;
}
inline bool NearlyEquals(double x1, double y1, double x2, double y2, double r)
{
	return ( ( (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) ) <  r*r ) ;
}