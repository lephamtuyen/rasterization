// -----------------------------------------------------------------------------
// Point-To-Line Relation Functions
//
// Functions:
//		PointOnLine - Checks whether a point is on a line (tolerance may be given)
//
//		PointProject2Line - Find a point on a line projected perpendicular
//							from a given point
//
//		Point2LineDistance - Distance from a point to a line
//
//		PointMirroredAgainstLine - Find a point mirrored with respect to a line
//
// Written in May 12, 2009 by Sungmin Kim
// -----------------------------------------------------------------------------
#include <stdlib.h>
#include <complex>
#include "gf_types.h"
#include "gfunc.h"

bool PointOnLine(double x, double y, double x1, double y1, double x2, double y2)
{
	if (x < x1 && x < x2) return false;
	if (x > x1 && x > x2) return false;
	if (y < y1 && y < y2) return false;
	if (y > y1 && y > y2) return false;
	if (x1 == x2) {
		if (x == x1)
			return true;
		else
			return false;
	}
	if (y1 == y2) {
		if (y == y1)
			return true;
		else
			return false;
	}

    double a, b;
    a = (y2 - y1) / (x2 - x1);
    b = y1 - a * x1;
    return (y == (a * x + b));
}

bool PointOnLine(gf_point p, gf_point v1, gf_point v2)
{
	return PointOnLine(p.x, p.y, v1.x, v1.y, v2.x, v2.y);
}

bool PointOnLine(double x, double y, double x1, double y1, double x2, double y2, double eps)
{
	if (x < x1 - eps && x < x2 - eps) return false;
	if (x > x1 + eps && x > x2 + eps) return false;
	if (y < y1 - eps && y < y2 - eps) return false;
	if (y > y1 + eps && y > y2 + eps) return false;

	if (ABS(x - x1) <= eps && ABS(y - y1) <= eps) return true;
	if (ABS(x - x2) <= eps && ABS(y - y2) <= eps) return true;

	if (ABS(x1 - x2) <= eps) {
		if (ABS(x - x1) <= eps)
			return true;
		else
			return false;
	}
	if (ABS(y1 - y2) <= eps) {
		if (ABS(y - y1) <= eps)
			return true;
		else
			return false;
	}

    double a1, b1;
    a1 = (y2 - y1) / (x2 - x1);
    b1 = y1 - a1 * x1;

	double a2, b2;
	a2 = -1 / a1;
	b2 = y - a2 * x;

	double xp = (b1 - b2) / (a2 - a1);
	double yp = a1 * xp + b1;

	double dx = xp - x;
	double dy = yp - y;
	double d = sqrt(dx * dx + dy * dy);

    return (d <= eps);
}

bool PointOnLine(gf_point p, gf_point v1, gf_point v2, double eps)
{
	return PointOnLine(p.x, p.y, v1.x, v1.y, v2.x, v2.y, eps);
}

bool PointProject2Line(double x, double y, double x1, double y1, double x2, double y2, double *px, double *py)
{
	double a, b, c, d, e;

	a = y1 - y2;
	b = x2 - x1;
	c = x1 * y2 - x2 * y1;

	// Projection to a horizontal line
	if (ABS(a) <= .0001 * ABS(b))
	{
		*px = x;
		*py = y1;
	}
	// Projection to non-horizontal lines
	else
	{
		d = b * x - a * y;
		e = a * a + b * b;
		*px = (-a * c + b * d) / e;
		*py = (-b * c - a * d) / e;
	}

	e = 1.0e-6 * AMAX(ABS(a), ABS(b));

	return PointOnLine(*px, *py, x1, y1, x2, y2, e);
}

bool PointProject2Line(gf_point p, gf_point v1, gf_point v2, gf_point *proj)
{
	return PointProject2Line(p.x, p.y, v1.x, v1.y, v2.x, v2.y, &(proj->x), &(proj->y));
}

double Point2LineDistance(double x, double y, double x1, double y1, double x2, double y2, bool projected_line, double *px, double *py)
{
	// If 'projected_line' variable is true, calculate the distance to the point projected perpendicularly
	// to an infinite line that connects (x1, y1) and (x2, y2)
	
	// If 'projected_line' variable if false, calculate the shortest distance from the point
	// to the line segment that connects (x1, y1) and (x2, y2)
	
	bool bProjected = PointProject2Line(x, y, x1, y1, x2, y2, px, py);
	
	double distance = (x - *px) * (x - *px) + (y - *py) * (y - *py);
	distance = sqrt(distance);
	
	if (!projected_line && !bProjected)
	{
		double a = (x - x1) * (x - x1) + (y - y1) * (y - y1);
		double b = (x - x2) * (x - x2) + (y - y2) * (y - y2);
		if (a < b)
		{
			distance = sqrt(a);
			*px = x1;
			*py = y1;
		}
		else
		{
			distance = sqrt(b);
			*px = x2;
			*py = y2;
		}
	}
	
	return distance;
}

double Point2LineDistance(double x, double y, double x1, double y1, double x2, double y2, bool projected_line)
{
	// If 'projected_line' variable is true, calculate the distance to the point projected perpendicularly
	// to an infinite line that connects (x1, y1) and (x2, y2)

	// If 'projected_line' variable if false, calculate the shortest distance from the point
	// to the line segment that connects (x1, y1) and (x2, y2)

	double px, py;

	double distance = Point2LineDistance(x, y, x1, y1, x2, y2, projected_line, &px, &py);

	return distance;
}

double Point2LineDistance(gf_point p, gf_point v1, gf_point v2, bool projected_line, gf_point *pxy)
{
	return Point2LineDistance(p.x, p.y, v1.x, v1.y, v2.x, v2.y, projected_line, &(pxy->x), &(pxy->y)); 
}

double Point2LineDistance(gf_point p, gf_point v1, gf_point v2, bool projected_line)
{
	return Point2LineDistance(p.x, p.y, v1.x, v1.y, v2.x, v2.y, projected_line); 
}


bool PointMirroredAgainstLine(double x, double y, double x1, double y1, double x2, double y2, double *px, double *py)
{
	bool bOnLine = PointProject2Line(x, y, x1, y1, x2, y2, px, py);
	*px -= (x - (*px));
	*py -= (y - (*py));
	return bOnLine;
}

bool PointMirroredAgainstLine(gf_point p, gf_point v1, gf_point v2, gf_point *mirroredPt)
{
	return PointMirroredAgainstLine(p.x, p.y, v1.x, v1.y, v2.x, v2.y, &(mirroredPt->x), &(mirroredPt->y));
}

int nPointNearStrip(double x, double y, double x1, double y1, double x2, double y2, double dWidth, double *dDistance, double *px, double *py)
{
	// Determine the distance from the reference point (x,y) to the line (x1,y1)-(x2,y2).
	// If the point (x,y) is within the strip with the width 'dWidth' along the (x1,y1)-(x2,y2) line, return true.
	// Otherwise, return false.
	// For either true or false, it returns the distance 'dDistance' and the closest point (px,py).

	double eps = 1.0e-9 * dWidth;

	if (fabs(x - x1) < eps && fabs(y - y1) < eps)
	{
		*px = x;
		*py = y;
		return GF_PNS_AT_BEGIN_POINT;	// The reference point coincides with the beginning point of the line.
	}
	else if (fabs(x - x2) < eps && fabs(y - y2) < eps)
	{
		*px = x;
		*py = y;
		return GF_PNS_AT_END_POINT;	// The reference point coincides with the ending point of the line.
	}

	double dTolerance = 0.5 * dWidth;

	*dDistance = Point2LineDistance(x, y, x1, y1, x2, y2, false, px, py);

	if (*dDistance >= dTolerance)
	{
		*dDistance = -fabs(*dDistance);
		return GF_PNS_FAR_FROM_STRIP;	// The reference point is farther than the 1/2 width from the center line
	}
	else if (*dDistance == 0)
	{
		return GF_PNS_ON_CENTER_LINE;	// The point is on the center line and between the two end points
	}

	int nPoints;
	gf_point *pRectPoints = pLine2Contour(x1, y1, x2, y2, dWidth, 0, 0, &nPoints);

	if (pRectPoints && nPoints > 3 && PointInContour(x, y, nPoints, pRectPoints))
	{
		gfFreePoints(pRectPoints);
		return GF_PNS_WITHIN_STRIP;		// The reference point is within the 1/2 width from the center line
	}
	gfFreePoints(pRectPoints);

	double dist1, dist2;
	double dX = x - x1;
	double dY = y - y1;
	dist1 = sqrt(dX * dX + dY * dY);

	dX = x - x2;
	dY = y - y2;
	dist2 = sqrt(dX * dX + dY * dY);


	if (dist1 < dist2)
	{
		return GF_PNS_NEAR_BEGIN_POINT;	// The reference point is within 1/2 width from the beginning point
	}
	else
	{
		return GF_PNS_NEAR_END_POINT;	// The reference point is within 1/2 width from the beginning point
	}
}
int nPointNearStrip(gf_point p, gf_point v1, gf_point v2, double dWidth, double *dDistance, gf_point *nearestPt)
{
	return nPointNearStrip(p.x, p.y, v1.x, v1.y, v2.x, v2.y, dWidth, dDistance, &(nearestPt->x), &(nearestPt->y));
}