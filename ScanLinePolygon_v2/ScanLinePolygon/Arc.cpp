// -----------------------------------------------------------------------------
// Arc Functions
//
// Functions:
//		pDivideArc - Divide arc to number of line segments
//		getArcCenter - Determines center of an arc
// -----------------------------------------------------------------------------
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "gf_types.h"
#include "gfunc.h"

gf_point *pDivideArc(double dXC, double dYC,			// Arc center
					 double dX1, double dY1,			// Starting arc point
					 double dX2, double dY2,			// Ending arc point
					 bool bCCW,							// Arc direction CCW (true), CW (false)
					 double dDivAngle,					// Approximate division angle
					 int nPointsMin, int nPointsMax,	// Min and max number of points
					 int *nRetDivPoints)				// outputs number of points in the result point list
{
	// Divide arc into line segments that are close to the specified division angle 'dDivAngle'

	gf_point *pResult;
	
	double dRadius;
	dRadius = (dX1 - dXC) * (dX1 - dXC) + (dY1 - dYC) * (dY1 - dYC);
	dRadius = sqrt(dRadius);

	double dAngle1, dAngle2;
	bool bCircle = false;

	if (dX1 == dX2 && dY1 == dY2)
	{
		// If the starting and ending points are the same, treat it as a circle.
		dAngle1 = 0;
		dAngle2 = 360;
		bCircle = true;
	}
	else
	{
		// Determine starting and ending dAngleles
		if (bCCW)
		{
			dAngle1 = dAngleFromEast(dXC, dYC, dX1, dY1);
			dAngle2 = dAngleFromEast(dXC, dYC, dX2, dY2);
		}
		else
		{
			dAngle1 = dAngleFromEast(dXC, dYC, dX2, dY2);
			dAngle2 = dAngleFromEast(dXC, dYC, dX1, dY1);
		}
	}

	double dAngle = dAngle2 - dAngle1;
	if (dAngle < 0) dAngle += 360;

	// Number of division points including both end points
	int nPoints = NINT(dAngle / dDivAngle) + 1;
	if (nPoints < 2) nPoints = 2;
	if (nPoints < nPointsMin) nPoints = nPointsMin;
	if (nPoints > nPointsMax) nPoints = nPointsMax;

	if (bCircle)
	{
		pResult = pCircle2Contour(dXC, dYC, 2 * dRadius, nPoints);
		*nRetDivPoints = nPoints;
		return pResult;
	}

	pResult = (gf_point *) malloc(nPoints * sizeof(gf_point));
	
	// The first and the last points must remain the same
	pResult[0].x = dX1;
	pResult[0].y = dY1;
	pResult[nPoints-1].x = dX2;
	pResult[nPoints-1].y = dY2;

	if (nPoints == 2)
	{
		// If only two points (the arc angle is smaller than the division angle),
		// no arc division is necessary.
		*nRetDivPoints = nPoints;
		return pResult;
	}

	// Angle increment to divide the arc
	dAngle = dAngle / (double) (nPoints - 1);
	dAngle = ONEPI * dAngle / 180;

	// Calculate all internal division points
	int n;
	dAngle1 = ONEPI * dAngle1 / 180;
	if (bCCW)
	{
		for (n = 1; n < nPoints-1; n++)
		{
			dAngle1 += dAngle;
			pResult[n].x = dXC + dRadius * cos(dAngle1);
			pResult[n].y = dYC + dRadius * sin(dAngle1);
		}
	}
	else
	{
		for (n = nPoints-2; n >= 0; n--)
		{
			dAngle1 += dAngle;
			pResult[n].x = dXC + dRadius * cos(dAngle1);
			pResult[n].y = dYC + dRadius * sin(dAngle1);
		}
	}

	*nRetDivPoints = nPoints;
	return pResult;
}

gf_point *pDivideArc(gf_point pCenter,					// Arc center
					 gf_point p1, gf_point p2,			// Starting and ending arc points
					 bool bCCW,							// Arc direction CCW (true), CW (false)
					 double dDivAngle,					// Approximate division angle
					 int nPointsMin, int nPointsMax,	// Min and max number of points
					 int *nRetDivPoints)				// outputs number of points in the result point list
{
	// Divide arc into line segments that are close to the specified division angle 'dDivAngle'

	return pDivideArc(pCenter.x, pCenter.y, p1.x, p1.y, p2.x, p2.y,
					  bCCW, dDivAngle, nPointsMin, nPointsMax, nRetDivPoints);
}

gf_point *pDivideArc(double dXC, double dYC,			// Arc center
					 double dX1, double dY1,			// Starting arc point
					 double dX2, double dY2,			// Ending arc point
					 bool bCCW,							// Arc direction CCW (true), CW (false)
					 int nDivPoints)					// number of division points
{
	// Divide arc into (nDivPoints - 1) line segments.
	// The nDivPoints must be greater than 1.
	
	gf_point *pResults;
	
	double dRadius;
	dRadius = (dX1 - dXC) * (dX1 - dXC) + (dY1 - dYC) * (dY1 - dYC);
	dRadius = sqrt(dRadius);

	double dAngle1, dAngle2;

	if (dX1 == dX2 && dY1 == dY2)
	{
		// If the starting and ending points are the same, treat it as a circle.
		pResults = pCircle2Contour(dXC, dYC, 2 * dRadius, nDivPoints);
		return pResults;
	}
	else
	{
		// Determine starting and ending dAngleles
		if (bCCW)
		{
			dAngle1 = dAngleFromEast(dXC, dYC, dX1, dY1);
			dAngle2 = dAngleFromEast(dXC, dYC, dX2, dY2);
		}
		else
		{
			dAngle1 = dAngleFromEast(dXC, dYC, dX2, dY2);
			dAngle2 = dAngleFromEast(dXC, dYC, dX1, dY1);
		}
	}
	double dAngle = dAngle2 - dAngle1;
	if (dAngle < 0) dAngle += 360;

	pResults = (gf_point *) malloc(nDivPoints * sizeof(gf_point));

	// The first and the last points must remain the same
	pResults[0].x = dX1;
	pResults[0].y = dY1;
	pResults[nDivPoints-1].x = dX2;
	pResults[nDivPoints-1].y = dY2;

	// If only two points, no arc division is necessary.
	if (nDivPoints == 2) return pResults;

	// Angle increment to divide the arc
	
	dAngle = dAngle / (double) (nDivPoints - 1);
	dAngle = ONEPI * dAngle / 180;

	// Calculate all internal division points
	int n;
	dAngle1 = ONEPI * dAngle1 / 180;
	if (bCCW)
	{
		for (n = 1; n < nDivPoints-1; n++)
		{
			dAngle1 += dAngle;
			pResults[n].x = dXC + dRadius * cos(dAngle1);
			pResults[n].y = dYC + dRadius * sin(dAngle1);
		}
	}
	else
	{
		for (n = nDivPoints-2; n >= 0; n--)
		{
			dAngle1 += dAngle;
			pResults[n].x = dXC + dRadius * cos(dAngle1);
			pResults[n].y = dYC + dRadius * sin(dAngle1);
		}
	}

	return pResults;
}

gf_point *pDivideArc(gf_point pCenter,					// Arc center
					 gf_point p1, gf_point p2,			// Starting and ending arc points
					 bool bCCW,							// Arc direction CCW (true), CW (false)
					 int nDivPoints)					// number of division points
{
	// Divide arc into (nDivPoints - 1) line segments.
	// The nDivPoints must be greater than 1.

	return pDivideArc(pCenter.x, pCenter.y, p1.x, p1.y, p2.x, p2.y, bCCW, nDivPoints);
}

gf_point getArcCenter(double dXstart, double dYstart,	// start point of the arc
					  double dXend, double dYend,		// end point of the arc
					  double dXmid, double dYmid)		// a point on the arc between start and end points
{
	// Returns a center coordinate of arc that is defined by three points.

	gf_point dCenter;

	double dSMALL = 1.0e-10;

	// Y = a1 * X + b1 for the line perpendicular to the line connecting start and mid points
	// Y = a2 * X + b2 for the line perpendicular to the line connecting mid and end points

	double a1, a2, b1, b2;

	if (ABS(dYmid - dYstart) > dSMALL)
	{
		a1 = - (dXmid - dXstart) / (dYmid - dYstart);
		b1 = 0.5 * ((dYstart + dYmid) - a1 * (dXstart + dXmid));
	}
	else
	{
		a1 = DBL_MAX;
		b1 = 0;
	}

	if (ABS(dYend - dYmid) > dSMALL)
	{
		a2 = - (dXend - dXmid) / (dYend - dYmid);
		b2 = 0.5 * ((dYmid+dYend) - a2 * (dXmid+dXend));
	}
	else
	{
		a2 = DBL_MAX;
		b2 = 0;
	}

	if ((ABS(dXend - dXstart) <= dSMALL && ABS(dYend - dYstart) < dSMALL) ||
		(ABS(dYmid - dYstart) <= dSMALL && ABS(dYend - dYmid) < dSMALL) ||
		ABS(a1 - a2) <= dSMALL)
	{
		// If the three points are collinear, output raidus = 0
		dCenter.x = dXstart;
		dCenter.y = dYstart;
	}
	else if (a1 > 1000)
	{
		// If the line connecting the start and mid points is near horizontal
		dCenter.x = 0.5 * (dXstart + dXmid);
		dCenter.y = a2 * dCenter.x + b2;
	}
	else if (a2 > 1000)
	{
		// If the line connecting the mid and end points is near horizontal
		dCenter.x = 0.5 * (dXmid + dXend);
		dCenter.y = a1 * dCenter.x + b1;
	}
	else
	{
		// For all other cases
		dCenter.x = (b1 - b2) / (a2 - a1);
		dCenter.y = a1 * dCenter.x + b1;
	}

	return dCenter;
}

gf_point getArcCenter(double dXstart, double dYstart,	// start point of the arc
					  double dXend, double dYend,		// end point of the arc
					  double dXoffset, double dYoffset,	// X and Y direction offsets of center from the start point
					  bool bEnable360)					// 360 degree enable flag
{
	// Returns a center coordinate of arc that is defined by two arc end points and offset distance to arc center
	// For more details, refer to the Circular Interpolation (G02, G03, G74, G75) in Gerber RS274X specification.

	// For 360 degree enabled arc, the center is located "dXoffset" amount in X-direction and
	// "dYoffset" amount in Y-direction from the starting point. "dXoffset" and "dYoffset" can
	// be either positive or negative.

	// For 360 degree disabled arc, the center is located "dXoffset" amount in X-direction and
	// "dYoffset" amount in Y-direction from the starting point. However "dXoffset" and "dYoffset"
	// are absolute values in this case, so the center is determined by finding the location that
	// is the same distance from both starting and ending arc points.

	gf_point dCenter;

	if (bEnable360)
	{
		// Multiquadrant (360? Circular Interpolation (G74, G75

		dCenter.x = dXstart + dXoffset;
		dCenter.y = dYstart + dYoffset;
	}
	else
	{
		// Single Quadrant Circular Interpolation (G74)

		// True radius^2 of the arc
		double dRadius = dXoffset * dXoffset + dYoffset * dYoffset;

		// Let's try all four quadrants and find the proper quadrant
		double dXabs = ABS(dXoffset);
		double dYabs = ABS(dYoffset);

		double dDist, dMinDiff;
		gf_point dTemp;

		// 1st quadrant
		dCenter.x = dXstart + dXabs;
		dCenter.y = dYstart + dYabs;
		dDist = (dCenter.x - dXend) * (dCenter.x - dXend) + (dCenter.y - dYend) * (dCenter.y - dYend);
		dMinDiff = ABS(dDist - dRadius);
		
		// 2nd quadrant
		dTemp.x = dXstart - dXabs;
		dTemp.y = dYstart + dYabs;
		dDist = (dTemp.x - dXend) * (dTemp.x - dXend) + (dTemp.y - dYend) * (dTemp.y - dYend);
		if (ABS(dDist - dRadius) < dMinDiff)
		{
			dMinDiff = ABS(dDist - dRadius);
			dCenter.x = dTemp.x;
			dCenter.y = dTemp.y;
		}
		
		// 3th quadrant
		dTemp.x = dXstart - dXabs;
		dTemp.y = dYstart - dYabs;
		dDist = (dTemp.x - dXend) * (dTemp.x - dXend) + (dTemp.y - dYend) * (dTemp.y - dYend);
		if (ABS(dDist - dRadius) < dMinDiff)
		{
			dMinDiff = ABS(dDist - dRadius);
			dCenter.x = dTemp.x;
			dCenter.y = dTemp.y;
		}
		
		// 4th quadrant
		dTemp.x = dXstart + dXabs;
		dTemp.y = dYstart - dYabs;
		dDist = (dTemp.x - dXend) * (dTemp.x - dXend) + (dTemp.y - dYend) * (dTemp.y - dYend);
		if (ABS(dDist - dRadius) < dMinDiff)
		{
			dCenter.x = dTemp.x;
			dCenter.y = dTemp.y;
		}
	}

	return dCenter;
}


gf_point getArcCenter(double x1s, double y1s, double x1e, double y1e,	// start and end points of the first tangent line. Tangent at the start point.
					  double x2s, double y2s, double x2e, double y2e)	// start and end points of the second tangent line. Tangent at the start point.
{
	// Returns a center coordinate of arc that is defined by two tangent lines.
	
	gf_point dCenter;
	
	double dSMALL = 1.0e-10;
	
	// Y = a1 * X + b1 for the line perpendicular to the first tangent line
	// Y = a2 * X + b2 for the line perpendicular to the second tangent line
	
	double a1, a2, b1, b2;
	
	if (ABS(y1e - y1s) > dSMALL)
	{
		a1 = - (x1e - x1s) / (y1e - y1s);
		b1 = y1s - a1 * x1s;
	}
	else
	{
		a1 = DBL_MAX;
		b1 = 0;
	}
	
	if (ABS(y2e - y2s) > dSMALL)
	{
		a2 = - (x2e - x2s) / (y2e - y2s);
		b2 = y2s - a2 * x2s;
	}
	else
	{
		a2 = DBL_MAX;
		b2 = 0;
	}
	
	if ((ABS(x1e - x1s) <= dSMALL && ABS(y1e - y1s) < dSMALL) ||
		(ABS(x2e - x2s) <= dSMALL && ABS(y2e - y2s) < dSMALL) ||
		ABS(a1 - a2) <= dSMALL)
	{
		// If the three points are collinear, output raidus = 0
		dCenter.x = x1s;
		dCenter.y = y1s;
	}
	else if (a1 > 1000)
	{
		// If the line connecting the start and mid points is near horizontal
		dCenter.x = x1s;
		dCenter.y = a2 * dCenter.x + b2;
	}
	else if (a2 > 1000)
	{
		// If the line connecting the mid and end points is near horizontal
		dCenter.x = x2s;
		dCenter.y = a1 * dCenter.x + b1;
	}
	else
	{
		// For all other cases
		dCenter.x = (b1 - b2) / (a2 - a1);
		dCenter.y = a1 * dCenter.x + b1;
	}
	
	return dCenter;
}

gf_box getArcBox(double dXC, double dYC,			// Arc center
				 double dX1, double dY1,			// Starting arc point
				 double dX2, double dY2,			// Ending arc point
				 bool bCCW)						// Arc direction CCW (true), CW (false)
{
	// Returns a bounding box of the arc
	gf_box box;

	double dRadius;;
	dRadius = (dX1 - dXC) * (dX1 - dXC) + (dY1 - dYC) * (dY1 - dYC);
	dRadius = sqrt(dRadius);

	// Box around a circle
	box.ll.x = dXC - dRadius;
	box.ll.y = dYC - dRadius;
	box.ur.x = dXC + dRadius;
	box.ur.y = dYC + dRadius;

	if (dX1 != dX2 || dY1 != dY2)
	{
		// Determine starting and ending dAngleles
		double dAngle1, dAngle2;

		if (bCCW)
		{
			dAngle1 = dAngleFromEast(dXC, dYC, dX1, dY1);
			dAngle2 = dAngleFromEast(dXC, dYC, dX2, dY2);
		}
		else
		{
			dAngle1 = dAngleFromEast(dXC, dYC, dX2, dY2);
			dAngle2 = dAngleFromEast(dXC, dYC, dX1, dY1);
		}
		if (dAngle1 <= 90)
		{
			if (dAngle2 <= 90)
			{
				if (dAngle2 > dAngle1)
				{
					box.ll.x = AMIN(dX1, dX2);
					box.ur.x = AMAX(dX1, dX2);
					box.ll.y = AMIN(dY1, dY2);
					box.ur.y = AMAX(dY1, dY2);
				}
			}
			else if (dAngle2 <= 180)
			{
				box.ll.x = AMIN(dX1, dX2);
				box.ur.x = AMAX(dX1, dX2);
				box.ll.y = AMIN(dY1, dY2);
			}
			else if (dAngle2 <= 270)
			{
				box.ur.x = AMAX(dX1, dX2);
				box.ll.y = AMIN(dY1, dY2);
			}
			else
			{
				box.ur.x = AMAX(dX1, dX2);
			}
		}
		else if (dAngle1 <= 180)
		{
			if (dAngle2 <= 90)
			{
				box.ur.y = AMAX(dY1, dY2);
			}
			else if (dAngle2 <= 180)
			{
				if (dAngle2 > dAngle1)
				{
					box.ll.x = AMIN(dX1, dX2);
					box.ur.x = AMAX(dX1, dX2);
					box.ll.y = AMIN(dY1, dY2);
					box.ur.y = AMAX(dY1, dY2);
				}
			}
			else if (dAngle2 <= 270)
			{
				box.ur.x = AMAX(dX1, dX2);
				box.ll.y = AMIN(dY1, dY2);
				box.ur.y = AMAX(dY1, dY2);
			}
			else
			{
				box.ur.x = AMAX(dX1, dX2);
				box.ur.y = AMAX(dY1, dY2);
			}
		}
		else if (dAngle1 <= 270)
		{
			if (dAngle2 <= 90)
			{
				box.ll.x = AMIN(dX1, dX2);
				box.ur.y = AMAX(dY1, dY2);
			}
			else if (dAngle2 <= 180)
			{
				box.ll.x = AMIN(dX1, dX2);
			}
			else if (dAngle2 <= 270)
			{
				if (dAngle2 > dAngle1)
				{
					box.ll.x = AMIN(dX1, dX2);
					box.ur.x = AMAX(dX1, dX2);
					box.ll.y = AMIN(dY1, dY2);
					box.ur.y = AMAX(dY1, dY2);
				}
			}
			else
			{
				box.ll.x = AMIN(dX1, dX2);
				box.ur.x = AMAX(dX1, dX2);
				box.ur.y = AMAX(dY1, dY2);
			}
		}
		else
		{
			if (dAngle2 <= 90)
			{
				box.ll.x = AMIN(dX1, dX2);
				box.ll.y = AMIN(dY1, dY2);
				box.ur.y = AMAX(dY1, dY2);
			}
			else if (dAngle2 <= 180)
			{
				box.ll.x = AMIN(dX1, dX2);
				box.ll.y = AMIN(dY1, dY2);
			}
			else if (dAngle2 <= 270)
			{
				box.ll.y = AMIN(dY1, dY2);
			}
			else
			{
				if (dAngle2 > dAngle1)
				{
					box.ll.x = AMIN(dX1, dX2);
					box.ur.x = AMAX(dX1, dX2);
					box.ll.y = AMIN(dY1, dY2);
					box.ur.y = AMAX(dY1, dY2);
				}
			}
		}
	}

	return box;
}
