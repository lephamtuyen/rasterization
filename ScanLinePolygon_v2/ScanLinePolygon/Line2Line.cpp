// -----------------------------------------------------------------------------
// Line-To-Line Relation Functions
//
// Functions:
//		dCrossProduct - Determine a cross product of two lines
//		dLineOverlap - Determine an overlapped length of two parallel lines
//		bTwoLinesParallel - Determine if two lines are parallel
//		bMergeTwoLines - Merge two lines if they are connected and parallel
//		bTwoLinesCross - Find a cross crossing point if extensions of two lines.
//						 Returns true or false.
//		nTwoLines -	Find a crossing point of extensions of two lines.
//					Returns flags -1, 0, 1, 2 and 3
//		dLine2LineDistance - Returns the closest distance between two lines.
//							 It returns a negative value if the two lines are crossing.
//						   - It also returns the two closest points.
//		vJoinLines - Join two lines with an arc or at the crossing point
//		bGetParallelLine - Determine a parallel line apart the given distance from
//						   the given line
//		pGetPathApart - Determine a path (polyline) from the reference path apart the
//                      specified distance on the right or left.

//		dAngleFromEast - Determine angle of a line from east
//		dAngleTwoLines - Determine angle between two separate lines
//		dAngle3Points - Determine angle between two lines formed by 3 consecutive points
// -----------------------------------------------------------------------------
#include <stdlib.h>
#include <float.h>
#include <complex>
#include "gf_types.h"
#include "gfunc.h"

double dCrossProduct(double x1s, double y1s, double x1e, double y1e,
					 double x2s, double y2s, double x2e, double y2e)
{
	double dx1, dy1, dx2, dy2, d;

	dx1 = x1e - x1s;
	dy1 = y1e - y1s;
	dx2 = x2e - x2s;
	dy2 = y2e - y2s;

	// Cross product
	d = dx1 * dy2 - dx2 * dy1;

	return d;
}

double dCrossProduct(gf_point *L1, gf_point *L2)
{
	return dCrossProduct(L1[0].x, L1[0].y, L1[1].x, L1[1].y,
						 L2[0].x, L2[0].y, L2[1].x, L2[1].y);
}

double dLineOverlap(double x1s, double y1s, double x1e, double y1e,		// Reference line
					double x2s, double y2s, double x2e, double y2e,		// Check line
					double dParallelCheckAngle,							// A small angle (degree) to check for parallel
					double *x1, double *y1, double *x2, double *y2)		// Overlapped section on the reference line
{
	// Returns zero if not parallel.
	// Returns negative if the two lines do not overlap.

	double dAng = dAngleTwoLines(x1s, y1s, x1e, y1e, x2s, y2s, x2e, y2e);
	dAng = ABS(dAng);
	if (dAng < dParallelCheckAngle || ABS(180-dAng) < dParallelCheckAngle)
	{
		double xs, ys, xe, ye;
		double s1, e1, s2, e2;
		double dX, dY;
		double sp, ep;
		double dOverlapLength;

		PointProject2Line(x2s, y2s, x1s, y1s, x1e, y1e, &xs, &ys);
		PointProject2Line(x2e, y2e, x1s, y1s, x1e, y1e, &xe, &ye);

		dX = x1e - x1s;
		dY = y1e - y1s;
		double d = sqrt(dX * dX + dY * dY);
		double c = dX / d;	// cosine
		double s = dY / d;	// sine
		
		s1 = 0;
		e1 = d;

		dX = xs - x1s;
		dY = ys - y1s;
		s2 = c * dX + s * dY;

		dX = xe - x1s;
		dY = ye - y1s;
		e2 = c * dX + s * dY;

		if (e2 < s2)
		{
			dX = s2;
			s2 = e2;
			e2 = dX;
		}

		// Distance between two lines if not overlapped
		if (s2 >= e1)
		{
			sp = e1;
			ep = s2;
			dOverlapLength = sp - ep;
		}
		else if (s1 >= e2)
		{
			sp = e2;
			ep = s1;
			dOverlapLength = sp - ep;
		}

		// Overlapped length
		
		else
		{
			if (s2 >= s1)
				sp = s2;
			else
				sp = s1;

			if (e2 >= e1)
				ep = e1;
			else
				ep = e2;

			dOverlapLength = ep - sp;
		}

		*x1 = c * sp + x1s;
		*y1 = s * sp + y1s;
		*x2 = c * ep + x1s;
		*y2 = s * ep + y1s;

		return dOverlapLength;
	}
	else
	{
		return 0;
	}
}

double dLineOverlap(gf_point p1s, gf_point p1e,			// Reference line
					gf_point p2s, gf_point p2e,			// Check line
					double dParallelCheckAngle,			// A small angle (degree) to check for parallel
					gf_point *pFrom, gf_point *pTo)		// Overlapped section on the reference line
{
	return dLineOverlap(p1s.x, p1s.y, p1e.x, p1e.y,
						p2s.x, p2s.y, p2e.x, p2e.y,
						dParallelCheckAngle,
						&(pFrom->x), &(pFrom->y), &(pTo->x), &(pTo->y));

}

bool bTwoLinesParallel(	double x1s, double y1s, double x1e, double y1e,
						double x2s, double y2s, double x2e, double y2e)
{
	double dx1, dy1, dx2, dy2, d;

	dx1 = x1e - x1s;
	dy1 = y1e - y1s;
	dx2 = x2e - x2s;
	dy2 = y2e - y2s;
	
	double dMax = fabs(dx1);
	if (fabs(dy1) > dMax) dMax = fabs(dy1);
	if (fabs(dx2) > dMax) dMax = fabs(dx2);
	if (fabs(dy2) > dMax) dMax = fabs(dy2);

	// Cross product
	d = dx1 * dy2 - dx2 * dy1;
	if (fabs(d) <= 1.0e-9 * dMax)
		return true;	// Parallel
	else
		return false;	// Not parallel

}

bool bTwoLinesParallel(gf_point *L1, gf_point *L2)
{
	return bTwoLinesParallel(L1[0].x, L1[0].y, L1[1].x, L1[1].y,
							 L2[0].x, L2[0].y, L2[1].x, L2[1].y);
}

int MergeTwoLines(double x1s, double y1s, double x1e, double y1e,
				  double x2s, double y2s, double x2e, double y2e, double tol /* = 0.0 */)
{
	// Merge two lines if they are connected and collinear

	int retFlag = 0;

	if (!bTwoLinesParallel(x1s, y1s, x1e, y1e, x2s, y2s, x2e, y2e)) return retFlag;
	
	double dx1, dy1, dx2, dy2, dMax;

	double tolerance = tol;
	if (tolerance == 0)
	{
		dx1 = x1e - x1s;
		dy1 = y1e - y1s;
		dx2 = x2e - x2s;
		dy2 = y2e - y2s;
		
		dMax = fabs(x1s);
		if (fabs(dy1) > dMax) dMax = fabs(dy1);
		if (fabs(dx2) > dMax) dMax = fabs(dx2);
		if (fabs(dy2) > dMax) dMax = fabs(dy2);
		tolerance = 1.0e-9 * dMax;
	}

	if ((fabs(x1s - x2s) <= tolerance && fabs(y1s - y2s) <= tolerance) ||
		(fabs(x1s - x2e) <= tolerance && fabs(y1s - y2e) <= tolerance) ||
		(fabs(x1e - x2s) <= tolerance && fabs(y1e - y2s) <= tolerance) ||
		(fabs(x1e - x2e) <= tolerance && fabs(y1e - y2e) <= tolerance))
	{
		dMax = (x1s - x2s) * (x1s - x2s) + (y1s - y2s) * (y1s - y2s);
		retFlag = 11;

		dx1 = (x1s - x2e) * (x1s - x2e) + (y1s - y2e) * (y1s - y2e);
		if (dx1 > dMax)
		{
			dMax = dx1;
			retFlag = 12;
		}
		dx1 = (x1e - x2s) * (x1e - x2s) + (y1e - y2s) * (y1e - y2s);
		if (dx1 > dMax)
		{
			dMax = dx1;
			retFlag = 21;
		}
		dx1 = (x1e - x2e) * (x1e - x2e) + (y1e - y2e) * (y1e - y2e);
		if (dx1 > dMax)
		{
			retFlag = 22;
		}
	}

	return retFlag;
}

int MergeTwoLines(gf_point *L1, gf_point *L2, double tol /* = 0.0 */)
{
	return MergeTwoLines(L1[0].x, L1[0].y, L1[1].x, L1[1].y,
						 L2[0].x, L2[0].y, L2[1].x, L2[1].y, tol);
}

bool bTwoLinesCross(double x1s, double y1s, double x1e, double y1e,
					double x2s, double y2s, double x2e, double y2e,
					double *x, double *y)
{
#if 0	// 2012/06/12 mjkim - 두 선분의 교점을 확인 하기 위하여 함수를 수정함. (관련 Task 확인 불가.)
		// 2012/08/30 mjkim - 이 함수는 선분에 대한 교차 검사가 아니고 선분이 형성하는 직선의 교점을 확인하기 위한 함수이기 때문에
		// #4477의 문제를 해결하기 위해 이전 Code로 돌림.
	enum
	{
		X = 0,
		Y = 1,
	};

	// make vector
	double v1s[] = {x1s, y1s}, v1e[] = {x1e, y1e}, 
		   v2s[] = {x2s, y2s}, v2e[] = {x2e, y2e};

	// make segment
	double s1c[] = {(v1s[X] + v1e[X]) * .5, (v1s[Y] + v1e[Y]) * .5};
	double s1d[] = {v1e[X] - v1s[X], v1e[Y] - v1s[Y]};
	double s1e = sqrt(pow(s1d[X], 2) + pow(s1d[Y], 2)) * .5;

	s1d[X] /= s1e * 2;
	s1d[Y] /= s1e * 2;

	double s2c[] = {(v2s[X] + v2e[X]) * .5, (v2s[Y] + v2e[Y]) * .5};
	double s2d[] = {v2e[X] - v2s[X], v2e[Y] - v2s[Y]};
	double s2e = sqrt(pow(s2d[X], 2) + pow(s2d[Y], 2)) * .5;

	s2d[X] /= s2e * 2;
	s2d[Y] /= s2e * 2;

	// find intersection
	double diff[] = {s2c[X] - s1c[X], s2c[Y] - s1c[Y]};
	double d1 = s1d[X] * s2d[Y] - s1d[Y] * s2d[X];			// dot prep

	if(abs(d1) > 0)	// intersects at point
	{
		*x = (diff[X] * s2d[Y] - diff[Y] * s2d[X]) / d1;	// dot prep / d
		*y = (diff[X] * s1d[Y] - diff[Y] * s1d[X]) / d1;	// dot prep / d

		// find point on segment
		if(abs(*x) <= s1e && abs(*y) <= s2e)
		{
			*x = s1c[X] + *x * s1d[X];
			*y = s1c[Y] + *x * s1d[Y];

			return true;
		}
		else
			return false;
	}

	double difflen = sqrt(pow(diff[X], 2) + pow(diff[Y], 2));

	diff[X] *= difflen;
	diff[Y] *= difflen;

	double d2 = diff[X] * s2d[Y] - diff[Y] * s2d[X];		// dot prep

	if(abs(d2) <= 0)	// same line
		return true;

	return false;	// parallel
#else
	double dx1, dy1, dx2, dy2;
	double dc1, dc2, d;

	dx1 = x1e - x1s;
	dy1 = y1e - y1s;
	dx2 = x2e - x2s;
	dy2 = y2e - y2s;

	double dMax = fabs(dx1);
	if (fabs(dy1) > dMax) dMax = fabs(dy1);
	if (fabs(dx2) > dMax) dMax = fabs(dx2);
	if (fabs(dy2) > dMax) dMax = fabs(dy2);

	// Cross product
	d = dx1 * dy2 - dx2 * dy1;
	if (fabs(d) <= 1.0e-9 * dMax) return false;  // Parallel

	// Find the intersecting point
	dc1 = dy1 * x1s - dx1 * y1s;
	dc2 = dy2 * x2s - dx2 * y2s;
	*x = (-dx2 * dc1 + dx1 * dc2) / d;
	*y = (-dy2 * dc1 + dy1 * dc2) / d;
	return true;
#endif
}

bool bTwoLinesCross(gf_point *L1, gf_point *L2, double *x, double *y)
{
	return bTwoLinesCross(L1[0].x, L1[0].y, L1[1].x, L1[1].y,
						  L2[0].x, L2[0].y, L2[1].x, L2[1].y,
						  x, y);
}

bool bTwoLinesCross(gf_line L1, gf_line L2, double *x, double *y)
{
	return bTwoLinesCross(L1.s.x, L1.s.y, L1.e.x, L1.e.y,
						  L2.s.x, L2.s.y, L2.e.x, L2.e.y,
						  x, y);
}

int nTwoLines(	double x1s, double y1s, double x1e, double y1e,
				double x2s, double y2s, double x2e, double y2e,
				double *x, double *y, double epsilonFactor)
{
	int flag = -1;

	if (bTwoLinesCross(x1s, y1s, x1e, y1e, x2s, y2s, x2e, y2e, x, y))
	{
		flag++;

		double dx1, dy1, dx2, dy2;
		dx1 = fabs(x1e - x1s);
		dy1 = fabs(y1e - y1s);
		dx2 = fabs(x2e - x2s);
		dy2 = fabs(y2e - y2s);

		double dist;
		dist = dx1;
		if (dy1 > dist) dist = dy1;
		if (dx2 > dist) dist = dx2;
		if (dy2 > dist) dist = dy2;
		dist *= epsilonFactor;

		if (PointInBetween(*x, *y, x1s, y1s, x1e, y1e, epsilonFactor) ||
			(fabs(*x - x1s) < dist && fabs(*y - y1s) < dist) ||
			(fabs(*x - x1e) < dist && fabs(*y - y1e) < dist))
		{
			flag += 1;
		}
		if (PointInBetween(*x, *y, x2s, y2s, x2e, y2e, epsilonFactor) ||
			(fabs(*x - x2s) < dist && fabs(*y - y2s) < dist) ||
			(fabs(*x - x2e) < dist && fabs(*y - y2e) < dist))
		{
			flag += 2;
		}
	}
	return flag;

}

int nTwoLines(gf_point *L1, gf_point *L2, double *x, double *y, double epsilonFactor)
{
	return nTwoLines(L1[0].x, L1[0].y, L1[1].x, L1[1].y,
					 L2[0].x, L2[0].y, L2[1].x, L2[1].y,
					 x, y, epsilonFactor);
}

int nTwoLines(gf_line L1, gf_line L2, double *x, double *y, double epsilonFactor)
{
	return nTwoLines(L1.s.x, L1.s.y, L1.e.x, L1.e.y,
					 L2.s.x, L2.s.y, L2.e.x, L2.e.y,
					 x, y, epsilonFactor);
}

int nTwoLines(double x1s, double y1s, double x1e, double y1e,
			  double x2s, double y2s, double x2e, double y2e,
			  double *x, double *y)
{
	int flag = -1;

	if (bTwoLinesCross(x1s, y1s, x1e, y1e, x2s, y2s, x2e, y2e, x, y))
	{
		flag++;
		if (PointInBetween(*x, *y, x1s, y1s, x1e, y1e)) flag += 1;
		if (PointInBetween(*x, *y, x2s, y2s, x2e, y2e)) flag += 2;
	}
	return flag;

}

int nTwoLines(gf_point *L1, gf_point *L2, double *x, double *y)
{
	return nTwoLines(L1[0].x, L1[0].y, L1[1].x, L1[1].y,
					 L2[0].x, L2[0].y, L2[1].x, L2[1].y,
					 x, y);
}

int nTwoLines(gf_line L1, gf_line L2, double *x, double *y)
{
	return nTwoLines(L1.s.x, L1.s.y, L1.e.x, L1.e.y,
					L2.s.x, L2.s.y, L2.e.x, L2.e.y,
					x, y);
}

double dLine2LineDistance(double dXA1, double dYA1, double dXA2, double dYA2,	// Start & end points of line A
						  double dXB1, double dYB1, double dXB2, double dYB2,	// Start & end points of line B
						  bool projected_line)									// if true, consider only the projected distance
{
	double closestDistance = DBL_MAX;
	double dist;
	
	dist = Point2LineDistance(dXA1, dYA1, dXB1, dYB1, dXB2, dYB2, projected_line);
	if (dist < closestDistance) closestDistance = dist;
	
	dist = Point2LineDistance(dXA2, dYA2, dXB1, dYB1, dXB2, dYB2, projected_line);
	if (dist < closestDistance) closestDistance = dist;
	
	dist = Point2LineDistance(dXB1, dYB1, dXA1, dYA1, dXA2, dYA2, projected_line);
	if (dist < closestDistance) closestDistance = dist;
	
	dist = Point2LineDistance(dXB2, dYB2, dXA1, dYA1, dXA2, dYA2, projected_line);
	if (dist < closestDistance) closestDistance = dist;
	
	double dX, dY;
	if (nTwoLines(dXA1, dYA1, dXA2, dYA2, dXB1, dYB1, dXB2, dYB2, &dX, &dY) == 3)
	{
		closestDistance *= -1;		// returns negative distance if the two lines are crossing
	}
	
	return closestDistance;
}

double dLine2LineDistance(double dXA1, double dYA1, double dXA2, double dYA2,	// Start & end points of line A
						  double dXB1, double dYB1, double dXB2, double dYB2,	// Start & end points of line B
						  bool projected_line,									// if true, consider only the projected distance
						  double *dXA, double *dYA,								// Closest point on line A
						  double *dXB, double *dYB)								// Closest point on line B
{
	// If the closest point can be from either start or end point of line A,
	// it will return the one closer to the start point of line A.

	double closestDistance = DBL_MAX;
	double dist;

	double distFromStartOfLineA = DBL_MAX;
	bool bFromStartOfLineA = false;

	double dX, dY;
	double deltaX, deltaY;

	dist = Point2LineDistance(dXA1, dYA1, dXB1, dYB1, dXB2, dYB2, projected_line, &dX, &dY);
	if (dist < closestDistance)
	{
		distFromStartOfLineA = 0;
		bFromStartOfLineA = true;
		*dXA = dXA1;
		*dYA = dYA1;
		*dXB = dX;
		*dYB = dY;
		closestDistance = dist;
	}
	
	dist = Point2LineDistance(dXA2, dYA2, dXB1, dYB1, dXB2, dYB2, projected_line, &dX, &dY);
	if ((dist + 1.0e-4 * dist) < closestDistance)
	{
		// Added a very small value to the distance to skip this distance
		// if the distance from the end point is about the same as the distance from the start point.
		deltaX = dXA2 - dXA1;
		deltaY = dYA2 - dYA1;
		distFromStartOfLineA = sqrt(deltaX * deltaX + deltaY * deltaY);
		bFromStartOfLineA = false;
		*dXA = dXA2;
		*dYA = dYA2;
		*dXB = dX;
		*dYB = dY;
		closestDistance = dist;
	}
	
	dist = Point2LineDistance(dXB1, dYB1, dXA1, dYA1, dXA2, dYA2, projected_line, &dX, &dY);
	if ((dist - 1.0e-4 * dist) < closestDistance)
	{
		// Subtracted a very small value from the distance to take this distance
		// if the distance from the this point is about the same as the previous closest point.
		deltaX = dX - dXA1;
		deltaY = dY - dYA1;
		distFromStartOfLineA = sqrt(deltaX * deltaX + deltaY * deltaY);
		*dXB = dXB1;
		*dYB = dYB1;
		*dXA = dX;
		*dYA = dY;
		closestDistance = dist;
	}
	
	dist = Point2LineDistance(dXB2, dYB2, dXA1, dYA1, dXA2, dYA2, projected_line, &dX, &dY);
	if ((dist - 1.1e-4 * dist) < closestDistance)
	{
		// Subtracted a very small value from the distance to take this distance
		// if the distance from the this point is about the same as the previous closest point.
		// Subtracted more than previous checking. (1.1e-4 vs 1.0e-4).
		deltaX = dX - dXA1;
		deltaY = dY - dYA1;
		double d = sqrt(deltaX * deltaX + deltaY * deltaY);
		if (closestDistance - dist > 1.0e-4 * dist ||		// This distance is closer enough than the previous closest point
			d < distFromStartOfLineA)						// Distance from this point to the start of line A is closer than the previous one
		{
			*dXB = dXB2;
			*dYB = dYB2;
			*dXA = dX;
			*dYA = dY;
			closestDistance = dist;
		}
	}
	
	if (nTwoLines(dXA1, dYA1, dXA2, dYA2, dXB1, dYB1, dXB2, dYB2, &dX, &dY) == 3)
	{
		*dXA = dX;
		*dYA = dY;
		*dXB = dX;
		*dYB = dY;
		closestDistance *= -1;		// returns negative distance if the two lines are crossing
	}
	
	return closestDistance;
}

bool bJoinLines(gf_point arcOrigin,					// Arc origin (unused if bArcJoin = false)
				gf_point *pLine1, gf_point *pLine2,	// two line vectors
				bool bArcJoin,						// join method (true: arc, false: straight line)
				gf_point *pPoint, int *nPoints)		// point array and index of the last new point
{
	// Join two lines with an arc or at the crossing point then add the join points to the existing point array 'pPoint'.
	// The arc origin is at 'arcOrigin'
	//
	// Returns:	true if the two lines are crossed, false if not crossed

	bool bCrossed = false;
	bool bDirectConnect = false;

	int p = *nPoints;

	double x, y;
	int nFlag = nTwoLines(pLine1, pLine2, &x, &y);

	if (nFlag < 0)			// Parallel lines
	{
		bDirectConnect = true;
	}
	else if (nFlag == 3)	// Crossing lines
	{
		pPoint[p].x = x;
		pPoint[p].y = y;
		p++;
		bCrossed = true;
	}
	else if (nFlag == 1 || nFlag == 2)	// Do nothing
	{
		return false;
	}
	else					// Neither parallel or crossing
	{
		if (bArcJoin)
		{
			double dCross = dCrossProduct(pLine1, pLine2);
			if (dCross == 0)
			{
				bDirectConnect = true;
			}
			else
			{
				// Join with an arc
				bool bCCW;
				if (dCross > 0)
					bCCW = true;
				else
					bCCW = false;
				int nArcPoints;
				gf_point *pArc = pDivideArc(arcOrigin, pLine1[1], pLine2[0], bCCW,
											LINE_JOINT_ARC_DIV_ANGLE, 2, LINE_JOINT_ARC_DIV_MAX_POINT, &nArcPoints);
				int i;
				for (i = 0; i < nArcPoints; i++)
				{
					pPoint[p].x = pArc[i].x;
					pPoint[p].y = pArc[i].y;
					p++;
				}
				FREE(pArc);
			}
		}
		else
		{
			bDirectConnect = true;
		}
	}

	if (bDirectConnect)
	{
		pPoint[p].x = pLine1[1].x;
		pPoint[p].y = pLine1[1].y;
		p++;
		pPoint[p].x = pLine2[0].x;
		pPoint[p].y = pLine2[0].y;
		p++;
	}

	*nPoints = p;

	return bCrossed;
}

bool bJoinLines(gf_point arcOrigin,					// Arc origin (unused if bArcJoin = false)
				gf_point *pLine1, gf_point *pLine2,	// two line vectors
				int nJoinType,						// join method (0: strait line, 1: arc, 2: crossing point)
				gf_point *pPoint, int *nPoints)		// point array and index of the last new point
{
	// Join two lines with an arc or at the crossing point then add the join points to the existing point array 'pPoint'.
	// The arc origin is at 'arcOrigin'
	//
	// Returns:	true if the two lines are crossed, false if not crossed
	
	bool bCrossed = false;
	bool bDirectConnect = false;
	
	int p = *nPoints;

	double x, y;
	int nFlag = nTwoLines(pLine1, pLine2, &x, &y);

	if (nFlag < 0)			// Parallel lines
	{
		bDirectConnect = true;
	}
	else if (nFlag == 3)	// Crossing lines
	{
		pPoint[p].x = x;
		pPoint[p].y = y;
		p++;
		bCrossed = true;
	}
	else if (nFlag == 1 || nFlag == 2)	// Do nothing
	{
		return false;
	}
	else					// Neither parallel or crossing
	{
		if (nJoinType == 0)				// Straight line (direct) join
		{
			bDirectConnect = true;
		}
		else if (nJoinType == 1)		// Arc join
		{
			double dCross = dCrossProduct(pLine1, pLine2);
			if (dCross == 0)
			{
				bDirectConnect = true;
			}
			else
			{
				// Join with an arc
				bool bCCW;
				if (dCross > 0)
					bCCW = true;
				else
					bCCW = false;
				int nArcPoints;
				gf_point *pArc = pDivideArc(arcOrigin, pLine1[1], pLine2[0], bCCW,
											LINE_JOINT_ARC_DIV_ANGLE, 2, LINE_JOINT_ARC_DIV_MAX_POINT, &nArcPoints);
				int i;
				for (i = 0; i < nArcPoints; i++)
				{
					pPoint[p].x = pArc[i].x;
					pPoint[p].y = pArc[i].y;
					p++;
				}
				FREE(pArc);
			}
		}
		else if (nJoinType == 2)		// Join at crossing point of line extensions
		{
			pPoint[p].x = pLine1[1].x;
			pPoint[p].y = pLine1[1].y;
			p++;
			pPoint[p].x = x;
			pPoint[p].y = y;
			p++;
			pPoint[p].x = pLine2[0].x;
			pPoint[p].y = pLine2[0].y;
			p++;
		}
	}
	
	if (bDirectConnect)
	{
		pPoint[p].x = pLine1[1].x;
		pPoint[p].y = pLine1[1].y;
		p++;
		pPoint[p].x = pLine2[0].x;
		pPoint[p].y = pLine2[0].y;
		p++;
	}
	
	*nPoints = p;

	return bCrossed;
}

int nCleanSelfCrossingLines(gf_point *pPath,		// List of points
							int nPoints)			// Number of points
{
	// returns a new number of points after clean up self crossing
	int nP = ABS(nPoints);
	int i, k, n;
	double x, y;
	
	int *skip = (int *) calloc(nP, sizeof(int));
	
	// Clean up self intersecting paths
	
	gf_point p[2], q[2];
	
	p[0].x = pPath[0].x;
	p[0].y = pPath[0].y;
	
	for (i = 1; i < nP; i++)
	{
		if (skip[i]) continue;
		p[1].x = pPath[i].x;
		p[1].y = pPath[i].y;
		
		k = i + 2;
		while (k < nP)
		{
			if (skip[k-1] == 0) break;
			k++;
		}
		
		q[0].x = pPath[k-1].x;
		q[0].y = pPath[k-1].y;
		for (; k < nP; k++)
		{
			if (skip[k]) continue;
			q[1].x = pPath[k].x;
			q[1].y = pPath[k].y;
			
			if (nTwoLines(p, q, &x, &y) == 3)
			{
				pPath[i].x = x;
				pPath[i].y = y;
				for (n = i + 1; n < k; n++)
				{
					skip[n] = 1;
				}
			}
			q[0].x = q[1].x;
			q[0].y = q[1].y;
		}
		p[0].x = p[1].x;
		p[0].y = p[1].y;
	}
	
	n = 0;
	for (i = 0; i < nP; i++)
	{
		if (!skip[i])
		{
			pPath[n].x = pPath[i].x;
			pPath[n].y = pPath[i].y;
			n++;
		}
	}
	
	return n;
}

bool bGetParallelLine(double xs, double ys, double xe, double ye,			// start end end points of the reference line
					  double *pxs, double *pys, double *pxe, double *pye,	// start and end points of the parallel line (output)
					  double distance,		// distance from the reference line
					  double extended,		// extend both ends of the parallel line
					  double epsilon)		// if the line length is shorter than this, return false with no parallel line output
{
	// Obtain a parallel line on the right hand side of the reference line vector
	double dx, dy, dd, c, s, xext, yext;

	dx = xe - xs;
	dy = ye - ys;
	dd = sqrt (dx * dx + dy * dy);
	if (dd < epsilon) return false;

	c = dx / dd;
	s = dy / dd;
	xext = extended * c;
	yext = extended * s;
	*pxs = (xs - xext) + distance * s;
	*pys = (ys - yext) - distance * c;
	*pxe = (xe + xext) + distance * s;
	*pye = (ye + yext) - distance * c;
	return true;
}

bool bGetParallelLine(gf_point *line, gf_point *rline,
					  double distance, double extended, double epsilon)
{
	return bGetParallelLine(line[0].x, line[0].y, line[1].x, line[1].y,
							&(rline[0].x), &(rline[0].y), &(rline[1].x), &(rline[1].y),
							distance, extended, epsilon);
}

bool bGetParallelLine(gf_line line, gf_line rline,
					  double distance, double extended, double epsilon)
{
	return bGetParallelLine(line.s.x, line.s.y, line.e.x, line.e.y,
							&(rline.s.x), &(rline.s.y), &(rline.e.x), &(rline.e.y),
							distance, extended, epsilon);
}


gf_point *pGetPathApart(gf_point *pLine,	// point list of the reference line
						int nPoints,		// number of points in the point list
						double dDistance,	// distance from the reference line (positive: on the right, negative: on the left)
						int nJoinType,		// line join method between two neighboring lines (0: straight line: 1: arc, 2: crossing point)
						int *nRetDivPoints)	// outputs number of points in the result point list
{
	// Determine a path (polyline) from the reference path apart a specified distance on the right or left.
	//

	double	dDIV_ANGLE = LINE_JOINT_ARC_DIV_ANGLE;
	int		nDIV_MAX = LINE_JOINT_ARC_DIV_MAX_POINT;

	gf_point prevSeg[2], newSeg[2];
	int i;

// 	double dEpsilon = fabs(dDistance) / 20;
	double dEpsilon = 1.0e-4 * fabs(dDistance);

	int p = 0;	// Actual point counter
	int q;

	// Allocate the result point array with enough number of points

	int nPointsMax = nPoints + ((nDIV_MAX-1) * nPoints);

	gf_point *pTemp = (gf_point *) malloc(nPointsMax * sizeof(gf_point));

	// Loop through all points in increasing index order
	for (i = 0; i < nPoints-1; i++)
	{
		// Ger a line segment parallel to the reference line
		if (!bGetParallelLine(&pLine[i], newSeg, dDistance, 0.0, dEpsilon)) continue;
		q = p;
		if (p == 0)
		{
			// keep the first point
			pTemp[p].x = newSeg[0].x;
			pTemp[p].y = newSeg[0].y;
			p++;
		}
		else
		{
			// joint the previous parallel line and the new parallel line
			bJoinLines(pLine[i], prevSeg, newSeg, nJoinType, pTemp, &p);
		}
		if (p != q)
		{
			prevSeg[0].x = newSeg[0].x;
			prevSeg[0].y = newSeg[0].y;
			prevSeg[1].x = newSeg[1].x;
			prevSeg[1].y = newSeg[1].y;
		}
	}
	if (p > 0)
	{
		// Keep the last point
		pTemp[p].x = prevSeg[1].x;
		pTemp[p].y = prevSeg[1].y;
		p++;
	}

	// Copy the results to a new point array

	gf_point *pResult = (gf_point *) malloc(p * sizeof(gf_point));

	for (i = 0; i < p; i++)
	{
		pResult[i].x = pTemp[i].x;
		pResult[i].y = pTemp[i].y;
	}
	FREE(pTemp);

	*nRetDivPoints = p;

	return pResult;
}

double dAngleFromEast(double dX1, double dY1,	// Starting point of the line
					  double dX2, double dY2)	// Ending point of the line
{
	// Determine angle of the line from east
	
	double dAngle = 0;
	
	// If the start and end points are the same, return the angle = 0
	if (dX1 == dX2 && dY1 == dY2) return dAngle;
	
	double dDeltaX = dX2 - dX1;
	double dDeltaY = dY2 - dY1;
	
	// vertical line
	if (fabs(dDeltaX) < 1.0e-20)
	{
		if (dDeltaY > 0)
			dAngle = 90;
		else
			dAngle = 270;
	}
	// non-vertical line
	else
	{
		dAngle = atan(dDeltaY / dDeltaX);
		dAngle *= 180.0 / ONEPI;
		if (dDeltaX < 0) dAngle += 180;
		while(dAngle < 0) {dAngle += 360;}		// Make sure the angle >= 0
		while(dAngle >= 360) {dAngle -= 360;}	// Make sure the angle < 360
	}
	
	return dAngle;
}

double safe_acos( double dVal )
{
	if( dVal < -1.0 )
		dVal = -1.0 ;
	else if( dVal > 1.0 )
		dVal = 1.0 ;

	return acos(dVal) ;
}

double dAngleTwoLines(double dXA1, double dYA1, double dXA2, double dYA2,	// Start & end points of line A
					  double dXB1, double dYB1, double dXB2, double dYB2)	// Start & end points of line B
{
	// Determine angle between two separate lines
	
	double dAngle = 0;
	
	// If the start and end points are the same, return the angle = 0
	if (dXA1 == dXA2 && dYA1 == dYA2) return dAngle;
	if (dXB1 == dXB2 && dYB1 == dYB2) return dAngle;
	
	// Get dot products and determinant of the two line vectors
	double dXA = dXA2 - dXA1;
	double dYA = dYA2 - dYA1;
	double dXB = dXB2 - dXB1;
	double dYB = dYB2 - dYB1;
	double dot = dXA * dXB + dYA * dYB;
	double det = sqrt(dXA * dXA + dYA * dYA) * sqrt(dXB * dXB + dYB * dYB);
	
	double dotdet = dot / det;

	dAngle = safe_acos(dotdet) * 180 / ONEPI;

	while(dAngle < 0) {dAngle += 360;}		// Make sure the angle >= 0
	while(dAngle >= 360) {dAngle -= 360;}	// Make sure the angle < 360

	return dAngle;
}

double dAngle3Points(double dX1, double dY1, double dX2, double dY2, double dX3, double dY3)
{
	// Determine angle between two lines formed by three consecutive points
	
	return dAngleTwoLines(dX1, dY1, dX2, dY2, dX3, dY3, dX2, dY2);
}