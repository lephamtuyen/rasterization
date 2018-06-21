//==============================================
//
// Functions Related to Contours and Point Lists
//
//==============================================
//
// dContourArea - calculates area of a contour defined in CW or CCW.
//				  The contour can be either convex or concave.
//
// bClockwise - returns true if clockwise, false if counter-clockwise
//
// pRemoveCollinearPoint - removes collinear points from the point list
// vRemoveCollinearPoint - removes collinear points from the point list
//
// pRemoveTinyHoles - removes tiny holes
//
// bContourInContour - returns true if a contour is inside of another contour
//
// nFindContour - finds an identical contour
//
// pRemoveRedundantContours - removes redundant contours
//
// vShiftContour  - translates (shifts) contour
// vRotateContour - rotates contour
// vMirrorContour - mirrors contour
//
// pBox2Contour		- converts a box to a contour
// pCircle2Contour	- converts a circle to a contour
// pEllipse2Contour - converts an ellipse to a contour
// pLine2Contour	- converts a line segment to a contour
// pTapered2Contour	- converts a tapered line segment to a contour
// pPath2Contour	- converts a line path to a contour
// pMoire2Contour	- converts a Moire circle to a set of contours
// pThermal2Contour	- converts a wagon wheel shape (thermal relief) to a set of contours
// pArc2Contour		- converts an arc to a contour
// pExpandContour	- expands or shrink contour
// -----------------------------------------------------------------------------
#include <stdlib.h>
#include <complex>
#include "gf_types.h"
#include "gfunc.h"

double dContourArea(double *x, double *y, int npoints)
{
	int    i;
	double area = 0;

	for (i = 0; i < npoints-1; i++)
	{
		area += x[i]   * y[i+1];
		area -= x[i+1] * y[i];
	}
	area += x[i] * y[0];
	area -= x[0] * y[i];
	if (area < 0) area = -area;
	area = area * 0.5;

	return area;
}

double dContourArea(gf_point *p, int npoints)
{
	int    i;
	double area = 0;

	for (i = 0; i < npoints-1; i++)
	{
		area += p[i].x   * p[i+1].y;
		area -= p[i+1].x * p[i].y;
	}
	area += p[i].x * p[0].y;
	area -= p[0].x * p[i].y;
	if (area < 0) area = -area;
	area = area * 0.5;

	return area;
}

bool bClockwise(double *x, double *y, int npoints)
{
	int    i;
	double sum = 0;

	for (i = 0; i < npoints-1; i++)
	{
		sum += x[i]   * y[i+1];
		sum -= x[i+1] * y[i];
	}
	sum += x[i] * y[0];
	sum -= x[0] * y[i];
	if (sum >= 0)
		return false;	// Counter clockwise
	else
		return true;	// Clockwise
}

bool bClockwise(gf_point *p, int npoints)
{
	int    i;
	double sum = 0;

	sum = 0;
	for (i = 0; i < npoints-1; i++) {
		sum += p[i].x   * p[i+1].y;
		sum -= p[i+1].x * p[i].y;
	}
	sum += p[i].x * p[0].y;
	sum -= p[0].x * p[i].y;
	if (sum >= 0)
		return false;	// Counter clockwise
	else
		return true;	// Clockwise
}

int nRemoveCollinearPoint(double *x, double *y, int npoints, bool bContour, double tolerance /* = 0.0 */)
{
	// bContour = true: the points define contour and the last line segment
	//					connects the first point and the last point. Therefore
	//					check collinearity between the first two points and
	//					the last two points.
	//			 false: the points define polyline. Therefore do not remove
	//					the first point and/or the last point.

	// returns a new number of points
	int nCountNew;;

	int    i, k, nKeep;
	bool   *keep;

	if (npoints < 2) return npoints ;

	// Allocate an array to save a flag to keep or discard
	keep = (bool *) malloc(npoints * sizeof(bool));

	// Remove duplicated points

	double eps = fabs(tolerance);

	keep[0] = true;
	for (i = 1; i < npoints; i++)
	{
		keep[i] = true;
		if (fabs(x[i] - x[i-1]) <= eps &&
			fabs(y[i] - y[i-1]) <= eps) keep[i] = false;
	}
	if (fabs(x[0] - x[i-1]) <= eps &&
		fabs(y[0] - y[i-1]) <= eps) keep[0] = false;

	// Find the first two points to keep
	nKeep = 0;
	k = -1;
	for (i = 0; i < npoints; i++)
	{
		if (keep[i])
		{
			nKeep++;
			if (k >= 0) break;	// When it breaks, nKeep = 2 and i = 2nd point to keep
			k = i;				// k is the first point to keep
		}
	}

	// check for collinear points only if the 2nd point is not the last point in the list
	if (i < npoints)
	{
		int first[2];	// starting and ending points of the first checked line segment
		int last[2];	// starting and ending points of the last checked line segment

		first[0] = k;
		first[1] = i;
		last[0] = i;
		last[1] = i;

		double dx1, dy1, dx2, dy2;

		dx1 = x[i] - x[k];
		dy1 = y[i] - y[k];
		k = i;
		i++;
		for (; i < npoints; i++)
		{
			if (keep[i])
			{
				dx2 = x[i] - x[k];
				dy2 = y[i] - y[k];

				if ((dx1 * dy2 - dx2 * dy1) == 0.0)		// Collinear checking
				{
					keep[k] = false;	// Remove collinear point
				}
				else					// Keep non-collinear point
				{
					dx1 = dx2;
					dy1 = dy2;
					last[0] = last[1];
					last[1] = i;
					nKeep++;
				}
				k = i;
			}
		}

		if (nKeep == 3)
		{
			// If there are only 3 points to keep, check the 2nd point
			dx1 = x[first[1]] - x[first[0]];
			dy1 = y[first[1]] - y[first[0]];
			dx2 = x[last[1]] - x[first[1]];
			dy2 = y[last[1]] - y[first[1]];
			if ((dx1 * dy2 - dx2 * dy1) == 0.0) keep[first[1]] = false;
		}
		else if (nKeep > 3 && bContour)
		{
			// Check the first point for collinearity
			dx1 = x[first[1]] - x[first[0]];
			dy1 = y[first[1]] - y[first[0]];
			dx2 = x[first[0]] - x[last[1]];
			dy2 = y[first[0]] - y[last[1]];
			if ((dx1 * dy2 - dx2 * dy1) == 0.0) keep[first[0]] = false;
			// Check the last point for collinearity
			dx1 = x[last[1]] - x[last[0]];
			dy1 = y[last[1]] - y[last[0]];
			dx2 = x[first[0]] - x[last[1]];
			dy2 = y[first[0]] - y[last[1]];
			if ((dx1 * dy2 - dx2 * dy1) == 0.0) keep[last[1]] = false;
		}
	}

	// Construct a new array to return
	nCountNew = 0;
	for (i = 0; i < npoints; i++)
	{
		if (keep[i])
		{
			x[nCountNew] = x[i];
			y[nCountNew] = y[i];
			nCountNew++;
		}
	}
	FREE(keep);
	return nCountNew;
}

int nRemoveCollinearPoint(gf_point *p, int npoints, bool bContour, double tolerance /* = 0.0 */)
{
	// bContour = true: the points define contour and the last line segment
	//					connects the first point and the last point. Therefore
	//					check collinearity between the first two points and
	//					the last two points.
	//			 false: the points define polyline. Therefore do not remove
	//					the first point and/or the last point.

	// returns a new number of points
	int nCountNew;;

	int    i, k, nKeep;
	bool   *keep;

	if (npoints < 2) return npoints ;

	// Allocate an array to save a flag to keep or discard
	keep = (bool *) malloc(npoints * sizeof(bool));

	// Remove duplicated points

	double eps = fabs(tolerance);

	keep[0] = true;
	for (i = 1; i < npoints; i++)
	{
		keep[i] = true;
		if (fabs(p[i].x - p[i-1].x) <= eps &&
			fabs(p[i].y - p[i-1].y) <= eps) keep[i] = false;
	}
	if (fabs(p[0].x - p[i-1].x) <= eps &&
		fabs(p[0].y - p[i-1].y) <= eps) keep[0] = false;

	// Find the first two points to keep
	nKeep = 0;
	k = -1;
	for (i = 0; i < npoints; i++)
	{
		if (keep[i])
		{
			nKeep++;
			if (k >= 0) break;	// When it breaks, nKeep = 2 and i = 2nd point to keep
			k = i;				// k is the first point to keep
		}
	}

	// check for collinear points only if the 2nd point is not the last point in the list
	if (i < npoints)
	{
		int first[2];	// starting and ending points of the first checked line segment
		int last[2];	// starting and ending points of the last checked line segment

		first[0] = k;
		first[1] = i;
		last[0] = i;
		last[1] = i;

		double dx1, dy1, dx2, dy2;

		dx1 = p[i].x - p[k].x;
		dy1 = p[i].y - p[k].y;
		k = i;
		i++;
		for (; i < npoints; i++)
		{
			if (keep[i])
			{
				dx2 = p[i].x - p[k].x;
				dy2 = p[i].y - p[k].y;

				if ((dx1 * dy2 - dx2 * dy1) == 0.0)		// Collinear checking
				{
					keep[k] = false;	// Remove collinear point
				}
				else					// Keep non-collinear point
				{
					dx1 = dx2;
					dy1 = dy2;
					last[0] = last[1];
					last[1] = i;
					nKeep++;
				}
				k = i;
			}
		}

		if (nKeep == 3)
		{
			// If there are only 3 points to keep, check the 2nd point
			dx1 = p[first[1]].x - p[first[0]].x;
			dy1 = p[first[1]].y - p[first[0]].y;
			dx2 = p[last[1]].x - p[first[1]].x;
			dy2 = p[last[1]].y - p[first[1]].y;
			if ((dx1 * dy2 - dx2 * dy1) == 0.0) keep[first[1]] = false;
		}
		else if (nKeep > 3 && bContour)
		{
			// Check the first point for collinearity
			dx1 = p[first[1]].x - p[first[0]].x;
			dy1 = p[first[1]].y - p[first[0]].y;
			dx2 = p[first[0]].x - p[last[1]].x;
			dy2 = p[first[0]].y - p[last[1]].y;
			if ((dx1 * dy2 - dx2 * dy1) == 0.0) keep[first[0]] = false;
			// Check the last point for collinearity
			dx1 = p[last[1]].x - p[last[0]].x;
			dy1 = p[last[1]].y - p[last[0]].y;
			dx2 = p[first[0]].x - p[last[1]].x;
			dy2 = p[first[0]].y - p[last[1]].y;
			if ((dx1 * dy2 - dx2 * dy1) == 0.0) keep[last[1]] = false;
		}
	}

	// Construct a new array to return
	nCountNew = 0;
	for (i = 0; i < npoints; i++)
	{
		if (keep[i])
		{
			p[nCountNew].x = p[i].x;
			p[nCountNew].y = p[i].y;
			nCountNew++;
		}
	}
	FREE(keep);
	return nCountNew;
}

gf_point *pRemoveCollinearPoint(gf_point *point, int *nPoints, bool bContour, double tolerance /* = 0.0 */)
{
	// bContour = true: the points define contour and the last line segment
	//					connects the first point and the last point. Therefore
	//					check collinearity between the first two points and
	//					the last two points.
	//			 false: the points define polyline. Therefore do not remove
	//					the first point and/or the last point.

	gf_point *result;

	int    nP, i, k, nKeep;
	bool   *keep;

	nP = ABS(*nPoints);
	if (nP < 2) return point ;

	// Allocate an array to save a flag to keep or discard
	keep = (bool *) malloc(nP * sizeof(bool));

	// Remove duplicated points

	double eps = fabs(tolerance);

	keep[0] = true;
	for (i = 1; i < nP; i++)
	{
		keep[i] = true;
		if (fabs(point[i].x - point[i-1].x) <= eps &&
			fabs(point[i].y - point[i-1].y) <= eps) keep[i] = false;
	}
	if (fabs(point[0].x - point[i-1].x) <= eps &&
		fabs(point[0].y - point[i-1].y) <= eps) keep[0] = false;

	// Find the first two points to keep
	nKeep = 0;
	k = -1;
	for (i = 0; i < nP; i++)
	{
		if (keep[i])
		{
			nKeep++;
			if (k >= 0) break;	// When it breaks, nKeep = 2 and i = 2nd point to keep
			k = i;				// k is the first point to keep
		}
	}

	// check for collinear points only if the 2nd point is not the last point in the list
	if (i < nP)
	{
		int first[2];	// starting and ending points of the first checked line segment
		int last[2];	// starting and ending points of the last checked line segment

		first[0] = k;
		first[1] = i;
		last[0] = i;
		last[1] = i;

		double dx1, dy1, dx2, dy2;

		dx1 = point[i].x - point[k].x;
		dy1 = point[i].y - point[k].y;
		k = i;
		i++;
		for (; i < nP; i++)
		{
			if (keep[i])
			{
				dx2 = point[i].x - point[k].x;
				dy2 = point[i].y - point[k].y;

				if ((dx1 * dy2 - dx2 * dy1) == 0.0)		// Collinear checking
				{
					keep[k] = false;	// Remove collinear point
				}
				else					// Keep non-collinear point
				{
					dx1 = dx2;
					dy1 = dy2;
					last[0] = last[1];
					last[1] = i;
					nKeep++;
				}
				k = i;
			}
		}

		if (nKeep == 3)
		{
			// If there are only 3 points to keep, check the 2nd point
			dx1 = point[first[1]].x - point[first[0]].x;
			dy1 = point[first[1]].y - point[first[0]].y;
			dx2 = point[last[1]].x - point[first[1]].x;
			dy2 = point[last[1]].y - point[first[1]].y;
			if ((dx1 * dy2 - dx2 * dy1) == 0.0) keep[first[1]] = false;
		}
		else if (nKeep > 3 && bContour)
		{
			// Check the first point for collinearity
			dx1 = point[first[1]].x - point[first[0]].x;
			dy1 = point[first[1]].y - point[first[0]].y;
			dx2 = point[first[0]].x - point[last[1]].x;
			dy2 = point[first[0]].y - point[last[1]].y;
			if ((dx1 * dy2 - dx2 * dy1) == 0.0) keep[first[0]] = false;
			// Check the last point for collinearity
			dx1 = point[last[1]].x - point[last[0]].x;
			dy1 = point[last[1]].y - point[last[0]].y;
			dx2 = point[first[0]].x - point[last[1]].x;
			dy2 = point[first[0]].y - point[last[1]].y;
			if ((dx1 * dy2 - dx2 * dy1) == 0.0) keep[last[1]] = false;
		}
	}

	// Count number of points to keep
	k = 0;
	for (i = 0; i < nP; i++) {
		if (keep[i]) k++;
	}

	// Output number of points and keep the +/- sign of the original number
	if (*nPoints > 0)
		*nPoints = k;
	else
		*nPoints = -k;

	// Construct a new array to return
	result = (gf_point *) malloc(k * sizeof(gf_point));
	k = 0;
	for (i = 0; i < nP; i++)
	{
		if (keep[i]) {
			result[k].x = point[i].x;
			result[k].y = point[i].y;
			k++;
		}
	}
	FREE(keep);
	FREE(point);
	return result;
}

void vRemoveCollinearPoint(gf_polygon *p, double tolerance /* = 0.0 */)
{
	// Remove collinear points within each contour in the polygon
	int c;
	for (c = 0; c < p->nContours; c++)
	{
		p->contour[c].point
			= pRemoveCollinearPoint(p->contour[c].point, &(p->contour[c].nPoints), true, tolerance);
	}
}

// Set and get criteria to determine tiny contour

double	TINY_HOLE_AREA = 0.0;
int		TINY_HOLE_POINTS = 36;

void vSetTinyAreaSize(double tinyArea, int nPoints)
{
	TINY_HOLE_AREA = tinyArea;
	TINY_HOLE_POINTS = nPoints;
}

int nGetTinyAreaPoints()
{
	return TINY_HOLE_POINTS;
}

double dGetTinyAreaSize()
{
	return TINY_HOLE_AREA;
}

gf_polygon *pRemoveTinyHoles(gf_polygon *p)
{
	// Remove tiny holes

	if (TINY_HOLE_AREA <= 0) return p;	//  Return the original polygon as it was

	int c;

	// Check if there is at least one hole
	bool bFound = false;
	for (c = 0; c < p->nContours; c++)
	{
		if (p->hole[c])
		{
			bFound = true;
			break;
		}
	}
	if (!bFound) return p;		// No holes. Return the original polygon as it was

	// Search for tiny holes
	double bigEnough = TINY_HOLE_AREA * 100;
	double dx, dy;
	int nRemove = 0;

	for (c = 0; c < p->nContours; c++)
	{
		if (p->hole[c] && p->contour[c].nPoints < TINY_HOLE_POINTS)
		{
			dx = p->box[c].ur.x - p->box[c].ll.x;
			dy = p->box[c].ur.y - p->box[c].ll.y;
			if (dx * dy > bigEnough) continue;

			if (dContourArea(p->contour[c].point, p->contour[c].nPoints) < TINY_HOLE_AREA)
			{
				// Found a tiny hole. Flag it to remove it later.
				p->contour[c].nPoints = -abs(p->contour[c].nPoints);
				nRemove++;

				// Check solid contours that immediately follow this hole
				int s = c + 1;
				for (; s < p->nContours; s++)
				{
					if (p->hole[s]) break;

					if (PointInContour(	p->contour[s].point[0].x, p->contour[s].point[0].y,
										p->contour[c].nPoints, p->contour[c].point))
					{
						// Found a tiny solid that is inside of the tiny hole. Flag it to remove it later.
						p->contour[c].nPoints = -abs(p->contour[c].nPoints);
						nRemove++;
					}
					else
					{
						s++;
						break;
					}
				}
				c = s - 1;
			}
		}
	}

	if (nRemove == 0) return p;		// No tiny holes. Just return the original polygon.

	// Reconstruct the polygon within the tiny holes

	int n, v;

	n = 0;
	for (c = 0; c < p->nContours; c++)
	{
		if (p->contour[c].nPoints > 2) n++;
	}

	gf_polygon *result = gfAllocPolygon(n);

	n = 0;
	for (c = 0; c < p->nContours; c++)
	{
		if (p->contour[c].nPoints <= 2) continue;

		result->hole[n] = p->hole[c];

		gfAllocContour(result, n, p->contour[c].nPoints);
		result->contour[n].nPoints = p->contour[c].nPoints;
		for (v = 0; v < result->contour[n].nPoints; v++)
		{
			result->contour[n].point[v].x = p->contour[c].point[v].x;
			result->contour[n].point[v].y = p->contour[c].point[v].y;
		}

		if (p->box)
		{
			result->box[n].ll.x = p->box[c].ll.x;
			result->box[n].ll.y = p->box[c].ll.y;
			result->box[n].ur.x = p->box[c].ur.x;
			result->box[n].ur.y = p->box[c].ur.y;
		}
		n++;
	}
	gfFreePolygon(p);

	return result;
}

bool bContoursIntersect(gf_point *pContour1, int nPoints1, gf_point *pContour2, int nPoints2, bool bCheckBoundingBox /* = true */)
{
	// Check whether the two contours intersect at all.

	if (pContour1 == NULL ||
		pContour2 == NULL) return false;

	if (bCheckBoundingBox)
	{
		gf_box box1 = gfGetBoundingBox(pContour1, nPoints1);
		gf_box box2 = gfGetBoundingBox(pContour2, nPoints2);
		if (!gfBoxesOverlap(box1, box2)) return false;
	}

	for (int v = 0; v < nPoints1; v++)
	{
		if (PointInContour(pContour1[v].x, pContour1[v].y, nPoints2, pContour2)) return true;
	}

	for (int v = 0; v < nPoints2; v++)
	{
		if (PointInContour(pContour2[v].x, pContour2[v].y, nPoints1, pContour1)) return true;
	}

	for (int i = 0; i < nPoints1; i++)
	{
		for (int j = 0; j < nPoints2; j++)
		{
			double x, y;

			if (bTwoLinesCross(	pContour1[i].x, pContour1[i].y, pContour1[(i + 1) % nPoints1].x, pContour1[(i + 1) % nPoints1].y,
								pContour2[j].x, pContour2[j].y, pContour2[(j + 1) % nPoints2].x, pContour2[(j + 1) % nPoints2].y,
								&x, &y)
				&& PointOnLine(x, y, pContour1[i].x, pContour1[i].y, pContour1[(i + 1) % nPoints1].x, pContour1[(i + 1) % nPoints1].y)
				&& PointOnLine(x, y, pContour2[j].x, pContour2[j].y, pContour2[(j + 1) % nPoints2].x, pContour2[(j + 1) % nPoints2].y))
				return true;
		}
	}
	return false;
}

bool bContoursIntersect(gf_point *pContour1, int nPoints1, gf_point *pContour2, int nPoints2,
						gf_box box1, gf_box box2)
{
	if (pContour1 == NULL ||
		pContour2 == NULL) return false;

	if (!gfBoxesOverlap(box1, box2)) return false;

	return bContoursIntersect(pContour1, nPoints1, pContour2, nPoints2, false);
}

bool bContourInContour(gf_point *pInside, int nInside, gf_point *pOutside, int nOutside)
{
	for(int i = 0; i < nInside; i++)
	{
		if(!PointInContour(pInside[i].x,
						   pInside[i].y,
						   nOutside,
						   pOutside))
			return false;
	}

	return true;
}

bool bContourInContour(gf_polygon *pInside, int nInID, gf_polygon *pOutside, int nOutID)
{
	// Check whether the (nInID)th contour of the polygon 'pInside' is completely
	// inside of the (nOutID)th contour of the polygon 'pOutside'

	//////////////////////////////////////////////////////////////////////////
	// pInside->contour[nInID].point 가 NULL 로 처리되 죽습니다.
	// 일단 예외 처리 하였습니다. chahn : 10.06.22
	if(NULL == pInside->contour[nInID].point 
		|| NULL == pOutside->contour[nOutID].point) return false;
	//////////////////////////////////////////////////////////////////////////

	gf_polygon *pResult;

	if (!PointInContour(pInside->contour[nInID].point[0].x,
						pInside->contour[nInID].point[0].y,
						pOutside->contour[nOutID].nPoints,
						pOutside->contour[nOutID].point)) return false;

	pResult = gfAllocPolygon(0);

	int nPointsInside = pInside->contour[0].nPoints;
	pInside->contour[0].nPoints = abs(pInside->contour[0].nPoints);

	int nPointsOutside = pOutside->contour[0].nPoints;
	pOutside->contour[0].nPoints = abs(pOutside->contour[0].nPoints);

	gfPolygonBool(GF_SUB, pInside, pOutside, pResult);

	pInside->contour[0].nPoints = nPointsInside;
	pOutside->contour[0].nPoints = nPointsOutside;

	int nContours = pResult->nContours;
	gfFreePolygon(pResult);

	if (nContours == 0)
		return true;
	else
		return false;
}

bool bContourInContour(gf_polygon *pInside, gf_polygon *pOutside)
{
	// Check whether the first contour of the polygon 'pInside' is completely
	// inside of the first contour of the polygon 'pOutside'

	return bContourInContour(pInside, 0, pOutside, 0);
}

int nFindContour(gf_contour contour, bool bHole,
				 gf_polygon *pPoly, int nStart, int nEnd,
				 bool bCheckType)
{
	// Find a contour in 'pPoly' that is identical to the contour 'contour'
	// then return the contour number if found. Otherwise, return -1.
	// Contours in 'pPoly' will be checked from contours number 'nStart' to 'nEnd'.
	// Check contour type (hole or dark) too if bCheckType is true.

	int nPoints = abs(contour.nPoints);

	int c, p;
	for (c = nStart; c <= nEnd; c++)
	{
		if (nPoints != pPoly->contour[c].nPoints) continue;		// point counts are different
		if (bCheckType && bHole != pPoly->hole[c]) continue;	// contour types are different

		// check coordinates of every point
		for (p = 0; p < nPoints; p++)
		{
			if (contour.point[p].x != pPoly->contour[c].point[p].x ||
				contour.point[p].y != pPoly->contour[c].point[p].y)
				break;
		}
		if (p == nPoints) return c;		// return the contour number if all points are the same
	}
	return -1;							// return -1 if there is no matching contour
}

int nFindContour(gf_polygon *pFind, int nContour,
				 gf_polygon *pPoly, int nStart, int nEnd,
				 bool bCheckType)
{
	// Find a contour in 'pPoly' that is identical to the contour number 'nContour' in 'pFind'
	// then return the contour number if found. Otherwise, return -1.
	// Contours in 'pPoly' will be checked from contours number 'nStart' to 'nEnd'.
	// Check contour type (hole or dark) too if bCheckType is true.

	return nFindContour(pFind->contour[nContour], pFind->hole[nContour],
						pPoly, nStart, nEnd, bCheckType);
}

gf_polygon *pRemoveRedundantContours(gf_polygon *pPoly)
{
	// Remove redundant contours that appear in the same polygon

	// Return the original polygon if only one contour in the polygon
	int nContours = abs(pPoly->nContours);
	if (nContours < 2) return pPoly;

	int c, p;
	int nKeep = 1;

	// Loop through all contours in the polygon
	for (c = 1; c < nContours; c++)
	{
		if (nFindContour(pPoly, c, pPoly, 0, c-1, true) >= 0)
		{
			// If the same contour is found, flag it to remove later
			pPoly->contour[c].nPoints = -abs(pPoly->contour[c].nPoints);
		}
		else
		{
			nKeep++;
		}
	}

	// Return the original polygon if no contour is duplicated
	if (nKeep == nContours) return(pPoly);

	gf_polygon *pResult = gfAllocPolygon(nKeep);

	// Copy the keep contours to the result polygon

	nKeep = 0;
	for (c = 0; c < nContours; c++)
	{
		if (pPoly->contour[c].nPoints > 0)
		{
			pResult->hole[nKeep] = pPoly->hole[c];
			gfAllocContour(pResult, nKeep, abs(pPoly->contour[c].nPoints));
			pResult->contour[nKeep].nPoints = pPoly->contour[c].nPoints;
			for (p = 0; p < abs(pPoly->contour[c].nPoints); p++)
			{
				pResult->contour[nKeep].point[p].x = pPoly->contour[c].point[p].x;
				pResult->contour[nKeep].point[p].y = pPoly->contour[c].point[p].y;
			}
			if (pPoly->box)
			{
				pResult->box[nKeep].ll.x = pPoly->box[c].ll.x;
				pResult->box[nKeep].ll.y = pPoly->box[c].ll.y;
				pResult->box[nKeep].ur.x = pPoly->box[c].ur.x;
				pResult->box[nKeep].ur.y = pPoly->box[c].ur.y;
			}
			nKeep++;
		}
	}

	gfFreePolygon(pPoly);

	return pResult;
}

void vShiftContour(gf_polygon *pPoly,
				   int iContour, int nContours,			// 'nContours' starting from the contour number 'iContour' of 'pPoly'
				   double dShiftX, double dShiftY)		// shift amount
{
	if (dShiftX == 0 && dShiftY == 0) return;
	
	int c, i;
	
	for (c = iContour; c < iContour + nContours; c++)
	{
		int nPoints = pPoly->contour[c].nPoints;
		
		for (i = 0; i < nPoints; i++)
		{
			pPoly->contour[c].point[i].x += dShiftX;
			pPoly->contour[c].point[i].y += dShiftY;
		}
	}
}

void vShiftContour(gf_point *pPoint, int nPoints,		// point array and number of points
				   double dShiftX, double dShiftY)		// shift amount
{
	if (dShiftX == 0 && dShiftY == 0) return;

	int i;
	for (i = 0; i < nPoints; i++)
	{
		pPoint[i].x += dShiftX;
		pPoint[i].y += dShiftY;
	}
}

void vRotateContour(gf_polygon *pPoly,
					int iContour, int nContours,	// 'nContours' starting from the contour number 'iContour' of 'pPoly'
					double dRotDegreeCCW,			// CCW rotation angle in degree
					double dRefX, double dRefY)		// rotates with respect to this reference point
{
	if (fmod(dRotDegreeCCW, 360) == 0.0) return;

	int c, i;
	double dAngle = dRotDegreeCCW * ONEPI / 180;
	double dCos = cos(dAngle);
	double dSin = sin(dAngle);
	double dX, dY;
	
	for (c = iContour; c < iContour + nContours; c++)
	{
		int nPoints = pPoly->contour[c].nPoints;
		
		for (i = 0; i < nPoints; i++)
		{
			dX = pPoly->contour[c].point[i].x - dRefX;
			dY = pPoly->contour[c].point[i].y - dRefY;
			pPoly->contour[c].point[i].x = dRefX + (dX * dCos - dY * dSin);
			pPoly->contour[c].point[i].y = dRefY + (dX * dSin + dY * dCos);
		}
	}
}

void vRotateContour(gf_point *pPoint,				// point coordinates
					int nPoints,					// number of points
					double dRotDegreeCCW,			// rotation angle in CCW degree
					double dXref, double dYref)		// rotation reference coordinates
{
	// Rotates the contour 'dRotDegree' CCW with respect to (dXref, dYref)
	
	if (nPoints <= 0) return;

	// If the rotation angle is a multiple of 360 degrees, no need to rotate
	if (fmod(dRotDegreeCCW, 360) == 0.0) return;

	int n;
	double dAngle = ONEPI * dRotDegreeCCW / 180.0;
	double dCos = cos(dAngle);
	double dSin = sin(dAngle);
	
	double dX, dY;
	
	for (n = 0; n < nPoints; n++)
	{
		dX = pPoint[n].x - dXref;
		dY = pPoint[n].y - dYref;
		pPoint[n].x = dXref + dCos * dX - dSin * dY;
		pPoint[n].y = dYref + dSin * dX + dCos * dY;
	}
}

void vRotateContour(gf_point *pPoint,				// point coordinates
					int nPoints,					// number of points
					double dRotDegreeCCW,			// rotation angle in CCW degree
					int nRefPoint)					// rotation reference point (if -1, use geometric center)
{
	// Rotates the contour 'dRotDegree' CCW with respect to the point number 'nRefPoint'
	// If the nRefPoint = -1, rotate with respect to the geometric center of the contour

	if (nPoints <= 0) return;
	
	// If the rotation angle is a multiple of 360 degrees, no need to rotate
	if (fmod(dRotDegreeCCW, 360) == 0.0) return;

	// Rotate with respect to the point 'nRefPoint'
	if (nRefPoint >= 0 && nRefPoint < nPoints)
	{
		vRotateContour(pPoint, nPoints, dRotDegreeCCW, pPoint[nRefPoint].x, pPoint[nRefPoint].y);
	}
	else if (nRefPoint == -1)
	{
		int n;
		double dX = 0;
		double dY = 0;
		for (n = 0; n < nPoints; n++)
		{
			dX += pPoint[n].x;
			dY += pPoint[n].y;
		}
		dX /= nPoints;
		dY /= nPoints;
		vRotateContour(pPoint, nPoints, dRotDegreeCCW, dX, dY);
	}
}

void vMirrorContour(gf_polygon *pPoly,				// point coordinates
					int iContour, int nContours,	// 'nContours' starting from the contour number 'iContour' of 'pPoly'
					int nMirror,					// mirror axis ('X' or 'Y')
					double dXref, double dYref)		// mirror reference coordinates
{
	// Mirrors the contour in 'X' or 'Y' with respect to (dXref, dYref)

	// If the mirror axis is not X or Y, no change is needed
	if (!(nMirror == 'X' || nMirror == 'x' || nMirror == 'Y' || nMirror == 'y')) return;
	
	int c;
	for (c = 0; c < pPoly->nContours; c++)
	{
		// Mirror each point with respect to the reference point
		int n;
		for (n = 0; n < pPoly->contour[c].nPoints; n++)
		{
			if (nMirror == 'X' || nMirror == 'x')
			{
				pPoly->contour[c].point[n].x = dXref - (pPoly->contour[c].point[n].x - dXref);
			}
			else if (nMirror == 'Y' || nMirror == 'y')
			{
				pPoly->contour[c].point[n].y = dYref - (pPoly->contour[c].point[n].y - dYref);
			}
		}
	}
}

void vMirrorContour(gf_point *pPoint,				// point coordinates
					int nPoints,					// number of points
					int nMirror,					// mirror axis ('X' or 'Y')
					double dXref, double dYref)		// mirror reference coordinates
{
	// Mirrors the contour in 'X' or 'Y' with respect to (dXref, dYref)
	
	if (nPoints <= 0) return;
	
	// If the mirror axis is not X or Y, no change is needed
	if (!(nMirror == 'X' || nMirror == 'x' || nMirror == 'Y' || nMirror == 'y')) return;
	
	// Mirror each point with respect to the reference point
	int n;
	for (n = 0; n < nPoints; n++)
	{
		if (nMirror == 'X' || nMirror == 'x')
		{
			pPoint[n].x = dXref - (pPoint[n].x - dXref);
		}
		else if (nMirror == 'Y' || nMirror == 'y')
		{
			pPoint[n].y = dYref - (pPoint[n].y - dYref);
		}
	}
}

void vMirrorContour(gf_point *pPoint,				// point coordinates
					int nPoints,					// number of points
					int nMirror,					// mirror axis ('X' or 'Y')
					int nRefPoint)					// mirror reference point (if -1, use geometric center)
{
	// Mirrors the contour in 'X' or 'Y' with respect to the reference point number 'nRefPoint'
	// If the nRefPoint = -1, mirror with respect to the geometric center of the contour
	
	if (nPoints <= 0) return;
	
	// If the mirror axis is not X or Y, no change is needed
	if (!(nMirror == 'X' || nMirror == 'x' || nMirror == 'Y' || nMirror == 'y')) return;
	
	// Mirror each point with respect to the reference point
	if (nRefPoint >= 0 && nRefPoint < nPoints)
	{
		vMirrorContour(pPoint, nPoints, nMirror, pPoint[nRefPoint].x, pPoint[nRefPoint].y);
	}
	else if (nRefPoint == -1)
	{
		int n;
		double dX = 0;
		double dY = 0;
		for (n = 0; n < nPoints; n++)
		{
			dX += pPoint[n].x;
			dY += pPoint[n].y;
		}
		dX /= nPoints;
		dY /= nPoints;
		vMirrorContour(pPoint, nPoints, nMirror, dX, dY);
	}
}

gf_point *pBox2Contour(double dXcenter, double dYcenter,	// box center
					   double dWidth, double dHeight,		// box width and height
					   double dRotDegree,					// rotation angle (degree in CCW)
					   int nMirror)							// mirror flag ('X', 'Y' or 'N')
{
	// Converts a box defined by center, width, height, rotation and mirror to a contour.
	// This function rotates the box then mirrors.
	
	gf_point *pResult = (gf_point *) malloc(4 * sizeof(gf_point));

	// convert to contour before any orientation
	pResult[0].x = 0.5 * dWidth;
	pResult[1].x = 0.5 * dWidth;
	pResult[2].x = -0.5 * dWidth;
	pResult[3].x = -0.5 * dWidth;
	pResult[0].y = -0.5 * dHeight;
	pResult[1].y = 0.5 * dHeight;
	pResult[2].y = 0.5 * dHeight;
	pResult[3].y = -0.5 * dHeight;

	int n;

	// If the rotation angle is a multiple of 180 degrees, no need to rotate or mirror
	if (fmod(dRotDegree, 180) == 0.0)
	{
		for (n = 0; n < 4; n++)
		{
			pResult[n].x += dXcenter;
			pResult[n].y += dYcenter;
		}
	}
	// Rotate and mirror
	else
	{
		double dAngle = ONEPI * dRotDegree / 180.0;
		double dCos = cos(dAngle);
		double dSin = sin(dAngle);
		
		double dX, dY;
		
		for (n = 0; n < 4; n++)
		{
			dX = pResult[n].x;
			dY = pResult[n].y;
			pResult[n].x = dCos * dX - dSin * dY;
			pResult[n].y = dSin * dX + dCos * dY;
			if (nMirror == 'X' || nMirror == 'x')
			{
				pResult[n].x = dXcenter - pResult[n].x;
				pResult[n].y = dYcenter + pResult[n].y;
			}
			else if (nMirror == 'Y' || nMirror == 'y')
			{
				pResult[n].x = dXcenter + pResult[n].x;
				pResult[n].y = dYcenter - pResult[n].y;
			}
			else
			{
				pResult[n].x = dXcenter + pResult[n].x;
				pResult[n].y = dYcenter + pResult[n].y;
			}
		}
	}
	return pResult;
}

gf_point *pBox2Contour(double dXcenter, double dYcenter,	// box center
					   double dWidth, double dHeight,		// box width and height
					   int nMirror,							// mirror flag ('X', 'Y' or 'N')
					   double dRotDegree)					// rotation angle (degree in CCW)
{
	// Converts a box defined by center, width, height, rotation and mirror to a contour.
	// This function mirrors the box then rotates.
	
	// Because the box is symmetric along X and Y axes before any orientation,
	// mirror first is the same as no mirror. Therefore, mirror->rotate is the same as rotate->no mirror.

	return pBox2Contour(dXcenter, dYcenter, dWidth, dHeight, dRotDegree, 'N');
}

gf_point *pBox2Contour(double dXcenter, double dYcenter,	// box center
					   double dWidth, double dHeight,		// box width and height
					   double dRotDegree)					// rotation angle (degree in CCW)
{
	// Converts a box defined by center, width, height, rotation to a contour.
	
	return pBox2Contour(dXcenter, dYcenter, dWidth, dHeight, dRotDegree, 'N');
}

gf_point *pBox2Contour(double dXcenter, double dYcenter,	// box center
					   double dWidth, double dHeight)		// box width and height
{
	// Converts a box defined by center, width and height to a contour.

	return pBox2Contour(dXcenter, dYcenter, dWidth, dHeight, 0.0, 'N');
}

gf_point *pRoundedBox2Contour(double dXcenter, double dYcenter,	// box center
							  double dWidth, double dHeight,	// box width and height
							  double dChamferRadius,			// Edge chamfer radius
							  int nChamferDivPoints,			// Number of division points on each chamfer arc
							  double dRotDegree,				// rotation angle (degree in CCW)
							  int nMirror)						// mirror flag ('X', 'Y' or 'N')
{
	// Converts a rounded box defined by center, width, height, chamfer radius, chamfer division points,
	// rotation and mirror to a contour.
	// This function rotates the rounded box then mirrors.

	gf_point *pResult = (gf_point *) malloc(4 * nChamferDivPoints * sizeof(gf_point));
	gf_point *pChamfer;
	int nPoints = 0;

	// convert to contour before any orientation

	double dX1, dY1, dX2, dY2;
	int i;

	dX1 = 0.5 * dWidth;
	dY1 = 0.5 * dHeight - dChamferRadius;
	dX2 = 0.5 * dWidth - dChamferRadius;
	dY2 = 0.5 * dHeight;
	pChamfer = pDivideArc(dX2, dY1, dX1, dY1, dX2, dY2, true, nChamferDivPoints);

	for (i = 0; i < nChamferDivPoints; i++)
	{
		pResult[nPoints].x = pChamfer[i].x;
		pResult[nPoints].y = pChamfer[i].y;
		nPoints++;
	}
	gfFreePoints(pChamfer);
	
	dX1 = -0.5 * dWidth + dChamferRadius;
	dY1 = 0.5 * dHeight;
	dX2 = -0.5 * dWidth;
	dY2 = 0.5 * dHeight - dChamferRadius;
	pChamfer = pDivideArc(dX1, dY2, dX1, dY1, dX2, dY2, true, nChamferDivPoints);
	
	for (i = 0; i < nChamferDivPoints; i++)
	{
		pResult[nPoints].x = pChamfer[i].x;
		pResult[nPoints].y = pChamfer[i].y;
		nPoints++;
	}
	gfFreePoints(pChamfer);
	
	dX1 = -0.5 * dWidth;
	dY1 = -0.5 * dHeight + dChamferRadius;
	dX2 = -0.5 * dWidth + dChamferRadius;
	dY2 = -0.5 * dHeight;
	pChamfer = pDivideArc(dX2, dY1, dX1, dY1, dX2, dY2, true, nChamferDivPoints);
	
	for (i = 0; i < nChamferDivPoints; i++)
	{
		pResult[nPoints].x = pChamfer[i].x;
		pResult[nPoints].y = pChamfer[i].y;
		nPoints++;
	}
	gfFreePoints(pChamfer);
	
	dX1 = 0.5 * dWidth - dChamferRadius;
	dY1 = -0.5 * dHeight;
	dX2 = 0.5 * dWidth;
	dY2 = -0.5 * dHeight + dChamferRadius;
	pChamfer = pDivideArc(dX1, dY2, dX1, dY1, dX2, dY2, true, nChamferDivPoints);
	
	for (i = 0; i < nChamferDivPoints; i++)
	{
		pResult[nPoints].x = pChamfer[i].x;
		pResult[nPoints].y = pChamfer[i].y;
		nPoints++;
	}
	gfFreePoints(pChamfer);
	
	// If the rotation angle is a multiple of 180 degrees, no need to rotate or mirror
	if (fmod(dRotDegree, 180) == 0.0)
	{
		for (i = 0; i < nPoints; i++)
		{
			pResult[i].x += dXcenter;
			pResult[i].y += dYcenter;
		}
	}
	// Rotate and mirror
	else
	{
		double dAngle = ONEPI * dRotDegree / 180.0;
		double dCos = cos(dAngle);
		double dSin = sin(dAngle);
		
		for (i = 0; i < nPoints; i++)
		{
			dX1 = pResult[i].x;
			dY1 = pResult[i].y;
			pResult[i].x = dCos * dX1 - dSin * dY1;
			pResult[i].y = dSin * dX1 + dCos * dY1;
			if (nMirror == 'X' || nMirror == 'x')
			{
				pResult[i].x = dXcenter - pResult[i].x;
				pResult[i].y = dYcenter + pResult[i].y;
			}
			else if (nMirror == 'Y' || nMirror == 'y')
			{
				pResult[i].x = dXcenter + pResult[i].x;
				pResult[i].y = dYcenter - pResult[i].y;
			}
			else
			{
				pResult[i].x = dXcenter + pResult[i].x;
				pResult[i].y = dYcenter + pResult[i].y;
			}
		}
	}
	return pResult;
}

gf_point *pRoundedBox2Contour(double dXcenter, double dYcenter,	// box center
							  double dWidth, double dHeight,	// box width and height
							  double dChamferRadius,			// Edge chamfer radius
							  int nChamferDivPoints,			// Number of division points on each chamfer arc
							  int nMirror,						// mirror flag ('X', 'Y' or 'N')
							  double dRotDegree)				// rotation angle (degree in CCW)
{
	// Converts a rounded box defined by center, width, height, chamfer radius, chamfer division points,
	// mirror and rotation to a contour.
	// This function mirrors the rounded box then rotates.
	
	// Because the rounded box is symmetric along X and Y axes before any orientation,
	// mirror first is the same as no mirror. Therefore, mirror->rotate is the same as rotate->no mirror.
	
	return pRoundedBox2Contour(dXcenter, dYcenter, dWidth, dHeight, dChamferRadius, nChamferDivPoints, dRotDegree, 'N');
}

gf_point *pRoundedBox2Contour(double dXcenter, double dYcenter,	// box center
							  double dWidth, double dHeight,	// box width and height
							  double dChamferRadius,			// Edge chamfer radius
							  int nChamferDivPoints,			// Number of division points on each chamfer arc
							  double dRotDegree)				// rotation angle (degree in CCW)
{
	// Converts a rounded box defined by center, width, height, chamfer radius, chamfer division points, rotation to a contour.
	
	return pRoundedBox2Contour(dXcenter, dYcenter, dWidth, dHeight, dChamferRadius, nChamferDivPoints, dRotDegree, 'N');
}

gf_point *pRoundedBox2Contour(double dXcenter, double dYcenter,	// box center
							  double dWidth, double dHeight,	// box width and height
							  double dChamferRadius,			// Edge chamfer radius
							  int nChamferDivPoints)			// Number of division points on each chamfer arc
{
	// Converts a box defined by center, width, height, chamfer radius, chamfer division points to a contour.
	
	return pRoundedBox2Contour(dXcenter, dYcenter, dWidth, dHeight, dChamferRadius, nChamferDivPoints, 0.0, 'N');
}

gf_point *pCircle2Contour(double dXcenter, double dYcenter,	// center
						  double dDia,						// diameter
						  int nSegments,					// number of segments
						  double dStartDegree)				// orientation angle of the 1st point from east 
{
	// Converts a circle to a contour with 'nSegments'.
	// The first point starts from 'dStartAngle' degrees from east
	
    double dRadius = dDia / 2;
    double dRadian = 2 * ONEPI / (double) nSegments;
	double dShift = ONEPI * dStartDegree / 180.0;
	
    gf_point *pResult = (gf_point *) malloc(nSegments * sizeof(gf_point));
	
	int n;
    for (n = 0; n < nSegments; n++)
	{
        double dAngle = dShift + dRadian * (double) n;
        pResult[n].x = dRadius * cos(dAngle) + dXcenter;
        pResult[n].y = dRadius * sin(dAngle) + dYcenter;
    }
    return pResult;
}

gf_point *pCircle2Contour(double dXcenter, double dYcenter,	// center
						  double dDia,						// diameter
						  int nSegments,					// number of segments
						  bool bHalfAngleShift)				// if true, shift 1/2 of one segment angle
{
	// Converts a circle to a contour with 'nSegments'
	// If the bHalfAgnleShift is true, start the first point 1/2 of one segment angle shifted from east.
	// Otherwise, start the first point from east.

	double dStartDegree = 0;
	if (bHalfAngleShift) dStartDegree = 0.5 * 360.0 / (double) nSegments;

    return pCircle2Contour(dXcenter, dYcenter, dDia, nSegments, dStartDegree);
}

gf_point *pCircle2Contour(double dXcenter, double dYcenter,	// center
						  double dDia,						// diameter
						  int nSegments)					// number of segments
{
	// Converts a circle to a contour with 'nSegments'
    return pCircle2Contour(dXcenter, dYcenter, dDia, nSegments, 0.0);
}

gf_point *pEllipse2Contour(double dXcenter, double dYcenter,	// center
						   double a, double b,					// major/2 and minor/2 axis
						   double rot,							// rotation
						   int nSegments)						// number of segments
{
	// Converts an ellipse to a contour with 'nSegments'.

    double dRadian = 2 * ONEPI / (double) nSegments;
	double angle = 0;

	double theta = ONEPI * rot / 180;
	double cosTheta = cos(theta);
	double sinTheta = sin(theta);

    gf_point *pResult = (gf_point *) malloc(nSegments * sizeof(gf_point));

	int n;
    for (n = 0; n < nSegments; n++)
	{
		double cosAlpha = cos(dRadian*n);
		double sinAlpha = sin(dRadian*n);
        pResult[n].x = dXcenter + a * cosAlpha * cosTheta - b * sinAlpha * sinTheta;
        pResult[n].y = dYcenter + a * cosAlpha * sinTheta + b * sinAlpha * cosTheta;
    }
    return pResult;
}

gf_point *pLine2Contour(gf_point *pLine,		// Starting and ending points
						double dWidth,			// Line width
						int nRoundCap,			// Round cap flag (0, 1, 2. 3)
						int nRoundCapPoints,	// Number of round cap discretization points. Unused if straight end cap.
						int *nRetDivPoints)		// outputs number of points in the result point list
{
	// Converts a line segment or a vector to a contour
	//
	// Round cap options (nRoundCap) are:
	//		1: Round cap at the start point only. It will be a straight end at the end point.
	//		2: Round cap at the end point only. It will be a straight end at the start point.
	//		3: Round cap at both ends of the line.
	//		0: No round cap at either end.
	//
	//		-------------       ------------      -----------        -----------
	//		|           |      (           |      |          )      (           )
	//		-------------       ------------      -----------        -----------
	//		nRoundCap = 0      nRoundCap = 1      RoundCap = 2      nRoundCap = 3

	double dEpsilson = dWidth / 200;

	int nPoints;

	switch (nRoundCap)
	{
	case 0:		// No round cap at either end
		nPoints = 4;
		break;
	case 1:		// Round cap at the start point
	case 2:		// ROund cap at the end point
		nPoints = nRoundCapPoints + 2;
		break;
	case 3:		// Round cap at both ends
		nPoints = 2 * nRoundCapPoints;
		break;
	}

	gf_point *pResult = NULL;

	gf_point rLine[2], tmpLine[2];

	// Get a line segment on the right hand side
	// and check if the line segment is too short
	if (!bGetParallelLine(pLine, rLine, 0.5 * dWidth, 0.0, dEpsilson))
	{
		if (nRoundCap == 3)		// Treat it as a circle if round end cap at both ends
		{
			nPoints = nPoints - 2;
			pResult = pCircle2Contour(0.5 * (pLine[0].x + pLine[1].x),
									  0.5 * (pLine[0].y + pLine[1].y),
									  dWidth, nPoints);
			*nRetDivPoints = nPoints;
			return pResult;
		}
		else					// Otherwise, do nothing (skip the short line).
		{
			*nRetDivPoints = 0;
			return pResult;
		}
	}

	pResult = (gf_point *) malloc(nPoints * sizeof(gf_point));

	int p = 0;

	// Keep the two points of the line segment on the right hand side
	pResult[p].x = rLine[0].x;
	pResult[p].y = rLine[0].y;
	p++;
	pResult[p].x = rLine[1].x;
	pResult[p].y = rLine[1].y;
	p++;

	// To get a line segment on the left hand side, switch the start and end points
	tmpLine[0].x = pLine[1].x;
	tmpLine[0].y = pLine[1].y;
	tmpLine[1].x = pLine[0].x;
	tmpLine[1].y = pLine[0].y;

	// Get a line segment on the left
	bGetParallelLine(tmpLine, rLine, 0.5 * dWidth, 0.0, dEpsilson);

	int n;
	gf_point *pArcPoint;

	if (nRoundCap == 2 || nRoundCap == 3)	// Round cap at the end point
	{
		// divide the first round cap arc at the end point
		pArcPoint = pDivideArc(pLine[1], pResult[1], rLine[0], true, nRoundCapPoints);
		for (n = 1; n < nRoundCapPoints-1; n++)
		{
			pResult[p].x = pArcPoint[n].x;
			pResult[p].y = pArcPoint[n].y;
			p++;
		}
		FREE(pArcPoint);
	}

	// Keep the two end points of the left hand side line segment
	pResult[p].x = rLine[0].x;
	pResult[p].y = rLine[0].y;
	p++;
	pResult[p].x = rLine[1].x;
	pResult[p].y = rLine[1].y;
	p++;
	if (nRoundCap == 1 || nRoundCap == 3)
	{
		// divide the second round cap arc at the start point
		pArcPoint = pDivideArc(pLine[0], rLine[1], pResult[0], true, nRoundCapPoints);
		for (n = 1; n < nRoundCapPoints-1; n++)
		{
			pResult[p].x = pArcPoint[n].x;
			pResult[p].y = pArcPoint[n].y;
			p++;
		}
		FREE(pArcPoint);
	}

	*nRetDivPoints = p;

	return pResult;
}

gf_point *pLine2Contour(double dX1, double dY1, double dX2, double dY2,	// Starting and ending points
						double dWidth,			// Line width
						int nRoundCap,			// Round cap flag (0, 1, 2, 3)
						int nRoundCapPoints,	// Number of round cap discretization points. Unused if straight end cap.
						int *nRetDivPoints)		// outputs number of points in the result point list
{
	// Converts a line segment or a vector to a contour.

	gf_point pLine[2];
	pLine[0].x = dX1;
	pLine[0].y = dY1;
	pLine[1].x = dX2;
	pLine[1].y = dY2;

	return pLine2Contour(pLine, dWidth, nRoundCap, nRoundCapPoints, nRetDivPoints);
}

gf_point *pTapered2Contour(gf_point *pLine,		// Starting and ending points
						   double dStartWidth,	// Start width
						   double dEndWidth,	// outputs number of points in the result point list
						   int nRoundCap,		// Round cap flag (0, 1, 2. 3)
						   int nRoundCapPoints,	// Number of round cap discretization points. Unused if straight end cap.
						   int *nRetDivPoints)	// outputs number of points in the result point list
{
	// Converts a tapered line segment to a contour
	//
	// Round cap options (nRoundCap) are:
	//		1: Round cap at the start point only. It will be a straight end at the end point.
	//		2: Round cap at the end point only. It will be a straight end at the start point.
	//		3: Round cap at both ends of the line.
	//		0: No round cap at either end.
	//
	//		-------------       ------------      -----------        -----------
	//		|           |      (           |      |          )      (           )
	//		-------------       ------------      -----------        -----------
	//		nRoundCap = 0      nRoundCap = 1      RoundCap = 2      nRoundCap = 3

	double dX = pLine[1].x - pLine[0].x;
	double dY = pLine[1].y - pLine[1].y;
	double dS = sqrt(dX * dX + dY * dY);

	*nRetDivPoints = 0;
	if (dS < 1.0e-6) return NULL;

	int k;
	int nPoints;
	
	switch (nRoundCap)
	{
	case 0:		// No round cap at either end
		nPoints = 4;
		break;
	case 1:		// Round cap at the start point
	case 2:		// ROund cap at the end point
		nPoints = nRoundCapPoints + 2;
		break;
	case 3:		// Round cap at both ends
		nPoints = 2 * nRoundCapPoints;
		break;
	}

	gf_point *pResult = (gf_point *) malloc(nPoints * sizeof(gf_point));
	gf_point *pPoints = NULL;

	int i = 0;
	double sxStart = dStartWidth * dY / dS;
	double syStart = dStartWidth * dX / dS;
	double sxEnd = dEndWidth * dY / dS;
	double syEnd = dEndWidth * dX / dS;

	double r = 0.5 * (dEndWidth - dStartWidth) / dS;
	double d;
	double x0, y0, x1, y1, x2, y2;

	x1 = pLine[0].x - sxStart;
	y1 = pLine[0].y + syStart;;
	x2 = pLine[0].x + sxStart;
	y2 = pLine[0].y - syStart;

	pResult[i].x = x1;
	pResult[i].y = y1;
	i++;
	if (nRoundCap == 1 || nRoundCap == 3)
	{
		d = 0.5 * dStartWidth * r;
		x0 = d * dX / dS + pLine[0].x;
		y0 = d * dY / dS + pLine[0].y;
		pPoints = pDivideArc(x0, y0, x1, y1, x2, y2, true, nRoundCapPoints);
		for (k = 1; k < nRoundCapPoints - 1; k++)
		{
			pResult[i].x = pPoints[k].x;
			pResult[i].y = pPoints[k].y;
			i++;
		}
		gfFreePoints(pPoints);
	}
	pResult[i].x = x2;
	pResult[i].y = y2;
	i++;
	
	x1 = pLine[1].x + sxEnd;
	y1 = pLine[1].y - syEnd;;
	x2 = pLine[1].x - sxEnd;
	y2 = pLine[1].y + syEnd;
	
	pResult[i].x = x1;
	pResult[i].y = y1;
	i++;
	if (nRoundCap == 2 || nRoundCap == 3)
	{
		d = 0.5 * dEndWidth * r;
		x0 = d * dX / dS + pLine[1].x;
		y0 = d * dY / dS + pLine[1].y;
		pPoints = pDivideArc(x0, y0, x1, y1, x2, y2, true, nRoundCapPoints);
		for (k = 1; k < nRoundCapPoints - 1; k++)
		{
			pResult[i].x = pPoints[k].x;
			pResult[i].y = pPoints[k].y;
			i++;
		}
		gfFreePoints(pPoints);
	}
	pResult[i].x = x2;
	pResult[i].y = y2;
	i++;

	*nRetDivPoints = i;

	return pResult;
}

gf_point *pTapered2Contour(double dX1, double dY1,
						   double dX2, double dY2,	// Starting and ending points
						   double dStartWidth,		// Start width
						   double dEndWidth,		// outputs number of points in the result point list
						   int nRoundCap,			// Round cap flag (0, 1, 2. 3)
						   int nRoundCapPoints,		// Number of round cap discretization points. Unused if straight end cap.
						   int *nRetDivPoints)		// outputs number of points in the result point list
{
	gf_point pLine[2];
	pLine[0].x = dX1;
	pLine[0].y = dY1;
	pLine[1].x = dX2;
	pLine[1].y = dY2;
	
	return pTapered2Contour(pLine, dStartWidth, dEndWidth, nRoundCap, nRoundCapPoints, nRetDivPoints);
}

gf_point *pPath2Contour(gf_point *pLine, int nPointsIn,	// point list and number of points in the point list
						double dWidth,		// width
						int nJoinType,		// line join method (0: strait line, 1: arc, 2: crossing point)
						bool bRoundCap,		// round cap flag (true or false)
						int nMethod,		//  1: generates contour while increasing the point index (for closed path)
											// -1: generates contour while decreasing the point index (for closed path)
											//  0: generates contour in both way (used for open path)
											//  2: generates contour while increasing the point index (used for open path)
											// -2: generates contour while decreasing the point index (used for open path)
						int *nRetDivPoints)	// outputs number of points in the result point list
{
	// Converts a path (polyline) to a contour.
	//
	// To use this function to expand or reduce closed path or contour,
	// first determine whether the closed path is CW or CCW by calling the 'bClockwise' function.
	//
	// To expand the closed path,
	//		use nMethod = 1 for CCW path or nMethod = -1 for CW path
	//		or make the point list ccw before calling this function then use nMethod = 1
	//
	// To shrink the closed path (abridgment),
	//		use nMethod = -1 for CCW path or nMethod = 1 for CW path
	//		or make the point list ccw before calling this function then use nMethod = -1

	double	dDIV_ANGLE = LINE_JOINT_ARC_DIV_ANGLE;
	int		nDIV_MAX = LINE_JOINT_ARC_DIV_MAX_POINT;

	gf_point firstSeg[2], prevSeg[2], newSeg[2], rectArea[4];
	int i, j, k, v;
	int p0;
	int nArcPoints;
	gf_point *pArc;
	double x, y;
	bool bCrossed;

// 	double dEpsilon = dWidth / 20;
	double dEpsilon = 1.0e-4 * dWidth;

	int p = 0;	// Actual point counter
	int q;

	int nPoints;
	if (pLine[0].x == pLine[nPointsIn-1].x && pLine[0].y == pLine[nPointsIn-1].y)
		nPoints = nPointsIn - 1;
	else
		nPoints = nPointsIn;

	// Allocate the result point array with enough number of points

	int nPointsMax = (2 * nPoints) + ((nDIV_MAX-1) * nPoints);
	if (nMethod == 1 || nMethod == -1) nPointsMax = nPointsMax - nPoints;

	gf_point *pTemp = (gf_point *) malloc(nPointsMax * sizeof(gf_point));

	// Go directly to the reverse order segment generation
	if (nMethod < 0) goto LAB_REVERSE;

	// Loop through all points in index increasing order
	if (nMethod == 1)
	{
		i = nPoints - 1;
		j = 0;
	}
	else
	{
		i = 0;
		j = i + 1;
	}
	for (; j < nPoints; j++)
	{
		// Get a line segment parallel to the reference line
		if (!bGetParallelLine(pLine[i].x, pLine[i].y, pLine[j].x, pLine[j].y,
							  &(newSeg[0].x), &(newSeg[0].y), &(newSeg[1].x), &(newSeg[1].y),
							  0.5 * dWidth, 0.0, dEpsilon)) continue;
		q = p;
		if (p == 0)
		{
			// keep the first segment
			firstSeg[0].x = newSeg[0].x;
			firstSeg[0].y = newSeg[0].y;
			firstSeg[1].x = newSeg[1].x;
			firstSeg[1].y = newSeg[1].y;

			// keep the first point
			pTemp[p].x = newSeg[0].x;
			pTemp[p].y = newSeg[0].y;
			p++;
		}
		else
		{
			rectArea[0].x = pLine[k].x;
			rectArea[0].y = pLine[k].y;
			rectArea[1].x = pLine[j].x;
			rectArea[1].y = pLine[j].y;
			rectArea[2].x = newSeg[1].x;
			rectArea[2].y = newSeg[1].y;
			rectArea[3].x = newSeg[0].x;
			rectArea[3].y = newSeg[0].y;
			if (nTwoLines(pLine[k].x, pLine[k].y, prevSeg[1].x, prevSeg[1].y,
						  pLine[j].x, pLine[j].y, newSeg[1].x, newSeg[1].y,
						  &x, &y) == 3)
			{
				// Do nothing
			}
			else if ((nPoints > 3 || nMethod != 1) &&
					  p > 1 && nTwoLines(&pTemp[p-2], newSeg, &x, &y) == 3)
			{
				pTemp[p-1].x = x;
				pTemp[p-1].y = y;
				q--;
			}
			else if (PointInContour(prevSeg[0].x, prevSeg[0].y, 4, rectArea) &&
					 PointInContour(prevSeg[1].x, prevSeg[1].y, 4, rectArea))
			{
				if (PointProject2Line(prevSeg[0], newSeg[0], newSeg[1], &pTemp[p]))
				{
					p++;
				}
			}
			else
			{
				// joint the previous parallel line and the new parallel line

				bCrossed = bJoinLines(pLine[i], prevSeg, newSeg, nJoinType, pTemp, &p);
				if (!bCrossed && p == q)
				{
					for (v = p - 3; v >= 0; v--)
					{
						if (nTwoLines(&pTemp[v], newSeg, &x, &y) == 3)
						{
							pTemp[v+1].x = x;
							pTemp[v+1].y = y;
							p = v + 2;
							break;
						}
					}
				}
			}
		}

		if (p != q)
		{
			prevSeg[0].x = newSeg[0].x;
			prevSeg[0].y = newSeg[0].y;
			prevSeg[1].x = newSeg[1].x;
			prevSeg[1].y = newSeg[1].y;
			k = j;
		}
		i = j;
	}
	
	if (p > 0)
	{
		// Keep the last point
		pTemp[p].x = prevSeg[1].x;
		pTemp[p].y = prevSeg[1].y;
		p++;
	}

	if (nMethod == 1)
	{
		// If only the increasing order was requested, join the first segment
		// and the last segment to close the contour
		if (p > 3)
		{
			if (bJoinLines(pLine[nPoints-1], prevSeg, firstSeg, nJoinType, pTemp, &p))
			{
				pTemp[0].x = pTemp[p-1].x;
				pTemp[0].y = pTemp[p-1].y;
				p = p - 2;
			}
		}
		if (p > 3)
		{
			prevSeg[0].x = pTemp[p-1].x;
			prevSeg[0].y = pTemp[p-1].y;
			prevSeg[1].x = pTemp[0].x;
			prevSeg[1].y = pTemp[0].y;
			for (j = p - 3; j >= 1; j--)
			{
				if (nTwoLines(prevSeg, &pTemp[j], &x, &y) == 3)
				{
					pTemp[0].x = x;
					pTemp[0].y = y;
					i = 1;
					j++;
					for ( ; j < p; j++)
					{
						pTemp[i].x = pTemp[j].x;
						pTemp[i].y = pTemp[j].y;
						i++;
					}
					p = i;
					break;
				}
			}
		}
		if (p > 3)
		{
			for (j = p - 2; j >= 2; j--)
			{
				if (nTwoLines(&pTemp[0], &pTemp[j], &x, &y) == 3)
				{
					pTemp[0].x = x;
					pTemp[0].y = y;
					p = j + 1;
					break;
				}
			}
		}
	}
	else if (nMethod == 2)
	{
		for (j = nPoints - 1; j >= 0; j--)
		{
			pTemp[p].x = pLine[j].x;
			pTemp[p].y = pLine[j].y;
			p++;
		}
	}

	if (nMethod != 0 || p < 2)	goto LAB_RETURN;	// No need to the reverse path

LAB_REVERSE:
	// build in reverse order
	p0 = p;

	if (nMethod == -1)
	{
		i = 0;
		j = nPoints - 1;
	}
	else
	{
		i = nPoints - 1;
		j = i - 1;
	}
	for (; j >= 0; j--)
	{
		if (!bGetParallelLine(pLine[i].x, pLine[i].y, pLine[j].x, pLine[j].y,
							  &(newSeg[0].x), &(newSeg[0].y), &(newSeg[1].x), &(newSeg[1].y),
							  0.5 * dWidth, 0.0, dEpsilon)) continue;
		q = p;
		if (p == p0)
		{
			if (nMethod == 0 && bRoundCap)
			{
				// put round cap at one end
				pArc = pDivideArc(pLine[i], prevSeg[1], newSeg[0], true, dDIV_ANGLE, 2, nDIV_MAX, &nArcPoints);
				for (j = 1; j < nArcPoints-1; j++)
				{
					pTemp[p].x = pArc[j].x;
					pTemp[p].y = pArc[j].y;
					p++;
				}
				FREE(pArc);
			}
			firstSeg[0].x = newSeg[0].x;
			firstSeg[0].y = newSeg[0].y;
			firstSeg[1].x = newSeg[1].x;
			firstSeg[1].y = newSeg[1].y;
			pTemp[p].x = newSeg[0].x;
			pTemp[p].y = newSeg[0].y;
			p++;
		}
		else
		{
			rectArea[0].x = pLine[k].x;
			rectArea[0].y = pLine[k].y;
			rectArea[1].x = pLine[j].x;
			rectArea[1].y = pLine[j].y;
			rectArea[2].x = newSeg[1].x;
			rectArea[2].y = newSeg[1].y;
			rectArea[3].x = newSeg[0].x;
			rectArea[3].y = newSeg[0].y;
			if (nTwoLines(pLine[k].x, pLine[k].y, prevSeg[1].x, prevSeg[1].y,
						  pLine[j].x, pLine[j].y, newSeg[1].x, newSeg[1].y,
						  &x, &y) == 3)
			{
				// Do nothing
			}
			else if ((nPoints > 3 || nMethod != -1) &&
					  p > 2 && p > p0 + 1 && nTwoLines(&pTemp[p-2], newSeg, &x, &y) == 3)
			{
				pTemp[p-1].x = x;
				pTemp[p-1].y = y;
				q--;
			}
			else if (p > 0 && nTwoLines(prevSeg, newSeg, &x, &y) == 3)
			{
				pTemp[p].x = x;
				pTemp[p].y = y;
				p++;
				if (j == 0)
				{
					if (nTwoLines(newSeg, &pTemp[0], &x, &y) == 3)
					{
						pTemp[0].x = x;
						pTemp[0].y = y;
					}
					else
					{
						pTemp[p].x = newSeg[1].x;
						pTemp[p].y = newSeg[1].y;
						p++;
					}
				}
			}
			else if (PointInContour(prevSeg[0].x, prevSeg[0].y, 4, rectArea) &&
					 PointInContour(prevSeg[1].x, prevSeg[1].y, 4, rectArea))
			{
				if (PointProject2Line(prevSeg[0], newSeg[0], newSeg[1], &pTemp[p]))
				{
					p++;
				}
			}
			else
			{
				bCrossed = bJoinLines(pLine[i], prevSeg, newSeg, nJoinType, pTemp, &p);
				if (!bCrossed && p == q)
				{
					for (v = p - 3; v >= 0; v--)
					{
						if (nTwoLines(&pTemp[v], newSeg, &x, &y) == 3)
						{
							pTemp[v+1].x = x;
							pTemp[v+1].y = y;
							p = v + 2;
							break;
						}
					}
					q--;
				}
			}
		}

		if (p != q)
		{
			prevSeg[0].x = newSeg[0].x;
			prevSeg[0].y = newSeg[0].y;
			prevSeg[1].x = newSeg[1].x;
			prevSeg[1].y = newSeg[1].y;
			k = j;
		}
		i = j;
	}
	
	if (p > 0)
	{
		pTemp[p].x = prevSeg[1].x;
		pTemp[p].y = prevSeg[1].y;
		p++;
		if (nMethod == 0 && bRoundCap)
		{
			// put round cap at another end
			pArc = pDivideArc(pLine[j], prevSeg[1], pTemp[0], true, dDIV_ANGLE, 2, nDIV_MAX, &nArcPoints);
			for (j = 1; j < nArcPoints-1; j++)
			{
				pTemp[p].x = pArc[j].x;
				pTemp[p].y = pArc[j].y;
				p++;
			}
			FREE(pArc);
		}
	}

	if (nMethod == -1)
	{
		if (p > 3)
		{
			if (bJoinLines(pLine[0], prevSeg, firstSeg, nJoinType, pTemp, &p))
			{
				pTemp[0].x = pTemp[p-1].x;
				pTemp[0].y = pTemp[p-1].y;
				p = p - 2;
			}
		}
		if (p > 3)
		{
			prevSeg[0].x = pTemp[p-1].x;
			prevSeg[0].y = pTemp[p-1].y;
			prevSeg[1].x = pTemp[0].x;
			prevSeg[1].y = pTemp[0].y;
			for (j = p - 3; j >= 1; j--)
			{
				if (nTwoLines(prevSeg, &pTemp[j], &x, &y) == 3)
				{
					pTemp[0].x = x;
					pTemp[0].y = y;
					i = 1;
					j++;
					for ( ; j < p; j++)
					{
						pTemp[i].x = pTemp[j].x;
						pTemp[i].y = pTemp[j].y;
						i++;
					}
					p = i;
					break;
				}
			}
		}
		if (p > 3)
		{
			for (j = p - 2; j >= 2; j--)
			{
				if (nTwoLines(&pTemp[0], &pTemp[j], &x, &y) == 3)
				{
					pTemp[0].x = x;
					pTemp[0].y = y;
					p = j + 1;
					break;
				}
			}
		}
	}
	else if (nMethod == -2)
	{
		for (j = 0; j < nPoints; j++)
		{
			pTemp[p].x = pLine[j].x;
			pTemp[p].y = pLine[j].y;
			p++;
		}
	}

LAB_RETURN:

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

gf_polygon *pMoire2Contour(double dXC, double dYC,					// Center point location
						   double dOutDia,							// iner diameter
						   double dRingThick, double dGapThick,		// ring and gap thickness
						   int nRings,								// number of rings
						   double dHairThick, double dHairLength,	// Cross hair thickness and length
						   double dRotDegree,						// rotation in CCW degree
						   int nDivPoints)							// number of circle division points
{
	// Converts a Moire circle that is defined as a Gerber Aperture Macro primitive number 6 to a set of contours
	
	gf_polygon *pRing, *pHair, *pResult;
	
	pRing = gfAllocPolygon(nRings * 2);
	
	// Build rings
	double dDia = dOutDia;
	int c = 0;
	int n;
	for (n = 0; n < nRings; n++)
	{
		pRing->hole[c] = false;
		pRing->contour[c].nPoints = nDivPoints;
		pRing->contour[c].point = pCircle2Contour(dXC, dYC, dDia, nDivPoints);
		c++;
		dDia -= 2 * dRingThick;
		
		pRing->hole[c] = true;
		pRing->contour[c].nPoints = nDivPoints;
		pRing->contour[c].point = pCircle2Contour(dXC, dYC, dDia, nDivPoints);
		c++;
		dDia -= 2 * dGapThick;
	}
	
	// Merge the rotated horizontal hair to the ring
	pHair = gfAllocPolygon(1);
	pHair->hole[0] = false;
	pHair->contour[0].nPoints = 4;
	pHair->contour[0].point = pBox2Contour(dXC, dYC, dHairLength, dHairThick, dRotDegree);
	
	pResult = gfAllocPolygon(0);
	gfPolygonBool(GF_UNION, pRing, pHair, pResult);
	gfFreePolygon(pRing);
	gfFreePolygon(pHair);
	pRing = pResult;
	
	// Merge the rotated vertical hair to the ring
	pHair = gfAllocPolygon(1);
	pHair->hole[0] = false;
	pHair->contour[0].nPoints = 4;
	pHair->contour[0].point = pBox2Contour(dXC, dYC, dHairThick, dHairLength, dRotDegree);
	
	pResult = gfAllocPolygon(0);
	gfPolygonBool(GF_UNION, pRing, pHair, pResult);
	gfFreePolygon(pRing);
	gfFreePolygon(pHair);
	
	pResult = gfContourHierarchy(pResult);
	
	return pResult;
}

gf_polygon *pThermal2Contour(double dXC, double dYC,		// Center point location
							 double dOutDia, double dInDia,	// Outer and inner diameters
							 double dHairWidth,				// cross hair width
							 double dRotDegree,				// rotation in CCW degree
							 int nDivPoints)				// number of circle division points
{
	// Convert a wagon wheel shape to a set of contours. The cross hairs in unrotated shape are
	// plus sign (+) shape, i.e. horizontal and vertial.
	// This is the shape defined as a Gerber Aperture Macro primitive number 7 that is Thermal.
	// In order to properly make a thermal relief shape, the result contours from this function
	// may need to be subtracted from a solid contour that enclose the thermal relief shape.

	gf_polygon *pDonut, *pHair, *pResult;

	// Construct a polygon that forms a donut shape
	pDonut = gfAllocPolygon(2);
	pDonut->hole[0] = false;
	pDonut->contour[0].nPoints = nDivPoints;
	pDonut->contour[0].point = pCircle2Contour(dXC, dYC, dOutDia, nDivPoints);

	pDonut->hole[1] = true;
	pDonut->contour[1].nPoints = nDivPoints;
	pDonut->contour[1].point = pCircle2Contour(dXC, dYC, dInDia, nDivPoints);

	// Subtract a horizontal box that is rotated 'dRotDegree'
	pHair = gfAllocPolygon(1);
	pHair->hole[0] = false;
	pHair->contour[0].nPoints = 4;
	pHair->contour[0].point = pBox2Contour(dXC, dYC, 1.1*dOutDia, dHairWidth, dRotDegree);

	pResult = gfAllocPolygon(0);
	gfPolygonBool(GF_SUB, pDonut, pHair, pResult);
	gfFreePolygon(pDonut);
	gfFreePolygon(pHair);
	pDonut = pResult;

	// Subtract a vertical box that is rotated 'dRotDegree'
	pHair = gfAllocPolygon(1);
	pHair->hole[0] = false;
	pHair->contour[0].nPoints = 4;
	pHair->contour[0].point = pBox2Contour(dXC, dYC, dHairWidth, 1.1*dOutDia, dRotDegree);

	pResult = gfAllocPolygon(0);
	gfPolygonBool(GF_SUB, pDonut, pHair, pResult);
	gfFreePolygon(pDonut);
	gfFreePolygon(pHair);

	pResult = gfContourHierarchy(pResult);

	return pResult;
}

gf_point *pArc2Contour(double dXC, double dYC,		// Arc center
					   double dX1, double dY1,		// Starting arc point
					   double dX2, double dY2,		// Ending arc point
					   bool bCCW,					// Arc direction CCW (true), CW (false)
					   int nDivPoints,				// number of arc division points
					   double dWidth,				// arc width
					   int nRoundCapPoints,			// round cap division points. If < 3, straight end
					   int *nRetDivPoints)			// outputs number of points in the result point list
{
	// Convert an arc to a contour.
	int nPoints = 0;
	gf_point *pResult = NULL;

	// Calculate radius
	double dRadius = sqrt((dXC - dX1) * (dXC - dX1) + (dYC - dY1) * (dYC - dY1));

	// Do not generate a contour if the radius is zero
	*nRetDivPoints = nPoints;
	if (dRadius <= 0) return pResult;

	gf_point *pArc;
	int i;
	int p = 0;

	// Determine start and end point of the inner arc then divide the inner arc
	double dRatio = (dRadius + 0.5 * dWidth) / dRadius;
	double dXout1 = dXC + (dX1 - dXC) * dRatio;
	double dYout1 = dYC + (dY1 - dYC) * dRatio;
	double dXout2 = dXC + (dX2 - dXC) * dRatio;
	double dYout2 = dYC + (dY2 - dYC) * dRatio;
	pArc = pDivideArc(dXC, dYC, dXout1, dYout1, dXout2, dYout2, bCCW, nDivPoints);

	// If the half of width is bigger than the radius, treat it as a pie.
	if (((0.5 * dWidth) - dRadius) > (-0.001 * dRadius))	// consider round-off error
	{
		nPoints = nDivPoints + 1;
		pResult = (gf_point *) malloc(nPoints * sizeof(gf_point));
		pResult[p].x = dXC;
		pResult[p].y = dYC;
		p++;
		for (i = 0; i < nDivPoints; i++)
		{
			pResult[p].x = pArc[i].x;
			pResult[p].y = pArc[i].y;
			p++;
		}
		FREE(pArc);
		*nRetDivPoints = nPoints;
		return pResult;
	}

	// Allocate the result point array
	if (nRoundCapPoints < 3)
		nPoints = 2 * nDivPoints;
	else
		nPoints = 2 * nDivPoints + 2 * (nRoundCapPoints - 2);

	pResult = (gf_point *) malloc(nPoints * sizeof(gf_point));

	// Copy the outer arc division points to the result point array
	for (i = 0; i < nDivPoints; i++)
	{
		pResult[p].x = pArc[i].x;
		pResult[p].y = pArc[i].y;
		p++;
	}
	FREE(pArc);

	// Determine start and end point of the inner arc
	dRatio = (dRadius - 0.5 * dWidth) / dRadius;
	double dXin1 = dXC + (dX2 - dXC) * dRatio;
	double dYin1 = dYC + (dY2 - dYC) * dRatio;
	double dXin2 = dXC + (dX1 - dXC) * dRatio;
	double dYin2 = dYC + (dY1 - dYC) * dRatio;

	// For the round cap
	if (nRoundCapPoints >= 3)
	{
		pArc = pDivideArc(dX2, dY2, dXout2, dYout2, dXin1, dYin1, bCCW, nRoundCapPoints);
		for (i = 1; i < nRoundCapPoints-1; i++)
		{
			pResult[p].x = pArc[i].x;
			pResult[p].y = pArc[i].y;
			p++;
		}
		FREE(pArc);
	}

	// divide inner arc
	pArc = pDivideArc(dXC, dYC, dXin1, dYin1, dXin2, dYin2, !bCCW, nDivPoints);
	for (i = 0; i < nDivPoints; i++)
	{
		pResult[p].x = pArc[i].x;
		pResult[p].y = pArc[i].y;
		p++;
	}
	FREE(pArc);

	// For the round cap
	if (nRoundCapPoints >= 3)
	{
		pArc = pDivideArc(dX1, dY1, dXin2, dYin2, dXout1, dYout1, bCCW, nRoundCapPoints);
		for (i = 1; i < nRoundCapPoints-1; i++)
		{
			pResult[p].x = pArc[i].x;
			pResult[p].y = pArc[i].y;
			p++;
		}
		FREE(pArc);
	}

	*nRetDivPoints = nPoints;
	return pResult;
}

#define DOES_NEEDED 1

gf_polygon *pExpandContour(gf_polygon *pPoly,
						   int iContour,
						   int nContours,			// 'nContours' starting from the contour number 'iContour' of 'pPoly'
						   double dExpand,			// expand amount. If negative, shrink the contour (abridgment)
						   bool bTopContourOnly,	// Expand or shrink the top (1st) contour only
						   int nJoinType)			// line join method (0: strait line, 1: arc, 2: crossing point of extended lines)
{
	// Expand (dExpand > 0) or shrink (dExpand < 0) a contour

	if (dExpand == 0) return NULL;				// return the original polygon if the expansion is zero.
	
	double dIgnoreSize = 2 * ABS(dExpand);
	int i, c, d, v;
	int ignored;
	int nKept = 0;

	gf_polygon *pSource;
	gf_polygon *pTarget;
	gf_polygon *pResult;
	if (bTopContourOnly && dExpand > 0)
	{
		if (pPoly->hole[iContour] &&
			((pPoly->box[iContour].ur.x - pPoly->box[iContour].ll.x) <= dIgnoreSize ||
			 (pPoly->box[iContour].ur.y - pPoly->box[iContour].ll.y) <= dIgnoreSize))
		{
			return NULL;
		}
		pResult = gfCopyContours2Polygon(pPoly, iContour, 1);
		nKept = 1;
	}
	else if (dExpand > 0)
	{
		pResult = gfAllocPolygon(nContours);
		i = 0;
		ignored = -1;
		
		for (c = 0, d = iContour; c < nContours; c++, d++)
		{
			if (ABS(pPoly->contour[d].nPoints) < 3) continue;

			if (pPoly->hole[d] &&
				((pPoly->box[d].ur.x - pPoly->box[d].ll.x) <= dIgnoreSize ||
				 (pPoly->box[d].ur.y - pPoly->box[d].ll.y) <= dIgnoreSize))
			{
				pPoly->contour[d].nPoints = -1 * ABS(pPoly->contour[d].nPoints);
				ignored = d;
				continue;
			}

			if (ignored >= 0)
			{
				if (PointInContour(pPoly->contour[d].point[0].x, pPoly->contour[d].point[0].y,
								   pPoly->contour[ignored].nPoints, pPoly->contour[ignored].point,
								   pPoly->box[ignored].ll.x, pPoly->box[ignored].ll.y,
								   pPoly->box[ignored].ur.x, pPoly->box[ignored].ur.y))
				{
					pPoly->contour[d].nPoints = -1 * ABS(pPoly->contour[d].nPoints);		// Flag it temporarily to ignore
					continue;
				}
			}

			pResult->hole[i] = pPoly->hole[d];
			
			gfAllocContour(pResult, i, abs(pPoly->contour[d].nPoints));
			
			for (v = 0; v < pPoly->contour[d].nPoints; v++)
			{
				pResult->contour[i].point[v].x = pPoly->contour[d].point[v].x;
				pResult->contour[i].point[v].y = pPoly->contour[d].point[v].y;
			}
			
			pResult->box[i].ll.x = pPoly->box[d].ll.x;
			pResult->box[i].ll.y = pPoly->box[d].ll.y;
			pResult->box[i].ur.x = pPoly->box[d].ur.x;
			pResult->box[i].ur.y = pPoly->box[d].ur.y;
			i++;
			ignored = -1;
		}
		nKept = i;
	}
	else
	{
		pResult = gfAllocPolygon(nContours);
		i = 0;
		ignored = -1;

		for (c = 0, d = iContour; c < nContours; c++, d++)
		{
			if (ABS(pPoly->contour[d].nPoints) < 3) continue;
			
			if (pPoly->hole[d] == false &&
				((pPoly->box[d].ur.x - pPoly->box[d].ll.x) <= dIgnoreSize ||
				(pPoly->box[d].ur.y - pPoly->box[d].ll.y) <= dIgnoreSize))
			{
				pPoly->contour[d].nPoints = -1 * ABS(pPoly->contour[d].nPoints);
				ignored = d;
				continue;
			}
			
			if (ignored >= 0)
			{
				if (pPoly->hole[d])
				{
					pPoly->contour[d].nPoints = -1 * ABS(pPoly->contour[d].nPoints);		// Flag it temporarily to ignore
					continue;
				}
			}
			
			pResult->hole[i] = pPoly->hole[d];
			
			gfAllocContour(pResult, i, abs(pPoly->contour[d].nPoints));
			
			for (v = 0; v < pPoly->contour[d].nPoints; v++)
			{
				pResult->contour[i].point[v].x = pPoly->contour[d].point[v].x;
				pResult->contour[i].point[v].y = pPoly->contour[d].point[v].y;
			}
			
			pResult->box[i].ll.x = pPoly->box[d].ll.x;
			pResult->box[i].ll.y = pPoly->box[d].ll.y;
			pResult->box[i].ur.x = pPoly->box[d].ur.x;
			pResult->box[i].ur.y = pPoly->box[d].ur.y;
			i++;
			ignored = -1;
		}
		nKept = i;
	}

	if (nKept == 0)
	{
		for (c = 0, d = iContour; c < nContours; c++, d++)
		{
			pPoly->contour[d].nPoints = ABS(pPoly->contour[d].nPoints);
		}
		gfFreePolygon(pResult);
		return NULL;
	}

	pResult->nContours = nKept;

	int iLastContour;
	if (bTopContourOnly)
		iLastContour = iContour;					// Expand/Shrink the 1st contour only
	else
		iLastContour = iContour + nContours - 1;	// Expand/Shrink all contours

	double eps = 0.001 * ABS(dExpand);
// 	eps = 0;

	// Loop through all contours
	for (c = iContour; c <= iLastContour; c++)
	{
		int nPoints = pPoly->contour[c].nPoints;
		
		if (nPoints < 3) continue;		// Less than 3 point is not a contour. This also skips the ignored contours.
		
		gf_point twoLines[3];
		double dX, dY;
		
		i = nPoints - 2;
		twoLines[0].x = pPoly->contour[c].point[i].x;
		twoLines[0].y = pPoly->contour[c].point[i].y;
		dX = twoLines[0].x - pPoly->contour[c].point[i-1].x;
		dY = twoLines[0].y - pPoly->contour[c].point[i-1].y;
		if (ABS(dX) < eps && ABS(dY) < eps)
		{
			twoLines[0].x = pPoly->contour[c].point[i-1].x;
			twoLines[0].y = pPoly->contour[c].point[i-1].y;
		}
		while (true)
		{
			i++;
			if (i >= nPoints) i = 0;
			twoLines[1].x = pPoly->contour[c].point[i].x;
			twoLines[1].y = pPoly->contour[c].point[i].y;
			dX = twoLines[1].x - twoLines[0].x;
			dY = twoLines[1].y - twoLines[0].y;
			if (ABS(dX) >= eps || ABS(dY) >= eps) break;
			if (i == nPoints - 2) break;
		}
		i++;
		if (i >= nPoints) i = 0;

		double dCrossP0, dCrossP1;
		dCrossP0 = 0;

		for (; i < nPoints; i++)
		{
			twoLines[2].x = pPoly->contour[c].point[i].x;
			twoLines[2].y = pPoly->contour[c].point[i].y;

			dX = twoLines[2].x - twoLines[1].x;
			dY = twoLines[2].y - twoLines[1].y;
			if (ABS(dX) < eps && ABS(dY) < eps) continue;
			
			dCrossP1 = dCrossProduct(twoLines, &twoLines[1]);
			if (dCrossP0 * dCrossP1 < 0) break;
			if (dCrossP1 != 0) dCrossP0 = dCrossP1;

			twoLines[0].x = twoLines[1].x;
			twoLines[0].y = twoLines[1].y;
			twoLines[1].x = twoLines[2].x;
			twoLines[1].y = twoLines[2].y;
		}

		if (i == nPoints)	// Convex contour
		{
			int nMethod = 1;
			if (bClockwise(pPoly->contour[c].point, pPoly->contour[c].nPoints)) nMethod = -1;

			if (
#if DOES_NEEDED
				dExpand >= 0 || 
#endif
				(pPoly->box[c].ur.x - pPoly->box[c].ll.x) <= ABS(dExpand) ||
				(pPoly->box[c].ur.y - pPoly->box[c].ll.y) <= ABS(dExpand))
			{
				pTarget = gfAllocPolygon(1);
			}
			else
			{
				pTarget = gfAllocPolygon(2);
			}

			pTarget->hole[0] = false;
			pTarget->contour[0].point
				= pPath2Contour(pPoly->contour[c].point, pPoly->contour[c].nPoints,
								2*ABS(dExpand), nJoinType, false, nMethod,
								&(pTarget->contour[0].nPoints));

			if (pTarget->nContours > 1)
			{
				nMethod *= -1;
				pTarget->hole[1] = true;
				pTarget->contour[1].point
					= pPath2Contour(pPoly->contour[c].point, pPoly->contour[c].nPoints,
									2*ABS(dExpand), nJoinType, false, nMethod,
									&(pTarget->contour[1].nPoints));
			}
		}

		else				// Non-Convex contour
		{
			pTarget = gfAllocPolygon(0);
		
			// Merge bands along the contour perimeter.
			// Do two line segments at a time to ensure proper joining of neighboring lines.

			i = nPoints - 2;
			twoLines[0].x = pPoly->contour[c].point[i].x;
			twoLines[0].y = pPoly->contour[c].point[i].y;
			dX = twoLines[0].x - pPoly->contour[c].point[i-1].x;
			dY = twoLines[0].y - pPoly->contour[c].point[i-1].y;
			if (ABS(dX) < eps && ABS(dY) < eps)
			{
				twoLines[0].x = pPoly->contour[c].point[i-1].x;
				twoLines[0].y = pPoly->contour[c].point[i-1].y;
			}
			while (true)
			{
				i++;
				if (i >= nPoints) i = 0;
				twoLines[1].x = pPoly->contour[c].point[i].x;
				twoLines[1].y = pPoly->contour[c].point[i].y;
				dX = twoLines[1].x - twoLines[0].x;
				dY = twoLines[1].y - twoLines[0].y;
				if (ABS(dX) >= eps || ABS(dY) >= eps) break;
				if (i == nPoints - 2) break;
			}
			i++;
			if (i >= nPoints) i = 0;

			for (; i < nPoints; i++)
			{
				twoLines[2].x = pPoly->contour[c].point[i].x;
				twoLines[2].y = pPoly->contour[c].point[i].y;

				dX = twoLines[2].x - twoLines[1].x;
				dY = twoLines[2].y - twoLines[1].y;
				if (ABS(dX) < eps && ABS(dY) < eps) continue;

				pSource = gfAllocPolygon(1);

				pSource->hole[0] = false;
				pSource->contour[0].point
					= pPath2Contour(twoLines, 3, 2*ABS(dExpand), nJoinType, false, 0, &(pSource->contour[0].nPoints));

				gf_polygon *pUnion = gfAllocPolygon(0);

				gfPolygonBool(GF_UNION, pTarget, pSource, pUnion);

				gfFreePolygon(pTarget);
				gfFreePolygon(pSource);
				pTarget = pUnion;

				twoLines[0].x = twoLines[1].x;
				twoLines[0].y = twoLines[1].y;
				twoLines[1].x = twoLines[2].x;
				twoLines[1].y = twoLines[2].y;
			}
		}

		// Union or subtract the bands on the original contours.
		gf_polygon *pMerged = gfAllocPolygon(0);

		if (dExpand > 0)
			gfPolygonBool(GF_UNION, pResult, pTarget, pMerged);
		else
			gfPolygonBool(GF_SUB, pResult, pTarget, pMerged);

		gfFreePolygon(pResult);
		gfFreePolygon(pTarget);
		pResult = pMerged;
	}

	for (c = 0, d = iContour; c < nContours; c++, d++)
	{
		pPoly->contour[d].nPoints = ABS(pPoly->contour[d].nPoints);
	}

	if (bTopContourOnly && dExpand > 0 && nContours > 1)
	{
		gfAppendContours2Polygon(pPoly, iContour+1, nContours-1, pResult);
		pTarget = gfPolygonMergeAll(pResult);
		gfFreePolygon(pResult);
		pResult = pTarget;
	}
	pResult = gfContourHierarchy(pResult);

	return pResult;
}

gf_polygon *pExpandContour(gf_polygon *pPoly,
						   int iContour,
						   int nContours,			// 'nContours' starting from the contour number 'iContour' of 'pPoly'
						   double dExpand,			// expand amount. If negative, shrink the contour (abridgment)
						   bool bTopContourOnly,	// Expand or shrink the top (1st) contour only
						   bool bArcJoin)			// join method between two neighboring nets(true: arc, false: straight line)
{
	int nJoinType = 0;
	if (bArcJoin) nJoinType = 1;
	return pExpandContour(pPoly, iContour, nContours, dExpand, bTopContourOnly,nJoinType);
}

gf_point *pContourBoxClip(gf_point *pPoints,		// Coordinates of input contour
						  int nPoints,				// Number of points in the input contour
						  gf_box box,				// A box that will clip the input contour
						  int *nReturnPoints)		// Number of points in a result contour
{
	double LFT = box.ll.x;
	double RHT = box.ur.x;
	double BTM = box.ll.y;
	double TOP = box.ur.y;

	gf_point *pTemp = (gf_point *) malloc (2 * nPoints * sizeof(gf_point)) ;
	gf_point *pClip = (gf_point *) malloc (2 * nPoints * sizeof(gf_point)) ;

	int nPointsNew = nPoints;

	int k;
	int u, v;
	bool inPrev, inCur;
	double dX, dY;

	// Check against the right edge
	k = 0;
	u = nPointsNew - 1;
	inPrev = pPoints[u].x > RHT ? false : true;
	for (v = 0; v < nPointsNew; v++)
	{
		inCur = pPoints[v].x > RHT ? false : true;

		if (inPrev != inCur)	// (inPrev && !inCur) || (!inPrev && inCur)
		{
			dX = pPoints[v].x - pPoints[u].x;
			dY = pPoints[v].y - pPoints[u].y;
			pTemp[k].x = RHT;
			pTemp[k].y = pPoints[u].y + (dY / dX) * (RHT - pPoints[u].x);
			k++;
		}
		if (inCur)
		{
			pTemp[k].x = pPoints[v].x;
			pTemp[k].y = pPoints[v].y;
			k++;
		}

		inPrev = inCur;
		u = v;
	}
	nPointsNew = k;
	if (nPointsNew < 3)
	{
		FREE(pTemp);
		FREE(pClip);
		*nReturnPoints = 0;
		return NULL;
	}

	// Check against the top edge
	k = 0;
	u = nPointsNew - 1;
	inPrev = pTemp[u].y > TOP ? false : true;
	for (v = 0; v < nPointsNew; v++)
	{
		inCur = pTemp[v].y > TOP ? false : true;

		if (inPrev != inCur)	// (inPrev && !inCur) || (!inPrev && inCur)
		{
			dX = pTemp[v].x - pTemp[u].x;
			dY = pTemp[v].y - pTemp[u].y;
			pClip[k].y = TOP;
			pClip[k].x = pTemp[u].x + (dX / dY) * (TOP - pTemp[u].y);
			k++;
		}
		if (inCur)
		{
			pClip[k].x = pTemp[v].x;
			pClip[k].y = pTemp[v].y;
			k++;
		}

		inPrev = inCur;
		u = v;
	}
	nPointsNew = k;
	if (nPointsNew < 3)
	{
		FREE(pTemp);
		FREE(pClip);
		*nReturnPoints = 0;
		return NULL;
	}

	// Check against the left edge
	k = 0;
	u = nPointsNew - 1;
	inPrev = pClip[u].x < LFT ? false : true;
	for (v = 0; v < nPointsNew; v++)
	{
		inCur = pClip[v].x < LFT ? false : true;

		if (inPrev != inCur)	// (inPrev && !inCur) || (!inPrev && inCur)
		{
			dX = pClip[v].x - pClip[u].x;
			dY = pClip[v].y - pClip[u].y;
			pTemp[k].x = LFT;
			pTemp[k].y = pClip[u].y + (dY / dX) * (LFT - pClip[u].x);
			k++;
		}
		if (inCur)
		{
			pTemp[k].x = pClip[v].x;
			pTemp[k].y = pClip[v].y;
			k++;
		}

		inPrev = inCur;
		u = v;
	}
	nPointsNew = k;
	if (nPointsNew < 3)
	{
		FREE(pTemp);
		FREE(pClip);
		*nReturnPoints = 0;
		return NULL;
	}

	// Check against the bottom edge
	k = 0;
	u = nPointsNew - 1;
	inPrev = pTemp[u].y < BTM ? false : true;
	for (v = 0; v < nPointsNew; v++)
	{
		inCur = pTemp[v].y < BTM ? false : true;

		if (inPrev != inCur)	// (inPrev && !inCur) || (!inPrev && inCur)
		{
			dX = pTemp[v].x - pTemp[u].x;
			dY = pTemp[v].y - pTemp[u].y;
			pClip[k].y = BTM;
			pClip[k].x = pTemp[u].x + (dX / dY) * (BTM - pTemp[u].y);
			k++;
		}
		if (inCur)
		{
			pClip[k].x = pTemp[v].x;
			pClip[k].y = pTemp[v].y;
			k++;
		}

		inPrev = inCur;
		u = v;
	}
	nPointsNew = k;

	if (nPointsNew < 3)
	{
		FREE(pTemp);
		FREE(pClip);
		*nReturnPoints = 0;
		return NULL;
	}

	FREE(pTemp);
	pTemp = (gf_point *) malloc(nPointsNew * sizeof(gf_point));
	for (v = 0; v < nPointsNew; v++)
	{
		pTemp[v].x = pClip[v].x;
		pTemp[v].y = pClip[v].y;
	}
	FREE(pClip);

	*nReturnPoints = nPointsNew;
	return pTemp;
}
