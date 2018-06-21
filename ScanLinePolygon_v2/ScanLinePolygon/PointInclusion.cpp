// -----------------------------------------------------------------------------
// Point Inclusion Checking
//
// Functions:
//		PointInBox     - Checks whether a point is in a rectangular box
//		PointInContour - Checks whether a point is in a contour
//		PointInPolygon - Checks whether a point is in a polygon that consists of
//						 one or more contours and some of them may be holes.
//
// The functions InContour and InPolygon works well for any shape of contours
// or polygons (convex, concave or spiral). The vertices can be defined in
// either clock-wise or counter-clock-wise.
//
// Written in Nay 08, 2009 by Sungmin Kim
// -----------------------------------------------------------------------------
#include <stdlib.h>
#include "gf_types.h"
#include "gfunc.h"

int CheckScanLineVsSegment(bool check_left, double tiny, double x, double y,
						   double x1, double y1, double x2, double y2);

// -----------------------------------------------------------------------------

bool PointInBox(double x, double y, gf_box box)
{
	if (x < box.ll.x || x > box.ur.x ||
		y < box.ll.y || y > box.ur.y)
		return false;
	else
		return true;
}

bool PointInBox(double x, double y, double xll, double yll, double xur, double yur)
{
	if (x < xll || x > xur || y < yll || y > yur)
		return false; 
	else
		return true;
}

bool PointInContour(double x, double y, int np, gf_point *p,
					double xll, double yll, double xur, double yur)
{
    int  n, k, m;

    n = abs(np);

	// Must have at least 3 points to be a contour

	if (n <= 2) return false;

	// Check with bounding box first

	if (x < xll || x > xur || y < yll || y > yur) return false; 

	// Determine whether to check line segment on the left or right of the check point

	bool check_left;
	if (x <= 0.5 * (xll + xur))
		check_left = true;
	else
		check_left = false;

	// Determine a tiny value to shift an end point slightly to avoid
	// the same end point getting counted twice.

	double tiny;
	tiny = 1.e-4 * (yur - yll);

	// Start with from the indicator flag "false" then toggle it as it passes line segment

	bool in = false;

	// Loop through all line segments

	double x1, y1, x2, y2;

	for (k = 0; k < n; k++)
	{
		x1 = p[k].x;
		y1 = p[k].y;

		if (k+1 == n)
		{
			x2 = p[0].x;
			y2 = p[0].y;
		}
		else
		{
			x2 = p[k+1].x;
			y2 = p[k+1].y;
		}

		// Check horizontal scan line vs the line segment

		m = CheckScanLineVsSegment(check_left, tiny, x, y, x1, y1, x2, y2);

		if (m == 1)
			in = !in;
		else if (m == -1)
			return true;
	};

	return in;
}

bool PointInContour(double x, double y, int np, double *px, double *py,
					double xll, double yll, double xur, double yur)
{
    int  n, k, m;

    n = abs(np);

	// Must have at least 3 points to be a contour

	if (n <= 2) return false;

	// Check with bounding box first

	if (x < xll || x > xur || y < yll || y > yur) return false; 

	// Determine whether to check line segment on the left or right of the check point

	bool check_left;
	if (x <= 0.5 * (xll + xur))
		check_left = true;
	else
		check_left = false;

	// Determine a tiny value to shift an end point slightly to avoid
	// the same end point getting counted twice.

	double tiny;
	tiny = 1.e-4 * (yur - yll);

	// Start with from the indicator flag "false" then toggle it as it passes line segment

	bool in = false;

	// Loop through all line segments

	double x1, y1, x2, y2;

	for (k = 0; k < n; k++)
	{
		x1 = px[k];
		y1 = py[k];

		if (k+1 == n)
		{
			x2 = px[0];
			y2 = py[0];
		}
		else
		{
			x2 = px[k+1];
			y2 = py[k+1];
		}

		// Check horizontal scan line vs the line segment

		m = CheckScanLineVsSegment(check_left, tiny, x, y, x1, y1, x2, y2);

		if (m == 1)
			in = !in;
		else if (m == -1)
			return true;
	};

	return in;
}

bool PointInContour(double x, double y, int np, double *p,
					double xll, double yll, double xur, double yur)
{
	// p:  x1, y1, x2, y2, ...

    int  n, i, k, m;

    n = abs(np);

	// Must have at least 3 points to be a contour

	if (n <= 2) return false;

	// Check with bounding box first

	if (x < xll || x > xur || y < yll || y > yur) return false; 

	// Determine whether to check line segment on the left or right of the check point

	bool check_left;
	if (x <= 0.5 * (xll + xur))
		check_left = true;
	else
		check_left = false;

	// Determine a tiny value to shift an end point slightly to avoid
	// the same end point getting counted twice.

	double tiny;
	tiny = 1.e-4 * (yur - yll);

	// Start with from the indicator flag "false" then toggle it as it passes line segment

	bool in = false;

	// Loop through all line segments

	double x1, y1, x2, y2;

	i = 0;
	for (k = 0; k < n; k++)
	{
		x1 = p[i];
		y1 = p[i+1];
		i += 2;
		if (k+1 == n)
		{
			x2 = p[0];
			y2 = p[1];
		}
		else
		{
			x2 = p[i];
			y2 = p[i+1];
		}

		// Check horizontal scan line vs the line segment

		m = CheckScanLineVsSegment(check_left, tiny, x, y, x1, y1, x2, y2);

		if (m == 1)
			in = !in;
		else if (m == -1)
			return true;
	};

	return in;
}

bool PointInContour(double x, double y, int np, gf_point *p)
{
    int  n, k;

    n = abs(np);

	// Must have at least 3 points to be a contour

	if (n <= 2) return false;

	// Determine X and Y boundary

	double xmin, xmax, ymin, ymax;
	xmin = p[0].x;
	xmax = p[0].x;
	ymin = p[0].y;
	ymax = p[0].y;
	for (k = 1; k < n; k++)
	{
		if (p[k].x < xmin) xmin = p[k].x;
		if (p[k].x > xmax) xmax = p[k].x;
		if (p[k].y < ymin) ymin = p[k].y;
		if (p[k].y > ymax) ymax = p[k].y;
	};

	return PointInContour(x, y, n, p, xmin, ymin, xmax, ymax);
}

bool PointInContour(double x, double y, int np, double *p)
{
	// p:  x1, y1, x2, y2, ...

    int  n, i, k;

    n = abs(np);

	// Must have at least 3 points to be a contour

	if (n <= 2) return false;

	// Determine X and Y boundary

	double xmin, xmax, ymin, ymax;
	i = 0; 
	xmin = p[i];
	xmax = p[i];
	i++;
	ymin = p[i];
	ymax = p[i];
	for (k = 1; k < n; k++)
	{
		i++;
		if (p[i] < xmin) xmin = p[i];
		if (p[i] > xmax) xmax = p[i];
		i++;
		if (p[i] < ymin) ymin = p[i];
		if (p[i] > ymax) ymax = p[i];
	};

	return PointInContour(x, y, n, p, xmin, ymin, xmax, ymax);
}

bool PointInContour(double x, double y, int np, double *px, double *py)
{
    int  n, k;

    n = abs(np);

	// Must have at least 3 points to be a contour

	if (n <= 2) return false;

	// Determine X and Y boundary

	double xmin, xmax, ymin, ymax;
	xmin = px[0];
	xmax = px[0];
	ymin = py[0];
	ymax = py[0];
	for (k = 1; k < n; k++)
	{
		if (px[k] < xmin) xmin = px[k];
		if (px[k] > xmax) xmax = px[k];
		if (py[k] < ymin) ymin = py[k];
		if (py[k] > ymax) ymax = py[k];
	};

	return PointInContour(x, y, n, px, py, xmin, ymin, xmax, ymax);
}

// bool PointInPolygon(double x, double y, gf_polygon *poly)
// {
// 	int i = 0;
// 	bool in = false;
// 
// 	while (i < poly->nContours)
// 	{
// 		int now = i;
// 		i++;
// 		in = false;
// 		if (!poly->hole[now])
// 		{
// 			in = PointInContour(x, y, poly->contour[now].nPoints, poly->contour[now].point);
// 			if (in)
// 			{
// 				while(in && i < poly->nContours && poly->hole[i])
// 				{
// 					if (PointInContour(x, y, poly->contour[i].nPoints, poly->contour[i].point))
// 					{
// 						in = false;
// 						i++;
// 						break;
// 					}
// 					else
// 					{
// 						i++;
// 					}
// 				}
// 			}
// 		}
// 		if (in) return true;
// 	}
// 	return false;
// }

bool PointInPolygon(double x, double y, gf_polygon *poly, int *ContourID)
{
	*ContourID = -1;
	int i;
	for (i = poly->nContours - 1; i >= 0; i--)
	{
		if (PointInContour(x, y, poly->contour[i].nPoints, poly->contour[i].point))
		{
			*ContourID = i;
			if (poly->hole[i])
				return false;
			else
				return true;
		}
	}
	return false;
}

bool PointInPolygon(double x, double y, gf_polygon *poly, int start, int nContours, int *ContourID)
{
	*ContourID = -1;
	int i;
	for (i = start + nContours - 1; i >= start; i--)
	{
		if (PointInContour(x, y, poly->contour[i].nPoints, poly->contour[i].point))
		{
			*ContourID = i;
			if (poly->hole[i])
				return false;
			else
				return true;
		}
	}
	return false;
}

bool PointInPolygon(double x, double y, gf_polygon *poly)
{
	int i;
	return PointInPolygon(x, y, poly, &i);
}

bool PointInPolygon(double x, double y, gf_polygon *poly, int start, int nContours)
{
	int i;
	return PointInPolygon(x, y, poly, start, nContours, &i);
}

int CheckScanLineVsSegment(bool check_left, double tiny, double x, double y,
						   double x1, double y1, double x2, double y2)
{
	// check_left:  True to check left hand side of the check point,
	//              False to check right hand side.
	// tiny:        A tiny value to shift an end point slightly to avoid the same end point
	//              getting counted twice.
	// x, y:        Check point.
	// x1,y1,x2,y2: Start and end points of the line segment.
	//
	// Returns: -1 if the check point is on the line segment
	//           0 if the scan line does not crosses the line segment
	//           1 if the scan line crosses the line segment

	// Check only one side of the check point
	if (check_left)
	{
		if (x1 > x && x2 > x) return 0;
	}
	else
	{
		if (x1 < x && x2 < x) return 0;
	}

	// The check point is exactly at one of the two end points of the line
	if ((x1 == x && y1 == y) || (x2 == x && y2 == y)) return -1;

	// Horizontal line
	if (y1 == y2 && y1 == y)
	{
		if (x1 > x && x2 > x) return 0;
		if (x1 < x && x2 < x) return 0;
		return -1;
	};

	// Vertical line
	if (x1 == x2 && x1 == x)
	{
		if (y1 > y && y2 > y) return 0;
		if (y1 < y && y2 < y) return 0;
		return -1;
	};

	// If y of either of two end points of the line is the same as y of the check point,
	// shift th end point slightly to avoid the same end point getting counted twice.
	double yy1 = y1;
	double yy2 = y2;
	if (yy1 == y) yy1 += tiny;
	if (yy2 == y) yy2 += tiny;

	// Check if both end points are above or below the check point
	if (yy1 > y && yy2 > y) return 0;
	if (yy1 < y && yy2 < y) return 0;

	// Calculate X coordinate of the point at which the horizontal line from the check point
	// and the line segment meet
	double xp = x1 + (x2 - x1) * (y - yy1) / (yy2 - yy1);

	// The check point is exactly on the line
	if (xp == x) return -1;

	// The horizontal scan line crosses the line segment
	if (check_left)
	{
		if (xp < x) return 1;
	}
	else
	{
		if (xp > x) return 1;
	}

	return 0;
}

bool PointInBetween(double x, double y, double x1, double y1, double x2, double y2, double epsilonFactor)
{
	gf_box box;
	double sx, sy;

	box.ll.x = AMIN(x1, x2);
	box.ll.y = AMIN(y1, y2);
	box.ur.x = AMAX(x1, x2);
	box.ur.y = AMAX(y1, y2);

	sx = box.ur.x - box.ll.x;
	sy = box.ur.y - box.ll.y;
	double smax = sx > sy ? sx : sy;

	if (smax == 0) return false;

	if (sy <= epsilonFactor * smax) {
		box.ll.y -= sx * epsilonFactor;
		box.ur.y += sx * epsilonFactor;
	}
	else if (sx <= epsilonFactor * smax) {
		box.ll.x -= sy * epsilonFactor;
		box.ur.x += sy * epsilonFactor;
	}

	return PointInBox(x, y, box);
}

bool PointInBetween(gf_point p, gf_point v1, gf_point v2, double epsilonFactor)
{
	return PointInBetween(p.x, p.y, v1.x, v1.y, v2.x, v2.y, epsilonFactor);
}

bool PointInBetween(double x, double y, double x1, double y1, double x2, double y2)
{
	double epsilonFactor = 1.0e-6;
	return PointInBetween(x, y, x1, y1, x2, y2, epsilonFactor);
}

bool PointInBetween(gf_point p, gf_point v1, gf_point v2)
{
	return PointInBetween(p.x, p.y, v1.x, v1.y, v2.x, v2.y);
}
