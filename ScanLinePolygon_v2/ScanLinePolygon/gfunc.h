#ifndef _GFUNC_H_
#define _GFUNC_H_

#include "gf_types.h"
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#define LINE_JOINT_ARC_DIV_ANGLE		15
#define LINE_JOINT_ARC_DIV_MAX_POINT	13

// PointNearStrip function return parameters
#define GF_PNS_FAR_FROM_STRIP		0
#define GF_PNS_AT_BEGIN_POINT		1
#define GF_PNS_AT_END_POINT			2
#define GF_PNS_ON_CENTER_LINE		3
#define GF_PNS_WITHIN_STRIP			4
#define GF_PNS_NEAR_BEGIN_POINT		5
#define GF_PNS_NEAR_END_POINT		6

//==============
// Point2Line
//==============
// PointOnLine
//		Arguments: a point to check, starting and ending points of a line.
//		If eps is given, the checking will be done with the tolerance of eps.
bool PointOnLine(double x, double y, double x1, double y1, double x2, double y2);
bool PointOnLine(gf_point p, gf_point v1, gf_point v2);
bool PointOnLine(double x, double y, double x1, double y1, double x2, double y2, double eps);
bool PointOnLine(gf_point p, gf_point v1, gf_point v2, double eps);

// PointProject2Line
//		Arguments: a point, starting and ending points of a line, projected coordinates (return)
//		It returns true if the projected point is on the line segment. Otherwise, return false.
bool PointProject2Line(double x, double y, double x1, double y1, double x2, double y2, double *px, double *py);
bool PointProject2Line(gf_point p, gf_point v1, gf_point v2, gf_point *proj);

// Point2LineDistance
//		Arguments: a point, starting and ending points of a line, projected distance flag.
//		If the projected_distance is true, it will be the projected distance to an infinite
//		line that connects the starting ending points.
//		If the projected_distance is false, it will be the closest (may not be a projected)
//		distance to the line segment that connects the starting and ending points.
double Point2LineDistance(double x, double y, double x1, double y1, double x2, double y2, bool projected_distance, double *px, double *py);
double Point2LineDistance(double x, double y, double x1, double y1, double x2, double y2, bool projected_distance);
double Point2LineDistance(gf_point p, gf_point v1, gf_point v2, bool projected_distance, gf_point *pxy);
double Point2LineDistance(gf_point p, gf_point v1, gf_point v2, bool projected_distance);

// PointMirroredAgainstLine
//		Arguments: a point, starting and ending points of a line, mirrored coordinates (return).
//		It returns true if the projected point is on the line segment. Otherwise, return false.
bool PointMirroredAgainstLine(double x, double y, double x1, double y1, double x2, double y2, double *px, double *py);
bool PointMirroredAgainstLine(gf_point p, gf_point v1, gf_point v2, gf_point *mirroredPt);

// nPointNearStrip
//		Arguments: reference point, starting and ending points of a line, line width, distance to the line, nearest point.
//		It returns GF_PNS_... parameters.
int nPointNearStrip(double x, double y, double x1, double y1, double x2, double y2, double dWidth, double *dDistance, double *px, double *py);
int nPointNearStrip(gf_point p, gf_point v1, gf_point v2, double dWidth, double *dDistance, gf_point *nearestPt);

//===============
// PointInclusion
//===============
bool PointInBox(double x, double y, gf_box box);
bool PointInBox(double x, double y, double xll, double yll, double xur, double yur);

 bool PointInContour(double x, double y, int np, gf_point *p,
								double xll, double yll, double xur, double yur);
 bool PointInContour(double x, double y, int np, double *p,
								double xll, double yll, double xur, double yur);
 bool PointInContour(double x, double y, int np, double *px, double *py,
								double xll, double yll, double xur, double yur);
 bool PointInContour(double x, double y, int np, gf_point *p);
 bool PointInContour(double x, double y, int np, double *p);
 bool PointInContour(double x, double y, int np, double *px, double *py);

 bool PointInPolygon(double x, double y, gf_polygon *poly, int *ContourID);
 bool PointInPolygon(double x, double y, gf_polygon *poly);
 bool PointInPolygon(double x, double y, gf_polygon *poly, int start, int nContours, int *ContourID);
 bool PointInPolygon(double x, double y, gf_polygon *poly, int start, int nContours);

 bool PointInBetween(double x, double y, double x1, double y1, double x2, double y2, double epsilonFactor);
 bool PointInBetween(gf_point p, gf_point v1, gf_point v2, double epsilonFactor);
 bool PointInBetween(double x, double y, double x1, double y1, double x2, double y2);
 bool PointInBetween(gf_point p, gf_point v1, gf_point v2);
	
//==========
// Line2Line
//==========

// dCrossProduct
//		Arguments: line1, line2
//      It returns cross product of the two vectors defined by line segments
 double dCrossProduct(double x1s, double y1s, double x1e, double y1e,
								 double x2s, double y2s, double x2e, double y2e);
 double dCrossProduct(gf_point *L1, gf_point *L2);

// dLineOverlap
//		Arguments: line1 (reference), line2 (check), parallel check angle,
//				   overlapped section on the reference line (output)
//      It returns cross product of the two vectors defined by line segments

 double dLineOverlap(double x1s, double y1s, double x1e, double y1e,
								double x2s, double y2s, double x2e, double y2e,
								double dParallelCheckAngle,
								double *x1, double *y1, double *x2, double *y2);
 double dLineOverlap(gf_point p1s, gf_point p1e,
								gf_point p2s, gf_point p2e,
								double dParallelCheckAngle,
								gf_point *pFrom, gf_point *pTo);

// bTwoLinesParallel
//		Arguments: line1, line2
//      It returns true, if they are parallel
//                 false, if they are not parallel
 bool bTwoLinesParallel(double x1s, double y1s, double x1e, double y1e,
								   double x2s, double y2s, double x2e, double y2e);
 bool bTwoLinesParallel(gf_point *L1, gf_point *L2);

// MergeTwoLines
//		Arguments: line1, line2
//      It returns 0, if they are not connected or collinear
//				 > 0, if they are connected and collinear
//					  11: new line connects from the 1st point of line1 to the 1st point of line2
//					  12: new line connects from the 1st point of line1 to the 2nd point of line2
//					  21: new line connects from the 2nd point of line1 to the 1st point of line2
//					  22: new line connects from the 2nd point of line1 to the 2nd point of line2
 int MergeTwoLines(double x1s, double y1s, double x1e, double y1e,
							  double x2s, double y2s, double x2e, double y2e, double tolerance = 0.0);
 int MergeTwoLines(gf_point *L1, gf_point *L2, double tolerance = 0.0);

// bTwoLinesCross
//		Arguments: line1, line2, crossing point (output)
//      It returns false, if they are parallel
//                 true, if they are not parallel. It outputs the crossing point.
 bool bTwoLinesCross(double x1s, double y1s, double x1e, double y1e,
								double x2s, double y2s, double x2e, double y2e,
								double *x, double *y);
 bool bTwoLinesCross(gf_point *L1, gf_point *L2, double *x, double *y);
 bool bTwoLinesCross(gf_line L1, gf_line L2, double *x, double *y);

// nTwoLines
//		Arguments: line1, line2, crossing point (output)
//		It returns -1: parallel
//					0: two lines meet (the crossing point is not on either line)
//					1: two lines meet (the crossing point is on the first line only)
//					2: two lines meet (the crossing point is on the second line only)
//					3: two lines meet (the crossing point is on both lines)
 int nTwoLines(double x1s, double y1s, double x1e, double y1e,
						  double x2s, double y2s, double x2e, double y2e,
						  double *x, double *y, double epsilonFactor);
 int nTwoLines(gf_point *L1, gf_point *L2, double *x, double *y, double epsilonFactor);
 int nTwoLines(double x1s, double y1s, double x1e, double y1e,
						  double x2s, double y2s, double x2e, double y2e,
						  double *x, double *y);
 int nTwoLines(gf_point *L1, gf_point *L2, double *x, double *y);

// dLine2LineDistance
//		Arguments: lineA, lineB, projected distance only (input)
//				   Coordinates of the two closest points (output)
//		It returns a distance between the two points.
//		If the two lines are crossing, it returns a negative value for the distance.
 double dLine2LineDistance(double dXA1, double dYA1, double dXA2, double dYA2,
									  double dXB1, double dYB1, double dXB2, double dYB2,
									  bool projected_line);
 double dLine2LineDistance(double dXA1, double dYA1, double dXA2, double dYA2,
									  double dXB1, double dYB1, double dXB2, double dYB2,
									  bool projected_line,
									  double *dXA, double *dYA,
									  double *dXB, double *dYB);

// vJoinLines
//		Join two lines with an arc or at the crossing point
//		Returns	true if the two lines were crossed.
//				false if the two lines were not crossed.
 bool bJoinLines(gf_point arcOrigin, gf_point *pLine1, gf_point *pLine2,
							bool bArcJoin, gf_point *pPoint, int *nPoints);
 bool bJoinLines(gf_point arcOrigin, gf_point *pLine1, gf_point *pLine2,
							int nJoinType, gf_point *pPoint, int *nPoints);

// bGetParallelLine
//		Arguments: line1, line2 (output), distance, extension, epsilon
//		It returns false if the length is less than epsilon.
//                  true if the length is greater than epsilon
//		Extends the line by the amount of extension
//		It makes a line2 on the right hand side of the line1 apart "distance"
 bool bGetParallelLine(double xs, double ys, double xe, double ye,
								  double *pxs, double *pys, double *pxe, double *pye,
								  double distance, double extended, double epsilon);
 bool bGetParallelLine(gf_point *line, gf_point *rline,
								  double distance, double extended, double epsilon);
 bool bGetParallelLine(gf_line line, gf_line rline,
								  double distance, double extended, double epsilon);

// nCleanSelfCrossingLines
//		Arguments: pPath, nPoints
//		It returns new number of points after clean up.
 int nCleanSelfCrossingLines(gf_point *pPath, int nPoints);

// pGetPathApart
//		Arguments: point list of the reference path, (+/-)distance, line join method
//		Outputs:   Number of points of the result line
//		Returns:   Point list of the result line
//		Makes a path from the reference line apart the given distance on the right or left
 gf_point *pGetPathApart(gf_point *pLine, int nPoints,
									double dDistance, int nJoinType,
									int *nRetDivPoints);
// Determine angle between lines
 double dAngleFromEast(double dX1, double dY1, double dX2, double dY2);
 double dAngleTwoLines(double dXA1, double dYA1, double dXA2, double dYA2,
								  double dXB1, double dYB1, double dXB2, double dYB2);
 double dAngle3Points(double dX1, double dY1, double dX2, double dY2, double dX3, double dY3);

//====
// Arc
//====

// Divide arc with a given division angle
 gf_point *pDivideArc(double dXC, double dYC, double dX1, double dY1, double dX2, double dY2,
								 bool bCCW, double dDivAngle, int nPointsMin, int nPointsMax, int *nRetDivPoints);
 gf_point *pDivideArc(gf_point pCenter, gf_point p1, gf_point p2,
								 bool bCCW, double dDivAngle, int nPointsMin, int nPointsMax, int *nRetDivPoints);

// Divide arc with a given number of division points
 gf_point *pDivideArc(double dXC, double dYC, double dX1, double dY1, double dX2, double dY2,
								 bool bCCW, int nDivPoints);
 gf_point *pDivideArc(gf_point pCenter, gf_point p1, gf_point p2, bool bCCW, int nDivPoints);

// Determines an arc center from 3 points one an arc. Also outputs radius and arc direction (CCW or CW)
 gf_point getArcCenter(double dXstart, double dYstart, double dXend, double dYend,	// arc defined by 3 points
								  double dXmid, double dYmid);
 gf_point getArcCenter(double dXstart, double dYstart, double dXend, double dYend,	// arc center defined by offsets
								  double dXoffset, double dYoffset, bool bEanble360);			// from the starting point
 gf_point getArcCenter(double x1s, double y1s, double x1e, double y1e,	// start and end points of the first tangent line. Tangent at the start point.
								  double x2s, double y2s, double x2e, double y2e);	// start and end points of the second tangent line. Tangent at the start point.

 gf_box getArcBox(double dXC, double dYC, double dX1, double dY1, double dX2, double dY2, bool bCCW);

//========
// Contour
//========

// dContourArea
 double dContourArea(double *x, double *y, int npoints);
 double dContourArea(gf_point *p, int npoints);

// bClockwise
//		returns true if clockwise contour, false if counter-clockwise
 bool bClockwise(double *x, double *y, int npoints);
 bool bClockwise(gf_point *p, int npoints);

// RemoveCollinearPoint
 int nRemoveCollinearPoint(double *x, double *y, int npoints, bool bContour, double tolerance = 0.0);
 int nRemoveCollinearPoint(gf_point *p, int npoints, bool bContour, double tolerance = 0.0);
//		returns a new number of points after collinear points are removed
//		x, y or p: arrays of x and y coordinates
//		npoints: number of points before collinear point removal
//		bContour: indicate the points are for contours instead of poly lines

 gf_point *pRemoveCollinearPoint(gf_point *point, int *nPoints, bool bContour, double tolerance = 0.0);
//		returns a new array after remove all collinear points
//		point: original point list
//		nPoint: number of points (this is updated with the new number of points)
//		bContour: indicate the points are for contours instead of poly lines

 void vRemoveCollinearPoint(gf_polygon *p, double tolerance = 0.0);
//		removes collinear points within each contour in the polygon

extern "C" double	TINY_HOLE_AREA;
extern "C" int		TINY_HOLE_POINTS;

 void	vSetTinyAreaSize(double tinyArea, int nPoints);
 int		nGetTinyAreaPoints();
 double	dGetTinyAreaSize();

 gf_polygon *pRemoveTinyHoles(gf_polygon *p);
//		removes tiny hole contours

// Check whether two contours intersect
 bool bContoursIntersect(gf_point *pContour1, int nPoints1, gf_point *pContour2, int nPoints2, bool bCheckBoundingBox = true);
 bool bContoursIntersect(gf_point *pContour1, int nPoints1, gf_point *pContour2, int nPoints2, gf_box box1, gf_box box2);

// Check whether a contour is completely inside of another contour
 bool bContourInContour(gf_point *pInside, int nInside, gf_point *pOutside, int nOutside);
 bool bContourInContour(gf_polygon *pInside, int nInID, gf_polygon *pOutside, int nOutID);
 bool bContourInContour(gf_polygon *pInside, gf_polygon *pOutside);

// Find an identical contour
 int nFindContour(gf_contour contour, bool bHole, gf_polygon *pPoly, int nStart, int nEnd, bool bCheckType);
 int nFindContour(gf_polygon *pFind, int nContour, gf_polygon *pPoly, int nStart, int nEnd, bool bCheckType);

// Remove redundant contours in the same polygon
 gf_polygon *pRemoveRedundantContours(gf_polygon *pPoly);

// Shifts (translates) contour
 void vShiftContour(gf_polygon *pPoly, int iContour, int nContours, double dShiftX, double dShiftY);
 void vShiftContour(gf_point *pPoint, int nPoints, double dShiftX, double dShiftY);

// Rotates contour
 void vRotateContour(gf_polygon *pPoly, int iContour, int nContours, double dRotDegreeCCW, double dRefX, double dRefY);
 void vRotateContour(gf_point *pPoint, int nPoints, double dRotDegreeCCW, double dXref, double dYref);
 void vRotateContour(gf_point *pPoint, int nPoints, double dRotDegreeCCW,int nRefPoint);

// Mirrors contour
 void vMirrorContour(gf_polygon *pPoly, int iContour, int nContours, int nMirror, double dXref, double dYref);
 void vMirrorContour(gf_point *pPoint, int nPoints, int nMirror,	double dXref, double dYref);
 void vMirrorContour(gf_point *pPoint, int nPoints, int nMirror,	int nRefPoint);

//===========================================
// Converts Geometric Description to Contours
//===========================================

// Converts a box to a contour.
// Can be used for the Gerber D-Code aperture type "R" (rectangle)
// and for the M-code primitive number "21" (center method line definition) and "22" (lower left method line definition).
 gf_point *pBox2Contour(double dXcenter, double dYcenter, double dWidth, double dHeight,	double dRotDegree, int mirror);
 gf_point *pBox2Contour(double dXcenter, double dYcenter, double dWidth, double dHeight, int mirror,	double dRotDegree);
 gf_point *pBox2Contour(double dXcenter, double dYcenter, double dWidth, double dHeight,	double dRotDegree);
 gf_point *pBox2Contour(double dXcenter, double dYcenter, double dWidth, double dHeight);

// Converts a rounded-box (edge-chamfered box) to a contour.
 gf_point *pRoundedBox2Contour(double dXcenter, double dYcenter, double dWidth, double dHeight,
										  double dChamferRadius, int nChamferDivPoints, double dRotDegree, int nMirror);
 gf_point *pRoundedBox2Contour(double dXcenter, double dYcenter, double dWidth, double dHeight,
										  double dChamferRadius, int nChamferDivPoints, int nMirror, double dRotDegree);
 gf_point *pRoundedBox2Contour(double dXcenter, double dYcenter, double dWidth, double dHeight,
										  double dChamferRadius, int nChamferDivPoints, double dRotDegree);
 gf_point *pRoundedBox2Contour(double dXcenter, double dYcenter, double dWidth, double dHeight,
										  double dChamferRadius, int nChamferDivPoints);

// Converts a circle to a contour.
// Can be used for the Gerber D-Code aperture type "C" (circle)
// and for the M-code primitive numbers "1" (circle) and "5" (polygon)
 gf_point *pCircle2Contour(double dXcenter, double dYcenter, double dDia, int nSegments,	double dStartDegree);
 gf_point *pCircle2Contour(double dXcenter, double dYcenter, double dDia, int nSegments,	bool bHalfAngleShift);
 gf_point *pCircle2Contour(double dXcenter, double dYcenter, double dDia, int nSegments);

// Converts an ellipse to a contour.
 gf_point *pEllipse2Contour(double dXcenter, double dYcenter, double a, double b, double rot, int nSegments);

// Converts a line segment to a contour.
// Can be used for the Gerber M-code primitive numbers "2" or "20" (vector format line definition)
 gf_point *pLine2Contour(gf_point *pLine, double dWidth, int nRoundCapType, int nRoundCapPoints, int *nRetDivPoints);
 gf_point *pLine2Contour(double dX1, double dY1, double dX2, double dY2,
									double dWidth, int nRoundCapType, int nRoundCapPoints, int *nRetDivPoints);

// Converts a tapered line segment to a contour.
 gf_point *pTapered2Contour(double dX1, double dY1, double dX2, double dY2, double dStartWidth, double dEndWidth,
									   int nRoundCapType, int nRoundCapPoints, int *nRetDivPoints);
 gf_point *pTapered2Contour(gf_point *pLine, double dStartWidth, double dEndWidth,
									   int nRoundCapType, int nRoundCapPoints, int *nRetDivPoints);

// Converts a path (polyline) to a contour.
// Can be used to expand or to reduce a closed polyline (contour)
 gf_point *pPath2Contour(gf_point *pLine, int nPoints, double dWidth,
									int nJoinType, bool bRoundCap, int nMethod, int *nRetDivPoints);

// Converts a Moire circle to a set of contours.
// Gerber Aperture Macro primitive number 6
 gf_polygon *pMoire2Contour(double dXC, double dYC, double dOutDia, double dRingThick, double dGapThick, int nRings,
									   double dHairThick, double dHairLength, double dRotDegree, int nDivPoints);

// Converts a wagon wheel shape (thermal relief) to a set of contours.
// Gerber Aperture Macro primitive number 7
 gf_polygon *pThermal2Contour(double dXC, double dYC, double dOutDia, double dInDia,
										 double dHairWidth, double dRotDegree, int nDivPoints);

// Converts an arc to a contour.
 gf_point *pArc2Contour(double dXC, double dYC, double dX1, double dY1, double dX2, double dY2,
								   bool bCCW, int nDivPoints, double dWidth, int nRoundCapPoints, int *nRetDivPoints);

//======
// Boxes
//======

 gf_box gfGetBoundingBox(gf_point *point, int nPoints);
 gf_box gfGetBoundingBox(gf_contour contour);

 void   gfSetBoundingBox(gf_polygon *p, int c);
 void   gfSetBoundingBox(gf_polygon *p);

 gf_box gfGetBoxOfBoxes(gf_box *boxes, int nBoxes);
 gf_box gfGetBoxOfBoxes(gf_box box1, gf_box box2);

 void   gfUpdateBoxOfBoxes(gf_box *box, gf_box *boxes, int nBoxes);

 bool   gfBoxesOverlap(gf_box b1, gf_box b2);
 double gfBoxesDistance(gf_box b1, gf_box b2);

 //===========================
 // Polygon boolean Operations
 //===========================

 typedef enum	// Polygon boolean operation types
 {
	 GF_UNION,	// Union
	 GF_SUB,	// Subtraction
	 GF_INT,	// Intersection
	 GF_XOR		// Exclusive OR
 } gf_boolean;

 bool gfPolygonBool(gf_boolean operation,
	 gf_polygon *target_polygon,
	 gf_polygon *source_polygon,
	 gf_polygon *result_polygon);

 bool gfPolygonBoolWithoutgfFreePolygon(gf_boolean operation,
	 gf_polygon *target_polygon,
	 gf_polygon *source_polygon,
	 gf_polygon *result_polygon);

 void gfFreePoints(gf_point *&point);

 void gfFreePolygon(gf_polygon *&polygon, bool bFreeContoursToo = true);

 void gfFreeContours(gf_polygon *&polygon);

 void gfAllocContour(gf_polygon *polygon, int contour_number, int nPoints);

 gf_polygon	*gfAllocPolygon (int nContours);

 void		gfAllocPolygon (gf_polygon *p, int nContours);

 void		gfReAllocPolygon (gf_polygon *p, int moreContours);

 gf_polygon *gfAssignContours2Polygon(gf_polygon *source, int c, int nContours);
 gf_polygon *gfCopyContours2Polygon(gf_polygon *source, int c, int nContours);

 void		gfAttachContours2Polygon(gf_polygon *source, int from, int nContours, gf_polygon *target);
 void		gfAppendContours2Polygon(gf_polygon *source, int from, int nContours, gf_polygon *target);

 gf_polygon	*gfPolygonMerge(gf_polygon *target, gf_polygon *source, bool bFreeSource = true);

 gf_polygon	*gfPolygonMergeAll(gf_polygon *p);

 gf_polygon	*SeparateHolesFromContour(gf_polygon *p);

 extern "C" int MEM_ALLOC_CHUNK;
 extern "C" int CONTOUR_MERGE_CHUNK;

 void gfSetChunkSize(int memChunk, int contourMergeChunk);
 void gfSetChunkSize (int nContours);
 int gfGetMemChunkSize ();
 int gfGetMergeChunkSize ();

 bool PolygonsIntersect(gf_polygon *pPoly1, gf_polygon *pPoly2);

 bool PolygonInPolygon(gf_polygon *pObjPoly, double xOffSet, double yOffSet, double dRot, double xLoc, double yLoc,
	 gf_polygon *pPoly, bool bCheckInside, double *dX, double *dY);

 double Polygon2PolygonDistance(gf_polygon *pPoly1,gf_polygon *pPoly2, double dClearance,
	 double *x1, double *y1, double *x2, double *y2);

//==================
// Contour Hierarchy
//==================

// Hierarchically orders all contours in a polygon
 gf_polygon *gfContourHierarchy(gf_polygon *pPoly);

// =====================================
// Contour Expand or Shrink (Abridgment)
// =====================================

// Expand or shrink contours
 gf_polygon *pExpandContour(gf_polygon *pPoly, int iContour, int nContours,
									   double dExpand, bool bTopContourOnly, int nJoinType);
 gf_polygon *pExpandContour(gf_polygon *pPoly, int iContour, int nContours,
									   double dExpand, bool bTopContourOnly, bool bArcJoin);

// =========================
// Contour Clipping by a Box
// =========================

 gf_point *pContourBoxClip(gf_point *pPoints, int nPoints, gf_box box, int *nReturnPoints);

// ===================
// Net Route Functions
// ===================

 int nBestRouteSearch(gf_point curPoint, int nActivePoints, gf_point *pActivePoints,
								 int iAnchor, int iConnectTo, bool bOffGridPin45,
								 double dGridSize, double dGridShiftX, double dGridShiftY,
								 int *nPointsBehind, gf_point *pPointsBehind,
								 int *nPointsAhead, gf_point *pPointsAhead);

 void vRatsConnections(int nPins, gf_point *Pins, int *pinStatus);

// ===================
// Sorting Functions
// ===================

 void vQuickSort(bool bAscend, int ibegin, int iend, double *dA, int *idx);
 void vQuickSort(bool bAscend, int ibegin, int iend, int *dA, int *idx);

 void vQuickSort(bool bAscend, int ibegin, int iend, double *dA);
 void vQuickSort(bool bAscend, int ibegin, int iend, int *dA);

// ===================
// B-spline Functions
// ===================

 void gfFreeBsplineKnots(double *&u);
 void gf3FreePoints(gf3_point *&point);

 gf3_point *pBspline(int nv, int p, int nc, gf3_point *v);
 gf_point *pBspline(int nv, int p, int nc, gf_point *v);

 gf3_point *pBsplineWithKnots(int nv, int p, int nu, int nc,
										 gf3_point *v, double *u, int *knotMultiples);
 gf3_point *pBsplineWithKnotsClosed(int nv, int p, int nu, int nc,
										 gf3_point *v, double *u, int *knotMultiples);

#endif