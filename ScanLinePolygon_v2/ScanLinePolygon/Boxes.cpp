//=========================
//
// Determine Bounding Boxes
//
//=========================

#include <float.h>
#include "gf_types.h"
#include "gfunc.h"

gf_box gfGetBoundingBox(gf_point *point, int nPoints)
{
	// Determine and return a bounding box of a set of point

	gf_box box;
	int v;
	box.ll.x = point[0].x;
	box.ll.y = point[0].y;
	box.ur.x = box.ll.x;
	box.ur.y = box.ll.y;
	for (v = 1; v < ABS(nPoints); v++)
	{
		box.ll.x = AMIN(box.ll.x, point[v].x); 
		box.ll.y = AMIN(box.ll.y, point[v].y); 
		box.ur.x = AMAX(box.ur.x, point[v].x); 
		box.ur.y = AMAX(box.ur.y, point[v].y); 
	}
	return box;
}

gf_box gfGetBoundingBox(gf_contour contour)
{
	// Determine and return a bounding box of a contour

	gf_box box;
	int v;
	
	//////////////////////////////////////////////////////////////////////////
	// contour.point 가 NULL 로 처리되 죽습니다.
	// 일단 예외 처리 하였습니다. chahn : 10.06.22
// 	if(NULL == contour.point)
// 	{
// 		box.ll.x = 0.0;
// 		box.ll.y = 0.0;
// 		box.ur.x = 0.0;
// 		box.ur.y = 0.0;
// 		return box;
// 	}
	//////////////////////////////////////////////////////////////////////////

	box.ll.x = contour.point[0].x;
	box.ll.y = contour.point[0].y;
	box.ur.x = box.ll.x;
	box.ur.y = box.ll.y;
	for (v = 1; v < ABS(contour.nPoints); v++)
	{
		box.ll.x = AMIN(box.ll.x, contour.point[v].x); 
		box.ll.y = AMIN(box.ll.y, contour.point[v].y); 
		box.ur.x = AMAX(box.ur.x, contour.point[v].x); 
		box.ur.y = AMAX(box.ur.y, contour.point[v].y); 
	}
	return box;
}

void gfSetBoundingBox(gf_polygon *p, int c)
{
	// Determine and set a bounding box of the contour no "c" in polygon

	int v;
	if (ABS(p->contour[c].nPoints) == 0) return;

	p->box[c].ll.x = p->contour[c].point[0].x;
	p->box[c].ll.y = p->contour[c].point[0].y;
	p->box[c].ur.x = p->box[c].ll.x;
	p->box[c].ur.y = p->box[c].ll.y;
	for (v = 1; v < ABS(p->contour[c].nPoints); v++)
	{
		p->box[c].ll.x = AMIN(p->box[c].ll.x, p->contour[c].point[v].x); 
		p->box[c].ll.y = AMIN(p->box[c].ll.y, p->contour[c].point[v].y); 
		p->box[c].ur.x = AMAX(p->box[c].ur.x, p->contour[c].point[v].x); 
		p->box[c].ur.y = AMAX(p->box[c].ur.y, p->contour[c].point[v].y); 
	}
}

void gfSetBoundingBox(gf_polygon *p)
{
	// Determine and set bounding boxes of each contour in polygon

	int c;
	for (c = 0; c < ABS(p->nContours); c++)
	{
		gfSetBoundingBox(p, c);
	}
}

gf_box gfGetBoxOfBoxes(gf_box *boxes, int nBoxes)
{
	// Determine and return a bounding box of boxes

	gf_box box;
	int b;

	box.ll.x = boxes[0].ll.x;
	box.ll.y = boxes[0].ll.y;
	box.ur.x = boxes[0].ur.x;
	box.ur.y = boxes[0].ur.y;

	for (b = 1; b < nBoxes; b++)
	{
		box.ll.x = AMIN(box.ll.x, boxes[b].ll.x); 
		box.ll.y = AMIN(box.ll.y, boxes[b].ll.y); 
		box.ur.x = AMAX(box.ur.x, boxes[b].ur.x); 
		box.ur.y = AMAX(box.ur.y, boxes[b].ur.y); 
	}
	return box;
}

gf_box gfGetBoxOfBoxes(gf_box box1, gf_box box2)
{
	// Determine and return a bounding box of two boxes

	gf_box box;

	box.ll.x = AMIN(box1.ll.x, box2.ll.x);
	box.ll.y = AMIN(box1.ll.y, box2.ll.y);
	box.ur.x = AMAX(box1.ur.x, box2.ur.x);
	box.ur.y = AMAX(box1.ur.y, box2.ur.y);

	return box;
}

void gfUpdateBoxOfBoxes(gf_box *box, gf_box *boxes, int nBoxes)
{
	// Update the bounding box "box" by checking new boxes "boxes"

	int b;

	for (b = 0; b < nBoxes; b++)
	{
		box->ll.x = AMIN(box->ll.x, boxes[b].ll.x); 
		box->ll.y = AMIN(box->ll.y, boxes[b].ll.y); 
		box->ur.x = AMAX(box->ur.x, boxes[b].ur.x); 
		box->ur.y = AMAX(box->ur.y, boxes[b].ur.y); 
	}
}

bool gfBoxesOverlap(gf_box b1, gf_box b2)
{
	// Check if two boxes are overlapping
	if (b1.ll.x > b2.ur.x) return false;
	if (b1.ll.y > b2.ur.y) return false;
	if (b2.ll.x > b1.ur.x) return false;
	if (b2.ll.y > b1.ur.y) return false;

	return true;
}

double gfBoxesDistance(gf_box b1, gf_box b2)
{
	// Determine closest distance between two boxes

	// If overlap, return -1;
	if (gfBoxesOverlap(b1, b2))
		return -1.0;

	double dMin = DBL_MAX;
	double d;

	d = b1.ll.x - b2.ur.x;
	if (d >= 0 && d < dMin) dMin = d;

	d = b1.ll.y - b2.ur.y;
	if (d >= 0 && d < dMin) dMin = d;

	d = b2.ll.x - b1.ur.x;
	if (d >= 0 && d < dMin) dMin = d;

	d = b2.ll.y - b1.ur.y;
	if (d >= 0 && d < dMin) dMin = d;

	return dMin;
}
