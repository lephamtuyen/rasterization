#include "TUYENPLE_Polygon.h"
#include "TUYENPLE_RectD.h"
#include "TUYENPLE_Math.h"
#include "TUYENPLE_Line.h"
#include "TUYENPLE_Arc.h"
#include <math.h>
#include <float.h>
#include "gf_types.h"
#include "gfunc.h"

Polygon::Polygon(void) {
	m_pBoundary = NULL;
	m_pCuttingArea = NULL;
}

Polygon::~Polygon(void) {
	delete m_pBoundary;
	delete m_pCuttingArea;
}

Polygon::Polygon(const Polygon& src)
{
	*this = src ;	
}

void Polygon::SetRect()
{
	double left = DBL_MAX;
	double top = DBL_MIN;
	double right = DBL_MIN;
	double bottom = DBL_MAX;

	for (unsigned int i = 0; i < m_pBoundary->size(); i++)
	{
		RectD rect = m_pBoundary->at(i)->GetRect();
		
		left = min(left, rect.GetLeft());
		top = max(top, rect.GetTop());
		right = max(right, rect.GetRight());
		bottom = min(bottom, rect.GetBottom());
	}

	DrawObj::SetRect(left, top, right, bottom);
}

bool Polygon::IsInner(PointD Pt, const bool& bPtOnBoundaryIsInner  )
{
	if(getBoundary()->size() == 0)
		return false;
	int nIsPolygonInner, nIsCuttingAreaInner;

	nIsPolygonInner = IsInner(Pt, m_pBoundary);
	if (0 == nIsPolygonInner) //Outer
		return false;

	if (-1 == nIsPolygonInner)//On the boundary
		return (true == bPtOnBoundaryIsInner ? true : false);

	//Inner
	int nCuttingAreSize = m_pCuttingArea->size();
	if (0 < nCuttingAreSize)
	{
		for (int i=0; i<nCuttingAreSize; i++)
		{	
			Polygon* pCuttingPolygon;
			pCuttingPolygon = m_pCuttingArea->at(i);

			nIsCuttingAreaInner = IsInner(Pt,pCuttingPolygon->getBoundary());
			if (0 != nIsCuttingAreaInner)
				break;
		}
		if (1 == nIsCuttingAreaInner)//Inside cutting areas
			return false;
		if (-1 == nIsCuttingAreaInner)//On cutting area boundary
			return (true == bPtOnBoundaryIsInner ? true : false);

	}	
	return true;	
}


// 0 : outer ; 1 : inner ; -1 : on boundary
int Polygon::IsInner(PointD Pt, const DrawObjArray* const pBoundaryObjs)
{
	//Point-location problems, "Computational Geometry", p41	
	int i;

	int nObjType;
	DrawObj* pCurrentObj;
	Line* pLineSeg;
	Line  LineSeg;
	Arc* pArc;

	for (i=0; i<pBoundaryObjs->size(); i++)
	{
		pCurrentObj = pBoundaryObjs->at(i);
		nObjType = pCurrentObj->GetType();

		if (CPGEOM_LINE == nObjType)
		{
			pLineSeg = (Line*)pCurrentObj;
			if(pLineSeg->Is_Pt_On_LineSeg_Math(Pt))
			{
				return -1;
			}
		}
		else if (CPGEOM_ARC == nObjType)
		{
			pArc = (Arc*)pCurrentObj;
			if (pArc->Is_Pt_On_Arc(Pt))
			{
				return -1;
			}	
		}
	}

	PointDArray	pointArray;
	ConvertToPointArray(pointArray, pBoundaryObjs);

	std::vector<gf_point> v;
	for(i = 0; i < pointArray.size();i++)
	{
		gf_point point;
		point.x = pointArray.at(i).GetX();
		point.y = pointArray.at(i).GetY();
		v.push_back(point);
	}

	if(PointInContour(Pt.GetX(), Pt.GetX(), v.size(), (double*)&v.front() ))
		return TRUE;

	return FALSE;
} 

void Polygon::ConvertToPointArray(PointDArray &PointArray, const DrawObjArray *pDrawObjArray, double dDivision /*= 10*/, int UnitType /*= 0*/)
{
	int		i = 0, j = 0;
	int		nDrawObjCount = 0;
	int		nBeforeType = -1;
	int		nFirstObjectType = -1;
	DrawObj	*pDrawObj = NULL;
	Arc		*pArc = NULL;
	Line	*pLineSeg = NULL;
	Line	*pArcLineSeg = NULL;
	LineArray	*pLineSegArray = NULL;

	PointD		StartPoint;
	PointD		EndPoint;

	for (i = 0; i < pDrawObjArray->size(); i++)
	{
		pDrawObj = pDrawObjArray->at(i);
		nDrawObjCount++;

		if (CPGEOM_LINE == pDrawObj->GetType())
		{
			pLineSeg = (Line*)pDrawObj;

			StartPoint = pLineSeg->GetStartPoint();
			EndPoint = pLineSeg->GetEndPoint();

			if(nBeforeType != -1)
			{
				PointArray.erase(PointArray.begin() + PointArray.size() - 1);
			}

			PointArray.push_back(StartPoint);
			PointArray.push_back(EndPoint);
		}
		else
		{
			pArc = (Arc	*)pDrawObj;
			double dAngleDivision;
			switch(UnitType)
			{
			case 0:	//		degree
				dAngleDivision = dDivision;
				break;

			case 1:	//		mm(length)
				dAngleDivision = 360.0/(2*pArc->GetRadius()*PI/dDivision);
				break;
			}

			if(pArc->GetStartAngle() == pArc->GetEndAngle()
				&& pArc->GetRadius() != 0.0)
			{
				Arc	*pTempArc = (Arc *)pArc->clone();
				pTempArc->SetEndAngle(pArc->GetEndAngle() + 360.0);

				if(dAngleDivision > 179)
				{
					pLineSegArray = pTempArc->GetLineSegArray(30);
				}
				else
				{
					int				nSeperateCount = 360.0 / dAngleDivision;
					double			dRemaindCount = fmod(360.0, dAngleDivision);

					if(dRemaindCount > 0.0)
					{
						nSeperateCount += 1;
					}
					dAngleDivision = 360.0/ nSeperateCount;

					pLineSegArray = pTempArc->GetLineSegArrayForCircle(dAngleDivision);
				}

				delete pTempArc;
			}
			else
			{
				if(dAngleDivision > pArc->GetSweepAngle())
				{
					int				nSeperateCount = pArc->GetSweepAngle() / 30.0;
					double			dRemaindCount = fmod(pArc->GetSweepAngle(), 30.0);

					if(dRemaindCount > 1.0)
					{
						nSeperateCount += 1;
					}

					dAngleDivision = pArc->GetSweepAngle()/ nSeperateCount;
					pLineSegArray = pArc->GetLineSegArray(dAngleDivision);
				}
				else
				{
					int				nSeperateCount = pArc->GetSweepAngle() / dAngleDivision;
					double			dRemaindCount = fmod(pArc->GetSweepAngle(), dAngleDivision);

					if(dRemaindCount > 0.0)
					{
						nSeperateCount += 1;
					}

					dAngleDivision = pArc->GetSweepAngle()/ nSeperateCount;
					pLineSegArray = pArc->GetLineSegArray(dAngleDivision);
				}
			}

			for(j = 0; j < pLineSegArray->size(); j++)
			{
				pArcLineSeg = pLineSegArray->at(j);
				StartPoint = pArcLineSeg->GetStartPoint();
				EndPoint = pArcLineSeg->GetEndPoint();

				if(j == 0)
				{
					if(nBeforeType != -1)
					{
						PointArray.push_back(EndPoint);
					}
					else
					{
						PointArray.push_back(StartPoint);
						PointArray.push_back(EndPoint);
					}
				}
				else
				{
					PointArray.push_back(EndPoint);
				}
			}
			delete pLineSegArray;
		}

		nBeforeType = pDrawObj->GetType();

		if(nDrawObjCount == 1)
			nFirstObjectType = nBeforeType;
	}

	if(PointArray.size()!= 0)
	{
		if(nFirstObjectType == CPGEOM_LINE)
		{
			PointArray[PointArray.size() - 1] =  PointArray[0];
		}
		else // ARC
		{
			if(nBeforeType == CPGEOM_ARC)
			{
				PointArray[PointArray.size() - 1] =  PointArray[0];
			}
			else//Current°¡ Line and First°¡ ARC
			{
				PointArray[0] =  PointArray[PointArray.size() - 1];
			}
		}
	}
}

gf_polygon * Polygon::Make_GFPolygon(double dDivision /*= 10*/, int UnitType /*= 0*/)
{
	int			i = 0, j = 0;
	int			nCount = 0;
	Polygon	*pCPCPolygon = NULL;
	gf_polygon	*pGFPolygon = NULL;
	PointD		pt;

	PointDArray	ptPointArray;ptPointArray.resize(0);

	ConvertToPointArray(ptPointArray, m_pBoundary, dDivision, UnitType);

	// The last point from ConvertToPointArray is always the same as the first point.
	int nPoints = ptPointArray.size() - 1;

	if (nPoints < 3) return NULL;	// A contour must have at least 3 points

	pGFPolygon = gfAllocPolygon( 1 + m_pCuttingArea->size());
	gfAllocContour(pGFPolygon, nCount, nPoints);

	pGFPolygon->hole[nCount] = FALSE;

	for (i = 0; i < nPoints; i++)
	{
		pt = ptPointArray.at(i);
		pGFPolygon->contour[nCount].point[i].x = pt.GetX();
		pGFPolygon->contour[nCount].point[i].y = pt.GetY();
	}
	pGFPolygon->box[nCount] = gfGetBoundingBox(pGFPolygon->contour[nCount]);
	nCount++;

	for(i = 0; i < m_pCuttingArea->size(); i++)
	{
		pCPCPolygon = m_pCuttingArea->at(i);
		ptPointArray.resize(0);
		ConvertToPointArray(ptPointArray, pCPCPolygon->m_pBoundary, dDivision, UnitType);

		nPoints = ptPointArray.size() - 1;
		if (nPoints < 3) continue;	// A contour must have at least 3 points

		gfAllocContour(pGFPolygon, nCount, nPoints);

		pGFPolygon->hole[nCount] = TRUE;

		for (j = 0; j < nPoints; j++)
		{
			pt = ptPointArray.at(j);
			pGFPolygon->contour[nCount].point[j].x = pt.GetX();
			pGFPolygon->contour[nCount].point[j].y = pt.GetY();
		}
		pGFPolygon->box[nCount] = gfGetBoundingBox(pGFPolygon->contour[nCount]);
		nCount++;
	}
	pGFPolygon->nContours = nCount;

	return pGFPolygon;
}