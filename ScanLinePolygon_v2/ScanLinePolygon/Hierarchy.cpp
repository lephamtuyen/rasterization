//===========================================
//
// Functions to order contours hierarchically
//
//===========================================
//
// Main Function:
//		gfContourHierarchy - Hierarchically orders all contours in a polygon
// -----------------------------------------------------------------------------

#include <stdlib.h>
#include "gf_types.h"
#include "gfunc.h"

//==================
// Private Functions
//==================

void vRegisterContour(int nContour, gf_polygon *pPoly, gf_polygon *pResult, int &nLatestIndex)
{
	// Register the specified contour in the hierarchical polygon family
	if (pPoly->contour[nContour].nPoints > 0)
	{
		// Copy the box, hole, and contour data if the contour was not registered
		pResult->box[nLatestIndex].ll.x = pPoly->box[nContour].ll.x;
		pResult->box[nLatestIndex].ll.y = pPoly->box[nContour].ll.y;
		pResult->box[nLatestIndex].ur.x = pPoly->box[nContour].ur.x;
		pResult->box[nLatestIndex].ur.y = pPoly->box[nContour].ur.y;
		pResult->hole[nLatestIndex] = pPoly->hole[nContour];

		gfAllocContour(pResult, nLatestIndex, pPoly->contour[nContour].nPoints);
		for (int nID = 0; nID < pResult->contour[nLatestIndex].nPoints; nID++)
		{
			pResult->contour[nLatestIndex].point[nID].x = pPoly->contour[nContour].point[nID].x;
			pResult->contour[nLatestIndex].point[nID].y = pPoly->contour[nContour].point[nID].y;
		}

		// Increase the hierarchical contour index for the next contour
		nLatestIndex++;

		// Set the number of point to negative to indicate the contour is registered
		pPoly->contour[nContour].nPoints = -ABS(pPoly->contour[nContour].nPoints);
	}
}

void vBuildDescendants(gf_tree *pRoot, gf_polygon *pPoly, gf_polygon *pResult, int &nLatestIndex)
{
	// Build descendant hierarchy of the trunk "pRoot"
	gf_tree *pTemp;

	// Register the root contour of the current trunk
	vRegisterContour(pRoot->id, pPoly, pResult, nLatestIndex);

	// Register all eldest descendants of the current trunk
	while(pRoot->child)
	{
		pRoot = pRoot->child;
		vRegisterContour(pRoot->id, pPoly, pResult, nLatestIndex);
	}

	// Register all siblings in each descent family
	while(pRoot->parent)
	{
		pTemp = pRoot;
		while(pTemp->next && pTemp->next->parent)
		{
			if(pRoot->parent->id == pTemp->next->parent->id)
			{
				pTemp = pTemp->next;
				vRegisterContour(pTemp->id, pPoly, pResult, nLatestIndex);

				while(pTemp->child)
				{
					pRoot = pTemp = pTemp->child;
					vRegisterContour(pTemp->id, pPoly, pResult, nLatestIndex);
				}
			}
			else
				pTemp = pTemp->next;
		}
		pRoot = pRoot->parent;
	}
}

void vBuildHierarchy(gf_tree *p1stHole, gf_polygon *pPoly, gf_polygon *pResult, int &nLatestIndex)
{
	gf_tree *pHole, *pRoot;

	pHole = p1stHole;

	// Loop through all holes
	while(pHole)
	{
		// Trace back to get a root contour of the current hole
		pRoot = pHole;
		while(pRoot->parent)
		{
			pRoot = pRoot->parent;
		}

		// Build descendant hierarchy of the current trunk
		vBuildDescendants(pRoot, pPoly, pResult, nLatestIndex);
		pHole = pHole->next;
	}
}

gf_tree *pGetTree(int nID, gf_tree *pStart)
{
	// Returns a tree whose ID is 'nID'

	gf_tree *pTree, *pCheck;

	pCheck = pStart;
	while(pCheck)
	{
		pTree = pCheck;
		while(pTree)
		{
			if (pTree->id == nID) return pTree;
			pTree = pTree->parent;
		}
		pCheck = pCheck->next;
	}
	return NULL;
}

gf_polygon *pContourBelongsTo(gf_polygon *pContour, gf_polygon *pPoly, bool bGetHole, int *nID)
{
	// Get a contour to which the 'pContour' belongs
	// Returns a hole contour if hGetHole is true.
	// Returns a dark contour if bGetHole is false.
	// Also outputs ID (nID) of the return contour(dark or hole)

	*nID = -1;

	gf_polygon *pReturn = NULL;
	gf_polygon *pCheck = NULL;

	int nContour;
	for (nContour = 0; nContour < pPoly->nContours; nContour++)
	{
		if (bGetHole)
		{
			// skip if requested to get a hole but the contour is not a hole
			if (!pPoly->hole[nContour]) continue;
		}
		else
		{
			// skip if requested to get a dark contour but the contour is a hole
			if (pPoly->hole[nContour]) continue;
		}

		// Construct a temporary polygon 
		if (pReturn == pCheck)
		{
			pCheck = gfAllocPolygon(1);
			pCheck->hole[0] = false;
		}
		pCheck->contour[0] = pPoly->contour[nContour];

		// Search for the smallest polygon that includes the 'pContour'
		if (bContourInContour(pContour, pCheck))
		{
			if (pReturn)
			{
				if (bContourInContour(pCheck, pReturn))
				{
					FREE(pReturn);
//					gfFreePolygon(pReturn);
					pReturn = pCheck;
					*nID = nContour;
				}
			}
			else
			{
				pReturn = pCheck;
				*nID = nContour;
			}
		}
	}

	if (*nID >= 0)
	{
		return pReturn;
	}
	else
	{
		gfFreeContours(pCheck);
		if (pReturn) FREE(pReturn);
//		if (pReturn) gfFreePolygon(pReturn);
		return NULL;
	}
}

gf_tree *pHoleHierarchy(gf_polygon *pPoly)
{
	// Construct hierarchies from independent holes

	gf_tree *p1stHole = NULL;
	gf_tree *pTreeContour = NULL;

	gf_polygon *pTempDark, *pTempHole;
	gf_tree    *pTree, *pTreeHole;

	int nContour, nID;

	for (nContour = 0; nContour < pPoly->nContours; nContour++)
	{
		if (!pPoly->hole[nContour]) continue;			// Not a hole
		if (pGetTree(nContour, p1stHole)) continue;		// Already processed

		gf_polygon *pContour = gfAllocPolygon(1);
		pContour->hole[0] = false;
		gfAllocContour(pContour, 0, pPoly->contour[nContour].nPoints);
		for (nID = 0; nID < pContour->contour[0].nPoints; nID++)
		{
			pContour->contour[0].point[nID].x = pPoly->contour[nContour].point[nID].x;
			pContour->contour[0].point[nID].y = pPoly->contour[nContour].point[nID].y;
		}

		// Created a new trunk
		pTreeHole = (gf_tree *) malloc(sizeof(gf_tree));
		pTreeHole->id     = nContour;
		pTreeHole->parent = NULL;
		pTreeHole->child  = NULL;
		pTreeHole->prev   = NULL;
		pTreeHole->next   = NULL;

		if (!p1stHole) p1stHole = pTreeHole;

		if (pTreeContour)
		{
			pTreeHole->prev = pTreeContour;
			pTreeContour->next = pTreeHole;
		}
		pTreeContour = pTreeHole;
		
		pTempHole = pContour;

		// Get a dark(solid) contour to which the hole 'pTempHole' belongs
		while ((pTempDark = pContourBelongsTo(pTempHole, pPoly, false, &nID)) != NULL)
		{
			pTree = pTreeHole;
			pTreeHole = (gf_tree *) malloc(sizeof(gf_tree));
			pTreeHole->id     = nID;
			pTreeHole->parent = NULL;
			pTreeHole->child  = pTree;
			pTreeHole->prev   = NULL;
			pTreeHole->next   = NULL;
			pTree->parent = pTreeHole;
			
			// Get a hole contour to which the dark contour 'pTempDark' belongs
			if (!(pTempHole = pContourBelongsTo(pTempDark, pPoly, true, &nID))) break;

			if ((pTree = pGetTree(nID, p1stHole)) != NULL)
			{
				while (pTree->next)
				{
					pTree = pTree->next;
				}
				pTree->next = pTreeHole;
				break;
			}
			pTree = pTreeHole;
			pTreeHole = (gf_tree *) malloc(sizeof(gf_tree));
			pTreeHole->id     = nID;
			pTreeHole->parent = NULL;
			pTreeHole->child  = pTree;
			pTreeHole->prev   = NULL;
			pTreeHole->next   = NULL;
			pTree->parent = pTreeHole;
		}
		
		gfFreePolygon(pContour);
	}

	return p1stHole;
}

//=================
// Global Functions
//=================

gf_polygon *gfContourHierarchy(gf_polygon *pPoly)
{
	// Construct a hierarchy of all contours in a polygon

	gf_polygon	*pResult;
	gf_tree		*p1stHole;
	int         nContour;

	// Construct a hole hierarchy
	p1stHole = pHoleHierarchy(pPoly);

	// Initialize the result polygon
	pResult = gfAllocPolygon(pPoly->nContours);

	// Initialize the contour index
	int nLatestIndex = 0;

	// Construct a contour hierarchy
	vBuildHierarchy(p1stHole, pPoly, pResult, nLatestIndex);

	// Register the rest of contours that do not belong to any trunk
	for (nContour = 0; nContour < pPoly->nContours; nContour++)
	{
		vRegisterContour(nContour, pPoly, pResult, nLatestIndex);
		pPoly->contour[nContour].nPoints = ABS(pPoly->contour[nContour].nPoints);
	}

	// Free memory allocated for the hole list
	gf_tree *pTemp, *pNext;
	pTemp = p1stHole;
	while (pTemp)
	{
		pNext = pTemp->next;
		FREE(pTemp);
		pTemp = pNext;
	}

	gfFreePolygon(pPoly);

	return pResult;
}
