// -----------------------------------------------------------------------------
// Sorting Functions
//
// Functions:
//		pQuickSort - Sort an array 
//                   then return sorted array and their original indices.
// -----------------------------------------------------------------------------
#include "gfunc.h"
#include <stdlib.h>

#define NISTACK 256

void vQuickSort(bool bAscend,		// ascending if true, descending if false
				int begin, int end,	// Beginning and ending array indices
				double *A,			// array to be sorted
				int *idx)			// indices of the sorted array
{
	// The array A is both input and output so that the result is output to A.
	// This function returns the original array index to an integer array.

	double X, TEMP;

	int i, j, iorder, is, l1, l2, m, itemp;
	int ist1[NISTACK], ist2[NISTACK];
	
	// initialize stacks
	for (i = 0; i < NISTACK; i++)
	{
		ist1[i] = 0;
		ist2[i] = 0;
	}

	// initialize the original array index
	for (i = 0; i < (end - begin + 1); i++)
	{
		idx[i] = i;
	}

	// set sorting order
	if (bAscend)
		iorder = 1;
	else
		iorder = -1;

	// initialize parameters
	is = 0;
	ist1[0] = begin;
	ist2[0] = end;

	// take top request from stack
	while (is >= 0)
	{
		l1 = ist1[is];
		l2 = ist2[is];
		is--;

		while (l1 < l2)
		{
			i = l1;
			j = l2;
			m = (l1 + l2) / 2;
			X = iorder * A[m];

			while(i <= j)
			{
				while(iorder * A[i] < X) i++;
				while(X < iorder * A[j]) j--;

				// swap
				if (i <= j)
				{
					itemp = idx[i - begin];
					idx[i - begin] = idx[j - begin];
					idx[j - begin] = itemp;
					TEMP = A[i];
					A[i] = A[j];
					A[j] = TEMP;
					i++;
					j--;
				}
			}

			// stack request to sort right partition
			if (i < l2)
			{
				is++;
				ist1[is] = i;
				ist2[is] = l2;
			}
			l2 = j;
		}
	}
	return;
}

void vQuickSort(bool bAscend,		// ascending if true, descending if false
				int begin, int end,	// Beginning and ending array indices
				int *A,				// array to be sorted
				int *idx)			// indices of the sorted array
{
	// The array A is both input and output so that the result is output to A.
	// This function returns the original array index to an integer array.
	
	int X, TEMP;
	
	int i, j, iorder, is, l1, l2, m, itemp;
	int ist1[NISTACK], ist2[NISTACK];
	
	// initialize stacks
	for (i = 0; i < NISTACK; i++)
	{
		ist1[i] = 0;
		ist2[i] = 0;
	}
	
	// initialize the original array index
	for (i = 0; i < (end - begin + 1); i++)
	{
		idx[i] = i;
	}
	
	// set sorting order
	if (bAscend)
		iorder = 1;
	else
		iorder = -1;
	
	// initialize parameters
	is = 0;
	ist1[0] = begin;
	ist2[0] = end;
	
	// take top request from stack
	while (is >= 0)
	{
		l1 = ist1[is];
		l2 = ist2[is];
		is--;
		
		while (l1 < l2)
		{
			i = l1;
			j = l2;
			m = (l1 + l2) / 2;
			X = iorder * A[m];
			
			while(i <= j)
			{
				while(iorder * A[i] < X) i++;
				while(X < iorder * A[j]) j--;
				
				// swap
				if (i <= j)
				{
					itemp = idx[i - begin];
					idx[i - begin] = idx[j - begin];
					idx[j - begin] = itemp;
					TEMP = A[i];
					A[i] = A[j];
					A[j] = TEMP;
					i++;
					j--;
				}
			}
			
			// stack request to sort right partition
			if (i < l2)
			{
				is++;
				ist1[is] = i;
				ist2[is] = l2;
			}
			l2 = j;
		}
	}
	return;
}

void vQuickSort(bool bAscend,		// ascending if true, descending if false
				int begin, int end,	// Beginning and ending array indices
				double *A)			// array to be sorted
{
	// The array A is both input and output so that the result is output to A.
	// This function returns the original array index to an integer array.

	double X, TEMP;

	int i, j, iorder, is, l1, l2, m;
	int ist1[NISTACK], ist2[NISTACK];
	
	// initialize stacks
	for (i = 0; i < NISTACK; i++)
	{
		ist1[i] = 0;
		ist2[i] = 0;
	}

	// set sorting order
	if (bAscend)
		iorder = 1;
	else
		iorder = -1;

	// initialize parameters
	is = 0;
	ist1[0] = begin;
	ist2[0] = end;

	// take top request from stack
	while (is >= 0)
	{
		l1 = ist1[is];
		l2 = ist2[is];
		is--;

		while (l1 < l2)
		{
			i = l1;
			j = l2;
			m = (l1 + l2) / 2;
			X = iorder * A[m];

			while(i <= j)
			{
				while(iorder * A[i] < X) i++;
				while(X < iorder * A[j]) j--;

				// swap
				if (i <= j)
				{
					TEMP = A[i];
					A[i] = A[j];
					A[j] = TEMP;
					i++;
					j--;
				}
			}

			// stack request to sort right partition
			if (i < l2)
			{
				is++;
				ist1[is] = i;
				ist2[is] = l2;
			}
			l2 = j;
		}
	}
	return;
}

void vQuickSort(bool bAscend,		// ascending if true, descending if false
				int begin, int end,	// Beginning and ending array indices
				int *A)				// array to be sorted
{
	// The array A is both input and output so that the result is output to A.
	// This function returns the original array index to an integer array.
	
	int X, TEMP;
	
	int i, j, iorder, is, l1, l2, m;
	int ist1[NISTACK], ist2[NISTACK];
	
	// initialize stacks
	for (i = 0; i < NISTACK; i++)
	{
		ist1[i] = 0;
		ist2[i] = 0;
	}
	
	// set sorting order
	if (bAscend)
		iorder = 1;
	else
		iorder = -1;
	
	// initialize parameters
	is = 0;
	ist1[0] = begin;
	ist2[0] = end;
	
	// take top request from stack
	while (is >= 0)
	{
		l1 = ist1[is];
		l2 = ist2[is];
		is--;
		
		while (l1 < l2)
		{
			i = l1;
			j = l2;
			m = (l1 + l2) / 2;
			X = iorder * A[m];
			
			while(i <= j)
			{
				while(iorder * A[i] < X) i++;
				while(X < iorder * A[j]) j--;
				
				// swap
				if (i <= j)
				{
					TEMP = A[i];
					A[i] = A[j];
					A[j] = TEMP;
					i++;
					j--;
				}
			}
			
			// stack request to sort right partition
			if (i < l2)
			{
				is++;
				ist1[is] = i;
				ist2[is] = l2;
			}
			l2 = j;
		}
	}
	return;
}
