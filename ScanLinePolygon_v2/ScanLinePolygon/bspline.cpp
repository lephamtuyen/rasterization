//===========================
//
// B-spline Functions
//
//===========================

#include "gf_types.h"
#include "gfunc.h"
#include <float.h>


double *basis(int p, double t, int nv, double *u);
double *pGetBsplineBasis(int p, double t, int nv, double *u);
double *pGetBsplineKnotsUniform(int nv, int p);
double *pGetBsplineKnotsPredefined(int& nu, int nw, double *w, int *knotMultiples);

gf3_point *pBspline(int nv, int p, int nc, gf3_point *v)
{
	//  Function to generate an open B-spline curve using an uniform open knot vector
	// 
	// 	INPUT:
	// 	   nv  = number of control points
	// 	   p   = order of the basis function
	// 	   nc  = number of points to be calculated on the curve
	// 	   v   = list of control points
	// 	RETURN:
	// 	   List of points on the B-spline curve

	int m = nv + p;	// last index of the knot vector array: m = n + p + 1 = nv + p (number of knot values = m + 1)
	double *u;		// knot vector
	double *N;		// basis for a single value of t
	double t;		// parameter value 0 <= t <= 1

	int i, k;

	double step;

	// generate the uniform open knot vector

	u = pGetBsplineKnotsUniform(nv, p);

	// calculate the points on the bspline curve

	gf3_point *c = (gf3_point *) malloc(nc * sizeof(gf3_point));

	t = 0;
	step = u[m] / (nc - 1);

	for (k = 0; k < nc; k++)
	{
		if (u[m] - t < 5.0e-6)
			t = u[m];

		N = pGetBsplineBasis(p, t, nv, u);		// generate the basis function for this value of t

		// generate a point on the curve for each axis (x, y, z)

		c[k].x = 0.;
		c[k].y = 0.;
		c[k].z = 0.;

		for (i = 0; i < nv; i++)				// Do local matrix multiplication
		{
			c[k].x += N[i] * v[i].x;
			c[k].y += N[i] * v[i].y;
			c[k].z += N[i] * v[i].z;
		}
		FREE(N);

		t += step;
	}
	FREE(u);

	return c;
}

gf3_point *pBsplineWithKnots(int nv, int p, int nw, int nc,
							 gf3_point *v, double *w, int *knotMultiples)
{
	//  Function to generate an open B-spline curve with a predefined knot vector
	// 
	// 	INPUT:
	// 	   nv  = number of control points
	// 	   p   = order of the basis function
	//	   nw  = number of knots
	// 	   nc  = number of points to be calculated on the curve
	// 	   v[] = list of control points
	//	   w[] = list of knots
	//	   knotMultiples[] = list of knot multiples
	// 	RETURN:
	// 	   List of points on the B-spline curve
	
	int nu;			// total number of knots
	int m;			// last index of the knot vector array ( = nu - 1 )
	double *u;		// knot vector
	double *N;		// basis for a single value of t
	double t;		// parameter value 0 <= t <= 1
	
	int i, k;
	
	double step;
	
	// convert multiple knots to single knot
	
	u = pGetBsplineKnotsPredefined(nu, nw, w, knotMultiples);
	m = nu - 1;
	
	// calculate the points on the bspline curve
	
	gf3_point *c = (gf3_point *) malloc(nc * sizeof(gf3_point));
	
	t = 0;
	step = u[m] / (nc - 1);
	
	for (k = 0; k < nc; k++)
	{
		if (u[m] - t < 5.0e-6)
			t = u[m];
		
		N = pGetBsplineBasis(p, t, nv, u);		// generate the basis function for this value of t
		
		// generate a point on the curve for each axis (x, y, z)
		
		c[k].x = 0.;
		c[k].y = 0.;
		c[k].z = 0.;
		
		for (i = 0; i < nv; i++)				// Do local matrix multiplication
		{
			c[k].x += N[i] * v[i].x;
			c[k].y += N[i] * v[i].y;
			c[k].z += N[i] * v[i].z;
		}
		FREE(N);
		
		t += step;
	}
	FREE(u);
	
	return c;
}

gf3_point *pBsplineWithKnotsClosed(int nv, int p, int nw, int nc,
							 gf3_point *v, double *w, int *knotMultiples)
{
	//  Function to generate a closed B-spline curve with a predefined knot vector
	// 
	// 	INPUT:
	// 	   nv  = number of control points
	// 	   p   = order of the basis function
	//	   nw  = number of knots
	// 	   nc  = number of points to be calculated on the curve
	// 	   v[] = list of control points
	//	   w[] = list of knots
	//	   knotMultiples[] = list of knot multiples
	// 	RETURN:
	// 	   List of points on the B-spline curve

	int nvClosed = nv + 1;
	gf3_point *vClosed = (gf3_point *) malloc(nvClosed * sizeof(gf3_point));

	int i;
	for (i = 0; i < nv; i++)
	{
		vClosed[i].x = v[i].x;
		vClosed[i].y = v[i].y;
		vClosed[i].z = v[i].z;
	}
	vClosed[i].x = v[0].x;
	vClosed[i].y = v[0].y;
	vClosed[i].z = v[0].z;

	if (v[nv-1].x == v[0].x && v[nv-1].y == v[0].y && v[nv-1].z == v[0].z)
		nvClosed--;

	// convert multiple knots to single knot
	
	int nu;
	double *u;
	u = pGetBsplineKnotsPredefined(nu, nw, w, knotMultiples);

	int nuClosed = nu + p + 2;
	double *uClosed = (double *) malloc(nuClosed * sizeof(double));
	
	int *multipliersClosed = (int *) malloc(nuClosed * sizeof(int));
	int k = 0;
	for (i = 0, k = 0; i < nu; i++, k++)
	{
		uClosed[k] = u[i];
		multipliersClosed[k] = 1;
	}
	for (i = 0; k < nuClosed; i++, k++)
	{
		uClosed[k] = u[i];
		multipliersClosed[k] = 1;
	}
	FREE(u);

	gf3_point *c = pBsplineWithKnots(nvClosed, p, nuClosed, nc, vClosed, uClosed, multipliersClosed);

	FREE(vClosed);
	FREE(uClosed);
	FREE(multipliersClosed);

	return c;
}

double *pGetBsplineKnotsUniform(int nv, int p)
{
	// 	Function to generate a B-spline open knot vector with multiplicity
	// 	equal to the (order + 1) at the ends in order to clamp at ends.
	// 
	// 	p	= order of the basis function
	// 	nv	= number of control points
	// 	u[]	= array containing the knot vector
	
	int m;			// last index of the knot vector array (number of know values = m + 1)
	int nClamp;
	int nMiddle;
	int i, j;

	m = nv + p;
	nClamp = p + 1;
	nMiddle = (m + 1) - 2 * nClamp;

	double *u = (double *) calloc(m + 1, sizeof(double));

	for (i = 0; i < nClamp; i++)
	{
		u[i] = 0;
	}
	for (j = 0; j <= nMiddle; j++, i++)
	{
		u[i] = u[i-1] + 1;
	}
	for (; i <= m; i++)
	{
		u[i] = u[i-1];
	}

	return u;
}

double *pGetBsplineKnotsPredefined(int& nu, int nw, double *w, int *knotMultiples)
{
	// 	Function to generate a B-spline open knot vector with multiplicity
	// 	equal to the (order + 1) at the ends in order to clamp at ends.
	// 
	// 	w	= order of the basis function
	// 	nv	= number of control points
	// 	u[]	= array containing the knot vector

	int i, j;

	nu = 0;
	for (i = 0; i < nw; i++)
	{
		nu += knotMultiples[i];
	}

	double *u = (double *) malloc(nu * sizeof(double));

	int inu = 0;
	for (i = 0; i < nw; i++)
	{
		for (j = 0; j < knotMultiples[i]; j++)
		{
			u[inu++] = w[i];
		}
	}

	return u;
}

double *pGetBsplineBasis(int p, double t, int nv, double *u)
{
	// 	Function to generate B-spline basis for open knot vectors
	// 
	// 	p        = order of the B-spline basis function
	// 	t        = parameter value
	// 	nv       = number of control points
	// 	u[]      = knot vector

	int i;
	int m = nv + p;		// last index of the knot vector array (number of know values = m + 1)

	double *N;			// Basis function
	N = (double *) calloc(nv, sizeof(double));

	// rule out special cases

	if (t == u[0])		
	{
		N[0] = 1.0;
		return N;
	}
	else if (t == u[m])
	{
		N[nv-1] = 1.0;
		return N;
	}

	double *Ntemp;
	Ntemp = basis(p, t, nv, u);

	for (i = 0; i < nv; i++)
	{
		N[i] = Ntemp[i];
	}
	for(int i=0; i< nv; i++)
	{
		double dValue = N[i];
		if(_isnan(dValue)!=0)
		{
			// isnan ... 
			N[i] = 0.0;
		}
	}
	FREE(Ntemp);

	return N;
}

double *basis(int p, double t, int nv, double *u)
{
	// 	Function to generate B-spline basis for open knot vectors
	//	This is a private function to supply results to pGetBsplineBasis and pGetBsplineBasisAll

	int i, k, d;
	int m = nv + p;			// last index of the knot vector array (number of know values = m + 1)

	double *N;				// Basis function
	N = (double *) calloc(m, sizeof(double));

	// calculate the basis functions of degree 0

	for (k = 0; k < m; k++)
	{
		if ( t >= u[k] && t < u[k+1] )	// t is within the range u[i] <= t < u[i+1]
		{
			N[k] = 1;
			break;
		}
	}

	// chahn : crash free(N)!
	if(k == m)
		k = k-1;
	//

	// calculate the basis functions for degree 1 and higher

	for (d = 1; d <= p; d++)
	{
		N[k-d] = N[k-d+1] * (u[k+1] - t) / (u[k+1] - u[k-d+1]);		// right (south-west corner) term only

		for (i = k-d+1; i < k; i++)									// compute internal terms
		{
			double a = N[i] * (t - u[i]) / (u[i+d] - u[i]);
			double b = N[i+1] * (u[i+d+1] - t) / (u[i+d+1] - u[i+1]);
			N[i] = a + b;
		}

		N[k] = N[k] * (t - u[k]) / (u[k+d] - u[k]);					// left (north-west corner) term only
	}

	return N;
}

gf3_point *pBsplineWithKnotsDeBoor(int nv, int p, int nu, int nc,
							 gf3_point *v, double *u, int *knotMultiples)
{
	//  Function to generate an open B-spline curve with a predefined knot vector
	// 
	// 	INPUT:
	// 	   nv  = number of control points
	// 	   p   = order of the basis function
	//	   nu  = number of knots
	// 	   nc  = number of points to be calculated on the curve
	// 	   v[] = list of control points
	//	   u[] = list of knots
	//	   knotMultiples[] = list of knot multiples
	// 	RETURN:
	// 	   List of points on the B-spline curve
	
	int m = nu - 1;	// last index of the knot vector array
	double t;		// parameter value 0 <= t <= 1
	
	int h;
	int s;

	int i, j, k, r;
	
	double step;
	
	// calculate the points on the bspline curve
	
	gf3_point *c = (gf3_point *) malloc(nc * sizeof(gf3_point));

	gf3_point *P = (gf3_point *) malloc(nv * sizeof(gf3_point));
	gf3_point *T = (gf3_point *) malloc(nv * sizeof(gf3_point));
	
	step = u[m] / (nc - 1);
	t = 0;
	
	for (j = 0; j < nc; j++)
	{

		if (u[m] - t < 5.0e-6)
			t = u[m];
		
		for (k = 0; k < m; k++)
		{
			if ( t >= u[k] && t < u[k+1] )	// t is within the range u[i] <= t < u[i+1]
				break;
		}

		if (t == u[k])
		{
			if (knotMultiples[k] > p)
				s = p;
			else
				s = knotMultiples[k];
			h = p - s;
		}
		else
		{
			s = 0;
			h = p;
		}

		for (i = k-p; i <= k-s; i++)
		{
			P[i].x = v[i].x;
			P[i].y = v[i].y;
			P[i].z = v[i].z;
		}
		
		for (r = 1; r <= h; r++)
		{
			for (i = k-p+r; i <= k-s; i++)
			{
				double a = (t - u[i]) / (u[i+p-r+1] - u[i]);
				T[i].x = (1 - a) * P[i-1].x + a * P[i].x;
				T[i].y = (1 - a) * P[i-1].y + a * P[i].y;
				T[i].z = (1 - a) * P[i-1].z + a * P[i].z;
			}
			for (i = k-p+r; i <= k-s; i++)
			{
				P[i].x = T[i].x;
				P[i].y = T[i].y;
				P[i].z = T[i].z;
			}
		}

		c[j].x = P[k-s].x;
		c[j].y = P[k-s].y;
		c[j].z = P[k-s].z;
	
		t += step;
	}
	FREE(P);
	FREE(T);
	
	return c;
}

gf3_point *pBsplineWithKnotsDeBoorClosed(int nv, int p, int nu, int nc,
							 gf3_point *v, double *u, int *knotMultiples)
{
	//  Function to generate a closed B-spline curve with a predefined knot vector
	// 
	// 	INPUT:
	// 	   nv  = number of control points
	// 	   p   = order of the basis function
	//	   nu  = number of knots
	// 	   nc  = number of points to be calculated on the curve
	// 	   v[] = list of control points
	//	   u[] = list of knots
	//	   knotMultiples[] = list of knot multiples
	// 	RETURN:
	// 	   List of points on the B-spline curve
	
	int m = nu - 1;	// last index of the knot vector array
	double t;		// parameter value 0 <= t <= 1
	
	int h;
	int s;
	
	int i, j, k, r;
	
	double step;
	
	// calculate the points on the bspline curve
	
	gf3_point *c = (gf3_point *) malloc(nc * sizeof(gf3_point));
	
	gf3_point *P = (gf3_point *) malloc(nv * sizeof(gf3_point));
	gf3_point *T = (gf3_point *) malloc(nv * sizeof(gf3_point));
	
	step = u[m] / (nc - 1);
	t = 0;
	
	for (j = 0; j < nc; j++)
	{
		
		if (u[m] - t < 5.0e-6)
			t = u[m];
		
		for (k = 0; k < m; k++)
		{
			if ( t >= u[k] && t < u[k+1] )	// t is within the range u[i] <= t < u[i+1]
				break;
		}
		
		if (t == u[k])
		{
			if (knotMultiples[k] > p)
				s = p;
			else
				s = knotMultiples[k];
			h = p - s;
		}
		else
		{
			s = 0;
			h = p;
		}
		
		for (i = k-p; i <= k-s; i++)
		{
			P[i].x = v[i].x;
			P[i].y = v[i].y;
			P[i].z = v[i].z;
		}
		
		for (r = 1; r <= h; r++)
		{
			for (i = k-p+r; i <= k-s; i++)
			{
				double a = (t - u[i]) / (u[i+p-r+1] - u[i]);
				T[i].x = (1 - a) * P[i-1].x + a * P[i].x;
				T[i].y = (1 - a) * P[i-1].y + a * P[i].y;
				T[i].z = (1 - a) * P[i-1].z + a * P[i].z;
			}
			for (i = k-p+r; i <= k-s; i++)
			{
				P[i].x = T[i].x;
				P[i].y = T[i].y;
				P[i].z = T[i].z;
			}
		}
		
		c[j].x = P[k-s].x;
		c[j].y = P[k-s].y;
		c[j].z = P[k-s].z;
		
		t += step;
	}
	FREE(u);
	FREE(P);
	FREE(T);
	
	return c;
}

void gf3FreePoints(gf3_point *&p)
{
	FREE(p);
}

void gfFreeBsplineKnots(double *&u)
{
	FREE(u);
}
