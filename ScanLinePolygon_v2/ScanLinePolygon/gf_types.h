#ifndef _GF_TYPES_H_
#define _GF_TYPES_H_

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

#if defined(__cplusplus) || \
	defined(bool) /* for C compilation with C99 bool (macro) */
typedef bool   bool_t;
#else
typedef BOOL   bool_t;
#endif /* __cplusplus */

#ifndef ONEPI
#define ONEPI 3.141592653589793238462
#endif

#ifndef AMIN
#define AMIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef AMAX
#define AMAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(a) ((a) > 0 ? (a) : -(a))
#endif

#ifndef NINT
#define NINT(a) ((a) > 0 ? ((int) ((a) + 0.5)) : ((int) ((a) - 0.5)))
#endif

#ifndef FREE
#define FREE(p) {if (p) {free(p); (p)= NULL;}}
#endif

typedef struct					// 3-D Point (double)
{
	double		x, y, z;
} gf3_point;

typedef struct					// Point (double)
{
	double		x, y;
} gf_point;

typedef struct					// Point (float)
{
	float		x, y;
} gf_point_float;

typedef struct					// Line
{
	gf_point	s, e;			// Line starting and ending points
} gf_line;

typedef struct					// Box
{
	gf_point	ll, ur;			// Lower-left and upper-right corners
} gf_box;

typedef struct					// Circle
{
	gf_point	c;				// Center point location
	double		r;				// Radius
} gf_circle;

typedef struct					// Point list for poly line or contour
{
	int			nPoints;		// Number of points in the list
	gf_point	*point;			// Pointer to the point array
} gf_point_list, gf_contour;

typedef struct					// Polygon (a set of contours: consists of zero or more contours)
{
	int			nContours;		// Number of contours belong to the polygon
	int			nAllocated;		// Allocated count (allocated >= nCoutours). Needed to allocate a chunk at a time.
	bool		*hole;			// Hole indicator of each contour (true if hole, false if dark)
	gf_box		*box;			// Pointer to the bounding box array
	gf_contour	*contour;		// Pointer to the contour array
} gf_polygon;

typedef struct _gf_tree			// Tree hierarchy
{
	int               id;		// ID number of itself
	struct _gf_tree  *parent;	// pointer to the parent
	struct _gf_tree  *child;	// pointer to the eldest child
	struct _gf_tree  *prev;		// pointer to the elder sibling
	struct _gf_tree  *next;		// pointer to the younger sibling
} gf_tree;

#endif
