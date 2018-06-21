//===========================
//
// Polygon boolean Operations
//
//===========================

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <float.h>
#include <math.h>
#include "gf_types.h"
#include "gfunc.h"
#include "comlib.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//==========
// Constants
//==========

#ifndef TRUE
#define FALSE              0
#define TRUE               1
#endif

#define LEFT   0
#define RIGHT  1

#define ABOVE  0
#define BELOW  1

#define CLIP   0
#define SUBJ   1

//=======
// Macros
//=======

#define GF_EQ(a, b)			(fabs((a) - (b)) <= DBL_EPSILON)

#define PREV_INDEX(i, n)	((i - 1 + n) % n)
#define NEXT_INDEX(i, n)	((i + 1    ) % n)

#define OPTIMAL(v, i, n)	((v[PREV_INDEX(i, n)].y != v[i].y) || \
							 (v[NEXT_INDEX(i, n)].y != v[i].y))

#define FWD_MIN(v, i, n)	((v[PREV_INDEX(i, n)].point.y >= v[i].point.y) \
						  && (v[NEXT_INDEX(i, n)].point.y > v[i].point.y))

#define NOT_FMAX(v, i, n)	(v[NEXT_INDEX(i, n)].point.y > v[i].point.y)

#define REV_MIN(v, i, n)	((v[PREV_INDEX(i, n)].point.y > v[i].point.y) \
						  && (v[NEXT_INDEX(i, n)].point.y >= v[i].point.y))

#define NOT_RMAX(v, i, n)	 (v[PREV_INDEX(i, n)].point.y > v[i].point.y)

#define P_EDGE(d,e,p,i,j)  {(d)= (e); \
	do {(d)= (d)->prev;} while (!(d)->outp[(p)]); \
	(i)= (d)->bot.x + (d)->dx * ((j)-(d)->bot.y);}

#define N_EDGE(d,e,p,i,j)  {(d)= (e); \
	do {(d)= (d)->next;} while (!(d)->outp[(p)]); \
	(i)= (d)->bot.x + (d)->dx * ((j)-(d)->bot.y);}


//===================
// Private Data Types
//===================

typedef enum	// Edge intersection classes
{
	NUL,		// Empty non-intersection
	EMX,		// External maximum
	ELI,		// External left intermediate
	TED,		// Top edge
	ERI,		// External right intermediate
	RED,		// Right edge
	IMM,		// Internal maximum and minimum
	IMN,		// Internal minimum
	EMN,		// External minimum
	EMM,		// External maximum and minimum
	LED,		// Left edge
	ILI,		// Internal left intermediate
	BED,		// Bottom edge
	IRI,		// Internal right intermediate
	IMX,		// Internal maximum
	FUL			// Full non-intersection
} vertex_type;

typedef enum	// Horizontal edge states
{
	NH,			// No horizontal edge
	BH,			// Bottom horizontal edge
	TH			// Top horizontal edge
} h_state;

typedef enum		// Edge bundle state
{
	UNBUNDLED,		// Isolated edge not within a bundle
	BUNDLE_HEAD,	// Bundle head node
	BUNDLE_TAIL		// Passive bundle tail node
} bundle_state;

typedef struct v_shape		// Internal vertex list datat ype
{
	double			x;		// X coordinate component
	double			y;		// Y coordinate component
	struct v_shape	*next;	// Pointer to next vertex in list
} vertex_node;

typedef struct p_shape		// Internal contour
{
	int				active;	// Active flag / vertex count
	bool			hole;	// Hole / external contour flag
	vertex_node		*v[2];	// Left and right vertex list ptrs
	struct p_shape	*next;	// Pointer to next polygon contour
	struct p_shape	*proxy;	// Pointer to actual structure used
} polygon_node;

typedef struct edge_shape				// Edge list data type
{
	gf_point			point;			// Piggy-backed contour vertex data
	gf_point			bot;			// Edge lower (x, y) coordinate
	gf_point			top;			// Edge upper (x, y) coordinate
	double				xb;				// Scanbeam bottom x coordinate
	double				xt;				// Scanbeam top x coordinate
	double				dx;				// Change in x for a unit y increase
	int					type;			// Suurce / target edge flag
	int					bundle[2][2];	// Bundle edge flags
	int					bside[2];		// Bundle left / right indicators
	bundle_state		bstate[2];		// Edge bundle state
	polygon_node		*outp[2];		// Output polygon / tristrip pointer
	struct edge_shape	*prev;			// Previous edge in the AET
	struct edge_shape	*next;			// Next edge in the AET
	struct edge_shape	*pred;			// Edge connected at the lower end
	struct edge_shape	*succ;			// Edge connected at the upper end
	struct edge_shape	*next_bound;	// Pointer to next bound in LMT
} edge_node;

typedef struct lmt_shape				// Local minima table
{
	double				y;				// Y coordinate at local minimum
	edge_node			*first_bound;	// Pointer to bound list
	struct lmt_shape	*next;			// Pointer to next local minimum
} lmt_node;

typedef struct sbt_t_shape			// Scanbeam tree
{
	double				y;			// Scanbeam node y value
	struct sbt_t_shape	*less;		// Pointer to nodes with lower y
	struct sbt_t_shape	*more;		// Pointer to nodes with higher y
	struct sbt_t_shape	*parent;	// Pointer to its parent nodes
	bool				id;			// Be used to traverse the whole tree
} sb_tree;

typedef struct it_shape			// Intersection table
{
	edge_node			*ie[2];	// Intersecting edge (bundle) pair
	gf_point			point;	// Point of intersection
	struct it_shape		*next;	// The next intersection table node
} it_node;

typedef struct st_shape		// Sorted edge table
{
	edge_node		*edge;	// Pointer to AET edge
	double			xb;		// Scanbeam bottom x coordinate
	double			xt;		// Scanbeam top x coordinate
	double			dx;		// Change in x for a unit y increase
	struct st_shape	*prev;	// Previous edge in sorted list
} st_node;


//============
// Global Data
//============

// Horizontal edge state transitions within scanbeam boundary
const h_state next_h_state[3][6]=
{
	//  ABOVE     BELOW     CROSS
	//  L   R     L   R     L   R
	{BH, TH,   TH, BH,   NH, NH}, // NH
	{NH, NH,   NH, NH,   TH, TH}, // BH
	{NH, NH,   NH, NH,   BH, BH}  // TH
};

bool bHoleExist;

//==================
// Private Functions
//==================

void *MALLOC(int size, char *str)
{
	void *NewAddress;

	if (size > 0)
	{
		NewAddress = (void*) malloc(size);
// 		if (!NewAddress)
// 		{
// 			exit(0);
// 		}
	}
	else
	{
		NewAddress = NULL;
	}
	return NewAddress;
}

void FreePolygon(gf_polygon *p)
{
	if (p)
	{
		int c;
		for (c = 0; c < p->nContours; c++)
		{
			if (p->contour[c].nPoints > 0) FREE(p->contour[c].point);
		}
		FREE(p->box);
		FREE(p->hole);
		FREE(p->contour);
		p->nContours= 0;
	}
}

static void reset_it(it_node **it)
{
	it_node *itn;

	while (*it)
	{
		itn = (*it)->next;
		FREE(*it);
		*it = itn;
	}
}


static void reset_lmt(lmt_node **lmt)
{
	lmt_node *lmtn;

	while (*lmt)
	{
		lmtn = (*lmt)->next;
		FREE(*lmt);
		*lmt = lmtn;
	}
}


static void insert_bound(edge_node **b, edge_node *e)
{
	edge_node *existing_bound;

	while (true)
	{
		if (!*b)
		{
			// Link node e to the tail of the list
			*b = e;
			return;
		}
		else
		{
			// Do primary sort on the x field
			if (e[0].bot.x < (*b)[0].bot.x)
			{
				// Insert a new node mid-list
				existing_bound = *b;
				*b = e;
				(*b)->next_bound = existing_bound;
				return;
			}
			else if (e[0].bot.x == (*b)[0].bot.x)
			{
				// Do secondary sort on the dx field
				if (e[0].dx < (*b)[0].dx)
				{
					// Insert a new node mid-list
					existing_bound = *b;
					*b= e;
					(*b)->next_bound= existing_bound;
					return;
				}
				else
				{
					b = &((*b)->next_bound);
				}
			}
			else
			{
				b = &((*b)->next_bound);
			}
		}
	}
}

static edge_node **bound_list(lmt_node **lmt, double y)
{
	lmt_node *existing_node;

	while (true)
	{
		if (!*lmt)
		{
			// Add node onto the tail end of the LMT
			*lmt = (lmt_node *) MALLOC(sizeof(lmt_node), "LMT insertion");
			(*lmt)->y = y;
			(*lmt)->first_bound = NULL;
			(*lmt)->next = NULL;
			return &((*lmt)->first_bound);
		}
		else
		{
			if (y < (*lmt)->y)
			{
				// Insert a new LMT node before the current node
				existing_node = *lmt;
				*lmt = (lmt_node *) MALLOC(sizeof(lmt_node), "LMT insertion");
				(*lmt)->y = y;
				(*lmt)->first_bound = NULL;
				(*lmt)->next = existing_node;
				return &((*lmt)->first_bound);
			}
			else if (y > (*lmt)->y)
			{
				lmt = &((*lmt)->next);
			}
			else
			{
				return &((*lmt)->first_bound);
			}
		}
	}
}


static void add_to_sbtree(int *entries, sb_tree **sbtree, double y)
{
	sb_tree *parent = NULL;

	while (true)
	{
		if (!*sbtree)
		{
			// Add a new tree node
			*sbtree = (sb_tree *) MALLOC(sizeof(sb_tree), "scanbeam tree insertion");
			(*sbtree)->y = y;
			(*sbtree)->less = NULL;
			(*sbtree)->more = NULL;
			(*sbtree)->parent = parent;
			(*sbtree)->id = false;
			(*entries)++;
			return;
		}
		else
		{
			parent = *sbtree;

			if ((*sbtree)->y > y)          
				sbtree = &((*sbtree)->less);
			else if ((*sbtree)->y < y)
				sbtree = &((*sbtree)->more);
			else
				return;
		}
	}
}

static void build_sbt(int *entries, double *sbt, sb_tree *sbtree)
{
	sb_tree *p = sbtree;

	while (p != NULL)
	{
		if (!p->id)
		{
			while (p->less!=NULL && !p->less->id) //traverse the left tree
				p = p->less;
			if (!p->id)
			{
				sbt[*entries] = p->y;
				(*entries)++;
				p->id = true;
			}
		}

		if (p->more!=NULL && !p->more->id)    //traverse the right tree
			p = p->more;
		else                            //back to the parent
			p = p->parent;
	}
}


static void free_sbtree(sb_tree **sbtree)
{
	sb_tree *p = *sbtree;
	while (p != NULL)
	{
		while (p->less != NULL) // the left tree
			p = p->less;
		if (p->more != NULL)    // the right tree
		{
			p = p->more;
			continue;
		}

		if (p->parent == NULL)
		{
			FREE(*sbtree);
			return;
		}
		else 
		{
			if (p->parent->less == p)
			{
				p = p->parent;
				FREE(p->less);
			}
			else
			{
				p = p->parent;
				FREE(p->more);
			}
		}
	}
}


static int count_optimal_vertices(gf_contour c)
{
	int result= 0, i;

	// Ignore non-contributing contours/
	if (c.nPoints > 0)
	{
		for (i= 0; i < c.nPoints; i++)
		{
			// Ignore superfluous vertices embedded in horizontal edges/
			if (OPTIMAL(c.point, i, c.nPoints))
				result++;
		}
	}
	return result;
}


static edge_node *build_lmt(lmt_node **lmt, sb_tree **sbtree,
							int *sbt_entries, gf_polygon *p, int type,
							gf_boolean op)
{
	int          c, i, min, max, num_edges, v, nPoints;
	int          total_vertices= 0, e_index=0;
	edge_node   *e, *edge_table;

	for (c = 0; c < p->nContours; c++)
		total_vertices += count_optimal_vertices(p->contour[c]);

	// Create the entire input polygon edge table in one go
	edge_table = (edge_node *) MALLOC(total_vertices * sizeof(edge_node),
		"edge table creation");

	for (c = 0; c < p->nContours; c++)
	{
		if (p->contour[c].nPoints < 0)
		{
			// Ignore the non-contributing contour and repair the vertex count
			p->contour[c].nPoints= -p->contour[c].nPoints;
		}
		else
		{
			// Perform contour optimization
			nPoints = 0;
			for (i = 0; i < p->contour[c].nPoints; i++)
				if (OPTIMAL(p->contour[c].point, i, p->contour[c].nPoints))
				{
					edge_table[nPoints].point.x = p->contour[c].point[i].x;
					edge_table[nPoints].point.y = p->contour[c].point[i].y;

					// Record vertex in the scanbeam table
					add_to_sbtree(sbt_entries, sbtree,
								  edge_table[nPoints].point.y);

					nPoints++;
				}

				// Do the contour forward pass
				for (min = 0; min < nPoints; min++)
				{
					// If a forward local minimum
					if (FWD_MIN(edge_table, min, nPoints))
					{
						// Search for the next local maximum
						num_edges = 1;
						max = NEXT_INDEX(min, nPoints);
						while (NOT_FMAX(edge_table, max, nPoints))
						{
							num_edges++;
							max = NEXT_INDEX(max, nPoints);
						}

						// Build the next edge list
						e = &edge_table[e_index];
						e_index += num_edges;
						v = min;
						e[0].bstate[BELOW] = UNBUNDLED;
						e[0].bundle[BELOW][CLIP] = FALSE;
						e[0].bundle[BELOW][SUBJ] = FALSE;
						for (i= 0; i < num_edges; i++)
						{
							e[i].xb = edge_table[v].point.x;
							e[i].bot.x = edge_table[v].point.x;
							e[i].bot.y = edge_table[v].point.y;

							v = NEXT_INDEX(v, nPoints);

							e[i].top.x = edge_table[v].point.x;
							e[i].top.y = edge_table[v].point.y;
							e[i].dx = (edge_table[v].point.x - e[i].bot.x) /
									  (e[i].top.y - e[i].bot.y);
							e[i].type= type;
							e[i].outp[ABOVE] = NULL;
							e[i].outp[BELOW] = NULL;
							e[i].next = NULL;
							e[i].prev = NULL;
							e[i].succ = ((num_edges > 1) && (i < (num_edges - 1))) ?
										&(e[i + 1]) : NULL;
							e[i].pred = ((num_edges > 1) && (i > 0)) ? &(e[i - 1]) : NULL;
							e[i].next_bound = NULL;
							e[i].bside[CLIP] = (op == GF_SUB) ? RIGHT : LEFT;
							e[i].bside[SUBJ] = LEFT;
						}
						insert_bound(bound_list(lmt, edge_table[min].point.y), e);
					}
				}

				// Do the contour reverse pass
				for (min = 0; min < nPoints; min++)
				{
					// If a reverse local minimum
					if (REV_MIN(edge_table, min, nPoints))
					{
						// Search for the previous local maximum
						num_edges = 1;
						max = PREV_INDEX(min, nPoints);
						while (NOT_RMAX(edge_table, max, nPoints))
						{
							num_edges++;
							max = PREV_INDEX(max, nPoints);
						}

						// Build the previous edge list
						e = &edge_table[e_index];
						e_index += num_edges;
						v = min;
						e[0].bstate[BELOW] = UNBUNDLED;
						e[0].bundle[BELOW][CLIP] = FALSE;
						e[0].bundle[BELOW][SUBJ] = FALSE;
						for (i= 0; i < num_edges; i++)
						{
							e[i].xb = edge_table[v].point.x;
							e[i].bot.x = edge_table[v].point.x;
							e[i].bot.y = edge_table[v].point.y;

							v = PREV_INDEX(v, nPoints);

							e[i].top.x = edge_table[v].point.x;
							e[i].top.y = edge_table[v].point.y;
							e[i].dx = (edge_table[v].point.x - e[i].bot.x) /
								(e[i].top.y - e[i].bot.y);
							e[i].type = type;
							e[i].outp[ABOVE] = NULL;
							e[i].outp[BELOW] = NULL;
							e[i].next = NULL;
							e[i].prev = NULL;
							e[i].succ = ((num_edges > 1) && (i < (num_edges - 1))) ?
								&(e[i + 1]) : NULL;
							e[i].pred = ((num_edges > 1) && (i > 0)) ? &(e[i - 1]) : NULL;
							e[i].next_bound = NULL;
							e[i].bside[CLIP] = (op == GF_SUB) ? RIGHT : LEFT;
							e[i].bside[SUBJ] = LEFT;
						}
						insert_bound(bound_list(lmt, edge_table[min].point.y), e);
					}
				}
		}
	}
	return edge_table;
}


static void add_edge_to_aet(edge_node **aet, edge_node *edge, edge_node *prev)
{
	while (true)
	{
		if (!*aet)
		{
			// Append edge onto the tail end of the AET
			*aet = edge;
			edge->prev = prev;
			edge->next = NULL;
			return;
		}
		else
		{
			// Do primary sort on the xb field
			if (edge->xb < (*aet)->xb)
			{
				// Insert edge here (before the AET edge)
				edge->prev = prev;
				edge->next = *aet;
				(*aet)->prev = edge;
				*aet= edge;
				return;
			}
			else if (edge->xb == (*aet)->xb)
			{
				// Do secondary sort on the dx field
				if (edge->dx < (*aet)->dx)
				{
					// Insert edge here (before the AET edge)
					edge->prev= prev;
					edge->next= *aet;
					(*aet)->prev= edge;
					*aet= edge;
					return;
				}
				else
				{
					// Head further into the AET
					prev = *aet;
					aet = &((*aet)->next);
				}
			}
			else
			{
				// Head further into the AET
				prev = *aet;
				aet = &((*aet)->next);
			}
		}
	}
}

static void add_intersection(it_node **it, edge_node *edge0, edge_node *edge1,
							 double x, double y)
{        
	it_node *existing_node;

	while (true)
	{
		if (!*it)
		{
			// Append a new node to the tail of the list
			*it = (it_node *) MALLOC(sizeof(it_node), "IT insertion");
			(*it)->ie[0] = edge0;
			(*it)->ie[1] = edge1;
			(*it)->point.x = x;
			(*it)->point.y = y;
			(*it)->next = NULL;
			return;
		}
		else if ((*it)->point.y > y)
		{
			// Insert a new node mid-list
			existing_node= *it;
			*it = (it_node *) MALLOC(sizeof(it_node), "IT insertion");
			(*it)->ie[0] = edge0;
			(*it)->ie[1] = edge1;
			(*it)->point.x = x;
			(*it)->point.y = y;
			(*it)->next = existing_node;
			return;
		}
		else
		{
			it = &((*it)->next);
		}
	}
}

static void add_st_edge(st_node **st, it_node **it, edge_node *edge,
						double dy)
{
	st_node *existing_node;
	double   den, r, x, y;

	while (true)
	{
		if (!*st)
		{
			// Append edge onto the tail end of the ST
			*st = (st_node *) MALLOC(sizeof(st_node), "ST insertion");
			(*st)->edge = edge;
			(*st)->xb = edge->xb;
			(*st)->xt = edge->xt;
			(*st)->dx = edge->dx;
			(*st)->prev = NULL;
			return;
		}
		else
		{
			den= ((*st)->xt - (*st)->xb) - (edge->xt - edge->xb);

			// If new edge and ST edge don't cross
			if ((edge->xt >= (*st)->xt) || (edge->dx == (*st)->dx)
				|| (fabs(den) <= DBL_EPSILON))
			{
				// No intersection - insert edge here (before the ST edge)
				existing_node = *st;
				*st = (st_node *) MALLOC(sizeof(st_node), "ST insertion");
				(*st)->edge = edge;
				(*st)->xb = edge->xb;
				(*st)->xt = edge->xt;
				(*st)->dx = edge->dx;
				(*st)->prev = existing_node;
				return;
			}
			else
			{
				// Compute intersection between new edge and ST edge
				r = (edge->xb - (*st)->xb) / den;
				x = (*st)->xb + r * ((*st)->xt - (*st)->xb);
				y = r * dy;

				// Insert the edge pointers and the intersection point in the IT
				add_intersection(it, (*st)->edge, edge, x, y);

				// Head further into the ST
				st = &((*st)->prev);
			}
		}
	}
}


static void build_intersection_table(it_node **it, edge_node *aet, double dy)
{
	st_node   *st, *stp;
	edge_node *edge;

	// Build intersection table for the current scanbeam
	reset_it(it);
	st = NULL;

	// Process each AET edge
	for (edge = aet; edge; edge = edge->next)
	{
		if ((edge->bstate[ABOVE] == BUNDLE_HEAD) ||
			edge->bundle[ABOVE][CLIP] || edge->bundle[ABOVE][SUBJ])
			add_st_edge(&st, it, edge, dy);
	}

	// Free the sorted edge table
	while (st)
	{
		stp = st->prev;
		FREE(st);
		st = stp;
	}
}

static int count_contours(polygon_node *polygon)
{
	int          nc, nv;
	vertex_node *v, *nextv;

	for (nc= 0; polygon; polygon= polygon->next)
		if (polygon->active)
		{
			// Count the vertices in the current contour
			nv= 0;
			for (v = polygon->proxy->v[LEFT]; v; v= v->next)
				nv++;

			// Record valid vertex counts in the active field
			if (nv > 2)
			{
				polygon->active= nv;
				nc++;
			}
			else
			{
				// Invalid contour: just free the heap
				for (v = polygon->proxy->v[LEFT]; v; v= nextv)
				{
					nextv = v->next;
					FREE(v);
				}
				polygon->active= 0;
			}
		}
		return nc;
}


static void add_left(polygon_node *p, double x, double y)
{
	if (p == NULL) return;

	vertex_node *nv;

	// Create a new vertex node and set its fields
	nv = (vertex_node *) MALLOC(sizeof(vertex_node), "vertex node creation");
	nv->x = x;
	nv->y = y;

	// Add vertex nv to the left end of the polygon's vertex list
	nv->next = p->proxy->v[LEFT];

	// Update proxy->[LEFT] to point to nv
	p->proxy->v[LEFT] = nv;
}


static void merge_left(polygon_node *p, polygon_node *q, polygon_node *list)
{
	if (p == NULL || q == NULL) return;

	polygon_node *target;

	// Label contour as a hole
	q->proxy->hole = true;

	if (p->proxy != q->proxy)
	{
		// Assign p's vertex list to the left end of q's list
		p->proxy->v[RIGHT]->next = q->proxy->v[LEFT];
		q->proxy->v[LEFT] = p->proxy->v[LEFT];

		// Redirect any p->proxy references to q->proxy

		for (target = p->proxy; list; list= list->next)
		{
			if (list->proxy == target)
			{
				list->active = FALSE;
				list->proxy = q->proxy;
			}
		}
	}
}

static void add_right(polygon_node *p, double x, double y)
{
	if (p == NULL) return;

	vertex_node *nv;

	// Create a new vertex node and set its fields
	nv = (vertex_node *) MALLOC(sizeof(vertex_node), "vertex node creation");
	nv->x = x;
	nv->y = y;
	nv->next = NULL;

	// Add vertex nv to the right end of the polygon's vertex list
	p->proxy->v[RIGHT]->next = nv;

	// Update proxy->v[RIGHT] to point to nv
	p->proxy->v[RIGHT] = nv;
}


static void merge_right(polygon_node *p, polygon_node *q, polygon_node *list)
{
	if (p == NULL || q == NULL) return;

	polygon_node *target;

	// Label contour as external
	q->proxy->hole = false;

	if (p->proxy != q->proxy)
	{
		// Assign p's vertex list to the right end of q's list
		q->proxy->v[RIGHT]->next = p->proxy->v[LEFT];
		q->proxy->v[RIGHT] = p->proxy->v[RIGHT];

		// Redirect any p->proxy references to q->proxy
		for (target= p->proxy; list; list= list->next)
		{
			if (list->proxy == target)
			{
				list->active = FALSE;
				list->proxy = q->proxy;
			}
		}
	}
}

static void add_local_min(polygon_node **p, edge_node *edge,
						  double x, double y)
{
	polygon_node *existing_min;
	vertex_node  *nv;

	existing_min= *p;

	*p = (polygon_node *) MALLOC(sizeof(polygon_node), "polygon node creation");

	// Create a new vertex node and set its fields
	nv = (vertex_node *) MALLOC(sizeof(vertex_node), "vertex node creation");
	nv->x = x;
	nv->y = y;
	nv->next = NULL;

	// Initialize proxy to point to p itself
	(*p)->proxy = (*p);
	(*p)->active = TRUE;
	(*p)->next = existing_min;

	// Make v[LEFT] and v[RIGHT] point to new vertex nv
	(*p)->v[LEFT] = nv;
	(*p)->v[RIGHT] = nv;

	// Assign polygon p to the edge
	edge->outp[ABOVE] = *p;
}

static gf_box *create_contour_bboxes(gf_polygon *p)
{
	gf_box	*box;
	int		c, v;

	box = (gf_box *) MALLOC(p->nContours * sizeof(gf_box), "Bounding box creation");

	// Construct contour bounding boxes
	for (c= 0; c < p->nContours; c++)
	{
		if(p->contour[c].nPoints == 0)
			continue;

		box[c].ll.x = p->contour[c].point[0].x;
		box[c].ll.y = p->contour[c].point[0].y;
		box[c].ur.x = p->contour[c].point[0].x;
		box[c].ur.y = p->contour[c].point[0].y;

		for (v= 1; v < p->contour[c].nPoints; v++)
		{
			box[c].ll.x = AMIN(box[c].ll.x, p->contour[c].point[v].x);
			box[c].ll.y = AMIN(box[c].ll.y, p->contour[c].point[v].y);
			box[c].ur.x = AMAX(box[c].ur.x, p->contour[c].point[v].x);
			box[c].ur.y = AMAX(box[c].ur.y, p->contour[c].point[v].y);
		}
	}
	return box;
}


static void minimax_test(gf_polygon *target, gf_polygon *source, gf_boolean op)
{
	gf_box *s_bbox, *c_bbox;
	int   s, c, *o_table, overlap;

	s_bbox = create_contour_bboxes(target);
	c_bbox = create_contour_bboxes(source);

	o_table = (int *) MALLOC(target->nContours * source->nContours * sizeof(int),
		"overlap table creation");

	// Check all target contour bounding boxes against source boxes
	for (s = 0; s < target->nContours; s++)
		for (c = 0; c < source->nContours; c++)
			o_table[c * target->nContours + s] =
				(!((s_bbox[s].ur.x < c_bbox[c].ll.x) ||
				   (s_bbox[s].ll.x > c_bbox[c].ur.x))) &&
				(!((s_bbox[s].ur.y < c_bbox[c].ll.y) ||
				   (s_bbox[s].ll.y > c_bbox[c].ur.y)));

	// For each source contour, search for any target contour overlaps
	for (c = 0; c < source->nContours; c++)
	{
		overlap = 0;
		for (s = 0; (!overlap) && (s < target->nContours); s++)
			overlap = o_table[c * target->nContours + s];

		if (!overlap)
			// Flag non contributing status by negating vertex count
			source->contour[c].nPoints = -source->contour[c].nPoints;
	}

	if (op == GF_INT)
	{
		// For each target contour, search for any source contour overlaps
		for (s = 0; s < target->nContours; s++)
		{
			overlap = 0;
			for (c = 0; (!overlap) && (c < source->nContours); c++)
				overlap = o_table[c * target->nContours + s];

			if (!overlap)
				// Flag non contributing status by negating vertex count
				target->contour[s].nPoints = -target->contour[s].nPoints;
		}
	}

	FREE(s_bbox);
	FREE(c_bbox);
	FREE(o_table);
}

gf_polygon *MergeOverlappedContours(gf_polygon *p, int from, int nContours)
{
	gf_polygon *target, *source, *result;
	gf_box		bound;

	target = gfCopyContours2Polygon(p, from, 1);

	bound.ll.x = target->box[0].ll.x;
	bound.ll.y = target->box[0].ll.y;
	bound.ur.x = target->box[0].ur.x;
	bound.ur.y = target->box[0].ur.y;

	int c;
	for (c = 1; c < nContours; c++)
	{
		int d;
		for (d = c; d < nContours; d++)
		{
		if (gfBoxesOverlap(bound, p->box[from+d])) break;
			bound = gfGetBoxOfBoxes(bound, p->box[from+d]);
		}
		if (d > c)
		{
			gfAttachContours2Polygon(p, from + c, d - c, target);
			int k;
			for (k = c; k < d; k++)
			{
				// prevents contour[j].vertex from being freed by gfFreePolygon
				p->contour[k+from].nPoints = 0;
			}
		}
		if (d < nContours)
		{
			source = gfAssignContours2Polygon(p, from + d, 1);
			result = gfAllocPolygon(0);
			gfPolygonBool(GF_UNION, target, source, result);
			gfFreePolygon(target);
			FREE(source);
			target = result;
			bound = gfGetBoxOfBoxes(bound, p->box[from+d]);
		}
		c = d;
	}
	return target;
}

gf_polygon *MergeContoursByChunk(gf_polygon *p, int from, int nContours)
{
	// Merge contours by chunk to improve th performance

	// Determine number of chunks
	int chunkCount = (nContours - 1) / CONTOUR_MERGE_CHUNK + 1;

	// Allocate an array of divided sets of contours
	gf_polygon **chunkPoly = (gf_polygon **) malloc(chunkCount * sizeof(gf_polygon *));
	gf_box *chunkBox = (gf_box *) malloc(chunkCount * sizeof(gf_box));

	// Divide the polygon into smaller sets of contours
	int chunk = CONTOUR_MERGE_CHUNK;
	int c = 0;
	int	k = 0;
	while (c < nContours)
	{
		if (c + chunk > nContours) chunk = nContours - c;

		if (chunk == 1)
		{
			// Just copy the original contour if there is only one contour in a chunk
			chunkPoly[k] = gfCopyContours2Polygon(p, from + c, 1);
		}
		else
		{
			// Merge all contours within a chunk
			chunkPoly[k] = MergeOverlappedContours(p, from + c, chunk);
			gfSetBoundingBox(chunkPoly[k]);
		}
		if (chunkPoly[k])
		{
			chunkBox[k] = gfGetBoxOfBoxes(chunkPoly[k]->box, chunkPoly[k]->nContours);
		}

		k++;
		c += chunk;
	}

	// Simply join those polygons if their bounding boxes are not overlapping

	for (k = 0; k < chunkCount; k++)
	{
		if (!chunkPoly[k]) continue;

		int m;
		for (m = k + 1; m < chunkCount; m++)
		{
			if (!chunkPoly[m]) continue;

			if (!gfBoxesOverlap(chunkBox[k], chunkBox[m]))
			{
				gfAttachContours2Polygon(chunkPoly[m], 0, chunkPoly[m]->nContours, chunkPoly[k]);
				gfUpdateBoxOfBoxes(&(chunkBox[k]), chunkPoly[m]->box, chunkPoly[m]->nContours);
				int i;
				for (i = 0; i < chunkPoly[m]->nContours; i++)
				{
					// prevents contour[j].vertex from being freed by gfFreePolygon
					chunkPoly[m]->contour[i].nPoints = 0;
				}
				gfFreeContours(chunkPoly[m]);
			}
		}
	}
	FREE (chunkBox);

	// Merge the result polygons of each chunk

	gf_polygon *temp;
	gf_polygon *result;

	result = chunkPoly[0];

	for (k = 1; k < chunkCount; k++)
	{
		if (!chunkPoly[k]) continue;

		temp = gfAllocPolygon(0);
		gfPolygonBool(GF_UNION, result, chunkPoly[k], temp);
		gfFreePolygon(result);
		gfFreePolygon(chunkPoly[k]);
		result = temp;
	}

	FREE(chunkPoly);

	return result;
}

gf_polygon *MergeContours(gf_polygon *p, int from, int nContours)
{
	// This is a private function that merges (nContours) solid contours
	// starting from (from) contour in polygon (p).

	gf_polygon *target, *source, *result;

	// Check for bounding box overlap
	bool *overlap = (bool *) malloc(nContours * sizeof(bool));
	int c;
	for (c = 0; c < nContours; c++)
	{
		overlap[c] = false;
	}

	for (c = 0; c < nContours - 1; c++)
	{
		if (overlap[c]) continue;

		int d;
		for (d = c + 1; d < nContours; d++)
		{
			if (overlap[d]) continue;

			if (gfBoxesOverlap(p->box[from+c], p->box[from+d]))
			{
				overlap[c] = true;
				overlap[d] = true;
			}
		}
	}

	// Count 
	int nOverlapped = 0;
	for (c = 0; c < nContours; c++)
	{
		if (overlap[c]) nOverlapped++;
	}

	target = gfAllocPolygon(nContours - nOverlapped);
	if (nOverlapped > 0)
		source = gfAllocPolygon(nOverlapped);
	else
		source = NULL;

	int s = 0;
	int t = 0;
	for (c = 0; c < nContours; c++)
	{
		int n = from + c;
		if (overlap[c])
		{
			source->hole[s] = false;
			source->contour[s].nPoints = p->contour[n].nPoints;
			source->contour[s].point = p->contour[n].point;
			source->box[s].ll.x = p->box[n].ll.x;
			source->box[s].ll.y = p->box[n].ll.y;
			source->box[s].ur.x = p->box[n].ur.x;
			source->box[s].ur.y = p->box[n].ur.y;
			s++;
		}
		else
		{
			target->hole[t] = false;
			target->contour[t].nPoints = p->contour[n].nPoints;
			target->contour[t].point = p->contour[n].point;
			target->box[t].ll.x = p->box[n].ll.x;
			target->box[t].ll.y = p->box[n].ll.y;
			target->box[t].ur.x = p->box[n].ur.x;
			target->box[t].ur.y = p->box[n].ur.y;
			t++;
		}
		// To prevent contour points get freed by gfFreePolygon
		p->contour[n].nPoints = 0;
	}

	FREE (overlap);

	if (source)
	{
		if (source->nContours <= CONTOUR_MERGE_CHUNK)
			result = MergeOverlappedContours(source, 0, source->nContours);
		else
			result = MergeContoursByChunk(source, 0, source->nContours);

		gfFreePolygon(source);
		source = result;

		result = gfAllocPolygon(0);
		gfPolygonBool(GF_UNION, target, source, result);
		gfFreePolygon(target);
		gfFreePolygon(source);
		target = result;
	}

	return target;
}

//=================
// Public Functions
//=================

void gfFreePoints(gf_point *&p)
{
	FREE(p);
}

void gfFreeContours(gf_polygon *&p)
{
	if (p)
	{
		FREE(p->hole);
		FREE(p->contour);
		FREE(p->box);
		FREE(p);
	}
}

void gfFreePolygon(gf_polygon *&p, bool bFreeContoursToo /* = true */)
{
	if (bFreeContoursToo)
		FreePolygon(p);

	FREE(p);
}

bool gfPolygonBoolWithoutgfFreePolygon(gf_boolean op, gf_polygon *target, gf_polygon *source,
				   gf_polygon *result)
{
	sb_tree       *sbtree= NULL;
	it_node       *it= NULL, *intersect;
	edge_node     *edge, *prev_edge, *next_edge, *succ_edge, *e0, *e1;
	edge_node     *aet= NULL, *c_heap= NULL, *s_heap= NULL;
	lmt_node      *lmt= NULL, *local_min;
	polygon_node  *out_poly= NULL, *p, *q, *poly, *npoly, *cf= NULL;
	vertex_node   *vtx, *nv;
	h_state        horiz[2];
	int            in[2], exists[2], parity[2] = {LEFT, LEFT};
	int            c, v, contributing = 0 , search, scanbeam = 0, sbt_entries = 0;
	int            vclass, bl = 0, br = 0, tl = 0, tr = 0;

	double        *sbt=NULL, xb, px, yb, yt=DBL_MAX, dy = 0 , ix, iy;

	bHoleExist = false;

	// Test for trivial NULL result cases
	if (((target->nContours == 0) && (source->nContours == 0))
	 || ((target->nContours == 0) && ((op == GF_INT) || (op == GF_SUB)))
	 || ((source->nContours == 0) &&  (op == GF_INT)))
	{
		result->nContours = 0;
		result->hole = NULL;
		result->contour = NULL;
		result->box = NULL;
		return false;
	}

	// Identify potentially contributing contours
	if (((op == GF_INT) || (op == GF_SUB))
		&& (target->nContours > 0) && (source->nContours > 0))
		minimax_test(target, source, op);

	// Build LMT
	if (target->nContours > 0)
		s_heap = build_lmt(&lmt, &sbtree, &sbt_entries, target, SUBJ, op);
	if (source->nContours > 0)
		c_heap = build_lmt(&lmt, &sbtree, &sbt_entries, source, CLIP, op);

	// Return a NULL result if no contours contribute
	if (lmt == NULL)
	{
		result->nContours = 0;
		result->hole = NULL;
		result->contour = NULL;
		result->box = NULL;
		reset_lmt(&lmt);
		FREE(s_heap);
		FREE(c_heap);
		return false;
	}

	// Build scanbeam table from scanbeam tree
	sbt = (double *) MALLOC(sbt_entries * sizeof(double), "sbt creation");
	build_sbt(&scanbeam, sbt, sbtree);
	scanbeam = 0;
	free_sbtree(&sbtree);


	// Invert source polygon for difference operation
	if (op == GF_SUB) parity[CLIP] = RIGHT;

	local_min= lmt;

	// ====== Process each scanbeam ======
	while (scanbeam < sbt_entries)
	{
		// Set yb and yt to the bottom and top of the scanbeam
		yb = sbt[scanbeam++];
		if (scanbeam < sbt_entries)
		{
			yt = sbt[scanbeam];
			dy = yt - yb;
		}

		// === SCANBEAM BOUNDARY PROCESSING ===

		// If LMT node corresponding to yb exists
		if (local_min)
		{
			if (local_min->y == yb)
			{
				// Add edges starting at this local minimum to the AET
				for (edge = local_min->first_bound; edge; edge = edge->next_bound)
					add_edge_to_aet(&aet, edge, NULL);

				local_min = local_min->next;
			}
		}

		// Set dummy previous x value
		px = -DBL_MAX;

		// Create bundles within AET
		e0 = aet;
		e1 = aet;

		// Set up bundle fields of first edge
		aet->bundle[ABOVE][ aet->type] = (aet->top.y != yb);
		aet->bundle[ABOVE][!aet->type] = FALSE;
		aet->bstate[ABOVE] = UNBUNDLED;

		for (next_edge = aet->next; next_edge; next_edge = next_edge->next)
		{
			// Set up bundle fields of next edge
			next_edge->bundle[ABOVE][ next_edge->type] = (next_edge->top.y != yb);
			next_edge->bundle[ABOVE][!next_edge->type] = FALSE;
			next_edge->bstate[ABOVE] = UNBUNDLED;

			// Bundle edges above the scanbeam boundary if they coincide
			if (next_edge->bundle[ABOVE][next_edge->type])
			{
				if (GF_EQ(e0->xb, next_edge->xb) && GF_EQ(e0->dx, next_edge->dx)
					&& (e0->top.y != yb))
				{
					next_edge->bundle[ABOVE][ next_edge->type] ^=
						e0->bundle[ABOVE][ next_edge->type];
					next_edge->bundle[ABOVE][!next_edge->type] =
						e0->bundle[ABOVE][!next_edge->type];
					next_edge->bstate[ABOVE] = BUNDLE_HEAD;
					e0->bundle[ABOVE][CLIP] = FALSE;
					e0->bundle[ABOVE][SUBJ] = FALSE;
					e0->bstate[ABOVE] = BUNDLE_TAIL;
				}
				e0 = next_edge;
			}
		}

		horiz[CLIP] = NH;
		horiz[SUBJ] = NH;

		// Process each edge at this scanbeam boundary
		for (edge = aet; edge; edge = edge->next)
		{
			exists[CLIP] = edge->bundle[ABOVE][CLIP] +
						  (edge->bundle[BELOW][CLIP] << 1);
			exists[SUBJ] = edge->bundle[ABOVE][SUBJ] +
						  (edge->bundle[BELOW][SUBJ] << 1);

			if (exists[CLIP] || exists[SUBJ])
			{
				// Set bundle side
				edge->bside[CLIP] = parity[CLIP];
				edge->bside[SUBJ] = parity[SUBJ];

				// Determine contributing status and quadrant occupancies
				switch (op)
				{
				case GF_SUB:
				case GF_INT:
					contributing = (exists[CLIP] && (parity[SUBJ] || horiz[SUBJ]))
						|| (exists[SUBJ] && (parity[CLIP] || horiz[CLIP]))
						|| (exists[CLIP] && exists[SUBJ]
					&& (parity[CLIP] == parity[SUBJ]));
					br = (parity[CLIP])
					  && (parity[SUBJ]);
					bl = (parity[CLIP] ^ edge->bundle[ABOVE][CLIP])
					  && (parity[SUBJ] ^ edge->bundle[ABOVE][SUBJ]);
					tr = (parity[CLIP] ^ (horiz[CLIP]!=NH))
					  && (parity[SUBJ] ^ (horiz[SUBJ]!=NH));
					tl = (parity[CLIP] ^ (horiz[CLIP]!=NH) ^ edge->bundle[BELOW][CLIP])
					  && (parity[SUBJ] ^ (horiz[SUBJ]!=NH) ^ edge->bundle[BELOW][SUBJ]);
					break;
				case GF_XOR:
					contributing = exists[CLIP] || exists[SUBJ];
					br = (parity[CLIP])
					   ^ (parity[SUBJ]);
					bl = (parity[CLIP] ^ edge->bundle[ABOVE][CLIP])
					   ^ (parity[SUBJ] ^ edge->bundle[ABOVE][SUBJ]);
					tr = (parity[CLIP] ^ (horiz[CLIP]!=NH))
					   ^ (parity[SUBJ] ^ (horiz[SUBJ]!=NH));
					tl = (parity[CLIP] ^ (horiz[CLIP]!=NH) ^ edge->bundle[BELOW][CLIP])
					   ^ (parity[SUBJ] ^ (horiz[SUBJ]!=NH) ^ edge->bundle[BELOW][SUBJ]);
					break;
				case GF_UNION:
					contributing = (exists[CLIP] && (!parity[SUBJ] || horiz[SUBJ]))
								|| (exists[SUBJ] && (!parity[CLIP] || horiz[CLIP]))
								|| (exists[CLIP] && exists[SUBJ]
								&& (parity[CLIP] == parity[SUBJ]));
					br = (parity[CLIP])
					  || (parity[SUBJ]);
					bl = (parity[CLIP] ^ edge->bundle[ABOVE][CLIP])
					  || (parity[SUBJ] ^ edge->bundle[ABOVE][SUBJ]);
					tr = (parity[CLIP] ^ (horiz[CLIP]!=NH))
					  || (parity[SUBJ] ^ (horiz[SUBJ]!=NH));
					tl = (parity[CLIP] ^ (horiz[CLIP]!=NH) ^ edge->bundle[BELOW][CLIP])
					  || (parity[SUBJ] ^ (horiz[SUBJ]!=NH) ^ edge->bundle[BELOW][SUBJ]);
					break;
				}

				// Update parity
				parity[CLIP] ^= edge->bundle[ABOVE][CLIP];
				parity[SUBJ] ^= edge->bundle[ABOVE][SUBJ];

				// Update horizontal state
				if (exists[CLIP])
					horiz[CLIP] = next_h_state[horiz[CLIP]]
				[((exists[CLIP] - 1) << 1) + parity[CLIP]];
				if (exists[SUBJ])
					horiz[SUBJ] = next_h_state[horiz[SUBJ]]
				[((exists[SUBJ] - 1) << 1) + parity[SUBJ]];

				vclass = tr + (tl << 1) + (br << 2) + (bl << 3);

				if (contributing)
				{
					xb = edge->xb;

					switch (vclass)
					{
					case EMN:
					case IMN:
						add_local_min(&out_poly, edge, xb, yb);
						px = xb;
						cf = edge->outp[ABOVE];
						break;
					case ERI:
						if (xb != px)
						{
							add_right(cf, xb, yb);
							px = xb;
						}
						edge->outp[ABOVE] = cf;
						cf = NULL;
						break;
					case ELI:
						add_left(edge->outp[BELOW], xb, yb);
						px = xb;
						cf = edge->outp[BELOW];
						break;
					case EMX:
						if (xb != px)
						{
							add_left(cf, xb, yb);
							px = xb;
						}
						merge_right(cf, edge->outp[BELOW], out_poly);
						cf = NULL;
						break;
					case ILI:
						if (xb != px)
						{
							add_left(cf, xb, yb);
							px = xb;
						}
						edge->outp[ABOVE] = cf;
						cf = NULL;
						break;
					case IRI:
						add_right(edge->outp[BELOW], xb, yb);
						px = xb;
						cf = edge->outp[BELOW];
						edge->outp[BELOW] = NULL;
						break;
					case IMX:
						if (xb != px)
						{
							add_right(cf, xb, yb);
							px = xb;
						}
						merge_left(cf, edge->outp[BELOW], out_poly);
						cf = NULL;
						edge->outp[BELOW] = NULL;
						break;
					case IMM:
						if (xb != px)
						{
							add_right(cf, xb, yb);
							px = xb;
						}
						merge_left(cf, edge->outp[BELOW], out_poly);
						edge->outp[BELOW] = NULL;
						add_local_min(&out_poly, edge, xb, yb);
						cf = edge->outp[ABOVE];
						break;
					case EMM:
						if (xb != px)
						{
							add_left(cf, xb, yb);
							px = xb;
						}
						merge_right(cf, edge->outp[BELOW], out_poly);
						edge->outp[BELOW] = NULL;
						add_local_min(&out_poly, edge, xb, yb);
						cf= edge->outp[ABOVE];
						break;
					case LED:
						if (edge->bot.y == yb)
							add_left(edge->outp[BELOW], xb, yb);
						edge->outp[ABOVE] = edge->outp[BELOW];
						px = xb;
						break;
					case RED:
						if (edge->bot.y == yb)
							add_right(edge->outp[BELOW], xb, yb);
						edge->outp[ABOVE] = edge->outp[BELOW];
						px = xb;
						break;
					default:
						break;
					} // End of switch
				} // End of contributing conditional
			} // End of edge exists conditional
		} // End of AET loop

		// Delete terminating edges from the AET, otherwise compute xt
		for (edge = aet; edge; edge = edge->next)
		{
			if (edge->top.y == yb)
			{
				prev_edge = edge->prev;
				next_edge = edge->next;
				if (prev_edge)
					prev_edge->next = next_edge;
				else
					aet = next_edge;
				if (next_edge)
					next_edge->prev = prev_edge;

				// Copy bundle head state to the adjacent tail edge if required
				if ((edge->bstate[BELOW] == BUNDLE_HEAD) && prev_edge)
				{
					if (prev_edge->bstate[BELOW] == BUNDLE_TAIL)
					{
						prev_edge->outp[BELOW] = edge->outp[BELOW];
						prev_edge->bstate[BELOW] = UNBUNDLED;
						if (prev_edge->prev)
							if (prev_edge->prev->bstate[BELOW] == BUNDLE_TAIL)
								prev_edge->bstate[BELOW] = BUNDLE_HEAD;
					}
				}
			}
			else
			{
				if (edge->top.y == yt)
					edge->xt = edge->top.x;
				else
					edge->xt = edge->bot.x + edge->dx * (yt - edge->bot.y);
			}
		}

		if (scanbeam < sbt_entries)
		{
			// === SCANBEAM INTERIOR PROCESSING ===

			build_intersection_table(&it, aet, dy);

			// Process each node in the intersection table
			for (intersect = it; intersect; intersect = intersect->next)
			{
				e0 = intersect->ie[0];
				e1 = intersect->ie[1];

				// Only generate output for contributing intersections
				if ((e0->bundle[ABOVE][CLIP] || e0->bundle[ABOVE][SUBJ])
					&& (e1->bundle[ABOVE][CLIP] || e1->bundle[ABOVE][SUBJ]))
				{
					p = e0->outp[ABOVE];
					q = e1->outp[ABOVE];
					ix = intersect->point.x;
					iy = intersect->point.y + yb;

					in[CLIP] = ( e0->bundle[ABOVE][CLIP] && !e0->bside[CLIP])
							|| ( e1->bundle[ABOVE][CLIP] &&  e1->bside[CLIP])
							|| (!e0->bundle[ABOVE][CLIP] && !e1->bundle[ABOVE][CLIP]
									  && e0->bside[CLIP] && e1->bside[CLIP]);
					in[SUBJ] = ( e0->bundle[ABOVE][SUBJ] && !e0->bside[SUBJ])
							|| ( e1->bundle[ABOVE][SUBJ] &&  e1->bside[SUBJ])
							|| (!e0->bundle[ABOVE][SUBJ] && !e1->bundle[ABOVE][SUBJ]
									  && e0->bside[SUBJ] && e1->bside[SUBJ]);

					// Determine quadrant occupancies
					switch (op)
					{
					case GF_SUB:
					case GF_INT:
						tr = (in[CLIP])
						  && (in[SUBJ]);
						tl = (in[CLIP] ^ e1->bundle[ABOVE][CLIP])
						  && (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ]);
						br = (in[CLIP] ^ e0->bundle[ABOVE][CLIP])
						  && (in[SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
						bl = (in[CLIP] ^ e1->bundle[ABOVE][CLIP] ^ e0->bundle[ABOVE][CLIP])
						  && (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
						break;
					case GF_XOR:
						tr = (in[CLIP])
						   ^ (in[SUBJ]);
						tl = (in[CLIP] ^ e1->bundle[ABOVE][CLIP])
						   ^ (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ]);
						br = (in[CLIP] ^ e0->bundle[ABOVE][CLIP])
						   ^ (in[SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
						bl = (in[CLIP] ^ e1->bundle[ABOVE][CLIP] ^ e0->bundle[ABOVE][CLIP])
						   ^ (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
						break;
					case GF_UNION:
						tr = (in[CLIP])
						  || (in[SUBJ]);
						tl = (in[CLIP] ^ e1->bundle[ABOVE][CLIP])
						  || (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ]);
						br = (in[CLIP] ^ e0->bundle[ABOVE][CLIP])
						  || (in[SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
						bl = (in[CLIP] ^ e1->bundle[ABOVE][CLIP] ^ e0->bundle[ABOVE][CLIP])
						  || (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
						break;
					}

					vclass = tr + (tl << 1) + (br << 2) + (bl << 3);

					switch (vclass)
					{
					case EMN:
						add_local_min(&out_poly, e0, ix, iy);
						e1->outp[ABOVE] = e0->outp[ABOVE];
						break;
					case ERI:
						if (p)
						{
							add_right(p, ix, iy);
							e1->outp[ABOVE] = p;
							e0->outp[ABOVE] = NULL;
						}
						break;
					case ELI:
						if (q)
						{
							add_left(q, ix, iy);
							e0->outp[ABOVE] = q;
							e1->outp[ABOVE] = NULL;
						}
						break;
					case EMX:
						if (p && q)
						{
							add_left(p, ix, iy);
							merge_right(p, q, out_poly);
							e0->outp[ABOVE] = NULL;
							e1->outp[ABOVE] = NULL;
						}
						break;
					case IMN:
						add_local_min(&out_poly, e0, ix, iy);
						e1->outp[ABOVE] = e0->outp[ABOVE];
						break;
					case ILI:
						if (p)
						{
							add_left(p, ix, iy);
							e1->outp[ABOVE] = p;
							e0->outp[ABOVE] = NULL;
						}
						break;
					case IRI:
						if (q)
						{
							add_right(q, ix, iy);
							e0->outp[ABOVE] = q;
							e1->outp[ABOVE] = NULL;
						}
						break;
					case IMX:
						if (p && q)
						{
							add_right(p, ix, iy);
							merge_left(p, q, out_poly);
							e0->outp[ABOVE] = NULL;
							e1->outp[ABOVE] = NULL;
						}
						break;
					case IMM:
						if (p && q)
						{
							add_right(p, ix, iy);
							merge_left(p, q, out_poly);
							add_local_min(&out_poly, e0, ix, iy);
							e1->outp[ABOVE] = e0->outp[ABOVE];
						}
						break;
					case EMM:
						if (p && q)
						{
							add_left(p, ix, iy);
							merge_right(p, q, out_poly);
							add_local_min(&out_poly, e0, ix, iy);
							e1->outp[ABOVE] = e0->outp[ABOVE];
						}
						break;
					default:
						break;
					} // End of switch
				} // End of contributing intersection conditional

				// Swap bundle sides in response to edge crossing
				if (e0->bundle[ABOVE][CLIP]) e1->bside[CLIP] = !e1->bside[CLIP];
				if (e1->bundle[ABOVE][CLIP]) e0->bside[CLIP] = !e0->bside[CLIP];
				if (e0->bundle[ABOVE][SUBJ]) e1->bside[SUBJ] = !e1->bside[SUBJ];
				if (e1->bundle[ABOVE][SUBJ]) e0->bside[SUBJ] = !e0->bside[SUBJ];

				// Swap e0 and e1 bundles in the AET
				prev_edge = e0->prev;
				next_edge = e1->next;
				if (next_edge) next_edge->prev = e0;

				if (e0->bstate[ABOVE] == BUNDLE_HEAD)
				{
					search = TRUE;
					while (search)
					{
						prev_edge= prev_edge->prev;
						if (prev_edge)
						{
							if (prev_edge->bstate[ABOVE] != BUNDLE_TAIL)
								search = FALSE;
						}
						else
							search = FALSE;
					}
				}
				if (!prev_edge)
				{
					aet->prev = e1;
					e1->next = aet;
					aet = e0->next;
				}
				else
				{
					prev_edge->next->prev = e1;
					e1->next = prev_edge->next;
					prev_edge->next = e0->next;
				}
				e0->next->prev = prev_edge;
				e1->next->prev = e1;
				e0->next = next_edge;
			} // End of IT loop

			// Prepare for next scanbeam
			for (edge = aet; edge; edge = next_edge)
			{
				next_edge = edge->next;
				succ_edge = edge->succ;

				if ((edge->top.y == yt) && succ_edge)
				{
					// Replace AET edge by its successor
					succ_edge->outp[BELOW] = edge->outp[ABOVE];
					succ_edge->bstate[BELOW] = edge->bstate[ABOVE];
					succ_edge->bundle[BELOW][CLIP] = edge->bundle[ABOVE][CLIP];
					succ_edge->bundle[BELOW][SUBJ] = edge->bundle[ABOVE][SUBJ];
					prev_edge = edge->prev;
					if (prev_edge)
						prev_edge->next = succ_edge;
					else
						aet = succ_edge;
					if (next_edge)
						next_edge->prev = succ_edge;
					succ_edge->prev = prev_edge;
					succ_edge->next = next_edge;
				}
				else
				{
					// Update this edge
					edge->outp[BELOW] = edge->outp[ABOVE];
					edge->bstate[BELOW] = edge->bstate[ABOVE];
					edge->bundle[BELOW][CLIP] = edge->bundle[ABOVE][CLIP];
					edge->bundle[BELOW][SUBJ] = edge->bundle[ABOVE][SUBJ];
					edge->xb= edge->xt;
				}
				edge->outp[ABOVE] = NULL;
			}
		}
	} // === END OF SCANBEAM PROCESSING ====

	// Generate result polygon from out_poly
	result->contour = NULL;
	result->hole = NULL;
	result->box = NULL;
	result->nContours = count_contours(out_poly);
	if (result->nContours > 0)
	{
		gfAllocPolygon(result, result->nContours);

		c = 0;
		for (poly = out_poly; poly; poly = npoly)
		{
			npoly = poly->next;
			if (poly->active)
			{
				if (poly->proxy->hole) bHoleExist = true;
				result->hole[c] = poly->proxy->hole;
				result->contour[c].nPoints = poly->active;
				result->contour[c].point = (gf_point *) MALLOC(
					result->contour[c].nPoints * sizeof(gf_point),
					"vertex creation");

				result->box[c].ll.x = DBL_MAX;
				result->box[c].ll.y = DBL_MAX;
				result->box[c].ur.x = -DBL_MAX;
				result->box[c].ur.y = -DBL_MAX;
				v = result->contour[c].nPoints - 1;
				for (vtx = poly->proxy->v[LEFT]; vtx; vtx = nv)
				{
					nv = vtx->next;
					result->contour[c].point[v].x = vtx->x;
					result->contour[c].point[v].y = vtx->y;
					if (vtx->x < result->box[c].ll.x) result->box[c].ll.x = vtx->x;
					if (vtx->x > result->box[c].ur.x) result->box[c].ur.x = vtx->x;
					if (vtx->y < result->box[c].ll.y) result->box[c].ll.y = vtx->y;
					if (vtx->y > result->box[c].ur.y) result->box[c].ur.y = vtx->y;
					FREE(vtx);
					v--;
				}
				c++;
			}
			FREE(poly);
		}
	}

	// Free allocated memories
	reset_it(&it);
	reset_lmt(&lmt);
	FREE(c_heap);
	FREE(s_heap);
	FREE(sbt);

	return bHoleExist;
}

bool gfPolygonBool(gf_boolean op, gf_polygon *target, gf_polygon *source,
				   gf_polygon *result)
{
	sb_tree       *sbtree= NULL;
	it_node       *it= NULL, *intersect;
	edge_node     *edge, *prev_edge, *next_edge, *succ_edge, *e0, *e1;
	edge_node     *aet= NULL, *c_heap= NULL, *s_heap= NULL;
	lmt_node      *lmt= NULL, *local_min;
	polygon_node  *out_poly= NULL, *p, *q, *poly, *npoly, *cf= NULL;
	vertex_node   *vtx, *nv;
	h_state        horiz[2];
	int            in[2], exists[2], parity[2] = {LEFT, LEFT};
	int            c, v, contributing = 0 , search, scanbeam = 0, sbt_entries = 0;
	int            vclass, bl = 0, br = 0, tl = 0, tr = 0;

	double        *sbt=NULL, xb, px, yb, yt=DBL_MAX, dy = 0 , ix, iy;

	bHoleExist = false;

	// Test for trivial NULL result cases
	if (((target->nContours == 0) && (source->nContours == 0))
	 || ((target->nContours == 0) && ((op == GF_INT) || (op == GF_SUB)))
	 || ((source->nContours == 0) &&  (op == GF_INT)))
	{
		result->nContours = 0;
		result->hole = NULL;
		result->contour = NULL;
		result->box = NULL;
		return false;
	}

	// Identify potentially contributing contours
	if (((op == GF_INT) || (op == GF_SUB))
		&& (target->nContours > 0) && (source->nContours > 0))
		minimax_test(target, source, op);

	// Build LMT
	if (target->nContours > 0)
		s_heap = build_lmt(&lmt, &sbtree, &sbt_entries, target, SUBJ, op);
	if (source->nContours > 0)
		c_heap = build_lmt(&lmt, &sbtree, &sbt_entries, source, CLIP, op);

	// Return a NULL result if no contours contribute
	if (lmt == NULL)
	{
		result->nContours = 0;
		result->hole = NULL;
		result->contour = NULL;
		result->box = NULL;
		reset_lmt(&lmt);
		FREE(s_heap);
		FREE(c_heap);
		return false;
	}

	// Build scanbeam table from scanbeam tree
	sbt = (double *) MALLOC(sbt_entries * sizeof(double), "sbt creation");
	build_sbt(&scanbeam, sbt, sbtree);
	scanbeam = 0;
	free_sbtree(&sbtree);

	// Allow pointer re-use without causing memory leak
 	if (target == result) gfFreePolygon(target);
 	if (source == result) gfFreePolygon(source);

	// Invert source polygon for difference operation
	if (op == GF_SUB) parity[CLIP] = RIGHT;

	local_min= lmt;

	// ====== Process each scanbeam ======
	while (scanbeam < sbt_entries)
	{
		// Set yb and yt to the bottom and top of the scanbeam
		yb = sbt[scanbeam++];
		if (scanbeam < sbt_entries)
		{
			yt = sbt[scanbeam];
			dy = yt - yb;
		}

		// === SCANBEAM BOUNDARY PROCESSING ===

		// If LMT node corresponding to yb exists
		if (local_min)
		{
			if (local_min->y == yb)
			{
				// Add edges starting at this local minimum to the AET
				for (edge = local_min->first_bound; edge; edge = edge->next_bound)
					add_edge_to_aet(&aet, edge, NULL);

				local_min = local_min->next;
			}
		}

		// Set dummy previous x value
		px = -DBL_MAX;

		// Create bundles within AET
		e0 = aet;
		e1 = aet;

		// Set up bundle fields of first edge
		aet->bundle[ABOVE][ aet->type] = (aet->top.y != yb);
		aet->bundle[ABOVE][!aet->type] = FALSE;
		aet->bstate[ABOVE] = UNBUNDLED;

		for (next_edge = aet->next; next_edge; next_edge = next_edge->next)
		{
			// Set up bundle fields of next edge
			next_edge->bundle[ABOVE][ next_edge->type] = (next_edge->top.y != yb);
			next_edge->bundle[ABOVE][!next_edge->type] = FALSE;
			next_edge->bstate[ABOVE] = UNBUNDLED;

			// Bundle edges above the scanbeam boundary if they coincide
			if (next_edge->bundle[ABOVE][next_edge->type])
			{
				if (GF_EQ(e0->xb, next_edge->xb) && GF_EQ(e0->dx, next_edge->dx)
					&& (e0->top.y != yb))
				{
					next_edge->bundle[ABOVE][ next_edge->type] ^=
						e0->bundle[ABOVE][ next_edge->type];
					next_edge->bundle[ABOVE][!next_edge->type] =
						e0->bundle[ABOVE][!next_edge->type];
					next_edge->bstate[ABOVE] = BUNDLE_HEAD;
					e0->bundle[ABOVE][CLIP] = FALSE;
					e0->bundle[ABOVE][SUBJ] = FALSE;
					e0->bstate[ABOVE] = BUNDLE_TAIL;
				}
				e0 = next_edge;
			}
		}

		horiz[CLIP] = NH;
		horiz[SUBJ] = NH;

		// Process each edge at this scanbeam boundary
		for (edge = aet; edge; edge = edge->next)
		{
			exists[CLIP] = edge->bundle[ABOVE][CLIP] +
						  (edge->bundle[BELOW][CLIP] << 1);
			exists[SUBJ] = edge->bundle[ABOVE][SUBJ] +
						  (edge->bundle[BELOW][SUBJ] << 1);

			if (exists[CLIP] || exists[SUBJ])
			{
				// Set bundle side
				edge->bside[CLIP] = parity[CLIP];
				edge->bside[SUBJ] = parity[SUBJ];

				// Determine contributing status and quadrant occupancies
				switch (op)
				{
				case GF_SUB:
				case GF_INT:
					contributing = (exists[CLIP] && (parity[SUBJ] || horiz[SUBJ]))
						|| (exists[SUBJ] && (parity[CLIP] || horiz[CLIP]))
						|| (exists[CLIP] && exists[SUBJ]
					&& (parity[CLIP] == parity[SUBJ]));
					br = (parity[CLIP])
					  && (parity[SUBJ]);
					bl = (parity[CLIP] ^ edge->bundle[ABOVE][CLIP])
					  && (parity[SUBJ] ^ edge->bundle[ABOVE][SUBJ]);
					tr = (parity[CLIP] ^ (horiz[CLIP]!=NH))
					  && (parity[SUBJ] ^ (horiz[SUBJ]!=NH));
					tl = (parity[CLIP] ^ (horiz[CLIP]!=NH) ^ edge->bundle[BELOW][CLIP])
					  && (parity[SUBJ] ^ (horiz[SUBJ]!=NH) ^ edge->bundle[BELOW][SUBJ]);
					break;
				case GF_XOR:
					contributing = exists[CLIP] || exists[SUBJ];
					br = (parity[CLIP])
					   ^ (parity[SUBJ]);
					bl = (parity[CLIP] ^ edge->bundle[ABOVE][CLIP])
					   ^ (parity[SUBJ] ^ edge->bundle[ABOVE][SUBJ]);
					tr = (parity[CLIP] ^ (horiz[CLIP]!=NH))
					   ^ (parity[SUBJ] ^ (horiz[SUBJ]!=NH));
					tl = (parity[CLIP] ^ (horiz[CLIP]!=NH) ^ edge->bundle[BELOW][CLIP])
					   ^ (parity[SUBJ] ^ (horiz[SUBJ]!=NH) ^ edge->bundle[BELOW][SUBJ]);
					break;
				case GF_UNION:
					contributing = (exists[CLIP] && (!parity[SUBJ] || horiz[SUBJ]))
								|| (exists[SUBJ] && (!parity[CLIP] || horiz[CLIP]))
								|| (exists[CLIP] && exists[SUBJ]
								&& (parity[CLIP] == parity[SUBJ]));
					br = (parity[CLIP])
					  || (parity[SUBJ]);
					bl = (parity[CLIP] ^ edge->bundle[ABOVE][CLIP])
					  || (parity[SUBJ] ^ edge->bundle[ABOVE][SUBJ]);
					tr = (parity[CLIP] ^ (horiz[CLIP]!=NH))
					  || (parity[SUBJ] ^ (horiz[SUBJ]!=NH));
					tl = (parity[CLIP] ^ (horiz[CLIP]!=NH) ^ edge->bundle[BELOW][CLIP])
					  || (parity[SUBJ] ^ (horiz[SUBJ]!=NH) ^ edge->bundle[BELOW][SUBJ]);
					break;
				}

				// Update parity
				parity[CLIP] ^= edge->bundle[ABOVE][CLIP];
				parity[SUBJ] ^= edge->bundle[ABOVE][SUBJ];

				// Update horizontal state
				if (exists[CLIP])
					horiz[CLIP] = next_h_state[horiz[CLIP]]
				[((exists[CLIP] - 1) << 1) + parity[CLIP]];
				if (exists[SUBJ])
					horiz[SUBJ] = next_h_state[horiz[SUBJ]]
				[((exists[SUBJ] - 1) << 1) + parity[SUBJ]];

				vclass = tr + (tl << 1) + (br << 2) + (bl << 3);

				if (contributing)
				{
					xb = edge->xb;

					switch (vclass)
					{
					case EMN:
					case IMN:
						add_local_min(&out_poly, edge, xb, yb);
						px = xb;
						cf = edge->outp[ABOVE];
						break;
					case ERI:
						if (xb != px)
						{
							add_right(cf, xb, yb);
							px = xb;
						}
						edge->outp[ABOVE] = cf;
						cf = NULL;
						break;
					case ELI:
						add_left(edge->outp[BELOW], xb, yb);
						px = xb;
						cf = edge->outp[BELOW];
						break;
					case EMX:
						if (xb != px)
						{
							add_left(cf, xb, yb);
							px = xb;
						}
						merge_right(cf, edge->outp[BELOW], out_poly);
						cf = NULL;
						break;
					case ILI:
						if (xb != px)
						{
							add_left(cf, xb, yb);
							px = xb;
						}
						edge->outp[ABOVE] = cf;
						cf = NULL;
						break;
					case IRI:
						add_right(edge->outp[BELOW], xb, yb);
						px = xb;
						cf = edge->outp[BELOW];
						edge->outp[BELOW] = NULL;
						break;
					case IMX:
						if (xb != px)
						{
							add_right(cf, xb, yb);
							px = xb;
						}
						merge_left(cf, edge->outp[BELOW], out_poly);
						cf = NULL;
						edge->outp[BELOW] = NULL;
						break;
					case IMM:
						if (xb != px)
						{
							add_right(cf, xb, yb);
							px = xb;
						}
						merge_left(cf, edge->outp[BELOW], out_poly);
						edge->outp[BELOW] = NULL;
						add_local_min(&out_poly, edge, xb, yb);
						cf = edge->outp[ABOVE];
						break;
					case EMM:
						if (xb != px)
						{
							add_left(cf, xb, yb);
							px = xb;
						}
						merge_right(cf, edge->outp[BELOW], out_poly);
						edge->outp[BELOW] = NULL;
						add_local_min(&out_poly, edge, xb, yb);
						cf= edge->outp[ABOVE];
						break;
					case LED:
						if (edge->bot.y == yb)
							add_left(edge->outp[BELOW], xb, yb);
						edge->outp[ABOVE] = edge->outp[BELOW];
						px = xb;
						break;
					case RED:
						if (edge->bot.y == yb)
							add_right(edge->outp[BELOW], xb, yb);
						edge->outp[ABOVE] = edge->outp[BELOW];
						px = xb;
						break;
					default:
						break;
					} // End of switch
				} // End of contributing conditional
			} // End of edge exists conditional
		} // End of AET loop

		// Delete terminating edges from the AET, otherwise compute xt
		for (edge = aet; edge; edge = edge->next)
		{
			if (edge->top.y == yb)
			{
				prev_edge = edge->prev;
				next_edge = edge->next;
				if (prev_edge)
					prev_edge->next = next_edge;
				else
					aet = next_edge;
				if (next_edge)
					next_edge->prev = prev_edge;

				// Copy bundle head state to the adjacent tail edge if required
				if ((edge->bstate[BELOW] == BUNDLE_HEAD) && prev_edge)
				{
					if (prev_edge->bstate[BELOW] == BUNDLE_TAIL)
					{
						prev_edge->outp[BELOW] = edge->outp[BELOW];
						prev_edge->bstate[BELOW] = UNBUNDLED;
						if (prev_edge->prev)
							if (prev_edge->prev->bstate[BELOW] == BUNDLE_TAIL)
								prev_edge->bstate[BELOW] = BUNDLE_HEAD;
					}
				}
			}
			else
			{
				if (edge->top.y == yt)
					edge->xt = edge->top.x;
				else
					edge->xt = edge->bot.x + edge->dx * (yt - edge->bot.y);
			}
		}

		if (scanbeam < sbt_entries)
		{
			// === SCANBEAM INTERIOR PROCESSING ===

			build_intersection_table(&it, aet, dy);

			// Process each node in the intersection table
			for (intersect = it; intersect; intersect = intersect->next)
			{
				e0 = intersect->ie[0];
				e1 = intersect->ie[1];

				// Only generate output for contributing intersections
				if ((e0->bundle[ABOVE][CLIP] || e0->bundle[ABOVE][SUBJ])
					&& (e1->bundle[ABOVE][CLIP] || e1->bundle[ABOVE][SUBJ]))
				{
					p = e0->outp[ABOVE];
					q = e1->outp[ABOVE];
					ix = intersect->point.x;
					iy = intersect->point.y + yb;

					in[CLIP] = ( e0->bundle[ABOVE][CLIP] && !e0->bside[CLIP])
							|| ( e1->bundle[ABOVE][CLIP] &&  e1->bside[CLIP])
							|| (!e0->bundle[ABOVE][CLIP] && !e1->bundle[ABOVE][CLIP]
									  && e0->bside[CLIP] && e1->bside[CLIP]);
					in[SUBJ] = ( e0->bundle[ABOVE][SUBJ] && !e0->bside[SUBJ])
							|| ( e1->bundle[ABOVE][SUBJ] &&  e1->bside[SUBJ])
							|| (!e0->bundle[ABOVE][SUBJ] && !e1->bundle[ABOVE][SUBJ]
									  && e0->bside[SUBJ] && e1->bside[SUBJ]);

					// Determine quadrant occupancies
					switch (op)
					{
					case GF_SUB:
					case GF_INT:
						tr = (in[CLIP])
						  && (in[SUBJ]);
						tl = (in[CLIP] ^ e1->bundle[ABOVE][CLIP])
						  && (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ]);
						br = (in[CLIP] ^ e0->bundle[ABOVE][CLIP])
						  && (in[SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
						bl = (in[CLIP] ^ e1->bundle[ABOVE][CLIP] ^ e0->bundle[ABOVE][CLIP])
						  && (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
						break;
					case GF_XOR:
						tr = (in[CLIP])
						   ^ (in[SUBJ]);
						tl = (in[CLIP] ^ e1->bundle[ABOVE][CLIP])
						   ^ (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ]);
						br = (in[CLIP] ^ e0->bundle[ABOVE][CLIP])
						   ^ (in[SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
						bl = (in[CLIP] ^ e1->bundle[ABOVE][CLIP] ^ e0->bundle[ABOVE][CLIP])
						   ^ (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
						break;
					case GF_UNION:
						tr = (in[CLIP])
						  || (in[SUBJ]);
						tl = (in[CLIP] ^ e1->bundle[ABOVE][CLIP])
						  || (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ]);
						br = (in[CLIP] ^ e0->bundle[ABOVE][CLIP])
						  || (in[SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
						bl = (in[CLIP] ^ e1->bundle[ABOVE][CLIP] ^ e0->bundle[ABOVE][CLIP])
						  || (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
						break;
					}

					vclass = tr + (tl << 1) + (br << 2) + (bl << 3);

					switch (vclass)
					{
					case EMN:
						add_local_min(&out_poly, e0, ix, iy);
						e1->outp[ABOVE] = e0->outp[ABOVE];
						break;
					case ERI:
						if (p)
						{
							add_right(p, ix, iy);
							e1->outp[ABOVE] = p;
							e0->outp[ABOVE] = NULL;
						}
						break;
					case ELI:
						if (q)
						{
							add_left(q, ix, iy);
							e0->outp[ABOVE] = q;
							e1->outp[ABOVE] = NULL;
						}
						break;
					case EMX:
						if (p && q)
						{
							add_left(p, ix, iy);
							merge_right(p, q, out_poly);
							e0->outp[ABOVE] = NULL;
							e1->outp[ABOVE] = NULL;
						}
						break;
					case IMN:
						add_local_min(&out_poly, e0, ix, iy);
						e1->outp[ABOVE] = e0->outp[ABOVE];
						break;
					case ILI:
						if (p)
						{
							add_left(p, ix, iy);
							e1->outp[ABOVE] = p;
							e0->outp[ABOVE] = NULL;
						}
						break;
					case IRI:
						if (q)
						{
							add_right(q, ix, iy);
							e0->outp[ABOVE] = q;
							e1->outp[ABOVE] = NULL;
						}
						break;
					case IMX:
						if (p && q)
						{
							add_right(p, ix, iy);
							merge_left(p, q, out_poly);
							e0->outp[ABOVE] = NULL;
							e1->outp[ABOVE] = NULL;
						}
						break;
					case IMM:
						if (p && q)
						{
							add_right(p, ix, iy);
							merge_left(p, q, out_poly);
							add_local_min(&out_poly, e0, ix, iy);
							e1->outp[ABOVE] = e0->outp[ABOVE];
						}
						break;
					case EMM:
						if (p && q)
						{
							add_left(p, ix, iy);
							merge_right(p, q, out_poly);
							add_local_min(&out_poly, e0, ix, iy);
							e1->outp[ABOVE] = e0->outp[ABOVE];
						}
						break;
					default:
						break;
					} // End of switch
				} // End of contributing intersection conditional

				// Swap bundle sides in response to edge crossing
				if (e0->bundle[ABOVE][CLIP]) e1->bside[CLIP] = !e1->bside[CLIP];
				if (e1->bundle[ABOVE][CLIP]) e0->bside[CLIP] = !e0->bside[CLIP];
				if (e0->bundle[ABOVE][SUBJ]) e1->bside[SUBJ] = !e1->bside[SUBJ];
				if (e1->bundle[ABOVE][SUBJ]) e0->bside[SUBJ] = !e0->bside[SUBJ];

				// Swap e0 and e1 bundles in the AET
				prev_edge = e0->prev;
				next_edge = e1->next;
				if (next_edge) next_edge->prev = e0;

				if (e0->bstate[ABOVE] == BUNDLE_HEAD)
				{
					search = TRUE;
					while (search)
					{
						prev_edge= prev_edge->prev;
						if (prev_edge)
						{
							if (prev_edge->bstate[ABOVE] != BUNDLE_TAIL)
								search = FALSE;
						}
						else
							search = FALSE;
					}
				}
				if (!prev_edge)
				{
					aet->prev = e1;
					e1->next = aet;
					aet = e0->next;
				}
				else
				{
					prev_edge->next->prev = e1;
					e1->next = prev_edge->next;
					prev_edge->next = e0->next;
				}
				if(e0->next != NULL)
				{
					e0->next->prev = prev_edge;
				}
				if(e1->next != NULL)
				{
					e1->next->prev = e1;
				}

// 				e1->next->prev = e1;
// 				e0->next->prev = prev_edge;
				e0->next = next_edge;
			} // End of IT loop

			// Prepare for next scanbeam
			for (edge = aet; edge; edge = next_edge)
			{
				next_edge = edge->next;
				succ_edge = edge->succ;

				if ((edge->top.y == yt) && succ_edge)
				{
					// Replace AET edge by its successor
					succ_edge->outp[BELOW] = edge->outp[ABOVE];
					succ_edge->bstate[BELOW] = edge->bstate[ABOVE];
					succ_edge->bundle[BELOW][CLIP] = edge->bundle[ABOVE][CLIP];
					succ_edge->bundle[BELOW][SUBJ] = edge->bundle[ABOVE][SUBJ];
					prev_edge = edge->prev;
					if (prev_edge)
						prev_edge->next = succ_edge;
					else
						aet = succ_edge;
					if (next_edge)
						next_edge->prev = succ_edge;
					succ_edge->prev = prev_edge;
					succ_edge->next = next_edge;
				}
				else
				{
					// Update this edge
					edge->outp[BELOW] = edge->outp[ABOVE];
					edge->bstate[BELOW] = edge->bstate[ABOVE];
					edge->bundle[BELOW][CLIP] = edge->bundle[ABOVE][CLIP];
					edge->bundle[BELOW][SUBJ] = edge->bundle[ABOVE][SUBJ];
					edge->xb= edge->xt;
				}
				edge->outp[ABOVE] = NULL;
			}
		}
	} // === END OF SCANBEAM PROCESSING ====

	// Generate result polygon from out_poly
	result->contour = NULL;
	result->hole = NULL;
	result->box = NULL;
	result->nContours = count_contours(out_poly);
	if (result->nContours > 0)
	{
		gfAllocPolygon(result, result->nContours);

		c = 0;
		for (poly = out_poly; poly; poly = npoly)
		{
			npoly = poly->next;
			if (poly->active)
			{
				if (poly->proxy->hole) bHoleExist = true;
				result->hole[c] = poly->proxy->hole;
				result->contour[c].nPoints = poly->active;
				result->contour[c].point = (gf_point *) MALLOC(
					result->contour[c].nPoints * sizeof(gf_point),
					"vertex creation");

				result->box[c].ll.x = DBL_MAX;
				result->box[c].ll.y = DBL_MAX;
				result->box[c].ur.x = -DBL_MAX;
				result->box[c].ur.y = -DBL_MAX;
				v = result->contour[c].nPoints - 1;
				for (vtx = poly->proxy->v[LEFT]; vtx; vtx = nv)
				{
					nv = vtx->next;
					result->contour[c].point[v].x = vtx->x;
					result->contour[c].point[v].y = vtx->y;
					if (vtx->x < result->box[c].ll.x) result->box[c].ll.x = vtx->x;
					if (vtx->x > result->box[c].ur.x) result->box[c].ur.x = vtx->x;
					if (vtx->y < result->box[c].ll.y) result->box[c].ll.y = vtx->y;
					if (vtx->y > result->box[c].ur.y) result->box[c].ur.y = vtx->y;
					FREE(vtx);
					v--;
				}
				c++;
			}
			FREE(poly);
		}
	}

	// Free allocated memories
	reset_it(&it);
	reset_lmt(&lmt);
	FREE(c_heap);
	FREE(s_heap);
	FREE(sbt);

	return bHoleExist;
}

void gfAllocPolygon(gf_polygon *p, int nContours)
{
	p->nContours = nContours;

	p->nAllocated = nContours;
	p->hole = (bool *) MALLOC(nContours * sizeof(bool), "hole flag array creation");
	p->box = (gf_box *) MALLOC(nContours * sizeof(gf_box), "box creation");
	p->contour = (gf_contour *) MALLOC(nContours * sizeof(gf_contour), "contour creation");
	int c;
	for (c = 0; c < nContours; c++)
	{
		p->contour[c].nPoints = 0;
		p->hole[c] = false;
	}
}

gf_polygon *gfAllocPolygon(int nContours)
{
	gf_polygon *p;
	p =  (gf_polygon *) MALLOC(sizeof(gf_polygon), "polygon MALLOC in gfAllocPolygon");
	gfAllocPolygon(p, nContours);

	return p;
}

void gfAllocContour(gf_polygon *p, int c, int nPoints)
{
	p->contour[c].nPoints = nPoints;
	p->contour[c].point = (gf_point *) MALLOC(nPoints * sizeof(gf_point),
		"vertex MALLOC in gfAllocContour");
}

gf_polygon *gfAssignContours2Polygon(gf_polygon *source, int from, int nContours)
{
	// Construct a polygon structure with a partial set of contours from the
	// existing polygon structure. ("nContours" from the "c" in the "source")
	// This function does not allocate memories for hole, contour or box. It
	// simply point to the existing data set.

	gf_polygon *result = gfAllocPolygon(0);

	result->nContours = nContours;
	result->nAllocated = nContours;
	result->hole = &(source->hole[from]);
	result->contour = &(source->contour[from]);
	result->box = &(source->box[from]);

	return result;
}

gf_polygon *gfCopyContours2Polygon(gf_polygon *source, int from, int nContours)
{
	// Copy a partial set of contours from the existing polygon structure
	// to a new polygon structure

	gf_polygon *result = gfAllocPolygon(nContours);

	int c, d, v;

	for (c = 0, d = from; c < nContours; c++, d++)
	{
		result->hole[c] = source->hole[d];

		gfAllocContour(result, c, abs(source->contour[d].nPoints));

		for (v = 0; v < source->contour[d].nPoints; v++)
		{
			result->contour[c].point[v].x = source->contour[d].point[v].x;
			result->contour[c].point[v].y = source->contour[d].point[v].y;
		}

		result->box[c].ll.x = source->box[d].ll.x;
		result->box[c].ll.y = source->box[d].ll.y;
		result->box[c].ur.x = source->box[d].ur.x;
		result->box[c].ur.y = source->box[d].ur.y;
	}

	return result;
}

void gfAttachContours2Polygon(gf_polygon *source, int from, int nContours, gf_polygon *target)
{
	// Attach a polygon structure with a partial set of contours from the
	// existing polygon structure ("nContours" from the "from" in the "source")
	// to the another existing structure ("target").
	// This function does not allocate memories for hole, contour or box. It
	// simply point to the existing "source" data set.

	int n, c;

	n = target->nContours;

	gfReAllocPolygon(target, nContours);

	for (c = 0; c < nContours; c++)
	{
		target->hole[n+c] = source->hole[from+c];
		target->contour[n+c].nPoints = source->contour[from+c].nPoints;
		target->contour[n+c].point = source->contour[from+c].point;
		target->box[n+c] = gfGetBoundingBox(target->contour[n+c]);
	}
}

void gfAppendContours2Polygon(gf_polygon *source, int from, int nContours, gf_polygon *target)
{
	// Append a partial set of contours from the existing polygon structure
	// ("nContours" from "from" of "source" to another existing polygon structure ("target").
	int n, c, v;

	n = target->nContours;

	gfReAllocPolygon(target, nContours);

	for (c = 0; c < nContours; c++)
	{
		target->hole[n+c] = source->hole[from+c];
		gfAllocContour(target, n + c, abs(source->contour[from+c].nPoints));
		target->contour[n+c].nPoints = source->contour[from+c].nPoints;

		if(target->contour[n+c].nPoints == 0)
			continue;

		for (v = 0; v < abs(target->contour[n+c].nPoints); v++)
		{
			target->contour[n+c].point[v].x = source->contour[from+c].point[v].x;
			target->contour[n+c].point[v].y = source->contour[from+c].point[v].y;
		}
		target->box[n+c] = gfGetBoundingBox(target->contour[n+c]);
	}
}

// Set polygon memory allocation size (MEM_ALLOC_CHUNK) and the
// number of contours to be boolean operated at a time (CONTOUR_MERGE_CHUNK).
//
// Use of MEM_ALLOC_CHUNK will make the memory allocation more efficient
// when the function "gfReAllocPolygon" used many times.
//
// Use of CONTOUR_MERGE_CHUNK will make the polygon boolean operation work
// faster than operating one at a time.

int MEM_ALLOC_CHUNK = 8;
int CONTOUR_MERGE_CHUNK = 100;

void gfSetChunkSize(int memChunk, int bopChunk)
{
	MEM_ALLOC_CHUNK = memChunk;
	CONTOUR_MERGE_CHUNK = bopChunk;
}

void gfSetChunkSize(int nContours)
{
	if (nContours <= 128)			MEM_ALLOC_CHUNK = 8;
	else if (nContours <= 1280)		MEM_ALLOC_CHUNK = 128;
	else if (nContours <= 40960)	MEM_ALLOC_CHUNK = 1024;
	else							MEM_ALLOC_CHUNK = 2048;

	CONTOUR_MERGE_CHUNK = 100;
	if (nContours > 10000)	CONTOUR_MERGE_CHUNK = 1000;
}

int gfGetMemChunkSize()
{
	return MEM_ALLOC_CHUNK;
}

int gfGetMergeChunkSize()
{
	return CONTOUR_MERGE_CHUNK;
}

void gfReAllocPolygon(gf_polygon *p, int moreContours)
{
	if (p == NULL || moreContours <= 0)
		return;

	// avoid memory fragments
	if ((p->nContours+moreContours) > p->nAllocated)
	{
// 		if (MEM_ALLOC_CHUNK < moreContours)
// 			p->nAllocated += moreContours;
// 		else
// 			p->nAllocated += MEM_ALLOC_CHUNK;

		int nChunks = (moreContours - 1) / MEM_ALLOC_CHUNK + 1;
		p->nAllocated += (nChunks * MEM_ALLOC_CHUNK);

		p->hole = (bool *) REALLOC(p->hole, p->nAllocated, sizeof(bool));
		p->box = (gf_box *) REALLOC(p->box, p->nAllocated, sizeof(gf_box));
		p->contour = (gf_contour *) REALLOC(p->contour, p->nAllocated,
			sizeof(gf_contour));
	}

	p->nContours += moreContours;
}

//===========================
//
// Polygon Merging Operations
//
//===========================

gf_polygon *gfPolygonMerge(gf_polygon *target, gf_polygon *source, bool bFreeSource /* = true */)
{
	// Merges solid contours in source to the target polygon and
	// subtract hole contours in source from the target polygon.
	//
	// This function is inefficient if the source polygon has
	// too many contours.

	gf_polygon *result, *temp;

	int c;
	for (c = 0; c < source->nContours; c++)
	{
		result = gfAllocPolygon(0);

		temp = gfAssignContours2Polygon(source, c, 1);

		if (source->hole[c])
			gfPolygonBool(GF_SUB, target, temp, result);
		else
			gfPolygonBool(GF_UNION, target, temp, result);

		FREE(temp);
		gfFreePolygon(target);
		target = result;
	}

	if (bFreeSource)
	{
		gfFreePolygon (source);
		source = NULL;
	}

	result = gfContourHierarchy(result);

	return result;
}


gf_polygon *gfPolygonMergeAll(gf_polygon *p)
{
	// Merges all contours within the polygon

	if (p->nContours <= 0) return p;	// return if there is no contour 

	gf_polygon *target, *source, *result;

	// Remove collinear points in all contours within the polygon
	vRemoveCollinearPoint(p);

	// Set bounding boxes of each contour within the polygon
	gfSetBoundingBox(p);

	// Merge contours by groups for efficiency
	target = NULL;

	int from;
	int c = 0;

	while (c < p->nContours)
	{
		// === Solid contours ===
		from = c;

		// collect consecutive solid contours
		for (; c < p->nContours; c++)
		{
			if (p->hole[c]) break;
		}

		// Merge the collected solid contours only
		if ((c - from) == 1)
		{
			source = gfCopyContours2Polygon(p, from, 1);
		}
		else
		{
			source = MergeContours(p, from, c - from);
		}

		// Merge the solid only merge result to the previously merged result
		if (target == NULL)
		{
			target = source;
		}
		else
		{
			result = gfAllocPolygon(0);
			gfPolygonBool(GF_UNION, target, source, result);
			gfFreePolygon(target);
			gfFreePolygon(source);
			target = result;
		}

		// === Hole contours ===
		from = c;

		// collect consecutive holes and change the type to solid for merging
		while (c < p->nContours && p->hole[c])
		{
			p->hole[c] = false;
			c++;
		}

		if (c > from && c <= p->nContours)
		{
			// Merge the collected hole contours only
			if ((c - from) == 1)
			{
				source = gfCopyContours2Polygon(p, from, 1);
			}
			else
			{
				source = MergeContours(p, from, c - from);
			}

			// Subtract the hole only merge result from the previously merged result
			result = gfAllocPolygon(0);
			gfPolygonBool(GF_SUB, target, source, result);
			gfFreePolygon(target);
			gfFreePolygon(source);
			target = result;
		}
	}

	// Remove collinear points in all contours in the result polygon
	vRemoveCollinearPoint(target);

	// Set bounding boxes of each contour in the result polygon
	gfSetBoundingBox(target);

	// Remove tiny holes
	target = pRemoveTinyHoles(target);

	// Organized the contours hierarchically
	target = gfContourHierarchy(target);

	return target;
}

bool PolygonsIntersect(gf_polygon *pPoly1, gf_polygon *pPoly2)
{
	if (!pPoly1 || pPoly1->nContours <= 0) return false;
	if (!pPoly2 || pPoly2->nContours <= 0) return false;

	gf_polygon *pResult = gfAllocPolygon(0);
	gfPolygonBool(GF_INT, pPoly1, pPoly2, pResult);

	bool bIntersect;
	if (pResult->nContours <= 0)
		bIntersect = false;
	else
		bIntersect = true;

	gfFreePolygon(pResult);

	return bIntersect;
}

bool PolygonInPolygon(gf_polygon *pObjPoly,
					  double xOffSet, double yOffSet, double dRot, double xLoc, double yLoc,	// pObjPoly transformation
					  gf_polygon *pPoly,
					  bool bCheckInside,
					  double *dX, double *dY)	// First violation point if FALSE
{
	// Checks to determine whether pObjPoly in completely inside (if bCheckInside = TRUE) or
	// completely outside (if bCheckInside = FALSE) of the pPoly.
	// 
	bool bInside;
	int c, v;
	
	if (xOffSet != 0.0 || yOffSet != 0)
	{
		for (c = 0; c < pObjPoly->nContours; c++)
		{
			vShiftContour(pObjPoly->contour[c].point, pObjPoly->contour[c].nPoints, xOffSet, yOffSet);
		}
	}
	if (dRot != 0.0)
	{
		for (c = 0; c < pObjPoly->nContours; c++)
		{
			vRotateContour(pObjPoly->contour[c].point, pObjPoly->contour[c].nPoints, dRot, 0.0, 0.0);
		}
	}
	if (xLoc != 0.0 || yLoc != 0)
	{
		for (c = 0; c < pObjPoly->nContours; c++)
		{
			vShiftContour(pObjPoly->contour[c].point, pObjPoly->contour[c].nPoints, xLoc, yLoc);
		}
	}
	
	for (c = 0; c < pObjPoly->nContours; c++)
	{
		for (v = 0; v < pObjPoly->contour[c].nPoints; v++)
		{
			bInside = PointInPolygon(pObjPoly->contour[c].point[v].x, pObjPoly->contour[c].point[v].y, pPoly);
			if ((bCheckInside && !bInside) ||
				(!bCheckInside && bInside))
			{
				*dX = pObjPoly->contour[c].point[v].x;
				*dY = pObjPoly->contour[c].point[v].y;
				return FALSE;
			}
		}
	}
	
	double x1, y1, x2, y2;
	double xa, ya, xb, yb;
	
	for (c = 0; c < pObjPoly->nContours; c++)
	{
		v = pObjPoly->contour[c].nPoints - 1;
		x1 = pObjPoly->contour[c].point[v].x;
		y1 = pObjPoly->contour[c].point[v].y;
		for (v = 0; v < pObjPoly->contour[c].nPoints; v++)
		{
			x2 = pObjPoly->contour[c].point[v].x;
			y2 = pObjPoly->contour[c].point[v].y;
			int cc, vv;
			for (cc = 0; cc < pPoly->nContours; cc++)
			{
				vv = pPoly->contour[cc].nPoints - 1;
				xa = pPoly->contour[cc].point[vv].x;
				ya = pPoly->contour[cc].point[vv].y;
				for (vv = 0; vv < pPoly->contour[cc].nPoints; vv++)
				{
					xb = pPoly->contour[cc].point[vv].x;
					yb = pPoly->contour[cc].point[vv].y;
					if (nTwoLines(x1, y1, x2, y2, xa, ya, xb, yb, dX, dY) == 3)
					{
						return FALSE;
					}
					xa = xb;
					ya = yb;
				}
			}
			x1 = x2;
			y1 = y2;
		}
	}
	
	return TRUE;
}

double Polygon2PolygonDistance(gf_polygon *pPoly1,gf_polygon *pPoly2,			// Two polygons
							   double dClearance,								// Stop calculating if a distance shorter than dClearance is found
							   double *x1, double *y1, double *x2, double *y2)	// Point coordinates between which the distance is calculated
{
	int c, v;
	
	double pX, pY, dX, dY;
	double xa, ya, xb, yb;
	double dist;
	double distMin = DBL_MAX;
	
	for (c = 0; c < pPoly1->nContours; c++)
	{
		for (v = 0; v < pPoly1->contour[c].nPoints; v++)
		{
			pX = pPoly1->contour[c].point[v].x;
			pY = pPoly1->contour[c].point[v].y;
			int cc, vv;
			for (cc = 0; cc < pPoly2->nContours; cc++)
			{
				vv = pPoly2->contour[cc].nPoints - 1;
				xa = pPoly2->contour[cc].point[vv].x;
				ya = pPoly2->contour[cc].point[vv].y;
				for (vv = 0; vv < pPoly1->contour[cc].nPoints; vv++)
				{
					xb = pPoly2->contour[cc].point[vv].x;
					yb = pPoly2->contour[cc].point[vv].y;
					dist = Point2LineDistance(pX, pY, xa, ya, xb, yb, false, &dX, &dY);
					if (dist < distMin)
					{
						distMin = dist;
						*x1 = pX;
						*y1 = pY;
						*x2 = dX;
						*y2 = dY;
						if (distMin < dClearance) return distMin;
					}
					xa = xb;
					ya = yb;
				}
			}
		}
	}
	return distMin;
}

bool SplitContour(gf_polygon *p, int c)
{
	bool bSplit = false;
	int u, v, i, j, k;
	int nPtOnLine = 0;
	int nPoints = p->contour[c].nPoints;
	if (nPoints <= 6) return bSplit;

	for (v = 0; v < nPoints; v++)
	{
		u = v + 3;
		if (u >= nPoints) u -= nPoints;
		int uEnd = v - 3;
		if (uEnd < 0) uEnd += nPoints;
		
		while (u != uEnd && u != v)
		{
			if (p->contour[c].point[u].x == p->contour[c].point[v].x &&
				p->contour[c].point[u].y == p->contour[c].point[v].y)
			{
				bSplit = true;
				break;
			}
// 			if (nPtOnLine < 2)
// 			{
// 				for (k = 0; k < v; k++)
// 				{
// 					if (PointOnLine(p->contour[c].point[u], p->contour[c].point[k], p->contour[c].point[k+1]))
// 					{
// 						nPtOnLine++;
// 					}
// 				}
// 			}
			u++;
			if (u >= nPoints) u -= nPoints;
		}
		if (bSplit) break;
	}
	if (!bSplit) return bSplit;

	gf_point *pList = (gf_point *) calloc(nPoints, sizeof(gf_point));
	for (k = 0; k < nPoints; k++)
	{
		pList[k].x = p->contour[c].point[k].x;
		pList[k].y = p->contour[c].point[k].y;
	}

	int nSplitPoints, nRemainedPoints;
	i = v + 1;
	if (i >= nPoints) i -= nPoints;
	j = u - 1;
	if (j < 0) j = nPoints - 1;
	if (pList[i].x == pList[j].x && pList[i].y == pList[j].y)
	{
		if (u > v)
			nSplitPoints = u - v - 2;
		else
			nSplitPoints = u - v - 2 + nPoints;
		
		nRemainedPoints = nPoints - nSplitPoints - 2;

		v++;
		if (v >= nPoints) v -= nPoints;
		nPtOnLine = 0;
	}
	else if (nPtOnLine > 1)
	{
		if (u > v)
			nSplitPoints = u - v - 1;
		else
			nSplitPoints = u - v - 1 + nPoints;

		nRemainedPoints = nPoints - nSplitPoints - 1;

		v++;
		if (v >= nPoints) v -= nPoints;
	}
	else
	{
		if (u > v)
			nSplitPoints = u - v - 1;
		else
			nSplitPoints = u - v - 1 + nPoints;

		nRemainedPoints = nPoints - nSplitPoints - 2;
		nPtOnLine = 0;
	}

	gfReAllocPolygon(p, 1);
	i = p->nContours - 1;
	p->hole[i] = true;
	gfAllocContour(p, i, nSplitPoints);

	gfFreePoints(p->contour[c].point);
	gfAllocContour(p, c, nRemainedPoints);
	
	for (j = 0; j < nSplitPoints; j++)
	{
		p->contour[i].point[j].x = pList[v].x;
		p->contour[i].point[j].y = pList[v].y;
		pList[v].x = DBL_MAX;
		v++;
		if (v >= nPoints) v -= nPoints;
	}
	pList[v].x = DBL_MAX;

	if (nPtOnLine == 0)
	{
		v++;
		if (v >= nPoints) v -= nPoints;
		pList[v].x = DBL_MAX;
	}

	j = 0;
	for (k = 0; k < nPoints; k++)
	{
		if (pList[k].x == DBL_MAX) continue;
		p->contour[c].point[j].x = pList[k].x;
		p->contour[c].point[j].y = pList[k].y;
		j++;
	}

	FREE(pList);
	return bSplit;
}

gf_polygon *SeparateHolesFromContour(gf_polygon *p)
{
	// Separate holes from a contour that defines holes by continuous line segments.
	// Example: Polygons or planes defined in Mentor Graphics designs.

	int k, n, c, u, v;

	// Pre-process each contour to add additional points if they are on line segments.

	for (c = 0; c < p->nContours; c++)
	{
		if (p->hole[c]) continue;

		gf_point *pList = (gf_point *) calloc(2 * p->contour[c].nPoints, sizeof(gf_point));
		u = p->contour[c].nPoints - 1;
		while (	p->contour[c].point[u].x == p->contour[c].point[0].x &&
				p->contour[c].point[u].y == p->contour[c].point[0].y)
		{
			p->contour[c].nPoints = u;
			u = p->contour[c].nPoints - 1;
		}
		n = 0;
		pList[n].x = p->contour[c].point[u].x;
		pList[n].y = p->contour[c].point[u].y;
		n++;
		for (v = 0; v < p->contour[c].nPoints - 1; v++)
		{
			if (v > 0)
			{
				if (p->contour[c].point[v].x == p->contour[c].point[v-1].x &&
					p->contour[c].point[v].y == p->contour[c].point[v-1].y) continue;
			}
			if (pList[n-1].x == p->contour[c].point[v].x)
			{
				for (k = p->contour[c].nPoints - 1; k > v; k--)
				{
					if (p->contour[c].point[k].x != pList[n-1].x) continue;
					if (p->contour[c].point[k].y <= pList[n-1].y &&
						p->contour[c].point[k].y <= p->contour[c].point[v].y) continue;
					if (p->contour[c].point[k].y >= pList[n-1].y &&
						p->contour[c].point[k].y >= p->contour[c].point[v].y) continue;
					if (k < p->contour[c].nPoints - 1)
					{
						if (p->contour[c].point[k].x == p->contour[c].point[k+1].x &&
							p->contour[c].point[k].y == p->contour[c].point[k+1].y) continue;
					}
					pList[n].x = p->contour[c].point[k].x;
					pList[n].y = p->contour[c].point[k].y;
					n++;
				}
			}
// 			for (k = p->contour[c].nPoints - 1; k > v; k--)
// 			{
// 				if (k < p->contour[c].nPoints - 1)
// 				{
// 					if (p->contour[c].point[k].x == p->contour[c].point[k+1].x &&
// 						p->contour[c].point[k].y == p->contour[c].point[k+1].y) continue;
// 				}
// 				if (p->contour[c].point[k].x == pList[n-1].x &&
// 					p->contour[c].point[k].y == pList[n-1].y) continue;
// 				if (p->contour[c].point[k].x == p->contour[c].point[v].x &&
// 					p->contour[c].point[k].y == p->contour[c].point[v].y) continue;
// 
// 				if (PointOnLine(p->contour[c].point[k], pList[n-1], p->contour[c].point[v]))
// 				{
// 					pList[n].x = p->contour[c].point[k].x;
// 					pList[n].y = p->contour[c].point[k].y;
// 					n++;
// 				}
// 			}
			pList[n].x = p->contour[c].point[v].x;
			pList[n].y = p->contour[c].point[v].y;
			n++;
		}

// 		if(n < 3)
// 		{
// 			gfFreePoints(pList);
// 			continue;
// 		}

		if (n != p->contour[c].nPoints)
		{
			gfFreePoints(p->contour[c].point);
			gfAllocContour(p, c, n);
			
			for (v = 0; v < p->contour[c].nPoints; v++)
			{
				p->contour[c].point[v].x = pList[v].x;
				p->contour[c].point[v].y = pList[v].y;
			}
		}
		gfFreePoints(pList);
	}

	// Search for holes in each contour

	int nContours = p->nContours;
	for (c = 0; c < nContours; c++)
	{
		while (SplitContour(p, c)) {};
		nContours = p->nContours;
	}

	return p;
}