#include <stdlib.h>
#include <vector>
#include <float.h>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "triangle.h"
#include <gl/glut.h>
#include <gl/GL.h>
#include <gl/GLU.h>

#include "TUYENPLE_Math.h"
#include "TUYENPLE_Polygon.h"
#include "TUYENPLE_Arc.h"
#include "TUYENPLE_Line.h"
#include "TUYENPLE_RectD.h"

#include "gfunc.h"

using namespace std;

Polygon *g_polygon = NULL;
DrawObjArray *g_scanlines = NULL;

// Only for OpenGL
double m_scale = 0.05;
double m_translateX = 0;
double m_translateY = 0;
double center_objectX = 0;
double center_objectY = 0;
int window_width = 1680;
int window_height = 1050;

struct DoublePoint1 
{
	double x;
	double y;

	DoublePoint1(double x_ = 0, double y_ = 0): x(x_), y(y_) {};

	friend inline bool operator== (const DoublePoint1& a, const DoublePoint1& b)
	{
		return abs(a.x-b.x) < 0.0001 && abs(a.y-b.y) < 0.0001;
	}
	friend inline bool operator!= (const DoublePoint1& a, const DoublePoint1& b)
	{
		return abs(a.x-b.x) > 0.0001 || abs(a.y-b.y) > 0.0001;
	}
	friend inline DoublePoint1 operator- (const DoublePoint1& a, const DoublePoint1& b)
	{
		return DoublePoint1(a.x-b.x, a.y-b.y);
	}
};

struct Triangle
{
	std::vector<int> apexes;

	std::vector<int> adjacentTriangles;
};

struct Edge
{
	int idxA;
	int idxB;

	Edge(int a_ = 0, int b_ = 0): idxA(a_), idxB(b_) {};

	friend inline bool operator== (const Edge& a, const Edge& b)
	{
		return (a.idxA==b.idxA && a.idxB==b.idxB) || (a.idxB==b.idxA && a.idxA==b.idxB);
	}
	friend inline bool operator!= (const Edge& a, const Edge& b)
	{
		return (a.idxA!=b.idxA && a.idxA!=b.idxB) || (a.idxB!=b.idxA && a.idxB!=b.idxB);
	}
};

struct Vertex
{
	DoublePoint1 point;
};

typedef std::vector< Edge *> Edges;
typedef std::vector< Vertex *> Vertices;
typedef std::vector<Triangle *> Triangles;

Vertices g_vertices;
Edges g_edges;
Triangles g_triangles;

void initializeOpengl(int width, int height)
{	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glViewport(0, 0, width, height); 

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// Calculate The Aspect Ratio Of The Window
	gluPerspective(45.0f,(GLfloat)width/(GLfloat)height, 0.1f, 100.0f);
}

void mykey(int key, int x, int y)
{
	double next_scale;
	switch (key)
	{
	case GLUT_KEY_PAGE_UP:
		m_scale += m_scale/2;
		break;
	case GLUT_KEY_PAGE_DOWN:
		next_scale = m_scale - m_scale/2;
		if (next_scale > 0)
			m_scale -= next_scale;
		break;
	case GLUT_KEY_UP:
		m_translateY -= 1/m_scale;
		break;
	case GLUT_KEY_DOWN:
		m_translateY += 1/m_scale;
		break;
	case GLUT_KEY_LEFT:
		m_translateX += 1/m_scale;
		break;
	case GLUT_KEY_RIGHT:
		m_translateX -= 1/m_scale;
		break;
	default:
		break;
	}

	glutPostRedisplay();
}

void drawObject(DrawObj *obj)
{
	if(obj->GetType() == CPGEOM_LINE)
	{
		Line *line = (Line *)obj;
		glVertex3f(line->GetStartPoint().GetX(), line->GetStartPoint().GetY(), -15.0f);
		glVertex3f(line->GetEndPoint().GetX(), line->GetEndPoint().GetY(), -15.0f);
	}
	else if(obj->GetType() == CPGEOM_ARC)
	{
		Arc *arc = (Arc *)obj;
		LineArray *lines = arc->GetLineSegArray(5);

		for (int i = 0; i < lines->size(); i++)
		{
			Line *line = lines->at(i);
			glVertex3f(line->GetStartPoint().GetX(), line->GetStartPoint().GetY(), -15.0f);
			glVertex3f(line->GetEndPoint().GetX(), line->GetEndPoint().GetY(), -15.0f);
		}
	}
}

void getIntersectionPoint(DrawObj *obj, double y, vector<double> &pointX /*output*/)
{
	if(obj->GetType() == CPGEOM_LINE)
	{
		Line *line = (Line *)obj;

		PointD startPoint = line->GetStartPoint();
		PointD endPoint = line->GetEndPoint();
		if ((y >= endPoint.GetY() && y <= startPoint.GetY()) || (y <= endPoint.GetY() && y >= startPoint.GetY())) 
		{
			if (abs(endPoint.GetY() - startPoint.GetY()) > 1e-10)
			{
				double x = (endPoint.GetX() + (y-endPoint.GetY())/(startPoint.GetY() - endPoint.GetY())*(startPoint.GetX() - endPoint.GetX())); 
				pointX.push_back(x);
			}
		}
	}
	else if(obj->GetType() == CPGEOM_ARC)
	{
		Arc *arc = (Arc *)obj;
		RectD rect = arc->GetRect();
		if ((y >= rect.GetBottom() && y <= rect.GetTop())) 
		{
			// Find Intersection point between horizontal line and circle
			double x1 = arc->GetCenter().GetX() + sqrt(pow(arc->GetRadius(),2) - pow(y-arc->GetCenter().GetY(),2));
			double x2 = arc->GetCenter().GetX() - sqrt(pow(arc->GetRadius(),2) - pow(y-arc->GetCenter().GetY(),2));

			// Check the candidate points are on the arc or not.
			if (arc->Is_Pt_On_Arc(PointD(x1,y)))
			{
				pointX.push_back(x1);
			}
			if (arc->Is_Pt_On_Arc(PointD(x2,y)))
			{
				pointX.push_back(x2);
			}
		}
	}
}

void displayFunc(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glScalef(m_scale, m_scale, 1.0f );
	glTranslatef(m_translateX-center_objectX, m_translateY-center_objectY, 1.0f);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	//// Draw polygon only
	//glColor3f(1.0, 1.0,1.0);

	//if (g_polygon != NULL)
	//{
	//	glColor3f(205.0/255.0, 192.0/255.0, 176.0/255.0);
	//	glBegin(GL_LINES);
	//	DrawObjArray *boundary = g_polygon->getBoundary();
	//	for (int i = 0; i < boundary->size(); i++)
	//	{
	//		drawObject(boundary->at(i));
	//	}
	//	glEnd();

	//	PolygonArray *cuttingArea = g_polygon->getCuttingArea();
	//	for (int i = 0; i < cuttingArea->size(); i++)
	//	{
	//		glBegin(GL_LINES);
	//		DrawObjArray *hole = cuttingArea->at(i)->getBoundary();
	//		for (int j = 0; j < hole->size(); j++)
	//		{
	//			drawObject(hole->at(j));
	//		}
	//		glEnd();
	//	}
	//}
	
	double x1, y1;
	double x2, y2;
	double x3, y3;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	// Draw triangle
	glColor3f(205.0/255.0, 192.0/255.0, 176.0/255.0);
	glBegin(GL_TRIANGLES);
	for (unsigned int i = 0; i < g_triangles.size(); i++)
	{
		x1 = g_vertices[g_triangles[i]->apexes[0]]->point.x;
		y1 = g_vertices[g_triangles[i]->apexes[0]]->point.y;
		x2 = g_vertices[g_triangles[i]->apexes[1]]->point.x;
		y2 = g_vertices[g_triangles[i]->apexes[1]]->point.y;
		x3 = g_vertices[g_triangles[i]->apexes[2]]->point.x;
		y3 = g_vertices[g_triangles[i]->apexes[2]]->point.y;

		glVertex3f(x1, y1, -15.0f);
		glVertex3f(x2, y2, -15.0f);
		glVertex3f(x3, y3, -15.0f);
	}
	glEnd();

	// Draw scanlines only
	glColor3f(1.0, 1.0,0.0);
	glBegin(GL_LINES);
	if (g_scanlines != NULL)
	{
		for (int i = 0; i < g_scanlines->size(); i++)
		{
			drawObject(g_scanlines->at(i));
		}
	}
	glEnd();

	glutSwapBuffers();
}

bool IsOnBoundary(PointD Pt, const DrawObjArray* const pBoundaryObjs)
{
	int i;

	int nObjType;
	DrawObj* pCurrentObj;
	Line* pLineSeg;
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
				return true;
			}
		}
		else if (CPGEOM_ARC == nObjType)
		{
			pArc = (Arc*)pCurrentObj;
			if (pArc->Is_Pt_On_Arc(Pt))
			{
				return true;
			}	
		}
	}

	return false;
}

DrawObjArray * Rasterization(Polygon *polygons, double offset)
{
	DrawObjArray *scanLines = new DrawObjArray;
	vector<double> pointX;

	// Step 1: Find the bottom point and the top point of polygons
	polygons->SetRect();
	double bottom_point_y = polygons->GetRect().GetBottom();
	double top_point_y = polygons->GetRect().GetTop();

	// Step 2: make gf_polygon
	gf_polygon *gf_poly = polygons->Make_GFPolygon();


	// Step 3: Loop from the bottom to the top.
	for (double pointY=bottom_point_y+offset; pointY<top_point_y; pointY=pointY+offset) 
	{
		// Step 3.1: For each scan line, Find all intersection points between the scan line and polygon
		pointX.resize(0);

		// Step 3.1.1: For each scan line, Find all intersection points between the scan line and boundary polygon
		DrawObjArray *boundary = g_polygon->getBoundary();
		for (int i = 0; i < boundary->size(); i++)
		{
			getIntersectionPoint(boundary->at(i), pointY, pointX);
		}

		// Step 3.1.2: For each scan line, Find all intersection points between the scan line and cutting areas (holes)
		PolygonArray *cuttingArea = g_polygon->getCuttingArea();
		for (int i = 0; i < cuttingArea->size(); i++)
		{
			DrawObjArray *hole = cuttingArea->at(i)->getBoundary();
			for (int j = 0; j < hole->size(); j++)
			{
				getIntersectionPoint(hole->at(j), pointY, pointX);
			}
		}

		//  Step 3.2: For each scan-line, sorting all intersection points from the left to the right.
		bool swapped = true;
		double tmp;
		while(swapped)
		{
			swapped = false;
			for (int i = 0; i < (int) pointX.size()-1; i++)
			{
				if (pointX[i] > pointX[i+1])
				{
					tmp = pointX[i];
					pointX[i]=pointX[i+1];
					pointX[i+1]=tmp;
					swapped = true;
				}
			}
		}

		// Step 3.3: For each two adjacent points in scan line, check if they can connect each other.
		for (int i = 0; i < (int)pointX.size()-1; i++)
		{
			if (abs(pointX[i] - pointX[i+1]) > 1e-10)
			{
				bool on_boundary_flag = false;
				double middle_point;

				// 3.3.1: Find the middle point between two adjacent points
				middle_point = (pointX[i] + pointX[i+1])/2;

				// 3.3.2: Check the middle point on the boundary or not
				DrawObjArray *boundary = g_polygon->getBoundary();
				on_boundary_flag = on_boundary_flag || IsOnBoundary(PointD(middle_point,pointY), boundary);

				if (on_boundary_flag)
				{
					continue;
				}

				// 3.3.3: Check the middle point on the cutting areas (holes) or not
				PolygonArray *cuttingArea = g_polygon->getCuttingArea();
				for (int j = 0; j < cuttingArea->size(); j++)
				{
					DrawObjArray *hole = cuttingArea->at(j)->getBoundary();
					on_boundary_flag = on_boundary_flag || IsOnBoundary(PointD(middle_point,pointY), hole);

					if (on_boundary_flag)
					{
						break;
					}
				}

				// 3.3.4: If the middle point is inside the polygon, add the line between two adjacent points to the scan-lines
				if (on_boundary_flag == false && PointInPolygon(middle_point, pointY, gf_poly))
				{
					Line *line = new Line;
					line->SetStartPoint(PointD(pointX[i], pointY));
					line->SetEndPoint(PointD(pointX[i+1], pointY));
					line->SetRect();

					//add line
					scanLines->push_back(line);
				}
			}
		}
	}

	gfFreePolygon(gf_poly);

	// return here.
	return scanLines;
}

Polygon *ReadPolygonFromFile(char *filename)
{
	int size = 1024, pos;
	int c;
	char *buffer = (char *)malloc(size);
	FILE *f = fopen(filename, "r");

	Polygon *polygon = new Polygon;

	bool void_flag = false;
	if(f) {
		do { // read all lines in file
			pos = 0;
			do{ // read one line
				c = fgetc(f);
				if(c != EOF) buffer[pos++] = (char)c;
				if(pos >= size - 1) { // increase buffer length - leave room for 0
					size *=2;
					buffer = (char*)realloc(buffer, size);
				}
			}while(c != EOF && c != '\n');
			buffer[pos] = 0;

			// line is now in buffer
			if (strstr(buffer, "POLYGON"))
			{
				if (void_flag == true)
				{
					if (polygon->getCuttingArea() == NULL)
					{
						polygon->setCuttingArea(new PolygonArray);
					}
					
					polygon->getCuttingArea()->push_back(new Polygon);

					int size = polygon->getCuttingArea()->size();

					polygon->getCuttingArea()->at(size-1)->setBoundary(new DrawObjArray);
				}
				else
				{
					polygon->setBoundary(new DrawObjArray);
				}
			}
			else if (!strstr(buffer, "VOID"))
			{
				char *result[11];
				char * temp;
				int cur_idx = 0;
				temp = strtok(buffer," )(");

				while (temp != NULL)
				{
					result[cur_idx] = temp;
					cur_idx++;
					temp = strtok(NULL," )(\n\r");
				}

				if (cur_idx > 4)
				{
					DrawObj *drawObj = NULL;
					if (strcmp(result[0], "DL") == 0)
					{
						double startX;
						double startY;

						double endX;
						double endY;

						sscanf(result[1],"%lf",&startX);
						sscanf(result[2],"%lf",&startY);

						sscanf(result[3],"%lf",&endX);
						sscanf(result[4],"%lf",&endY);

						Line *line = new Line;
						line->SetStartPoint(PointD(startX, startY));
						line->SetEndPoint(PointD(endX, endY));
						line->SetRect();

						drawObj = line;
					}
					else
					{
						double startX;
						double startY;
						double endX;
						double endY;
						double centerX;
						double centerY;
						double radius;
						double startDegree;
						double endDegree;

						sscanf(result[1],"%lf",&radius);
						sscanf(result[2],"%lf",&startDegree);
						sscanf(result[3],"%lf",&endDegree);
						sscanf(result[5],"%lf",&startX);
						sscanf(result[6],"%lf",&startY);
						sscanf(result[7],"%lf",&endX);
						sscanf(result[8],"%lf",&endY);
						sscanf(result[9],"%lf",&centerX);
						sscanf(result[10],"%lf",&centerY);

						Arc *arc = new Arc;
						arc->SetCenter(PointD(centerX, centerY));
						arc->SetStartPoint(PointD(startX, startY));
						arc->SetEndPoint(PointD(endX,endY));
						arc->SetAngle(startDegree,endDegree);
						arc->SetRadius(radius);

						if (strcmp(result[4], "CCW") == 0)
						{
							arc->setDirection(DIRECTION_COUNTER_CLOCKWISE);
						}
						else
						{
							arc->setDirection(DIRECTION_CLOCKWISE);
						}

						arc->SetRect();
						drawObj = arc;
					}
					if (void_flag == true)
					{
						int size = polygon->getCuttingArea()->size();

						polygon->getCuttingArea()->at(size-1)->getBoundary()->push_back(drawObj);
					}
					else
					{
						polygon->getBoundary()->push_back(drawObj);
					}
				}
			}
			else
			{
				void_flag = true;
			}

		} while(c != EOF); 
		fclose(f);           
	}
	free(buffer);

	return polygon;
}

void releaseTriangles(Triangles & triangles)
{
	while (!triangles.empty())
	{
		delete triangles[0];
		triangles.erase(triangles.begin());
	}
	triangles.clear();
}

void releaseEdges(Edges & edges)
{
	while (!edges.empty())
	{
		delete edges[0];
		edges.erase(edges.begin());
	}
	edges.clear();
}
void releaseVertices(Vertices & vertices)
{
	while (!vertices.empty())
	{
		delete vertices[0];
		vertices.erase(vertices.begin());
	}
	vertices.clear();
}

double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

void GetAPointInsidePolygon(double minX, double minY, double maxX, double maxY, gf_contour contour,double& resultX,double &resultY)
{
	int i, j;
	bool c = false;
	double testx , testy;

	while (c == false)
	{
		testx = fRand(minX, maxX);
		testy = fRand(minY, maxY);
		for (i = 0, j = contour.nPoints-1; i < contour.nPoints; j = i++) {
			if ( ((contour.point[i].y>testy) != (contour.point[j].y>testy)) &&
				(testx < (contour.point[j].x-contour.point[i].x) * (testy-contour.point[i].y) / (contour.point[j].y-contour.point[i].y) + contour.point[i].x) )
				c = !c;
		}
	}

	resultX = testx;
	resultY = testy;
}

void GetTriangulate(gf_polygon *poly, Vertices & vertices, Edges & edges, Triangles & triangles)
{
	// Declare input of triangulate function
	struct triangulateio in, out;
	memset(&in, 0, sizeof(triangulateio));
	memset(&out, 0, sizeof(triangulateio));

	// Loop through all vertex of polygon
	for (int i = 0; i < poly->nContours; i++)
	{
		int startVertexIdx = vertices.size();
		double curX,curY;
		int numberOfPoint = 0;
		for (int j = 0; j < poly->contour[i].nPoints; j++)
		{
			if (j == 0)
			{
				curX = poly->contour[i].point[j].x;
				curY = poly->contour[i].point[j].y;

				// Add vertex to vertices list
				Vertex *vertex = new Vertex;
				vertex->point = DoublePoint1(poly->contour[i].point[j].x, poly->contour[i].point[j].y);
				vertices.push_back(vertex);
				numberOfPoint++;
			}
			else
			{
				if (fabs(poly->contour[i].point[j].x - curX) <= 0.001 && 
					fabs(poly->contour[i].point[j].y - curY) <= 0.001)
				{
					continue;
				}

				curX = poly->contour[i].point[j].x;
				curY = poly->contour[i].point[j].y;

				// Add vertex to vertices list
				Vertex *vertex = new Vertex;
				vertex->point = DoublePoint1(poly->contour[i].point[j].x, poly->contour[i].point[j].y);
				vertices.push_back(vertex);
				numberOfPoint++;

				edges.push_back(new Edge(vertices.size()-2, vertices.size()-1));
			}
		}
		edges.push_back(new Edge(vertices.size()-1, startVertexIdx));

		// Number of vertex should be greater than 2
		if (numberOfPoint <= 2)
		{
			releaseVertices(vertices);
			releaseEdges(edges);
			return;
		}
	}

	// Make point list with structure of input's triangulate function
	double *pointlist = new double[vertices.size()*2]();
	for (unsigned int i = 0; i < vertices.size(); i++)
	{
		pointlist[2*i]=vertices[i]->point.x;
		pointlist[2*i+1]=vertices[i]->point.y;
	}

	// Make segment list with structure of input's triangulate function
	int *segmentlist = new int[edges.size()*2]();
	for (unsigned int i = 0; i < edges.size(); i++)
	{
		segmentlist[2*i]=edges[i]->idxA;
		segmentlist[2*i+1]=edges[i]->idxB;
	}

	/// Input points
	in.numberofpoints = vertices.size();
	in.pointlist = pointlist;

	/// Input segments
	in.numberofsegments = edges.size();
	in.segmentlist = segmentlist;

	// Input holes
	double *holes = new double[(poly->nContours-1)*2]();
	in.numberofholes = poly->nContours-1;

	// With each hole, find a point inside it
	double dResultX;
	double dResultY;
	for (int i = 1; i <= in.numberofholes; i++)
	{
		GetAPointInsidePolygon(poly->box[i].ll.x, poly->box[i].ll.y, poly->box[i].ur.x, poly->box[i].ur.y, poly->contour[i],dResultX,dResultY);
		holes[(i-1)*2] = dResultX;
		holes[(i-1)*2+1] = dResultY;
	}
	in.holelist = holes;

	triangulate("pDzen", &in/*input*/, &out/*output*/, 0/*don't care*/);

	// Get vertices list after triangulation
	releaseVertices(vertices);
	int numberofvertices = out.numberofpoints;
	for (int i = 0; i < numberofvertices; ++i)
	{
		Vertex *vertex = new Vertex;
		vertex->point= DoublePoint1(out.pointlist[2*i], out.pointlist[2*i+1]);
		vertices.push_back(vertex);
	}

	// Get edges list after triangulation
	releaseEdges(edges);
	int numberofedges = out.numberofedges;
	for (int i = 0; i < numberofedges; ++i)
	{
		Edge *edge = new Edge(out.edgelist[2*i], out.edgelist[2*i+1]);
		edges.push_back(edge);
	}

	// Get triangles list after triangulation
	int k = 0;
	for (int i = 0; i < out.numberoftriangles; i++)
	{
		Triangle * triangle = new Triangle;
		vector<int> boudary_vertices;
		for (int j = 0; j < 3; j++)
		{
			triangle->apexes.push_back(out.trianglelist[3*i+j]);

			if (out.neighborlist[3*i+j] != -1)
			{
				triangle->adjacentTriangles.push_back(out.neighborlist[3*i+j]);
			}
		}
		triangles.push_back(triangle);
	}

	// Free memory
	if(in.pointlist)
		free(in.pointlist);
	if(in.pointattributelist)
		free(in.pointattributelist);
	if(in.pointmarkerlist)
		free(in.pointmarkerlist);
	if(in.trianglelist)
		free(in.trianglelist);
	if(in.triangleattributelist)
		free(in.triangleattributelist);
	if(in.trianglearealist)
		free(in.trianglearealist);
	if(in.segmentmarkerlist)
		free(in.segmentmarkerlist);
	if(in.holelist)
		free(in.holelist);
	if(in.segmentlist)
		free(in.segmentlist);

	if(out.edgelist)
		free(out.edgelist);
	if(out.pointattributelist)
		free(out.pointattributelist);
	if(out.pointmarkerlist)
		free(out.pointmarkerlist);
	if(out.trianglelist)
		free(out.trianglelist);
	if(out.triangleattributelist)
		free(out.triangleattributelist);
	if(out.neighborlist)
		free(out.neighborlist);
	if(out.segmentlist)
		free(out.segmentlist);
	if(out.segmentmarkerlist)
		free(out.segmentmarkerlist);
	if(out.regionlist)
		free(out.regionlist);
	if(out.edgemarkerlist)
		free(out.edgemarkerlist);
	if(out.normlist)
		free(out.normlist);
}

int main(int argc, char** argv)
{
	const char* windowName = "Raster";
	glutInit(&argc, argv);

	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(window_width, window_height);
	glutCreateWindow(windowName);
	glutDisplayFunc(displayFunc);
	glutSpecialFunc(mykey);
	initializeOpengl(window_width, window_height);

	

	// Read polygon from file or create a polygon
	g_polygon = ReadPolygonFromFile("sample6.txt");

	gf_polygon *gf_poly = g_polygon->Make_GFPolygon();
	// Step 2: Make triangulate
	GetTriangulate(gf_poly/*input*/, g_vertices/*output*/, g_edges/*output*/, g_triangles/*output*/);

	gfFreePolygon(gf_poly);

	time_t time_start = clock();

	// Offset between 2 scanlines
	double offset = 0.3;

	g_scanlines = Rasterization(g_polygon, offset);

	double time_elapsed_total = double(clock() - time_start)/CLOCKS_PER_SEC;
	cout << "\n time_elapsed_total in " << time_elapsed_total << " secs \n";

	glutMainLoop();

	return 0;
}