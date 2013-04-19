#pragma once

#include "includes.h"
#include "Mesh.h"
#include "QuickArray.h"

#define MAX_GEO 1000
//assume all faces are triangles.
//each vertex incident with 3 other vertices.
//just vertices and faces for now

struct edge;
struct vertex;
struct face;

struct edge
{
	vertex *v1, *v2;
	face *f1, *f2;
	edge *prev1, *prev2, *next1, *next2;

};

struct vertex
{
	vector3<double> coords;
	vector3<double> norm;
	//edge* edges[3];
};

struct face
{
	vertex* verts[4];
	//edge* edges[3];
};

struct bary
{
	int element_no;
	vector3<double> barycoords;
};

typedef QuickArray<vertex,MAX_GEO> vertArray;
typedef QuickArray<edge,MAX_GEO> edgeArray;
typedef QuickArray<face,MAX_GEO> faceArray;
typedef QuickArray<bary,MAX_GEO> baryArray;

class Model
{
public:
	vertArray verts;
	edgeArray edges;
	faceArray faces;
	baryArray barys;

	Mesh* mesh;

	void computeBarycentricCoords();

	Model(void (*Modelfunc)(vertArray*,	edgeArray*,	faceArray*), void (*Meshfunc)(Mesh**));
	void interpolateVerts();
	void render();
	~Model(void);
};

