#pragma once

#include "includes.h"
#include "Mesh.h"
#include "QuadTetMesh.h"
#include "QuickArray.h"

#define MAX_GEO 3000
//assume all faces are triangles.
//each vertex incident with 3 other vertices.
//just vertices and faces for now

enum FACE_TYPE
{
	TRIANGLE=0, QUAD=1
};

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
	vector3<float> coords;
	vector3<float> norm;
	//edge* edges[3];
};

struct face
{
	FACE_TYPE type;
	int vindex[4];
	//edge* edges[3];
};

struct bary
{
	int element_no;
	vector3<float> barycoords;
};

typedef QuickArray<vertex,MAX_GEO> vertArray;
typedef QuickArray<edge,MAX_GEO> edgeArray;
typedef QuickArray<face,MAX_GEO> faceArray;
typedef QuickArray<bary,MAX_GEO> baryArray;

class Model
{
	float fact(int n);
public:
	vertArray verts;
	edgeArray edges;
	faceArray faces;
	baryArray barys;

	MESH* mesh;

	void computeBarycentricCoords();

	Model(void (*Modelfunc)(vertArray*,	edgeArray*,	faceArray*), void (*Meshfunc)(MESH**));
	Model(void (*Modelfunc)(vertArray*,	edgeArray*,	faceArray*), MESH* mesh);
	void interpolateVerts();
	void render();
	~Model(void);
};

