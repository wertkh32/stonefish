#pragma once
#include "defines.h"
#include "QuickArray.h"
#include "QuadTetElement.h"


class QuadTetMesh
{
public:
	int numnodes;

	QuickArray<QuadTetElement*,MAX_ELEMENTS> elements;
	int nodeIndices[MAX_ELEMENTS][10];
	QuickArray<Node*,MAX_NODES> nodes;
	float dt;

	int getNoOfElements(){return elements.size();}
	int getNoOfNodes(){return numnodes;}

	void addElement(int ind[10], float _E, float _v, float density);

	QuadTetMesh(Node node_list[], int n);
	~QuadTetMesh(void);
};

