#pragma once
#include "includes.h"
#include "TetElement.h"
#include "QuickArray.h"
#include "PolarDecompose.h"

class Mesh
{
	int numnodes;

	//////////////////////////////////////////////
public:
	QuickArray<TetElement*,MAX_ELEMENTS> elements;
	int nodeIndices[MAX_ELEMENTS][4];
	QuickArray<Node*,MAX_NODES> nodes;
	float **globalStiffness, **globalMass;
	float dt;

	int getNoOfElements(){return elements.size();}
	int getNoOfNodes(){return numnodes;}
	
	void addElement(int ind[4], float _E, float _v, float density);
	void resetGlobalStiffness();
	void resetGlobalMass();
	float** assembleGlobalMass();
	float** assembleGlobalStiffness();

	Mesh(Node nodes[], int n);
	~Mesh(void);
};

