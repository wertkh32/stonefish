#pragma once
#include "includes.h"
#include "Element.h"
#include "QuickArray.h"

class Mesh
{
	int numnodes;
public:
	QuickArray<Element*,MAX_ELEMENTS> elements;
	int nodeIndices[MAX_ELEMENTS][4];
	QuickArray<Node*,MAX_NODES> nodes;
	double **globalStiffness, **globalMass;

	int getNoOfElements(){return elements.size();}
	int getNoOfNodes(){return numnodes;}
	
	void addElement(int ind0, int ind1, int ind2, int ind3, double _E, double _v, double density);
	void resetGlobalStiffness();
	void resetGlobalMass();
	double** assembleGlobalMass();
	double** assembleGlobalStiffness();
	Mesh(Node nodes[], int n);
	~Mesh(void);
};

