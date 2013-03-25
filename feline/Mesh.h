#pragma once
#include "includes.h"
#include "Element.h"
#include "QuickArray.h"

class Mesh
{
	int max_index;
public:
	QuickArray<Element*,MAX_ELEMENTS> elements;
	int nodeIndices[MAX_ELEMENTS][4];
	QuickArray<Node*,MAX_NODES> nodes;
	float **globalStiffness, **globalMass;


	int getNoOfElements(){return elements.size();}
	int getNoOfNodes(){return max_index;}
	
	void addElement(int ind0, int ind1, int ind2, int ind3, float _E, float _v, float density);
	void resetGlobalStiffness();
	void resetGlobalMass();
	void assembleGlobalMass();
	void assembleGlobalStiffness();

	Mesh(Node nodes[], int n);
	~Mesh(void);
};

