#include "Mesh.h"


Mesh::Mesh(Node node_list[], int n)
{
	for(int i=0;i<n;i++)
		nodes.push(node_list+i);
	max_index = n;

	globalStiffness = (float**)malloc(sizeof(float*) * n * 3);
	globalMass = (float**)malloc(sizeof(float*) * n * 3);
	for(int i=0;i<n;i++)
	{
		globalStiffness[i] = (float*)malloc(sizeof(float) * n * 3);
		globalMass[i] = (float*)malloc(sizeof(float) * n * 3);
	}

	resetGlobalStiffness();
	resetGlobalMass();
}

void
Mesh::assembleMassMat()
{
	for(int i=0;i<max_index;i++)
		for(int j=0;j<max_index;j++)
		{
		}
}

void
Mesh::resetGlobalStiffness()
{
	for(int i=0;i<max_index;i++)
		for(int j=0;j<max_index;j++)
			globalStiffness[i][j] = 0.0;
}

void
Mesh::resetGlobalMass()
{
	for(int i=0;i<max_index;i++)
		for(int j=0;j<max_index;j++)
			globalMass[i][j] = 0.0;
}

void
Mesh::addElement(int ind0, int ind1, int ind2, int ind3, float _E, float _v, float density)
{
	elements.push(new Element(nodes[ind0], nodes[ind1], nodes[ind2], nodes[ind3], _E, _v, density));
}

Mesh::~Mesh(void)
{
}
