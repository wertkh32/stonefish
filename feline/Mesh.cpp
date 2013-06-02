#include "Mesh.h"


Mesh::Mesh(Node node_list[], int n)
{
	for(int i=0;i<n;i++)
		nodes.push(node_list+i);
	numnodes = n;

	globalStiffness = (float**)malloc(sizeof(float*) * n * 3);
	globalMass = (float**)malloc(sizeof(float*) * n * 3);
	for(int i=0;i<n * 3;i++)
	{
		globalStiffness[i] = (float*)malloc(sizeof(float) * n * 3);
		globalMass[i] = (float*)malloc(sizeof(float) * n * 3);
	}

	resetGlobalStiffness();
	resetGlobalMass();
}

float**
Mesh::assembleGlobalMass()
{
	resetGlobalMass();

	for(int i=0;i<elements.size();i++)
	{
		GenMatrix<float,12,12>* masse = elements[i]->getMassMat();

		for(int a=0;a<4;a++)
			for(int b=0;b<4;b++)
				for(int j=0;j<3;j++)
					for(int k=0;k<3;k++)
					{
						globalMass[nodeIndices[i][a] * 3 + j][nodeIndices[i][b] * 3 + k] += (*masse)(a + j, b + k);
					}
			
	}

	return globalMass;
}

float**
Mesh::assembleGlobalStiffness()
{
	resetGlobalStiffness();

	for(int i=0;i<elements.size();i++)
	{
		GenMatrix<float,12,12>* stiffe = elements[i]->getStiffnessMat();

		for(int a=0;a<4;a++)
			for(int b=0;b<4;b++)
				for(int j=0;j<3;j++)
					for(int k=0;k<3;k++)
					{
						globalStiffness[nodeIndices[i][a] * 3+j][nodeIndices[i][b] * 3+k] += (*stiffe)(a + j, b + k);
						//printf("%lf ",(*stiffe)(a+j,b+k));
					}
			
	}

	return globalStiffness;
}

void
Mesh::resetGlobalStiffness()
{
	for(int i=0;i<numnodes * 3;i++)
		for(int j=0;j<numnodes * 3;j++)
			globalStiffness[i][j] = 0.0;
}

void
Mesh::resetGlobalMass()
{
	for(int i=0;i<numnodes * 3;i++)
		for(int j=0;j<numnodes * 3;j++)
			globalMass[i][j] = 0.0;
}

void
Mesh::addElement(int ind0, int ind1, int ind2, int ind3, float _E, float _v, float density)
{
	nodeIndices[elements.size()][0] = ind0;
	nodeIndices[elements.size()][1] = ind1;
	nodeIndices[elements.size()][2] = ind2;
	nodeIndices[elements.size()][3] = ind3;
	elements.push(new Element(nodes[ind0], nodes[ind1], nodes[ind2], nodes[ind3], _E, _v, density));
}

Mesh::~Mesh(void)
{
}
