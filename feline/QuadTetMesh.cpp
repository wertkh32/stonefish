#include "QuadTetMesh.h"


QuadTetMesh::QuadTetMesh(Node node_list[], int n)
{
	for(int i=0;i<n;i++)
		nodes.push(node_list+i);
	numnodes = n;

	dt = 1.0/FPS;
}

void
QuadTetMesh::addElement(int ind[10], float _E, float _v, float density)
{
	for(int i=0;i<10;i++)
		nodeIndices[elements.size()][i] = ind[i];

	Node* elenodes[10];

	for(int i=0;i<10;i++)
		elenodes[i] = nodes[ind[i]];

	elements.push(new QuadTetElement(elenodes, _E, _v, density));
}


QuadTetMesh::~QuadTetMesh(void)
{
}
