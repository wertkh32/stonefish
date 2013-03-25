#include "Mesh.h"


Mesh::Mesh(Node node_list[], int n)
{
	for(int i=0;i<n;i++)
		nodes.push(node_list+i);
	max_index = n;
}

void
Mesh::addElement(int ind0, int ind1, int ind2, int ind3, float _E, float _v, float density)
{
	elements.push(new Element(nodes[ind0], nodes[ind1], nodes[ind2], nodes[ind3], _E, _v, density));
}

Mesh::~Mesh(void)
{
}
