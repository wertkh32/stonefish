#include "Element.h"


Element::Element(Node* n1, Node* n2, Node* n3, Node* n4)
{
	
}

void Element::preCompute()
{
	undeformShapeMat = 
			Matrix3d(nodes[0]->pos.x - nodes[3]->pos.x,nodes[1]->pos.x - nodes[3]->pos.x,nodes[2]->pos.x - nodes[3]->pos.x,
					nodes[0]->pos.y - nodes[3]->pos.y,nodes[1]->pos.y - nodes[3]->pos.y,nodes[2]->pos.y - nodes[3]->pos.y,
					nodes[0]->pos.z - nodes[3]->pos.z,nodes[1]->pos.z - nodes[3]->pos.z,nodes[2]->pos.z - nodes[3]->pos.z);
	undeformShapeMatInv = undeformShapeMat.inverse();
	undeformVolume = (1.0/6.0) * undeformShapeMat.determinant();
}

void Element::preComputeShapeFuncDeriv()
{
	
}

Matrix3d Element::computeDeformationMat()
{
	Matrix3d deformShapeMat = 
			Matrix3d(nodes[0]->pos_t.x - nodes[3]->pos_t.x,nodes[1]->pos_t.x - nodes[3]->pos_t.x,nodes[2]->pos_t.x - nodes[3]->pos_t.x,
					nodes[0]->pos_t.y - nodes[3]->pos_t.y,nodes[1]->pos_t.y - nodes[3]->pos_t.y,nodes[2]->pos_t.y - nodes[3]->pos_t.y,
					nodes[0]->pos_t.z - nodes[3]->pos_t.z,nodes[1]->pos_t.z - nodes[3]->pos_t.z,nodes[2]->pos_t.z - nodes[3]->pos_t.z);
	return deformShapeMat * undeformShapeMatInv;
}


Element::~Element(void)
{
}
