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
	//inv is the inverse of the matrix relating volume coords to cartesian coords
	Matrix4d inv =	  Matrix4d
					  ( 1.0, 1.0, 1.0, 1.0,
						nodes[0]->pos.x, nodes[1]->pos.x, nodes[2]->pos.x, nodes[3]->pos.x,
						nodes[0]->pos.y, nodes[1]->pos.y, nodes[2]->pos.y, nodes[3]->pos.y,
						nodes[0]->pos.z, nodes[1]->pos.z, nodes[2]->pos.z, nodes[3]->pos.z).inverse();

	//strain matrix B = LN = dN/dx
	float strainMatrix[6][12] = 
	{
		{ inv(0,1), 0, 0, inv(1,1), 0, 0, inv(2,1), 0, 0, inv(3,1), 0, 0 },
		{ 0, inv(0,2), 0, 0, inv(1,2), 0, 0, inv(2,2), 0, 0, inv(3,2), 0 },
		{ 0, 0, inv(0,3), 0, 0, inv(1,3), 0, 0, inv(2,3), 0, 0, inv(3,3) },
		{ inv(0,2), inv(0,1), 0, inv(1,2), inv(1,1), 0, inv(2,2), inv(2,1), 0, inv(3,2), inv(3,1), 0 },
		{ 0, inv(0,3), inv(0,2), 0, inv(1,3), inv(1,2), 0, inv(2,3), inv(2,2), 0, inv(3,3), inv(3,2) },
		{ inv(0,3), 0, inv(0,1), inv(1,3), 0, inv(1,1), inv(2,3), 0, inv(2,1), inv(3,3), 0, inv(3,1) }
	};

	strainMat = GenMatrix<float,6,12>(strainMatrix);

	//material constants
	float C[6][6] = 
	{
	};
}

Matrix3d Element::computeDeformationMat()
{
	Matrix3d deformShapeMat 
				   (nodes[0]->pos_t.x - nodes[3]->pos_t.x,nodes[1]->pos_t.x - nodes[3]->pos_t.x,nodes[2]->pos_t.x - nodes[3]->pos_t.x,
					nodes[0]->pos_t.y - nodes[3]->pos_t.y,nodes[1]->pos_t.y - nodes[3]->pos_t.y,nodes[2]->pos_t.y - nodes[3]->pos_t.y,
					nodes[0]->pos_t.z - nodes[3]->pos_t.z,nodes[1]->pos_t.z - nodes[3]->pos_t.z,nodes[2]->pos_t.z - nodes[3]->pos_t.z);
	return deformShapeMat * undeformShapeMatInv;
}


Element::~Element(void)
{
}
