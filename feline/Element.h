#pragma once
#include "includes.h"
#include "Node.h"
#include "Matrix3d.h"

class Element
{
	//implements a tetrahedra element
	Node* nodes[4];
	Matrix3d stiffnessMat;
	Matrix3d undeformShapeMat,
		     undeformShapeMatInv;
	float undeformVolume;
	void preCompute();
	void preComputeShapeFuncDeriv();

public:
	Element(Node* n1, Node* n2, Node* n3, Node* n4);
	Matrix3d computeDeformationMat();
	//computeundeformstiffnessmat B^T * c * B where B = LN
	//computestiffnessmat dfe/dx |x = x*
	~Element(void);
};

