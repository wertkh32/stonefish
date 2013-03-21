#pragma once
#include "includes.h"
#include "Node.h"
#include "Matrix3d.h"
#include "Matrix4d.h"
#include "GenMatrix.h"

class Element
{
	//implements a tetrahedra element
	Node* nodes[4];
	Matrix3d stiffnessMat;
	Matrix3d undeformShapeMat,
		     undeformShapeMatInv;
	GenMatrix<float,6,12> strainMat;
	GenMatrix<float,6,6> matConstantsMat;
	GenMatrix<float,12,12> undeformStiffnessMat;
	float undeformVolume;
	float E, v; // E is Young's Modulus, v is Poisson's Ratio
	void preCompute();
	void preComputeUndeformedStiffnessMat();

public:
	Element(Node* n1, Node* n2, Node* n3, Node* n4);
	Matrix3d computeDeformationMat();
	//computeundeformstiffnessmat B^T * c * B where B = LN
	//computestiffnessmat dfe/dx |x = x*
	~Element(void);
};

