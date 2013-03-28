#pragma once
#include "includes.h"
#include "Node.h"
#include "Matrix3d.h"
#include "Matrix4d.h"
#include "GenMatrix.h"
#include "PolarDecompose.h"

class Element
{
	//implements a tetrahedra element
	Node* nodes[4];
	Matrix3d stiffnessMat;
	Matrix3d undeformShapeMat,
		     undeformShapeMatInv;

	GenMatrix<double,6,12> strainMat;
	GenMatrix<double,6,6> matConstantsMat;
	GenMatrix<double,12,12> undeformStiffnessMat,
						   massMat, RK, RKRT;
	double undeformVolume;
	double E, v; // E is Young's Modulus, v is Poisson's Ratio
	double density;
	void preCompute();
	void preComputeUndeformedStiffnessMat();
	void preComputeMassMat();

public:
	Element(Node* n1, Node* n2, Node* n3, Node* n4, double _E, double _v, double _density);
	Matrix3d computeDeformationMat();

	//computeundeformstiffnessmat B^T * c * B where B = LN
	//computestiffnessmat dfe/dx |x = x*

	GenMatrix<double,12,12>* getStiffnessMat(){return &undeformStiffnessMat;}
	GenMatrix<double,12,12>* getMassMat(){return &massMat;}
	void getRKRTandRK(GenMatrix<double,12,12>& RK, GenMatrix<double,12,12>& RKRT);

	void renderElement();

	~Element(void);
};

