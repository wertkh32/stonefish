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

	GenMatrix<float,6,12> strainMat;
	GenMatrix<float,6,6> matConstantsMat;
	GenMatrix<float,12,12> undeformStiffnessMat,
						   massMat, RK, RKRT;
	float undeformVolume;
	float E, v; // E is Young's Modulus, v is Poisson's Ratio
	float density;
	void preCompute();
	void preComputeUndeformedStiffnessMat();
	void preComputeMassMat();

public:
	Element(Node* n1, Node* n2, Node* n3, Node* n4, float _E, float _v, float _density);
	Matrix3d computeDeformationMat();
	Matrix3d computeDeformationMatDeriv();
	Matrix3d computeDeformShapeMat();
	Matrix3d computeDeformShapeMatDeriv();


	Matrix3d& getUndeformShapeMat(){return undeformShapeMat;}
	Matrix3d& getUndeformShapeMatInv(){return undeformShapeMatInv;}
	Node** getNodes(){return nodes;}
	//computeundeformstiffnessmat B^T * c * B where B = LN
	//computestiffnessmat dfe/dx |x = x*

	GenMatrix<float,12,12>* getStiffnessMat(){return &undeformStiffnessMat;}
	GenMatrix<float,12,12>* getMassMat(){return &massMat;}
	void getRKRTandRK(GenMatrix<float,12,12>& RK, GenMatrix<float,12,12>& RKRT);

	void renderElement();

	~Element(void);
};

