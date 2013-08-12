#pragma once
#include "includes.h"
#include "Node.h"
#include "Matrix3d.h"
#include "Matrix4d.h"
#include "GenMatrix.h"
#include "SymSparseMatrix.h"
#include "PolarDecompose.h"

class Element
{
public:
	//implements a tetrahedra element
	Node* nodes[4];
	Matrix3d stiffnessMat;
	Matrix3d undeformShapeMat,
		     undeformShapeMatInv;


	GenMatrix<float,6,12> strainMat;
	GenMatrix<float,6,6> matConstantsMat;
	GenMatrix<float,12,12> undeformStiffnessMat,
						   massMat, RK, RKRT, A;
	SparseMatrix *sparseStiff;

	float undeformVolume;
	float mass, nodalMass;
	float E, v; // E is Young's Modulus, v is Poisson's Ratio
	float density;
	void preCompute();
	void preComputeUndeformedStiffnessMat();
	void preComputeMassMat();
	float dt;
	//float B[3][3];

	Element(Node* n1, Node* n2, Node* n3, Node* n4, float _E, float _v, float _density);
	Matrix3d computeDeformationMat();
	Matrix3d computeDeformShapeMat();

	Matrix3d& getUndeformShapeMat(){return undeformShapeMat;}
	Matrix3d& getUndeformShapeMatInv(){return undeformShapeMatInv;}
	
	
	Node** getNodes(){return nodes;}
	float getDensity(){return density;}
	float getVolume(){return undeformVolume;}
	float getYoungMod(){return E;}
	float getPoissonRatio(){return v;}
	//computeundeformstiffnessmat B^T * c * B where B = LN
	//computestiffnessmat dfe/dx |x = x*

	GenMatrix<float,12,12>* getStiffnessMat(){return &undeformStiffnessMat;}
	GenMatrix<float,12,12>* getMassMat(){return &massMat;}
	void computeRKRTandRK();
	void getRKRTandRK(GenMatrix<float,12,12>*& RK, GenMatrix<float,12,12>*& RKRT);
	void getRKRTandRK(SparseMatrix& RK, SparseMatrix& RKRT);

	GenMatrix<float,12,12>* getRK(){return &RK;}
	GenMatrix<float,12,12>* getRKRT(){return &RKRT;}
	GenMatrix<float,12,12>* getA(){return &A;}
	
	//mat free vars for timestep
	//Matrix3d undeformShapeMatInvT;
	//Matrix3d Ft, Rt, St, trSI_Sinv;
	//float S_Itrace;
	//float miu, lambda;
	///////////////////////////
	//compute and store for mat products
	//void computeMatFreeVars();
	//Matrix3d& getUndeformShapeMatInvT(){return undeformShapeMatInvT;}
	//Matrix3d computePiolaStressTensor();
	//Matrix3d computeRKtensor(Matrix3d& dF);
	//Matrix3d computeDiffPiolaStressTensor(Matrix3d& dF);

	void renderElement();

	~Element(void);
};

