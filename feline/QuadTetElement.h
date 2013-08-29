#pragma once
#include "Node.h"
#include "Matrix4d.h"
#include "Matrix3d.h"
#include "GenMatrix.h"
#include "PolarDecompose.h"

class QuadTetElement
{
public:
	Node* nodes[10];
	float x[10][3];
	float volume;
	float mass;
	float nodemass[10];
	float E, v, density, dt;
	GenMatrix<float, 30, 30> K;
	Matrix3d undeformShapeMatInv, R;
	/*Node* corner1, Node* corner2, Node* corner3, Node* corner4,
				   Node* mid12, Node* mid23, Node* mid31, Node* mid41,
				   Node* mid42, Node* mid43*/

	QuadTetElement(Node* nodess[10],
				   float _E, float _v, float _density);

	void precompute();
	void computeLumpedMasses();
	Matrix3d computeDeformationMat();

	void computeB(float s[4], GenMatrix<float, 6, 30>* B, float* Jdet);
	void computeStiffness();
	void renderElement();

	GenMatrix<float, 30, 30>& getStiffnessMat(){return K;}
	//timestep
	void computeRotation();
	Matrix3d& getRotation(){return R;}

	~QuadTetElement(void);
};