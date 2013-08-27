#pragma once
#include "Node.h"
#include "Matrix4d.h"
#include "Matrix3d.h"

class QuadTetElement
{
public:
	Node* nodes[10];
	float E, v, density, dt;
	float b[4][10];
	
	QuadTetElement(Node* corner1, Node* corner2, Node* corner3, Node* corner4,
				   Node* mid12, Node* mid23, Node* mid31, Node* mid41,
				   Node* mid42, Node* mid43,
				   float _E, float _v, float _density);

	void computeB(float point[4]);


	~QuadTetElement(void);
};

