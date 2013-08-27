#include "QuadTetElement.h"


QuadTetElement::QuadTetElement(Node* corner1, Node* corner2, Node* corner3, Node* corner4,
							   Node* mid12, Node* mid23, Node* mid31, Node* mid41,
							   Node* mid42, Node* mid43,
							   float _E, float _v, float _density)
{
		E = _E;
		v = _v;
		density = _density;
		nodes[0] = corner1;
		nodes[1] = corner2;
		nodes[2] = corner3;
		nodes[3] = corner4;
		nodes[4] = mid12;
		nodes[5] = mid23;
		nodes[6] = mid31;
		nodes[7] = mid41;
		nodes[8] = mid42;
		nodes[9] = mid43;
}

#define JAC(d,x,y) (J[x][d] - J[y][d])

void
QuadTetElement::computeB(float s[4])
{
	// Jx1 Jy1 Jz1
	// Jx2 Jy2 Jz2
	// Jx3 Jy3 Jz3
	// Jx4 Jy4 Jz4

	float x[10][3];

	for(int i=0;i<10;i++)
	{
		x[i][0] = nodes[i]->pos.x;
		x[i][1] = nodes[i]->pos.y;
		x[i][2] = nodes[i]->pos.z;
	}

	Matrix4d J = Matrix4d(0.25, 0.25, 0.25, 0.25,
						  x[0][0] * (s[0] - 0.25) + x[4][0] * s[1] + x[6][0] * s[2] + x[7][0] * s[3],   x[4][0] * s[0] + x[1][0] * (s[1] - 0.25) + x[5][0] * s[2] + x[0][0] * s[3],   x[6][0] * s[0] + x[5][0] * s[1] + x[2][0] * (s[2] - 0.25) + x[9][0] * s[3],   x[7][0] * s[0] + x[8][0] * s[1] + x[9][0] * s[2] + x[3][0] * (s[3] - 0.25),
						  x[0][1] * (s[0] - 0.25) + x[4][1] * s[1] + x[6][1] * s[2] + x[7][1] * s[3],   x[4][1] * s[0] + x[1][1] * (s[1] - 0.25) + x[5][1] * s[2] + x[0][1] * s[3],   x[6][1] * s[0] + x[5][1] * s[1] + x[2][1] * (s[2] - 0.25) + x[9][1] * s[3],   x[7][1] * s[0] + x[8][1] * s[1] + x[9][1] * s[2] + x[3][1] * (s[3] - 0.25),
						  x[0][2] * (s[0] - 0.25) + x[4][2] * s[1] + x[6][2] * s[2] + x[7][2] * s[3],   x[4][2] * s[0] + x[1][2] * (s[1] - 0.25) + x[5][2] * s[2] + x[0][2] * s[3],   x[6][2] * s[0] + x[5][2] * s[1] + x[2][2] * (s[2] - 0.25) + x[9][2] * s[3],   x[7][2] * s[0] + x[8][2] * s[1] + x[9][2] * s[2] + x[3][2] * (s[3] - 0.25)) * 4;
	
	Matrix4d Jinv = J.inverse();
	
}

QuadTetElement::~QuadTetElement(void)
{
}
