#include "StdAfx.h"
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

void
QuadTetElement::computeB(float point[4])
{
	// Jx1 Jy1 Jz1
	// Jx2 Jy2 Jz2
	// Jx3 Jy3 Jz3
	// Jx4 Jy4 Jz4

	float J[4][3];
	float x[10][3];

	for(int i=0;i<10;i++)
	{
		x[i][0] = nodes[i]->pos.x;
		x[i][1] = nodes[i]->pos.y;
		x[i][2] = nodes[i]->pos.z;
	}

	J[0][0] = 4 * (x[0][0] * (point[0] - 0.25) + x[4][0] * point[1] + x[6][0] * point[2] + x[7][0] * point[3]);
	J[0][1] = 4 * (x[0][1] * (point[0] - 0.25) + x[4][1] * point[1] + x[6][1] * point[2] + x[7][1] * point[3]);
	J[0][2] = 4 * (x[0][2] * (point[0] - 0.25) + x[4][2] * point[1] + x[6][2] * point[2] + x[7][2] * point[3]);

	J[1][0] = 4 * (x[4][0] * point[0] + x[1][0] * (point[1] - 0.25) + x[5][0] * point[2] + x[8][0] * point[3]);
	J[1][1] = 4 * (x[4][1] * point[0] + x[1][1] * (point[1] - 0.25) + x[5][1] * point[2] + x[8][1] * point[3]);
	J[1][2] = 4 * (x[4][2] * point[0] + x[1][2] * (point[1] - 0.25) + x[5][2] * point[2] + x[8][2] * point[3]);

	J[2][0] = 4 * (x[6][0] * point[0] + x[5][0] * point[1] + x[2][0] * (point[2] - 0.25) + x[9][0] * point[3]);
	J[2][1] = 4 * (x[6][1] * point[0] + x[5][1] * point[1] + x[2][1] * (point[2] - 0.25) + x[9][1] * point[3]);
	J[2][2] = 4 * (x[6][2] * point[0] + x[5][2] * point[1] + x[2][2] * (point[2] - 0.25) + x[9][2] * point[3]);

	J[3][0] = 4 * (x[7][0] * point[0] + x[8][0] * point[1] + x[9][0] * point[2] + x[3][0] * (point[3] - 0.25) );
	J[3][1] = 4 * (x[7][1] * point[0] + x[8][1] * point[1] + x[9][1] * point[2] + x[3][1] * (point[3] - 0.25) );
	J[3][2] = 4 * (x[7][2] * point[0] + x[8][2] * point[1] + x[9][2] * point[2] + x[3][2] * (point[3] - 0.25) );

	

}

QuadTetElement::~QuadTetElement(void)
{
}
