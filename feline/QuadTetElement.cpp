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

	//PT
	float P[3][4] = { {Jinv(0,1), Jinv(1,1), Jinv(2,1), Jinv(3,1)},
					  {Jinv(0,2), Jinv(1,2), Jinv(2,2), Jinv(3,2)},
					  {Jinv(0,3), Jinv(1,3), Jinv(2,3), Jinv(3,3)} };


	float dN[4][10] = { {4 * s[0] - 1, 0, 0, 0, 4 * s[1], 0, 4 * s[2], 4 * s[3], 0, 0},

						{0, 4 * s[1] - 1, 0, 0, 4 * s[0], 4 * s[2], 0, 0, 4 * s[3], 0},
						
						{0, 0, 4 * s[2] - 1, 0, 0, 4 * s[1], 4 * s[0], 0, 0, 4 * s[3]},
						
						{0, 0, 0, 4 * s[0] - 1, 0, 0, 0, 4 * s[0], 4 * s[1], 4 * s[2]} };

	float dNdX[3][10];

	//PT * dN
	for(int i=0;i<3;i++)
		for(int j=0;j<10;j++)
		{
			dNdX[i][j] = 0;
			for(int k=0;k<4;k++)
				dNdX[i][j] += P[i][k] * dN[k][j];
		}

		float b[6][30] = { {dNdX[0][0], 0, 0, dNdX[0][1], 0, 0, dNdX[0][2], 0, 0, dNdX[0][3], 0, 0, dNdX[0][4], 0, 0, dNdX[0][5], 0, 0, dNdX[0][6], 0, 0, dNdX[0][7], 0, 0, dNdX[0][8], 0, 0, dNdX[0][9], 0, 0},
						   {0, dNdX[1][0], 0, 0, dNdX[1][1], 0, 0, dNdX[1][2], 0, 0, dNdX[1][3], 0, 0, dNdX[1][4], 0, 0, dNdX[1][5], 0, 0, dNdX[1][6], 0, 0, dNdX[1][7], 0, 0, dNdX[1][8], 0, 0, dNdX[1][9], 0},
						   {0, 0, dNdX[2][0], 0, 0, dNdX[2][1], 0, 0, dNdX[2][2], 0, 0, dNdX[2][3], 0, 0, dNdX[2][4], 0, 0, dNdX[2][5], 0, 0, dNdX[2][6], 0, 0, dNdX[2][7], 0, 0, dNdX[2][8], 0, 0, dNdX[2][9]},
						   {dNdX[1][0], dNdX[0][0], 0, dNdX[1][1], dNdX[0][1], 0, dNdX[1][2], dNdX[0][2], 0, dNdX[1][3], dNdX[0][3], 0, dNdX[1][4], dNdX[0][4], 0, dNdX[1][5], dNdX[0][5], 0, dNdX[1][6], dNdX[0][6], 0, dNdX[1][7], dNdX[0][7], 0, dNdX[1][8], dNdX[0][8], 0, dNdX[1][9], dNdX[0][9], 0},
						   {0, dNdX[2][0], dNdX[1][0], 0, dNdX[2][1], dNdX[1][1], 0, dNdX[2][2], dNdX[1][2], 0, dNdX[2][3], dNdX[1][3], 0, dNdX[2][4], dNdX[1][4], 0, dNdX[2][5], dNdX[1][5], 0, dNdX[2][6], dNdX[1][6], 0, dNdX[2][7], dNdX[1][7], 0, dNdX[2][8], dNdX[1][8], 0, dNdX[2][9], dNdX[1][9]},
						   {dNdX[2][0], 0, dNdX[0][0], dNdX[2][1], 0, dNdX[0][1], dNdX[2][2], 0, dNdX[0][2], dNdX[2][3], 0, dNdX[0][3], dNdX[2][4], 0, dNdX[0][4], dNdX[2][5], 0, dNdX[0][5], dNdX[2][6], 0, dNdX[0][6], dNdX[2][7], 0, dNdX[0][7], dNdX[2][8], 0, dNdX[0][8], dNdX[2][9], 0, dNdX[0][9]} };

		B = GenMatrix<float,6,30>(b);
}

void 
QuadTetElement::computeStiffness()
{
	//TO ADD: 4 point Gaussian Quadrature
	//material constants
	float c1 = (E*(1-v))/((1.0-2.0*v)*(1.0+v)),
		c2 = (E*v)/((1.0-2.0*v)*(1.0+v)),
		c3 = (c1 - c2)/2.0;
	//printf("%f %f %f\n",c1,c2,c3);

	float C[6][6] =
	{
		{ c1, c2, c2, 0, 0, 0 },
		{c2, c1, c2, 0, 0, 0 },
		{c2, c2, c1, 0, 0, 0 },
		{0, 0, 0, c3, 0, 0 },
		{0, 0, 0, 0, c3, 0 },
		{0, 0, 0, 0, 0, c3 }
	};

	GenMatrix<float,6,6> matConstantsMat = GenMatrix<float,6,6>(C);
	
	K = B.transpose() * matConstantsMat * B;

}

QuadTetElement::~QuadTetElement(void)
{
}
