#include "QuadTetElement.h"

//first 4 are corner vertices
//next six are points on the edges
//1,2,3,4, 12, 23, 31, 41, 42, 43
QuadTetElement::QuadTetElement(Node* nodess[10],
							   float _E, float _v, float _density)
{
		E = _E;
		v = _v;
		density = _density;

		for(int i=0;i<10;i++)
			nodes[i] = nodess[i];
	
		for(int i=0;i<10;i++)
		{
			x[i][0] = nodes[i]->pos.x;
			x[i][1] = nodes[i]->pos.y;
			x[i][2] = nodes[i]->pos.z;

			printf("%f %f %f\n",x[i][0],x[i][1],x[i][2]);
		}

		precompute();
}

void QuadTetElement::renderElement()
{

	//1,2,3,4, 5 = 12, 6 = 23, 7 = 31, 8 = 41, 9 = 42, 10 = 43
	glPushMatrix();

	glBegin(GL_LINE_STRIP);

		glVertex3fv(nodes[3]->pos_t.coords);
		glVertex3fv(nodes[7]->pos_t.coords);
		glVertex3fv(nodes[0]->pos_t.coords);

	glEnd();

	glBegin(GL_LINE_STRIP);

		glVertex3fv(nodes[3]->pos_t.coords);
		glVertex3fv(nodes[8]->pos_t.coords);
		glVertex3fv(nodes[1]->pos_t.coords);

	glEnd();

	glBegin(GL_LINE_STRIP);

		glVertex3fv(nodes[3]->pos_t.coords);
		glVertex3fv(nodes[9]->pos_t.coords);
		glVertex3fv(nodes[2]->pos_t.coords);

	glEnd();

	glBegin(GL_LINE_STRIP);

		glVertex3fv(nodes[0]->pos_t.coords);
		glVertex3fv(nodes[4]->pos_t.coords);
		glVertex3fv(nodes[1]->pos_t.coords);

	glEnd();

	glBegin(GL_LINE_STRIP);
		
		glVertex3fv(nodes[1]->pos_t.coords);
		glVertex3fv(nodes[5]->pos_t.coords);
		glVertex3fv(nodes[2]->pos_t.coords);
	
	glEnd();

	glBegin(GL_LINE_STRIP);

		glVertex3fv(nodes[2]->pos_t.coords);
		glVertex3fv(nodes[6]->pos_t.coords);
		glVertex3fv(nodes[0]->pos_t.coords);

	glEnd();


	glPopMatrix();
}

void QuadTetElement::precompute()
{
	//from characteristic polynomial
	undeformShapeMatInv =
		Matrix3d(nodes[0]->pos.x - nodes[3]->pos.x,nodes[1]->pos.x - nodes[3]->pos.x,nodes[2]->pos.x - nodes[3]->pos.x,
				 nodes[0]->pos.y - nodes[3]->pos.y,nodes[1]->pos.y - nodes[3]->pos.y,nodes[2]->pos.y - nodes[3]->pos.y,
				 nodes[0]->pos.z - nodes[3]->pos.z,nodes[1]->pos.z - nodes[3]->pos.z,nodes[2]->pos.z - nodes[3]->pos.z).inverse();

	computeLumpedMasses();
	computeStiffness();

	//for(int i=0;i<30;i++, putchar('\n'))
	//	for(int j=0;j<30;j++)
	//		printf("%f ",K(i,j));
}

void QuadTetElement::computeLumpedMasses()
{
	//from characteristic polynomial
	volume = fabs(Matrix3d(x[0][0] - x[3][0], x[1][0] - x[3][0], x[2][0] - x[3][0],
	        			   x[0][1] - x[3][1], x[1][1] - x[3][1], x[2][1] - x[3][1],
						   x[0][2] - x[3][2], x[1][2] - x[3][2], x[2][2] - x[3][2]).determinant());

	mass = volume * density;

	nodemass[0] = nodemass[1] = nodemass[2] = nodemass[3] = (1.0/32.0) * mass;
	nodemass[4] = nodemass[5] = nodemass[6] = nodemass[7] = nodemass[8] = nodemass[9] = (7.0/48.0) * mass; 
}

Matrix3d QuadTetElement::computeDeformationMat()
{
	Matrix3d deformShapeMat
		(nodes[0]->pos_t.x - nodes[3]->pos_t.x,nodes[1]->pos_t.x - nodes[3]->pos_t.x,nodes[2]->pos_t.x - nodes[3]->pos_t.x,
		nodes[0]->pos_t.y - nodes[3]->pos_t.y,nodes[1]->pos_t.y - nodes[3]->pos_t.y,nodes[2]->pos_t.y - nodes[3]->pos_t.y,
		nodes[0]->pos_t.z - nodes[3]->pos_t.z,nodes[1]->pos_t.z - nodes[3]->pos_t.z,nodes[2]->pos_t.z - nodes[3]->pos_t.z);
	return deformShapeMat * undeformShapeMatInv;
}

void QuadTetElement::computeRotation()
{
	Matrix3d F,S;
	F = computeDeformationMat();
	PolarDecompose::compute(F,R,S);
}

void
QuadTetElement::computeB(float s[4], GenMatrix<float, 6, 30>* B, float* Jdet)
{
	// Jx1 Jy1 Jz1
	// Jx2 Jy2 Jz2
	// Jx3 Jy3 Jz3
	// Jx4 Jy4 Jz4


	Matrix4d J = Matrix4d(0.25, 0.25, 0.25, 0.25,
						  x[0][0] * (s[0] - 0.25) + x[4][0] * s[1] + x[6][0] * s[2] + x[7][0] * s[3],   x[4][0] * s[0] + x[1][0] * (s[1] - 0.25) + x[5][0] * s[2] + x[8][0] * s[3],   x[6][0] * s[0] + x[5][0] * s[1] + x[2][0] * (s[2] - 0.25) + x[9][0] * s[3],   x[7][0] * s[0] + x[8][0] * s[1] + x[9][0] * s[2] + x[3][0] * (s[3] - 0.25),
						  x[0][1] * (s[0] - 0.25) + x[4][1] * s[1] + x[6][1] * s[2] + x[7][1] * s[3],   x[4][1] * s[0] + x[1][1] * (s[1] - 0.25) + x[5][1] * s[2] + x[8][1] * s[3],   x[6][1] * s[0] + x[5][1] * s[1] + x[2][1] * (s[2] - 0.25) + x[9][1] * s[3],   x[7][1] * s[0] + x[8][1] * s[1] + x[9][1] * s[2] + x[3][1] * (s[3] - 0.25),
						  x[0][2] * (s[0] - 0.25) + x[4][2] * s[1] + x[6][2] * s[2] + x[7][2] * s[3],   x[4][2] * s[0] + x[1][2] * (s[1] - 0.25) + x[5][2] * s[2] + x[8][2] * s[3],   x[6][2] * s[0] + x[5][2] * s[1] + x[2][2] * (s[2] - 0.25) + x[9][2] * s[3],   x[7][2] * s[0] + x[8][2] * s[1] + x[9][2] * s[2] + x[3][2] * (s[3] - 0.25)) * 4;
	
	Matrix4d Jinv = J.inverse();

	//PT
	float P[3][4] = { {Jinv(0,1), Jinv(1,1), Jinv(2,1), Jinv(3,1)},
					  {Jinv(0,2), Jinv(1,2), Jinv(2,2), Jinv(3,2)},
					  {Jinv(0,3), Jinv(1,3), Jinv(2,3), Jinv(3,3)} };


	float dN[4][10] = { {4 * s[0] - 1, 0, 0, 0, 4 * s[1], 0, 4 * s[2], 4 * s[3], 0, 0},

						{0, 4 * s[1] - 1, 0, 0, 4 * s[0], 4 * s[2], 0, 0, 4 * s[3], 0},
						
						{0, 0, 4 * s[2] - 1, 0, 0, 4 * s[1], 4 * s[0], 0, 0, 4 * s[3]},
						
						{0, 0, 0, 4 * s[3] - 1, 0, 0, 0, 4 * s[0], 4 * s[1], 4 * s[2]} };

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

		*B = GenMatrix<float,6,30>(b);
		*Jdet = J.determinant();
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

	GenMatrix<float,6,6> E(C);


	K.zeroOut();

	float a = (5.0 + 3.0 * sqrt(5.0))/20.;
	float b = (5.0 - sqrt(5.0))/20.;

	float S[4][4] = { {a,b,b,b},
					  {b,a,b,b},
					  {b,b,a,b},
					  {b,b,b,a} };

	float weight = 0.25;

	for(int i=0;i<4;i++)
	{
		GenMatrix<float, 6, 30> B;
		float det;

		computeB(S[i],&B, &det);
		K = K + (B.transpose() * E * B * (det * weight));
	}

	K = K * (1.0/6.0);
}

QuadTetElement::~QuadTetElement(void)
{
}
