#include "Element.h"


Element::Element(Node* n1, Node* n2, Node* n3, Node* n4, float _E, float _v, float _density)
{
	nodes[0] = n1;
	nodes[1] = n2;
	nodes[2] = n3;
	nodes[3] = n4;
	E =_E;
	v = _v;
	density = _density;
	
	preCompute();
}

void Element::preCompute()
{
	undeformShapeMat = 
			Matrix3d(nodes[0]->pos.x - nodes[3]->pos.x,nodes[1]->pos.x - nodes[3]->pos.x,nodes[2]->pos.x - nodes[3]->pos.x,
					nodes[0]->pos.y - nodes[3]->pos.y,nodes[1]->pos.y - nodes[3]->pos.y,nodes[2]->pos.y - nodes[3]->pos.y,
					nodes[0]->pos.z - nodes[3]->pos.z,nodes[1]->pos.z - nodes[3]->pos.z,nodes[2]->pos.z - nodes[3]->pos.z);
	undeformShapeMatInv = undeformShapeMat.inverse();
	undeformVolume = (1.0/6.0) * fabs(undeformShapeMat.determinant());

	preComputeUndeformedStiffnessMat();
	preComputeMassMat();
}

void Element::preComputeUndeformedStiffnessMat()
{
	
	//inv is the inverse of the matrix relating volume coords to cartesian coords
	Matrix4d inv =	  Matrix4d
					  ( 1.0, 1.0, 1.0, 1.0,
						nodes[0]->pos.x, nodes[1]->pos.x, nodes[2]->pos.x, nodes[3]->pos.x,
						nodes[0]->pos.y, nodes[1]->pos.y, nodes[2]->pos.y, nodes[3]->pos.y,
						nodes[0]->pos.z, nodes[1]->pos.z, nodes[2]->pos.z, nodes[3]->pos.z).inverse();

	//strain matrix B = LN = dN/dx
	//checked correct.
	float strainMatrix[6][12] = 
	{
		{ inv(0,1), 0, 0, inv(1,1), 0, 0, inv(2,1), 0, 0, inv(3,1), 0, 0 },
		{ 0, inv(0,2), 0, 0, inv(1,2), 0, 0, inv(2,2), 0, 0, inv(3,2), 0 },
		{ 0, 0, inv(0,3), 0, 0, inv(1,3), 0, 0, inv(2,3), 0, 0, inv(3,3) },
		{ inv(0,2), inv(0,1), 0, inv(1,2), inv(1,1), 0, inv(2,2), inv(2,1), 0, inv(3,2), inv(3,1), 0 },
		{ 0, inv(0,3), inv(0,2), 0, inv(1,3), inv(1,2), 0, inv(2,3), inv(2,2), 0, inv(3,3), inv(3,2) },
		{ inv(0,3), 0, inv(0,1), inv(1,3), 0, inv(1,1), inv(2,3), 0, inv(2,1), inv(3,3), 0, inv(3,1) }
	};

	strainMat = GenMatrix<float,6,12>(strainMatrix);
	strainMat.scalarMul(1/(2 * undeformVolume));

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

	matConstantsMat = GenMatrix<float,6,6>(C);

	undeformStiffnessMat = strainMat.transpose() * matConstantsMat * strainMat;
	undeformStiffnessMat.scalarMul(undeformVolume);
	/*
	for(int i=0;i<12;i++)
		for(int j=0;j<12;j++)
		{
			printf("%lf ",undeformStiffnessMat(i,j));
		}
	*/
}

void Element::preComputeMassMat()
{
	float mass[12][12] = 
	{
		{ 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0 },
		{ 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0 },
		{ 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1 },
		{ 1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0 },
		{ 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0 },
		{ 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 1 },
		{ 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0 },
		{ 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0 },
		{ 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1 },
		{ 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0 },
		{ 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0 },
		{ 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2 }
	};

	massMat = GenMatrix<float,12,12>(mass);
	//massMat.scalarMul( (density * undeformVolume) / 20.);
	
	/*
	for(int i=0;i<12;i++)
	{
		for(int j=0;j<12;j++)
			printf("%f ", massMat(i,j));
		printf("\n");
	}
	*/
}

Matrix3d Element::computeDeformationMat()
{
	Matrix3d deformShapeMat 
				   (nodes[0]->pos_t.x - nodes[3]->pos_t.x,nodes[1]->pos_t.x - nodes[3]->pos_t.x,nodes[2]->pos_t.x - nodes[3]->pos_t.x,
					nodes[0]->pos_t.y - nodes[3]->pos_t.y,nodes[1]->pos_t.y - nodes[3]->pos_t.y,nodes[2]->pos_t.y - nodes[3]->pos_t.y,
					nodes[0]->pos_t.z - nodes[3]->pos_t.z,nodes[1]->pos_t.z - nodes[3]->pos_t.z,nodes[2]->pos_t.z - nodes[3]->pos_t.z);
	return deformShapeMat * undeformShapeMatInv;
}

Matrix3d Element::computeDeformShapeMat()
{
	Matrix3d deformShapeMat 
				   (nodes[0]->pos_t.x - nodes[3]->pos_t.x,nodes[1]->pos_t.x - nodes[3]->pos_t.x,nodes[2]->pos_t.x - nodes[3]->pos_t.x,
					nodes[0]->pos_t.y - nodes[3]->pos_t.y,nodes[1]->pos_t.y - nodes[3]->pos_t.y,nodes[2]->pos_t.y - nodes[3]->pos_t.y,
					nodes[0]->pos_t.z - nodes[3]->pos_t.z,nodes[1]->pos_t.z - nodes[3]->pos_t.z,nodes[2]->pos_t.z - nodes[3]->pos_t.z);
	return deformShapeMat;
}

void Element::getRKRTandRK(GenMatrix<float,12,12>& RK, GenMatrix<float,12,12>& RKRT)
{
	Matrix3d F,R,S;
	GenMatrix<float,12,12> Rot;
	F = computeDeformationMat();
	PolarDecompose::compute(F,R,S);
	
	for(int i=0;i<4;i++)
		for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
			{
				Rot(i * 3 + j, i * 3 + k) = R(j,k);
			}

	RK = Rot * undeformStiffnessMat;
	RKRT = RK * Rot.transpose();
}

void
Element::renderElement()
{
	glPushMatrix();
	glBegin(GL_LINES);
	glVertex3fv(nodes[3]->pos_t.coords);
	glVertex3fv(nodes[0]->pos_t.coords);

	glVertex3fv(nodes[3]->pos_t.coords);
	glVertex3fv(nodes[1]->pos_t.coords);

	glVertex3fv(nodes[3]->pos_t.coords);
	glVertex3fv(nodes[2]->pos_t.coords);


	glVertex3fv(nodes[0]->pos_t.coords);
	glVertex3fv(nodes[1]->pos_t.coords);

	glVertex3fv(nodes[1]->pos_t.coords);
	glVertex3fv(nodes[2]->pos_t.coords);

	glVertex3fv(nodes[2]->pos_t.coords);
	glVertex3fv(nodes[0]->pos_t.coords);

	glEnd();

	
	glPopMatrix();
}

Element::~Element(void)
{
}
