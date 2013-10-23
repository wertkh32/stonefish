#include "TetElement.h"

TetElement::TetElement(Node* n1, Node* n2, Node* n3, Node* n4, float _E, float _v, float _density)
{
	nodes[0] = n1;
	nodes[1] = n2;
	nodes[2] = n3;
	nodes[3] = n4;
	E =_E;
	v = _v;
	density = _density;
	dt = 1.0/FPS;

	//miu = E / (2.0*(1+v));
	//lambda = (E*v)/((1.0+v)*(1.0-2.0*v));

	preCompute();

	mass = density * undeformVolume;
	nodalMass = mass / 4.0;

	nodemass[0] = nodemass[1] = nodemass[2] = nodemass[3] = nodalMass;
}

void TetElement::preCompute()
{
	undeformShapeMat =
		Matrix3d(nodes[0]->pos.x - nodes[3]->pos.x,nodes[1]->pos.x - nodes[3]->pos.x,nodes[2]->pos.x - nodes[3]->pos.x,
		nodes[0]->pos.y - nodes[3]->pos.y,nodes[1]->pos.y - nodes[3]->pos.y,nodes[2]->pos.y - nodes[3]->pos.y,
		nodes[0]->pos.z - nodes[3]->pos.z,nodes[1]->pos.z - nodes[3]->pos.z,nodes[2]->pos.z - nodes[3]->pos.z);
	undeformShapeMatInv = undeformShapeMat.inverse();
	//undeformShapeMatInv.print();
	//undeformShapeMatInvT = undeformShapeMatInv.transpose();
	undeformVolume = (1.0/6.0) * fabs(undeformShapeMat.determinant());

	preComputeUndeformedStiffnessMat();
	preComputeMassMat();
}

void TetElement::preComputeUndeformedStiffnessMat()
{

	//inv is the inverse of the matrix relating volume coords to cartesian coords
	Matrix4d inv =	Matrix4d
		( 1.0, 1.0, 1.0, 1.0,
		nodes[0]->pos.x, nodes[1]->pos.x, nodes[2]->pos.x, nodes[3]->pos.x,
		nodes[0]->pos.y, nodes[1]->pos.y, nodes[2]->pos.y, nodes[3]->pos.y,
		nodes[0]->pos.z, nodes[1]->pos.z, nodes[2]->pos.z, nodes[3]->pos.z).inverse();

	//inv.print();
	//system("pause");

	//for(int i=0;i<4;i++,putchar('\n'))
	//	for(int j=0;j<4;j++)
	//		printf("%f ",inv(i,j));
	//system("pause");

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


	//for(int i=0;i<6;i++,putchar('\n'))
	//	for(int j=0;j<12;j++)
	//		printf("%f ",strainMatrix[i][j]);
	//system("pause");

	GenMatrix<float,6,12> strainMat = GenMatrix<float,6,12>(strainMatrix);
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

	GenMatrix<float,6,6> matConstantsMat = GenMatrix<float,6,6>(C);

	undeformStiffnessMat = strainMat.transpose() * matConstantsMat * strainMat;
	undeformStiffnessMat.scalarMul(undeformVolume);

	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			B[i][j] = inv(i,j+1);

	con1 = c1 * undeformVolume * (1.0/(2.0 * undeformVolume)) * (1.0/(2.0 * undeformVolume));
	con2 = c2 * undeformVolume * (1.0/(2.0 * undeformVolume)) * (1.0/(2.0 * undeformVolume));
	/*
	printf("start");
	for(int i=0;i<12;i++,printf("\n"))
	for(int j=0;j<12;j++)
	{
	printf("%f ",undeformStiffnessMat(i,j));
	}
	*/

	/*
	for(int i=0;i<4;i++)
	for(int j=0;j<4;j++)
	sparseStiff->setBlockFilled(i,j,true);

	for(int i=0;i<12;i++)
	for(int j=0;j<12;j++)
	{
	sparseStiff->setValue(i,j,undeformStiffnessMat(i,j));
	}
	*/

}

void TetElement::preComputeMassMat()
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
	massMat.scalarMul( (density * undeformVolume) / 20.);

	/*
	for(int i=0;i<12;i++)
	{
	for(int j=0;j<12;j++)
	printf("%f ", massMat(i,j));
	printf("\n");
	}
	*/
}

Matrix3d TetElement::computeDeformationMat()
{
	Matrix3d deformShapeMat
		(nodes[0]->pos_t.x - nodes[3]->pos_t.x,nodes[1]->pos_t.x - nodes[3]->pos_t.x,nodes[2]->pos_t.x - nodes[3]->pos_t.x,
		nodes[0]->pos_t.y - nodes[3]->pos_t.y,nodes[1]->pos_t.y - nodes[3]->pos_t.y,nodes[2]->pos_t.y - nodes[3]->pos_t.y,
		nodes[0]->pos_t.z - nodes[3]->pos_t.z,nodes[1]->pos_t.z - nodes[3]->pos_t.z,nodes[2]->pos_t.z - nodes[3]->pos_t.z);
	return deformShapeMat * undeformShapeMatInv;
}



Matrix3d TetElement::computeDeformShapeMat()
{
	Matrix3d deformShapeMat
		(nodes[0]->pos_t.x - nodes[3]->pos_t.x,nodes[1]->pos_t.x - nodes[3]->pos_t.x,nodes[2]->pos_t.x - nodes[3]->pos_t.x,
		nodes[0]->pos_t.y - nodes[3]->pos_t.y,nodes[1]->pos_t.y - nodes[3]->pos_t.y,nodes[2]->pos_t.y - nodes[3]->pos_t.y,
		nodes[0]->pos_t.z - nodes[3]->pos_t.z,nodes[1]->pos_t.z - nodes[3]->pos_t.z,nodes[2]->pos_t.z - nodes[3]->pos_t.z);
	return deformShapeMat;
}

void TetElement::computeRotation()
{
	Matrix3d F,RR,S;
	F = computeDeformationMat();
	PolarDecompose::compute(F,RR,S);

	R = RR;
}

void
TetElement::computeRKRTandRK()
{
	static const float alpha = 0.1, beta = 0.1;
	static const float coeffK = dt * beta + dt * dt, coeffM = 1 + dt * alpha;

	Matrix3d F,R,S;
	F = computeDeformationMat();
	PolarDecompose::compute(F,R,S);

	//for(int i=0;i<4;i++)
	// for(int j=0;j<3;j++)
	// for(int k=0;k<3;k++)
	// {
	// Rot(i * 3 + j, i * 3 + k) = R(j,k);
	// }


	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			for(int a=0;a<3;a++)
				for(int b=0;b<3;b++)
				{
					RK(a + i * 3, b + j * 3) = 0;
					for(int c=0;c<3;c++)
						RK(a + i * 3, b + j * 3) += R(a,c) * undeformStiffnessMat(c + i * 3, b + j * 3);
				}
		}
	}

	//Matrix3d RT = R.transpose();
	//upper triangle
	for(int i=0;i<4;i++)
	{
		for(int j=i;j<4;j++)
		{
			for(int a=0;a<3;a++)
				for(int b=0;b<3;b++)
				{
					RKRT(a + i * 3, b + j * 3) = 0;
					for(int c=0;c<3;c++)
						RKRT(a + i * 3, b + j * 3) += RK(a + i * 3, c + j * 3) * R(b , c); //R(c, b) but its RT so, R(b , c)
				}
		}
	}
	//lower triangle
	for(int i=1;i<4;i++)
		for(int j=0;j<i;j++)
		{
			for(int a=0;a<3;a++)
				for(int b=0;b<3;b++)
					RKRT(a + i * 3, b + j * 3) = RKRT(b + j * 3, a + i * 3);
		}

		for(int i=0;i<12;i++)
			for(int j=0;j<12;j++)
			{
				A(i,j) = coeffK * RKRT(i,j);
				if(i==j)
					A(i,j) += coeffM * nodalMass;
			}

}

void TetElement::getRKRTandRK(GenMatrix<float,12,12>*& ptrRK, GenMatrix<float,12,12>*& ptrRKRT)
{
	computeRKRTandRK();
	//RK = Rot * undeformStiffnessMat;
	//RKRT = RK * Rot.transpose();
	ptrRK = &RK;
	ptrRKRT = &RKRT;
}

void TetElement::getRKRTandRK(SparseMatrix& RK, SparseMatrix& RKRT)
{
	/*
	Matrix3d F,R,S;
	F = computeDeformationMat();
	PolarDecompose::compute(F,R,S);

	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			RK.setBlock(i,j,R * sparseStiff->getBlock(i,j),true);
		}
	}


	Matrix3d RT = R.transpose();

	//upper triangle
	for(int i=0;i<4;i++)
	{
		for(int j=i;j<4;j++)
		{
			RKRT.setBlock(i,j,RK.getBlock(i,j) * RT,true);
		}
	}
	//lower triangle
	for(int i=1;i<4;i++)
	{
		for(int j=0;j<i;j++)
		{
			RKRT.setBlock(i,j,RKRT.getBlock(j,i).transpose(),true);
		}
	}
	*/
}
/*
void Element::computeMatFreeVars()
{
Matrix3d deformShapeMat
(nodes[0]->pos_t.x - nodes[3]->pos_t.x,nodes[1]->pos_t.x - nodes[3]->pos_t.x,nodes[2]->pos_t.x - nodes[3]->pos_t.x,
nodes[0]->pos_t.y - nodes[3]->pos_t.y,nodes[1]->pos_t.y - nodes[3]->pos_t.y,nodes[2]->pos_t.y - nodes[3]->pos_t.y,
nodes[0]->pos_t.z - nodes[3]->pos_t.z,nodes[1]->pos_t.z - nodes[3]->pos_t.z,nodes[2]->pos_t.z - nodes[3]->pos_t.z);
Ft = deformShapeMat * undeformShapeMatInv;
PolarDecompose::computeFull(Ft,Rt,St);
float traceS = St.trace();
Matrix3d trSI = Matrix3d(traceS,0,0,
0,traceS,0,
0,0,traceS);
trSI_Sinv = (trSI - St).inverse();
//trSI_Sinv.print();
S_Itrace = (St - Matrix3d(1,0,0,
0,1,0,
0,0,1)).trace();
}

Matrix3d
Element::computePiolaStressTensor()
{
//Matrix3d I = Matrix3d(1,0,0,
// 0,1,0,
// 0,0,1);
//float dF_Itrace = (Ft-I).trace();
return //(Ft + Ft.transpose() - I * 2) * miu + I * dF_Itrace;
(Ft - Rt) * (2 * miu) + Rt * (S_Itrace * lambda);
}

Matrix3d
Element::computeRKtensor(Matrix3d& dF)
{
Matrix3d I = Matrix3d(1,0,0,
0,1,0,
0,0,1);
//Assume compute Ft Rt St has been computed
Matrix3d RtT = Rt.transpose();
Matrix3d W = RtT * dF;
vector3<float> w = vector3<float>(W(1,2) - W(2,1), W(2,0) - W(0,2),W(0,1)-W(1,0));
Matrix3d dR = Rt * Matrix3d::skew(trSI_Sinv * w);

float F_Itrace = (Ft-I).trace();
float dFtrace = dF.trace();

return //dR * ((Ft + Ft.transpose() - I * 2) * miu + (I *(F_Itrace * lambda))) + Rt *
((dF + dF.transpose()) * miu + (I * (dFtrace * lambda)));

}

Matrix3d
Element::computeDiffPiolaStressTensor(Matrix3d& dF)
{
Matrix3d I = Matrix3d(1,0,0,
0,1,0,
0,0,1);
//Assume compute Ft Rt St has been computed
Matrix3d RtT = Rt.transpose();
Matrix3d W = RtT * dF;
vector3<float> w = vector3<float>(W(1,2) - W(2,1), W(2,0) - W(0,2),W(0,1)-W(1,0));

Matrix3d dR = Rt * Matrix3d::skew(trSI_Sinv * w);
// 2EdF + v * tr(W) * R + {v * tr(S - I) - 2 * E} * dR
float dFtrace = dF.trace();
//printf("[miu: %f, lambda: %f, dftrace:%f]\n",miu,lambda,dFtrace);
//(I * (dFtrace * lambda)).print();

return ((dF + dF.transpose()) * miu + (I * (dFtrace * lambda)));
//dF * (2 * miu) + Rt * (W.trace() * lambda) + dR * (v * S_Itrace - 2 * miu);
}
*/

void
TetElement::renderElement()
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

TetElement::~TetElement(void)
{
}