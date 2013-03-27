#include "Matrix4d.h"


Matrix4d::Matrix4d(void)
{
	mat[0][0]=1;mat[0][1]=0;mat[0][2]=0;mat[0][3]=0;
	mat[1][0]=0;mat[1][1]=1;mat[1][2]=0;mat[1][3]=0;
	mat[2][0]=0;mat[2][1]=0;mat[2][2]=1;mat[2][3]=0;
	mat[3][0]=0;mat[3][1]=0;mat[3][2]=0;mat[3][3]=1;
}

Matrix4d::Matrix4d(const Matrix4d& mm)
{
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			mat[i][j]=mm.mat[i][j];
}

Matrix4d::Matrix4d(double m[4][4]){
for(int i=0;i<4;i++)
	for(int j=0;j<4;j++)
		mat[i][j]=m[i][j];
}

Matrix4d::Matrix4d(double a00,double a01,double a02,double a03,
					double a10,double a11,double a12,double a13,
					double a20,double a21,double a22,double a23,
					double a30,double a31,double a32,double a33)
{
	mat[0][0]=a00;mat[0][1]=a01;mat[0][2]=a02;mat[0][3]=a03;
	mat[1][0]=a10;mat[1][1]=a11;mat[1][2]=a12;mat[1][3]=a13;
	mat[2][0]=a20;mat[2][1]=a21;mat[2][2]=a22;mat[2][3]=a23;
	mat[3][0]=a30;mat[3][1]=a31;mat[3][2]=a32;mat[3][3]=a33;
	
}

Matrix4d Matrix4d::inverse()
{
	//from Mathematica
	double result[4][4];

	result[0][0] = -mat[2][3] * mat[3][2] * mat[1][1] + mat[2][2] * mat[3][3] * mat[1][1] + mat[2][3] * mat[3][1] * mat[1][2] - mat[2][2] * mat[3][1] * mat[1][3] - mat[3][3] * mat[1][2] * mat[2][1] + mat[3][2] * mat[1][3] * mat[2][1];
	result[0][1] = mat[0][1] * mat[2][3] * mat[3][2] - mat[0][1] * mat[2][2] * mat[3][3] - mat[2][3] * mat[3][1] * mat[0][2] + mat[2][2] * mat[3][1] * mat[0][3] + mat[3][3] * mat[0][2] * mat[2][1] - mat[3][2] * mat[0][3] * mat[2][1];
	result[0][2] = -mat[3][3] * mat[0][2] * mat[1][1] + mat[3][2] * mat[0][3] * mat[1][1] + mat[0][1] * mat[3][3] * mat[1][2] - mat[3][1] * mat[0][3] * mat[1][2] - mat[0][1] * mat[3][2] * mat[1][3] + mat[3][1] * mat[0][2] * mat[1][3];
	result[0][3] = mat[2][3] * mat[0][2] * mat[1][1] - mat[2][2] * mat[0][3] * mat[1][1] - mat[0][1] * mat[2][3] * mat[1][2] + mat[0][1] * mat[2][2] * mat[1][3] + mat[0][3] * mat[1][2] * mat[2][1] - mat[0][2] * mat[1][3] * mat[2][1];
	result[1][0] = mat[2][3] * mat[3][2] * mat[1][0] - mat[2][2] * mat[3][3] * mat[1][0] - mat[2][3] * mat[3][0] * mat[1][2] + mat[2][2] * mat[3][0] * mat[1][3] + mat[3][3] * mat[1][2] * mat[2][0] - mat[3][2] * mat[1][3] * mat[2][0];
	result[1][1] = -mat[0][0] * mat[2][3] * mat[3][2] + mat[0][0] * mat[2][2] * mat[3][3] + mat[2][3] * mat[3][0] * mat[0][2] - mat[2][2] * mat[3][0] * mat[0][3] - mat[3][3] * mat[0][2] * mat[2][0] + mat[3][2] * mat[0][3] * mat[2][0];
	result[1][2] = mat[3][3] * mat[0][2] * mat[1][0] - mat[3][2] * mat[0][3] * mat[1][0] - mat[0][0] * mat[3][3] * mat[1][2] + mat[3][0] * mat[0][3] * mat[1][2] + mat[0][0] * mat[3][2] * mat[1][3] - mat[3][0] * mat[0][2] * mat[1][3];
	result[1][3] = -mat[2][3] * mat[0][2] * mat[1][0] + mat[2][2] * mat[0][3] * mat[1][0] + mat[0][0] * mat[2][3] * mat[1][2] - mat[0][0] * mat[2][2] * mat[1][3] - mat[0][3] * mat[1][2] * mat[2][0] + mat[0][2] * mat[1][3] * mat[2][0];
	result[2][0] = -mat[2][3] * mat[3][1] * mat[1][0] + mat[2][3] * mat[3][0] * mat[1][1] - mat[3][3] * mat[1][1] * mat[2][0] + mat[3][1] * mat[1][3] * mat[2][0] + mat[3][3] * mat[1][0] * mat[2][1] - mat[3][0] * mat[1][3] * mat[2][1];
	result[2][1] = -mat[0][1] * mat[2][3] * mat[3][0] + mat[0][0] * mat[2][3] * mat[3][1] + mat[0][1] * mat[3][3] * mat[2][0] - mat[3][1] * mat[0][3] * mat[2][0] - mat[0][0] * mat[3][3] * mat[2][1] + mat[3][0] * mat[0][3] * mat[2][1];
	result[2][2] = -mat[0][1] * mat[3][3] * mat[1][0] + mat[3][1] * mat[0][3] * mat[1][0] + mat[0][0] * mat[3][3] * mat[1][1] - mat[3][0] * mat[0][3] * mat[1][1] + mat[0][1] * mat[3][0] * mat[1][3] - mat[0][0] * mat[3][1] * mat[1][3];
	result[2][3] = mat[0][1] * mat[2][3] * mat[1][0] - mat[0][0] * mat[2][3] * mat[1][1] + mat[0][3] * mat[1][1] * mat[2][0] - mat[0][1] * mat[1][3] * mat[2][0] - mat[0][3] * mat[1][0] * mat[2][1] + mat[0][0] * mat[1][3] * mat[2][1]; 
	result[3][0] = mat[2][2] * mat[3][1] * mat[1][0] - mat[2][2] * mat[3][0] * mat[1][1] + mat[3][2] * mat[1][1] * mat[2][0] - mat[3][1] * mat[1][2] * mat[2][0] - mat[3][2] * mat[1][0] * mat[2][1] + mat[3][0] * mat[1][2] * mat[2][1];
	result[3][1] = mat[0][1] * mat[2][2] * mat[3][0] - mat[0][0] * mat[2][2] * mat[3][1] - mat[0][1] * mat[3][2] * mat[2][0] + mat[3][1] * mat[0][2] * mat[2][0] + mat[0][0] * mat[3][2] * mat[2][1] - mat[3][0] * mat[0][2] * mat[2][1]; 
	result[3][2] = mat[0][1] * mat[3][2] * mat[1][0] - mat[3][1] * mat[0][2] * mat[1][0] - mat[0][0] * mat[3][2] * mat[1][1] + mat[3][0] * mat[0][2] * mat[1][1] - mat[0][1] * mat[3][0] * mat[1][2] + mat[0][0] * mat[3][1] * mat[1][2];
	result[3][3] = -mat[0][1] * mat[2][2] * mat[1][0] + mat[0][0] * mat[2][2] * mat[1][1] - mat[0][2] * mat[1][1] * mat[2][0] + mat[0][1] * mat[1][2] * mat[2][0] + mat[0][2] * mat[1][0] *mat[2][1] - mat[0][0] * mat[1][2] * mat[2][1];

	double d = 1.0 / (mat[0][0] * result[0][0] + mat[0][1] * result[1][0] + mat[0][2] * result[2][0] + mat[0][3] * result[3][0]);

	for(int i=0; i<4; i++)
	  for(int j=0; j<4; j++)
			result[i][j] *= d;
	return Matrix4d(result);
}

Matrix4d::~Matrix4d(void)
{
}
