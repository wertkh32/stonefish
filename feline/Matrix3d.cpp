#include "Matrix3d.h"


Matrix3d::Matrix3d(void)
{
	mat[0][0]=1;mat[0][1]=0;mat[0][2]=0;
	mat[1][0]=0;mat[1][1]=1;mat[1][2]=0;
	mat[2][0]=0;mat[2][1]=0;mat[2][2]=1;
}

Matrix3d::Matrix3d(const Matrix3d& mm)
{
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			mat[i][j]=mm.mat[i][j];
}

Matrix3d::Matrix3d(double m[3][3]){
for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
		mat[i][j]=m[i][j];
}

Matrix3d::Matrix3d(double a00,double a01,double a02,
				   double a10,double a11,double a12,
				   double a20,double a21,double a22)
{
	mat[0][0]=a00;mat[0][1]=a01;mat[0][2]=a02;
	mat[1][0]=a10;mat[1][1]=a11;mat[1][2]=a12;
	mat[2][0]=a20;mat[2][1]=a21;mat[2][2]=a22;
	
}

Matrix3d Matrix3d::skew(vector3<double>& v){
	return Matrix3d(0, -v.z, v.y,
					v.z, 0, -v.x,
					-v.y, v.x, 0);
}

Matrix3d Matrix3d::inverse()
{
	double d = determinant();

	return Matrix3d( (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])/d , -(mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1])/d ,  (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1])/d,
					-(mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])/d ,  (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0])/d , -(mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0])/d,
					 (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0])/d , -(mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0])/d ,  (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0])/d );
}

Matrix3d::~Matrix3d(void)
{
}
