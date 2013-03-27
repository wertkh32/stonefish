#pragma once
#include "vector3.h"
class Matrix4d
{
public:
	double mat[4][4];
	Matrix4d(void);
	Matrix4d(double[4][4]);
	Matrix4d(const Matrix4d& m);
	Matrix4d(double a00,double a01,double a02,double a03,
			 double a10,double a11,double a12,double a13,
			 double a20,double a21,double a22,double a23,
			 double a30,double a31,double a32,double a33);

	inline Matrix4d transpose();
	inline double determinant();
	inline Matrix4d operator*(Matrix4d&);
	inline Matrix4d operator*(double);
	inline Matrix4d operator+(Matrix4d&);
	inline Matrix4d operator-(Matrix4d&);
	double& operator()(int i,int j){return mat[i][j];}

	Matrix4d inverse();

	~Matrix4d(void);
};

inline
double Matrix4d::determinant()
{
	return (mat[0][0] * (-mat[2][3] * mat[3][2] * mat[1][1] + mat[2][2] * mat[3][3] * mat[1][1] + mat[2][3] * mat[3][1] * mat[1][2] - mat[2][2] * mat[3][1] * mat[1][3] - mat[3][3] * mat[1][2] * mat[2][1] + mat[3][2] * mat[1][3] * mat[2][1]) + 
			mat[0][1] * (mat[2][3] * mat[3][2] * mat[1][0] - mat[2][2] * mat[3][3] * mat[1][0] - mat[2][3] * mat[3][0] * mat[1][2] + mat[2][2] * mat[3][0] * mat[1][3] + mat[3][3] * mat[1][2] * mat[2][0] - mat[3][2] * mat[1][3] * mat[2][0]) + 
			mat[0][2] * (-mat[2][3] * mat[3][1] * mat[1][0] + mat[2][3] * mat[3][0] * mat[1][1] - mat[3][3] * mat[1][1] * mat[2][0] + mat[3][1] * mat[1][3] * mat[2][0] + mat[3][3] * mat[1][0] * mat[2][1] - mat[3][0] * mat[1][3] * mat[2][1]) + 
			mat[0][3] * (mat[2][2] * mat[3][1] * mat[1][0] - mat[2][2] * mat[3][0] * mat[1][1] + mat[3][2] * mat[1][1] * mat[2][0] - mat[3][1] * mat[1][2] * mat[2][0] - mat[3][2] * mat[1][0] * mat[2][1] + mat[3][0] * mat[1][2] * mat[2][1]));
}

inline
Matrix4d Matrix4d::transpose(){
return Matrix4d(mat[0][0],mat[1][0],mat[2][0],mat[3][0],
				mat[0][1],mat[1][1],mat[2][1],mat[3][1],
				mat[0][2],mat[1][2],mat[2][2],mat[3][2],
				mat[0][3],mat[1][3],mat[2][3],mat[3][3]);
}

inline
Matrix4d Matrix4d::operator*(Matrix4d& mm){
	double m[4][4]={0};
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			for(int k=0;k<4;k++)
				m[i][j]+= mat[i][k]*mm.mat[k][j];
	return Matrix4d(m);
}

inline
Matrix4d Matrix4d::operator*(double k){
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
				mat[i][j] *= k;
	return Matrix4d(mat);
}

inline
Matrix4d Matrix4d::operator+(Matrix4d& mm){
	double m[4][4];
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			m[i][j] = mat[i][j]+mm.mat[i][j];
	return Matrix4d(m);
}

inline
Matrix4d Matrix4d::operator-(Matrix4d& mm){
	double m[4][4];
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			m[i][j] = mat[i][j]-mm.mat[i][j];
	return Matrix4d(m);
}

