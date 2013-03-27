#pragma once
#include "vector3.h"
class Matrix3d
{
public:
	double mat[3][3];
	Matrix3d(void);
	Matrix3d(double[3][3]);
	Matrix3d(const Matrix3d& m);
	Matrix3d(double a00,double a01,double a02,
			 double a10,double a11,double a12,
			 double a20,double a21,double a22);
	inline Matrix3d transpose();
	inline double determinant();
	inline Matrix3d operator*(Matrix3d&);
	inline Matrix3d operator*(double);
	inline vector3<double> operator*(vector3<double>&);
	inline Matrix3d operator+(Matrix3d&);
	inline Matrix3d operator-(Matrix3d&);
	double& operator()(int i,int j){return mat[i][j];}

	Matrix3d inverse();
	static Matrix3d skew(vector3<double>& v);

	~Matrix3d(void);
};

inline
double Matrix3d::determinant()
{
	return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) 
		 - mat[0][1] * (mat[2][2] * mat[1][0] - mat[1][2] * mat[2][0]) 
		 + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
}

inline
Matrix3d Matrix3d::transpose(){
return Matrix3d(mat[0][0],mat[1][0],mat[2][0],
				mat[0][1],mat[1][1],mat[2][1],
				mat[0][2],mat[1][2],mat[2][2]);
}

inline
Matrix3d Matrix3d::operator*(Matrix3d& mm){
	double m[3][3]={0};
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
				m[i][j]+= mat[i][k]*mm.mat[k][j];
	return Matrix3d(m);
}

inline
vector3<double> Matrix3d::operator*(vector3<double>& vv){
	double v[3]={0};
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			v[i]+=mat[i][j] * vv.coords[j];
	return vector3<double>(v);
}

inline
Matrix3d Matrix3d::operator*(double k){
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
				mat[i][j] *= k;
	return Matrix3d(mat);
}

inline
Matrix3d Matrix3d::operator+(Matrix3d& mm){
	double m[3][3];
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			m[i][j] = mat[i][j]+mm.mat[i][j];
	return Matrix3d(m);
}

inline
Matrix3d Matrix3d::operator-(Matrix3d& mm){
	double m[3][3];
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			m[i][j] = mat[i][j]-mm.mat[i][j];
	return Matrix3d(m);
}

