#pragma once
#include "vector3.h"
class Matrix3d
{
public:
	float mat[3][3];
	Matrix3d(void);
	Matrix3d(float[3][3]);
	Matrix3d(const Matrix3d& m);
	Matrix3d(float a00,float a01,float a02,
			 float a10,float a11,float a12,
			 float a20,float a21,float a22);
	inline Matrix3d transpose();
	inline float determinant();
	inline Matrix3d operator*(Matrix3d&);
	inline Matrix3d operator*(float);
	inline vector3<float> operator*(vector3<float>&);
	inline Matrix3d operator+(Matrix3d&);
	inline Matrix3d operator-(Matrix3d&);
	float& operator()(int i,int j){return mat[i][j];}

	Matrix3d inverse();
	static Matrix3d skew(vector3<float>& v);

	~Matrix3d(void);
};

inline
float Matrix3d::determinant()
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
	float m[3][3]={0};
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
				m[i][j]+= mat[i][k]*mm.mat[k][j];
	return Matrix3d(m);
}

inline
vector3<float> Matrix3d::operator*(vector3<float>& vv){
	float v[3]={0};
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			v[i]+=mat[i][j] * vv.coords[j];
	return vector3<float>(v);
}

inline
Matrix3d Matrix3d::operator*(float k){
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
				mat[i][j] *= k;
	return Matrix3d(mat);
}

inline
Matrix3d Matrix3d::operator+(Matrix3d& mm){
	float m[3][3];
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			m[i][j] = mat[i][j]+mm.mat[i][j];
	return Matrix3d(m);
}

inline
Matrix3d Matrix3d::operator-(Matrix3d& mm){
	float m[3][3];
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			m[i][j] = mat[i][j]-mm.mat[i][j];
	return Matrix3d(m);
}

