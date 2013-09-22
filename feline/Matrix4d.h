#pragma once
#include "vector3.h"
class Matrix4d
{
public:
	float mat[4][4];
	Matrix4d(void);
	Matrix4d(float[4][4]);
	Matrix4d(const Matrix4d& m);
	Matrix4d(float a00,float a01,float a02,float a03,
			 float a10,float a11,float a12,float a13,
			 float a20,float a21,float a22,float a23,
			 float a30,float a31,float a32,float a33);

	inline Matrix4d transpose();
	inline float determinant();
	inline Matrix4d operator*(Matrix4d&);
	inline Matrix4d operator*(float);
	inline Matrix4d operator+(Matrix4d&);
	inline Matrix4d operator-(Matrix4d&);
	float& operator()(int i,int j){return mat[i][j];}

	Matrix4d inverse();

	~Matrix4d(void);
};

inline
float Matrix4d::determinant()
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
	float m[4][4]={0};
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			for(int k=0;k<4;k++)
				m[i][j]+= mat[i][k]*mm.mat[k][j];
	return Matrix4d(m);
}

inline
Matrix4d Matrix4d::operator*(float k){
	float m[4][4] = {0};
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
				m[i][j] = mat[i][j] * k;
	return Matrix4d(m);
}

inline
Matrix4d Matrix4d::operator+(Matrix4d& mm){
	float m[4][4];
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			m[i][j] = mat[i][j]+mm.mat[i][j];
	return Matrix4d(m);
}

inline
Matrix4d Matrix4d::operator-(Matrix4d& mm){
	float m[4][4];
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			m[i][j] = mat[i][j]-mm.mat[i][j];
	return Matrix4d(m);
}

