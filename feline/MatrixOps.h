#pragma once
/////
//matrix lib for T**s 
/////
template <typename T>
void MatrixAdd(T** mat1, T** mat2, int n, int m)
{
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			mat1[i][j] += mat[i][j];
}

template <typename T>
void MatrixSub(T** mat1, T** mat2, int n, int m)
{
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			mat1[i][j] -= mat[i][j];
}


template <typename T>
void MatrixScalarMul(T** mat1, int k, int n, int m)
{
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			mat1[i][j] *= k;
}

template <typename T>
void MatrixVectorMul(T** mat, T* vec, T* outvec, int n, int m)
{
	for(int i=0;i<n;i++)
	{
		outvec[i] = 0;
		for(int j=0;j<m;j++)
		{
			outvec[i] += mat[i][j] * vec[j];
		}
	}
}

