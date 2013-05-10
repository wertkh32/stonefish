#pragma once
#include <stdio.h>
#include <stdlib.h>
#include "GenVector.h"

template <class T,int N, int M>
class GenMatrix
{
public:
	T mat[N][M];
	GenMatrix(void)
	{
		for(int i=0;i<N;i++)
			for(int j=0;j<M;j++)
				mat[i][j] = 0;
	}

	GenMatrix(T m[N][M])
	{
		for(int i=0;i<N;i++)
			for(int j=0;j<M;j++)
				mat[i][j] = m[i][j];
	}

	GenMatrix(GenMatrix<T,N,M>& m)
	{
		for(int i=0;i<N;i++)
			for(int j=0;j<M;j++)
				mat[i][j] = m(i,j);
	}

	T& operator()(int x, int y){return mat[x][y];}
	
	template<int K>
	GenMatrix<T,N,K> operator*(GenMatrix<T,M,K>& mat2);

	GenMatrix<T,M,N> transpose();
	GenMatrix<T,M,N> inverse();

	void scalarMul(T k)
	{
		for(int i=0;i<N;i++)
			for(int j=0;j<M;j++)
				mat[i][j] *= k;
	}

};

template<class T,int N, int M>
template<int K>
GenMatrix<T,N,K> GenMatrix<T,N,M>::operator*(GenMatrix<T,M,K>& mat2)
{
	GenMatrix<T,N,K> result;
	for(int i = 0; i<N; i++)
		for(int j = 0; j<K; j++)
			for(int k = 0; k<M; k++)
				result(i,j) += mat[i][k] * mat2(k,j);
	return result;
}

template<class T,int N, int M>
GenMatrix<T,M,N> GenMatrix<T,N,M>::transpose()
{
	GenMatrix<T,M,N> result;
	for(int i=0; i<N; i++)
		for(int j=0;j<M;j++)
			result(j,i) = mat[i][j];

	return result;
}

template<class T,int N, int M>
GenMatrix<T,M,N> GenMatrix<T,N,M>::inverse()
{
	//////////////////////////////////////////
	//simulate find of inverse using a float width matrix in form
	//[A | I] -> [I | A^-1]
	//idea from http://programming-technique.blogspot.sg/2011/09/numerical-methods-inverse-of-nxn-matrix.html
	//////////////////////////////////////////

	if(N!=M) return GenMatrix<T,M,N>();
	else
	{
		T result[N][N * 2];

		//build [ A | I ]
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				result[i][j] = mat[i][j];

		for(int i=0;i<M;i++)
			for(int j=N;j<N * 2;j++)
				if(i==(j-N))
					result[i][j] = 1;
				else
					result[i][j] = 0;

		for(int i = 0; i < N; i++)
		{
			for(int j = 0; j < N; j++)
			{
				if(i!=j)
				{
					T ratio = result[j][i]/result[i][i];
					for(int k = 0; k < 2*n; k++)
						result[j][k] -= ratio * result[i][k];
				}
			}
		}

		for(int i = 0; i < n; i++)
		{
			T a = result[i][i];
			for(int j = 0; j < 2*n; j++)
				result[i][j] /= a;
		}

		return GenMatrix<T,M,N>(result);
	}
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



template<class T>
class GenMatrixDyn
{
public:
	T** mat;
	int n,m;
	GenMatrixDyn(int _n, int _m)
	{
		n = _n;
		m = _m;
		mat = (T**)malloc(sizeof(T*) * n);
		for(int i=0;i<n;i++)
		{
			mat[i] = (T*)malloc(sizeof(T) * m);
			for(int j=0;j<m;j++)
				if(i==j) mat[i][j] = 1;
				else mat[i][j] = 0;
		}
	}

	GenMatrixDyn(T** m, int _n, int _m)
	{
		n = _n;
		m = _m;
		mat = m;
	}

	T& operator()(int x, int y){return mat[x][y];}
	
	void matMul(GenMatrixDyn<T>* mat2, GenMatrixDyn<T>* out);

	void transpose(GenMatrixDyn<T>* out);
	
	//GenMatrixDyn<T> inverse(GenMatrixDyn<T>* out);

	void matAddToThis(GenMatrixDyn<T>* mat2)
	{
		if(n == mat2->n && m == mat2->m)
		{
			for(int i=0;i<n;i++)
				for(int j=0;j<m;j++)
					mat[i][j] += mat2(i,j);
		}
	}

	void scalarMul(T k)
	{
		for(int i=0;i<N;i++)
			for(int j=0;j<M;j++)
				mat[i][j] *= k;
	}

	void vectorMul(GenVectorDyn<T>* vec, GenVectorDyn<T>* out)
	{
		if(m == vec->n && out->n == n)
		{
			
			for(int i=0;i<n;i++)
			{
				(*out)[i] = 0;
				for(int j=0;j<m;j++)
				{
					(*out)[i] += mat[i][j] * (*vec)[j];
				}
			}
		}
	}



};


template<class T>
void GenMatrixDyn<T>::matMul(GenMatrixDyn<T>* mat2, GenMatrixDyn<T>* out)
{
	if(out->n == n && out->m == mat2->m && m == mat2->n)
	{
		for(int i = 0; i<n; i++)
			for(int j = 0; j<mat2->m; j++)
			{
				out(i,j) = 0;
				for(int k = 0; k<m; k++)
					out(i,j) += mat[i][k] * mat2(k,j);
			}
	}
}


template<class T>
void GenMatrixDyn<T>::transpose(GenMatrixDyn<T>* out)
{
	if(out->m == n && out->n == m)
	{
		for(int i=0; i<n; i++)
			for(int j=0;j<m;j++)
				out(j,i) = mat[i][j];

	}
}

/*
template<class T>
void GenMatrixDyn<T>::inverse(GenMatrixDyn<T>* out)
{
	//////////////////////////////////////////
	//simulate find of inverse using a float width matrix in form
	//[A | I] -> [I | A^-1]
	//idea from http://programming-technique.blogspot.sg/2011/09/numerical-methods-inverse-of-nxn-matrix.html
	//////////////////////////////////////////

	if(N!=M) return GenMatrix<T,M,N>();
	else
	{
		T result[N][N * 2];

		//build [ A | I ]
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				result[i][j] = mat[i][j];

		for(int i=0;i<M;i++)
			for(int j=N;j<N * 2;j++)
				if(i==(j-N))
					result[i][j] = 1;
				else
					result[i][j] = 0;

		for(int i = 0; i < N; i++)
		{
			for(int j = 0; j < N; j++)
			{
				if(i!=j)
				{
					float ratio = result[j][i]/result[i][i];
					for(int k = 0; k < 2*n; k++)
						result[j][k] -= ratio * result[i][k];
				}
			}
		}

		for(int i = 0; i < n; i++)
		{
			float a = result[i][i];
			for(int j = 0; j < 2*n; j++)
				result[i][j] /= a;
		}

		return GenMatrix<T,M,N>(result);
	}
}
*/