#pragma once

#include "includes.h"
#include "QuickArray.h"

#define MAX_ITER 1000
#define EPSILON 0.001

#define MAX_ROWS_CONSTRAINED 50

enum DOF
{
	X = 0, Y = 1, Z = 2
};

struct ConstrainedRows
{
	QuickArray<int,MAX_ROWS_CONSTRAINED> list;
	ConstrainedRows()
	{
	}
	
	void add(int node, DOF dof)
	{
		list.push(node * 3 + dof);
	}
	
	void add(int node)
	{
		list.push(node * 3);
		list.push(node * 3 + 1);
		list.push(node * 3 + 2);
	}


};

//use pure arrays
class ConjugateGradientSolver
{
	float *r, *x, *d, *q, *tempo;
	float *flatA;
	float** A;
	int n, nb;

	void sysMul(float* in, float* out, float** A, bool** matmap, bool* allowed)
	{
		for(int i=0;i<n;i++)
			out[i] = 0;

		for(int i=0;i<nb;i++)
			for(int j=0;j<nb;j++)
			{
				if(matmap[i][j])
				{
					for(int a=0;a<3;a++)
					{
						int xx = i * 3 + a;
						if(allowed[xx])
							for(int b=0;b<3;b++)
							{
								int yy = j * 3 + b;
								if(allowed[yy])
									out[xx] += A[xx][yy] * in[yy];
							}
					}
				}
			}
	}
public:
	ConjugateGradientSolver();
	~ConjugateGradientSolver(void);
	void initSolver(int _n, float** A);

	void solve(float* x, float* b);
	void solveWithConstraints(float* x, float* b, bool* allowed);
	void solveWithConstraints(float* x, float* b, bool* allowed, bool** matmap);

	float dot(float* a, float* b, int k)
	{
		float r=0;
		for(int i=0;i<k;i++)
			r+=a[i]*b[i];
		return r;
	}
	
	void removeRows(int r);
};

