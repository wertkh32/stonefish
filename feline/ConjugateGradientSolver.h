#pragma once
#include "includes.h"

#define MAX_ITER 1000
#define EPSILON 0.001
//use pure arrays
class ConjugateGradientSolver
{
	float *r, *x, *d, *q;
	float** A;
	int n;
public:
	ConjugateGradientSolver(int _n, float** A);
	~ConjugateGradientSolver(void);
	void solve(float* x, float* b);
	float dot(float* a, float* b, int k)
	{
		float r=0;
		for(int i=0;i<k;i++)
			r+=a[i]*b[i];
		return r;
	}
	void initSolver();
};

