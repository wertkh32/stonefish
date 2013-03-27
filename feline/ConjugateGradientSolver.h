#pragma once
#include "includes.h"

#define MAX_ITER 1000
#define EPSILON 0.001
//use pure arrays
class ConjugateGradientSolver
{
	double *r, *x, *d, *q;
	double** A;
	int n;
public:
	ConjugateGradientSolver(int _n, double** A);
	~ConjugateGradientSolver(void);
	void solve(double* x, double* b);
	double dot(double* a, double* b, int k)
	{
		double r=0;
		for(int i=0;i<k;i++)
			r+=a[i]*b[i];
		return r;
	}
	void initSolver();
};

