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
	double *r, *x, *d, *q;
	double** A;
	int n;
public:
	ConjugateGradientSolver(int _n, double** A);
	~ConjugateGradientSolver(void);
	void solve(double* x, double* b);
	void solveWithConstraints(double* x, double* b, ConstrainedRows* rowSet);

	double dot(double* a, double* b, int k)
	{
		double r=0;
		for(int i=0;i<k;i++)
			r+=a[i]*b[i];
		return r;
	}
	void initSolver();
	void removeRows(int r);
};

