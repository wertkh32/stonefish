#include "ConjugateGradientSolver.h"


ConjugateGradientSolver::ConjugateGradientSolver(int _n, float** _A)
{
	n = _n;
	A = _A;
	r = new float[n];
	x = new float[n];
	d = new float[n];
}

void
ConjugateGradientSolver::initSolver()
{
	
}

void
ConjugateGradientSolver::solve(float* x, float* b)
{
	float deltaOld, deltaNew, delta0,alpha,beta;
	float* q = new float[n];
	int it;
	
	for(int i=0;i<n;i++)
	{
		r[i]=0;
		d[i]=0;
		x[i]=0;
	}

	for(int i=0; i<n;i++)
	{
		d[i] = r[i] = b[i] - dot(A[i],x,n);
	}

	it=0;
	deltaNew = dot(r,r,n);
	delta0 = deltaNew;
	alpha = beta = 0;

	while(it < MAX_ITER && deltaNew > EPSILON*EPSILON*delta0)
	{
		for(int i=0;i<n;i++)
			q[i] = dot(A[i],d,n);

		alpha = deltaNew/dot(d,q,n);
		
		for(int i=0;i<n;i++)
			x[i] = x[i] + alpha*d[i];

		if(it%50==0)
		{
			//refresh r of its horrible floating point errors
			for(int i=0;i<n;i++)
				r[i] = b[i] - dot(A[i],x,n);
		}
		else
		{
			for(int i=0;i<n;i++)
				r[i] = r[i] - alpha * q[i];
		}

		deltaOld = deltaNew;
		deltaNew = dot(r,r,n);
		beta = deltaNew/deltaOld;

		for(int i=0;i<n;i++)
			d[i] = r[i] + beta * d[i];

		it++;
	}


}

ConjugateGradientSolver::~ConjugateGradientSolver(void)
{
}
