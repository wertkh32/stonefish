#include "ConjugateGradientSolver.h"


ConjugateGradientSolver::ConjugateGradientSolver(int _n, float** _A)
{
	n = _n;
	A = _A;
	r = (float*)malloc(sizeof(float) * n); 
	x = (float*)malloc(sizeof(float) * n); 
	d = (float*)malloc(sizeof(float) * n); 
	q = (float*)malloc(sizeof(float) * n); 
}

void
ConjugateGradientSolver::initSolver()
{
	
}

void
ConjugateGradientSolver::solve(float* x, float* b)
{
	float deltaOld, deltaNew, delta0,alpha,beta;
	int it;
	
	for(int i=0;i<n;i++)
	{
		r[i]=0;
		d[i]=0;
		//x[i]=0;
		q[i]=0;
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

		if(it%10==0)
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
