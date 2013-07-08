#pragma once

#include "includes.h"
#include "Mesh.h"
#include "QuickArray.h"
#include "ConstrainedRows.h"

#define MAX_ITER 10
#define EPSILON 0.01



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


void sysMulMatFree(float* in, float* out, bool* allowed, Mesh* mesh)
{	
		for(int i=0;i<n;i++)
			out[i] = 0;

		//mesh->mulA(in,out);
		//mesh->mulRKRT(in,out);
		
		for(int i=0;i<mesh->elements.size();i++)
		{
			Element& ele = *(mesh->elements[i]);


			GenMatrix<float,12,12>& A = *(ele.getA());

			for(int a=0;a<4;a++)
			{
				int x = mesh->nodeIndices[i][a];
				for(int b=0;b<3;b++)
				{
					if(allowed[x*3+b])
					for(int c=0;c<4;c++)
					{
						    int y = mesh->nodeIndices[i][c];
							for(int d=0;d<3;d++)
							{
								if(allowed[y*3+d])
									out[x*3 + b] += A(a*3+b,c*3+d) * in[y*3 + d];
		
							}
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
	void solveWithConstraints(float* x, float* b, bool* allowed, Mesh* mesh);

	float dot(float* a, float* b, int k)
	{
		float r=0;
		for(int i=0;i<k;i++)
			r+=a[i]*b[i];
		return r;
	}
	
	void removeRows(int r);
};

