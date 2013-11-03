#pragma once

#include "includes.h"
#include "Mesh.h"
#include "QuadTetMesh.h"
#include "QuickArray.h"
#include "ConstrainedRows.h"

#define MAX_ITER 20
#define EPSILON 0.05


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


void sysMulMatFree(float* in, float* out, bool* allowed, MESH* mesh)
{	
	static const float alpha = 0.1, beta = 0.1;
	static const float coeffK = (1.0/FPS) * beta + (1.0/FPS) * (1.0/FPS), coeffM = 1 + (1.0/FPS) * alpha;

	for(int i=0;i<n;i++)
	{
		out[i] = 0;
		if(!allowed[i])
			in[i] = 0;
	}
		//mesh->mulA(in,out);
		//mesh->mulRKRT(in,out);
		
	for(int i=0;i<mesh->elements.size();i++)
		{
			ELEMENT& ele = *(mesh->elements[i]);
			GenMatrix<float,NUM_NODES_PER_ELE * 3,NUM_NODES_PER_ELE * 3>& A = ele.getStiffnessMat();
			Matrix3d& R = ele.getRotation();

			//float mass = ele.nodalMass * coeffM;

			float pos[NUM_NODES_PER_ELE * 3] = {0};
			float pos2[NUM_NODES_PER_ELE * 3] = {0};
			float pos3[NUM_NODES_PER_ELE * 3] = {0};

			for(int a=0;a<NUM_NODES_PER_ELE;a++)
			{
				int x = mesh->nodeIndices[i][a];

				pos[a * 3] = in[x * 3];
				pos[a * 3 + 1] = in[x * 3 + 1];
				pos[a * 3 + 2] = in[x * 3 + 2];
			}

			for(int a=0;a<NUM_NODES_PER_ELE;a++)
			{
				for(int b=0;b<3;b++)
					for(int c=0;c<3;c++)
						pos2[a*3+b] += R(c,b) * pos[a*3+c];
			}

			for(int a=0;a<NUM_NODES_PER_ELE * 3;a++)
				for(int b=0;b<NUM_NODES_PER_ELE * 3;b++)
					pos3[a] += A(a,b) * pos2[b];

			for(int a=0;a<NUM_NODES_PER_ELE;a++)
			{
				pos2[a*3] = 0;
				pos2[a*3+1] = 0;
				pos2[a*3+2] = 0;

				for(int b=0;b<3;b++)
					for(int c=0;c<3;c++)
						pos2[a*3+b] += R(b,c) * pos3[a*3+c];
			}

			for(int a=0;a<NUM_NODES_PER_ELE;a++)
			{
				int x = mesh->nodeIndices[i][a];

				out[x * 3] += pos2[a * 3] * coeffK + ele.nodemass[a] * coeffM * pos[a * 3];
				out[x * 3 + 1] += pos2[a * 3 + 1] * coeffK + ele.nodemass[a] * coeffM * pos[a * 3 + 1];
				out[x * 3 + 2] += pos2[a * 3 + 2] * coeffK + ele.nodemass[a] * coeffM * pos[a * 3 + 2];
			}

		}
		
		for(int i=0;i<n;i++)
		{
			if(!allowed[i])
				out[i] = 0;
		}
		
}

public:
	ConjugateGradientSolver();
	~ConjugateGradientSolver(void);
	void initSolver(int _n, float** A);

	void solve(float* x, float* b);
	void solveWithConstraints(float* x, float* b, bool* allowed);
	void solveWithConstraints(float* x, float* b, bool* allowed, bool** matmap);
	void solveWithConstraints(float* x, float* b, bool* allowed, MESH* mesh);

	float dot(float* a, float* b, int k)
	{
		float r=0;
		for(int i=0;i<k;i++)
			r+=a[i]*b[i];
		return r;
	}
	
	void removeRows(int r);
};

