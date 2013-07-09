#pragma once
#include <cuda.h>
#include <stdio.h>
#include "defines.h"
#include "GPUDataStructs.cuh"
#include "GPUPolarDecompose.cu"

#define BLOCK_SIZE 256
#define ALPHA 0.1
#define BETA 0.1

__constant__ float COEFFK, COEFFM, dt;

GPUElement* gpuptrElements;
GPUNode*   gpuptrNodes;
mulData*	gpuptrMulData;
float*   gpuptr_x0;//const
float*   gpuptr_xt;//dynamic
float*   gpuptr_vt;//dynamic
float*	 gpuptr_extforces;//dynamic


void
gpuInitVars(int numele, int numnodes)
{
	cudaMalloc(&gpuptrElements, numele * sizeof(GPUElement));
	cudaMalloc(&gpuptrNodes, numnodes * sizeof(GPUNode));
	cudaMalloc(&gpuptrMulData, numele * sizeof(mulData));
	cudaMalloc(&gpuptr_xt, numnodes * 3 * sizeof(float));
	cudaMalloc(&gpuptr_vt, numnodes * 3 * sizeof(float));
	cudaMalloc(&gpuptr_extforces, numnodes * 3 * sizeof(float));

	float coeffK = dt * BETA + dt * dt, coeffM = 1 + dt * ALPHA;
	float dt = DT;

	cudaMemcpyToSymbol("COEFFK", &coeffK, sizeof(float));
	cudaMemcpyToSymbol("COEFFM", &coeffM, sizeof(float));
	cudaMemcpyToSymbol("dt", &dt, sizeof(float));
}

void
gpuUploadExtForces(float* extforces, int numnodes)
{
	cudaMemcpy(gpuptr_extforces, extforces, numnodes*3*sizeof(float),cudaMemcpyHostToDevice);
}

__device__
void makeRK(float mat[12][12], float R[3][3])
{
	float RK[12][12];
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			for(int a=0;a<3;a++)
			{
				for(int b=0;b<3;b++)
				{
					RK[a + i * 3][b + j * 3]=0;
						
					for(int c=0;c<3;c++)
						RK[a + i * 3][b + j * 3] += R[a][c] * mat[c + i * 3][b + j * 3];
				}
			}
		}
	}

	for(int i=0;i<12;i++)
		for(int j=0;j<12;j++)
			mat[i][j] = RK[i][j];
}


__device__
void makeRKRT(float mat[12][12], float R[3][3])
{
	float RKRT[12][12];

	for(int i=0;i<4;i++)
	{
		for(int j=i;j<4;j++)
		{
			for(int a=0;a<3;a++)
				for(int b=0;b<3;b++)
				{
					RKRT[a + i * 3][b + j * 3] = 0;
					for(int c=0;c<3;c++)
						RKRT[a + i * 3][b + j * 3] += mat[a + i * 3][c + j * 3] * R[b][c]; //R(c, b) but its RT so, R(b , c)
				}
		}
	}
	//lower triangle
	for(int i=1;i<4;i++)
		for(int j=0;j<i;j++)
		{
			for(int a=0;a<3;a++)
				for(int b=0;b<3;b++)
						RKRT[a + i * 3][b + j * 3] = RKRT[b + j * 3][a + i * 3];
		}

	for(int i=0;i<12;i++)
		for(int j=0;j<12;j++)
			mat[i][j] = RKRT[i][j];
}


__global__
void precompute(GPUElement* elements, mulData* solverData, float* xt, float* vt, float* extforces)
{
	int tid = threadIdx.x + blockIdx.x * BLOCK_SIZE;
	GPUElement* t_ele = &(elements[tid]);
	mulData* t_solvedata = &(solverData[tid]);

	float nodalmass = t_ele->nodalmass;

	float nodes[12], vel[12], F[3][3], R[3][3];

	for(int i=0;i<4;i++)
	{
		nodes[i * 3] = xt[t_ele->nodeindex[i] * 3];
		nodes[i * 3 + 1] = xt[t_ele->nodeindex[i] * 3 + 1];
		nodes[i * 3 + 2] = xt[t_ele->nodeindex[i] * 3 + 2];
	}

	for(int i=0;i<4;i++)
	{
		vel[i * 3] = vt[t_ele->nodeindex[i] * 3];
		vel[i * 3 + 1] = vt[t_ele->nodeindex[i] * 3 + 1];
		vel[i * 3 + 2] = vt[t_ele->nodeindex[i] * 3 + 2];
	}

	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
		{
			F[i][j] = nodes[j*3 + i] - nodes[9 + i];
		}

	gpuComputePolarDecomposition(F,R);


	for(int i=0;i<12;i++)
		for(int j=0;j<12;j++)
			t_solvedata->system[i][j] = t_ele->unwarpK[i][j];
	
	for(int i=0;i<4;i++)
	{
		t_ele->b[i * 3] = extforces[t_ele->nodeindex[i] * 3];
		t_ele->b[i * 3 + 1] = extforces[t_ele->nodeindex[i] * 3 + 1];
		t_ele->b[i * 3 + 2] = extforces[t_ele->nodeindex[i] * 3 + 2];
	}

	makeRK(t_solvedata->system, R);

	for(int i=0;i<12;i++)
		for(int j=0;j<12;j++)
			t_ele->b[i] += t_solvedata->system[i][j] * t_ele->x0[j];
	
	makeRKRT(t_solvedata->system, R);

	for(int i=0;i<12;i++)
		for(int j=0;j<12;j++)
			t_ele->b[i] -= t_solvedata->system[i][j] * nodes[j];

	for(int i=0;i<12;i++)
		t_ele->b[i] = t_ele->b[i] * dt + nodalmass * vel[i];

	//final system matrix
	for(int i=0;i<12;i++)
		for(int j=0;j<12;j++)
		{
			t_solvedata->system[i][j] *= COEFFK;
			if(i==j)
				t_solvedata->system[i][i] += COEFFM * nodalmass;
		}
}

__global__
void collateB(GPUNode* nodes, GPUElement* elements, float* b)
{
	int tid = threadIdx.x + blockIdx.x * BLOCK_SIZE;
	GPUNode* node = &(nodes[tid]);

	int n = node->n;

	b[tid * 3] = 0;
	b[tid * 3 + 1] = 0;
	b[tid * 3 + 2] = 0;

	for(int i=0;i<n;i++)
	{
		int tetindex = node->elementindex[i][0];
		int nodeindex = node->elementindex[i][1];

		b[tid * 3] += elements[tetindex].b[nodeindex * 3];
		b[tid * 3 + 1] += elements[tetindex].b[nodeindex * 3 + 1];
		b[tid * 3 + 2] += elements[tetindex].b[nodeindex * 3 + 2];
	}
}

__global__
void mulSystem(GPUElement* elements, mulData* solverData, float* x)
{
	int tid = threadIdx.x + blockIdx.x * BLOCK_SIZE;
	GPUElement* t_ele = &(elements[tid]);
	mulData* t_solvedata = &(solverData[tid]);

	float nodes[12];

	for(int i=0;i<4;i++)
	{
		nodes[i * 3] = x[t_ele->nodeindex[i] * 3];
		nodes[i * 3 + 1] = x[t_ele->nodeindex[i] * 3 + 1];
		nodes[i * 3 + 2] = x[t_ele->nodeindex[i] * 3 + 2];
	}

	for(int i=0;i<12;i++)
	{
		t_solvedata->product[i] = 0;
		for(int j=0;j<12;j++)
			t_solvedata->product[i] += t_solvedata->system[i][j] * nodes[j];
	}
}

__global__ 
void mulSystemGather(GPUNode* nodes, mulData* solverData, float* x)
{
	int tid = threadIdx.x + blockIdx.x * BLOCK_SIZE;
	GPUNode* node = &(nodes[tid]);

	int n = node->n;
	
	x[tid * 3] = 0;
	x[tid * 3 + 1] = 0;
	x[tid * 3 + 2] = 0;

	for(int i=0;i<n;i++)
	{
		int tetindex = node->elementindex[i][0];
		int nodeindex = node->elementindex[i][1];

		x[tid * 3] += solverData[tetindex].product[nodeindex * 3];
		x[tid * 3 + 1] += solverData[tetindex].product[nodeindex * 3 + 1];
		x[tid * 3 + 2] += solverData[tetindex].product[nodeindex * 3 + 2];
	}
}

void
timestep()
{
	
}

